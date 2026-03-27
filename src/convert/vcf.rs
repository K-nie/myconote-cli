/// VCF conversion utilities
///
/// VCF → BED intervals, TSV table, consensus FASTA, ANNOVAR input, MAF
/// Handles plain .vcf and gzip-compressed .vcf.gz files.

use crate::utils::error::{MycoNoteError, Result};
use std::io::{BufRead, Write};
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// Shared VCF line reader (handles plain and .gz)
// ─────────────────────────────────────────────────────────────────────────────

fn open_vcf(path: &Path) -> Result<Box<dyn BufRead>> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;

    if ext == "gz" || ext == "bgz" {
        use std::io::BufReader;
        // Use flate2 to decompress
        let decoder = flate2::read::GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(std::io::BufReader::new(file)))
    }
}

/// A single parsed VCF variant record.
#[derive(Debug)]
pub struct VcfRecord {
    pub chrom:  String,
    pub pos:    u64,        // 1-based
    pub id:     String,
    pub ref_:   String,
    pub alt:    String,     // first ALT allele
    pub qual:   Option<f64>,
    pub filter: String,
    pub info:   std::collections::HashMap<String, String>,
    pub samples: Vec<String>,
}

impl VcfRecord {
    fn from_line(line: &str) -> Option<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 { return None; }

        let pos: u64 = f[1].parse().ok()?;
        let qual = if f[5] == "." { None } else { f[5].parse().ok() };

        // Parse INFO field into a HashMap
        let mut info = std::collections::HashMap::new();
        if f[7] != "." {
            for pair in f[7].split(';') {
                let mut kv = pair.splitn(2, '=');
                let k = kv.next().unwrap_or("").to_string();
                let v = kv.next().unwrap_or("true").to_string();
                info.insert(k, v);
            }
        }

        let samples: Vec<String> = f[9..].iter().map(|s| s.to_string()).collect();

        Some(VcfRecord {
            chrom:   f[0].to_string(),
            pos,
            id:      f[2].to_string(),
            ref_:    f[3].to_string(),
            alt:     f[4].split(',').next().unwrap_or(".").to_string(),
            qual,
            filter:  f[6].to_string(),
            info,
            samples,
        })
    }

    /// End position in VCF coordinates (1-based)
    pub fn end_pos(&self) -> u64 {
        self.pos + self.ref_.len() as u64 - 1
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VCF → BED
// ─────────────────────────────────────────────────────────────────────────────

/// Convert VCF variant positions to BED intervals.
/// SNPs become 1-bp intervals; indels span the REF allele length.
pub fn vcf_to_bed(input: &Path, output: &Path) -> Result<usize> {
    let reader = open_vcf(input)?;
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0;

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') { continue; }
        let rec = match VcfRecord::from_line(&line) { Some(r) => r, None => continue };

        // VCF POS is 1-based → BED chromStart is 0-based
        let bed_start = rec.pos - 1;
        let bed_end   = rec.pos + rec.ref_.len() as u64 - 1;
        let name = if rec.id == "." {
            format!("{}_{}_{}/{}", rec.chrom, rec.pos, rec.ref_, rec.alt)
        } else {
            rec.id.clone()
        };

        writeln!(out, "{}\t{}\t{}\t{}\t.\t+",
            rec.chrom, bed_start, bed_end, name
        ).map_err(MycoNoteError::Io)?;
        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// VCF → TSV table
// ─────────────────────────────────────────────────────────────────────────────

/// Flatten a VCF file into a TSV table.
/// Standard columns + all INFO key=value pairs as named columns.
pub fn vcf_to_table(input: &Path, output: &Path) -> Result<usize> {
    // First pass: collect all INFO keys to build header
    let mut info_keys: Vec<String> = Vec::new();
    {
        let reader = open_vcf(input)?;
        for line in reader.lines() {
            let line = line.map_err(MycoNoteError::Io)?;
            if line.starts_with("##INFO=<ID=") {
                if let Some(id_start) = line.find("ID=") {
                    let rest = &line[id_start + 3..];
                    let id: String = rest.chars().take_while(|&c| c != ',').collect();
                    if !info_keys.contains(&id) { info_keys.push(id); }
                }
            }
            if !line.starts_with('#') { break; }
        }
    }

    let reader = open_vcf(input)?;
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;

    // Header
    let mut header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER".to_string();
    for k in &info_keys { header.push('\t'); header.push_str(k); }
    writeln!(out, "{}", header).map_err(MycoNoteError::Io)?;

    let mut count = 0;
    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') { continue; }
        let rec = match VcfRecord::from_line(&line) { Some(r) => r, None => continue };

        let mut row = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.chrom, rec.pos, rec.id, rec.ref_, rec.alt,
            rec.qual.map(|q| format!("{:.2}", q)).unwrap_or_else(|| ".".into()),
            rec.filter,
        );
        for k in &info_keys {
            row.push('\t');
            row.push_str(rec.info.get(k).map(|s| s.as_str()).unwrap_or("."));
        }
        writeln!(out, "{}", row).map_err(MycoNoteError::Io)?;
        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// VCF → Consensus FASTA
// ─────────────────────────────────────────────────────────────────────────────

/// Apply SNP/indel variants to a reference FASTA to produce a consensus sequence.
/// Only handles haploid calls (takes the first ALT allele for any variant).
pub fn vcf_to_consensus(vcf_path: &Path, ref_fasta: &Path, output: &Path) -> Result<usize> {
    use crate::parser::fasta::read_fasta;

    let references = read_fasta(ref_fasta)?;
    let reader = open_vcf(vcf_path)?;

    // Load all variants per chromosome
    let mut variants: std::collections::HashMap<String, Vec<(u64, String, String)>> =
        std::collections::HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') { continue; }
        if let Some(rec) = VcfRecord::from_line(&line) {
            if rec.alt != "." && rec.alt != rec.ref_ {
                variants.entry(rec.chrom.clone())
                    .or_default()
                    .push((rec.pos, rec.ref_.clone(), rec.alt.clone()));
            }
        }
    }

    // Sort variants by position (descending so substitutions don't shift coordinates)
    for v in variants.values_mut() {
        v.sort_by_key(|(pos, _, _)| *pos);
        v.reverse();
    }

    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0;

    for rec in &references {
        let mut seq = rec.sequence.clone();

        if let Some(vars) = variants.get(&rec.id) {
            for (pos, ref_allele, alt_allele) in vars {
                let s = (*pos as usize).saturating_sub(1);  // 0-based
                let e = s + ref_allele.len();
                if e <= seq.len() {
                    seq.replace_range(s..e, alt_allele.to_uppercase().as_str());
                }
            }
        }

        writeln!(out, ">{} consensus", rec.id).map_err(MycoNoteError::Io)?;
        for chunk in seq.as_bytes().chunks(60) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                .map_err(MycoNoteError::Io)?;
        }
        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// VCF → ANNOVAR input format
// ─────────────────────────────────────────────────────────────────────────────

/// Convert VCF to ANNOVAR input format (avinput).
/// Columns: Chr  Start  End  Ref  Alt
pub fn vcf_to_annovar(input: &Path, output: &Path) -> Result<usize> {
    let reader = open_vcf(input)?;
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0;

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') { continue; }
        let rec = match VcfRecord::from_line(&line) { Some(r) => r, None => continue };

        // ANNOVAR uses 1-based coordinates (same as VCF)
        let _end = rec.pos + rec.ref_.len() as u64 - 1;

        // Handle indels: ANNOVAR convention
        let (start, end_out, ref_out, alt_out) = if rec.ref_.len() == 1 && rec.alt.len() == 1 {
            // SNV
            (rec.pos, rec.pos, rec.ref_.clone(), rec.alt.clone())
        } else if rec.ref_.len() > rec.alt.len() {
            // Deletion: trim leading common base
            let r = rec.ref_[1..].to_string();
            let a = if rec.alt.len() > 1 { rec.alt[1..].to_string() } else { "-".to_string() };
            (rec.pos + 1, rec.pos + r.len() as u64, r, a)
        } else {
            // Insertion
            let a = rec.alt[1..].to_string();
            (rec.pos, rec.pos + 1, "-".to_string(), a)
        };

        writeln!(out, "{}\t{}\t{}\t{}\t{}", rec.chrom, start, end_out, ref_out, alt_out)
            .map_err(MycoNoteError::Io)?;
        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// VCF → MAF (Mutation Annotation Format)
// ─────────────────────────────────────────────────────────────────────────────

/// Convert VCF to MAF format (used by TCGA and cBioPortal).
/// Produces a simplified MAF with the core required columns.
pub fn vcf_to_maf(input: &Path, output: &Path, tumor_sample: &str) -> Result<usize> {
    let reader = open_vcf(input)?;
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;

    // MAF header
    writeln!(out,
        "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\t\
         Start_Position\tEnd_Position\tStrand\tVariant_Classification\t\
         Variant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\t\
         Tumor_Sample_Barcode"
    ).map_err(MycoNoteError::Io)?;

    let mut count = 0;
    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') { continue; }
        let rec = match VcfRecord::from_line(&line) { Some(r) => r, None => continue };

        let variant_type = classify_variant(&rec.ref_, &rec.alt);
        let end_pos = rec.pos + rec.ref_.len() as u64 - 1;

        writeln!(out,
            "Unknown\t0\t.\tGRCh38\t{}\t{}\t{}\t+\tMissense_Mutation\t{}\t{}\t{}\t{}\t{}",
            rec.chrom, rec.pos, end_pos,
            variant_type,
            rec.ref_, rec.ref_, rec.alt,
            tumor_sample,
        ).map_err(MycoNoteError::Io)?;
        count += 1;
    }

    Ok(count)
}

fn classify_variant(ref_: &str, alt: &str) -> &'static str {
    match (ref_.len(), alt.len()) {
        (1, 1) => "SNP",
        (1, n) if n > 1 => "INS",
        (n, 1) if n > 1 => "DEL",
        _ => "DNP",
    }
}
