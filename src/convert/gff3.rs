use crate::parser::gff::{GFFReader, GFFRecord};
/// GFF3 conversion utilities
///
/// Supported output formats:
///   GTF        — GENCODE/Ensembl gene transfer format
///   BED6       — 6-column BED for genome browsers
///   BED12      — 12-column BED with exon block structure per transcript
///   BEDGraph   — per-feature coverage depth track
///   TSV table  — flat feature table (all attributes expanded)
///   Protein    — translated CDS sequences as FASTA (.faa)
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 → GTF
// ─────────────────────────────────────────────────────────────────────────────

/// Convert a GFF3 file to GTF (GENCODE-style).
///
/// Handles the two key differences between formats:
///   1. Feature type renaming (`mRNA` → `transcript`)
///   2. Attribute serialisation (`key=value` → `key "value";`)
///      with mandatory `gene_id` and `transcript_id` fields on every line
pub fn gff3_to_gtf(input: &Path, output: &Path) -> Result<usize> {
    let records: Vec<GFFRecord> = GFFReader::from_path(input)?
        .filter_map(|r| r.ok())
        .collect();

    // Build ID → Parent map so we can resolve gene_id / transcript_id
    // for grandchild features (exon, CDS, UTR)
    let mut id_to_parent: HashMap<String, String> = HashMap::new();
    for rec in &records {
        if let (Some(id), Some(parent)) = (rec.id(), rec.parent()) {
            id_to_parent.insert(id.clone(), parent.clone());
        }
    }

    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(out, "##gtf-version 2").map_err(MycoNoteError::Io)?;
    let mut count = 0;

    for rec in &records {
        let gtf_feature = match rec.feature_type.as_str() {
            "mRNA" | "transcript" => "transcript",
            "five_prime_UTR" => "5UTR",
            "three_prime_UTR" => "3UTR",
            "region" | "chromosome" => continue, // skip assembly regions
            other => other,
        };

        let id = rec.id().map(|s| s.as_str()).unwrap_or("");
        let parent = rec.parent().map(|s| s.as_str()).unwrap_or("");

        // Resolve gene_id and transcript_id for every feature type
        let (gene_id, transcript_id) = match rec.feature_type.as_str() {
            "gene" => (id, ""),
            "mRNA" | "transcript" => (parent, id),
            _ => {
                // exon / CDS / UTR: parent = mRNA, grandparent = gene
                let tid = parent;
                let gid = id_to_parent.get(parent).map(|s| s.as_str()).unwrap_or("");
                (gid, tid)
            }
        };

        // Build attribute string
        let mut attrs = format!(
            "gene_id \"{}\"; transcript_id \"{}\";",
            gene_id, transcript_id
        );
        if let Some(name) = rec.attributes.get("Name") {
            attrs.push_str(&format!(" gene_name \"{}\";", name));
        }
        if let Some(biotype) = rec
            .attributes
            .get("gene_biotype")
            .or_else(|| rec.attributes.get("biotype"))
        {
            attrs.push_str(&format!(" gene_biotype \"{}\";", biotype));
        }
        if let Some(product) = rec.attributes.get("product") {
            attrs.push_str(&format!(" product \"{}\";", product));
        }

        let score = rec
            .score
            .map(|s| format!("{}", s))
            .unwrap_or_else(|| ".".into());
        let phase = rec
            .phase
            .map(|p| format!("{}", p))
            .unwrap_or_else(|| ".".into());

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.seqid, rec.source, gtf_feature, rec.start, rec.end, score, rec.strand, phase, attrs
        )
        .map_err(MycoNoteError::Io)?;

        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 → BED6
// ─────────────────────────────────────────────────────────────────────────────

/// Convert GFF3 features to 6-column BED.
///
/// If `feature_types` is empty, all features are converted.
/// Coordinates are converted from GFF3 1-based inclusive to BED 0-based half-open.
pub fn gff3_to_bed(input: &Path, output: &Path, feature_types: &[&str]) -> Result<usize> {
    let filter: std::collections::HashSet<&str> = feature_types.iter().copied().collect();
    let filter_all = filter.is_empty();

    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0;

    for rec_res in GFFReader::from_path(input)? {
        let rec = match rec_res {
            Ok(r) => r,
            Err(_) => continue,
        };

        if !filter_all && !filter.contains(rec.feature_type.as_str()) {
            continue;
        }

        let name = rec
            .id()
            .or_else(|| rec.attributes.get("Name"))
            .map(|s| s.as_str())
            .unwrap_or(".");
        let score = rec
            .score
            .map(|s| format!("{:.0}", s))
            .unwrap_or_else(|| "0".into());

        // GFF3: 1-based inclusive → BED: 0-based half-open
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            rec.seqid,
            rec.start - 1,
            rec.end,
            name,
            score,
            rec.strand,
        )
        .map_err(MycoNoteError::Io)?;

        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 → BED12
// ─────────────────────────────────────────────────────────────────────────────

/// Convert GFF3 transcripts + their exons to BED12 format.
/// One output line per mRNA/transcript, blocks derived from exon children.
pub fn gff3_to_bed12(input: &Path, output: &Path) -> Result<usize> {
    let records: Vec<GFFRecord> = GFFReader::from_path(input)?
        .filter_map(|r| r.ok())
        .collect();

    let mut transcripts: HashMap<String, GFFRecord> = HashMap::new();
    let mut tx_exons: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    let mut tx_cds: HashMap<String, (u64, u64)> = HashMap::new();

    for rec in &records {
        match rec.feature_type.as_str() {
            "mRNA" | "transcript" => {
                if let Some(id) = rec.id() {
                    transcripts.insert(id.clone(), rec.clone());
                }
            }
            "exon" => {
                if let Some(p) = rec.parent() {
                    tx_exons
                        .entry(p.clone())
                        .or_default()
                        .push((rec.start, rec.end));
                }
            }
            "CDS" => {
                if let Some(p) = rec.parent() {
                    let e = tx_cds.entry(p.clone()).or_insert((rec.start, rec.end));
                    e.0 = e.0.min(rec.start);
                    e.1 = e.1.max(rec.end);
                }
            }
            _ => {}
        }
    }

    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0;

    let mut ids: Vec<&String> = transcripts.keys().collect();
    ids.sort();

    for id in ids {
        let tx = &transcripts[id];
        let chrom_start = tx.start - 1; // 0-based
        let chrom_end = tx.end;

        let (thick_start, thick_end) = tx_cds
            .get(id)
            .map(|&(s, e)| (s - 1, e))
            .unwrap_or((chrom_start, chrom_start)); // no CDS → thickStart = thickEnd

        let name = tx
            .attributes
            .get("Name")
            .map(|s| s.as_str())
            .unwrap_or(id.as_str());
        let score = tx
            .score
            .map(|s| format!("{:.0}", s))
            .unwrap_or_else(|| "0".into());

        // Build exon blocks (0-based coordinates, sorted)
        let blocks: Vec<(u64, u64)> = tx_exons
            .get(id)
            .map(|exons| {
                let mut b: Vec<(u64, u64)> = exons.iter().map(|&(s, e)| (s - 1, e)).collect();
                b.sort_by_key(|&(s, _)| s);
                b
            })
            .unwrap_or_else(|| vec![(chrom_start, chrom_end)]);

        let block_count = blocks.len();
        let block_sizes: Vec<String> = blocks.iter().map(|(s, e)| (e - s).to_string()).collect();
        let block_starts: Vec<String> = blocks
            .iter()
            .map(|(s, _)| (s - chrom_start).to_string())
            .collect();

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t{},\t{},",
            tx.seqid,
            chrom_start,
            chrom_end,
            name,
            score,
            tx.strand,
            thick_start,
            thick_end,
            block_count,
            block_sizes.join(","),
            block_starts.join(","),
        )
        .map_err(MycoNoteError::Io)?;

        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 → BEDGraph (feature density track)
// ─────────────────────────────────────────────────────────────────────────────

/// Produce a BEDGraph track with value = 1 for each feature interval.
/// Useful for visualising gene density in IGV or UCSC browser.
pub fn gff3_to_bedgraph(input: &Path, output: &Path, feature_type: &str) -> Result<usize> {
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(
        out,
        "track type=bedGraph name=\"{}\" visibility=full",
        feature_type
    )
    .map_err(MycoNoteError::Io)?;
    let mut count = 0;

    for rec_res in GFFReader::from_path(input)? {
        let rec = match rec_res {
            Ok(r) => r,
            Err(_) => continue,
        };
        if rec.feature_type != feature_type {
            continue;
        }

        writeln!(out, "{}\t{}\t{}\t1", rec.seqid, rec.start - 1, rec.end)
            .map_err(MycoNoteError::Io)?;
        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 → Feature table (TSV)
// ─────────────────────────────────────────────────────────────────────────────

/// Flatten all GFF3 features into a tab-separated table.
/// Standard columns + all attribute key=value pairs as the last column.
pub fn gff3_to_table(input: &Path, output: &Path) -> Result<usize> {
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(
        out,
        "seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tID\tParent\tName\tattributes"
    )
    .map_err(MycoNoteError::Io)?;

    let mut count = 0;
    for rec_res in GFFReader::from_path(input)? {
        let rec = match rec_res {
            Ok(r) => r,
            Err(_) => continue,
        };

        let id = rec.attributes.get("ID").map(|s| s.as_str()).unwrap_or("");
        let parent = rec
            .attributes
            .get("Parent")
            .map(|s| s.as_str())
            .unwrap_or("");
        let name = rec.attributes.get("Name").map(|s| s.as_str()).unwrap_or("");

        let extra: Vec<String> = rec
            .attributes
            .iter()
            .filter(|(k, _)| !matches!(k.as_str(), "ID" | "Parent" | "Name"))
            .map(|(k, v)| format!("{}={}", k, v))
            .collect();

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.seqid,
            rec.source,
            rec.feature_type,
            rec.start,
            rec.end,
            rec.score
                .map(|s| format!("{}", s))
                .unwrap_or_else(|| ".".into()),
            rec.strand,
            rec.phase
                .map(|p| format!("{}", p))
                .unwrap_or_else(|| ".".into()),
            id,
            parent,
            name,
            extra.join(";"),
        )
        .map_err(MycoNoteError::Io)?;

        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 → Protein FASTA
// ─────────────────────────────────────────────────────────────────────────────

/// Extract translated CDS sequences for all mRNA features.
/// Requires the paired genome FASTA.
pub fn gff3_to_protein(gff_path: &Path, fasta_path: &Path, output: &Path) -> Result<usize> {
    use crate::parser::fasta::{read_fasta_index, reverse_complement};

    let fasta_index = read_fasta_index(fasta_path)?;
    let records: Vec<GFFRecord> = GFFReader::from_path(gff_path)?
        .filter_map(|r| r.ok())
        .collect();

    let mut mrna_cds: HashMap<String, Vec<GFFRecord>> = HashMap::new();
    let mut mrna_to_gene: HashMap<String, String> = HashMap::new();

    for rec in &records {
        match rec.feature_type.as_str() {
            "CDS" => {
                if let Some(p) = rec.parent() {
                    mrna_cds.entry(p.clone()).or_default().push(rec.clone());
                }
            }
            "mRNA" | "transcript" => {
                if let (Some(id), Some(p)) = (rec.id(), rec.parent()) {
                    mrna_to_gene.insert(id.clone(), p.clone());
                }
            }
            _ => {}
        }
    }

    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0;

    for (mrna_id, mut cds_list) in mrna_cds {
        cds_list.sort_by_key(|r| r.start);

        let seq_rec = match fasta_index.get(&cds_list[0].seqid) {
            Some(s) => s,
            None => continue,
        };

        let strand = cds_list[0].strand;
        let mut cds_seq = String::new();
        for cds in &cds_list {
            cds_seq.push_str(seq_rec.subsequence(cds.start, cds.end));
        }
        if strand == '-' {
            cds_seq = reverse_complement(&cds_seq);
        }

        let protein = translate_dna(&cds_seq);
        if protein.len() < 10 {
            continue;
        }

        let gene_id = mrna_to_gene
            .get(&mrna_id)
            .map(|s| s.as_str())
            .unwrap_or(mrna_id.as_str());
        writeln!(out, ">{} gene={}", mrna_id, gene_id).map_err(MycoNoteError::Io)?;
        for chunk in protein.as_bytes().chunks(60) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                .map_err(MycoNoteError::Io)?;
        }
        count += 1;
    }

    if count == 0 {
        return Err(MycoNoteError::InvalidFormat(
            "No proteins written: check that GFF3 seqids match FASTA headers \
             (seqid must equal the first whitespace-delimited token of the '>' line)"
                .to_string(),
        ));
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// Translation helper (shared with annotate module)
// ─────────────────────────────────────────────────────────────────────────────

pub fn translate_dna(dna: &str) -> String {
    let bytes = dna.as_bytes();
    let mut prot = String::with_capacity(bytes.len() / 3);
    let mut i = 0;
    while i + 2 < bytes.len() {
        let codon = [
            bytes[i].to_ascii_uppercase(),
            bytes[i + 1].to_ascii_uppercase(),
            bytes[i + 2].to_ascii_uppercase(),
        ];
        prot.push(codon_to_aa(&codon));
        i += 3;
    }
    if prot.ends_with('*') {
        prot.pop();
    }
    prot
}

fn codon_to_aa(c: &[u8; 3]) -> char {
    match c {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}
