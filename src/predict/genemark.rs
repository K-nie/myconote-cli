/// GeneMark-ES / ET / EP+ / ETP+ wrapper
///
/// Four self-training modes, all driven by `gmes_petap.pl`:
///
///   ES   self-training only (no extrinsic evidence)        — `--ES`
///   ET   RNA-seq intron hints (HISAT2/STAR splice-GFF)     — `--ET`
///   EP+  protein-evidence hints (built by ProtHint)        — `--EP`
///   ETP+ both protein and RNA-seq hints                    — `--ETP`
///
/// EP+ and ETP+ require ProtHint to convert a protein FASTA + genome into
/// the GFF hint files that GeneMark consumes.  ProtHint ships bundled with
/// the GeneMark-ES installer tarball (`<install>/ProtHint/bin/prothint.py`)
/// and with the bioconda `braker3` package; there is no standalone bioconda
/// recipe.  See `predict.md` for installation guidance.
///
/// GeneMark itself requires a free academic licence from Georgia Tech:
///   http://topaz.gatech.edu/GeneMark/
use crate::utils::error::{MycoNoteError, Result};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Mode enum
// ─────────────────────────────────────────────────────────────────────────────

/// Which GeneMark training algorithm to run.  Default is `Es`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneMarkMode {
    /// Self-training, no extrinsic evidence.
    Es,
    /// RNA-seq intron hints (`--ET`).
    Et,
    /// Protein-evidence hints from ProtHint (`--EP`, EP+ algorithm).
    Ep,
    /// Combined protein + RNA-seq hints (`--ETP`, ETP+ algorithm).
    Etp,
}

impl GeneMarkMode {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "es" => Some(Self::Es),
            "et" => Some(Self::Et),
            "ep" | "ep+" => Some(Self::Ep),
            "etp" | "etp+" => Some(Self::Etp),
            _ => None,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Es => "es",
            Self::Et => "et",
            Self::Ep => "ep",
            Self::Etp => "etp",
        }
    }

    pub fn long_name(&self) -> &'static str {
        match self {
            Self::Es => "GeneMark-ES",
            Self::Et => "GeneMark-ET",
            Self::Ep => "GeneMark-EP+",
            Self::Etp => "GeneMark-ETP+",
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Availability
// ─────────────────────────────────────────────────────────────────────────────

pub fn genemark_available() -> bool {
    for name in &["gmes_petap.pl", "gmes_linux_64", "gmhmme3"] {
        if Command::new("which")
            .arg(name)
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
        {
            return true;
        }
    }
    false
}

/// True when ProtHint's `prothint.py` is reachable on PATH.
pub fn prothint_available() -> bool {
    Command::new("which")
        .arg("prothint.py")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

fn genemark_binary() -> Option<String> {
    for name in &["gmes_petap.pl", "gmes_linux_64"] {
        if Command::new("which")
            .arg(name)
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
        {
            return Some(name.to_string());
        }
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Run GeneMark-ES
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_genemark_es(
    genome_fasta: &Path,
    out_dir: &Path,
    is_fungus: bool,
    threads: usize,
) -> Result<PathBuf> {
    let binary = genemark_binary().ok_or_else(|| {
        MycoNoteError::ExternalTool(
            "GeneMark-ES (gmes_petap.pl) not found.\n  \
             GeneMark requires a license: http://topaz.gatech.edu/GeneMark/\n  \
             After obtaining a license: conda install -c bioconda genemark-es"
                .to_string(),
        )
    })?;

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let mut cmd = Command::new(&binary);
    cmd.arg("--ES")
        .arg("--sequence")
        .arg(genome_fasta)
        .arg("--cores")
        .arg(threads.to_string())
        .current_dir(out_dir);

    if is_fungus {
        cmd.arg("--fungus");
    }

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("gmes_petap.pl: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "GeneMark-ES exited with non-zero status".to_string(),
        ));
    }

    // GeneMark-ES outputs genemark.gtf in the current dir
    let gtf_out = out_dir.join("genemark.gtf");
    if !gtf_out.exists() {
        return Err(MycoNoteError::ExternalTool(
            "GeneMark-ES output genemark.gtf not found".to_string(),
        ));
    }

    // Convert to GFF3
    let gff3_out = out_dir.join("genemark.gff3");
    convert_genemark_gtf_to_gff3(&gtf_out, &gff3_out)?;

    Ok(gff3_out)
}

/// Run GeneMark-ET with RNA-seq splice junction hints.
pub fn run_genemark_et(
    genome_fasta: &Path,
    hint_file: &Path, // intron hints in GFF format from HISAT2/STAR
    out_dir: &Path,
    is_fungus: bool,
    threads: usize,
) -> Result<PathBuf> {
    let binary = genemark_binary()
        .ok_or_else(|| MycoNoteError::ExternalTool("gmes_petap.pl not found".to_string()))?;

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let mut cmd = Command::new(&binary);
    cmd.arg("--ET")
        .arg(hint_file)
        .arg("--sequence")
        .arg(genome_fasta)
        .arg("--cores")
        .arg(threads.to_string())
        .current_dir(out_dir);

    if is_fungus {
        cmd.arg("--fungus");
    }

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("gmes_petap.pl --ET: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "GeneMark-ET failed".to_string(),
        ));
    }

    let gtf_out = out_dir.join("genemark.gtf");
    let gff3_out = out_dir.join("genemark.gff3");
    convert_genemark_gtf_to_gff3(&gtf_out, &gff3_out)?;

    Ok(gff3_out)
}

// ─────────────────────────────────────────────────────────────────────────────
// ProtHint — convert a protein FASTA into splice-site hints for EP/ETP
// ─────────────────────────────────────────────────────────────────────────────

/// Run ProtHint to generate `prothint.gff` + `evidence.gff` from a genome
/// and a protein database.  Returns the path to `prothint.gff` (the full
/// hint set; `evidence.gff` is the high-confidence subset).
///
/// `out_dir` is created if missing; ProtHint runs there and leaves all of
/// its intermediate files behind for inspection / re-use.
pub fn run_prothint(
    genome_fasta: &Path,
    protein_fasta: &Path,
    out_dir: &Path,
    is_fungus: bool,
    threads: usize,
) -> Result<PathBuf> {
    if !prothint_available() {
        return Err(MycoNoteError::ExternalTool(
            "ProtHint (prothint.py) not found on PATH.\n  \
             ProtHint ships bundled with the GeneMark-ES installer\n  \
             (<install>/ProtHint/bin) and with the bioconda `braker3` package.\n  \
             Run `myconote-cli install predict` for installation guidance."
                .to_string(),
        ));
    }

    if !protein_fasta.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Protein FASTA for ProtHint not found: {}",
            protein_fasta.display()
        )));
    }

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let mut cmd = Command::new("prothint.py");
    cmd.arg("--workdir")
        .arg(out_dir)
        .arg("--threads")
        .arg(threads.to_string());

    if is_fungus {
        cmd.arg("--fungus");
    }

    cmd.arg(genome_fasta).arg(protein_fasta);

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("prothint.py: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "ProtHint exited with non-zero status".to_string(),
        ));
    }

    let hints = out_dir.join("prothint.gff");
    if !hints.exists() {
        return Err(MycoNoteError::ExternalTool(format!(
            "ProtHint output {} not found",
            hints.display()
        )));
    }
    Ok(hints)
}

// ─────────────────────────────────────────────────────────────────────────────
// Run GeneMark-EP+ (protein-evidence-guided)
// ─────────────────────────────────────────────────────────────────────────────

/// Run GeneMark-EP+ end-to-end: ProtHint converts the genome + protein FASTA
/// into hints, then `gmes_petap.pl --EP <prothint.gff> --evidence <evidence.gff>`
/// trains and predicts.
///
/// Particularly useful for novel CTG-clade fungi without RNA-seq, since the
/// alternative-yeast nuclear code makes published predictors error-prone.
pub fn run_genemark_ep(
    genome_fasta: &Path,
    protein_fasta: &Path,
    out_dir: &Path,
    is_fungus: bool,
    threads: usize,
) -> Result<PathBuf> {
    let binary = genemark_binary()
        .ok_or_else(|| MycoNoteError::ExternalTool("gmes_petap.pl not found".to_string()))?;

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    // ── Step 1: ProtHint ──────────────────────────────────────────────────────
    let prothint_dir = out_dir.join("prothint");
    let hints = run_prothint(
        genome_fasta,
        protein_fasta,
        &prothint_dir,
        is_fungus,
        threads,
    )?;
    let evidence = prothint_dir.join("evidence.gff");

    // ── Step 2: GeneMark-EP+ ──────────────────────────────────────────────────
    let mut cmd = Command::new(&binary);
    cmd.arg("--EP")
        .arg(&hints)
        .arg("--sequence")
        .arg(genome_fasta)
        .arg("--cores")
        .arg(threads.to_string())
        .current_dir(out_dir);

    // The EP+ algorithm uses the high-confidence evidence as anchors.
    if evidence.exists() {
        cmd.arg("--evidence").arg(&evidence);
    }

    if is_fungus {
        cmd.arg("--fungus");
    }

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("gmes_petap.pl --EP: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "GeneMark-EP+ exited with non-zero status".to_string(),
        ));
    }

    let gtf_out = out_dir.join("genemark.gtf");
    if !gtf_out.exists() {
        return Err(MycoNoteError::ExternalTool(
            "GeneMark-EP+ output genemark.gtf not found".to_string(),
        ));
    }

    let gff3_out = out_dir.join("genemark.gff3");
    convert_genemark_gtf_to_gff3(&gtf_out, &gff3_out)?;
    Ok(gff3_out)
}

// ─────────────────────────────────────────────────────────────────────────────
// Run GeneMark-ETP+ (combined RNA-seq + protein evidence)
// ─────────────────────────────────────────────────────────────────────────────

/// Run GeneMark-ETP+: requires *both* an RNA-seq intron-hints GFF (--ET-style)
/// and a protein FASTA that ProtHint can convert into protein hints (--EP).
/// gmes_petap.pl is invoked with `--ETP --ET <rna_hints> --EP <prothint.gff>`.
///
/// This is the highest-quality GeneMark mode — the recommended path for fungi
/// with both RNA-seq coverage and a curated protein database (e.g. OrthoDB
/// fungi).
pub fn run_genemark_etp(
    genome_fasta: &Path,
    rna_hints: &Path,
    protein_fasta: &Path,
    out_dir: &Path,
    is_fungus: bool,
    threads: usize,
) -> Result<PathBuf> {
    let binary = genemark_binary()
        .ok_or_else(|| MycoNoteError::ExternalTool("gmes_petap.pl not found".to_string()))?;

    if !rna_hints.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "RNA-seq intron hints not found: {}",
            rna_hints.display()
        )));
    }

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    // ── Step 1: ProtHint ──────────────────────────────────────────────────────
    let prothint_dir = out_dir.join("prothint");
    let prot_hints = run_prothint(
        genome_fasta,
        protein_fasta,
        &prothint_dir,
        is_fungus,
        threads,
    )?;

    // ── Step 2: GeneMark-ETP+ ─────────────────────────────────────────────────
    let mut cmd = Command::new(&binary);
    cmd.arg("--ETP")
        .arg("--ET")
        .arg(rna_hints)
        .arg("--EP")
        .arg(&prot_hints)
        .arg("--sequence")
        .arg(genome_fasta)
        .arg("--cores")
        .arg(threads.to_string())
        .current_dir(out_dir);

    if is_fungus {
        cmd.arg("--fungus");
    }

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("gmes_petap.pl --ETP: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "GeneMark-ETP+ exited with non-zero status".to_string(),
        ));
    }

    let gtf_out = out_dir.join("genemark.gtf");
    if !gtf_out.exists() {
        return Err(MycoNoteError::ExternalTool(
            "GeneMark-ETP+ output genemark.gtf not found".to_string(),
        ));
    }

    let gff3_out = out_dir.join("genemark.gff3");
    convert_genemark_gtf_to_gff3(&gtf_out, &gff3_out)?;
    Ok(gff3_out)
}

// ─────────────────────────────────────────────────────────────────────────────
// GTF → GFF3 conversion (GeneMark-specific)
// ─────────────────────────────────────────────────────────────────────────────

/// Convert GeneMark's GTF output to GFF3 with gene/mRNA/CDS/exon hierarchy.
pub fn convert_genemark_gtf_to_gff3(gtf: &Path, gff3: &Path) -> Result<usize> {
    use std::collections::HashMap;

    let file = std::fs::File::open(gtf).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut out = std::fs::File::create(gff3).map_err(MycoNoteError::Io)?;

    writeln!(out, "##gff-version 3").map_err(MycoNoteError::Io)?;
    writeln!(out, "# Source: GeneMark-ES/ET").map_err(MycoNoteError::Io)?;

    // Collect all CDS features keyed by transcript_id
    // GTF line: seqname source feature start end score strand frame attributes
    // attributes: gene_id "xxx"; transcript_id "yyy";
    struct TxRecord {
        seqid: String,
        strand: char,
        cdss: Vec<(u64, u64, u8)>, // start, end, phase
    }

    let mut transcripts: HashMap<String, TxRecord> = HashMap::new();
    let mut tx_order: Vec<String> = vec![]; // insertion order

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let t = line.trim();
        if t.starts_with('#') || t.is_empty() {
            continue;
        }

        let fields: Vec<&str> = t.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        if fields[2] != "CDS" && fields[2] != "exon" {
            continue;
        }

        let seqid = fields[0].to_string();
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand: char = fields[6].chars().next().unwrap_or('+');
        let phase: u8 = fields[7].parse().unwrap_or(0);
        let attrs = fields[8];

        // Parse transcript_id from attributes
        let tx_id = parse_gtf_attr(attrs, "transcript_id")
            .unwrap_or_else(|| format!("{}_{}_{}", seqid, start, end));

        let entry = transcripts.entry(tx_id.clone()).or_insert_with(|| {
            tx_order.push(tx_id.clone());
            TxRecord {
                seqid: seqid.clone(),
                strand,
                cdss: vec![],
            }
        });
        entry.cdss.push((start, end, phase));
    }

    let mut gene_count = 0usize;
    for (idx, tx_id) in tx_order.iter().enumerate() {
        let rec = &transcripts[tx_id];
        if rec.cdss.is_empty() {
            continue;
        }

        let g_start = rec.cdss.iter().map(|c| c.0).min().unwrap_or(0);
        let g_end = rec.cdss.iter().map(|c| c.1).max().unwrap_or(0);

        let gene_id = format!("GM_{:06}", idx + 1);
        let mrna_id = format!("{}.mRNA1", gene_id);

        writeln!(
            out,
            "{}\tGeneMark\tgene\t{}\t{}\t.\t{}\t.\tID={}",
            rec.seqid, g_start, g_end, rec.strand, gene_id
        )
        .map_err(MycoNoteError::Io)?;
        writeln!(
            out,
            "{}\tGeneMark\tmRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={}",
            rec.seqid, g_start, g_end, rec.strand, mrna_id, gene_id
        )
        .map_err(MycoNoteError::Io)?;

        let mut sorted_cds = rec.cdss.clone();
        sorted_cds.sort_by_key(|c| c.0);

        for (i, (start, end, phase)) in sorted_cds.iter().enumerate() {
            let exon_id = format!("{}.exon{}", mrna_id, i + 1);
            let cds_id = format!("{}.CDS{}", mrna_id, i + 1);
            writeln!(
                out,
                "{}\tGeneMark\texon\t{}\t{}\t.\t{}\t.\tID={};Parent={}",
                rec.seqid, start, end, rec.strand, exon_id, mrna_id
            )
            .map_err(MycoNoteError::Io)?;
            writeln!(
                out,
                "{}\tGeneMark\tCDS\t{}\t{}\t.\t{}\t{}\tID={};Parent={}",
                rec.seqid, start, end, rec.strand, phase, cds_id, mrna_id
            )
            .map_err(MycoNoteError::Io)?;
        }

        gene_count += 1;
    }

    Ok(gene_count)
}

fn parse_gtf_attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        if part.starts_with(key) {
            let value = part[key.len()..].trim().trim_matches('"').to_string();
            if !value.is_empty() {
                return Some(value);
            }
        }
    }
    None
}
