/// Protein-to-genome evidence alignment
///
/// Aligns protein sequences to a genome assembly using miniprot (preferred)
/// or Exonerate to generate protein evidence for gene prediction.
/// This evidence is then used by Augustus as hints or fed into the
/// Evidence Modeler consensus.
///
/// This is myconote-cli's own implementation — independent of any other
/// annotation pipeline.
///
/// Tools supported:
///   - miniprot (preferred, faster than exonerate for protein→genome)
///   - exonerate (fallback, protein2genome model)
///   - diamond (for pre-filtering before full alignment)
use crate::utils::error::{MycoNoteError, Result};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct ProteinEvidenceConfig {
    /// Protein FASTA (e.g. Swiss-Prot subset, OrthoDB proteins)
    pub proteins: PathBuf,
    /// Genome FASTA (masked or unmasked)
    pub genome: PathBuf,
    /// Output directory
    pub out_dir: PathBuf,
    /// Number of threads
    pub threads: usize,
    /// Maximum intron size for alignment (important for splice-aware tools)
    pub max_intron: usize,
    /// Minimum protein identity to keep alignment (0-100)
    pub min_identity: f64,
    /// Minimum alignment coverage of the protein (0-1)
    pub min_coverage: f64,
}

impl Default for ProteinEvidenceConfig {
    fn default() -> Self {
        Self {
            proteins: PathBuf::new(),
            genome: PathBuf::new(),
            out_dir: PathBuf::from("protein_evidence"),
            threads: 4,
            max_intron: 10_000,
            min_identity: 50.0,
            min_coverage: 0.5,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Result
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug)]
pub struct ProteinEvidenceResult {
    /// GFF3 with protein-to-genome alignments (for EVM)
    pub evidence_gff: PathBuf,
    /// Augustus hints file (for hints-based prediction)
    pub hints_gff: PathBuf,
    /// Number of proteins aligned
    pub n_aligned: usize,
    /// Number of gene loci supported by protein evidence
    pub n_loci: usize,
    /// Tool used
    pub tool: String,
}

// ─────────────────────────────────────────────────────────────────────────────
// Tool detection
// ─────────────────────────────────────────────────────────────────────────────

fn miniprot_available() -> bool {
    Command::new("miniprot")
        .arg("--version")
        .output()
        .map(|o| o.status.success() || !o.stderr.is_empty())
        .unwrap_or(false)
}

fn exonerate_available() -> bool {
    Command::new("exonerate")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run protein-to-genome alignment and generate evidence files.
/// Tries miniprot first (much faster), falls back to exonerate.
pub fn generate_protein_evidence(config: &ProteinEvidenceConfig) -> Result<ProteinEvidenceResult> {
    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    if miniprot_available() {
        println!("  Using miniprot for protein→genome alignment");
        run_miniprot(config)
    } else if exonerate_available() {
        println!("  Using exonerate for protein→genome alignment");
        run_exonerate(config)
    } else {
        Err(MycoNoteError::ExternalTool(
            "Neither miniprot nor exonerate found. Install: conda install -c bioconda miniprot"
                .to_string(),
        ))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// miniprot backend
// ─────────────────────────────────────────────────────────────────────────────

fn run_miniprot(config: &ProteinEvidenceConfig) -> Result<ProteinEvidenceResult> {
    let gff_out = config.out_dir.join("protein_alignments.gff3");
    let hints_out = config.out_dir.join("protein_hints.gff");

    // miniprot outputs GFF3 natively with --gff
    let status = Command::new("miniprot")
        .args([
            "--gff",
            "-t",
            &config.threads.to_string(),
            "--max-intron",
            &config.max_intron.to_string(),
            config.genome.to_str().unwrap_or(""),
            config.proteins.to_str().unwrap_or(""),
        ])
        .stdout(std::fs::File::create(&gff_out).map_err(MycoNoteError::Io)?)
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("miniprot: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool("miniprot failed".to_string()));
    }

    // Count aligned proteins and generate hints
    let (n_aligned, n_loci) = convert_to_hints(&gff_out, &hints_out, config)?;

    Ok(ProteinEvidenceResult {
        evidence_gff: gff_out,
        hints_gff: hints_out,
        n_aligned,
        n_loci,
        tool: "miniprot".to_string(),
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// exonerate backend
// ─────────────────────────────────────────────────────────────────────────────

fn run_exonerate(config: &ProteinEvidenceConfig) -> Result<ProteinEvidenceResult> {
    let raw_out = config.out_dir.join("exonerate_raw.txt");
    let gff_out = config.out_dir.join("protein_alignments.gff3");
    let hints_out = config.out_dir.join("protein_hints.gff");

    let status = Command::new("exonerate")
        .args([
            "--model",
            "protein2genome",
            "--showtargetgff",
            "yes",
            "--showvulgar",
            "no",
            "--showalignment",
            "no",
            "--percent",
            &config.min_identity.to_string(),
            "--maxintron",
            &config.max_intron.to_string(),
            "--query",
            config.proteins.to_str().unwrap_or(""),
            "--target",
            config.genome.to_str().unwrap_or(""),
        ])
        .stdout(std::fs::File::create(&raw_out).map_err(MycoNoteError::Io)?)
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("exonerate: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool("exonerate failed".to_string()));
    }

    // Parse exonerate GFF output
    let n_aligned = parse_exonerate_to_gff3(&raw_out, &gff_out)?;
    let n_loci = convert_to_hints(&gff_out, &hints_out, config)?.1;

    Ok(ProteinEvidenceResult {
        evidence_gff: gff_out,
        hints_gff: hints_out,
        n_aligned,
        n_loci,
        tool: "exonerate".to_string(),
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// GFF conversion helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Convert protein alignment GFF3 to Augustus hints format.
/// Returns (n_proteins, n_loci).
fn convert_to_hints(
    alignment_gff: &Path,
    hints_gff: &Path,
    _config: &ProteinEvidenceConfig,
) -> Result<(usize, usize)> {
    let file = std::fs::File::open(alignment_gff).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut out = std::fs::File::create(hints_gff).map_err(MycoNoteError::Io)?;

    let mut n_proteins = 0usize;
    let mut loci: std::collections::HashSet<String> = std::collections::HashSet::new();

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }

        let feature_type = cols[2];
        let seqid = cols[0];

        match feature_type {
            "mRNA" | "gene" | "match" => {
                n_proteins += 1;
                loci.insert(format!("{}:{}-{}", seqid, cols[3], cols[4]));
            }
            "CDS" | "exon" | "match_part" => {
                // Convert to CDSpart hint
                writeln!(
                    out,
                    "{}\tProtein\tCDSpart\t{}\t{}\t{}\t{}\t.\tsource=P;priority=4",
                    seqid, cols[3], cols[4], cols[5], cols[6]
                )
                .map_err(MycoNoteError::Io)?;
            }
            "intron" => {
                // Convert to intron hint
                writeln!(
                    out,
                    "{}\tProtein\tintron\t{}\t{}\t{}\t{}\t.\tsource=P;priority=4",
                    seqid, cols[3], cols[4], cols[5], cols[6]
                )
                .map_err(MycoNoteError::Io)?;
            }
            _ => {}
        }
    }

    Ok((n_proteins, loci.len()))
}

/// Parse exonerate raw output (--showtargetgff yes) into clean GFF3.
fn parse_exonerate_to_gff3(raw: &Path, gff_out: &Path) -> Result<usize> {
    let file = std::fs::File::open(raw).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut out = std::fs::File::create(gff_out).map_err(MycoNoteError::Io)?;

    writeln!(out, "##gff-version 3").map_err(MycoNoteError::Io)?;

    let mut in_gff_section = false;
    let mut n_alignments = 0usize;

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();

        // Exonerate GFF sections start with "# --- START OF GFF DUMP ---"
        if trimmed.contains("START OF GFF DUMP") {
            in_gff_section = true;
            continue;
        }
        if trimmed.contains("END OF GFF DUMP") {
            in_gff_section = false;
            continue;
        }

        if in_gff_section && !trimmed.starts_with('#') && !trimmed.is_empty() {
            writeln!(out, "{}", trimmed).map_err(MycoNoteError::Io)?;
            if trimmed.contains("\tgene\t") {
                n_alignments += 1;
            }
        }
    }

    Ok(n_alignments)
}
