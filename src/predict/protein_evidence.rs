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

/// Is `diamond` available for the protein-DB prefilter?
pub fn diamond_available() -> bool {
    Command::new("diamond")
        .arg("version")
        .output()
        .map(|o| o.status.success() || !o.stdout.is_empty())
        .unwrap_or(false)
}

// ─────────────────────────────────────────────────────────────────────────────
// Protein-DB prefilter (FIX B, v0.7.10) — BRAKER-style reduction before miniprot
// ─────────────────────────────────────────────────────────────────────────────
//
// Feeding an entire OrthoDB partition (e.g. Fungi.fa, ~3.9 GB) straight to
// miniprot never finishes in a per-genome budget (an A/B run was still aligning
// after 12.5 h). BRAKER/ProtHint avoid this by prefiltering the protein set to
// the genome's plausible homologs before any splice-aware alignment. We do the
// same: diamond blastx the genome against the protein DB, keep only proteins
// with a hit, and cap the survivors. With no diamond, we cap the DB to a bounded
// subset so the run still completes rather than aligning the whole thing.

/// Tunables for `prefilter_protein_db`.
#[derive(Debug, Clone)]
pub struct PrefilterParams {
    /// diamond `--evalue` threshold for keeping a protein hit.
    pub evalue: f64,
    /// Maximum proteins fed to miniprot after prefiltering. `0` = unbounded.
    pub max_proteins: usize,
    /// Threads for diamond.
    pub threads: usize,
}

impl Default for PrefilterParams {
    fn default() -> Self {
        Self {
            evalue: 1e-5,
            max_proteins: 50_000,
            threads: 4,
        }
    }
}

/// Outcome of a prefilter pass.
#[derive(Debug)]
pub struct PrefilterResult {
    /// Reduced protein FASTA to feed to miniprot.
    pub proteins: PathBuf,
    /// Number of proteins kept.
    pub n_kept: usize,
    /// How the set was reduced: `"diamond"`, `"cap"`, or `"passthrough"`.
    pub method: &'static str,
    /// Whether the `max_proteins` cap actually truncated the set.
    pub capped: bool,
}

/// Write a single FASTA record (60-char wrapped) to `out`.
fn write_fasta_record(out: &mut impl Write, rec: &crate::parser::fasta::FastaRecord) -> Result<()> {
    writeln!(out, ">{}", rec.header).map_err(MycoNoteError::Io)?;
    for chunk in rec.sequence.as_bytes().chunks(60) {
        out.write_all(chunk).map_err(MycoNoteError::Io)?;
        writeln!(out).map_err(MycoNoteError::Io)?;
    }
    Ok(())
}

/// Parse diamond `--outfmt 6` output and collect the set of subject (protein)
/// IDs with a hit. Column 2 (0-based index 1) is `sseqid`.
fn parse_diamond_hit_ids(tsv: &Path) -> Result<std::collections::HashSet<String>> {
    let file = std::fs::File::open(tsv).map_err(MycoNoteError::Io)?;
    let mut ids = std::collections::HashSet::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        if let Some(sseqid) = t.split('\t').nth(1) {
            if !sseqid.is_empty() {
                ids.insert(sseqid.to_string());
            }
        }
    }
    Ok(ids)
}

/// Stream `in_fa`, writing only records whose bare ID is in `ids`, up to `max`
/// records (`0` = unbounded). Streams one record at a time so a multi-GB DB is
/// never loaded into memory. Returns `(n_kept, capped)`.
fn subset_fasta_by_ids(
    in_fa: &Path,
    out_fa: &Path,
    ids: &std::collections::HashSet<String>,
    max: usize,
) -> Result<(usize, bool)> {
    let reader = crate::parser::fasta::FastaReader::from_path(in_fa)?;
    let mut out = std::fs::File::create(out_fa).map_err(MycoNoteError::Io)?;
    let mut kept = 0usize;
    let mut capped = false;
    for rec in reader {
        let rec = rec?;
        if !ids.contains(&rec.id) {
            continue;
        }
        if max != 0 && kept >= max {
            capped = true;
            break;
        }
        write_fasta_record(&mut out, &rec)?;
        kept += 1;
    }
    Ok((kept, capped))
}

/// Stream `in_fa`, writing the first `max` records (`0` = unbounded, copies all).
/// Returns `(n_written, capped)`.
fn cap_fasta(in_fa: &Path, out_fa: &Path, max: usize) -> Result<(usize, bool)> {
    let reader = crate::parser::fasta::FastaReader::from_path(in_fa)?;
    let mut out = std::fs::File::create(out_fa).map_err(MycoNoteError::Io)?;
    let mut n = 0usize;
    let mut capped = false;
    for rec in reader {
        let rec = rec?;
        if max != 0 && n >= max {
            capped = true;
            break;
        }
        write_fasta_record(&mut out, &rec)?;
        n += 1;
    }
    Ok((n, capped))
}

/// Run `diamond makedb` + `diamond blastx` (genome query vs protein DB) and
/// return the set of protein IDs with a hit at or below `params.evalue`. A tool
/// failure bubbles up so the caller can fall back to the bounded-subset path.
fn run_diamond_prefilter(
    proteins: &Path,
    genome: &Path,
    out_dir: &Path,
    params: &PrefilterParams,
) -> Result<std::collections::HashSet<String>> {
    let db = out_dir.join("prefilter_db");
    let hits = out_dir.join("prefilter_hits.tsv");

    let makedb = Command::new("diamond")
        .arg("makedb")
        .arg("--in")
        .arg(proteins)
        .arg("--db")
        .arg(&db)
        .args(["--threads", &params.threads.to_string()])
        .arg("--quiet")
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("diamond makedb: {}", e)))?;
    if !makedb.success() {
        return Err(MycoNoteError::ExternalTool(
            "diamond makedb failed".to_string(),
        ));
    }

    let search = Command::new("diamond")
        .arg("blastx")
        .arg("--db")
        .arg(&db)
        .arg("--query")
        .arg(genome)
        .arg("--out")
        .arg(&hits)
        .args(["--outfmt", "6", "qseqid", "sseqid", "pident", "evalue"])
        .args(["--evalue", &format!("{:e}", params.evalue)])
        .args(["--max-target-seqs", "25"])
        .args(["--threads", &params.threads.to_string()])
        .arg("--quiet")
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("diamond blastx: {}", e)))?;
    if !search.success() {
        return Err(MycoNoteError::ExternalTool(
            "diamond blastx failed".to_string(),
        ));
    }

    parse_diamond_hit_ids(&hits)
}

/// Reduce a (potentially huge) protein DB to a miniprot-sized set before
/// alignment, BRAKER-style:
///   1. If diamond is available, blastx the genome against the DB and keep only
///      proteins with a hit (≤ `evalue`), capped at `max_proteins`.
///   2. Otherwise — or if diamond fails or finds nothing — fall back to the
///      first `max_proteins` records of the DB.
///
/// Always returns a usable reduced FASTA; never fails the run merely because
/// diamond is absent — the whole point is to avoid handing 3.9 GB to miniprot.
pub fn prefilter_protein_db(
    proteins: &Path,
    genome: &Path,
    out_dir: &Path,
    params: &PrefilterParams,
) -> Result<PrefilterResult> {
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;
    let reduced = out_dir.join("proteins_prefiltered.faa");

    if diamond_available() {
        println!(
            "  Prefiltering protein DB with diamond blastx (evalue ≤ {:e})…",
            params.evalue
        );
        match run_diamond_prefilter(proteins, genome, out_dir, params) {
            Ok(ids) if !ids.is_empty() => {
                let (kept, capped) =
                    subset_fasta_by_ids(proteins, &reduced, &ids, params.max_proteins)?;
                if kept > 0 {
                    println!(
                        "      diamond kept {} proteins{} → {}",
                        kept,
                        if capped {
                            format!(" (capped at {})", params.max_proteins)
                        } else {
                            String::new()
                        },
                        reduced.display()
                    );
                    return Ok(PrefilterResult {
                        proteins: reduced,
                        n_kept: kept,
                        method: "diamond",
                        capped,
                    });
                }
                eprintln!(
                    "  ⚠  diamond hits matched no DB record IDs — falling back to a bounded subset."
                );
            }
            Ok(_) => {
                eprintln!(
                    "  ⚠  diamond found no protein hits — falling back to a bounded subset of {} proteins.",
                    params.max_proteins
                );
            }
            Err(e) => {
                eprintln!(
                    "  ⚠  diamond prefilter failed ({}) — falling back to a bounded subset.",
                    e
                );
            }
        }
    } else {
        eprintln!(
            "  ⚠  diamond not found — capping protein DB to {} records before miniprot \
             (install diamond for a genome-specific prefilter: conda install -c bioconda diamond).",
            params.max_proteins
        );
    }

    // Fallback: bounded subset (first N records).
    let (kept, capped) = cap_fasta(proteins, &reduced, params.max_proteins)?;
    let method = if capped { "cap" } else { "passthrough" };
    println!(
        "      using {} proteins{} → {}",
        kept,
        if capped { " (capped)" } else { "" },
        reduced.display()
    );
    Ok(PrefilterResult {
        proteins: reduced,
        n_kept: kept,
        method,
        capped,
    })
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

// ─────────────────────────────────────────────────────────────────────────────
// Tests — protein-DB prefilter (FIX B)
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn write(path: &Path, body: &str) {
        std::fs::write(path, body).unwrap();
    }

    const PROTEINS: &str = "\
>p1 first
MKVLAA
>p2 second
MARNDC
>p3 third
MQEGHI
>p4 fourth
MKLMNP
>p5 fifth
MQRSTV
";

    // A mock diamond outfmt-6 table → only the hit proteins are passed on.
    #[test]
    fn parse_diamond_hits_collects_subject_ids() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("hits.tsv");
        write(
            &tsv,
            "contig1\tp2\t88.0\t1e-20\ncontig1\tp4\t72.0\t3e-10\ncontig2\tp2\t91.0\t1e-30\n",
        );
        let ids = parse_diamond_hit_ids(&tsv).unwrap();
        let expected: HashSet<String> = ["p2".to_string(), "p4".to_string()].into_iter().collect();
        assert_eq!(ids, expected);
    }

    // Prefilter reduces the protein set: a mock diamond hit set → only the hit
    // proteins survive into the reduced FASTA.
    #[test]
    fn subset_keeps_only_hit_proteins() {
        let dir = tempfile::tempdir().unwrap();
        let in_fa = dir.path().join("db.faa");
        let out_fa = dir.path().join("reduced.faa");
        write(&in_fa, PROTEINS);

        let ids: HashSet<String> = ["p1".to_string(), "p3".to_string()].into_iter().collect();
        let (kept, capped) = subset_fasta_by_ids(&in_fa, &out_fa, &ids, 0).unwrap();
        assert_eq!(kept, 2);
        assert!(!capped);

        let body = std::fs::read_to_string(&out_fa).unwrap();
        assert!(body.contains(">p1"));
        assert!(body.contains(">p3"));
        assert!(!body.contains(">p2"), "non-hit protein must be dropped");
        assert!(!body.contains(">p5"));
    }

    // The max cap bounds the subset even when more proteins would hit.
    #[test]
    fn subset_respects_max_cap() {
        let dir = tempfile::tempdir().unwrap();
        let in_fa = dir.path().join("db.faa");
        let out_fa = dir.path().join("reduced.faa");
        write(&in_fa, PROTEINS);

        let ids: HashSet<String> = ["p1", "p2", "p3", "p4", "p5"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let (kept, capped) = subset_fasta_by_ids(&in_fa, &out_fa, &ids, 2).unwrap();
        assert_eq!(kept, 2);
        assert!(capped);
    }

    // The cap-only fallback path (no diamond) bounds the input set.
    #[test]
    fn cap_fasta_bounds_input() {
        let dir = tempfile::tempdir().unwrap();
        let in_fa = dir.path().join("db.faa");
        write(&in_fa, PROTEINS);

        let capped_out = dir.path().join("capped.faa");
        let (n, capped) = cap_fasta(&in_fa, &capped_out, 3).unwrap();
        assert_eq!(n, 3);
        assert!(capped);
        let body = std::fs::read_to_string(&capped_out).unwrap();
        assert_eq!(body.matches('>').count(), 3);

        // max = 0 copies everything, no truncation.
        let all_out = dir.path().join("all.faa");
        let (n_all, capped_all) = cap_fasta(&in_fa, &all_out, 0).unwrap();
        assert_eq!(n_all, 5);
        assert!(!capped_all);
    }

    // End-to-end prefilter: with no diamond hits (or no diamond), the orchestrator
    // always returns a usable, bounded reduced FASTA — never the whole DB, never
    // an error for a missing tool. Robust whether or not diamond is installed:
    // both the diamond-empty and no-diamond paths converge on the cap fallback.
    #[test]
    fn prefilter_db_bounds_input_and_is_nonfatal() {
        let dir = tempfile::tempdir().unwrap();
        let proteins = dir.path().join("db.faa");
        write(&proteins, PROTEINS);
        // A tiny genome with no homology to the mock proteins → zero diamond hits.
        let genome = dir.path().join("genome.fa");
        write(&genome, ">contig1\nACGTACGTACGTACGTACGTACGT\n");
        let out_dir = dir.path().join("prefilter");

        let params = PrefilterParams {
            evalue: 1e-5,
            max_proteins: 2,
            threads: 1,
        };
        let res = prefilter_protein_db(&proteins, &genome, &out_dir, &params).unwrap();

        assert!(res.proteins.exists(), "reduced FASTA must be written");
        assert!(res.n_kept >= 1, "at least one protein survives");
        assert!(
            res.n_kept <= 2,
            "the max cap bounds the set regardless of path (got {})",
            res.n_kept
        );
    }
}
