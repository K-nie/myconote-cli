/// Gene model update / UTR refinement
///
/// The `update` command refines an initial gene annotation by:
///   1. Aligning assembled transcripts (Trinity / existing RNA-seq) to the genome
///      with minimap2 (splice-aware, long-read mode for Trinity assemblies)
///   2. Loading transcript alignments into a PASA database (if PASA is available)
///   3. Updating gene model boundaries (UTR extension, intron corrections, alt-isoforms)
///   4. Writing a new, refined GFF3 and summary of changes
///
/// If PASA is not available, a lightweight fallback updates UTRs using only
/// coverage information from the minimap2 BAM/PAF output.
///
/// Equivalent to: funannotate update --input <gff3> --fasta <genome> --rna-bam <bam>
pub mod pasa_update;

use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct UpdateConfig {
    /// Input GFF3 from `myconote predict`
    pub gff: PathBuf,
    /// Genome FASTA (must match GFF3 seqids)
    pub fasta: PathBuf,
    /// Output directory
    pub out_dir: PathBuf,
    /// RNA-seq reads FASTQ (paired or single) — assembled with Trinity internally
    pub rna_r1: Option<PathBuf>,
    pub rna_r2: Option<PathBuf>,
    /// Pre-assembled transcript FASTA (Trinity output or any cDNA FASTA)
    pub transcripts: Option<PathBuf>,
    /// Pre-aligned BAM (minimap2 / HISAT2 / STAR vs genome)
    pub rna_bam: Option<PathBuf>,
    /// PASA database name (sqlite db created in out_dir if not given)
    pub pasa_db: Option<String>,
    /// Organism name for PASA config
    pub organism: Option<String>,
    /// Locus tag prefix (must match predict step)
    pub locus_prefix: String,
    /// Threads
    pub threads: usize,
    /// Maximum UTR extension in bp (prevents runaway UTRs)
    pub max_utr_extension: u64,
    /// Minimum transcript alignment identity (0–1)
    pub min_identity: f64,
    /// Minimum transcript alignment coverage (0–1)
    pub min_coverage: f64,
}

impl Default for UpdateConfig {
    fn default() -> Self {
        Self {
            gff: PathBuf::new(),
            fasta: PathBuf::new(),
            out_dir: PathBuf::from("update_out"),
            rna_r1: None,
            rna_r2: None,
            transcripts: None,
            rna_bam: None,
            pasa_db: None,
            organism: None,
            locus_prefix: "GENE".to_string(),
            threads: 4,
            max_utr_extension: 2000,
            min_identity: 0.95,
            min_coverage: 0.90,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Update result
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
pub struct UpdateResult {
    pub output_gff: PathBuf,
    pub n_genes_total: usize,
    pub n_genes_updated: usize,
    pub n_utr5_added: usize,
    pub n_utr3_added: usize,
    pub n_isoforms_added: usize,
    pub pasa_used: bool,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_update(config: &UpdateConfig) -> Result<UpdateResult> {
    if !config.gff.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "GFF3 not found: {}",
            config.gff.display()
        )));
    }
    if !config.fasta.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Genome FASTA not found: {}",
            config.fasta.display()
        )));
    }

    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    println!("── Gene model update ────────────────────────────────────────");
    println!("  GFF3     : {}", config.gff.display());
    println!("  Genome   : {}", config.fasta.display());
    println!("  Output   : {}", config.out_dir.display());

    // ── Step 1: assemble / gather transcripts ────────────────────────────────
    let transcript_fa = resolve_transcripts(config)?;

    // ── Step 2: align transcripts → genome ───────────────────────────────────
    let bam_path = if let Some(ref bam) = config.rna_bam {
        println!("  Using pre-aligned BAM: {}", bam.display());
        bam.clone()
    } else if let Some(ref tx) = transcript_fa {
        println!("  Aligning transcripts to genome (minimap2)…");
        align_transcripts_to_genome(tx, &config.fasta, &config.out_dir, config.threads)?
    } else {
        return Err(MycoNoteError::InvalidFormat(
            "No RNA-seq input provided. Supply --rna-r1/r2, --transcripts, or --rna-bam."
                .to_string(),
        ));
    };

    // ── Step 3: PASA update or lightweight fallback ──────────────────────────
    let use_pasa = pasa_update::pasa_available();

    let result = if use_pasa {
        println!("  Running PASA gene model update…");
        pasa_update::run_pasa_update(config, &bam_path, transcript_fa.as_deref())?
    } else {
        println!("  PASA not found — using lightweight UTR extension fallback.");
        println!("  (Install PASA for full isoform-aware gene model updates)");
        run_lightweight_update(config, &bam_path)?
    };

    // ── Summary ──────────────────────────────────────────────────────────────
    println!(
        "  ✓  Update complete: {} / {} genes updated",
        result.n_genes_updated, result.n_genes_total
    );
    if result.n_utr5_added > 0 || result.n_utr3_added > 0 {
        println!(
            "     5' UTRs added: {}   3' UTRs added: {}",
            result.n_utr5_added, result.n_utr3_added
        );
    }
    if result.n_isoforms_added > 0 {
        println!(
            "     Alternative isoforms added: {}",
            result.n_isoforms_added
        );
    }
    println!("  ✓  Updated GFF3 → {}", result.output_gff.display());

    Ok(result)
}

// ─────────────────────────────────────────────────────────────────────────────
// Transcript resolution / Trinity assembly
// ─────────────────────────────────────────────────────────────────────────────

fn resolve_transcripts(config: &UpdateConfig) -> Result<Option<PathBuf>> {
    if let Some(ref tx) = config.transcripts {
        if tx.exists() {
            return Ok(Some(tx.clone()));
        }
        return Err(MycoNoteError::InvalidFormat(format!(
            "Transcript FASTA not found: {}",
            tx.display()
        )));
    }

    if let Some(ref r1) = config.rna_r1 {
        use crate::train::{trinity, TrainConfig};
        let trinity_dir = config.out_dir.join("trinity_update");
        println!("  Assembling RNA-seq with Trinity…");
        let mut train_cfg = TrainConfig::default();
        train_cfg.left_reads = vec![r1.clone()];
        if let Some(ref r2_path) = config.rna_r2 {
            train_cfg.right_reads = vec![r2_path.clone()];
        }
        train_cfg.threads = config.threads;
        train_cfg.out_dir = trinity_dir.clone();
        let trinity_fa = trinity::run_trinity(&train_cfg, &trinity_dir)?;
        return Ok(Some(trinity_fa));
    }

    Ok(None)
}

// ─────────────────────────────────────────────────────────────────────────────
// Transcript → genome alignment (minimap2, splice mode)
// ─────────────────────────────────────────────────────────────────────────────

fn align_transcripts_to_genome(
    transcripts: &Path,
    genome: &Path,
    out_dir: &Path,
    threads: usize,
) -> Result<PathBuf> {
    let bam_path = out_dir.join("transcripts_vs_genome.bam");

    // Check minimap2
    let mm2_ok = Command::new("minimap2")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !mm2_ok {
        return Err(MycoNoteError::ExternalTool(
            "minimap2 not found. Install: conda install -c bioconda minimap2".to_string(),
        ));
    }

    // minimap2 -ax splice:hq for assembled transcripts (Trinity cDNA)
    // Pipe through samtools sort → BAM
    let sam_path = out_dir.join("transcripts_vs_genome.sam");

    let genome_s = genome.to_string_lossy();
    let tx_s = transcripts.to_string_lossy();
    let sam_s = sam_path.to_string_lossy();

    let mm2_status = Command::new("minimap2")
        .args([
            "-ax",
            "splice:hq",
            "--secondary=no",
            "-t",
            &threads.to_string(),
            genome_s.as_ref(),
            tx_s.as_ref(),
            "-o",
            sam_s.as_ref(),
        ])
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("minimap2: {}", e)))?;

    if !mm2_status.success() {
        return Err(MycoNoteError::ExternalTool(
            "minimap2 alignment failed".to_string(),
        ));
    }

    // samtools sort
    let bam_s = bam_path.to_string_lossy();
    let sam_ok = Command::new("samtools")
        .args([
            "sort",
            "-o",
            bam_s.as_ref(),
            "-@",
            &threads.to_string(),
            sam_s.as_ref(),
        ])
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("samtools sort: {}", e)));

    match sam_ok {
        Ok(s) if s.success() => {
            let _ = std::fs::remove_file(&sam_path);
            // Index BAM
            let _ = Command::new("samtools")
                .args(["index", bam_s.as_ref()])
                .status();
        }
        _ => {
            // samtools unavailable — keep SAM, return as "BAM" path
            return Ok(sam_path);
        }
    }

    Ok(bam_path)
}

// ─────────────────────────────────────────────────────────────────────────────
// Lightweight UTR extension fallback (no PASA required)
// ─────────────────────────────────────────────────────────────────────────────

/// Parse coverage from a BED/PAF/BAM to extend UTRs, without a full PASA run.
/// This is a simplified approach: parse the BAM with `samtools view` into
/// transcript blocks, then extend gene 5'/3' ends to match transcript extents.
fn run_lightweight_update(config: &UpdateConfig, bam_path: &Path) -> Result<UpdateResult> {
    use crate::parser::gff::GFFReader;

    // ── Read existing GFF3 ────────────────────────────────────────────────────
    let records: Vec<_> = GFFReader::from_path(&config.gff)?
        .filter_map(|r| r.ok())
        .collect();

    // Build per-gene structures: gene_id → (seqid, strand, start, end)
    let mut genes: HashMap<String, (String, char, u64, u64)> = HashMap::new();
    for rec in &records {
        if rec.feature_type == "gene" {
            if let Some(id) = rec.id() {
                genes.insert(
                    id.clone(),
                    (rec.seqid.clone(), rec.strand, rec.start, rec.end),
                );
            }
        }
    }

    // ── Build transcript coverage blocks from BAM (samtools view) ─────────────
    // transcript_id → (seqid, start, end)
    let mut tx_extents: HashMap<String, (String, u64, u64)> = HashMap::new();
    let bam_str = bam_path.to_string_lossy();
    if let Ok(out) = Command::new("samtools")
        .args(["view", "-F", "4", &bam_str])
        .output()
    {
        for line in out.stdout.split(|&b| b == b'\n') {
            let line = std::str::from_utf8(line).unwrap_or("").trim();
            if line.is_empty() || line.starts_with('@') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 10 {
                continue;
            }
            let tx_id = cols[0].to_string();
            let seqid = cols[2].to_string();
            let pos: u64 = cols[3].parse().unwrap_or(0);
            let cigar = cols[5];
            let end = cigar_end(pos, cigar);
            let entry = tx_extents
                .entry(tx_id)
                .or_insert_with(|| (seqid.clone(), pos, end));
            if pos < entry.1 {
                entry.1 = pos;
            }
            if end > entry.2 {
                entry.2 = end;
            }
        }
    }

    // ── Write updated GFF3 ───────────────────────────────────────────────────
    let out_gff = config.out_dir.join("updated.gff3");
    let mut f = std::fs::File::create(&out_gff).map_err(MycoNoteError::Io)?;
    writeln!(f, "##gff-version 3").map_err(MycoNoteError::Io)?;

    let mut n_updated = 0usize;
    let mut n_utr5 = 0usize;
    let mut n_utr3 = 0usize;
    let n_total = genes.len();

    // Track updated gene boundaries
    let mut updated_genes: HashMap<String, (u64, u64)> = HashMap::new();

    for (gene_id, (seqid, strand, g_start, g_end)) in &genes {
        let mut new_start = *g_start;
        let mut new_end = *g_end;
        let max_ext = config.max_utr_extension;

        // Find transcript alignments overlapping this gene
        for (tx_id, (tx_seqid, tx_start, tx_end)) in &tx_extents {
            if tx_seqid != seqid {
                continue;
            }
            if *tx_end < *g_start || *tx_start > *g_end {
                continue;
            } // no overlap
            let _ = tx_id; // suppress unused warning

            // Extend UTR (capped at max_utr_extension)
            if *tx_start < new_start && new_start.saturating_sub(*tx_start) <= max_ext {
                new_start = *tx_start;
            }
            if *tx_end > new_end && tx_end.saturating_sub(new_end) <= max_ext {
                new_end = *tx_end;
            }
        }

        let changed = new_start != *g_start || new_end != *g_end;
        if changed {
            n_updated += 1;
            if *strand == '+' {
                if new_start < *g_start {
                    n_utr5 += 1;
                }
                if new_end > *g_end {
                    n_utr3 += 1;
                }
            } else {
                if new_end > *g_end {
                    n_utr5 += 1;
                }
                if new_start < *g_start {
                    n_utr3 += 1;
                }
            }
            updated_genes.insert(gene_id.clone(), (new_start, new_end));
        }
    }

    // Write all records, updating gene/mRNA boundaries where changed
    for mut rec in records {
        if rec.feature_type == "gene" {
            if let Some(id) = rec.id().map(|s| s.to_string()) {
                if let Some(&(ns, ne)) = updated_genes.get(&id) {
                    rec.start = ns;
                    rec.end = ne;
                }
            }
        } else if rec.feature_type == "mRNA" {
            // Find parent gene
            if let Some(parent) = rec.parent().map(|s| s.to_string()) {
                if let Some(&(ns, ne)) = updated_genes.get(&parent) {
                    if ns < rec.start {
                        rec.start = ns;
                    }
                    if ne > rec.end {
                        rec.end = ne;
                    }
                }
            }
        }
        writeln!(f, "{}", rec.to_gff3_line()).map_err(MycoNoteError::Io)?;
    }

    Ok(UpdateResult {
        output_gff: out_gff,
        n_genes_total: n_total,
        n_genes_updated: n_updated,
        n_utr5_added: n_utr5,
        n_utr3_added: n_utr3,
        n_isoforms_added: 0,
        pasa_used: false,
    })
}

/// Compute alignment end position from CIGAR string.
fn cigar_end(start: u64, cigar: &str) -> u64 {
    let mut pos = start;
    let mut num = 0u64;
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch as u64 - '0' as u64);
        } else {
            match ch {
                'M' | 'D' | 'N' | '=' | 'X' => pos += num,
                _ => {}
            }
            num = 0;
        }
    }
    pos
}
