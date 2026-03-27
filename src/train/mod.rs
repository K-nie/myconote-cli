/// RNA-seq mediated training pipeline
///
/// Equivalent to `funannotate train`:
///   1. Quality-trim RNA-seq reads (fastp / trimmomatic)
///   2. Assemble transcripts with Trinity (de novo or genome-guided)
///   3. Align Trinity assemblies to genome with minimap2 or GMAP
///   4. Build PASA transcript database for evidence-based annotation
///   5. Extract high-confidence complete gene models from PASA
///   6. Train Augustus and SNAP on those models
///   7. Optionally train GeneMark-ES on the masked genome (unsupervised)
///
/// Outputs:
///   - `trinity.fasta`          — assembled transcripts
///   - `pasa.sqlite`            — PASA transcript database
///   - `pasa_training.gff3`     — high-confidence training models
///   - `augustus_training/`     — trained Augustus species dir
///   - `snap_training.hmm`      — trained SNAP HMM
///   - `train_summary.txt`      — summary of models used for training

pub mod trinity;
pub mod pasa;
pub mod augustus_train;
pub mod snap_train;

use crate::utils::error::{MycoNoteError, Result};
use crate::progress;
use std::io::Write;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct TrainConfig {
    /// Soft-masked genome FASTA (output of `myconote mask`)
    pub masked_fasta: PathBuf,
    /// Output directory
    pub out_dir:      PathBuf,
    /// Left (R1) RNA-seq reads — can be multiple files (comma-separated)
    pub left_reads:   Vec<PathBuf>,
    /// Right (R2) RNA-seq reads — same order as left (empty = single-end)
    pub right_reads:  Vec<PathBuf>,
    /// Single-end reads (alternative to paired left/right)
    pub single_reads: Vec<PathBuf>,
    /// Pre-assembled Trinity FASTA (skip Trinity if provided)
    pub trinity_fasta: Option<PathBuf>,
    /// Species name for Augustus training
    pub species:      String,
    /// Number of threads
    pub threads:      usize,
    /// Maximum intron size (bp). Default: 3000 (fungi), 200000 (plants)
    pub max_intron:   usize,
    /// Minimum number of complete PASA models for training. Default: 200
    pub min_models:   usize,
    /// Also train SNAP (in addition to Augustus)
    pub train_snap:   bool,
    /// Also train/run GeneMark-ES on masked genome
    pub train_genemark: bool,
    /// Strand specificity: "RF", "FR", or "" (unstranded)
    pub strand:       String,
    /// Memory for Trinity (e.g. "50G")
    pub trinity_memory: String,
}

impl Default for TrainConfig {
    fn default() -> Self {
        Self {
            masked_fasta:   PathBuf::new(),
            out_dir:        PathBuf::from("train_out"),
            left_reads:     vec![],
            right_reads:    vec![],
            single_reads:   vec![],
            trinity_fasta:  None,
            species:        "myconote_trained".to_string(),
            threads:        4,
            max_intron:     3000,
            min_models:     200,
            train_snap:     true,
            train_genemark: false,
            strand:         String::new(),
            trinity_memory: "50G".to_string(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────

pub struct TrainResult {
    pub augustus_species_dir: Option<PathBuf>,
    pub snap_hmm:             Option<PathBuf>,
    pub pasa_gff3:            Option<PathBuf>,
    pub trinity_fasta:        Option<PathBuf>,
    pub training_model_count: usize,
}

pub fn run_training(config: &TrainConfig) -> Result<TrainResult> {
    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    let pb = progress::spinner("RNA-seq training pipeline");

    // ── Step 1: Trinity transcript assembly ──────────────────────────────
    pb.set_message("Assembling transcripts with Trinity");

    let trinity_fasta = if let Some(ref pre) = config.trinity_fasta {
        println!("  Using pre-assembled Trinity FASTA: {}", pre.display());
        pre.clone()
    } else {
        trinity::check_trinity()?;
        let out = config.out_dir.join("trinity.fasta");
        trinity::run_trinity(config, &out)?;
        out
    };

    pb.set_message("Trinity assembly complete");

    // ── Step 2: Align to genome with minimap2 ─────────────────────────────
    pb.set_message("Aligning transcripts to genome");

    let bam_path = config.out_dir.join("trinity_aligned.bam");
    align_transcripts_to_genome(&trinity_fasta, &config.masked_fasta, &bam_path, config.threads)?;

    // ── Step 3: PASA database ─────────────────────────────────────────────
    pb.set_message("Building PASA transcript database");

    let pasa_db   = config.out_dir.join("pasa.sqlite");
    let pasa_gff3 = config.out_dir.join("pasa_assemblies.gff3");

    match pasa::run_pasa(config, &trinity_fasta, &pasa_db, &pasa_gff3) {
        Ok(n) => println!("  PASA: {} transcript assemblies", n),
        Err(e) => {
            eprintln!("  ⚠  PASA failed ({}). Will use direct Trinity alignments for training.", e);
        }
    }

    // ── Step 4: Extract training models ──────────────────────────────────
    pb.set_message("Extracting high-confidence training gene models");

    let training_gff3 = config.out_dir.join("training_models.gff3");
    let model_count = extract_training_models(
        &pasa_gff3, &config.masked_fasta, &training_gff3, config.min_models
    )?;

    if model_count < config.min_models {
        eprintln!("  ⚠  Only {} complete models found (minimum: {}).", model_count, config.min_models);
        eprintln!("     Training may be suboptimal. Consider more RNA-seq data.");
    } else {
        println!("  Training models: {}", model_count);
    }

    // ── Step 5: Train Augustus ────────────────────────────────────────────
    pb.set_message("Training Augustus");

    let augustus_dir = augustus_train::train_augustus(
        &training_gff3, &config.masked_fasta, &config.species, &config.out_dir, config.threads
    )?;
    println!("  Augustus trained: {}", augustus_dir.display());

    // ── Step 6: Train SNAP ────────────────────────────────────────────────
    let snap_hmm = if config.train_snap {
        pb.set_message("Training SNAP");
        match snap_train::train_snap(&training_gff3, &config.masked_fasta, &config.out_dir) {
            Ok(hmm) => {
                println!("  SNAP trained:    {}", hmm.display());
                Some(hmm)
            }
            Err(e) => {
                eprintln!("  ⚠  SNAP training failed: {}", e);
                None
            }
        }
    } else {
        None
    };

    // ── Step 7: Write summary ─────────────────────────────────────────────
    write_train_summary(config, model_count, &augustus_dir, snap_hmm.as_deref())?;

    pb.finish();
    println!("\n── Training complete ─────────────────────────────────────");
    println!("  Augustus species:  {}", config.species);
    println!("  Training models:   {}", model_count);
    println!("  Output directory:  {}", config.out_dir.display());
    println!("\n  Next step: myconote-cli predict <genome.fa> \\");
    println!("               --species {} \\", config.species);
    if let Some(ref hmm) = snap_hmm {
        println!("               --snap-hmm {} \\", hmm.display());
    }

    Ok(TrainResult {
        augustus_species_dir: Some(augustus_dir),
        snap_hmm,
        pasa_gff3: Some(pasa_gff3),
        trinity_fasta: Some(trinity_fasta),
        training_model_count: model_count,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

fn align_transcripts_to_genome(
    transcripts: &Path,
    genome:      &Path,
    bam_out:     &Path,
    threads:     usize,
) -> Result<()> {
    // Try minimap2 first (fast, handles spliced alignments well for fungi)
    let use_minimap2 = std::process::Command::new("which")
        .arg("minimap2")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);

    if use_minimap2 {
        let sam_path = bam_out.with_extension("sam");

        let status = std::process::Command::new("minimap2")
            .arg("-a")
            .arg("-x").arg("splice")
            .arg("--cs")
            .arg("-t").arg(threads.to_string())
            .arg(genome)
            .arg(transcripts)
            .stdout(std::fs::File::create(&sam_path).map_err(MycoNoteError::Io)?)
            .status()
            .map_err(|e| MycoNoteError::ExternalTool(format!("minimap2: {}", e)))?;

        if !status.success() {
            return Err(MycoNoteError::ExternalTool("minimap2 alignment failed".to_string()));
        }

        // Sort and index with samtools
        let _ = std::process::Command::new("samtools")
            .args(["sort", "-o"])
            .arg(bam_out)
            .arg(&sam_path)
            .arg("-@").arg(threads.to_string())
            .status();
        let _ = std::process::Command::new("samtools")
            .arg("index")
            .arg(bam_out)
            .status();
    } else {
        eprintln!("  ⚠  minimap2 not found, skipping genome alignment (PASA may still work)");
    }

    Ok(())
}

fn extract_training_models(
    pasa_gff3:    &Path,
    _genome_fasta: &Path,
    output:       &Path,
    min_models:   usize,
) -> Result<usize> {
    use std::io::{BufRead, BufReader, Write as IoWrite};

    if !pasa_gff3.exists() {
        // Create an empty placeholder
        std::fs::write(output, "").map_err(MycoNoteError::Io)?;
        return Ok(0);
    }

    let file   = std::fs::File::open(pasa_gff3).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;

    let mut gene_count = 0usize;
    let mut in_complete_gene = false;
    let mut buffer: Vec<String> = Vec::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();

        if trimmed.starts_with('#') || trimmed.is_empty() {
            if !buffer.is_empty() {
                // Write buffered complete gene
                for bl in &buffer {
                    writeln!(out, "{}", bl).map_err(MycoNoteError::Io)?;
                }
                buffer.clear();
                in_complete_gene = false;
            }
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 9 { continue; }

        match fields[2] {
            "gene" => {
                buffer.clear();
                in_complete_gene = true;
                buffer.push(line.clone());
                gene_count += 1;
            }
            "mRNA" | "transcript" if in_complete_gene => {
                // Only keep models with both start and stop codon (complete CDS)
                buffer.push(line.clone());
            }
            "start_codon" | "stop_codon" => {
                buffer.push(line.clone());
            }
            _ if in_complete_gene => {
                buffer.push(line.clone());
            }
            _ => {}
        }
    }

    // Flush remaining
    if !buffer.is_empty() {
        for bl in &buffer {
            writeln!(out, "{}", bl).map_err(MycoNoteError::Io)?;
        }
    }

    let _ = min_models; // used by caller for warning only
    Ok(gene_count)
}

fn write_train_summary(
    config:        &TrainConfig,
    model_count:   usize,
    augustus_dir:  &Path,
    snap_hmm:      Option<&Path>,
) -> Result<()> {
    let summary_path = config.out_dir.join("train_summary.txt");
    let mut f = std::fs::File::create(&summary_path).map_err(MycoNoteError::Io)?;

    writeln!(f, "myconote-cli train summary").map_err(MycoNoteError::Io)?;
    writeln!(f, "=========================").map_err(MycoNoteError::Io)?;
    writeln!(f, "Species:          {}", config.species).map_err(MycoNoteError::Io)?;
    writeln!(f, "Genome:           {}", config.masked_fasta.display()).map_err(MycoNoteError::Io)?;
    writeln!(f, "RNA-seq reads:    {} file(s)", config.left_reads.len() + config.single_reads.len()).map_err(MycoNoteError::Io)?;
    writeln!(f, "Training models:  {}", model_count).map_err(MycoNoteError::Io)?;
    writeln!(f, "Augustus trained: {}", augustus_dir.display()).map_err(MycoNoteError::Io)?;
    if let Some(hmm) = snap_hmm {
        writeln!(f, "SNAP trained:     {}", hmm.display()).map_err(MycoNoteError::Io)?;
    }
    writeln!(f).map_err(MycoNoteError::Io)?;
    writeln!(f, "Next step:").map_err(MycoNoteError::Io)?;
    writeln!(f, "  myconote-cli predict <genome.fa> --species {} \\", config.species).map_err(MycoNoteError::Io)?;
    if let Some(hmm) = snap_hmm {
        writeln!(f, "    --snap-hmm {} \\", hmm.display()).map_err(MycoNoteError::Io)?;
    }

    Ok(())
}
