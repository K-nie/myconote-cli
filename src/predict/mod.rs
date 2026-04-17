pub mod augustus;
pub mod evidence;
pub mod genemark;
pub mod glimmer;
/// Gene prediction pipeline
///
/// Orchestrates the full predict workflow:
///   1. Run Augustus (ab initio, kingdom-aware)
///   2. Run SNAP (secondary ab initio, optional)
///   3. Merge predictions with Evidence Modeler-style scoring
///   4. Write consensus GFF3 with sequential locus tags
///
/// External tools used (all optional — graceful fallback):
///   - augustus   (conda install -c bioconda augustus)
///   - snap        (conda install -c bioconda snap)
pub mod kingdom;
pub mod ploidy;
pub mod protein_evidence;
pub mod snap;
pub mod train;

use crate::progress;
use crate::utils::error::{MycoNoteError, Result};
use kingdom::Kingdom;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct PredictConfig {
    /// Soft-masked FASTA (output of `myconote mask`)
    pub masked_fasta: PathBuf,
    /// Output directory for all predict files
    pub out_dir: PathBuf,
    /// Organism kingdom — drives species defaults and weights
    pub kingdom: Kingdom,
    /// Augustus species override (None = use kingdom default)
    pub augustus_species: Option<String>,
    /// Use SNAP in addition to Augustus
    pub use_snap: bool,
    /// SNAP HMM parameter file (None = use kingdom default)
    pub snap_hmm: Option<String>,
    /// Optional protein BLAST hints (tabular fmt6)
    pub protein_evidence: Option<PathBuf>,
    /// Gene ID prefix for locus tags (e.g. "MYCO" → "MYCO_000001")
    pub locus_prefix: String,
    /// Parallel threads for each tool
    pub threads: usize,
    /// Evidence weights (None = use defaults: augustus=10, snap=3, protein=20)
    pub weights: Option<evidence::EvidenceWeights>,
    /// Auto-train Augustus on a high-confidence subset before predicting
    pub self_train: bool,
    /// Species name for trained model (default: "<locus_prefix>_trained")
    pub train_species: Option<String>,
    /// Run GlimmerHMM (adds a third ab initio predictor)
    pub use_glimmerhmm: bool,
    /// GlimmerHMM training directory (None = auto-detect from kingdom)
    pub glimmer_dir: Option<PathBuf>,
    /// Run GeneMark-ES (self-training, no pre-existing model needed)
    pub use_genemark: bool,
    /// Use GeneMark-ET (RNA-seq intron hints file in GFF format)
    pub genemark_hints: Option<PathBuf>,
    /// Protein FASTA for protein→genome evidence (miniprot/exonerate)
    pub protein_fasta: Option<PathBuf>,
    /// Maximum intron size for protein alignment
    pub max_intron: usize,
    /// Ploidy level (None = auto-detect or haploid)
    pub ploidy: Option<u8>,
    /// Evidence weights TOML file (None = use defaults)
    pub weights_file: Option<PathBuf>,
}

impl Default for PredictConfig {
    fn default() -> Self {
        Self {
            masked_fasta: PathBuf::new(),
            out_dir: PathBuf::from("predict_out"),
            kingdom: Kingdom::Fungi,
            augustus_species: None,
            use_snap: true,
            snap_hmm: None,
            protein_evidence: None,
            locus_prefix: "GENE".to_string(),
            threads: 4,
            weights: None,
            self_train: false,
            train_species: None,
            use_glimmerhmm: false,
            glimmer_dir: None,
            use_genemark: false,
            genemark_hints: None,
            protein_fasta: None,
            max_intron: 10_000,
            ploidy: None,
            weights_file: None,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SNAP HMM defaults per kingdom
// ─────────────────────────────────────────────────────────────────────────────

fn default_snap_hmm(k: &Kingdom) -> &'static str {
    match k {
        Kingdom::Fungi => "fungal",
        Kingdom::Plant => "worm", // closest available for plants
        Kingdom::Animal => "human",
        Kingdom::Insect => "fly",
        Kingdom::Protist => "fungal",
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Pipeline entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run the full gene prediction pipeline.
/// Returns the path to the final consensus GFF3 and the number of genes called.
pub fn run_prediction(config: &PredictConfig) -> Result<(PathBuf, usize)> {
    // ── Validate input ────────────────────────────────────────────────────────
    if !config.masked_fasta.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Masked FASTA not found: {}",
            config.masked_fasta.display()
        )));
    }

    // ── Create output directory ───────────────────────────────────────────────
    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    println!("── Gene prediction ──────────────────────────────────────────");
    println!("  Kingdom : {}", config.kingdom.display_name());
    println!("  Input   : {}", config.masked_fasta.display());
    println!("  Output  : {}", config.out_dir.display());

    let mut prediction_inputs: Vec<(PathBuf, &'static str, f64)> = Vec::new();

    // Load weights from file if specified, otherwise use config or defaults
    let weights = if let Some(ref wf) = config.weights_file {
        match evidence::EvidenceWeights::from_toml(wf) {
            Ok(w) => {
                println!("  Weights loaded from: {}", wf.display());
                w
            }
            Err(e) => {
                eprintln!(
                    "  Warning: failed to load weights file: {}. Using defaults.",
                    e
                );
                config.weights.clone().unwrap_or_default()
            }
        }
    } else {
        config.weights.clone().unwrap_or_default()
    };

    // Save weights used for reproducibility
    let weights_out = config.out_dir.join("evidence_weights.toml");
    let _ = weights.write_toml(&weights_out);

    // ── 0. Optional self-training ─────────────────────────────────────────────
    let trained_species: Option<String> = if config.self_train {
        progress::step(1, 4, "Augustus self-training (first-pass prediction)…");
        match run_self_training(config) {
            Ok(name) => {
                println!("  ✓ Trained species: {}", name);
                Some(name)
            }
            Err(e) => {
                eprintln!("  ⚠  Self-training failed: {}", e);
                eprintln!("     Continuing with pre-trained species model.");
                None
            }
        }
    } else {
        None
    };

    let step_offset = if config.self_train { 1 } else { 0 };

    // ── 1. Augustus ───────────────────────────────────────────────────────────
    let aug_gff = config.out_dir.join("augustus.gff3");
    // Use trained species if available, otherwise configured/default species
    let aug_species_owned = trained_species
        .clone()
        .or_else(|| config.augustus_species.clone())
        .unwrap_or_else(|| config.kingdom.default_augustus_species().to_string());
    let aug_species = aug_species_owned.as_str();

    // Optionally build protein hints first
    let hints_path = if let Some(ref blast_tsv) = config.protein_evidence {
        let hp = config.out_dir.join("protein_hints.gff");
        println!("  Building protein hints from {}…", blast_tsv.display());
        augustus::make_protein_hints(blast_tsv, &hp, 4)?;
        Some(hp)
    } else {
        None
    };

    let aug_cfg = augustus::AugustusConfig {
        species: aug_species.to_string(),
        threads: config.threads,
        utr: config.kingdom.augustus_utr(),
        hints_file: hints_path,
        extra_args: Vec::new(),
    };

    progress::step(
        step_offset + 1,
        step_offset + 3,
        &format!("Augustus ({})…", aug_species),
    );
    let pb_aug = progress::spinner(format!("Running augustus --species={}…", aug_species));
    match augustus::run(&config.masked_fasta, &aug_gff, &aug_cfg) {
        Ok(()) => {
            progress::finish_spinner(
                &pb_aug,
                format!("Augustus complete → {}", aug_gff.display()),
            );
            prediction_inputs.push((aug_gff, "Augustus", weights.augustus));
        }
        Err(e) => {
            progress::warn_spinner(&pb_aug, format!("Augustus failed: {}", e));
            eprintln!("     At least one predictor must succeed.");
        }
    }

    // ── 2. SNAP (optional) ────────────────────────────────────────────────────
    if config.use_snap {
        let snap_gff = config.out_dir.join("snap.gff3");
        let hmm = config
            .snap_hmm
            .as_deref()
            .unwrap_or_else(|| default_snap_hmm(&config.kingdom));

        let snap_cfg = snap::SnapConfig {
            hmm: hmm.to_string(),
            threads: config.threads,
        };

        progress::step(step_offset + 2, step_offset + 3, "SNAP gene prediction…");
        let pb_snap = progress::spinner(format!("Running SNAP (hmm: {})…", hmm));
        match snap::run(&config.masked_fasta, &snap_gff, &snap_cfg) {
            Ok(()) => {
                progress::finish_spinner(&pb_snap, "SNAP complete");
                prediction_inputs.push((snap_gff, "SNAP", weights.snap));
            }
            Err(e) => {
                progress::warn_spinner(&pb_snap, format!("SNAP failed (non-fatal): {}", e));
            }
        }
    }

    // ── 2a. GlimmerHMM (optional third ab-initio predictor) ─────────────────
    // The audit found that `use_glimmerhmm` was parsed from the CLI and set
    // on PredictConfig but never read here — the flag was a silent no-op.
    // Now dispatched like SNAP: non-fatal on failure, merges into EVM as
    // its own weighted source.
    if config.use_glimmerhmm {
        let glimmer_gff = config.out_dir.join("glimmerhmm.gff3");
        let pb = progress::spinner("Running GlimmerHMM…");
        match config.glimmer_dir.as_ref() {
            Some(train_dir) if train_dir.exists() => {
                match glimmer::run_glimmerhmm(
                    &config.masked_fasta,
                    train_dir,
                    &glimmer_gff,
                    config.threads,
                ) {
                    Ok(n) => {
                        progress::finish_spinner(&pb, format!("GlimmerHMM: {} genes", n));
                        prediction_inputs.push((glimmer_gff, "GlimmerHMM", weights.glimmerhmm));
                    }
                    Err(e) => {
                        progress::warn_spinner(&pb, format!("GlimmerHMM failed (non-fatal): {}", e));
                    }
                }
            }
            Some(train_dir) => {
                progress::warn_spinner(
                    &pb,
                    format!("GlimmerHMM training dir not found: {}", train_dir.display()),
                );
            }
            None => {
                progress::warn_spinner(
                    &pb,
                    "GlimmerHMM requested but no training dir. Pass --glimmer-dir <trained_dir> \
                     (glimmerhmm ships example trainings under /usr/share/glimmerhmm/trained_dir)",
                );
            }
        }
    }

    // ── 2b. GeneMark-ES / ET (optional self-training predictor) ─────────────
    // Same silent-no-op bug: `use_genemark` and `genemark_hints` were set but
    // never dispatched. GeneMark-ES is self-training; GeneMark-ET adds RNA-seq
    // splice hints when `--genemark-hints <intron_hints.gff>` is supplied.
    if config.use_genemark {
        let gm_dir = config.out_dir.join("genemark");
        let is_fungus = matches!(config.kingdom, Kingdom::Fungi);
        let pb = progress::spinner("Running GeneMark-ES/ET…");
        let run_result = if let Some(hints) = config.genemark_hints.as_ref() {
            if hints.exists() {
                genemark::run_genemark_et(
                    &config.masked_fasta,
                    hints,
                    &gm_dir,
                    is_fungus,
                    config.threads,
                )
            } else {
                progress::warn_spinner(
                    &pb,
                    format!(
                        "GeneMark hints file not found: {} (falling back to --ES)",
                        hints.display()
                    ),
                );
                genemark::run_genemark_es(&config.masked_fasta, &gm_dir, is_fungus, config.threads)
            }
        } else {
            genemark::run_genemark_es(&config.masked_fasta, &gm_dir, is_fungus, config.threads)
        };
        match run_result {
            Ok(gm_gff) => {
                progress::finish_spinner(&pb, format!("GeneMark complete → {}", gm_gff.display()));
                prediction_inputs.push((gm_gff, "GeneMark", weights.genemark));
            }
            Err(e) => {
                progress::warn_spinner(&pb, format!("GeneMark failed (non-fatal): {}", e));
            }
        }
    }

    // ── 2c. Protein evidence (miniprot / exonerate) ─────────────────────────
    if let Some(ref prot_fa) = config.protein_fasta {
        if prot_fa.exists() {
            let prot_cfg = protein_evidence::ProteinEvidenceConfig {
                proteins: prot_fa.clone(),
                genome: config.masked_fasta.clone(),
                out_dir: config.out_dir.join("protein_evidence"),
                threads: config.threads,
                max_intron: config.max_intron,
                ..protein_evidence::ProteinEvidenceConfig::default()
            };

            println!(
                "  Running protein-to-genome alignment ({})...",
                prot_fa.display()
            );
            match protein_evidence::generate_protein_evidence(&prot_cfg) {
                Ok(result) => {
                    println!(
                        "  {} proteins aligned to {} loci ({})",
                        result.n_aligned, result.n_loci, result.tool
                    );
                    prediction_inputs.push((result.evidence_gff, "Protein", weights.protein));
                }
                Err(e) => {
                    eprintln!("  Warning: protein evidence alignment failed: {}", e);
                }
            }
        } else {
            eprintln!("  Warning: protein FASTA not found: {}", prot_fa.display());
        }
    }

    // ── Guard: at least one predictor must have succeeded ─────────────────────
    if prediction_inputs.is_empty() {
        return Err(MycoNoteError::InvalidFormat(
            "All gene predictors failed. Check that Augustus (and optionally SNAP) \
             are installed and that the species model is correct."
                .to_string(),
        ));
    }

    // ── 3. Evidence Modeler consensus ─────────────────────────────────────────
    progress::step(
        step_offset + 3,
        step_offset + 3,
        "Merging predictions (Evidence Modeler)…",
    );
    let consensus_gff = config.out_dir.join("consensus.gff3");

    let inputs_ref: Vec<(&Path, &str, f64)> = prediction_inputs
        .iter()
        .map(|(p, s, w)| (p.as_path(), *s, *w))
        .collect();

    let gene_count =
        evidence::merge_predictions(&inputs_ref, &consensus_gff, &config.locus_prefix)?;

    println!("  ✓  Consensus GFF3 → {}", consensus_gff.display());
    println!("     {} genes called", gene_count);

    // ── 4. Write a short summary ──────────────────────────────────────────────
    let summary_path = config.out_dir.join("predict_summary.txt");
    write_summary(&summary_path, config, gene_count, &prediction_inputs)?;

    Ok((consensus_gff, gene_count))
}

// ─────────────────────────────────────────────────────────────────────────────
// Summary report
// ─────────────────────────────────────────────────────────────────────────────

fn write_summary(
    path: &Path,
    config: &PredictConfig,
    gene_count: usize,
    inputs: &[(PathBuf, &'static str, f64)],
) -> Result<()> {
    use std::io::Write;

    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "myconote predict — Summary").map_err(MycoNoteError::Io)?;
    writeln!(f, "==========================").map_err(MycoNoteError::Io)?;
    writeln!(f, "Input FASTA  : {}", config.masked_fasta.display()).map_err(MycoNoteError::Io)?;
    writeln!(f, "Kingdom      : {}", config.kingdom.display_name()).map_err(MycoNoteError::Io)?;
    writeln!(f, "Locus prefix : {}", config.locus_prefix).map_err(MycoNoteError::Io)?;
    writeln!(f, "").map_err(MycoNoteError::Io)?;
    writeln!(f, "Predictors used:").map_err(MycoNoteError::Io)?;
    for (path, source, weight) in inputs {
        writeln!(
            f,
            "  {:<12} weight={:.0}  → {}",
            source,
            weight,
            path.display()
        )
        .map_err(MycoNoteError::Io)?;
    }
    writeln!(f, "").map_err(MycoNoteError::Io)?;
    writeln!(f, "Final gene count : {}", gene_count).map_err(MycoNoteError::Io)?;

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Self-training helper
// ─────────────────────────────────────────────────────────────────────────────

/// Run a first-pass prediction with the default model, then use those gene
/// models to train a custom Augustus species for the actual prediction.
/// Returns the trained species name.
fn run_self_training(config: &PredictConfig) -> Result<String> {
    let train_dir = config.out_dir.join("training");
    std::fs::create_dir_all(&train_dir).map_err(MycoNoteError::Io)?;

    let default_species = config
        .augustus_species
        .as_deref()
        .unwrap_or_else(|| config.kingdom.default_augustus_species());

    // First-pass Augustus prediction for training material
    let first_pass_gff = train_dir.join("first_pass.gff3");
    let aug_cfg = augustus::AugustusConfig {
        species: default_species.to_string(),
        threads: config.threads,
        utr: false, // no UTR for training pass
        hints_file: None,
        extra_args: Vec::new(),
    };
    augustus::run(&config.masked_fasta, &first_pass_gff, &aug_cfg)?;

    // Train a custom species from the first-pass predictions
    let species_name = config
        .train_species
        .clone()
        .unwrap_or_else(|| format!("{}_trained", config.locus_prefix.to_lowercase()));

    let train_config = train::TrainConfig {
        gff: first_pass_gff,
        fasta: config.masked_fasta.clone(),
        species_name: species_name.clone(),
        out_dir: train_dir,
        optimize: false, // skip optimization for speed
        threads: config.threads,
        ..train::TrainConfig::default()
    };

    train::run_training(&train_config)?;
    Ok(species_name)
}
