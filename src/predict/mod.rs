pub mod augustus;
pub mod braker;
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
use genemark::GeneMarkMode;
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
    /// Run GeneMark in any of its self-training modes.  When `None`, GeneMark
    /// is skipped.  Modes: ES (no evidence), ET (RNA-seq introns), EP+
    /// (ProtHint-derived protein hints), ETP+ (both).
    pub genemark_mode: Option<GeneMarkMode>,
    /// Legacy alias kept for backward compatibility — `--genemark` sets this
    /// and is translated to `genemark_mode = Some(GeneMarkMode::Es)` if no
    /// explicit `--genemark-mode` is given.
    pub use_genemark: bool,
    /// RNA-seq intron hints in GFF format (HISAT2 / STAR splice output).
    /// Required for GeneMark-ET and GeneMark-ETP+.  When supplied without
    /// `--genemark-mode`, mode defaults to ET (legacy `--genemark-hints`).
    pub genemark_hints: Option<PathBuf>,
    /// Protein FASTA for protein→genome evidence (miniprot / exonerate) and
    /// for ProtHint when running GeneMark-EP+ / ETP+.
    pub protein_fasta: Option<PathBuf>,
    /// Maximum intron size for protein alignment
    pub max_intron: usize,
    /// Ploidy level (None = auto-detect or haploid)
    pub ploidy: Option<u8>,
    /// Evidence weights TOML file (None = use defaults)
    pub weights_file: Option<PathBuf>,
    // ── BRAKER (v0.6.0) ──────────────────────────────────────────────────────
    /// When true, hand the entire prediction over to BRAKER and skip the
    /// standard ab-initio + EVM consensus stack.  Mutually exclusive with the
    /// per-tool predictor flags (Augustus / SNAP / GlimmerHMM / GeneMark /
    /// miniprot are not invoked when BRAKER is the engine).
    pub use_braker: bool,
    /// BRAKER1/2/3 mode override.  `None` means auto-detect from the inputs:
    /// RNA-only -> 1, protein-only -> 2, both -> 3.
    pub braker_mode: Option<braker::BrakerMode>,
    /// RNA-seq BAM(s) for BRAKER (`--bam`, repeatable).
    pub braker_rna_bams: Vec<PathBuf>,
    /// Protein FASTA for BRAKER (`--prot_seq`, typically OrthoDB fungi).
    pub braker_proteins: Option<PathBuf>,
    /// NCBI translation table forwarded to BRAKER as `--translation_table`
    /// and through Augustus / GeneMark.  Only meaningful when `use_braker`.
    pub genetic_code: u8,
    // ── Consensus false-positive filter (v0.7.7) ─────────────────────────────
    /// Minimum number of distinct predictor sources that must corroborate a
    /// consensus locus for it to be emitted.  `1` (default) keeps every locus
    /// and reproduces pre-filter behaviour; `2` drops loci called by a single
    /// predictor alone.  See `evidence::ConsensusFilter`.
    pub consensus_min_support: usize,
    /// Optional rescue threshold on the agreement-weighted consensus score: a
    /// locus failing the support test is still kept when its score reaches this
    /// value.  `None` (default) disables the rescue.
    pub consensus_min_score: Option<f64>,
    // ── Augustus training mode (v0.7.7) ──────────────────────────────────────
    // (self-training itself is driven by `self_train`; `--augustus-training
    // self` is an alias that sets it.  No extra field needed.)
    // ── Protein → Augustus hints, BRAKER-style (v0.7.7) ──────────────────────
    /// Protein FASTA aligned with miniprot to build Augustus extrinsic hints
    /// that steer Augustus's own ab-initio calls (distinct from
    /// `protein_fasta`, which adds protein alignments as a separate EVM track,
    /// and from `protein_evidence`, which consumes a precomputed BLAST table).
    /// `None` (default) → ab-initio only, unchanged behaviour.
    pub protein_hints: Option<PathBuf>,
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
            genemark_mode: None,
            use_genemark: false,
            genemark_hints: None,
            protein_fasta: None,
            max_intron: 10_000,
            ploidy: None,
            weights_file: None,
            use_braker: false,
            braker_mode: None,
            braker_rna_bams: Vec::new(),
            braker_proteins: None,
            genetic_code: 1,
            consensus_min_support: 1,
            consensus_min_score: None,
            protein_hints: None,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GeneMark-mode resolution (handles legacy `--genemark` / `--genemark-hints`)
// ─────────────────────────────────────────────────────────────────────────────

/// Translate the deprecated boolean / hints-only flags into the new
/// `GeneMarkMode` enum.  Resolution rules, in order:
///
///   1. If `genemark_mode` is set explicitly, use it as-is.
///   2. Else if `use_genemark` is true and `genemark_hints` is `Some`, ET.
///   3. Else if `use_genemark` is true alone, ES.
///   4. Else if `genemark_hints` is set without `use_genemark`, ET (preserves
///      the old behaviour where supplying hints implied GeneMark-ET).
///   5. Else `None` — GeneMark is skipped.
///
/// EP / ETP modes can only be requested via `--genemark-mode`; there is no
/// legacy alias for them since they did not exist in earlier releases.
pub fn resolve_genemark_mode(
    explicit: Option<GeneMarkMode>,
    use_genemark: bool,
    genemark_hints: Option<&Path>,
) -> Option<GeneMarkMode> {
    if let Some(m) = explicit {
        return Some(m);
    }
    if use_genemark {
        if genemark_hints.is_some() {
            return Some(GeneMarkMode::Et);
        }
        return Some(GeneMarkMode::Es);
    }
    if genemark_hints.is_some() {
        return Some(GeneMarkMode::Et);
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// BRAKER conflict check
// ─────────────────────────────────────────────────────────────────────────────

/// Reject predict configs where `--use-braker` is combined with the standard
/// per-predictor flags.  BRAKER is a complete alternative engine; the EVM
/// stack does not run when it's on, so any per-predictor flag is at best
/// dead weight and at worst confusing (e.g. `--genemark-mode etp` while
/// BRAKER also runs its own GeneMark-ETP+ internally).
///
/// This is a pure function so it can be unit-tested without spinning up
/// the whole pipeline.
pub fn check_braker_conflicts(config: &PredictConfig) -> Result<()> {
    if !config.use_braker {
        return Ok(());
    }
    let mut conflicts: Vec<&'static str> = Vec::new();
    if config.genemark_mode.is_some() {
        conflicts.push("--genemark-mode");
    }
    if config.use_genemark {
        conflicts.push("--genemark");
    }
    if config.genemark_hints.is_some() {
        conflicts.push("--genemark-hints");
    }
    if config.protein_evidence.is_some() {
        conflicts.push("--protein-evidence");
    }
    if config.protein_hints.is_some() {
        conflicts.push("--protein-hints");
    }
    if config.protein_fasta.is_some() {
        conflicts.push("--protein-fasta");
    }
    if config.use_glimmerhmm {
        conflicts.push("--glimmerhmm");
    }
    // Note: --no-snap is NOT a conflict — it's a way to disable a predictor
    // we wouldn't run anyway, so silently allow it.
    if conflicts.is_empty() {
        Ok(())
    } else {
        Err(MycoNoteError::InvalidFormat(format!(
            "--use-braker is mutually exclusive with the standard predictor \
             flags. Drop these and re-run: {}. \
             Pass RNA-seq + protein evidence via --braker-rna-bam and \
             --braker-proteins instead.",
            conflicts.join(", ")
        )))
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

    // ── Reject predict configs that mix BRAKER with the standard stack ───────
    check_braker_conflicts(config)?;

    // ── Create output directory ───────────────────────────────────────────────
    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    // ── BRAKER short-circuit ─────────────────────────────────────────────────
    // When `--use-braker` is on, hand the entire prediction over to braker.pl
    // and skip the EVM consensus stack entirely.  The downstream stages
    // (`update`, `annotate`, `submit`) consume `braker.gff3` exactly the way
    // they consume the EVM-generated `consensus.gff3`, so we just publish
    // braker.gff3 as the canonical predict output and return.
    if config.use_braker {
        println!("── Gene prediction (BRAKER) ─────────────────────────────────");
        println!("  Kingdom : {}", config.kingdom.display_name());
        println!("  Input   : {}", config.masked_fasta.display());
        println!("  Output  : {}", config.out_dir.display());

        let species = config
            .augustus_species
            .clone()
            .unwrap_or_else(|| format!("{}_braker", config.locus_prefix.to_lowercase()));

        let braker_cfg = braker::BrakerConfig {
            genome: config.masked_fasta.clone(),
            rna_bams: config.braker_rna_bams.clone(),
            proteins: config.braker_proteins.clone(),
            out_dir: config.out_dir.join("braker"),
            species,
            threads: config.threads,
            genetic_code: config.genetic_code,
            mode: config.braker_mode,
            fungus: matches!(config.kingdom, Kingdom::Fungi),
        };

        let braker_gff = braker::run_braker(&braker_cfg)?;
        // Publish under the canonical name downstream stages expect, so the
        // BRAKER and EVM paths are interchangeable from the user's POV.
        let consensus_gff = config.out_dir.join("consensus.gff3");
        std::fs::copy(&braker_gff, &consensus_gff).map_err(MycoNoteError::Io)?;

        let gene_count = count_genes_in_gff3(&consensus_gff)?;
        println!("  ✓  BRAKER consensus → {}", consensus_gff.display());
        println!("     {} genes called", gene_count);

        // Reuse the existing summary writer with a single-source list so the
        // output shape is identical to the EVM path.
        let summary_path = config.out_dir.join("predict_summary.txt");
        let inputs = vec![(braker_gff.clone(), "BRAKER", 0.0)];
        write_summary(&summary_path, config, gene_count, &inputs)?;

        return Ok((consensus_gff, gene_count));
    }

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

    // ── Augustus extrinsic hints ─────────────────────────────────────────────
    // Two independent hint sources can feed Augustus, and both are optional:
    //   * --protein-evidence <blast.tsv>  (existing)  → CDSpart hints from a
    //     precomputed BLAST table (augustus::make_protein_hints).
    //   * --protein-hints   <proteins.fa> (v0.7.7, BRAKER-style) → run miniprot
    //     protein→genome and convert the alignment into Augustus hints.
    // Whatever is available is concatenated into one hints file. The presence of
    // protein-FASTA hints also selects a P-aware extrinsic config so Augustus
    // actually weights the `src=P` lines. With neither flag, Augustus runs
    // ab-initio exactly as before (hints_path = None, extrinsic_cfg = None).
    let mut hint_files: Vec<PathBuf> = Vec::new();

    if let Some(ref blast_tsv) = config.protein_evidence {
        let hp = config.out_dir.join("protein_hints.gff");
        println!("  Building protein hints from {}…", blast_tsv.display());
        augustus::make_protein_hints(blast_tsv, &hp, 4)?;
        hint_files.push(hp);
    }

    if let Some(ref prot_fa) = config.protein_hints {
        // Non-fatal: a missing protein file or missing aligner logs and leaves
        // Augustus ab-initio rather than failing the whole run.
        if let Some(hp) = build_protein_hints_miniprot(config, prot_fa)? {
            hint_files.push(hp);
        }
    }

    let use_protein_extrinsic = config.protein_hints.is_some() && !hint_files.is_empty();
    let hints_path = combine_hint_files(&hint_files, &config.out_dir)?;
    let extrinsic_cfg = augustus_extrinsic_cfg(use_protein_extrinsic && hints_path.is_some());

    let aug_cfg = augustus::AugustusConfig {
        species: aug_species.to_string(),
        threads: config.threads,
        utr: config.kingdom.augustus_utr(),
        hints_file: hints_path,
        extrinsic_cfg,
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
            prediction_inputs.push((aug_gff.clone(), "Augustus", weights.augustus));
        }
        Err(e) => {
            progress::warn_spinner(&pb_aug, format!("Augustus failed: {}", e));
            eprintln!("     At least one predictor must succeed.");
        }
    }

    // ── 2. SNAP (optional) ────────────────────────────────────────────────────
    //
    // SNAP needs an HMM parameter file and there is NO stock fungal HMM shipped
    // with the `snap` conda package (its share/snap/HMM dir carries worm / fly /
    // human / plant models, not a "fungal" one). The old code fell back to the
    // literal string "fungal", which does not resolve, so SNAP errored and the
    // error was swallowed by `warn_spinner` — invisible in non-TTY cluster logs.
    // The whole benchmark therefore ran Augustus-only.
    //
    // New behaviour: an explicit --snap-hmm still wins. Otherwise self-train a
    // genome-specific HMM from the Augustus first-pass calls (this is exactly
    // what funannotate does) and feed that into SNAP. Every branch logs via
    // println! so non-TTY logs show whether SNAP ran and why not.
    if config.use_snap {
        println!("  [+] SNAP (secondary ab-initio predictor)…");
        let snap_gff = config.out_dir.join("snap.gff3");

        let resolved_hmm: Option<String> = match config.snap_hmm.as_deref() {
            Some(h) => {
                println!("      using provided HMM: {}", h);
                Some(h.to_string())
            }
            None => {
                if !aug_gff.exists() {
                    println!("      ⚠ SNAP skipped: no Augustus first-pass GFF to self-train on.");
                    None
                } else if !crate::train::snap_train::snap_available() {
                    println!("      ⚠ SNAP skipped: `snap` not found on PATH.");
                    None
                } else {
                    println!("      self-training SNAP HMM on Augustus first-pass calls…");
                    match crate::train::snap_train::train_snap(
                        &aug_gff,
                        &config.masked_fasta,
                        &config.out_dir,
                    ) {
                        Ok(hmm_path) => {
                            println!("      SNAP HMM trained → {}", hmm_path.display());
                            Some(hmm_path.to_string_lossy().into_owned())
                        }
                        Err(e) => {
                            println!("      ⚠ SNAP self-training failed (non-fatal): {}", e);
                            None
                        }
                    }
                }
            }
        };

        if let Some(hmm) = resolved_hmm {
            let snap_cfg = snap::SnapConfig {
                hmm,
                threads: config.threads,
            };
            match snap::run(&config.masked_fasta, &snap_gff, &snap_cfg) {
                Ok(()) => {
                    println!("      SNAP complete → {}", snap_gff.display());
                    prediction_inputs.push((snap_gff, "SNAP", weights.snap));
                }
                Err(e) => {
                    println!("      ⚠ SNAP failed (non-fatal): {}", e);
                }
            }
        }
    }

    // ── 2a. GlimmerHMM (optional third ab-initio predictor) ─────────────────
    //
    // Prior audit found that `use_glimmerhmm` was parsed from the CLI and
    // set on PredictConfig but never read here — a silent no-op. The next
    // audit (v0.7.2) found that even after the dispatch was wired, the
    // outcome went through `progress::warn_spinner`, which is invisible in
    // non-TTY (batch / cluster) runs. Callers piping stdout to a log file
    // therefore had no way to see that GlimmerHMM was skipped or why.
    //
    // Current behavior: always emit a `println!` line for the step and its
    // outcome so non-TTY logs make the decision visible, and auto-detect
    // the training dir when the user did not pass `--glimmer-dir`.
    if config.use_glimmerhmm {
        println!("  [+] GlimmerHMM (third ab-initio predictor)…");
        let glimmer_gff = config.out_dir.join("glimmerhmm.gff3");

        // Resolve the training dir: prefer explicit --glimmer-dir; otherwise
        // fall back to `find_training_dir(default_species_for_kingdom)`.
        let resolved_dir: Option<std::path::PathBuf> = match config.glimmer_dir.as_ref() {
            Some(p) if p.exists() => Some(p.clone()),
            Some(p) => {
                println!(
                    "      ⚠ --glimmer-dir path does not exist: {}. \
                     Trying auto-detect for kingdom={}.",
                    p.display(),
                    config.kingdom.display_name()
                );
                let species = glimmer::default_training_species(config.kingdom.display_name());
                glimmer::find_training_dir(species)
            }
            None => {
                let species = glimmer::default_training_species(config.kingdom.display_name());
                glimmer::find_training_dir(species)
            }
        };

        match resolved_dir {
            Some(train_dir) => {
                println!("      training dir: {}", train_dir.display());
                match glimmer::run_glimmerhmm(
                    &config.masked_fasta,
                    &train_dir,
                    &glimmer_gff,
                    config.threads,
                ) {
                    Ok(n) => {
                        println!("      GlimmerHMM finished → {} genes", n);
                        prediction_inputs.push((glimmer_gff, "GlimmerHMM", weights.glimmerhmm));
                    }
                    Err(e) => {
                        println!("      ⚠ GlimmerHMM failed (non-fatal): {}", e);
                    }
                }
            }
            None => {
                println!(
                    "      ⚠ GlimmerHMM training dir not found. Pass \
                     --glimmer-dir <trained_dir> explicitly. Common locations \
                     that were checked: /usr/share/glimmerhmm/trained_dir, \
                     /opt/conda/share/glimmerhmm/trained_dir, \
                     /usr/local/share/glimmerhmm/trained_dir, ~/.myconote/glimmerhmm."
                );
            }
        }
    }

    // ── 2b. GeneMark-ES / ET / EP+ / ETP+ (optional) ─────────────────────────
    //
    // The actual mode is resolved from three flag sources by
    // `resolve_genemark_mode`: the new `--genemark-mode <es|et|ep|etp>`, the
    // deprecated `--genemark` boolean, and the `--genemark-hints <gff>` path.
    // EP / ETP additionally require `--protein-fasta` because ProtHint needs
    // a protein database to convert into hints.
    // Resolve the GeneMark mode from the explicit flags, then fall back to a
    // fungi default: GeneMark-ES self-trains with no external evidence and is a
    // genuinely independent predictor (unlike SNAP, which we train off Augustus
    // above), so for fungi we enable it by default whenever the GeneMark-ES
    // suite is on PATH. It stays off (with an info line) when gmes is absent or
    // the kingdom is non-fungal, so unprovisioned environments are unaffected.
    let gm_mode = resolve_genemark_mode(
        config.genemark_mode,
        config.use_genemark,
        config.genemark_hints.as_deref(),
    )
    .or_else(|| {
        if !matches!(config.kingdom, Kingdom::Fungi) {
            None
        } else if genemark::genemark_available() {
            println!("  [+] GeneMark-ES enabled by default for fungi (gmes found on PATH).");
            Some(GeneMarkMode::Es)
        } else {
            println!("  [i] GeneMark-ES not run: GeneMark-ES suite not found on PATH (optional).");
            None
        }
    });
    if let Some(mode) = gm_mode {
        let gm_dir = config.out_dir.join("genemark");
        let is_fungus = matches!(config.kingdom, Kingdom::Fungi);
        let pb = progress::spinner(format!("Running {}…", mode.long_name()));
        let run_result = match mode {
            GeneMarkMode::Es => {
                genemark::run_genemark_es(&config.masked_fasta, &gm_dir, is_fungus, config.threads)
            }
            GeneMarkMode::Et => match config.genemark_hints.as_ref() {
                Some(hints) if hints.exists() => genemark::run_genemark_et(
                    &config.masked_fasta,
                    hints,
                    &gm_dir,
                    is_fungus,
                    config.threads,
                ),
                Some(hints) => Err(MycoNoteError::InvalidFormat(format!(
                    "GeneMark-ET hints file not found: {}. \
                     Pass --genemark-hints <intron_hints.gff> or drop to \
                     --genemark-mode es.",
                    hints.display()
                ))),
                None => Err(MycoNoteError::InvalidFormat(
                    "GeneMark-ET requires --genemark-hints <intron_hints.gff>. \
                     Use --genemark-mode es for self-training only."
                        .to_string(),
                )),
            },
            GeneMarkMode::Ep => match config.protein_fasta.as_ref() {
                Some(prot) => genemark::run_genemark_ep(
                    &config.masked_fasta,
                    prot,
                    &gm_dir,
                    is_fungus,
                    config.threads,
                ),
                None => Err(MycoNoteError::InvalidFormat(
                    "GeneMark-EP+ requires --protein-fasta <proteins.fa> \
                     (e.g. OrthoDB fungi). \
                     Use --genemark-mode es to skip protein evidence."
                        .to_string(),
                )),
            },
            GeneMarkMode::Etp => match (
                config.genemark_hints.as_ref(),
                config.protein_fasta.as_ref(),
            ) {
                (Some(rna), Some(prot)) => genemark::run_genemark_etp(
                    &config.masked_fasta,
                    rna,
                    prot,
                    &gm_dir,
                    is_fungus,
                    config.threads,
                ),
                _ => Err(MycoNoteError::InvalidFormat(
                    "GeneMark-ETP+ requires *both* --genemark-hints <intron_hints.gff> \
                     and --protein-fasta <proteins.fa>. \
                     Use --genemark-mode et or ep to drop one input."
                        .to_string(),
                )),
            },
        };
        match run_result {
            Ok(gm_gff) => {
                progress::finish_spinner(
                    &pb,
                    format!("{} complete → {}", mode.long_name(), gm_gff.display()),
                );
                prediction_inputs.push((gm_gff, "GeneMark", weights.genemark));
            }
            Err(e) => {
                progress::warn_spinner(
                    &pb,
                    format!("{} failed (non-fatal): {}", mode.long_name(), e),
                );
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

    // Consensus false-positive filter. Defaults (support=1, score=None) make
    // this a no-op that leaves the output byte-identical to prior releases.
    let consensus_filter = evidence::ConsensusFilter {
        min_predictor_support: config.consensus_min_support,
        min_consensus_score: config.consensus_min_score,
    };
    let gene_count = evidence::merge_predictions_filtered(
        &inputs_ref,
        &consensus_gff,
        &config.locus_prefix,
        &consensus_filter,
    )?;

    println!("  ✓  Consensus GFF3 → {}", consensus_gff.display());
    println!("     {} genes called", gene_count);

    // ── 4. Write a short summary ──────────────────────────────────────────────
    let summary_path = config.out_dir.join("predict_summary.txt");
    write_summary(&summary_path, config, gene_count, &prediction_inputs)?;

    Ok((consensus_gff, gene_count))
}

// ─────────────────────────────────────────────────────────────────────────────
// Gene counter (used on the BRAKER short-circuit path where EVM doesn't run)
// ─────────────────────────────────────────────────────────────────────────────

/// Count `gene` features in a GFF3 file.  Used after a BRAKER run where we
/// don't go through the EVM consensus stage (which already returned a count).
fn count_genes_in_gff3(gff: &Path) -> Result<usize> {
    use std::io::{BufRead, BufReader};
    let file = std::fs::File::open(gff).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut n = 0usize;
    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() >= 3 && cols[2] == "gene" {
            n += 1;
        }
    }
    Ok(n)
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
        extrinsic_cfg: None,
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

    let report = train::run_training(&train_config)?;

    // A run that technically "succeeds" can still leave an unusable model — no
    // registered species directory, or one that recognises zero genes on the
    // held-out test set. Shipping such a model collapses prediction on divergent
    // fungi (the exact failure mode of the stock S288C species), so treat it as
    // a training failure and let `run_prediction` fall back to the stock species.
    if !trained_model_usable(&report) {
        return Err(MycoNoteError::InvalidFormat(format!(
            "trained species '{}' is unusable (species_path={:?}, gene_sensitivity={:?}); \
             falling back to stock/default species",
            report.species_name, report.species_path, report.gene_sensitivity
        )));
    }

    Ok(species_name)
}

/// Decide whether a freshly trained Augustus model is fit to predict with.
///
/// A model is usable only when Augustus actually registered a species directory
/// for it AND, where a held-out accuracy estimate exists, gene-level sensitivity
/// is above zero. A model that identifies no genes on its own test set is worse
/// than the stock species, so the caller falls back rather than shipping it.
/// When no evaluation is available (e.g. the test split was too small to score)
/// we trust etraining and keep the model.
fn trained_model_usable(report: &train::TrainReport) -> bool {
    if report.species_path.is_none() {
        return false;
    }
    match report.gene_sensitivity {
        Some(sens) => sens > 0.0,
        None => true,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Protein → Augustus hints (Feature 3, BRAKER-style)
// ─────────────────────────────────────────────────────────────────────────────

/// Run miniprot (or exonerate) on a protein FASTA and return the Augustus hints
/// file it produces. Robust by design: a missing protein file or a missing
/// aligner logs a warning and returns `Ok(None)` so prediction proceeds
/// ab-initio. Only a genuine I/O error (directory creation) propagates.
fn build_protein_hints_miniprot(config: &PredictConfig, prot_fa: &Path) -> Result<Option<PathBuf>> {
    if !prot_fa.exists() {
        eprintln!(
            "  ⚠  --protein-hints file not found: {} — continuing ab-initio (non-fatal).",
            prot_fa.display()
        );
        return Ok(None);
    }

    let prot_cfg = protein_evidence::ProteinEvidenceConfig {
        proteins: prot_fa.to_path_buf(),
        genome: config.masked_fasta.clone(),
        out_dir: config.out_dir.join("protein_hints"),
        threads: config.threads,
        max_intron: config.max_intron,
        ..protein_evidence::ProteinEvidenceConfig::default()
    };

    println!(
        "  Building Augustus protein hints via miniprot ({})…",
        prot_fa.display()
    );
    match protein_evidence::generate_protein_evidence(&prot_cfg) {
        Ok(result) => {
            println!(
                "      {} proteins aligned to {} loci ({}) → {}",
                result.n_aligned,
                result.n_loci,
                result.tool,
                result.hints_gff.display()
            );
            Ok(Some(result.hints_gff))
        }
        Err(e) => {
            // miniprot/exonerate absent or alignment failed — not fatal.
            eprintln!(
                "  ⚠  protein-hints alignment unavailable: {} — continuing ab-initio (non-fatal).",
                e
            );
            Ok(None)
        }
    }
}

/// Collapse zero/one/many Augustus hint files into a single `Option<PathBuf>`.
/// Zero files → `None` (ab-initio). One file → that file as-is. Many → their
/// concatenation written to `<out_dir>/combined_hints.gff` so Augustus reads a
/// single `--hintsfile`. Concatenation keeps every source's lines (no silent
/// drop of either protein-BLAST or miniprot hints when both flags are given).
fn combine_hint_files(files: &[PathBuf], out_dir: &Path) -> Result<Option<PathBuf>> {
    match files.len() {
        0 => Ok(None),
        1 => Ok(Some(files[0].clone())),
        _ => {
            use std::io::Write;
            let combined = out_dir.join("combined_hints.gff");
            let mut out = std::fs::File::create(&combined).map_err(MycoNoteError::Io)?;
            for f in files {
                let body = std::fs::read_to_string(f).map_err(MycoNoteError::Io)?;
                out.write_all(body.as_bytes()).map_err(MycoNoteError::Io)?;
                if !body.ends_with('\n') {
                    writeln!(out).map_err(MycoNoteError::Io)?;
                }
            }
            Ok(Some(combined))
        }
    }
}

/// Choose the Augustus `--extrinsicCfgFile` basename. When protein hints are in
/// play we need a config whose `[SOURCES]` block lists `P`; Augustus ships
/// `extrinsic.M.RM.E.W.P.cfg` for exactly this. Otherwise `None` lets
/// `AugustusConfig` fall back to its historical default, leaving the existing
/// EST/BLAST hints path untouched.
fn augustus_extrinsic_cfg(use_protein: bool) -> Option<String> {
    if use_protein {
        Some("extrinsic.M.RM.E.W.P.cfg".to_string())
    } else {
        None
    }
}

#[cfg(test)]
mod predict_feature_tests {
    use super::*;

    fn report(species_path: Option<PathBuf>, sens: Option<f64>) -> train::TrainReport {
        train::TrainReport {
            species_name: "test_trained".to_string(),
            n_training_genes: 200,
            n_test_genes: 20,
            gene_sensitivity: sens,
            gene_specificity: sens, // mirror; not inspected by trained_model_usable
            species_path,
        }
    }

    // ── Feature 2: self-trained Augustus usability + fallback ─────────────────

    #[test]
    fn trained_model_usable_requires_species_path() {
        // No registered species dir → unusable, caller must fall back.
        assert!(!trained_model_usable(&report(None, Some(80.0))));
    }

    #[test]
    fn trained_model_unusable_at_zero_sensitivity() {
        let p = Some(PathBuf::from("/tmp/augustus/species/test_trained"));
        assert!(!trained_model_usable(&report(p, Some(0.0))));
    }

    #[test]
    fn trained_model_usable_with_positive_sensitivity() {
        let p = Some(PathBuf::from("/tmp/augustus/species/test_trained"));
        assert!(trained_model_usable(&report(p, Some(62.5))));
    }

    #[test]
    fn trained_model_usable_when_no_eval_available() {
        // Test set too small to score → trust etraining, keep the model.
        let p = Some(PathBuf::from("/tmp/augustus/species/test_trained"));
        assert!(trained_model_usable(&report(p, None)));
    }

    // ── Feature 3: protein → Augustus hints plumbing ──────────────────────────

    #[test]
    fn extrinsic_cfg_selected_only_with_protein_hints() {
        assert_eq!(augustus_extrinsic_cfg(false), None);
        assert_eq!(
            augustus_extrinsic_cfg(true).as_deref(),
            Some("extrinsic.M.RM.E.W.P.cfg")
        );
    }

    #[test]
    fn combine_hint_files_zero_one_many() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path();

        // Zero → None (ab-initio).
        assert!(combine_hint_files(&[], out).unwrap().is_none());

        // One → passed through unchanged (no copy).
        let a = out.join("a.gff");
        std::fs::write(&a, "chr1\tP\tCDSpart\t1\t9\t.\t+\t.\tsrc=P\n").unwrap();
        let one = combine_hint_files(&[a.clone()], out).unwrap().unwrap();
        assert_eq!(one, a);

        // Many → concatenated into combined_hints.gff with every line kept.
        let b = out.join("b.gff");
        std::fs::write(&b, "chr1\tP\tintron\t20\t40\t.\t+\t.\tsrc=P\n").unwrap();
        let many = combine_hint_files(&[a, b], out).unwrap().unwrap();
        assert_eq!(many.file_name().unwrap(), "combined_hints.gff");
        let body = std::fs::read_to_string(&many).unwrap();
        assert!(body.contains("CDSpart"), "first source kept");
        assert!(body.contains("intron"), "second source kept");
    }

    #[test]
    fn protein_hints_missing_file_is_nonfatal() {
        // A protein FASTA path that does not exist must not fail the run; it
        // logs and returns Ok(None) so Augustus proceeds ab-initio.
        let dir = tempfile::tempdir().unwrap();
        let cfg = PredictConfig {
            masked_fasta: dir.path().join("genome.fa"),
            out_dir: dir.path().to_path_buf(),
            ..PredictConfig::default()
        };
        let missing = dir.path().join("does_not_exist.faa");
        let res = build_protein_hints_miniprot(&cfg, &missing).unwrap();
        assert!(res.is_none(), "missing protein file → ab-initio, non-fatal");
    }
}
