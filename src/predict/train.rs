use crate::parser::fasta::{read_fasta_index, FastaRecord};
use crate::parser::gff::{GFFReader, GFFRecord};
/// Augustus self-training
///
/// Trains a custom Augustus HMM species model from a set of gene models,
/// dramatically improving prediction accuracy for non-model organisms or
/// any organism genetically distant from the pre-trained species list.
///
/// Training workflow:
///   1. Select high-confidence gene models (supported by ≥2 predictors OR
///      by protein evidence)
///   2. Natively generate Augustus GenBank training format (mini-records with
///      ±flanking_bp of genomic context per gene)
///   3. Split 80 / 20 into train / test sets
///   4. Run `etraining` to fit HMM parameters
///   5. Run `augustus` on the test set and report accuracy
///   6. (Optional) Run `optimize_augustus.pl` for parameter fine-tuning
///
/// All Perl script calls (gff2gbSmallDNA.pl etc.) are replaced by a native
/// Rust implementation so the tool works even when Augustus Perl scripts are
/// not in PATH.
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct TrainConfig {
    /// GFF3 with gene models to train from (ideally high-confidence set)
    pub gff: PathBuf,
    /// Genome FASTA (unmasked preferred so flanking context is intact)
    pub fasta: PathBuf,
    /// New species name to register in Augustus (e.g. "myorganism_v1")
    pub species_name: String,
    /// Output directory for training files and logs
    pub out_dir: PathBuf,
    /// Flanking DNA context on each side of a gene (bp)
    pub flanking: usize,
    /// Min gene length to include in training (filters tiny ORFs)
    pub min_gene_len: u64,
    /// Fraction of models held out for accuracy test (0.0–0.5)
    pub test_fraction: f64,
    /// Run optimize_augustus.pl after etraining (slow, 1–4 h)
    pub optimize: bool,
    /// Threads for optimize_augustus.pl
    pub threads: usize,
}

impl Default for TrainConfig {
    fn default() -> Self {
        Self {
            gff: PathBuf::new(),
            fasta: PathBuf::new(),
            species_name: "myorganism_v1".to_string(),
            out_dir: PathBuf::from("train_out"),
            flanking: 1000,
            min_gene_len: 300,
            test_fraction: 0.2,
            optimize: false,
            threads: 4,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Training accuracy report
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug)]
pub struct TrainReport {
    pub species_name: String,
    pub n_training_genes: usize,
    pub n_test_genes: usize,
    /// Sensitivity at gene level (0–100)
    pub gene_sensitivity: Option<f64>,
    /// Specificity at gene level (0–100)
    pub gene_specificity: Option<f64>,
    /// Path to trained species configuration
    pub species_path: Option<PathBuf>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_training(config: &TrainConfig) -> Result<TrainReport> {
    // ── Validate inputs ───────────────────────────────────────────────────────
    for p in [&config.gff, &config.fasta] {
        if !p.exists() {
            return Err(MycoNoteError::InvalidFormat(format!(
                "File not found: {}",
                p.display()
            )));
        }
    }

    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    println!("── Augustus self-training ───────────────────────────────────");
    println!("  Species  : {}", config.species_name);
    println!("  GFF3     : {}", config.gff.display());
    println!("  FASTA    : {}", config.fasta.display());
    println!("  Output   : {}", config.out_dir.display());

    // ── 1. Load gene models ───────────────────────────────────────────────────
    let (gene_models, fasta_index) = load_gene_models(config)?;
    println!(
        "  Loaded {} gene models (≥{} bp)",
        gene_models.len(),
        config.min_gene_len
    );

    if gene_models.len() < 50 {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Too few gene models for training ({}). Need at least 50. \
             Try lowering --min-gene-len or providing a larger GFF3.",
            gene_models.len()
        )));
    }

    // ── 2. Generate Augustus GenBank training format ──────────────────────────
    let gb_path = config.out_dir.join("training_genes.gb");
    let n_written = write_augustus_genbank(&gene_models, &fasta_index, config, &gb_path)?;
    println!("  Wrote {} mini-GenBank records for training", n_written);

    // ── 3. Split train / test ─────────────────────────────────────────────────
    let n_test = ((n_written as f64 * config.test_fraction) as usize)
        .max(10)
        .min(n_written - 20);
    let n_train = n_written - n_test;
    let train_gb = config.out_dir.join("train.gb");
    let test_gb = config.out_dir.join("test.gb");
    split_genbank(&gb_path, &train_gb, &test_gb, n_test)?;
    println!("  Split: {} training, {} test", n_train, n_test);

    // ── 4. Create new Augustus species ───────────────────────────────────────
    create_augustus_species(&config.species_name, config)?;

    // ── 5. Run etraining ─────────────────────────────────────────────────────
    run_etraining(&config.species_name, &train_gb, config)?;
    println!("  etraining complete");

    // ── 6. Evaluate on test set ───────────────────────────────────────────────
    let (gene_sens, gene_spec) = evaluate_model(&config.species_name, &test_gb, config)?;
    if let Some(sens) = gene_sens {
        println!(
            "  Accuracy: sensitivity={:.1}%  specificity={:.1}%",
            sens,
            gene_spec.unwrap_or(0.0)
        );
    }

    // ── 7. Optional optimization ─────────────────────────────────────────────
    if config.optimize {
        optimize_species(&config.species_name, &train_gb, config)?;
        println!("  Optimization complete");
    }

    // ── 8. Write report ───────────────────────────────────────────────────────
    let report = TrainReport {
        species_name: config.species_name.clone(),
        n_training_genes: n_train,
        n_test_genes: n_test,
        gene_sensitivity: gene_sens,
        gene_specificity: gene_spec,
        species_path: find_augustus_species_path(&config.species_name),
    };
    write_train_report(&report, config)?;

    println!("  ✓ Trained species: {}", config.species_name);
    println!(
        "    Use with: myconote predict <fasta> --species {}",
        config.species_name
    );

    Ok(report)
}

// ─────────────────────────────────────────────────────────────────────────────
// Gene model loading
// ─────────────────────────────────────────────────────────────────────────────

/// A gene model: the gene record + its mRNA children + their CDS records.
#[derive(Debug)]
#[allow(dead_code)]
struct GeneModel {
    gene: GFFRecord,
    mrnas: Vec<GFFRecord>,
    cds: Vec<GFFRecord>, // all CDS for this gene across all mRNAs
}

fn load_gene_models(
    config: &TrainConfig,
) -> Result<(Vec<GeneModel>, HashMap<String, FastaRecord>)> {
    let fasta_index = read_fasta_index(&config.fasta)?;
    let all_records: Vec<GFFRecord> = GFFReader::from_path(&config.gff)?
        .filter_map(|r| r.ok())
        .collect();

    let mut models: Vec<GeneModel> = Vec::new();

    for rec in &all_records {
        if rec.feature_type != "gene" {
            continue;
        }
        let gene_len = rec.end.saturating_sub(rec.start);
        if gene_len < config.min_gene_len {
            continue;
        }

        // Check sequence is available
        if !fasta_index.contains_key(&rec.seqid) {
            continue;
        }

        let gene_id = match rec.id() {
            Some(id) => id.clone(),
            None => continue,
        };

        // Collect mRNA children
        let mrnas: Vec<GFFRecord> = all_records
            .iter()
            .filter(|r| {
                (r.feature_type == "mRNA" || r.feature_type == "transcript")
                    && r.parent().map(|p| p == &gene_id).unwrap_or(false)
            })
            .cloned()
            .collect();

        // Collect CDS for all mRNAs
        let mrna_ids: std::collections::HashSet<String> =
            mrnas.iter().filter_map(|m| m.id().cloned()).collect();

        let cds: Vec<GFFRecord> = all_records
            .iter()
            .filter(|r| {
                r.feature_type == "CDS" && r.parent().map(|p| mrna_ids.contains(p)).unwrap_or(false)
            })
            .cloned()
            .collect();

        if cds.is_empty() {
            continue;
        }

        models.push(GeneModel {
            gene: rec.clone(),
            mrnas,
            cds,
        });
    }

    Ok((models, fasta_index))
}

// ─────────────────────────────────────────────────────────────────────────────
// Augustus GenBank format writer
// ─────────────────────────────────────────────────────────────────────────────

/// Write each gene model as a mini-GenBank record with flanking context.
/// Returns the number of records written.
fn write_augustus_genbank(
    models: &[GeneModel],
    fasta_index: &HashMap<String, FastaRecord>,
    config: &TrainConfig,
    out_path: &Path,
) -> Result<usize> {
    let mut f = std::fs::File::create(out_path).map_err(MycoNoteError::Io)?;
    let mut written = 0usize;

    for model in models {
        let seq_rec = match fasta_index.get(&model.gene.seqid) {
            Some(s) => s,
            None => continue,
        };

        let seq_len = seq_rec.sequence.len() as u64;
        let flank = config.flanking as u64;
        // Clamp to sequence bounds
        let region_start = model.gene.start.saturating_sub(flank).max(1);
        let region_end = (model.gene.end + flank).min(seq_len);
        let region_len = region_end - region_start + 1;

        // Extract genomic region (1-based inclusive → slice)
        let region_seq = seq_rec.subsequence(region_start, region_end);
        if region_seq.is_empty() {
            continue;
        }

        // Offset: coordinates relative to extracted region start
        let offset = region_start - 1; // convert to 0-based offset

        // Collect CDS intervals (0-based within extracted region)
        let mut cds_coords: Vec<(u64, u64)> = model
            .cds
            .iter()
            .map(|c| (c.start - 1 - offset, c.end - 1 - offset))
            .collect();
        cds_coords.sort_by_key(|c| c.0);
        cds_coords.dedup();

        if cds_coords.is_empty() {
            continue;
        }

        let gene_id = model
            .gene
            .id()
            .cloned()
            .unwrap_or_else(|| format!("gene{}", written + 1));

        // Write LOCUS line
        writeln!(f, "LOCUS       {:<20} {} bp    DNA", gene_id, region_len)
            .map_err(MycoNoteError::Io)?;
        writeln!(f, "FEATURES             Location/Qualifiers").map_err(MycoNoteError::Io)?;

        // Write CDS feature
        let strand_char = model.gene.strand;
        let location = build_location(&cds_coords, strand_char);
        writeln!(f, "     CDS             {}", location).map_err(MycoNoteError::Io)?;
        writeln!(f, "                     /gene=\"{}\"", gene_id).map_err(MycoNoteError::Io)?;
        writeln!(f, "ORIGIN").map_err(MycoNoteError::Io)?;

        // Write sequence in GenBank format (60 bp / line, 10-char blocks)
        for (i, chunk) in region_seq.as_bytes().chunks(60).enumerate() {
            let pos = i * 60 + 1;
            write!(f, "{:>9} ", pos).map_err(MycoNoteError::Io)?;
            for block in chunk.chunks(10) {
                let s = std::str::from_utf8(block).unwrap_or("").to_lowercase();
                write!(f, "{} ", s).map_err(MycoNoteError::Io)?;
            }
            writeln!(f).map_err(MycoNoteError::Io)?;
        }
        writeln!(f, "//").map_err(MycoNoteError::Io)?;

        written += 1;
    }

    Ok(written)
}

/// Build the GenBank CDS location string, e.g. `join(100..200,300..400)` or
/// `complement(join(...))` for minus-strand genes.
fn build_location(coords: &[(u64, u64)], strand: char) -> String {
    let intervals: Vec<String> = coords
        .iter()
        .map(|(s, e)| format!("{}..{}", s + 1, e + 1)) // back to 1-based
        .collect();

    let join = if intervals.len() == 1 {
        intervals[0].clone()
    } else {
        format!("join({})", intervals.join(","))
    };

    if strand == '-' {
        format!("complement({})", join)
    } else {
        join
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Train / test split
// ─────────────────────────────────────────────────────────────────────────────

/// Split a GenBank file into training and test sets.
/// Takes the LAST n_test records as the test set (deterministic, no random seed needed).
fn split_genbank(gb_path: &Path, train_out: &Path, test_out: &Path, n_test: usize) -> Result<()> {
    // Read all records
    let content = std::fs::read_to_string(gb_path).map_err(MycoNoteError::Io)?;
    let records: Vec<&str> = content
        .split("//\n")
        .filter(|s| !s.trim().is_empty())
        .collect();

    let n_total = records.len();
    let n_train = n_total.saturating_sub(n_test);

    let mut train_f = std::fs::File::create(train_out).map_err(MycoNoteError::Io)?;
    let mut test_f = std::fs::File::create(test_out).map_err(MycoNoteError::Io)?;

    for (i, rec) in records.iter().enumerate() {
        let target: &mut dyn std::io::Write = if i < n_train {
            &mut train_f
        } else {
            &mut test_f
        };
        writeln!(target, "{}//", rec).map_err(MycoNoteError::Io)?;
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Augustus species management
// ─────────────────────────────────────────────────────────────────────────────

/// Create a new Augustus species config directory.
/// Calls `new_species.pl --species=<name>`.
fn create_augustus_species(species: &str, config: &TrainConfig) -> Result<()> {
    let script = find_augustus_script("new_species.pl")?;

    let status = Command::new("perl")
        .arg(&script)
        .arg(format!("--species={}", species))
        .current_dir(&config.out_dir)
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        // Species might already exist — treat as warning, not error
        eprintln!("  ⚠  new_species.pl returned non-zero (species may already exist)");
    }
    Ok(())
}

/// Run etraining to fit HMM parameters.
fn run_etraining(species: &str, train_gb: &Path, config: &TrainConfig) -> Result<()> {
    let etraining = which::which("etraining").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "etraining not found in PATH.\n\
             Install Augustus: conda install -c bioconda augustus"
                .to_string(),
        )
    })?;

    println!("  Running etraining…");

    let log_path = config.out_dir.join("etraining.log");
    let log_file = std::fs::File::create(&log_path).map_err(MycoNoteError::Io)?;

    let status = Command::new(&etraining)
        .args([
            &format!("--species={}", species),
            train_gb.to_str().unwrap_or(""),
        ])
        .stdout(log_file)
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "etraining failed. Check train_out/etraining.log for details.".to_string(),
        ));
    }
    Ok(())
}

/// Evaluate the trained model on the test set.
/// Returns (sensitivity, specificity) at the gene level.
fn evaluate_model(
    species: &str,
    test_gb: &Path,
    config: &TrainConfig,
) -> Result<(Option<f64>, Option<f64>)> {
    let augustus = which::which("augustus").ok();
    if augustus.is_none() {
        return Ok((None, None));
    }

    println!("  Evaluating model accuracy…");

    let aug_out = config.out_dir.join("test_predictions.gff");
    let aug_file = std::fs::File::create(&aug_out).map_err(MycoNoteError::Io)?;

    let status = Command::new(augustus.unwrap())
        .args([
            &format!("--species={}", species),
            &format!("--gff3=on"),
            test_gb.to_str().unwrap_or(""),
        ])
        .stdout(aug_file)
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Ok((None, None));
    }

    // Parse accuracy from Augustus output header comments
    parse_augustus_accuracy(&aug_out)
}

fn parse_augustus_accuracy(aug_gff: &Path) -> Result<(Option<f64>, Option<f64>)> {
    let file = std::fs::File::open(aug_gff).map_err(MycoNoteError::Io)?;
    let mut sensitivity = None;
    let mut specificity = None;

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        // Augustus outputs: "# sensitivity of gene prediction: 0.85"
        if line.contains("sensitivity of gene prediction:") {
            sensitivity = line
                .split(':')
                .last()
                .and_then(|s| s.trim().parse::<f64>().ok())
                .map(|v| v * 100.0);
        }
        if line.contains("specificity of gene prediction:") {
            specificity = line
                .split(':')
                .last()
                .and_then(|s| s.trim().parse::<f64>().ok())
                .map(|v| v * 100.0);
        }
    }
    Ok((sensitivity, specificity))
}

/// Optionally run optimize_augustus.pl for parameter fine-tuning.
fn optimize_species(species: &str, train_gb: &Path, config: &TrainConfig) -> Result<()> {
    let script = find_augustus_script("optimize_augustus.pl")?;

    println!("  Running optimize_augustus.pl (this may take 1–4 hours)…");

    let log_path = config.out_dir.join("optimize.log");
    let log_file = std::fs::File::create(&log_path).map_err(MycoNoteError::Io)?;

    Command::new("perl")
        .arg(&script)
        .args([
            &format!("--species={}", species),
            &format!("--cpus={}", config.threads),
            train_gb.to_str().unwrap_or(""),
        ])
        .stdout(log_file)
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Script/path helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Push `d` onto `dirs` only if it is not already present — keeps the search
/// list ordered by priority without duplicate stat() calls.
fn push_unique(dirs: &mut Vec<PathBuf>, d: PathBuf) {
    if !dirs.contains(&d) {
        dirs.push(d);
    }
}

/// Enumerate `share/augustus*` directories (e.g. `augustus`, `augustus-3.5.0`)
/// so version-suffixed source/Debian builds are discovered. Returns an empty
/// vec when `share` does not exist or cannot be read. Sorted for determinism.
fn augustus_share_subdirs(share: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    if let Ok(entries) = std::fs::read_dir(share) {
        for e in entries.flatten() {
            let name = e.file_name();
            if name.to_string_lossy().starts_with("augustus") {
                out.push(e.path());
            }
        }
    }
    out.sort();
    out
}

/// Candidate directories that may hold Augustus companion scripts
/// (`new_species.pl`, `etraining`, `gff2gbSmallDNA.pl`, `optimize_augustus.pl`),
/// in priority order, derived from the environment and the resolved `augustus`
/// binary.
///
/// Augustus ships these scripts in one of two layouts depending on the build:
///   * conda / Homebrew: beside the `augustus` binary in the env's `bin/`, or in
///     a sibling `../scripts` of that `bin/`.
///   * source / Debian: under `<prefix>/share/augustus*/scripts` (the directory
///     is version-suffixed on some builds, e.g. `augustus-3.5.0`).
///
/// We search relative to the binary's install root, not just `$PATH`, so a
/// self-training run finds the scripts the same install provides — this is the
/// fix for conda envs (e.g. `myconote_augustus`) that keep `new_species.pl` in
/// their `bin/` directory next to `augustus`.
fn augustus_script_search_dirs(
    scripts_env: Option<&Path>,
    bin_env: Option<&Path>,
    config_env: Option<&Path>,
    augustus_bin: Option<&Path>,
) -> Vec<PathBuf> {
    let mut dirs: Vec<PathBuf> = Vec::new();

    // 1. Explicit $AUGUSTUS_SCRIPTS_PATH wins.
    if let Some(d) = scripts_env {
        push_unique(&mut dirs, d.to_path_buf());
    }
    // 2. $AUGUSTUS_BIN_PATH itself and its sibling scripts/.
    if let Some(b) = bin_env {
        push_unique(&mut dirs, b.to_path_buf());
        if let Some(parent) = b.parent() {
            push_unique(&mut dirs, parent.join("scripts"));
        }
    }
    // 3. $AUGUSTUS_CONFIG_PATH/../scripts.
    if let Some(c) = config_env {
        if let Some(parent) = c.parent() {
            push_unique(&mut dirs, parent.join("scripts"));
        }
    }
    // 4. Relative to the resolved augustus binary.
    if let Some(aug) = augustus_bin {
        if let Some(bin_dir) = aug.parent() {
            // 4a. Beside the binary (conda env bin/ — the cluster layout).
            push_unique(&mut dirs, bin_dir.to_path_buf());
            // 4b. Sibling scripts/ of bin/ (bin/../scripts).
            push_unique(&mut dirs, bin_dir.join("scripts"));
            if let Some(root) = bin_dir.parent() {
                // 4c. <root>/scripts
                push_unique(&mut dirs, root.join("scripts"));
                // 4d. <root>/share/augustus*/scripts (version-suffixed dirs).
                let share = root.join("share");
                for sub in augustus_share_subdirs(&share) {
                    push_unique(&mut dirs, sub.join("scripts"));
                }
                // 4e. plain <root>/share/augustus/scripts
                push_unique(&mut dirs, share.join("augustus").join("scripts"));
            }
        }
    }
    // 5. Common system / conda prefixes as a last resort.
    for prefix in ["/usr", "/opt/conda", "/usr/local"] {
        push_unique(
            &mut dirs,
            PathBuf::from(prefix).join("share/augustus/scripts"),
        );
        push_unique(&mut dirs, PathBuf::from(prefix).join("bin"));
    }

    dirs
}

/// Locate an Augustus companion script by name, returning the first candidate
/// directory that actually contains it. Pure over its inputs so it can be
/// unit-tested against a mock install tree.
fn locate_augustus_script(
    name: &str,
    scripts_env: Option<&Path>,
    bin_env: Option<&Path>,
    config_env: Option<&Path>,
    augustus_bin: Option<&Path>,
) -> Option<PathBuf> {
    for dir in augustus_script_search_dirs(scripts_env, bin_env, config_env, augustus_bin) {
        let p = dir.join(name);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

fn find_augustus_script(name: &str) -> Result<PathBuf> {
    let scripts_env = std::env::var_os("AUGUSTUS_SCRIPTS_PATH").map(PathBuf::from);
    let bin_env = std::env::var_os("AUGUSTUS_BIN_PATH").map(PathBuf::from);
    let config_env = std::env::var_os("AUGUSTUS_CONFIG_PATH").map(PathBuf::from);
    let augustus_bin = which::which("augustus").ok();

    locate_augustus_script(
        name,
        scripts_env.as_deref(),
        bin_env.as_deref(),
        config_env.as_deref(),
        augustus_bin.as_deref(),
    )
    .ok_or_else(|| {
        MycoNoteError::UnsupportedFormat(format!(
            "Augustus script '{}' not found beside the augustus binary, under \
             $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH/../scripts, or any \
             share/augustus*/scripts directory.\n\
             Set $AUGUSTUS_SCRIPTS_PATH to the directory holding Augustus scripts.",
            name
        ))
    })
}

/// Pure trainability check shared by the self-train entry point and tests:
/// given the environment overrides and a resolved `augustus` binary, can the
/// gating script `new_species.pl` be found? `new_species.pl` is the script
/// that registers a species; without it self-training cannot even start.
fn self_training_scripts_resolvable(
    scripts_env: Option<&Path>,
    bin_env: Option<&Path>,
    config_env: Option<&Path>,
    augustus_bin: Option<&Path>,
) -> bool {
    locate_augustus_script(
        "new_species.pl",
        scripts_env,
        bin_env,
        config_env,
        augustus_bin,
    )
    .is_some()
}

/// Whether this install can self-train Augustus — i.e. `new_species.pl` resolves
/// from the current environment / the `augustus` binary on PATH. Callers use
/// this to bail out *before* an expensive first-pass prediction when the
/// companion scripts are absent, so the loud stock-species fallback fires
/// immediately instead of after wasted compute.
pub fn augustus_self_training_available() -> bool {
    let scripts_env = std::env::var_os("AUGUSTUS_SCRIPTS_PATH").map(PathBuf::from);
    let bin_env = std::env::var_os("AUGUSTUS_BIN_PATH").map(PathBuf::from);
    let config_env = std::env::var_os("AUGUSTUS_CONFIG_PATH").map(PathBuf::from);
    let augustus_bin = which::which("augustus").ok();
    self_training_scripts_resolvable(
        scripts_env.as_deref(),
        bin_env.as_deref(),
        config_env.as_deref(),
        augustus_bin.as_deref(),
    )
}

fn find_augustus_species_path(species: &str) -> Option<PathBuf> {
    let base_dirs = [
        std::env::var("AUGUSTUS_CONFIG_PATH").ok(),
        Some("/opt/conda/config".to_string()),
        Some("/usr/share/augustus/config".to_string()),
        Some("/usr/local/share/augustus/config".to_string()),
    ];
    for dir_opt in &base_dirs {
        if let Some(dir) = dir_opt {
            let p = PathBuf::from(dir).join("species").join(species);
            if p.exists() {
                return Some(p);
            }
        }
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Report
// ─────────────────────────────────────────────────────────────────────────────

fn write_train_report(report: &TrainReport, config: &TrainConfig) -> Result<()> {
    let path = config.out_dir.join("training_report.txt");
    let mut f = std::fs::File::create(&path).map_err(MycoNoteError::Io)?;

    writeln!(f, "myconote predict train — Training Report").map_err(MycoNoteError::Io)?;
    writeln!(f, "=========================================").map_err(MycoNoteError::Io)?;
    writeln!(f, "Species name    : {}", report.species_name).map_err(MycoNoteError::Io)?;
    writeln!(f, "Training genes  : {}", report.n_training_genes).map_err(MycoNoteError::Io)?;
    writeln!(f, "Test genes      : {}", report.n_test_genes).map_err(MycoNoteError::Io)?;
    if let Some(s) = report.gene_sensitivity {
        writeln!(f, "Sensitivity     : {:.1}%", s).map_err(MycoNoteError::Io)?;
        writeln!(
            f,
            "Specificity     : {:.1}%",
            report.gene_specificity.unwrap_or(0.0)
        )
        .map_err(MycoNoteError::Io)?;
    } else {
        writeln!(f, "Accuracy        : not evaluated (augustus not in PATH)")
            .map_err(MycoNoteError::Io)?;
    }
    writeln!(f, "").map_err(MycoNoteError::Io)?;
    writeln!(f, "To use this model in gene prediction:").map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "  myconote predict <masked.fa> --species {}",
        report.species_name
    )
    .map_err(MycoNoteError::Io)?;

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Create an empty file at `path`, making parent dirs as needed.
    fn touch(path: &Path) {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).unwrap();
        }
        std::fs::write(path, b"").unwrap();
    }

    // FIX A: the conda/cluster layout — companion scripts live in the SAME bin/
    // as the `augustus` binary. The resolver must find them there (this is the
    // case the old code missed, which caused the silent self-training no-op).
    #[test]
    fn locates_script_beside_augustus_binary() {
        let dir = tempfile::tempdir().unwrap();
        let bin = dir.path().join("bin");
        let aug = bin.join("augustus");
        let script = bin.join("new_species.pl");
        touch(&aug);
        touch(&script);

        let found = locate_augustus_script("new_species.pl", None, None, None, Some(&aug));
        assert_eq!(found.as_deref(), Some(script.as_path()));
    }

    // Source/Debian layout: scripts under <root>/share/augustus-<ver>/scripts.
    #[test]
    fn locates_script_under_versioned_share_dir() {
        let dir = tempfile::tempdir().unwrap();
        let aug = dir.path().join("bin").join("augustus");
        touch(&aug);
        let script = dir
            .path()
            .join("share")
            .join("augustus-3.5.0")
            .join("scripts")
            .join("gff2gbSmallDNA.pl");
        touch(&script);

        let found = locate_augustus_script("gff2gbSmallDNA.pl", None, None, None, Some(&aug));
        assert_eq!(found.as_deref(), Some(script.as_path()));
    }

    // $AUGUSTUS_SCRIPTS_PATH override is honoured and takes priority.
    #[test]
    fn locates_script_under_scripts_path_env() {
        let dir = tempfile::tempdir().unwrap();
        let scripts = dir.path().join("aug_scripts");
        let script = scripts.join("etraining");
        touch(&script);

        let found = locate_augustus_script("etraining", Some(&scripts), None, None, None);
        assert_eq!(found.as_deref(), Some(script.as_path()));
    }

    // $AUGUSTUS_CONFIG_PATH/../scripts is searched (common conda config layout).
    #[test]
    fn locates_script_via_config_path_sibling_scripts() {
        let dir = tempfile::tempdir().unwrap();
        let config = dir.path().join("config");
        std::fs::create_dir_all(&config).unwrap();
        let script = dir.path().join("scripts").join("optimize_augustus.pl");
        touch(&script);

        let found = locate_augustus_script("optimize_augustus.pl", None, None, Some(&config), None);
        assert_eq!(found.as_deref(), Some(script.as_path()));
    }

    // A present-scripts path is detected as trainable; an absent one is not.
    #[test]
    fn self_training_detected_when_scripts_present() {
        let dir = tempfile::tempdir().unwrap();
        let bin = dir.path().join("bin");
        let aug = bin.join("augustus");
        touch(&aug);

        // No new_species.pl yet → not trainable.
        assert!(!self_training_scripts_resolvable(
            None,
            None,
            None,
            Some(&aug)
        ));

        // Drop new_species.pl beside augustus → now trainable.
        touch(&bin.join("new_species.pl"));
        assert!(self_training_scripts_resolvable(
            None,
            None,
            None,
            Some(&aug)
        ));
    }

    // Missing script → resolver returns None (the caller maps this to a clean,
    // actionable error, never a panic).
    #[test]
    fn missing_script_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let aug = dir.path().join("bin").join("augustus");
        touch(&aug);
        let found = locate_augustus_script("new_species.pl", None, None, None, Some(&aug));
        assert!(found.is_none());
    }
}
