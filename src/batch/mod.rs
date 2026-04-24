/// Batch genome annotation — run the full pipeline across multiple genomes.
///
/// Accepts a directory of FASTA files or a sample sheet (TSV), runs
/// sort → mask → predict → annotate → submit per genome, tracks progress
/// with a TTY-aware dashboard, and optionally generates HTCondor submit files.
pub mod condor;
pub mod dashboard;
pub mod state;

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::utils::error::{MycoNoteError, Result};
use dashboard::Dashboard;
use state::{BatchSettings, BatchState};

/// Default pipeline stages in order.
pub const ALL_STAGES: &[&str] = &["sort", "mask", "predict", "annotate", "submit"];

/// A single genome entry parsed from a sample sheet or directory scan.
#[derive(Debug, Clone)]
pub struct GenomeEntry {
    /// Display name (stem of FASTA filename).
    pub name: String,
    /// Path to the input FASTA.
    pub fasta: PathBuf,
    /// Kingdom for prediction (default: fungi).
    pub kingdom: String,
    /// Augustus species override.
    pub species: Option<String>,
    /// Translation table (default: 1 = standard).
    pub genetic_code: u8,
    /// Locus tag prefix.
    pub locus_prefix: String,
}

/// Batch run configuration.
#[derive(Debug, Clone)]
pub struct BatchConfig {
    /// Input: either a directory of FASTAs or a sample sheet TSV.
    pub input: PathBuf,
    /// Output batch directory.
    pub output_dir: PathBuf,
    /// Which stages to run.
    pub stages: Vec<String>,
    /// Kingdom default for genomes without a sample sheet entry.
    pub kingdom: String,
    /// Threads per genome.
    pub threads: usize,
    /// Max genomes to run in parallel (local mode).
    pub max_parallel: usize,
    /// Resume from a previous run.
    pub resume: Option<PathBuf>,
    /// Generate HTCondor submit files instead of running locally.
    pub condor: bool,
    /// HTCondor CPUs per job.
    pub condor_cpus: usize,
    /// HTCondor memory per job.
    pub condor_mem: String,
    /// HTCondor disk per job.
    pub condor_disk: String,
    /// HTCondor accounting group.
    pub condor_queue: Option<String>,
    /// Extra HTCondor submit directives file.
    pub condor_extra: Option<PathBuf>,
    /// Minimum contig length for sort stage.
    pub min_length: usize,
    /// Masking engine.
    pub mask_engine: String,
    /// Genetic code default.
    pub genetic_code: u8,
    /// Locus prefix default.
    pub locus_prefix: String,
}

impl Default for BatchConfig {
    fn default() -> Self {
        Self {
            input: PathBuf::new(),
            output_dir: PathBuf::from("batch_out"),
            stages: ALL_STAGES.iter().map(|s| s.to_string()).collect(),
            kingdom: "fungi".to_string(),
            threads: 4,
            max_parallel: 2,
            resume: None,
            condor: false,
            condor_cpus: 8,
            condor_mem: "32G".to_string(),
            condor_disk: "50G".to_string(),
            condor_queue: None,
            condor_extra: None,
            min_length: 500,
            mask_engine: "repeatmodeler".to_string(),
            genetic_code: 1,
            locus_prefix: "GENE".to_string(),
        }
    }
}

/// Entry point for the batch subcommand.
pub fn run_batch(config: &BatchConfig) -> Result<()> {
    // ── 1. Discover genomes ──
    let genomes = discover_genomes(config)?;

    if genomes.is_empty() {
        return Err(MycoNoteError::BatchError(
            "no FASTA files found. Provide a directory of .fa/.fasta files or a sample sheet TSV."
                .to_string(),
        ));
    }

    println!("  Found {} genome(s) to annotate", genomes.len());

    // ── 2. HTCondor mode — generate submit files and exit ──
    if config.condor {
        return run_condor_mode(config, &genomes);
    }

    // ── 3. Local mode — run pipeline per genome ──
    run_local_mode(config, &genomes)
}

/// Discover genomes from a directory or sample sheet.
fn discover_genomes(config: &BatchConfig) -> Result<Vec<GenomeEntry>> {
    let input = &config.input;

    if input.is_file() {
        // Treat as sample sheet if it's a .tsv or .txt file
        parse_sample_sheet(input, config)
    } else if input.is_dir() {
        scan_directory(input, config)
    } else {
        Err(MycoNoteError::BatchError(format!(
            "'{}' is not a file or directory",
            input.display()
        )))
    }
}

/// Parse a TSV sample sheet.
///
/// Format (header required):
///   name\tfasta\tkingdom\tspecies\tgenetic_code\tlocus_prefix
///
/// Only `name` and `fasta` are required; others use defaults.
fn parse_sample_sheet(path: &Path, config: &BatchConfig) -> Result<Vec<GenomeEntry>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| MycoNoteError::BatchError(format!("cannot read sample sheet: {}", e)))?;

    let mut entries = Vec::new();
    let mut lines = text.lines();

    // Skip header
    let header = lines.next().unwrap_or("");
    let cols: Vec<&str> = header.split('\t').collect();

    // Find column indices
    let idx = |name: &str| {
        cols.iter()
            .position(|c| c.trim().eq_ignore_ascii_case(name))
    };
    let i_name = idx("name");
    let i_fasta = idx("fasta")
        .or_else(|| idx("path"))
        .or_else(|| idx("genome"));
    let i_kingdom = idx("kingdom");
    let i_species = idx("species");
    let i_code = idx("genetic_code").or_else(|| idx("code"));
    let i_prefix = idx("locus_prefix").or_else(|| idx("prefix"));

    let i_fasta = i_fasta.ok_or_else(|| {
        MycoNoteError::BatchError(
            "sample sheet must have a 'fasta' (or 'path' or 'genome') column".to_string(),
        )
    })?;

    for (line_no, line) in lines.enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();

        let fasta = fields.get(i_fasta).map(|s| s.trim()).unwrap_or("");
        if fasta.is_empty() {
            continue;
        }

        let name = i_name
            .and_then(|i| fields.get(i))
            .map(|s| s.trim().to_string())
            .unwrap_or_else(|| {
                Path::new(fasta)
                    .file_stem()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .to_string()
            });

        let kingdom = i_kingdom
            .and_then(|i| fields.get(i))
            .map(|s| s.trim().to_string())
            .unwrap_or_else(|| config.kingdom.clone());

        let species = i_species
            .and_then(|i| fields.get(i))
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty() && s != "auto");

        let genetic_code = i_code
            .and_then(|i| fields.get(i))
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(config.genetic_code);

        let locus_prefix = i_prefix
            .and_then(|i| fields.get(i))
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| config.locus_prefix.clone());

        // Validate FASTA exists
        if !Path::new(fasta).exists() {
            eprintln!(
                "  Warning: line {}: FASTA '{}' not found, skipping",
                line_no + 2,
                fasta
            );
            continue;
        }

        entries.push(GenomeEntry {
            name,
            fasta: PathBuf::from(fasta),
            kingdom,
            species,
            genetic_code,
            locus_prefix,
        });
    }

    Ok(entries)
}

/// Scan a directory for FASTA files.
fn scan_directory(dir: &Path, config: &BatchConfig) -> Result<Vec<GenomeEntry>> {
    let mut entries = Vec::new();

    let read_dir = std::fs::read_dir(dir)
        .map_err(|e| MycoNoteError::BatchError(format!("cannot read directory: {}", e)))?;

    let mut paths: Vec<PathBuf> = read_dir
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && matches!(
                    p.extension().and_then(|e| e.to_str()),
                    Some("fa") | Some("fasta") | Some("fas") | Some("fna")
                )
        })
        .collect();

    paths.sort();

    for path in paths {
        let name = path
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();

        entries.push(GenomeEntry {
            name,
            fasta: path,
            kingdom: config.kingdom.clone(),
            species: None,
            genetic_code: config.genetic_code,
            locus_prefix: config.locus_prefix.clone(),
        });
    }

    Ok(entries)
}

/// Generate HTCondor submit files and exit.
fn run_condor_mode(config: &BatchConfig, genomes: &[GenomeEntry]) -> Result<()> {
    let condor_genomes: Vec<condor::CondorGenome> = genomes
        .iter()
        .map(|g| condor::CondorGenome {
            name: g.name.clone(),
            fasta_path: g.fasta.display().to_string(),
            kingdom: g.kingdom.clone(),
            species: g.species.clone(),
            genetic_code: g.genetic_code,
            locus_prefix: Some(g.locus_prefix.clone()),
        })
        .collect();

    let condor_cfg = condor::CondorConfig {
        cpus: config.condor_cpus,
        memory: config.condor_mem.clone(),
        disk: config.condor_disk.clone(),
        queue: config.condor_queue.clone(),
        extra_file: config.condor_extra.clone(),
        max_idle: None,
    };

    let stage_refs: Vec<&str> = config.stages.iter().map(|s| s.as_str()).collect();

    let sub_path = condor::generate_condor_submit(
        &config.output_dir,
        &condor_genomes,
        &condor_cfg,
        &stage_refs,
        config.threads,
    )?;

    println!("\n  \x1b[32m\u{2713}\x1b[0m HTCondor submit files generated:");
    println!("    Submit file:   {}", sub_path.display());
    println!(
        "    Genome list:   {}",
        config.output_dir.join("condor_genomes.txt").display()
    );
    println!(
        "    Wrapper script: {}",
        config.output_dir.join("run_genome.sh").display()
    );
    println!("\n  To submit:");
    println!("    \x1b[33mcondor_submit {}\x1b[0m", sub_path.display());
    println!("\n  To monitor:");
    println!("    condor_q");
    println!(
        "    tail -f {}/condor_logs/job_0.out",
        config.output_dir.display()
    );

    Ok(())
}

/// Run the pipeline locally for each genome.
fn run_local_mode(config: &BatchConfig, genomes: &[GenomeEntry]) -> Result<()> {
    // Create output directory
    std::fs::create_dir_all(&config.output_dir).map_err(|e| {
        MycoNoteError::BatchError(format!(
            "cannot create {}: {}",
            config.output_dir.display(),
            e
        ))
    })?;

    // Load or create state
    let mut batch_state = if let Some(ref resume_path) = config.resume {
        let state_file = resume_path.join("status.json");
        if state_file.exists() {
            println!("  Resuming from {}", state_file.display());
            BatchState::load(&state_file)?
        } else {
            return Err(MycoNoteError::BatchError(format!(
                "no status.json found in {}",
                resume_path.display()
            )));
        }
    } else {
        let settings = BatchSettings {
            kingdom: config.kingdom.clone(),
            threads: config.threads,
            max_parallel: config.max_parallel,
            stages: config.stages.clone(),
        };
        BatchState::new(&config.output_dir, settings)
    };

    // Initialize dashboard
    let mut dash = Dashboard::new(genomes.len());
    for g in genomes {
        dash.add_genome(&g.name);
    }

    // Run genomes sequentially (parallel via rayon deferred to avoid
    // complexity with external tool subprocesses competing for resources)
    for genome in genomes {
        let gs = batch_state.genome_mut(&genome.name);
        if gs.failed {
            continue; // Skip previously failed genomes on resume
        }

        // Check if all stages already done (resume case)
        let all_done = config.stages.iter().all(|s| gs.is_stage_done(s));
        if all_done {
            dash.genome_done(&genome.name, gs.elapsed_s);
            continue;
        }

        let genome_start = Instant::now();
        let genome_dir = config.output_dir.join(&genome.name);
        std::fs::create_dir_all(&genome_dir).ok();

        let mut current_fasta = genome.fasta.clone();
        let mut current_gff: Option<PathBuf> = None;
        let mut failed = false;

        for stage_name in &config.stages {
            let gs = batch_state.genome_mut(&genome.name);
            if gs.is_stage_done(stage_name) {
                // Recover output path from previous run
                if let Some(p) = gs.outputs.get(stage_name.as_str()) {
                    let p = PathBuf::from(p);
                    if p.exists() {
                        match stage_name.as_str() {
                            "sort" | "mask" => current_fasta = p,
                            "predict" | "update" | "annotate" => current_gff = Some(p),
                            _ => {}
                        }
                    }
                }
                continue;
            }

            dash.stage_start(&genome.name, stage_name);
            batch_state.genome_mut(&genome.name).start_stage(stage_name);
            batch_state.save().ok();

            let result = run_stage(
                stage_name,
                genome,
                &genome_dir,
                &current_fasta,
                current_gff.as_deref(),
                config,
            );

            match result {
                Ok(output_path) => {
                    let completed: Vec<&str> = batch_state
                        .genome_mut(&genome.name)
                        .completed_stages
                        .iter()
                        .map(|s| s.as_str())
                        .collect();
                    dash.stage_done(&genome.name, stage_name, &completed);

                    let output_str = output_path.as_ref().map(|p| p.display().to_string());
                    batch_state
                        .genome_mut(&genome.name)
                        .complete_stage(stage_name, output_str.as_deref());

                    // Chain outputs
                    if let Some(ref p) = output_path {
                        match stage_name.as_str() {
                            "sort" | "mask" => current_fasta = p.clone(),
                            "predict" | "update" | "annotate" => current_gff = Some(p.clone()),
                            _ => {}
                        }
                    }

                    batch_state.save().ok();
                }
                Err(e) => {
                    let msg = format!("{}", e);
                    dash.genome_failed(&genome.name, stage_name, &msg);
                    batch_state
                        .genome_mut(&genome.name)
                        .mark_failed(stage_name, &msg);
                    batch_state.save().ok();
                    failed = true;
                    break;
                }
            }
        }

        if !failed {
            let elapsed = genome_start.elapsed().as_secs_f64();
            batch_state.genome_mut(&genome.name).elapsed_s = elapsed;
            dash.genome_done(&genome.name, elapsed);
            batch_state.save().ok();
        }
    }

    // Final summary
    dash.print_summary(&batch_state);

    // Point users at real downstream commands once at least two genomes
    // finished annotation.
    let (done, _, _) = batch_state.summary();
    if done >= 2 {
        println!("  \x1b[2mNext step:\x1b[0m");
        println!(
            "  \x1b[33mmyconote-cli compare {}/<A>/annotate_out/annotated.gff3 <A.fa> \\\n                       {}/<B>/annotate_out/annotated.gff3 <B.fa>\x1b[0m\n",
            config.output_dir.display(),
            config.output_dir.display()
        );
    }

    Ok(())
}

/// Run a single pipeline stage for a genome. Returns the primary output path.
fn run_stage(
    stage: &str,
    genome: &GenomeEntry,
    genome_dir: &Path,
    current_fasta: &Path,
    current_gff: Option<&Path>,
    config: &BatchConfig,
) -> Result<Option<PathBuf>> {
    match stage {
        "sort" => {
            let output = genome_dir.join(format!("{}_sorted.fa", genome.name));
            let args = vec![
                current_fasta.display().to_string(),
                "--output".to_string(),
                output.display().to_string(),
                "--min-length".to_string(),
                config.min_length.to_string(),
            ];
            run_myconote_stage("sort", &args)?;
            Ok(Some(output))
        }
        "mask" => {
            let output = genome_dir.join(format!("{}_masked.fa", genome.name));
            let args = vec![
                current_fasta.display().to_string(),
                "--output".to_string(),
                output.display().to_string(),
                "--engine".to_string(),
                config.mask_engine.clone(),
                "--threads".to_string(),
                config.threads.to_string(),
            ];
            run_myconote_stage("mask", &args)?;
            Ok(Some(output))
        }
        "train" => {
            let out_dir = genome_dir.join("train_out");
            let mut args = vec![
                current_fasta.display().to_string(),
                "--output".to_string(),
                out_dir.display().to_string(),
                "--threads".to_string(),
                config.threads.to_string(),
            ];
            if let Some(ref sp) = genome.species {
                args.push("--species".to_string());
                args.push(sp.clone());
            }
            run_myconote_stage("train", &args)?;
            Ok(Some(out_dir))
        }
        "predict" => {
            let out_dir = genome_dir.join("predict_out");
            let mut args = vec![
                current_fasta.display().to_string(),
                "--output".to_string(),
                out_dir.display().to_string(),
                "--kingdom".to_string(),
                genome.kingdom.clone(),
                "--locus-prefix".to_string(),
                genome.locus_prefix.clone(),
                "--threads".to_string(),
                config.threads.to_string(),
            ];
            if let Some(ref sp) = genome.species {
                args.push("--species".to_string());
                args.push(sp.clone());
            }
            run_myconote_stage("predict", &args)?;
            let gff = out_dir.join("consensus.gff3");
            Ok(Some(gff))
        }
        "update" => {
            let gff = current_gff.ok_or_else(|| {
                MycoNoteError::BatchError("update requires a GFF3 from predict".to_string())
            })?;
            let out_dir = genome_dir.join("update_out");
            let args = vec![
                gff.display().to_string(),
                "--fasta".to_string(),
                current_fasta.display().to_string(),
                "--output".to_string(),
                out_dir.display().to_string(),
                "--locus-prefix".to_string(),
                genome.locus_prefix.clone(),
                "--threads".to_string(),
                config.threads.to_string(),
            ];
            run_myconote_stage("update", &args)?;
            let updated_gff = out_dir.join("updated.gff3");
            Ok(Some(updated_gff))
        }
        "annotate" => {
            let gff = current_gff.ok_or_else(|| {
                MycoNoteError::BatchError(
                    "annotate requires a GFF3 from predict/update".to_string(),
                )
            })?;
            let out_dir = genome_dir.join("annotate_out");
            let mut args = vec![
                gff.display().to_string(),
                "--fasta".to_string(),
                current_fasta.display().to_string(),
                "--output".to_string(),
                out_dir.display().to_string(),
                "--kingdom".to_string(),
                genome.kingdom.clone(),
                "--locus-prefix".to_string(),
                genome.locus_prefix.clone(),
                "--threads".to_string(),
                config.threads.to_string(),
            ];
            if genome.genetic_code != 1 {
                args.push("--genetic-code".to_string());
                args.push(genome.genetic_code.to_string());
            }
            run_myconote_stage("annotate", &args)?;
            let annotated_gff = out_dir.join("annotated.gff3");
            Ok(Some(annotated_gff))
        }
        "submit" => {
            let gff = current_gff.ok_or_else(|| {
                MycoNoteError::BatchError("submit requires a GFF3 from annotate".to_string())
            })?;
            let out_dir = genome_dir.join("submit_out");
            let mut args = vec![
                gff.display().to_string(),
                "--fasta".to_string(),
                current_fasta.display().to_string(),
                "--output".to_string(),
                out_dir.display().to_string(),
                "--locus-prefix".to_string(),
                genome.locus_prefix.clone(),
            ];
            if genome.genetic_code != 1 {
                args.push("--genetic-code".to_string());
                args.push(genome.genetic_code.to_string());
            }
            run_myconote_stage("submit", &args)?;
            Ok(Some(out_dir))
        }
        other => {
            eprintln!("  Warning: unknown stage '{}', skipping", other);
            Ok(None)
        }
    }
}

/// Run a myconote-cli subcommand as a child process.
fn run_myconote_stage(subcommand: &str, args: &[String]) -> Result<()> {
    let exe = std::env::current_exe().map_err(|e| {
        MycoNoteError::BatchError(format!("cannot determine myconote-cli path: {}", e))
    })?;

    let mut cmd = std::process::Command::new(&exe);
    cmd.arg(subcommand);
    for arg in args {
        cmd.arg(arg);
    }

    let output = cmd
        .output()
        .map_err(|e| MycoNoteError::BatchError(format!("{} failed to start: {}", subcommand, e)))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        let msg = if stderr.is_empty() {
            stdout.to_string()
        } else {
            stderr.to_string()
        };
        return Err(MycoNoteError::BatchError(format!(
            "{} exited with {}: {}",
            subcommand,
            output.status,
            msg.trim()
        )));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_scan_directory() {
        let tmp = std::env::temp_dir().join("myconote_batch_scan");
        let _ = std::fs::remove_dir_all(&tmp);
        std::fs::create_dir_all(&tmp).unwrap();

        // Create test FASTA files
        std::fs::write(tmp.join("genome_a.fa"), ">seq1\nACGT\n").unwrap();
        std::fs::write(tmp.join("genome_b.fasta"), ">seq1\nACGT\n").unwrap();
        std::fs::write(tmp.join("readme.txt"), "not a fasta").unwrap();

        let config = BatchConfig {
            input: tmp.clone(),
            ..Default::default()
        };

        let entries = scan_directory(&tmp, &config).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].name, "genome_a");
        assert_eq!(entries[1].name, "genome_b");
        assert_eq!(entries[0].kingdom, "fungi");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn test_parse_sample_sheet() {
        let tmp = std::env::temp_dir().join("myconote_batch_sheet");
        let _ = std::fs::remove_dir_all(&tmp);
        std::fs::create_dir_all(&tmp).unwrap();

        // Create a dummy FASTA so validation passes
        let fasta_a = tmp.join("a.fa");
        let fasta_b = tmp.join("b.fa");
        std::fs::write(&fasta_a, ">s\nACGT\n").unwrap();
        std::fs::write(&fasta_b, ">s\nACGT\n").unwrap();

        let sheet = tmp.join("samples.tsv");
        let mut f = std::fs::File::create(&sheet).unwrap();
        writeln!(
            f,
            "name\tfasta\tkingdom\tspecies\tgenetic_code\tlocus_prefix"
        )
        .unwrap();
        writeln!(
            f,
            "isolate_A\t{}\tfungi\tsaccharomyces\t1\tISOA",
            fasta_a.display()
        )
        .unwrap();
        writeln!(f, "isolate_B\t{}\tplant\tauto\t1\tGENE", fasta_b.display()).unwrap();

        let config = BatchConfig::default();
        let entries = parse_sample_sheet(&sheet, &config).unwrap();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].name, "isolate_A");
        assert_eq!(entries[0].kingdom, "fungi");
        assert_eq!(entries[0].species, Some("saccharomyces".to_string()));
        assert_eq!(entries[0].locus_prefix, "ISOA");

        assert_eq!(entries[1].name, "isolate_B");
        assert_eq!(entries[1].kingdom, "plant");
        assert!(entries[1].species.is_none()); // "auto" is filtered

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn test_empty_directory() {
        let tmp = std::env::temp_dir().join("myconote_batch_empty");
        let _ = std::fs::remove_dir_all(&tmp);
        std::fs::create_dir_all(&tmp).unwrap();

        let config = BatchConfig {
            input: tmp.clone(),
            ..Default::default()
        };

        let entries = scan_directory(&tmp, &config).unwrap();
        assert!(entries.is_empty());

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn test_batch_config_default() {
        let cfg = BatchConfig::default();
        assert_eq!(cfg.stages.len(), 5);
        assert_eq!(cfg.kingdom, "fungi");
        assert_eq!(cfg.threads, 4);
        assert_eq!(cfg.max_parallel, 2);
        assert!(!cfg.condor);
    }

    #[test]
    fn test_all_stages_order() {
        assert_eq!(
            ALL_STAGES,
            &["sort", "mask", "predict", "annotate", "submit"]
        );
    }
}
