//! RNA-seq expression quantification for fungal genomes.
//!
//! Full pipeline: parse a sample sheet → build (or reuse) a decoy-aware
//! salmon index → per sample, fastp QC/trim then `salmon quant` → merge
//! the per-sample `quant.sf` files into a wide count matrix → write a
//! reproducibility manifest. DE analysis stays in R.
//!
//! The pipeline is additive to the existing annotation flow: users run
//! `predict → annotate`, then `convert --to cds` to derive the
//! transcript FASTA salmon indexes, then `quant` here.
//!
//! See `scratch/rnaseq_spec.md` and `scratch/rnaseq_spec_decisions.md`
//! for the design that drives this module.

pub mod bundle;
pub mod fastp;
pub mod index;
pub mod merge;
pub mod salmon;
pub mod sample_sheet;

use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

/// Parsed command-line configuration for a `quant` run.
#[derive(Debug, Clone)]
pub struct QuantConfig {
    pub samples_sheet: PathBuf,
    pub cds_fa: PathBuf,
    pub genome_fa: PathBuf,
    pub output_dir: PathBuf,
    pub threads: usize,
    /// Outer sample concurrency. Per `scratch/rnaseq_spec_decisions.md`
    /// §5 we ship 0.3.0 with a default of 1 (no outer parallelism)
    /// until the A. niger benchmark picks a better value. Users who
    /// want to experiment can override with `--jobs N` and pay the
    /// "no supported default yet" note in the docs.
    pub jobs: usize,
    pub k: usize,
    pub tmpdir: Option<PathBuf>,
    pub index_cache: Option<PathBuf>,
    /// `Some(dir)` → persist trimmed FASTQs under `dir/` rather than
    /// cleaning them up after each sample's salmon quant finishes.
    pub keep_trimmed: Option<PathBuf>,
    pub skip_fastp: bool,
    pub rebuild_index: bool,
    pub seed: u64,
    pub fastp_bin: Option<String>,
    pub salmon_bin: Option<String>,
    /// Raw argv for the bundle's `command` field. Preserved verbatim.
    pub raw_cmdline: String,
}

impl QuantConfig {
    fn threads_default() -> usize {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
    }

    fn parse(args: &[String]) -> Result<Self> {
        let mut samples_sheet: Option<PathBuf> = None;
        let mut genome_fa: Option<PathBuf> = None;
        let mut output_dir: PathBuf = PathBuf::from("quant_out");
        let mut threads = Self::threads_default();
        let mut jobs: usize = 1;
        let mut k: usize = 31;
        let mut tmpdir: Option<PathBuf> = None;
        let mut index_cache: Option<PathBuf> = None;
        let mut keep_trimmed: Option<PathBuf> = None;
        let mut skip_fastp = false;
        let mut rebuild_index = false;
        let mut seed: u64 = 42;
        let mut fastp_bin: Option<String> = None;
        let mut salmon_bin: Option<String> = None;
        let mut positional: Vec<String> = Vec::new();

        let mut i = 0usize;
        while i < args.len() {
            match args[i].as_str() {
                "--samples" | "-s" if i + 1 < args.len() => {
                    samples_sheet = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--genome" | "-g" if i + 1 < args.len() => {
                    genome_fa = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--output" | "-o" if i + 1 < args.len() => {
                    output_dir = PathBuf::from(&args[i + 1]);
                    i += 2;
                }
                "--threads" | "-t" if i + 1 < args.len() => {
                    threads = args[i + 1].parse().unwrap_or(threads).max(1);
                    i += 2;
                }
                "--jobs" | "-j" if i + 1 < args.len() => {
                    jobs = args[i + 1].parse().unwrap_or(jobs).max(1);
                    i += 2;
                }
                "-k" if i + 1 < args.len() => {
                    k = args[i + 1].parse().unwrap_or(k);
                    i += 2;
                }
                "--tmpdir" if i + 1 < args.len() => {
                    tmpdir = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--index-cache" if i + 1 < args.len() => {
                    index_cache = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--keep-trimmed" if i + 1 < args.len() => {
                    keep_trimmed = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--skip-fastp" => {
                    skip_fastp = true;
                    i += 1;
                }
                "--rebuild-index" => {
                    rebuild_index = true;
                    i += 1;
                }
                "--seed" if i + 1 < args.len() => {
                    seed = args[i + 1].parse().unwrap_or(seed);
                    i += 2;
                }
                "--fastp" if i + 1 < args.len() => {
                    fastp_bin = Some(args[i + 1].clone());
                    i += 2;
                }
                "--salmon" if i + 1 < args.len() => {
                    salmon_bin = Some(args[i + 1].clone());
                    i += 2;
                }
                other if other.starts_with("--") => {
                    return Err(MycoNoteError::QuantSheet(format!(
                        "unknown flag '{other}' for quant"
                    )));
                }
                _ => {
                    positional.push(args[i].clone());
                    i += 1;
                }
            }
        }

        // First positional is the CDS FASTA (transcripts). We require
        // it as a positional so the common case is short:
        //   myconote-cli quant cds.fa --samples s.tsv --genome g.fa
        let cds_fa = if let Some(cds) = positional.into_iter().next() {
            PathBuf::from(cds)
        } else {
            return Err(MycoNoteError::QuantSheet(
                "quant: missing CDS FASTA. Usage:\n  \
                 myconote-cli quant <cds.fa> --samples <sheet.tsv> --genome <genome.fa>"
                    .to_string(),
            ));
        };

        let samples_sheet = samples_sheet.ok_or_else(|| {
            MycoNoteError::QuantSheet("quant: --samples <sheet.tsv> is required".to_string())
        })?;
        let genome_fa = genome_fa.ok_or_else(|| {
            MycoNoteError::QuantSheet(
                "quant: --genome <genome.fa> is required (used as salmon decoy)".to_string(),
            )
        })?;

        let raw_cmdline = std::iter::once("myconote-cli quant".to_string())
            .chain(args.iter().cloned())
            .collect::<Vec<_>>()
            .join(" ");

        Ok(Self {
            samples_sheet,
            cds_fa,
            genome_fa,
            output_dir,
            threads,
            jobs,
            k,
            tmpdir,
            index_cache,
            keep_trimmed,
            skip_fastp,
            rebuild_index,
            seed,
            fastp_bin,
            salmon_bin,
            raw_cmdline,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Dispatcher
// ─────────────────────────────────────────────────────────────────────────────

/// CLI entry point: `myconote-cli quant …`. Top-level help is handled
/// by `src/main.rs`; this path assumes the caller already checked for
/// the `--help` flag and chose to dispatch here.
pub fn run_quant(args: &[String]) -> Result<()> {
    let cfg = QuantConfig::parse(args)?;
    dispatch(&cfg)
}

fn dispatch(cfg: &QuantConfig) -> Result<()> {
    // ── 1. Input validation (fail fast) ──────────────────────────────────
    require_file(&cfg.cds_fa, "--to cds output (CDS FASTA)")?;
    require_file(&cfg.genome_fa, "--genome genome FASTA")?;
    require_file(&cfg.samples_sheet, "--samples sheet")?;

    // ── 2. Load the sample sheet. Errors carry line numbers. ─────────────
    let sheet = sample_sheet::parse_sheet(&cfg.samples_sheet)?;
    eprintln!(
        "  ✓ loaded {} samples from {}",
        sheet.samples.len(),
        cfg.samples_sheet.display()
    );

    // Validate every FASTQ exists up front — surfacing all missing
    // files in one shot is much friendlier than discovering them
    // sample-by-sample at salmon time.
    let mut missing: Vec<String> = Vec::new();
    for s in &sheet.samples {
        if !s.fastq_r1.is_file() {
            missing.push(format!(
                "sample {}: fastq_r1 not found at {}",
                s.sample_id,
                s.fastq_r1.display()
            ));
        }
        if let Some(ref r2) = s.fastq_r2 {
            if !r2.is_file() {
                missing.push(format!(
                    "sample {}: fastq_r2 not found at {}",
                    s.sample_id,
                    r2.display()
                ));
            }
        }
    }
    if !missing.is_empty() {
        return Err(MycoNoteError::QuantSheet(format!(
            "{} missing FASTQ file(s):\n  - {}",
            missing.len(),
            missing.join("\n  - ")
        )));
    }

    std::fs::create_dir_all(&cfg.output_dir)?;

    // ── 3. Resolve and build / reuse the salmon index. ──────────────────
    let cache_env = index::CacheEnv::from_process(cfg.index_cache.clone());
    let cache_root = index::resolve_cache_root(&cache_env)?;
    if cfg.rebuild_index {
        // User explicitly asked for a rebuild — delete only the entry
        // for this input tuple's hash so other cached indices are
        // preserved.
        let salmon = cfg.salmon_bin.as_deref().unwrap_or("salmon");
        let version = index::detect_salmon_version(salmon)?;
        let key = index::compute_cache_key(&cfg.cds_fa, &cfg.genome_fa, &version, cfg.k)?;
        let doomed = cache_root.join(&key);
        if doomed.exists() {
            let _ = std::fs::remove_dir_all(&doomed);
            eprintln!(
                "  ↪ --rebuild-index: cleared cache entry {}",
                doomed.display()
            );
        }
    }
    let ispec = index::IndexSpec {
        cds_fa: cfg.cds_fa.clone(),
        genome_fa: cfg.genome_fa.clone(),
        k: cfg.k,
        threads: cfg.threads,
        salmon_bin: cfg.salmon_bin.clone(),
    };
    let idx = index::build_or_reuse_index(&ispec, &cache_root)?;
    eprintln!(
        "  ✓ salmon index {} (k={}, targets={}, decoys={}) {}",
        idx.path.display(),
        idx.k,
        idx.targets,
        idx.decoys,
        if idx.cached { "[cached]" } else { "[fresh]" }
    );

    // ── 4. Per-sample pipeline. jobs=1 for 0.3.0 (benchmark pending). ───
    let tmpdir = resolve_tmpdir(cfg)?;
    let fastp_subdir = cfg.output_dir.join("fastp");
    let salmon_subdir = cfg.output_dir.join("salmon");
    std::fs::create_dir_all(&fastp_subdir)?;
    std::fs::create_dir_all(&salmon_subdir)?;

    let mut sample_summaries: Vec<bundle::SampleSummary> = Vec::new();
    let mut merge_inputs: Vec<(String, PathBuf)> = Vec::new();

    if cfg.jobs != 1 {
        eprintln!(
            "  ⚠ --jobs {} requested; 0.3.0 serial path is the only one benchmarked — \
             running samples sequentially regardless.",
            cfg.jobs
        );
    }

    for sample in &sheet.samples {
        eprintln!("  → {}", sample.sample_id);

        // Per-sample tmpdir so filename collisions between parallel
        // samples are impossible even when the outer loop goes
        // parallel in a later release.
        let sample_tmp = tmpdir.join(format!("myconote_quant_{}", sample.sample_id));
        let _ = std::fs::remove_dir_all(&sample_tmp);
        std::fs::create_dir_all(&sample_tmp)?;

        let trimmed_r1: PathBuf;
        let trimmed_r2: Option<PathBuf>;
        let fastp_json: PathBuf;
        let fastp_html: Option<PathBuf>;
        let _fastp_output_guard: Option<fastp::FastpOutput>;

        if cfg.skip_fastp {
            // Trust the user — feed FASTQs directly to salmon.
            trimmed_r1 = sample.fastq_r1.clone();
            trimmed_r2 = sample.fastq_r2.clone();
            fastp_json = PathBuf::new();
            fastp_html = None;
            _fastp_output_guard = None;
        } else {
            let fspec = fastp::FastpSpec {
                sample,
                tmpdir: &sample_tmp,
                threads: cfg.threads,
                fastp_bin: cfg.fastp_bin.as_deref(),
                persist_trimmed: cfg.keep_trimmed.is_some(),
            };
            let f_out = fastp::run_fastp(&fspec)?;

            // Copy fastp JSON into quant_out/fastp/<sample>.json so
            // the permanent artifact layout doesn't depend on tmpdir
            // survival. HTML likewise, for human browsing.
            let persisted_json = fastp_subdir.join(format!("{}.json", sample.sample_id));
            std::fs::copy(&f_out.json_report, &persisted_json)?;
            let persisted_html = fastp_subdir.join(format!("{}.html", sample.sample_id));
            let _ = std::fs::copy(&f_out.html_report, &persisted_html);

            trimmed_r1 = f_out.trimmed_r1.clone();
            trimmed_r2 = f_out.trimmed_r2.clone();
            fastp_json = persisted_json;
            fastp_html = Some(persisted_html);
            _fastp_output_guard = Some(f_out); // keep the drop guard alive through salmon quant
        }

        // Run salmon quant.
        let sample_quant_dir = salmon_subdir.join(&sample.sample_id);
        let sspec = salmon::SalmonQuantSpec {
            sample,
            index_path: &idx.path,
            trimmed_r1: &trimmed_r1,
            trimmed_r2: trimmed_r2.as_deref(),
            output_dir: &sample_quant_dir,
            threads: cfg.threads,
            salmon_bin: cfg.salmon_bin.as_deref(),
        };
        let qout = salmon::run_salmon_quant(&sspec)?;
        merge_inputs.push((sample.sample_id.clone(), qout.quant_sf.clone()));

        // Per-sample bundle entry.
        let fsum = if cfg.skip_fastp {
            bundle::FastpSummary {
                reads_before_filtering: 0,
                reads_after_filtering: 0,
                reads_passed_pct: 0.0,
                q30_rate_before: 0.0,
                q30_rate_after: 0.0,
                adapter_trimmed_reads: None,
                adapter_trimmed_bases: None,
                duplication_rate: 0.0,
                insert_size_peak: None,
            }
        } else {
            bundle::extract_fastp_summary(&fastp_json)?
        };

        sample_summaries.push(bundle::SampleSummary {
            sample_id: sample.sample_id.clone(),
            mapping_rate: qout.mapping_rate,
            library_size: qout.library_size,
            reads_before_filtering: fsum.reads_before_filtering,
            reads_after_filtering: fsum.reads_after_filtering,
            reads_passed_pct: fsum.reads_passed_pct,
            q30_rate_before: fsum.q30_rate_before,
            q30_rate_after: fsum.q30_rate_after,
            adapter_trimmed_reads: fsum.adapter_trimmed_reads,
            adapter_trimmed_bases: fsum.adapter_trimmed_bases,
            duplication_rate: fsum.duplication_rate,
            insert_size_peak: fsum.insert_size_peak,
        });

        // Persist trimmed FASTQs to --keep-trimmed target if requested.
        if let (Some(keep_dir), Some(ref fout)) = (&cfg.keep_trimmed, _fastp_output_guard.as_ref())
        {
            std::fs::create_dir_all(keep_dir)?;
            let dest = keep_dir.join(format!("{}_R1.fq.gz", sample.sample_id));
            let _ = std::fs::copy(&fout.trimmed_r1, &dest);
            if let Some(ref r2) = fout.trimmed_r2 {
                let dest2 = keep_dir.join(format!("{}_R2.fq.gz", sample.sample_id));
                let _ = std::fs::copy(r2, &dest2);
            }
        }
        let _ = fastp_html; // unused past logging; keep var for future UI work.
    }

    // ── 5. Merge → counts.tsv + tpm.tsv ──────────────────────────────────
    let merged = merge::merge(&merge_inputs)?;
    let counts_path = cfg.output_dir.join("counts.tsv");
    let tpm_path = cfg.output_dir.join("tpm.tsv");
    merged.write_counts(&counts_path)?;
    merged.write_tpm(&tpm_path)?;
    eprintln!(
        "  ✓ merged {} transcripts × {} samples",
        merged.transcripts.len(),
        merged.samples.len()
    );

    // ── 6. Write bundle JSON + normalized sample sheet ───────────────────
    let inputs = bundle::BundleInputs {
        transcripts: bundle::HashedPath {
            path: cfg.cds_fa.to_string_lossy().into_owned(),
            sha256: bundle::sha256_file(&cfg.cds_fa)?,
        },
        genome: bundle::HashedPath {
            path: cfg.genome_fa.to_string_lossy().into_owned(),
            sha256: bundle::sha256_file(&cfg.genome_fa)?,
        },
        sample_sheet: bundle::HashedPath {
            path: cfg.samples_sheet.to_string_lossy().into_owned(),
            sha256: bundle::sha256_file(&cfg.samples_sheet)?,
        },
    };
    let bundle_index = bundle::BundleIndex {
        path: idx.path.to_string_lossy().into_owned(),
        sha256: idx.cache_key.clone(),
        k: idx.k,
        decoys: idx.decoys,
        targets: idx.targets,
    };
    let fastp_version = detect_fastp_version(cfg.fastp_bin.as_deref().unwrap_or("fastp"))
        .unwrap_or_else(|_| "unknown".to_string());
    let bundle_obj = bundle::QuantBundle {
        version: env!("CARGO_PKG_VERSION").to_string(),
        run_timestamp: rfc3339_utc_now(),
        command: cfg.raw_cmdline.clone(),
        external_tools: bundle::ExternalTools {
            fastp: fastp_version,
            salmon: idx.salmon_version.clone(),
        },
        inputs,
        index: bundle_index,
        samples: sample_summaries,
        seed: cfg.seed,
    };
    let bundle_path = cfg.output_dir.join("quant_bundle.json");
    bundle_obj.write(&bundle_path)?;

    // Normalized sample sheet — for now, just copy the input verbatim
    // so downstream runs have a reproducible pointer. A full
    // normalization pass (absolute paths + derived QC stats appended)
    // can land in a later revision without breaking the bundle shape.
    let normalized_sheet = cfg.output_dir.join("sample_sheet.tsv");
    std::fs::copy(&cfg.samples_sheet, &normalized_sheet)?;

    eprintln!(
        "  ✓ wrote {} and {}",
        counts_path.display(),
        tpm_path.display()
    );
    eprintln!("  ✓ bundle {}", bundle_path.display());
    eprintln!(
        "\nNext step: load the matrix into R via tximport on \
         {}/salmon/<sample>/quant.sf and run DESeq2.",
        cfg.output_dir.display()
    );

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn require_file(p: &Path, role: &str) -> Result<()> {
    if !p.is_file() {
        return Err(MycoNoteError::QuantSheet(format!(
            "{role}: file not found at {}",
            p.display()
        )));
    }
    Ok(())
}

fn resolve_tmpdir(cfg: &QuantConfig) -> Result<PathBuf> {
    let root = cfg.tmpdir.clone().unwrap_or_else(std::env::temp_dir);
    std::fs::create_dir_all(&root)?;
    Ok(root)
}

/// Probe fastp's `--version` output (format: `fastp 1.3.2`). Kept here
/// rather than in `fastp.rs` because the version string only matters
/// for the bundle manifest — the runtime wrapper doesn't key on it.
fn detect_fastp_version(fastp_bin: &str) -> Result<String> {
    let output = duct::cmd!(fastp_bin, "--version")
        .stderr_to_stdout()
        .read()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "fastp".to_string(),
            message: format!("running '{fastp_bin} --version': {e}"),
        })?;
    for line in output.lines() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.first().map(|t| t.eq_ignore_ascii_case("fastp")) == Some(true) && toks.len() >= 2 {
            // fastp prints "fastp v1.3.2" on stderr; strip the leading 'v'.
            let v = toks[1].trim_start_matches('v').to_string();
            return Ok(v);
        }
    }
    Err(MycoNoteError::QuantTool {
        tool: "fastp".to_string(),
        message: format!("unexpected '--version' output: {}", output.trim()),
    })
}

/// RFC 3339 timestamp in UTC. Rolled inline rather than pulling a
/// full date crate — the bundle is the only consumer and the format
/// is fixed by the spec.
fn rfc3339_utc_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    // Convert to broken-down UTC. Pure arithmetic; avoids libc's
    // gmtime (which touches env vars).
    format_rfc3339(secs)
}

fn format_rfc3339(secs: i64) -> String {
    // Days from epoch + seconds of day.
    let days = secs.div_euclid(86400);
    let sod = secs.rem_euclid(86400);
    let hour = sod / 3600;
    let minute = (sod / 60) % 60;
    let second = sod % 60;

    // Civil-from-days algorithm (Howard Hinnant, public domain).
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };

    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        y, m, d, hour, minute, second
    )
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(s: &str) -> Vec<String> {
        s.split_whitespace().map(String::from).collect()
    }

    #[test]
    fn parse_minimal_ok() {
        let a = argv("cds.fa --samples s.tsv --genome g.fa");
        let cfg = QuantConfig::parse(&a).unwrap();
        assert_eq!(cfg.cds_fa, PathBuf::from("cds.fa"));
        assert_eq!(cfg.samples_sheet, PathBuf::from("s.tsv"));
        assert_eq!(cfg.genome_fa, PathBuf::from("g.fa"));
        assert_eq!(cfg.output_dir, PathBuf::from("quant_out"));
        assert_eq!(cfg.jobs, 1);
        assert_eq!(cfg.k, 31);
        assert_eq!(cfg.seed, 42);
    }

    #[test]
    fn parse_missing_cds_errors() {
        let a = argv("--samples s.tsv --genome g.fa");
        let err = QuantConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("missing CDS FASTA"));
    }

    #[test]
    fn parse_missing_samples_errors() {
        let a = argv("cds.fa --genome g.fa");
        let err = QuantConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--samples"));
    }

    #[test]
    fn parse_missing_genome_errors() {
        let a = argv("cds.fa --samples s.tsv");
        let err = QuantConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--genome"));
    }

    #[test]
    fn parse_all_flags() {
        let a = argv(
            "cds.fa --samples s.tsv --genome g.fa \
             --output out -t 8 -j 2 -k 25 \
             --tmpdir /tmp/mc --index-cache /opt/cache \
             --keep-trimmed /keep --skip-fastp --rebuild-index \
             --seed 7 --fastp /bin/fastp --salmon /bin/salmon",
        );
        let cfg = QuantConfig::parse(&a).unwrap();
        assert_eq!(cfg.output_dir, PathBuf::from("out"));
        assert_eq!(cfg.threads, 8);
        assert_eq!(cfg.jobs, 2);
        assert_eq!(cfg.k, 25);
        assert_eq!(cfg.tmpdir.as_deref(), Some(Path::new("/tmp/mc")));
        assert_eq!(cfg.index_cache.as_deref(), Some(Path::new("/opt/cache")));
        assert_eq!(cfg.keep_trimmed.as_deref(), Some(Path::new("/keep")));
        assert!(cfg.skip_fastp);
        assert!(cfg.rebuild_index);
        assert_eq!(cfg.seed, 7);
        assert_eq!(cfg.fastp_bin.as_deref(), Some("/bin/fastp"));
        assert_eq!(cfg.salmon_bin.as_deref(), Some("/bin/salmon"));
    }

    #[test]
    fn parse_rejects_unknown_flag() {
        let a = argv("cds.fa --samples s.tsv --genome g.fa --bogus-flag");
        let err = QuantConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("unknown flag"));
    }

    #[test]
    fn rfc3339_at_unix_epoch() {
        assert_eq!(format_rfc3339(0), "1970-01-01T00:00:00Z");
    }

    #[test]
    fn rfc3339_known_timestamp() {
        // Verified via python:
        //   datetime(2026,4,23,17,30,0,tzinfo=timezone.utc).timestamp() == 1776965400
        assert_eq!(format_rfc3339(1776965400), "2026-04-23T17:30:00Z");
    }
}
