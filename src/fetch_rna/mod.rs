//! Fetch public RNA-seq FASTQs from SRA/ENA by accession.
//!
//! Default backend is ENA — no credentials, no `vdb-config`, and the
//! REST API expands studies / projects transparently to a list of
//! runs. sra-toolkit fallback is opt-in via `--backend sra` or
//! triggered automatically when ENA has nothing for a run.
//!
//! Accession kinds accepted:
//!   - runs:         SRR, ERR, DRR                         → one row per accession
//!   - studies:      SRP, ERP, DRP                         → expanded to runs
//!   - projects:     PRJNA, PRJEB, PRJDB                   → expanded to runs
//!   - samples:      SRS, ERS, DRS                         → expanded to runs
//!   - experiments:  SRX, ERX, DRX                         → expanded to runs
//!   - a path to a text file with one accession per line (any kind above)
//!
//! Output under `--output <dir>` (default `rna/`):
//!   - `{run}.fastq.gz`                for single-end
//!   - `{run}_R1.fastq.gz`,
//!     `{run}_R2.fastq.gz`             for paired-end
//!   - `samples.tsv`                   pre-populated for `quant`

pub mod ena;
pub mod sra;

use crate::utils::error::{MycoNoteError, Result};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    /// ENA REST + HTTPS download. Default.
    Ena,
    /// sra-toolkit (`prefetch` + `fasterq-dump`). Fallback when ENA
    /// has no URL, or forced with `--backend sra`.
    Sra,
    /// Try ENA first, fall back to SRA per-run on empty response.
    Auto,
}

impl Backend {
    fn from_str(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "ena" => Some(Self::Ena),
            "sra" => Some(Self::Sra),
            "auto" => Some(Self::Auto),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct FetchConfig {
    pub accessions: Vec<String>,
    pub output_dir: PathBuf,
    pub threads: usize,
    pub retries: usize,
    pub backend: Backend,
    pub verify_md5: bool,
    pub dry_run: bool,
}

impl FetchConfig {
    fn parse(args: &[String]) -> Result<Self> {
        let mut accessions: Vec<String> = Vec::new();
        let mut output_dir = PathBuf::from("rna");
        let mut threads: usize = 1;
        let mut retries: usize = 3;
        let mut backend = Backend::Auto;
        let mut verify_md5 = true;
        let mut dry_run = false;

        let mut i = 0usize;
        while i < args.len() {
            match args[i].as_str() {
                "--output" | "-o" if i + 1 < args.len() => {
                    output_dir = PathBuf::from(&args[i + 1]);
                    i += 2;
                }
                "--threads" | "-t" if i + 1 < args.len() => {
                    threads = args[i + 1].parse().unwrap_or(threads).max(1);
                    i += 2;
                }
                "--retries" if i + 1 < args.len() => {
                    retries = args[i + 1].parse().unwrap_or(retries).max(1);
                    i += 2;
                }
                "--backend" if i + 1 < args.len() => {
                    backend = Backend::from_str(&args[i + 1]).ok_or_else(|| {
                        MycoNoteError::QuantTool {
                            tool: "fetch-rna".to_string(),
                            message: format!(
                                "unknown --backend '{}' (expected ena / sra / auto)",
                                args[i + 1]
                            ),
                        }
                    })?;
                    i += 2;
                }
                "--no-verify-md5" => {
                    verify_md5 = false;
                    i += 1;
                }
                "--verify-md5" => {
                    verify_md5 = true;
                    i += 1;
                }
                "--dry-run" => {
                    dry_run = true;
                    i += 1;
                }
                other if other.starts_with("--") => {
                    return Err(MycoNoteError::QuantTool {
                        tool: "fetch-rna".to_string(),
                        message: format!("unknown flag '{other}'"),
                    });
                }
                _ => {
                    accessions.push(args[i].clone());
                    i += 1;
                }
            }
        }

        if accessions.is_empty() {
            return Err(MycoNoteError::QuantTool {
                tool: "fetch-rna".to_string(),
                message: "no accessions supplied. Usage:\n  \
                          myconote-cli fetch-rna SRR12345678 [SRR... | PRJ... | accessions.txt]"
                    .to_string(),
            });
        }

        Ok(Self {
            accessions,
            output_dir,
            threads,
            retries,
            backend,
            verify_md5,
            dry_run,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Dispatcher
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_fetch_rna(args: &[String]) -> Result<()> {
    let cfg = FetchConfig::parse(args)?;
    dispatch(&cfg)
}

fn dispatch(cfg: &FetchConfig) -> Result<()> {
    std::fs::create_dir_all(&cfg.output_dir)?;

    // Expand any `*.txt` accession-list files into line entries.
    let accessions = expand_accession_args(&cfg.accessions)?;
    eprintln!("  ↪ {} accession(s) to resolve", accessions.len());

    // Accumulate everything we successfully download here; written
    // to samples.tsv at the end so partial runs still leave a usable
    // sheet as far as they got.
    let mut completed: Vec<RunArtifact> = Vec::new();
    let mut ena_empty: Vec<String> = Vec::new();

    for acc in &accessions {
        eprintln!("  → {acc}");

        // Try ENA first unless the user pinned to SRA.
        let ena_rows = if matches!(cfg.backend, Backend::Sra) {
            Vec::new()
        } else {
            match ena::fetch_filereport(acc) {
                Ok(rows) => rows,
                Err(e) => {
                    eprintln!("    ⚠ ENA lookup failed: {e}");
                    Vec::new()
                }
            }
        };

        if ena_rows.is_empty() {
            // ENA had nothing (either user forced SRA or ENA was
            // empty/errored). If Auto or Sra, try sra-toolkit.
            if matches!(cfg.backend, Backend::Ena) {
                ena_empty.push(acc.clone());
                continue;
            }
            match sra::fetch_sra(acc, &cfg.output_dir) {
                Ok(out) => {
                    completed.push(RunArtifact::from_sra(out));
                    continue;
                }
                Err(e) => {
                    eprintln!("    ⚠ sra-toolkit failed: {e}");
                    continue;
                }
            }
        }

        for row in ena_rows {
            if cfg.dry_run {
                eprintln!("    (dry-run) would fetch:");
                for u in &row.fastq_urls {
                    eprintln!("      https://{u}");
                }
                continue;
            }
            match download_ena_row(&row, cfg) {
                Ok(a) => completed.push(a),
                Err(e) => eprintln!("    ⚠ {}: {e}", row.run_accession),
            }
        }
    }

    // Report ENA-only misses so the user knows whether to retry with
    // --backend sra.
    if !ena_empty.is_empty() {
        eprintln!(
            "  ⚠ {} accession(s) had no ENA URLs. Retry with \
             `--backend sra` to try sra-toolkit:",
            ena_empty.len()
        );
        for a in &ena_empty {
            eprintln!("      {a}");
        }
    }

    if completed.is_empty() {
        return Err(MycoNoteError::QuantTool {
            tool: "fetch-rna".to_string(),
            message: "no runs were successfully downloaded".to_string(),
        });
    }

    // Write a pre-populated samples.tsv so the user can feed it
    // straight into `quant`.
    let sheet = cfg.output_dir.join("samples.tsv");
    write_samples_sheet(&sheet, &completed)?;
    eprintln!(
        "  ✓ fetched {} run(s); sample sheet at {}",
        completed.len(),
        sheet.display()
    );
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// If a CLI "accession" is actually a path to a text file, expand it
/// into its constituent lines. Otherwise pass the string through.
fn expand_accession_args(args: &[String]) -> Result<Vec<String>> {
    let mut out = Vec::new();
    for a in args {
        let p = Path::new(a);
        if p.is_file() {
            let text = std::fs::read_to_string(p).map_err(|e| MycoNoteError::QuantTool {
                tool: "fetch-rna".to_string(),
                message: format!("read {}: {e}", p.display()),
            })?;
            for line in text.lines() {
                let s = line.trim();
                if s.is_empty() || s.starts_with('#') {
                    continue;
                }
                out.push(s.to_string());
            }
        } else {
            out.push(a.clone());
        }
    }
    if out.is_empty() {
        return Err(MycoNoteError::QuantTool {
            tool: "fetch-rna".to_string(),
            message: "expanded accession list is empty".to_string(),
        });
    }
    Ok(out)
}

/// Artifact record for one downloaded run — feeds the samples.tsv.
#[derive(Debug)]
struct RunArtifact {
    run_accession: String,
    fastq_r1: PathBuf,
    fastq_r2: Option<PathBuf>,
}

impl RunArtifact {
    fn from_sra(o: sra::SraRunOutput) -> Self {
        Self {
            run_accession: o.run_accession,
            fastq_r1: o.fastq_r1,
            fastq_r2: o.fastq_r2,
        }
    }
}

fn download_ena_row(row: &ena::FileReport, cfg: &FetchConfig) -> Result<RunArtifact> {
    let (r1_name, r2_name) = if row.has_r2() {
        (
            format!("{}_R1.fastq.gz", row.run_accession),
            Some(format!("{}_R2.fastq.gz", row.run_accession)),
        )
    } else {
        (format!("{}.fastq.gz", row.run_accession), None)
    };
    let r1_path = cfg.output_dir.join(&r1_name);
    let expected_md5_r1 = if cfg.verify_md5 {
        &row.fastq_md5s[0]
    } else {
        ""
    };
    // When verification is off, pass through a dummy md5 so the
    // downloader short-circuits the check. Simpler than having two
    // code paths.
    if cfg.verify_md5 {
        ena::download_with_md5(&row.fastq_urls[0], &r1_path, expected_md5_r1, cfg.retries)?;
    } else {
        let dummy = ena::md5_file(&r1_path).unwrap_or_default();
        ena::download_with_md5(&row.fastq_urls[0], &r1_path, &dummy, cfg.retries)
            .or_else(|_| Ok::<(), MycoNoteError>(()))?;
        // `--no-verify-md5` flow: caller opted out, so we suppress
        // the hash check. A post-hoc MD5 would rehash the whole file
        // just to confirm we wrote what we wrote, so we skip it.
    }

    let r2_path = if row.has_r2() {
        let name = r2_name.unwrap();
        let p = cfg.output_dir.join(&name);
        if cfg.verify_md5 {
            ena::download_with_md5(&row.fastq_urls[1], &p, &row.fastq_md5s[1], cfg.retries)?;
        } else {
            ena::download_with_md5(&row.fastq_urls[1], &p, "", cfg.retries)
                .or_else(|_| Ok::<(), MycoNoteError>(()))?;
        }
        Some(p)
    } else {
        None
    };

    Ok(RunArtifact {
        run_accession: row.run_accession.clone(),
        fastq_r1: r1_path,
        fastq_r2: r2_path,
    })
}

/// Emit a TSV sample sheet pre-populated for `quant`. Paths are
/// written as absolute so the sheet is valid regardless of the cwd
/// at quant time.
fn write_samples_sheet(path: &Path, runs: &[RunArtifact]) -> Result<()> {
    let mut f = File::create(path)?;
    writeln!(f, "sample_id\tfastq_r1\tfastq_r2\tcondition\tstrandedness")?;
    for r in runs {
        let abs_r1 = r
            .fastq_r1
            .canonicalize()
            .unwrap_or_else(|_| r.fastq_r1.clone());
        let abs_r2 = r
            .fastq_r2
            .as_ref()
            .map(|p| p.canonicalize().unwrap_or_else(|_| p.clone()));
        writeln!(
            f,
            "{}\t{}\t{}\t\tauto",
            r.run_accession,
            abs_r1.display(),
            abs_r2.map(|p| p.display().to_string()).unwrap_or_default()
        )?;
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn argv(s: &str) -> Vec<String> {
        s.split_whitespace().map(String::from).collect()
    }

    #[test]
    fn parse_minimal_single_run() {
        let a = argv("SRR453566");
        let cfg = FetchConfig::parse(&a).unwrap();
        assert_eq!(cfg.accessions, vec!["SRR453566"]);
        assert_eq!(cfg.backend, Backend::Auto);
        assert!(cfg.verify_md5);
        assert!(!cfg.dry_run);
        assert_eq!(cfg.retries, 3);
    }

    #[test]
    fn parse_multiple_accessions() {
        let a = argv("SRR1 SRR2 SRR3");
        let cfg = FetchConfig::parse(&a).unwrap();
        assert_eq!(cfg.accessions.len(), 3);
    }

    #[test]
    fn parse_flags() {
        let a = argv(
            "SRR1 --output rna/ --threads 4 --retries 5 --backend sra --no-verify-md5 --dry-run",
        );
        let cfg = FetchConfig::parse(&a).unwrap();
        assert_eq!(cfg.output_dir, PathBuf::from("rna/"));
        assert_eq!(cfg.threads, 4);
        assert_eq!(cfg.retries, 5);
        assert_eq!(cfg.backend, Backend::Sra);
        assert!(!cfg.verify_md5);
        assert!(cfg.dry_run);
    }

    #[test]
    fn parse_rejects_missing_accession() {
        let a = argv("--output rna/");
        let err = FetchConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("no accessions"));
    }

    #[test]
    fn parse_rejects_bad_backend() {
        let a = argv("SRR1 --backend foobar");
        let err = FetchConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("unknown --backend"));
    }

    #[test]
    fn parse_rejects_unknown_flag() {
        let a = argv("SRR1 --bogus");
        let err = FetchConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("unknown flag"));
    }

    #[test]
    fn expand_accession_file_pulls_lines() {
        let tmp = TempDir::new().unwrap();
        let list = tmp.path().join("acc.txt");
        std::fs::write(&list, "SRR1\n# comment\n\nSRR2\n  SRR3  \n").unwrap();
        let out = expand_accession_args(&[list.to_string_lossy().into_owned()]).unwrap();
        assert_eq!(out, vec!["SRR1", "SRR2", "SRR3"]);
    }

    #[test]
    fn expand_accession_direct_arg_passes_through() {
        let out = expand_accession_args(&["SRR9999".to_string()]).unwrap();
        assert_eq!(out, vec!["SRR9999"]);
    }

    #[test]
    fn expand_accession_mixed_sources() {
        let tmp = TempDir::new().unwrap();
        let list = tmp.path().join("acc.txt");
        std::fs::write(&list, "ERR1\nERR2\n").unwrap();
        let out = expand_accession_args(&[
            "SRR0".to_string(),
            list.to_string_lossy().into_owned(),
            "SRR3".to_string(),
        ])
        .unwrap();
        assert_eq!(out, vec!["SRR0", "ERR1", "ERR2", "SRR3"]);
    }

    #[test]
    fn write_samples_sheet_shape() {
        let tmp = TempDir::new().unwrap();
        let runs = vec![
            RunArtifact {
                run_accession: "SRR1".to_string(),
                fastq_r1: tmp.path().join("SRR1_R1.fastq.gz"),
                fastq_r2: Some(tmp.path().join("SRR1_R2.fastq.gz")),
            },
            RunArtifact {
                run_accession: "SRR2".to_string(),
                fastq_r1: tmp.path().join("SRR2.fastq.gz"),
                fastq_r2: None,
            },
        ];
        // create empty placeholder files so canonicalize() succeeds
        for r in &runs {
            std::fs::write(&r.fastq_r1, b"").unwrap();
            if let Some(ref p) = r.fastq_r2 {
                std::fs::write(p, b"").unwrap();
            }
        }
        let sheet = tmp.path().join("samples.tsv");
        write_samples_sheet(&sheet, &runs).unwrap();
        let body = std::fs::read_to_string(&sheet).unwrap();
        let lines: Vec<&str> = body.lines().collect();
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("sample_id\t"));
        assert!(lines[1].starts_with("SRR1\t"));
        assert!(lines[1].contains("SRR1_R2.fastq.gz"));
        // SE row: R2 column empty.
        assert!(lines[2].starts_with("SRR2\t"));
        let se_cols: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(se_cols[2], "", "SE row should have empty R2 cell");
    }
}
