//! Kallisto-based transcript abundance filtering for the `update` command.
//!
//! funannotate's `update` runs Kallisto to estimate per-transcript TPM from
//! RNA-seq reads, then drops transcripts whose abundance falls below a
//! threshold (default 1.0 TPM) before PASA's UTR-extension pass. This
//! eliminates spurious low-coverage Trinity assemblies and assembly
//! artifacts that would otherwise feed PASA noise and inflate isoform
//! counts.
//!
//! MycoNote-CLI exposes the step under `update --kallisto`. When the flag
//! is set, this module:
//!   1. Builds a kallisto index over the transcript FASTA.
//!   2. Runs `kallisto quant` for each provided RNA-seq sample.
//!   3. Parses each sample's `abundance.tsv` into per-transcript TPM.
//!   4. Aggregates across samples (mean TPM per transcript) and writes a
//!      filtered `passing_transcripts.txt` containing transcript IDs that
//!      cleared the threshold in any sample (more permissive than
//!      requiring all samples).
//!
//! The filtered transcript ID list is then passed to PASA's update step
//! as a hard restriction, so PASA only updates gene models with
//! abundance-supported transcripts.
//!
//! The kallisto CLI surface is verified against
//! https://pachterlab.github.io/kallisto/manual (kallisto 0.50.0,
//! checked 2026-04-24):
//!
//!   kallisto index -i <out>.idx <transcripts.fa>
//!   kallisto quant -i <idx> -o <out_dir> <r1> <r2>          # paired
//!   kallisto quant -i <idx> -o <out_dir> --single -l 200 -s 20 <r1>
//!
//! abundance.tsv columns (header row, then data):
//!   target_id  length  eff_length  est_counts  tpm

use crate::utils::error::{MycoNoteError, Result};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

/// One row of a parsed kallisto `abundance.tsv` file.
#[derive(Debug, Clone, PartialEq)]
pub struct AbundanceRow {
    pub target_id: String,
    pub length: u64,
    pub eff_length: f64,
    pub est_counts: f64,
    pub tpm: f64,
}

/// Check whether `kallisto` is on PATH.
pub fn kallisto_available() -> bool {
    which::which("kallisto").is_ok()
}

/// Parse a kallisto `abundance.tsv` file. Returns the row vector in the
/// order kallisto wrote it; downstream callers can index by `target_id`.
///
/// The header line is required and must contain at least the columns
/// `target_id`, `length`, `eff_length`, `est_counts`, `tpm`. Anything
/// else (tab-delimited extras) is ignored.
pub fn parse_abundance_tsv(path: &Path) -> Result<Vec<AbundanceRow>> {
    let f = std::fs::File::open(path).map_err(|e| {
        MycoNoteError::Io(std::io::Error::new(
            e.kind(),
            format!("opening {}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(f);
    let mut lines = reader.lines();

    let header_line = match lines.next() {
        Some(Ok(h)) => h,
        Some(Err(e)) => return Err(MycoNoteError::Io(e)),
        None => {
            return Err(MycoNoteError::InvalidFormat(format!(
                "{}: empty abundance.tsv (kallisto produced no output)",
                path.display()
            )));
        }
    };

    let cols: Vec<&str> = header_line.split('\t').collect();
    let idx_target = col_index(&cols, "target_id", path)?;
    let idx_length = col_index(&cols, "length", path)?;
    let idx_eff = col_index(&cols, "eff_length", path)?;
    let idx_counts = col_index(&cols, "est_counts", path)?;
    let idx_tpm = col_index(&cols, "tpm", path)?;

    let mut rows = Vec::new();
    for (lineno, line) in lines.enumerate() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let need = idx_tpm
            .max(idx_counts)
            .max(idx_eff)
            .max(idx_length)
            .max(idx_target)
            + 1;
        if fields.len() < need {
            return Err(MycoNoteError::InvalidFormat(format!(
                "{}: line {}: expected at least {} tab-separated fields, found {}",
                path.display(),
                lineno + 2, // +1 for header, +1 for 1-based
                need,
                fields.len()
            )));
        }
        let target_id = fields[idx_target].to_string();
        let length: u64 = fields[idx_length].parse().map_err(|_| {
            MycoNoteError::InvalidFormat(format!(
                "{}: line {}: column 'length' is not an integer: {:?}",
                path.display(),
                lineno + 2,
                fields[idx_length]
            ))
        })?;
        let eff_length: f64 = fields[idx_eff].parse().map_err(|_| {
            MycoNoteError::InvalidFormat(format!(
                "{}: line {}: column 'eff_length' is not a number: {:?}",
                path.display(),
                lineno + 2,
                fields[idx_eff]
            ))
        })?;
        let est_counts: f64 = fields[idx_counts].parse().map_err(|_| {
            MycoNoteError::InvalidFormat(format!(
                "{}: line {}: column 'est_counts' is not a number: {:?}",
                path.display(),
                lineno + 2,
                fields[idx_counts]
            ))
        })?;
        let tpm: f64 = fields[idx_tpm].parse().map_err(|_| {
            MycoNoteError::InvalidFormat(format!(
                "{}: line {}: column 'tpm' is not a number: {:?}",
                path.display(),
                lineno + 2,
                fields[idx_tpm]
            ))
        })?;

        rows.push(AbundanceRow {
            target_id,
            length,
            eff_length,
            est_counts,
            tpm,
        });
    }

    Ok(rows)
}

fn col_index(cols: &[&str], wanted: &str, path: &Path) -> Result<usize> {
    cols.iter().position(|c| *c == wanted).ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!(
            "{}: header is missing required column '{}' (got: {})",
            path.display(),
            wanted,
            cols.join(", ")
        ))
    })
}

/// Filter abundance rows to transcripts whose TPM is at least
/// `min_tpm`. Returns the passing `target_id`s in the order they
/// appeared (kallisto orders by index, which is stable across runs).
pub fn filter_by_min_tpm(rows: &[AbundanceRow], min_tpm: f64) -> Vec<String> {
    rows.iter()
        .filter(|r| r.tpm >= min_tpm)
        .map(|r| r.target_id.clone())
        .collect()
}

/// Aggregate per-sample abundance rows into a per-transcript max TPM.
/// "Pass in any sample" semantics: a transcript clears the threshold if
/// its TPM in at least one sample is >= `min_tpm`. Returns the sorted
/// list of passing transcript IDs (BTreeMap → deterministic order).
pub fn aggregate_passing_transcripts(samples: &[Vec<AbundanceRow>], min_tpm: f64) -> Vec<String> {
    let mut max_tpm: BTreeMap<String, f64> = BTreeMap::new();
    for sample in samples {
        for row in sample {
            let entry = max_tpm.entry(row.target_id.clone()).or_insert(0.0);
            if row.tpm > *entry {
                *entry = row.tpm;
            }
        }
    }
    max_tpm
        .into_iter()
        .filter(|(_id, tpm)| *tpm >= min_tpm)
        .map(|(id, _)| id)
        .collect()
}

/// Configuration for one kallisto-based filtering pass.
#[derive(Debug, Clone)]
pub struct KallistoFilter {
    /// Transcript FASTA to index.
    pub transcripts: PathBuf,
    /// RNA-seq input pairs. Each entry is (R1, optional R2). Single-end
    /// samples set R2 to None.
    pub samples: Vec<(PathBuf, Option<PathBuf>)>,
    /// Filter threshold in TPM. funannotate's default is 1.0; users
    /// with low-coverage RNA-seq should drop this with `--kallisto-min-tpm`.
    pub min_tpm: f64,
    /// Output directory for the kallisto index + per-sample quant runs.
    pub out_dir: PathBuf,
    /// Threads for `kallisto quant -t`.
    pub threads: usize,
    /// Single-end fragment-length mean (kallisto requires `-l` for
    /// single-end). Default 200.
    pub single_frag_len: f64,
    /// Single-end fragment-length sd (kallisto requires `-s` for
    /// single-end). Default 20.
    pub single_frag_sd: f64,
}

impl Default for KallistoFilter {
    fn default() -> Self {
        Self {
            transcripts: PathBuf::new(),
            samples: Vec::new(),
            min_tpm: 1.0,
            out_dir: PathBuf::from("kallisto_out"),
            threads: 4,
            single_frag_len: 200.0,
            single_frag_sd: 20.0,
        }
    }
}

/// Run the full kallisto pipeline: index → per-sample quant → aggregate.
/// Returns the list of passing transcript IDs (sorted, deterministic).
///
/// Caller is responsible for having already verified `kallisto_available()`.
pub fn run_kallisto_filter(cfg: &KallistoFilter) -> Result<Vec<String>> {
    if cfg.samples.is_empty() {
        return Err(MycoNoteError::InvalidFormat(
            "kallisto filter: no RNA-seq samples provided. Pass --rna-r1 (and optionally --rna-r2) or supply --transcripts plus a separate --rna-bam path.".to_string(),
        ));
    }

    if !cfg.transcripts.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "kallisto filter: transcript FASTA not found: {}",
            cfg.transcripts.display()
        )));
    }

    std::fs::create_dir_all(&cfg.out_dir).map_err(MycoNoteError::Io)?;

    // ── 1. Build index ───────────────────────────────────────────────────
    let idx_path = cfg.out_dir.join("transcripts.idx");
    println!(
        "    kallisto index → {}  ({})",
        idx_path.display(),
        cfg.transcripts.display()
    );
    let idx_status = Command::new("kallisto")
        .args([
            "index",
            "-i",
            idx_path.to_str().unwrap_or(""),
            cfg.transcripts.to_str().unwrap_or(""),
        ])
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("kallisto index: {}", e)))?;

    if !idx_status.success() {
        return Err(MycoNoteError::ExternalTool(
            "kallisto index exited with non-zero status".to_string(),
        ));
    }

    // ── 2. Per-sample quant ──────────────────────────────────────────────
    let mut all_rows: Vec<Vec<AbundanceRow>> = Vec::with_capacity(cfg.samples.len());
    for (i, (r1, r2_opt)) in cfg.samples.iter().enumerate() {
        let sample_dir = cfg.out_dir.join(format!("sample_{:02}", i + 1));
        std::fs::create_dir_all(&sample_dir).map_err(MycoNoteError::Io)?;

        println!(
            "    kallisto quant sample {:02} → {}",
            i + 1,
            sample_dir.display()
        );

        let mut cmd = Command::new("kallisto");
        cmd.arg("quant");
        cmd.args(["-i", idx_path.to_str().unwrap_or("")]);
        cmd.args(["-o", sample_dir.to_str().unwrap_or("")]);
        cmd.args(["-t", &cfg.threads.to_string()]);
        if let Some(r2) = r2_opt {
            cmd.arg(r1.to_str().unwrap_or(""));
            cmd.arg(r2.to_str().unwrap_or(""));
        } else {
            // Single-end: kallisto requires -l and -s
            cmd.arg("--single");
            cmd.args(["-l", &format!("{}", cfg.single_frag_len)]);
            cmd.args(["-s", &format!("{}", cfg.single_frag_sd)]);
            cmd.arg(r1.to_str().unwrap_or(""));
        }

        let st = cmd
            .status()
            .map_err(|e| MycoNoteError::ExternalTool(format!("kallisto quant: {}", e)))?;
        if !st.success() {
            return Err(MycoNoteError::ExternalTool(format!(
                "kallisto quant failed for sample {} ({})",
                i + 1,
                r1.display()
            )));
        }

        let abu = sample_dir.join("abundance.tsv");
        let rows = parse_abundance_tsv(&abu)?;
        all_rows.push(rows);
    }

    // ── 3. Aggregate + write list ────────────────────────────────────────
    let passing = aggregate_passing_transcripts(&all_rows, cfg.min_tpm);
    let passing_path = cfg.out_dir.join("passing_transcripts.txt");
    let mut f = std::fs::File::create(&passing_path).map_err(MycoNoteError::Io)?;
    for id in &passing {
        writeln!(f, "{}", id).map_err(MycoNoteError::Io)?;
    }

    let total_targets: usize = all_rows.first().map(|r| r.len()).unwrap_or(0);
    println!(
        "    kallisto filter: {} / {} transcripts passed (min_tpm={})",
        passing.len(),
        total_targets,
        cfg.min_tpm
    );
    println!("    Passing list:    {}", passing_path.display());

    Ok(passing)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_fixture(body: &str) -> tempfile::TempDir {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("abundance.tsv");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        dir
    }

    #[test]
    fn parse_minimal_abundance_tsv() {
        // Three transcripts: low / medium / high TPM. Tab-separated,
        // trailing newline like kallisto writes.
        let body = "target_id\tlength\teff_length\test_counts\ttpm\n\
                    TX_001\t1500\t1300.5\t12.4\t0.5\n\
                    TX_002\t800\t650.0\t250.0\t15.7\n\
                    TX_003\t2200\t2050.0\t1000.0\t250.3\n";
        let dir = write_fixture(body);
        let path = dir.path().join("abundance.tsv");
        let rows = parse_abundance_tsv(&path).unwrap();

        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].target_id, "TX_001");
        assert_eq!(rows[0].length, 1500);
        assert!((rows[0].eff_length - 1300.5).abs() < 1e-9);
        assert!((rows[0].est_counts - 12.4).abs() < 1e-9);
        assert!((rows[0].tpm - 0.5).abs() < 1e-9);
        assert_eq!(rows[2].target_id, "TX_003");
        assert!((rows[2].tpm - 250.3).abs() < 1e-9);
    }

    #[test]
    fn parse_tolerates_extra_columns_and_reorders() {
        // kallisto's column order is fixed but defensive parsing should
        // tolerate any order so long as the required columns are present.
        let body = "tpm\ttarget_id\test_counts\teff_length\tlength\n\
                    1.0\tTX_A\t5.0\t900.0\t1000\n";
        let dir = write_fixture(body);
        let path = dir.path().join("abundance.tsv");
        let rows = parse_abundance_tsv(&path).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].target_id, "TX_A");
        assert_eq!(rows[0].length, 1000);
        assert!((rows[0].tpm - 1.0).abs() < 1e-9);
    }

    #[test]
    fn parse_rejects_missing_required_column() {
        let body = "target_id\tlength\teff_length\test_counts\n\
                    TX\t100\t90.0\t1.0\n";
        let dir = write_fixture(body);
        let path = dir.path().join("abundance.tsv");
        let err = parse_abundance_tsv(&path).unwrap_err();
        let msg = format!("{}", err);
        assert!(
            msg.contains("tpm") && msg.contains("missing"),
            "error must name the missing column: {}",
            msg
        );
    }

    #[test]
    fn parse_rejects_non_numeric_tpm() {
        let body = "target_id\tlength\teff_length\test_counts\ttpm\n\
                    TX\t100\t90.0\t1.0\tnot_a_number\n";
        let dir = write_fixture(body);
        let path = dir.path().join("abundance.tsv");
        let err = parse_abundance_tsv(&path).unwrap_err();
        let msg = format!("{}", err);
        assert!(
            msg.contains("tpm") && msg.contains("number"),
            "error must name the bad column: {}",
            msg
        );
    }

    #[test]
    fn parse_empty_file_errors_clearly() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("abundance.tsv");
        std::fs::write(&path, "").unwrap();
        let err = parse_abundance_tsv(&path).unwrap_err();
        assert!(format!("{}", err).contains("empty"));
    }

    #[test]
    fn filter_at_threshold_is_inclusive() {
        // 1.0 TPM should be kept when threshold is 1.0 (inclusive).
        let rows = vec![
            AbundanceRow {
                target_id: "low".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.5,
                tpm: 0.5,
            },
            AbundanceRow {
                target_id: "edge".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 1.0,
                tpm: 1.0,
            },
            AbundanceRow {
                target_id: "high".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 100.0,
                tpm: 100.0,
            },
        ];
        let kept = filter_by_min_tpm(&rows, 1.0);
        assert_eq!(kept, vec!["edge".to_string(), "high".to_string()]);
    }

    #[test]
    fn filter_drops_zero_tpm() {
        let rows = vec![AbundanceRow {
            target_id: "ghost".into(),
            length: 100,
            eff_length: 100.0,
            est_counts: 0.0,
            tpm: 0.0,
        }];
        let kept = filter_by_min_tpm(&rows, 1.0);
        assert!(kept.is_empty());
    }

    #[test]
    fn aggregate_uses_pass_in_any_sample_semantics() {
        // TX_A: 0.5 / 2.0 / 0.0 → max 2.0, passes at min_tpm=1.0
        // TX_B: 0.5 / 0.5 / 0.5 → max 0.5, fails
        // TX_C: 0.0 / 0.0 / 5.0 → max 5.0, passes
        let s1 = vec![
            AbundanceRow {
                target_id: "TX_A".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 0.5,
            },
            AbundanceRow {
                target_id: "TX_B".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 0.5,
            },
            AbundanceRow {
                target_id: "TX_C".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 0.0,
            },
        ];
        let s2 = vec![
            AbundanceRow {
                target_id: "TX_A".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 2.0,
            },
            AbundanceRow {
                target_id: "TX_B".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 0.5,
            },
            AbundanceRow {
                target_id: "TX_C".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 0.0,
            },
        ];
        let s3 = vec![
            AbundanceRow {
                target_id: "TX_A".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 0.0,
            },
            AbundanceRow {
                target_id: "TX_B".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 0.5,
            },
            AbundanceRow {
                target_id: "TX_C".into(),
                length: 100,
                eff_length: 100.0,
                est_counts: 0.0,
                tpm: 5.0,
            },
        ];
        let passing = aggregate_passing_transcripts(&[s1, s2, s3], 1.0);
        assert_eq!(passing, vec!["TX_A".to_string(), "TX_C".to_string()]);
    }

    #[test]
    fn run_kallisto_filter_errors_when_no_samples() {
        let dir = tempfile::tempdir().unwrap();
        let tx = dir.path().join("tx.fa");
        std::fs::write(&tx, ">a\nACGT\n").unwrap();
        let cfg = KallistoFilter {
            transcripts: tx,
            samples: vec![],
            min_tpm: 1.0,
            out_dir: dir.path().join("out"),
            threads: 1,
            single_frag_len: 200.0,
            single_frag_sd: 20.0,
        };
        let err = run_kallisto_filter(&cfg).unwrap_err();
        let msg = format!("{}", err);
        assert!(
            msg.to_lowercase().contains("no rna-seq samples"),
            "error must call out the missing samples: {}",
            msg
        );
    }

    /// End-to-end smoke test against a real `kallisto` binary. Gated
    /// under `#[ignore]` because it's a multi-MB fixture with a
    /// bioconda-installable dependency, matching the same convention
    /// used by `predict::braker` and `de-template`.
    #[test]
    #[ignore]
    fn kallisto_e2e_smoke() {
        if !kallisto_available() {
            eprintln!("skipping: kallisto not on PATH");
            return;
        }
        // A real test would build a tiny transcript FASTA, simulate
        // reads with wgsim, run quant, and assert the high-coverage
        // transcript passes the 1.0-TPM filter. Left as a placeholder
        // so the test list signals the gating semantics.
    }
}
