//! fastp wrapper — QC + adapter trim for one sample.
//!
//! Per the design decision in `scratch/rnaseq_spec_decisions.md` §1,
//! fastp output is written to temp files rather than streamed into
//! salmon via stdin. Salmon needs both `-1` and `-2` on separate file
//! args for paired-end input, which is incompatible with a single
//! child stdin. Temp files trade ~1-2% wall-clock for debuggability
//! and clean interrupt handling; the dispatcher cleans them up after
//! salmon quant finishes for that sample.
//!
//! Each call produces:
//!   - `<tmpdir>/<sample_id>_R1.trimmed.fq.gz`  (always)
//!   - `<tmpdir>/<sample_id>_R2.trimmed.fq.gz`  (paired only)
//!   - `<tmpdir>/<sample_id>.fastp.json`        (summary for bundle)
//!   - `<tmpdir>/<sample_id>.fastp.html`        (persisted when
//!                                               --keep-trimmed is set)
//!
//! The trimmed FASTQs are returned via a `FastpOutput` struct that
//! holds a drop-guard cleanup hook. The dispatcher releases the
//! guard by calling `.persist()` when `--keep-trimmed <dir>` was
//! passed, at which point the trimmed files are moved to the
//! persistence directory instead of deleted.

use crate::quant::sample_sheet::Sample;
use crate::utils::error::{MycoNoteError, Result};
use std::fs;
use std::path::{Path, PathBuf};

/// Inputs to a single-sample fastp run.
pub struct FastpSpec<'a> {
    pub sample: &'a Sample,
    /// Directory to write trimmed + JSON files into. Cleaned at end
    /// of run by the dispatcher; tests pass a TempDir.
    pub tmpdir: &'a Path,
    pub threads: usize,
    /// Path to fastp binary. `None` → "fastp" (PATH lookup).
    pub fastp_bin: Option<&'a str>,
    /// When true, suppress the drop guard — caller takes ownership of
    /// the trimmed files. Used by the `--keep-trimmed` CLI flag.
    pub persist_trimmed: bool,
}

impl<'a> FastpSpec<'a> {
    fn bin(&self) -> &str {
        self.fastp_bin.unwrap_or("fastp")
    }
}

/// Paths + cleanup hook returned by a successful fastp run.
///
/// The `_cleanup` field is a drop guard; when `FastpOutput` goes out
/// of scope the trimmed FASTQs are deleted (unless `persist_trimmed`
/// was set on the spec). The dispatcher either reads the outputs
/// while the struct is alive or explicitly moves them with
/// `forget_cleanup()` after copying to a --keep-trimmed target.
pub struct FastpOutput {
    pub trimmed_r1: PathBuf,
    pub trimmed_r2: Option<PathBuf>,
    pub json_report: PathBuf,
    /// HTML report — kept for users who want to eyeball per-base
    /// quality plots; the dispatcher does not read this file.
    pub html_report: PathBuf,
    _cleanup: TrimmedCleanup,
}

impl FastpOutput {
    /// Disarm the cleanup guard. After this, the trimmed FASTQs stay
    /// on disk until the caller (or process exit) removes them.
    pub fn forget_cleanup(mut self) -> Self {
        self._cleanup.armed = false;
        self
    }
}

/// Drop guard that unlinks the trimmed FASTQs when the `FastpOutput`
/// goes out of scope. Keeps the JSON + HTML reports — those live
/// inside `quant_out/fastp/` after the dispatcher copies them.
pub struct TrimmedCleanup {
    files: Vec<PathBuf>,
    armed: bool,
}

impl Drop for TrimmedCleanup {
    fn drop(&mut self) {
        if !self.armed {
            return;
        }
        for f in &self.files {
            let _ = fs::remove_file(f);
        }
    }
}

/// Run fastp on a single sample. Returns paths to the trimmed FASTQs
/// and the JSON report; errors carry enough context (sample_id,
/// stdout+stderr from fastp) that a user can diagnose without a
/// re-run.
pub fn run_fastp(spec: &FastpSpec) -> Result<FastpOutput> {
    let sid = &spec.sample.sample_id;
    let trimmed_r1 = spec.tmpdir.join(format!("{sid}_R1.trimmed.fq.gz"));
    let json_report = spec.tmpdir.join(format!("{sid}.fastp.json"));
    let html_report = spec.tmpdir.join(format!("{sid}.fastp.html"));

    // Common args: threads, json/html output paths.
    let mut args: Vec<String> = vec![
        "--thread".to_string(),
        spec.threads.to_string(),
        "--json".to_string(),
        json_report.to_string_lossy().to_string(),
        "--html".to_string(),
        html_report.to_string_lossy().to_string(),
        // Mark HTML report with the sample_id so users browsing a
        // report directory can tell samples apart at a glance.
        "--report_title".to_string(),
        format!("MycoNote-CLI quant: {sid}"),
    ];

    // Input + output FASTQ args. Salmon will read from the `-o` /
    // `-O` paths after we return.
    args.extend_from_slice(&[
        "--in1".to_string(),
        spec.sample.fastq_r1.to_string_lossy().to_string(),
        "--out1".to_string(),
        trimmed_r1.to_string_lossy().to_string(),
    ]);

    let mut trimmed_r2: Option<PathBuf> = None;
    if let Some(ref r2_in) = spec.sample.fastq_r2 {
        let r2_out = spec.tmpdir.join(format!("{sid}_R2.trimmed.fq.gz"));
        args.extend_from_slice(&[
            "--in2".to_string(),
            r2_in.to_string_lossy().to_string(),
            "--out2".to_string(),
            r2_out.to_string_lossy().to_string(),
        ]);
        trimmed_r2 = Some(r2_out);
    }

    let output = duct::cmd(spec.bin(), &args)
        .stderr_to_stdout()
        .unchecked()
        .run()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "fastp".to_string(),
            message: format!("sample {sid}: spawn failed: {e}"),
        })?;

    if !output.status.success() {
        return Err(MycoNoteError::QuantTool {
            tool: "fastp".to_string(),
            message: format!(
                "sample {sid}: fastp exited {}:\n{}",
                output.status,
                String::from_utf8_lossy(&output.stdout)
            ),
        });
    }

    // Defensive existence check — fastp is conservative about emitting
    // empty outputs on malformed input, but we want a clear error
    // rather than "file not found" arriving from salmon later.
    for (label, p) in [("trimmed_r1", &trimmed_r1), ("fastp.json", &json_report)] {
        if !p.is_file() {
            return Err(MycoNoteError::QuantTool {
                tool: "fastp".to_string(),
                message: format!(
                    "sample {sid}: expected output {label} missing at {}",
                    p.display()
                ),
            });
        }
    }
    if let Some(ref p) = trimmed_r2 {
        if !p.is_file() {
            return Err(MycoNoteError::QuantTool {
                tool: "fastp".to_string(),
                message: format!(
                    "sample {sid}: expected trimmed R2 missing at {}",
                    p.display()
                ),
            });
        }
    }

    // Arm cleanup over the trimmed FASTQs only. JSON + HTML stay; the
    // dispatcher copies them into quant_out/fastp/ after the sample
    // succeeds, and the JSON is later consumed by extract_fastp_summary.
    let mut cleanup_files = vec![trimmed_r1.clone()];
    if let Some(ref r2) = trimmed_r2 {
        cleanup_files.push(r2.clone());
    }

    Ok(FastpOutput {
        trimmed_r1,
        trimmed_r2,
        json_report,
        html_report,
        _cleanup: TrimmedCleanup {
            files: cleanup_files,
            armed: !spec.persist_trimmed,
        },
    })
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quant::sample_sheet::Strandedness;
    use std::collections::BTreeMap;
    use std::io::Write;
    use tempfile::TempDir;

    /// Write a tiny gzipped FASTQ with N reads. fastp reads gzipped
    /// input natively, so we don't need plain-text FASTQ test fixtures.
    fn write_tiny_fq_gz(path: &Path, n_reads: usize, read_name_prefix: &str) {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let f = std::fs::File::create(path).unwrap();
        let mut enc = GzEncoder::new(f, Compression::fast());
        for i in 0..n_reads {
            // 50 bp of the same 4-mer repeated — fastp doesn't care
            // about sequence content, only about format validity.
            writeln!(enc, "@{read_name_prefix}_{i}").unwrap();
            writeln!(enc, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC").unwrap();
            writeln!(enc, "+").unwrap();
            writeln!(enc, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();
        }
        enc.finish().unwrap();
    }

    fn make_sample(id: &str, r1: PathBuf, r2: Option<PathBuf>) -> Sample {
        Sample {
            sample_id: id.to_string(),
            fastq_r1: r1,
            fastq_r2: r2,
            condition: None,
            strandedness: Strandedness::Auto,
            batch: None,
            source_line: 2,
            extras: BTreeMap::new(),
        }
    }

    #[test]
    #[ignore]
    fn fastp_live_paired_end_smoke() {
        if which::which("fastp").is_err() {
            eprintln!("fastp not on PATH — skipping");
            return;
        }

        let tmp = TempDir::new().unwrap();
        let r1 = tmp.path().join("in_R1.fq.gz");
        let r2 = tmp.path().join("in_R2.fq.gz");
        write_tiny_fq_gz(&r1, 100, "r1");
        write_tiny_fq_gz(&r2, 100, "r2");

        let sample = make_sample("smoke_PE", r1, Some(r2));
        let spec = FastpSpec {
            sample: &sample,
            tmpdir: tmp.path(),
            threads: 2,
            fastp_bin: None,
            persist_trimmed: true, // keep files so we can inspect them
        };

        let out = run_fastp(&spec).unwrap();
        assert!(out.trimmed_r1.is_file(), "trimmed R1 missing");
        assert!(
            out.trimmed_r2.as_ref().unwrap().is_file(),
            "trimmed R2 missing"
        );
        assert!(out.json_report.is_file(), "JSON report missing");
        // JSON should parse and contain the top-level 'summary' key.
        let json = std::fs::read_to_string(&out.json_report).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v.get("summary").is_some());
        assert!(v.get("filtering_result").is_some());

        // Confirm our bundle extractor copes with this live JSON.
        let summary = crate::quant::bundle::extract_fastp_summary(&out.json_report).unwrap();
        assert!(summary.reads_before_filtering > 0);
        assert!(summary.reads_after_filtering > 0);
    }

    #[test]
    #[ignore]
    fn fastp_live_single_end_smoke() {
        if which::which("fastp").is_err() {
            eprintln!("fastp not on PATH — skipping");
            return;
        }

        let tmp = TempDir::new().unwrap();
        let r1 = tmp.path().join("single.fq.gz");
        write_tiny_fq_gz(&r1, 50, "se");

        let sample = make_sample("smoke_SE", r1, None);
        let spec = FastpSpec {
            sample: &sample,
            tmpdir: tmp.path(),
            threads: 2,
            fastp_bin: None,
            persist_trimmed: true,
        };

        let out = run_fastp(&spec).unwrap();
        assert!(out.trimmed_r1.is_file());
        assert!(out.trimmed_r2.is_none());
        assert!(out.json_report.is_file());
    }

    #[test]
    #[ignore]
    fn fastp_drop_guard_removes_trimmed_files() {
        if which::which("fastp").is_err() {
            eprintln!("fastp not on PATH — skipping");
            return;
        }

        let tmp = TempDir::new().unwrap();
        let r1 = tmp.path().join("in_R1.fq.gz");
        write_tiny_fq_gz(&r1, 20, "r1");
        let sample = make_sample("drop_test", r1, None);

        let (saved_path, saved_json_path) = {
            let spec = FastpSpec {
                sample: &sample,
                tmpdir: tmp.path(),
                threads: 1,
                fastp_bin: None,
                persist_trimmed: false, // guard armed; should clean on drop
            };
            let out = run_fastp(&spec).unwrap();
            (out.trimmed_r1.clone(), out.json_report.clone())
            // `out` drops here; cleanup should fire
        };

        assert!(
            !saved_path.exists(),
            "trimmed file should have been removed by drop guard, still at {}",
            saved_path.display()
        );
        // JSON report is intentionally NOT cleaned up — the
        // dispatcher reads it after fastp returns.
        assert!(
            saved_json_path.exists(),
            "JSON report should NOT be removed by drop guard"
        );
    }

    #[test]
    #[ignore]
    fn fastp_forget_cleanup_persists_trimmed() {
        if which::which("fastp").is_err() {
            eprintln!("fastp not on PATH — skipping");
            return;
        }

        let tmp = TempDir::new().unwrap();
        let r1 = tmp.path().join("in_R1.fq.gz");
        write_tiny_fq_gz(&r1, 20, "r1");
        let sample = make_sample("forget_test", r1, None);

        let saved_path = {
            let spec = FastpSpec {
                sample: &sample,
                tmpdir: tmp.path(),
                threads: 1,
                fastp_bin: None,
                persist_trimmed: false,
            };
            let out = run_fastp(&spec).unwrap().forget_cleanup();
            out.trimmed_r1.clone()
        };

        assert!(
            saved_path.exists(),
            "forget_cleanup() should prevent deletion, file at {}",
            saved_path.display()
        );
    }
}
