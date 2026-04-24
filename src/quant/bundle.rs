//! Reproducibility manifest for a `quant` run.
//!
//! Writes `quant_out/quant_bundle.json` — a fixed-schema JSON summary
//! that records exactly what ran: external-tool versions, input file
//! SHA256s, salmon index hash + k + target/decoy counts, and per-sample
//! QC stats pulled from each fastp JSON. The full fastp JSON per
//! sample is persisted under `quant_out/fastp/{sample_id}.json` —
//! this bundle keeps the one-screen summary while the long-form
//! artifact stays one `ls` away.
//!
//! The JSON-path mapping from fastp 1.3.2's output to our per-sample
//! summary was empirically verified against a live fastp run. Field
//! locations:
//!
//!   reads_before_filtering   ← summary.before_filtering.total_reads
//!   reads_after_filtering    ← summary.after_filtering.total_reads
//!   q30_rate_before          ← summary.before_filtering.q30_rate
//!   q30_rate_after           ← summary.after_filtering.q30_rate
//!   reads_passed_pct         derived from filtering_result.passed_filter_reads
//!   adapter_trimmed_reads    ← adapter_cutting.adapter_trimmed_reads (Option)
//!   adapter_trimmed_bases    ← adapter_cutting.adapter_trimmed_bases (Option)
//!   duplication_rate         ← duplication.rate
//!   insert_size_peak         ← insert_size.peak (Option; absent for SE)
//!
//! The `adapter_cutting` section is missing entirely when fastp runs
//! with `--disable_adapter_trimming`, which is why those two fields
//! are `Option<u64>` here.

use crate::utils::error::{MycoNoteError, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// Top-level manifest structure
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Serialize, Deserialize)]
pub struct QuantBundle {
    pub version: String,
    pub run_timestamp: String,
    pub command: String,
    pub external_tools: ExternalTools,
    pub inputs: BundleInputs,
    pub index: BundleIndex,
    pub samples: Vec<SampleSummary>,
    pub seed: u64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ExternalTools {
    pub fastp: String,
    pub salmon: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BundleInputs {
    pub transcripts: HashedPath,
    pub genome: HashedPath,
    pub sample_sheet: HashedPath,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct HashedPath {
    pub path: String,
    pub sha256: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BundleIndex {
    pub path: String,
    pub sha256: String,
    pub k: usize,
    pub decoys: usize,
    pub targets: usize,
}

/// Per-sample QC + quant summary. Salmon-derived fields
/// (`mapping_rate`, `library_size`) come from the sample's
/// `logs/salmon_quant.log` + `quant.sf`; fastp-derived fields come
/// from `fastp.json`. The dispatcher assembles this struct after both
/// tools finish for the sample.
#[derive(Debug, Serialize, Deserialize)]
pub struct SampleSummary {
    pub sample_id: String,
    pub mapping_rate: f64,
    pub library_size: u64,
    pub reads_before_filtering: u64,
    pub reads_after_filtering: u64,
    pub reads_passed_pct: f64,
    pub q30_rate_before: f64,
    pub q30_rate_after: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub adapter_trimmed_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub adapter_trimmed_bases: Option<u64>,
    pub duplication_rate: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub insert_size_peak: Option<u64>,
}

impl QuantBundle {
    /// Serialize the bundle as pretty JSON and write to `path`.
    pub fn write(&self, path: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(self)?;
        let mut f = File::create(path)?;
        f.write_all(json.as_bytes())?;
        f.write_all(b"\n")?;
        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers — file hashing + fastp JSON extraction
// ─────────────────────────────────────────────────────────────────────────────

/// Streaming SHA256 of a file. Used to stamp input provenance into
/// the bundle so reruns against the "same" genome.fa that has been
/// silently edited can be detected.
pub fn sha256_file(path: &Path) -> Result<String> {
    let f = File::open(path)?;
    let mut reader = BufReader::new(f);
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 8192];
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

/// The summary fields we lift out of a fastp `--json` report. Kept
/// separate from `SampleSummary` because mapping_rate + library_size
/// come from salmon, not fastp — the dispatcher combines them.
#[derive(Debug, Clone, PartialEq)]
pub struct FastpSummary {
    pub reads_before_filtering: u64,
    pub reads_after_filtering: u64,
    pub reads_passed_pct: f64,
    pub q30_rate_before: f64,
    pub q30_rate_after: f64,
    pub adapter_trimmed_reads: Option<u64>,
    pub adapter_trimmed_bases: Option<u64>,
    pub duplication_rate: f64,
    pub insert_size_peak: Option<u64>,
}

/// Parse a fastp JSON report and extract the bundle-relevant fields.
/// Returns `QuantTool { tool: "fastp", .. }` on any structural issue
/// so the error message names which sample's JSON was malformed.
pub fn extract_fastp_summary(json_path: &Path) -> Result<FastpSummary> {
    let text = std::fs::read_to_string(json_path).map_err(|e| MycoNoteError::QuantTool {
        tool: "fastp".to_string(),
        message: format!("read {}: {}", json_path.display(), e),
    })?;
    parse_fastp_summary(&text).map_err(|msg| MycoNoteError::QuantTool {
        tool: "fastp".to_string(),
        message: format!("{}: {}", json_path.display(), msg),
    })
}

fn parse_fastp_summary(text: &str) -> std::result::Result<FastpSummary, String> {
    let v: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("malformed JSON: {e}"))?;

    let get_u64 = |path: &[&str]| -> std::result::Result<u64, String> {
        let mut cur = &v;
        for key in path {
            cur = cur
                .get(*key)
                .ok_or_else(|| format!("missing field '{}' (path {:?})", key, path))?;
        }
        cur.as_u64()
            .ok_or_else(|| format!("field {:?} is not u64: {:?}", path, cur))
    };
    let get_f64 = |path: &[&str]| -> std::result::Result<f64, String> {
        let mut cur = &v;
        for key in path {
            cur = cur
                .get(*key)
                .ok_or_else(|| format!("missing field '{}' (path {:?})", key, path))?;
        }
        cur.as_f64()
            .ok_or_else(|| format!("field {:?} is not f64: {:?}", path, cur))
    };
    let try_u64 = |path: &[&str]| -> Option<u64> {
        let mut cur = &v;
        for key in path {
            cur = cur.get(*key)?;
        }
        cur.as_u64()
    };

    let reads_before = get_u64(&["summary", "before_filtering", "total_reads"])?;
    let reads_after = get_u64(&["summary", "after_filtering", "total_reads"])?;
    let q30_before = get_f64(&["summary", "before_filtering", "q30_rate"])?;
    let q30_after = get_f64(&["summary", "after_filtering", "q30_rate"])?;
    let passed = get_u64(&["filtering_result", "passed_filter_reads"])?;
    let dup_rate = get_f64(&["duplication", "rate"])?;

    // Derived percentage. Guard against fastp emitting 0 total_reads
    // on an empty input — division by zero would produce NaN and
    // serde_json::to_string_pretty would then emit `null`, breaking
    // round-trip parsing of the bundle.
    let passed_pct = if reads_before == 0 {
        0.0
    } else {
        passed as f64 / reads_before as f64 * 100.0
    };

    // Optional sections. adapter_cutting is absent when fastp ran
    // with --disable_adapter_trimming. insert_size.peak is 0 for
    // single-end — treat 0 as None so the bundle doesn't claim
    // an insert size for unpaired reads.
    let adapter_reads = try_u64(&["adapter_cutting", "adapter_trimmed_reads"]);
    let adapter_bases = try_u64(&["adapter_cutting", "adapter_trimmed_bases"]);
    let insert_peak = match try_u64(&["insert_size", "peak"]) {
        Some(0) | None => None,
        Some(n) => Some(n),
    };

    Ok(FastpSummary {
        reads_before_filtering: reads_before,
        reads_after_filtering: reads_after,
        reads_passed_pct: passed_pct,
        q30_rate_before: q30_before,
        q30_rate_after: q30_after,
        adapter_trimmed_reads: adapter_reads,
        adapter_trimmed_bases: adapter_bases,
        duplication_rate: dup_rate,
        insert_size_peak: insert_peak,
    })
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    // Trimmed real fastp 1.3.2 JSON — paired-end, adapter-cutting
    // present, nonzero inputs. Captured from a live run during the
    // 0.3.0 design phase; kept inline rather than as a fixture file
    // to avoid committing generated artifacts.
    const FASTP_PAIRED_JSON: &str = r#"{
      "summary": {
        "fastp_version": "1.3.2",
        "sequencing": "paired end (53 cycles + 53 cycles)",
        "before_filtering": {
          "total_reads": 1000,
          "total_bases": 53000,
          "q20_bases": 52000,
          "q30_bases": 51500,
          "q20_rate": 0.981132,
          "q30_rate": 0.971698,
          "read1_mean_length": 53,
          "read2_mean_length": 53,
          "gc_content": 0.45
        },
        "after_filtering": {
          "total_reads": 980,
          "total_bases": 51940,
          "q20_bases": 51000,
          "q30_bases": 50600,
          "q20_rate": 0.981902,
          "q30_rate": 0.974200,
          "read1_mean_length": 53,
          "read2_mean_length": 53,
          "gc_content": 0.4502
        }
      },
      "filtering_result": {
        "passed_filter_reads": 980,
        "low_quality_reads": 12,
        "too_many_N_reads": 0,
        "too_short_reads": 8,
        "too_long_reads": 0
      },
      "duplication": { "rate": 0.02 },
      "insert_size": { "peak": 180, "unknown": 5, "histogram": [] },
      "adapter_cutting": {
        "adapter_trimmed_reads": 40,
        "adapter_trimmed_bases": 520
      }
    }"#;

    const FASTP_SE_NO_ADAPTER_JSON: &str = r#"{
      "summary": {
        "before_filtering": {
          "total_reads": 500,
          "q30_rate": 0.95
        },
        "after_filtering": {
          "total_reads": 495,
          "q30_rate": 0.96
        }
      },
      "filtering_result": { "passed_filter_reads": 495 },
      "duplication": { "rate": 0.01 },
      "insert_size": { "peak": 0 }
    }"#;

    #[test]
    fn extract_paired_fastp() {
        let s = parse_fastp_summary(FASTP_PAIRED_JSON).unwrap();
        assert_eq!(s.reads_before_filtering, 1000);
        assert_eq!(s.reads_after_filtering, 980);
        assert!((s.reads_passed_pct - 98.0).abs() < 1e-9);
        assert!((s.q30_rate_before - 0.971698).abs() < 1e-6);
        assert!((s.q30_rate_after - 0.974200).abs() < 1e-6);
        assert_eq!(s.adapter_trimmed_reads, Some(40));
        assert_eq!(s.adapter_trimmed_bases, Some(520));
        assert!((s.duplication_rate - 0.02).abs() < 1e-9);
        assert_eq!(s.insert_size_peak, Some(180));
    }

    #[test]
    fn extract_single_end_missing_adapter_section() {
        let s = parse_fastp_summary(FASTP_SE_NO_ADAPTER_JSON).unwrap();
        assert_eq!(s.reads_before_filtering, 500);
        assert_eq!(s.reads_after_filtering, 495);
        assert_eq!(s.adapter_trimmed_reads, None);
        assert_eq!(s.adapter_trimmed_bases, None);
        assert_eq!(s.insert_size_peak, None, "peak=0 must collapse to None");
    }

    #[test]
    fn reject_malformed_fastp_json() {
        let err = parse_fastp_summary("not json").unwrap_err();
        assert!(err.contains("malformed JSON"));
    }

    #[test]
    fn reject_fastp_missing_required_field() {
        let missing = r#"{"summary": {"before_filtering": {}}}"#;
        let err = parse_fastp_summary(missing).unwrap_err();
        assert!(err.contains("missing field"), "got: {err}");
    }

    #[test]
    fn sha256_file_matches_known_value() {
        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("hello.txt");
        std::fs::write(&p, b"hello").unwrap();
        let h = sha256_file(&p).unwrap();
        // SHA256("hello") — precomputed, locks our hash function in
        // place against accidental truncation / encoding bugs.
        assert_eq!(
            h,
            "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
        );
    }

    #[test]
    fn sha256_file_matches_on_longer_than_buffer() {
        // Default read buffer is 8 KiB — exercise the multi-chunk path.
        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("big.bin");
        let payload = vec![0u8; 20_000];
        std::fs::write(&p, &payload).unwrap();
        let h = sha256_file(&p).unwrap();
        // SHA256 of 20000 zero bytes — verified via
        // `python3 -c 'import hashlib; print(hashlib.sha256(b"\0" * 20000).hexdigest())'`
        assert_eq!(
            h,
            "28b4f41a7f3ee6d8cc87272db6e09c6d3566551fd4d18702b041a21658272a85"
        );
    }

    #[test]
    fn bundle_roundtrip_serde() {
        let b = QuantBundle {
            version: "0.3.0".to_string(),
            run_timestamp: "2026-04-23T17:30:00Z".to_string(),
            command: "myconote-cli quant cds.fa --samples s.tsv".to_string(),
            external_tools: ExternalTools {
                fastp: "1.3.2".to_string(),
                salmon: "1.11.4".to_string(),
            },
            inputs: BundleInputs {
                transcripts: HashedPath {
                    path: "cds.fa".to_string(),
                    sha256: "a".repeat(64),
                },
                genome: HashedPath {
                    path: "genome.fa".to_string(),
                    sha256: "b".repeat(64),
                },
                sample_sheet: HashedPath {
                    path: "s.tsv".to_string(),
                    sha256: "c".repeat(64),
                },
            },
            index: BundleIndex {
                path: "/cache/salmon_index/xyz".to_string(),
                sha256: "d".repeat(64),
                k: 31,
                decoys: 12,
                targets: 10547,
            },
            samples: vec![SampleSummary {
                sample_id: "WT_rep1".to_string(),
                mapping_rate: 0.89,
                library_size: 12_345_678,
                reads_before_filtering: 10_000_000,
                reads_after_filtering: 9_800_000,
                reads_passed_pct: 98.0,
                q30_rate_before: 0.95,
                q30_rate_after: 0.97,
                adapter_trimmed_reads: Some(1_000),
                adapter_trimmed_bases: Some(20_000),
                duplication_rate: 0.05,
                insert_size_peak: Some(200),
            }],
            seed: 42,
        };

        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("quant_bundle.json");
        b.write(&p).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let parsed: QuantBundle = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed.version, "0.3.0");
        assert_eq!(parsed.samples.len(), 1);
        assert_eq!(parsed.samples[0].sample_id, "WT_rep1");
        assert_eq!(parsed.samples[0].adapter_trimmed_reads, Some(1_000));
    }

    #[test]
    fn bundle_omits_none_adapter_fields() {
        // When adapter_trimmed_{reads,bases} is None, the serialized
        // JSON should not include those keys — `#[serde(skip_serializing_if)]`
        // controls this. Tested so a future refactor can't silently
        // regress to emitting `null`, which confuses downstream parsers.
        let b = QuantBundle {
            version: "0.3.0".to_string(),
            run_timestamp: "2026-04-23T17:30:00Z".to_string(),
            command: "x".to_string(),
            external_tools: ExternalTools {
                fastp: "1.3.2".to_string(),
                salmon: "1.11.4".to_string(),
            },
            inputs: BundleInputs {
                transcripts: HashedPath {
                    path: "t".to_string(),
                    sha256: "0".repeat(64),
                },
                genome: HashedPath {
                    path: "g".to_string(),
                    sha256: "0".repeat(64),
                },
                sample_sheet: HashedPath {
                    path: "s".to_string(),
                    sha256: "0".repeat(64),
                },
            },
            index: BundleIndex {
                path: "i".to_string(),
                sha256: "0".repeat(64),
                k: 31,
                decoys: 0,
                targets: 0,
            },
            samples: vec![SampleSummary {
                sample_id: "SE".to_string(),
                mapping_rate: 0.5,
                library_size: 100,
                reads_before_filtering: 100,
                reads_after_filtering: 100,
                reads_passed_pct: 100.0,
                q30_rate_before: 0.9,
                q30_rate_after: 0.9,
                adapter_trimmed_reads: None,
                adapter_trimmed_bases: None,
                duplication_rate: 0.0,
                insert_size_peak: None,
            }],
            seed: 0,
        };

        let json = serde_json::to_string_pretty(&b).unwrap();
        assert!(
            !json.contains("adapter_trimmed_reads"),
            "None field should not be serialized: {json}"
        );
        assert!(
            !json.contains("insert_size_peak"),
            "None field should not be serialized: {json}"
        );
    }
}
