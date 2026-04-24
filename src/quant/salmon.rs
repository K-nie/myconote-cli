//! salmon quant wrapper — map + quantify one sample against the
//! pre-built decoy-aware index.
//!
//! Spec defaults (`scratch/rnaseq_spec.md` §3): libType from the
//! sample sheet's `strandedness` column (Auto when unset),
//! `--validateMappings`, `--gcBias`, `--seqBias`, `--threads N`.
//! Fragment-length distribution, equivalence classes, and bootstrap
//! replicates stay at salmon defaults — the bundle records the full
//! command so a user who wants bootstraps can re-run with
//! `--salmon-extra "--numBootstraps 30"` once that flag lands (not
//! in 0.3.0).
//!
//! Output tree per sample (all under `output_dir/<sample_id>/`):
//!   - `quant.sf`                 tximport-native format
//!   - `lib_format_counts.json`   compatible library-type breakdown
//!   - `logs/salmon_quant.log`    contains the "Mapping rate" line
//!                                 we lift into the bundle
//!   - plus other salmon artifacts we don't parse here
//!
//! The caller (dispatcher) is responsible for reading `quant.sf` via
//! `merge::parse_quant_sf`. This module is a thin wrapper around the
//! subprocess boundary.

use crate::quant::sample_sheet::Sample;
use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};

pub struct SalmonQuantSpec<'a> {
    pub sample: &'a Sample,
    pub index_path: &'a Path,
    pub trimmed_r1: &'a Path,
    pub trimmed_r2: Option<&'a Path>,
    /// Root output directory for this sample's salmon run. The
    /// wrapper does NOT create a per-sample sub-directory — the
    /// caller is responsible for picking
    /// `quant_out/salmon/<sample_id>/` (or equivalent) and passing
    /// that path here so the sample_id does not re-appear twice in
    /// the path.
    pub output_dir: &'a Path,
    pub threads: usize,
    pub salmon_bin: Option<&'a str>,
}

impl<'a> SalmonQuantSpec<'a> {
    fn bin(&self) -> &str {
        self.salmon_bin.unwrap_or("salmon")
    }
}

#[derive(Debug)]
pub struct SalmonQuantOutput {
    pub quant_sf: PathBuf,
    pub mapping_rate: f64,
    pub library_size: u64,
}

/// Run `salmon quant` for one sample. Parses the mapping rate out of
/// salmon's log and the library size out of `quant.sf` so the
/// dispatcher can populate the bundle's per-sample summary.
pub fn run_salmon_quant(spec: &SalmonQuantSpec) -> Result<SalmonQuantOutput> {
    std::fs::create_dir_all(spec.output_dir)?;

    let paired = spec.trimmed_r2.is_some();
    let libtype = spec.sample.strandedness.salmon_libtype(paired);

    let mut args: Vec<String> = vec![
        "quant".to_string(),
        "--index".to_string(),
        spec.index_path.to_string_lossy().to_string(),
        "--libType".to_string(),
        libtype.to_string(),
        "--validateMappings".to_string(),
        "--gcBias".to_string(),
        "--seqBias".to_string(),
        "--threads".to_string(),
        spec.threads.to_string(),
        "-o".to_string(),
        spec.output_dir.to_string_lossy().to_string(),
    ];
    if let Some(r2) = spec.trimmed_r2 {
        args.extend_from_slice(&[
            "-1".to_string(),
            spec.trimmed_r1.to_string_lossy().to_string(),
            "-2".to_string(),
            r2.to_string_lossy().to_string(),
        ]);
    } else {
        args.extend_from_slice(&[
            "-r".to_string(),
            spec.trimmed_r1.to_string_lossy().to_string(),
        ]);
    }

    let output = duct::cmd(spec.bin(), &args)
        .stderr_to_stdout()
        .unchecked()
        .run()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "salmon".to_string(),
            message: format!("sample {}: salmon spawn failed: {e}", spec.sample.sample_id),
        })?;

    if !output.status.success() {
        return Err(MycoNoteError::QuantTool {
            tool: "salmon".to_string(),
            message: format!(
                "sample {}: salmon quant exited {}:\n{}",
                spec.sample.sample_id,
                output.status,
                String::from_utf8_lossy(&output.stdout)
            ),
        });
    }

    let quant_sf = spec.output_dir.join("quant.sf");
    if !quant_sf.is_file() {
        return Err(MycoNoteError::QuantTool {
            tool: "salmon".to_string(),
            message: format!(
                "sample {}: quant.sf missing after successful salmon run at {}",
                spec.sample.sample_id,
                quant_sf.display()
            ),
        });
    }

    let log_path = spec.output_dir.join("logs").join("salmon_quant.log");
    let mapping_rate = parse_mapping_rate_from_log(&log_path).unwrap_or(0.0);
    let library_size = sum_num_reads_in_quant_sf(&quant_sf)?;

    Ok(SalmonQuantOutput {
        quant_sf,
        mapping_rate,
        library_size,
    })
}

/// Lift the mapping rate from salmon's log. Log line format (salmon
/// 1.10+): `[jointLog] [info] Mapping rate = 85.3821%`. We return the
/// rate as a fraction (0.0–1.0) — the bundle renders it that way.
/// Returns `None` if the log is missing or the line cannot be found;
/// callers fall back to `0.0` which is honest (unknown).
fn parse_mapping_rate_from_log(log_path: &Path) -> Option<f64> {
    let text = std::fs::read_to_string(log_path).ok()?;
    parse_mapping_rate(&text)
}

fn parse_mapping_rate(log_text: &str) -> Option<f64> {
    for line in log_text.lines() {
        // Accept either "Mapping rate = X%" or "Mapping rate: X%".
        // Salmon has shipped both historically.
        let idx_rate = line.find("Mapping rate");
        if idx_rate.is_none() {
            continue;
        }
        // Find the number between "Mapping rate" and the trailing '%'.
        let sep_idx = line[idx_rate.unwrap()..].find(|c: char| c == '=' || c == ':')?;
        let after_sep = &line[idx_rate.unwrap() + sep_idx + 1..];
        let pct_idx = after_sep.find('%')?;
        let num_str = after_sep[..pct_idx].trim();
        let pct: f64 = num_str.parse().ok()?;
        return Some(pct / 100.0);
    }
    None
}

/// Sum the `NumReads` column of `quant.sf`. Represents the assigned
/// (not raw) library size — what DESeq2 will see via tximport.
fn sum_num_reads_in_quant_sf(quant_sf: &Path) -> Result<u64> {
    let text = std::fs::read_to_string(quant_sf).map_err(|e| MycoNoteError::QuantTool {
        tool: "salmon".to_string(),
        message: format!("reading {}: {e}", quant_sf.display()),
    })?;
    let mut total: f64 = 0.0;
    for (i, line) in text.lines().enumerate() {
        if i == 0 {
            continue; // header
        }
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            continue;
        }
        if let Ok(n) = cols[4].parse::<f64>() {
            total += n;
        }
    }
    // Round rather than truncate — salmon's NumReads is already a
    // floating-point estimate, so nearest-integer is the honest
    // representation of library size.
    Ok(total.round() as u64)
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quant::fastp::{run_fastp, FastpSpec};
    use crate::quant::index::{build_or_reuse_index, IndexSpec};
    use crate::quant::sample_sheet::{Sample, Strandedness};
    use std::collections::BTreeMap;
    use std::io::Write;
    use tempfile::TempDir;

    // ── mapping-rate log parser (pure) ───────────────────────────────────────

    #[test]
    fn parse_mapping_rate_equals_form() {
        let log = "[2026-04-23 10:00:00.000] [jointLog] [info] Mapping rate = 85.3821%\n";
        let r = parse_mapping_rate(log).unwrap();
        assert!((r - 0.853821).abs() < 1e-9);
    }

    #[test]
    fn parse_mapping_rate_colon_form() {
        let log = "[info] Mapping rate: 50.0%\n";
        assert_eq!(parse_mapping_rate(log), Some(0.5));
    }

    #[test]
    fn parse_mapping_rate_zero_percent() {
        let log = "Mapping rate = 0%\n";
        assert_eq!(parse_mapping_rate(log), Some(0.0));
    }

    #[test]
    fn parse_mapping_rate_missing_returns_none() {
        let log = "[info] Indexing done.\n[info] Done.\n";
        assert_eq!(parse_mapping_rate(log), None);
    }

    #[test]
    fn parse_mapping_rate_handles_surrounding_lines() {
        let log = "[info] Start\n\
                   [info] Mapping rate = 12.5%\n\
                   [info] End\n";
        assert_eq!(parse_mapping_rate(log), Some(0.125));
    }

    // ── NumReads sum (pure) ──────────────────────────────────────────────────

    #[test]
    fn sum_num_reads_rounds_fractional() {
        let tmp = TempDir::new().unwrap();
        let sf = tmp.path().join("quant.sf");
        std::fs::write(
            &sf,
            "Name\tLength\tEffectiveLength\tTPM\tNumReads\n\
             t1\t100\t80\t10\t123.4\n\
             t2\t100\t80\t20\t456.8\n",
        )
        .unwrap();
        assert_eq!(sum_num_reads_in_quant_sf(&sf).unwrap(), 580); // 123.4 + 456.8 = 580.2 → 580
    }

    #[test]
    fn sum_num_reads_empty_sf_is_zero() {
        let tmp = TempDir::new().unwrap();
        let sf = tmp.path().join("quant.sf");
        std::fs::write(&sf, "Name\tLength\tEffectiveLength\tTPM\tNumReads\n").unwrap();
        assert_eq!(sum_num_reads_in_quant_sf(&sf).unwrap(), 0);
    }

    // ── live end-to-end smoke: index → fastp → salmon quant ─────────────────
    //
    // Verifies the full subprocess chain works against salmon 1.11.4
    // and fastp 1.3.2. Uses synthetic FASTQ reads drawn DIRECTLY
    // from the indexed transcripts, so mapping rate should be >0
    // even on this tiny fixture.

    fn write_tiny_fq_gz(path: &Path, reads: &[(&str, &str)]) {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let f = std::fs::File::create(path).unwrap();
        let mut enc = GzEncoder::new(f, Compression::fast());
        for (name, seq) in reads {
            writeln!(enc, "@{name}").unwrap();
            writeln!(enc, "{seq}").unwrap();
            writeln!(enc, "+").unwrap();
            writeln!(enc, "{}", "I".repeat(seq.len())).unwrap();
        }
        enc.finish().unwrap();
    }

    fn make_sample(id: &str, r1: PathBuf, r2: Option<PathBuf>) -> Sample {
        Sample {
            sample_id: id.to_string(),
            fastq_r1: r1,
            fastq_r2: r2,
            condition: None,
            // Explicit unstranded — auto-detection needs more reads
            // than our synthetic fixture provides.
            strandedness: Strandedness::Unstranded,
            batch: None,
            source_line: 2,
            extras: BTreeMap::new(),
        }
    }

    #[test]
    #[ignore]
    fn salmon_live_end_to_end_single_end() {
        if which::which("salmon").is_err() || which::which("fastp").is_err() {
            eprintln!("salmon or fastp missing from PATH — skipping");
            return;
        }

        let tmp = TempDir::new().unwrap();

        // Transcripts: two 60-bp sequences. k=21 → salmon needs each
        // transcript ≥ k bp; 60 gives comfortable margin.
        let t1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let t2 = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG";
        let cds = tmp.path().join("cds.fa");
        std::fs::write(&cds, format!(">t1\n{t1}\n>t2\n{t2}\n")).unwrap();

        // Decoy genome: longer so salmon is happy with the decoy set.
        // Content is independent from the CDS so reads don't match it.
        let genome = tmp.path().join("genome.fa");
        std::fs::write(
            &genome,
            ">chr1\n\
             GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC\n\
             >chr2\n\
             TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT\n",
        )
        .unwrap();

        // Build the index.
        let cache = tmp.path().join("cache");
        let ispec = IndexSpec {
            cds_fa: cds,
            genome_fa: genome,
            k: 21,
            threads: 2,
            salmon_bin: None,
        };
        let idx = build_or_reuse_index(&ispec, &cache).expect("index build");

        // Build a tiny single-end FASTQ of reads drawn from t1 + t2.
        // 30 bp reads → 10 reads each → 20 reads total.
        let mut reads = Vec::new();
        for i in 0..10 {
            let start = i * 3;
            let end = start + 30;
            if end <= t1.len() {
                reads.push((format!("r1_{i}"), &t1[start..end]));
            }
        }
        for i in 0..10 {
            let start = i * 3;
            let end = start + 30;
            if end <= t2.len() {
                reads.push((format!("r2_{i}"), &t2[start..end]));
            }
        }
        let reads_slices: Vec<(&str, &str)> = reads.iter().map(|(n, s)| (n.as_str(), *s)).collect();
        let r1 = tmp.path().join("reads.fq.gz");
        write_tiny_fq_gz(&r1, &reads_slices);

        // Trim with fastp.
        let sample = make_sample("live_SE", r1, None);
        let trim_dir = tmp.path().join("trim");
        std::fs::create_dir_all(&trim_dir).unwrap();
        let fspec = FastpSpec {
            sample: &sample,
            tmpdir: &trim_dir,
            threads: 2,
            fastp_bin: None,
            persist_trimmed: true,
        };
        let fastp_out = run_fastp(&fspec).expect("fastp ran");

        // Quantify.
        let quant_dir = tmp.path().join("quant_out").join(&sample.sample_id);
        let sspec = SalmonQuantSpec {
            sample: &sample,
            index_path: &idx.path,
            trimmed_r1: &fastp_out.trimmed_r1,
            trimmed_r2: None,
            output_dir: &quant_dir,
            threads: 2,
            salmon_bin: None,
        };
        let qout = run_salmon_quant(&sspec).expect("salmon quant ran");

        assert!(qout.quant_sf.is_file(), "quant.sf missing");
        assert!(
            qout.library_size > 0,
            "expected at least some reads to map, got library_size={}",
            qout.library_size
        );
        // Mapping rate should be non-trivial on this fixture — reads
        // are exact substrings of the indexed transcripts. Relaxed
        // bound (>30%) guards against salmon heuristics occasionally
        // rejecting short reads.
        assert!(
            qout.mapping_rate > 0.3,
            "mapping rate too low: {}",
            qout.mapping_rate
        );
    }
}
