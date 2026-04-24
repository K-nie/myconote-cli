//! Per-sample × per-haplotype salmon driver for ASE.
//!
//! For each sample in the sheet:
//!   1. Run `fastp` once (haplotype-independent QC/trim).
//!   2. For each haplotype: run `salmon quant` against that
//!      haplotype's index. Two salmon runs per sample.
//!
//! Salmon indices per haplotype are built (or reused from cache) by
//! the existing `src/quant/index.rs` infrastructure. Each haplotype's
//! personalized CDS FASTA has a different SHA256 → different cache
//! slot automatically; we don't have to thread haplotype identity
//! through the cache key.
//!
//! Mapping-rate asymmetry detection: if the two haplotypes' mapping
//! rates for the same sample differ by more than the configured
//! threshold (default 5 %), we flag the sample. Strong asymmetry
//! usually signals a phasing error or an assembly issue, not real
//! biology.

use crate::quant::fastp::{run_fastp, FastpOutput, FastpSpec};
use crate::quant::index::{build_or_reuse_index, IndexResult, IndexSpec};
use crate::quant::salmon::{run_salmon_quant, SalmonQuantOutput, SalmonQuantSpec};
use crate::quant::sample_sheet::Sample;
use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration for the per-sample-per-haplotype run
// ─────────────────────────────────────────────────────────────────────────────

/// Inputs needed to quantify all samples against all haplotypes.
/// The `AseQuantSpec` is assembled once by the dispatcher (`mod.rs`)
/// and then consumed for each sample.
pub struct AseQuantSpec<'a> {
    /// Per-haplotype CDS FASTA paths on disk (emitted by
    /// `src/ase/personalize.rs` + `write_cds_fasta`). The
    /// index is built per haplotype.
    pub haplotype_cds: &'a [(String, PathBuf)],
    /// Shared genome FASTA for the decoy set — the same for all
    /// haplotypes. Point variants don't change which contigs serve
    /// as decoys.
    pub genome_fa: &'a Path,
    /// Output directory root (typically `ase_out/`).
    pub output_dir: &'a Path,
    /// Index cache root (resolved by the CLI; honors
    /// `--index-cache` / `MYCONOTE_INDEX_CACHE` / XDG fallbacks via
    /// the caller).
    pub index_cache_root: &'a Path,
    /// salmon `-k` value. Default 31 per the existing `quant` spec.
    pub k: usize,
    /// Per-sample salmon / fastp thread count.
    pub threads: usize,
    /// Temp-dir root for fastp-trimmed FASTQs (typically `$TMPDIR`
    /// or a user-supplied `--tmpdir`).
    pub tmpdir: &'a Path,
    /// Optional override for the fastp binary (falls back to PATH).
    pub fastp_bin: Option<&'a str>,
    /// Optional override for the salmon binary.
    pub salmon_bin: Option<&'a str>,
    /// Mapping-rate asymmetry threshold (fraction). Any sample
    /// whose per-haplotype mapping rates differ by more than this
    /// is flagged. Default 0.05 (5 %).
    pub asymmetry_threshold: f64,
}

/// One salmon run's result plus its haplotype label. The dispatcher
/// collects two of these per sample and writes them into the bundle.
#[derive(Debug)]
pub struct SampleHapResult {
    pub haplotype: String,
    pub quant: SalmonQuantOutput,
}

/// Everything we produced for a single sample: fastp output (shared
/// across haplotypes) plus one `SampleHapResult` per haplotype.
pub struct SampleAseResult {
    pub sample_id: String,
    pub fastp: FastpOutput,
    pub hap_results: Vec<SampleHapResult>,
    /// `true` when |mapping_rate_hap0 − mapping_rate_hap1| > threshold.
    pub asymmetry_flag: bool,
}

// ─────────────────────────────────────────────────────────────────────────────
// Index build — once per haplotype, reuses the quant index cache
// ─────────────────────────────────────────────────────────────────────────────

/// Build (or reuse) a salmon index for each haplotype. Each call
/// goes through `src/quant/index.rs::build_or_reuse_index`, so:
///   - identical personalized CDS + genome + salmon version + k
///     → same cache slot across runs
///   - different haplotypes → different cache slots automatically
///     (because each haplotype's CDS FASTA hashes differently)
///
/// Returns `(haplotype_name, IndexResult)` pairs in the same order
/// as the input slice.
pub fn build_haplotype_indices(spec: &AseQuantSpec) -> Result<Vec<(String, IndexResult)>> {
    let mut out = Vec::with_capacity(spec.haplotype_cds.len());
    for (hap_name, cds_fa) in spec.haplotype_cds {
        let ispec = IndexSpec {
            cds_fa: cds_fa.clone(),
            genome_fa: spec.genome_fa.to_path_buf(),
            k: spec.k,
            threads: spec.threads,
            salmon_bin: spec.salmon_bin.map(|s| s.to_string()),
        };
        let result = build_or_reuse_index(&ispec, spec.index_cache_root)?;
        eprintln!(
            "  ✓ haplotype {hap_name}: index at {} ({})",
            result.path.display(),
            if result.cached { "cached" } else { "fresh" }
        );
        out.push((hap_name.clone(), result));
    }
    Ok(out)
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-sample quantification
// ─────────────────────────────────────────────────────────────────────────────

/// Run the full ASE pipeline for one sample: fastp once, salmon
/// twice. Writes outputs into
/// `<output_dir>/salmon/<sample_id>.<hap_name>/`, keeping the naming
/// convention spelled out in `scratch/ase_spec.md`.
pub fn quantify_sample(
    sample: &Sample,
    indices: &[(String, IndexResult)],
    spec: &AseQuantSpec,
) -> Result<SampleAseResult> {
    let sample_tmpdir = spec
        .tmpdir
        .join(format!("myconote_ase_{}", sample.sample_id));
    let _ = std::fs::remove_dir_all(&sample_tmpdir);
    std::fs::create_dir_all(&sample_tmpdir)?;

    // fastp once, reuse the output across both haplotype salmon runs.
    let fastp_spec = FastpSpec {
        sample,
        tmpdir: &sample_tmpdir,
        threads: spec.threads,
        fastp_bin: spec.fastp_bin,
        // `persist_trimmed = true` so the trimmed files survive
        // beyond the FastpOutput's drop guard — the second salmon
        // run needs them. We manually clean up after both finish.
        persist_trimmed: true,
    };
    let fastp_out = run_fastp(&fastp_spec)?;

    let salmon_root = spec.output_dir.join("salmon");
    std::fs::create_dir_all(&salmon_root)?;

    let mut hap_results: Vec<SampleHapResult> = Vec::with_capacity(indices.len());
    for (hap_name, idx) in indices {
        let out_dir = salmon_root.join(format!("{}.{}", sample.sample_id, hap_name));
        let spec_q = SalmonQuantSpec {
            sample,
            index_path: &idx.path,
            trimmed_r1: &fastp_out.trimmed_r1,
            trimmed_r2: fastp_out.trimmed_r2.as_deref(),
            output_dir: &out_dir,
            threads: spec.threads,
            salmon_bin: spec.salmon_bin,
        };
        let qout = run_salmon_quant(&spec_q)?;
        eprintln!(
            "  ✓ {}.{}: mapping_rate={:.3}, library_size={}",
            sample.sample_id, hap_name, qout.mapping_rate, qout.library_size
        );
        hap_results.push(SampleHapResult {
            haplotype: hap_name.clone(),
            quant: qout,
        });
    }

    // Clean up the trimmed FASTQs now that both salmon runs are
    // done. The JSON reports (fastp_out.json_report + html_report)
    // are NOT touched — callers (bundle.rs / dispatcher) persist
    // them under ase_out/fastp/.
    cleanup_trimmed(&fastp_out);

    let asymmetry_flag = detect_asymmetry(&hap_results, spec.asymmetry_threshold);
    if asymmetry_flag {
        eprintln!(
            "  ⚠  {}: mapping-rate asymmetry between haplotypes exceeds {:.1}%",
            sample.sample_id,
            spec.asymmetry_threshold * 100.0
        );
    }

    Ok(SampleAseResult {
        sample_id: sample.sample_id.clone(),
        fastp: fastp_out,
        hap_results,
        asymmetry_flag,
    })
}

/// Max pairwise absolute difference in mapping rate across
/// haplotypes. For the usual diploid case this is a single
/// subtraction; we write it more generally so triploid / higher-ploidy
/// extensions (0.5.x+) reuse the check unchanged.
fn detect_asymmetry(results: &[SampleHapResult], threshold: f64) -> bool {
    let mut max_diff = 0.0_f64;
    for i in 0..results.len() {
        for j in (i + 1)..results.len() {
            let diff = (results[i].quant.mapping_rate - results[j].quant.mapping_rate).abs();
            if diff > max_diff {
                max_diff = diff;
            }
        }
    }
    max_diff > threshold
}

fn cleanup_trimmed(out: &FastpOutput) {
    let _ = std::fs::remove_file(&out.trimmed_r1);
    if let Some(ref r2) = out.trimmed_r2 {
        let _ = std::fs::remove_file(r2);
    }
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn asymmetry_under_threshold_returns_false() {
        use crate::quant::salmon::SalmonQuantOutput;
        let rs = vec![
            SampleHapResult {
                haplotype: "hap0".to_string(),
                quant: SalmonQuantOutput {
                    quant_sf: PathBuf::new(),
                    mapping_rate: 0.89,
                    library_size: 1_000_000,
                },
            },
            SampleHapResult {
                haplotype: "hap1".to_string(),
                quant: SalmonQuantOutput {
                    quant_sf: PathBuf::new(),
                    mapping_rate: 0.88,
                    library_size: 1_000_000,
                },
            },
        ];
        // 0.01 diff vs 0.05 threshold → no flag.
        assert!(!detect_asymmetry(&rs, 0.05));
    }

    #[test]
    fn asymmetry_over_threshold_returns_true() {
        use crate::quant::salmon::SalmonQuantOutput;
        let rs = vec![
            SampleHapResult {
                haplotype: "hap0".to_string(),
                quant: SalmonQuantOutput {
                    quant_sf: PathBuf::new(),
                    mapping_rate: 0.89,
                    library_size: 1_000_000,
                },
            },
            SampleHapResult {
                haplotype: "hap1".to_string(),
                quant: SalmonQuantOutput {
                    quant_sf: PathBuf::new(),
                    mapping_rate: 0.72,
                    library_size: 1_000_000,
                },
            },
        ];
        // 0.17 diff → flagged at 0.05 threshold.
        assert!(detect_asymmetry(&rs, 0.05));
    }

    #[test]
    fn asymmetry_single_haplotype_never_flagged() {
        use crate::quant::salmon::SalmonQuantOutput;
        // Degenerate: only one haplotype quant. No pairs to compare,
        // so the function returns false regardless of mapping rate.
        let rs = vec![SampleHapResult {
            haplotype: "hap0".to_string(),
            quant: SalmonQuantOutput {
                quant_sf: PathBuf::new(),
                mapping_rate: 0.50,
                library_size: 100,
            },
        }];
        assert!(!detect_asymmetry(&rs, 0.05));
    }
}
