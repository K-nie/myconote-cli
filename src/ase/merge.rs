//! Combine per-sample-per-haplotype salmon `quant.sf` output into
//! wide allele matrices for ASE.
//!
//! Column layout: `transcript` + one column per `sample_id.hap_name`
//! combination (e.g. `WT1.hap0`, `WT1.hap1`, `KO1.hap0`, `KO1.hap1`).
//! Rows = transcripts, sorted lexicographically for determinism.
//!
//! Unlike `src/quant/merge.rs` (bulk-expression), ASE keeps
//! haplotypes as adjacent columns so downstream R can `pivot_longer`
//! or group by sample trivially. We also emit `ase_summary.tsv`
//! flagging per-transcript informativeness — critical for the
//! binomial ASE tests that should skip transcripts where the two
//! haplotypes carry identical sequences.

use crate::quant::merge::parse_quant_sf;
use crate::utils::error::{MycoNoteError, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// One (sample, haplotype) pair and the path to its `quant.sf`.
#[derive(Debug, Clone)]
pub struct HapQuantInput {
    pub sample_id: String,
    pub haplotype: String,
    pub quant_sf: PathBuf,
}

/// Merged allele matrices. `columns` holds the `sample_id.hap_name`
/// labels in the order they were passed in; `transcripts` is sorted;
/// `counts` and `tpm` are both row-major with the same dimensions.
#[derive(Debug)]
pub struct AseMatrix {
    pub columns: Vec<String>,
    pub transcripts: Vec<String>,
    pub counts: Vec<Vec<f64>>,
    pub tpm: Vec<Vec<f64>>,
}

/// Merge per-sample-per-haplotype quant.sf files into the wide
/// allele matrix. Columns appear in the order given by `inputs` so
/// the dispatcher controls grouping (typical pattern: emit sample A
/// hap0, sample A hap1, sample B hap0, sample B hap1, …).
pub fn merge_allele_quant(inputs: &[HapQuantInput]) -> Result<AseMatrix> {
    if inputs.is_empty() {
        return Err(MycoNoteError::QuantSheet(
            "ase merge: no (sample, haplotype) inputs provided".to_string(),
        ));
    }

    let mut per_col: Vec<BTreeMap<String, crate::quant::merge::QuantRecord>> =
        Vec::with_capacity(inputs.len());
    for inp in inputs {
        per_col.push(parse_quant_sf(&inp.quant_sf)?);
    }

    // Union of transcript IDs across every column, sorted.
    let mut tx_set: BTreeSet<String> = BTreeSet::new();
    for col in &per_col {
        tx_set.extend(col.keys().cloned());
    }
    let transcripts: Vec<String> = tx_set.into_iter().collect();
    let n_tx = transcripts.len();
    let n_cols = inputs.len();

    let mut counts = vec![vec![0.0_f64; n_cols]; n_tx];
    let mut tpm = vec![vec![0.0_f64; n_cols]; n_tx];

    for (j, col) in per_col.iter().enumerate() {
        for (i, tx) in transcripts.iter().enumerate() {
            if let Some(rec) = col.get(tx) {
                counts[i][j] = rec.num_reads;
                tpm[i][j] = rec.tpm;
            }
        }
    }

    let columns: Vec<String> = inputs
        .iter()
        .map(|i| format!("{}.{}", i.sample_id, i.haplotype))
        .collect();

    Ok(AseMatrix {
        columns,
        transcripts,
        counts,
        tpm,
    })
}

impl AseMatrix {
    /// Write the counts matrix as a wide TSV. Integer-valued cells
    /// format without decimals; fractional cells use six decimals —
    /// matches the convention already established in
    /// `src/quant/merge.rs` so downstream R consumers see consistent
    /// precision regardless of which path emitted the matrix.
    pub fn write_counts(&self, path: &Path) -> Result<()> {
        self.write_matrix(path, &self.counts)
    }

    pub fn write_tpm(&self, path: &Path) -> Result<()> {
        self.write_matrix(path, &self.tpm)
    }

    fn write_matrix(&self, path: &Path, data: &[Vec<f64>]) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        write!(w, "transcript")?;
        for col in &self.columns {
            write!(w, "\t{col}")?;
        }
        writeln!(w)?;
        for (i, tx) in self.transcripts.iter().enumerate() {
            write!(w, "{tx}")?;
            for v in &data[i] {
                if v.fract() == 0.0 && v.is_finite() {
                    write!(w, "\t{}", *v as u64)?;
                } else {
                    write!(w, "\t{:.6}", v)?;
                }
            }
            writeln!(w)?;
        }
        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-transcript summary: informativeness, variant counts, asymmetry
// ─────────────────────────────────────────────────────────────────────────────

/// Per-transcript ASE-level metadata. One row per transcript in
/// `ase_summary.tsv`. The dispatcher populates this by joining
/// personalize.rs output (variants applied per haplotype per
/// transcript) with quant-level read counts.
#[derive(Debug, Clone)]
pub struct TranscriptSummaryRow {
    pub transcript_id: String,
    /// Non-silent variants on at least one haplotype (i.e. the two
    /// haplotype sequences differ somewhere in this transcript).
    /// Transcripts with `false` here should be filtered out of the
    /// binomial ASE test downstream.
    pub informative: bool,
    /// Variant counts per haplotype.
    pub n_variants_hap0: usize,
    pub n_variants_hap1: usize,
    /// Max read-count asymmetry across samples for this transcript.
    /// Computed as `max_i |count_{i,hap0} - count_{i,hap1}| /
    /// (count_{i,hap0} + count_{i,hap1} + 1)`. Small denominators
    /// are stabilized by the `+1`.
    pub max_asymmetry: f64,
}

pub fn write_summary(path: &Path, rows: &[TranscriptSummaryRow]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "transcript\tinformative\tn_variants_hap0\tn_variants_hap1\tmax_asymmetry"
    )?;
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.4}",
            r.transcript_id, r.informative, r.n_variants_hap0, r.n_variants_hap1, r.max_asymmetry
        )?;
    }
    Ok(())
}

/// Build a summary row per transcript from the merged matrix. The
/// caller supplies per-transcript variant counts (precomputed from
/// the personalize output) as a pair of HashMaps. Missing
/// transcripts in the maps default to zero variants.
pub fn build_summary(
    matrix: &AseMatrix,
    variants_hap0: &std::collections::HashMap<String, usize>,
    variants_hap1: &std::collections::HashMap<String, usize>,
    hap_names: &[String; 2],
) -> Vec<TranscriptSummaryRow> {
    let mut out = Vec::with_capacity(matrix.transcripts.len());
    // Locate column indices for each haplotype per sample. Columns
    // are labelled `sample_id.hap_name`; we split and group.
    let mut hap_col_indices: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
    for (idx, col) in matrix.columns.iter().enumerate() {
        if let Some(dot) = col.rfind('.') {
            let hap = &col[dot + 1..];
            if hap == hap_names[0] {
                hap_col_indices[0].push(idx);
            } else if hap == hap_names[1] {
                hap_col_indices[1].push(idx);
            }
        }
    }

    for (i, tx) in matrix.transcripts.iter().enumerate() {
        let n0 = variants_hap0.get(tx).copied().unwrap_or(0);
        let n1 = variants_hap1.get(tx).copied().unwrap_or(0);
        let informative = (n0 + n1) > 0 && n0 != n1 || (n0 > 0 && n1 > 0);
        // `informative` = at least one haplotype differs from REF at
        // some position AND the two haplotypes are not trivially
        // identical. For the common ASE case both haps will have
        // ≥ 1 variant each (n0 > 0 && n1 > 0), or one will be zero
        // (het where the other haplotype carries REF).
        let informative = informative || (n0 > 0) != (n1 > 0);

        let mut max_asym = 0.0_f64;
        let n_samples = hap_col_indices[0].len().min(hap_col_indices[1].len());
        for s in 0..n_samples {
            let c0 = matrix.counts[i][hap_col_indices[0][s]];
            let c1 = matrix.counts[i][hap_col_indices[1][s]];
            let denom = c0 + c1 + 1.0;
            let asym = (c0 - c1).abs() / denom;
            if asym > max_asym {
                max_asym = asym;
            }
        }

        out.push(TranscriptSummaryRow {
            transcript_id: tx.clone(),
            informative,
            n_variants_hap0: n0,
            n_variants_hap1: n1,
            max_asymmetry: max_asym,
        });
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use tempfile::TempDir;

    /// Write a minimal quant.sf with the given (transcript, tpm,
    /// num_reads) triples. Used by the merge tests — keeps the
    /// fixture cost in-file.
    fn write_sf(dir: &Path, name: &str, rows: &[(&str, f64, f64)]) -> PathBuf {
        let p = dir.join(name);
        let mut s = String::from("Name\tLength\tEffectiveLength\tTPM\tNumReads\n");
        for (tx, tpm, reads) in rows {
            s.push_str(&format!("{tx}\t1000\t800\t{tpm}\t{reads}\n"));
        }
        std::fs::write(&p, s).unwrap();
        p
    }

    #[test]
    fn merge_two_samples_two_haps() {
        let tmp = TempDir::new().unwrap();
        let inputs = vec![
            HapQuantInput {
                sample_id: "s1".into(),
                haplotype: "hap0".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s1.hap0.sf",
                    &[("tx1", 10.0, 60.0), ("tx2", 20.0, 120.0)],
                ),
            },
            HapQuantInput {
                sample_id: "s1".into(),
                haplotype: "hap1".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s1.hap1.sf",
                    &[("tx1", 11.0, 40.0), ("tx2", 22.0, 110.0)],
                ),
            },
            HapQuantInput {
                sample_id: "s2".into(),
                haplotype: "hap0".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s2.hap0.sf",
                    &[("tx1", 5.0, 100.0), ("tx2", 30.0, 200.0)],
                ),
            },
            HapQuantInput {
                sample_id: "s2".into(),
                haplotype: "hap1".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s2.hap1.sf",
                    &[("tx1", 5.5, 90.0), ("tx2", 33.0, 190.0)],
                ),
            },
        ];
        let m = merge_allele_quant(&inputs).unwrap();
        assert_eq!(m.columns, vec!["s1.hap0", "s1.hap1", "s2.hap0", "s2.hap1"]);
        assert_eq!(m.transcripts, vec!["tx1", "tx2"]);
        // Rows are sorted alphabetically; tx1 is row 0.
        assert_eq!(m.counts[0], vec![60.0, 40.0, 100.0, 90.0]);
        assert_eq!(m.counts[1], vec![120.0, 110.0, 200.0, 190.0]);
    }

    #[test]
    fn merge_fills_missing_transcripts_with_zero() {
        let tmp = TempDir::new().unwrap();
        // hap0 has tx1+tx2; hap1 has tx2+tx3. Union = {tx1, tx2, tx3}.
        let inputs = vec![
            HapQuantInput {
                sample_id: "s1".into(),
                haplotype: "hap0".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s1.hap0.sf",
                    &[("tx1", 5.0, 50.0), ("tx2", 10.0, 100.0)],
                ),
            },
            HapQuantInput {
                sample_id: "s1".into(),
                haplotype: "hap1".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s1.hap1.sf",
                    &[("tx2", 11.0, 110.0), ("tx3", 3.0, 30.0)],
                ),
            },
        ];
        let m = merge_allele_quant(&inputs).unwrap();
        assert_eq!(m.transcripts, vec!["tx1", "tx2", "tx3"]);
        // tx1 is missing from hap1 → 0 in column 1.
        assert_eq!(m.counts[0], vec![50.0, 0.0]);
        // tx2 present in both.
        assert_eq!(m.counts[1], vec![100.0, 110.0]);
        // tx3 missing from hap0 → 0 in column 0.
        assert_eq!(m.counts[2], vec![0.0, 30.0]);
    }

    #[test]
    fn merge_empty_inputs_errors() {
        let err = merge_allele_quant(&[]).unwrap_err();
        assert!(format!("{err}").contains("no (sample, haplotype)"));
    }

    #[test]
    fn write_counts_produces_deterministic_tsv() {
        let tmp = TempDir::new().unwrap();
        let inputs = vec![
            HapQuantInput {
                sample_id: "s1".into(),
                haplotype: "hap0".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s1.hap0.sf",
                    &[("tx2", 20.0, 100.0), ("tx1", 10.0, 5.0)],
                ),
            },
            HapQuantInput {
                sample_id: "s1".into(),
                haplotype: "hap1".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s1.hap1.sf",
                    &[("tx1", 11.0, 6.5), ("tx2", 22.0, 110.0)],
                ),
            },
        ];
        let m = merge_allele_quant(&inputs).unwrap();
        let out = tmp.path().join("counts.tsv");
        m.write_counts(&out).unwrap();
        let written = std::fs::read_to_string(&out).unwrap();
        let expected = "transcript\ts1.hap0\ts1.hap1\n\
                        tx1\t5\t6.500000\n\
                        tx2\t100\t110\n";
        assert_eq!(written, expected);
    }

    #[test]
    fn summary_flags_informative_transcripts() {
        let tmp = TempDir::new().unwrap();
        let inputs = vec![
            HapQuantInput {
                sample_id: "s1".into(),
                haplotype: "hap0".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s1.hap0.sf",
                    &[
                        ("tx1", 10.0, 100.0),
                        ("tx2", 20.0, 50.0),
                        ("tx3", 5.0, 20.0),
                    ],
                ),
            },
            HapQuantInput {
                sample_id: "s1".into(),
                haplotype: "hap1".into(),
                quant_sf: write_sf(
                    tmp.path(),
                    "s1.hap1.sf",
                    &[
                        ("tx1", 10.0, 100.0),
                        ("tx2", 22.0, 150.0),
                        ("tx3", 5.0, 20.0),
                    ],
                ),
            },
        ];
        let m = merge_allele_quant(&inputs).unwrap();

        // tx1: no variants on either haplotype (uninformative).
        // tx2: variants on both haplotypes (informative — ASE-testable).
        // tx3: variant only on hap0 (still informative — hap1 carries REF).
        let mut v0: HashMap<String, usize> = HashMap::new();
        v0.insert("tx2".into(), 3);
        v0.insert("tx3".into(), 1);
        let mut v1: HashMap<String, usize> = HashMap::new();
        v1.insert("tx2".into(), 3);

        let haps = ["hap0".to_string(), "hap1".to_string()];
        let rows = build_summary(&m, &v0, &v1, &haps);
        let by_tx: HashMap<&str, &TranscriptSummaryRow> =
            rows.iter().map(|r| (r.transcript_id.as_str(), r)).collect();

        assert!(!by_tx["tx1"].informative, "tx1 has no variants anywhere");
        assert!(by_tx["tx2"].informative, "tx2 has variants on both haps");
        assert!(
            by_tx["tx3"].informative,
            "tx3 has variants on hap0 only — still informative"
        );
        assert_eq!(by_tx["tx2"].n_variants_hap0, 3);
        assert_eq!(by_tx["tx2"].n_variants_hap1, 3);
    }
}
