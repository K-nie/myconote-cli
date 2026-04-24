//! Merge per-sample `quant.sf` files into wide count + TPM matrices.
//!
//! Salmon's `quant.sf` is the tximport input format — one file per
//! sample, tab-separated, with columns:
//!
//! ```text
//! Name    Length    EffectiveLength    TPM    NumReads
//! ```
//!
//! We keep the per-sample files untouched so tximport/DESeq2 works
//! directly, and additionally emit two wide matrices (transcripts ×
//! samples) for rapid triage:
//!
//!   - `counts.tsv` — estimated counts (`NumReads`), the matrix DESeq2
//!     will construct from tximport anyway
//!   - `tpm.tsv` — TPM-normalized expression
//!
//! Transcript rows are sorted ASCII-ascending for determinism
//! (HashMap iteration order is otherwise not reproducible, per the
//! project's determinism rule).
//!
//! Missing transcripts across samples are impossible when all samples
//! share the same salmon index — which is the pipeline's invariant —
//! but we still tolerate sparse rows defensively: any missing
//! (transcript, sample) cell is written as `0`.

use crate::utils::error::{MycoNoteError, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// One record from a salmon `quant.sf`. `tpm` and `num_reads` are the
/// only columns we expose downstream; `length` / `effective_length`
/// are preserved from the per-sample file for tximport to consume.
#[derive(Debug, Clone, PartialEq)]
pub struct QuantRecord {
    pub tpm: f64,
    pub num_reads: f64,
}

/// Parse a single `quant.sf` into a `transcript_id → QuantRecord` map.
/// Uses `BTreeMap` so callers iterating on the result see deterministic
/// order; the merger below relies on this for stable output.
pub fn parse_quant_sf(path: &Path) -> Result<BTreeMap<String, QuantRecord>> {
    let f = File::open(path).map_err(|e| MycoNoteError::QuantTool {
        tool: "salmon".to_string(),
        message: format!("cannot open {}: {}", path.display(), e),
    })?;
    let reader = BufReader::new(f);
    parse_quant_sf_reader(reader, path)
}

fn parse_quant_sf_reader<R: BufRead>(
    reader: R,
    origin: &Path,
) -> Result<BTreeMap<String, QuantRecord>> {
    let mut out = BTreeMap::new();
    let mut header_checked = false;

    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| MycoNoteError::QuantTool {
            tool: "salmon".to_string(),
            message: format!("read error on {} line {}: {}", origin.display(), i + 1, e),
        })?;
        if line.is_empty() {
            continue;
        }

        if !header_checked {
            // Salmon emits the header as `Name\tLength\tEffectiveLength\tTPM\tNumReads`.
            // We check the first three column names rather than all five so minor
            // salmon-version column additions don't break parsing — extra columns
            // after NumReads are ignored.
            let cols: Vec<&str> = line.split('\t').collect();
            let ok = cols.len() >= 5
                && cols[0] == "Name"
                && cols[1] == "Length"
                && cols[3] == "TPM"
                && cols[4] == "NumReads";
            if !ok {
                return Err(MycoNoteError::QuantTool {
                    tool: "salmon".to_string(),
                    message: format!(
                        "{}: unexpected quant.sf header (expected 'Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads'): {}",
                        origin.display(),
                        line
                    ),
                });
            }
            header_checked = true;
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            return Err(MycoNoteError::QuantTool {
                tool: "salmon".to_string(),
                message: format!(
                    "{}: malformed row (expected 5 columns, got {}): {}",
                    origin.display(),
                    cols.len(),
                    line
                ),
            });
        }
        let name = cols[0].to_string();
        let tpm: f64 = cols[3].parse().map_err(|_| MycoNoteError::QuantTool {
            tool: "salmon".to_string(),
            message: format!(
                "{}: non-numeric TPM for transcript {}: {}",
                origin.display(),
                name,
                cols[3]
            ),
        })?;
        let num_reads: f64 = cols[4].parse().map_err(|_| MycoNoteError::QuantTool {
            tool: "salmon".to_string(),
            message: format!(
                "{}: non-numeric NumReads for transcript {}: {}",
                origin.display(),
                name,
                cols[4]
            ),
        })?;

        out.insert(name, QuantRecord { tpm, num_reads });
    }

    if !header_checked {
        return Err(MycoNoteError::QuantTool {
            tool: "salmon".to_string(),
            message: format!("{}: quant.sf is empty", origin.display()),
        });
    }

    Ok(out)
}

/// Combined per-transcript × per-sample matrix. `samples` is the
/// column order (preserved from input), `transcripts` is the sorted
/// row order, and the two vectors of vectors carry the data in the
/// same row/column order.
#[derive(Debug)]
pub struct MergedMatrix {
    pub samples: Vec<String>,
    pub transcripts: Vec<String>,
    pub counts: Vec<Vec<f64>>,
    pub tpm: Vec<Vec<f64>>,
}

impl MergedMatrix {
    /// Write the estimated-counts matrix as a wide TSV to `path`.
    pub fn write_counts(&self, path: &Path) -> Result<()> {
        self.write_matrix(path, &self.counts)
    }

    /// Write the TPM matrix as a wide TSV to `path`.
    pub fn write_tpm(&self, path: &Path) -> Result<()> {
        self.write_matrix(path, &self.tpm)
    }

    fn write_matrix(&self, path: &Path, data: &[Vec<f64>]) -> Result<()> {
        let mut f = File::create(path)?;
        // Header: leading 'transcript' column, then sample IDs.
        write!(f, "transcript")?;
        for s in &self.samples {
            write!(f, "\t{}", s)?;
        }
        writeln!(f)?;

        for (i, tx) in self.transcripts.iter().enumerate() {
            write!(f, "{}", tx)?;
            for v in &data[i] {
                // Integer-valued counts stay integer-looking; anything
                // fractional gets 6 decimals. Matches the precision
                // R/tximport readers expect without introducing
                // parse-time rounding errors.
                if v.fract() == 0.0 && v.is_finite() {
                    write!(f, "\t{}", *v as u64)?;
                } else {
                    write!(f, "\t{:.6}", v)?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

/// Merge a list of `(sample_id, quant.sf path)` pairs into a wide
/// matrix. Output row order is sorted ASCII-ascending over the union
/// of transcript IDs seen across samples. Output column order matches
/// the input vector (callers preserve the sample-sheet order).
pub fn merge(samples: &[(String, std::path::PathBuf)]) -> Result<MergedMatrix> {
    if samples.is_empty() {
        return Err(MycoNoteError::QuantSheet(
            "merge: no samples provided".to_string(),
        ));
    }

    // Parse each sample's quant.sf first so we can fail fast on any
    // bad input before allocating the (potentially large) matrix.
    let mut per_sample: Vec<BTreeMap<String, QuantRecord>> = Vec::with_capacity(samples.len());
    for (_, path) in samples {
        per_sample.push(parse_quant_sf(path)?);
    }

    // Union of transcript IDs across samples. BTreeSet gives sorted
    // iteration; the project's determinism rule forbids HashMap here.
    let mut tx_union: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for s in &per_sample {
        tx_union.extend(s.keys().cloned());
    }
    let transcripts: Vec<String> = tx_union.into_iter().collect();

    let n_tx = transcripts.len();
    let n_samples = samples.len();
    let mut counts = vec![vec![0.0_f64; n_samples]; n_tx];
    let mut tpm = vec![vec![0.0_f64; n_samples]; n_tx];

    for (j, sample_map) in per_sample.iter().enumerate() {
        for (i, tx) in transcripts.iter().enumerate() {
            if let Some(rec) = sample_map.get(tx) {
                counts[i][j] = rec.num_reads;
                tpm[i][j] = rec.tpm;
            }
        }
    }

    Ok(MergedMatrix {
        samples: samples.iter().map(|(id, _)| id.clone()).collect(),
        transcripts,
        counts,
        tpm,
    })
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::path::PathBuf;
    use tempfile::TempDir;

    fn sf_text() -> &'static str {
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n\
         tx001\t1200\t950.5\t123.456\t42.0\n\
         tx002\t800\t650.2\t0\t0\n\
         tx003\t1500\t1250.0\t7.5\t100.2\n"
    }

    #[test]
    fn parse_quant_sf_basic() {
        let rdr = Cursor::new(sf_text());
        let m = parse_quant_sf_reader(rdr, Path::new("test.sf")).unwrap();
        assert_eq!(m.len(), 3);
        assert_eq!(m["tx001"].tpm, 123.456);
        assert_eq!(m["tx001"].num_reads, 42.0);
        assert_eq!(m["tx002"].num_reads, 0.0);
        assert_eq!(m["tx003"].num_reads, 100.2);
    }

    #[test]
    fn parse_quant_sf_rejects_bad_header() {
        let bad = "FOO\tBAR\tBAZ\n";
        let rdr = Cursor::new(bad);
        let err = parse_quant_sf_reader(rdr, Path::new("bad.sf")).unwrap_err();
        assert!(format!("{err}").contains("unexpected quant.sf header"));
    }

    #[test]
    fn parse_quant_sf_rejects_empty_file() {
        let rdr = Cursor::new("");
        let err = parse_quant_sf_reader(rdr, Path::new("empty.sf")).unwrap_err();
        assert!(format!("{err}").contains("empty"));
    }

    #[test]
    fn parse_quant_sf_rejects_non_numeric_tpm() {
        let bad = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n\
                   tx001\t100\t80\tNaN-ish\t5\n";
        let rdr = Cursor::new(bad);
        let err = parse_quant_sf_reader(rdr, Path::new("bad.sf")).unwrap_err();
        assert!(format!("{err}").contains("non-numeric TPM"));
    }

    /// Write a fake quant.sf file and return its path. Used by merge tests.
    fn write_sf(dir: &Path, name: &str, rows: &[(&str, f64, f64)]) -> PathBuf {
        let p = dir.join(name);
        let mut contents = String::from("Name\tLength\tEffectiveLength\tTPM\tNumReads\n");
        for (tx, tpm, reads) in rows {
            contents.push_str(&format!("{}\t1000\t800\t{}\t{}\n", tx, tpm, reads));
        }
        std::fs::write(&p, contents).unwrap();
        p
    }

    #[test]
    fn merge_two_samples_same_transcripts() {
        let tmp = TempDir::new().unwrap();
        let a = write_sf(
            tmp.path(),
            "a.sf",
            &[("tx01", 10.0, 5.0), ("tx02", 20.0, 100.0)],
        );
        let b = write_sf(
            tmp.path(),
            "b.sf",
            &[("tx01", 11.0, 6.0), ("tx02", 22.0, 110.0)],
        );

        let m = merge(&[("s1".to_string(), a), ("s2".to_string(), b)]).unwrap();

        assert_eq!(m.samples, vec!["s1", "s2"]);
        assert_eq!(m.transcripts, vec!["tx01", "tx02"]);
        assert_eq!(m.counts[0], vec![5.0, 6.0]);
        assert_eq!(m.counts[1], vec![100.0, 110.0]);
        assert_eq!(m.tpm[0], vec![10.0, 11.0]);
        assert_eq!(m.tpm[1], vec![20.0, 22.0]);
    }

    #[test]
    fn merge_fills_missing_transcripts_with_zero() {
        let tmp = TempDir::new().unwrap();
        // s1 has tx01 + tx02; s2 has tx02 + tx03. Union = {tx01, tx02, tx03}.
        // Missing cells should be 0, and row order must be sorted.
        let a = write_sf(
            tmp.path(),
            "a.sf",
            &[("tx01", 10.0, 5.0), ("tx02", 20.0, 100.0)],
        );
        let b = write_sf(
            tmp.path(),
            "b.sf",
            &[("tx02", 22.0, 110.0), ("tx03", 30.0, 7.0)],
        );

        let m = merge(&[("s1".to_string(), a), ("s2".to_string(), b)]).unwrap();

        assert_eq!(m.transcripts, vec!["tx01", "tx02", "tx03"]);
        assert_eq!(m.counts[0], vec![5.0, 0.0]); // tx01: s1=5, s2=missing
        assert_eq!(m.counts[1], vec![100.0, 110.0]); // tx02: both present
        assert_eq!(m.counts[2], vec![0.0, 7.0]); // tx03: s1=missing, s2=7
    }

    #[test]
    fn merge_writes_deterministic_tsv() {
        let tmp = TempDir::new().unwrap();
        let a = write_sf(
            tmp.path(),
            "a.sf",
            &[("tx02", 20.0, 100.0), ("tx01", 10.0, 5.0)],
        );
        let b = write_sf(
            tmp.path(),
            "b.sf",
            &[("tx01", 11.0, 6.5), ("tx02", 22.0, 110.0)],
        );

        let m = merge(&[("s1".to_string(), a), ("s2".to_string(), b)]).unwrap();

        let out = tmp.path().join("counts.tsv");
        m.write_counts(&out).unwrap();
        let written = std::fs::read_to_string(&out).unwrap();
        // Row order must be alpha-sorted; sample column order preserved.
        // Integer-valued cells (5, 100, 110) format without decimals;
        // the fractional 6.5 keeps its precision.
        let expected = "transcript\ts1\ts2\n\
                        tx01\t5\t6.500000\n\
                        tx02\t100\t110\n";
        assert_eq!(written, expected);
    }

    #[test]
    fn merge_empty_input_errors() {
        let err = merge(&[]).unwrap_err();
        assert!(format!("{err}").contains("no samples"));
    }
}
