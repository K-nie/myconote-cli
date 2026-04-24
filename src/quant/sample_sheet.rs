//! Sample-sheet parser for the `quant` subcommand.
//!
//! TSV format, header required, UTF-8. Columns (column order does not
//! matter; unknown columns are preserved so users can keep custom
//! metadata alongside the sheet):
//!
//! | column        | required | meaning                                          |
//! |---------------|----------|--------------------------------------------------|
//! | sample_id     | yes      | unique, filename-safe `[A-Za-z0-9_-]+`           |
//! | fastq_r1      | yes      | R1 FASTQ path, gzipped or plain                  |
//! | fastq_r2      | no       | R2 for paired-end; empty / missing → single-end  |
//! | condition     | no       | free-text DE grouping label, passed through      |
//! | strandedness  | no       | `unstranded` / `forward` / `reverse` / `auto` (default `auto`) |
//! | batch         | no       | free-text technical batch label                  |
//!
//! Rules:
//!   - Header row required, tab-separated.
//!   - Lines starting with `#` are ignored.
//!   - Blank lines are ignored.
//!   - Relative `fastq_r1` / `fastq_r2` paths are resolved against the
//!     directory that contains the sheet.
//!   - Errors carry a 1-based line number from the original file.

use crate::utils::error::{MycoNoteError, Result};
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

/// Strandedness label from the sheet, mapped to salmon `--libType` at
/// quant time. Using an enum rather than a string so typos surface at
/// parse time, not buried inside a salmon failure message.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strandedness {
    Unstranded,
    Forward,
    Reverse,
    Auto,
}

impl Strandedness {
    /// Translate to the salmon `--libType` token. `paired = true` selects
    /// the `I*` (inward paired) variants; single-end uses bare letters.
    pub fn salmon_libtype(self, paired: bool) -> &'static str {
        match (self, paired) {
            (Strandedness::Unstranded, true) => "IU",
            (Strandedness::Unstranded, false) => "U",
            (Strandedness::Forward, true) => "ISF",
            (Strandedness::Forward, false) => "SF",
            (Strandedness::Reverse, true) => "ISR",
            (Strandedness::Reverse, false) => "SR",
            (Strandedness::Auto, _) => "A",
        }
    }

    fn from_str(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "unstranded" | "u" | "iu" => Some(Self::Unstranded),
            "forward" | "sf" | "isf" | "fr" => Some(Self::Forward),
            "reverse" | "sr" | "isr" | "rf" => Some(Self::Reverse),
            "auto" | "a" | "" => Some(Self::Auto),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Sample {
    pub sample_id: String,
    pub fastq_r1: PathBuf,
    pub fastq_r2: Option<PathBuf>,
    pub condition: Option<String>,
    pub strandedness: Strandedness,
    pub batch: Option<String>,
    /// Line number in the original sheet (1-based, counting headers +
    /// comments), surfaced in error messages and in the normalized
    /// output for downstream bookkeeping.
    pub source_line: usize,
    /// Columns present in the sheet that are not part of the known
    /// schema, preserved verbatim so a normalized output round-trip
    /// keeps user metadata. Key = header, value = cell.
    pub extras: BTreeMap<String, String>,
}

impl Sample {
    pub fn is_paired(&self) -> bool {
        self.fastq_r2.is_some()
    }
}

#[derive(Debug, Clone)]
pub struct SampleSheet {
    pub samples: Vec<Sample>,
    /// Header row as it appeared in the input, in column order. Needed
    /// to round-trip extras in the same order when emitting the
    /// normalized sheet.
    pub input_header: Vec<String>,
    /// Directory the sheet was loaded from. Relative paths resolve
    /// against this; absolute paths are passed through untouched.
    pub root_dir: PathBuf,
}

/// Parse a TSV sample sheet from disk. All fastq paths are resolved to
/// absolute but *not* canonicalized (symlinks are preserved) and *not*
/// checked for existence — the caller (quant dispatcher) runs the
/// existence check after parsing so users see every missing file at
/// once instead of one per re-run.
pub fn parse_sheet(path: &Path) -> Result<SampleSheet> {
    let text = fs::read_to_string(path).map_err(|e| {
        MycoNoteError::QuantSheet(format!("failed to read {}: {}", path.display(), e))
    })?;
    let root_dir = path
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));

    parse_sheet_from_str(&text, &root_dir)
}

/// Parse a sheet from an in-memory TSV. Factored out so unit tests can
/// exercise the parser without touching disk.
pub fn parse_sheet_from_str(text: &str, root_dir: &Path) -> Result<SampleSheet> {
    let mut header: Option<Vec<String>> = None;
    let mut samples: Vec<Sample> = Vec::new();
    let mut seen_ids: HashSet<String> = HashSet::new();

    for (idx, raw_line) in text.lines().enumerate() {
        let line_num = idx + 1;
        let trimmed = raw_line.trim_end_matches('\r');

        // Skip blank and comment lines. A comment starts with `#` as
        // the first non-whitespace character — allowing leading
        // whitespace would be surprising in a TSV but we keep it
        // permissive for hand-edited sheets.
        let first_non_ws = trimmed.trim_start();
        if first_non_ws.is_empty() || first_non_ws.starts_with('#') {
            continue;
        }

        let cells: Vec<&str> = trimmed.split('\t').collect();

        if header.is_none() {
            // First non-blank / non-comment line is the header.
            let hdr: Vec<String> = cells.iter().map(|s| s.trim().to_string()).collect();
            require_header(&hdr, line_num)?;
            header = Some(hdr);
            continue;
        }

        let hdr = header.as_ref().unwrap();
        if cells.len() > hdr.len() {
            return Err(MycoNoteError::QuantSheet(format!(
                "line {line_num}: row has {} columns but header has {}",
                cells.len(),
                hdr.len()
            )));
        }

        let sample = parse_row(hdr, &cells, root_dir, line_num)?;

        if !seen_ids.insert(sample.sample_id.clone()) {
            return Err(MycoNoteError::QuantSheet(format!(
                "line {line_num}: duplicate sample_id '{}'",
                sample.sample_id
            )));
        }
        samples.push(sample);
    }

    let header = header.ok_or_else(|| {
        MycoNoteError::QuantSheet("sheet is empty (no header row found)".to_string())
    })?;

    if samples.is_empty() {
        return Err(MycoNoteError::QuantSheet(
            "sheet has a header but no sample rows".to_string(),
        ));
    }

    Ok(SampleSheet {
        samples,
        input_header: header,
        root_dir: root_dir.to_path_buf(),
    })
}

fn require_header(hdr: &[String], line_num: usize) -> Result<()> {
    let mut seen = HashSet::new();
    for h in hdr {
        if !seen.insert(h.as_str()) {
            return Err(MycoNoteError::QuantSheet(format!(
                "line {line_num}: duplicate column '{h}' in header"
            )));
        }
    }
    for required in ["sample_id", "fastq_r1"] {
        if !hdr.iter().any(|h| h == required) {
            return Err(MycoNoteError::QuantSheet(format!(
                "header missing required column '{required}'"
            )));
        }
    }
    Ok(())
}

fn parse_row(hdr: &[String], cells: &[&str], root_dir: &Path, line_num: usize) -> Result<Sample> {
    let get = |name: &str| -> Option<&str> {
        hdr.iter()
            .position(|h| h == name)
            .and_then(|i| cells.get(i))
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
    };

    let sample_id = get("sample_id")
        .ok_or_else(|| MycoNoteError::QuantSheet(format!("line {line_num}: sample_id is empty")))?
        .to_string();
    require_filename_safe(&sample_id, line_num)?;

    let fastq_r1_raw = get("fastq_r1")
        .ok_or_else(|| MycoNoteError::QuantSheet(format!("line {line_num}: fastq_r1 is empty")))?;
    let fastq_r1 = resolve_path(fastq_r1_raw, root_dir);

    let fastq_r2 = get("fastq_r2").map(|p| resolve_path(p, root_dir));

    let condition = get("condition").map(|s| s.to_string());
    let batch = get("batch").map(|s| s.to_string());

    let strand_raw = get("strandedness").unwrap_or("auto");
    let strandedness = Strandedness::from_str(strand_raw).ok_or_else(|| {
        MycoNoteError::QuantSheet(format!(
            "line {line_num}: unknown strandedness '{strand_raw}' \
             (expected unstranded / forward / reverse / auto)"
        ))
    })?;

    // Collect any extras — any header column that is not one of the
    // known schema columns. Stored in a BTreeMap so output order is
    // deterministic regardless of input column order.
    let known = [
        "sample_id",
        "fastq_r1",
        "fastq_r2",
        "condition",
        "strandedness",
        "batch",
    ];
    let mut extras = BTreeMap::new();
    for (i, col) in hdr.iter().enumerate() {
        if known.contains(&col.as_str()) {
            continue;
        }
        let value = cells.get(i).map(|s| s.trim()).unwrap_or("").to_string();
        extras.insert(col.clone(), value);
    }

    Ok(Sample {
        sample_id,
        fastq_r1,
        fastq_r2,
        condition,
        strandedness,
        batch,
        source_line: line_num,
        extras,
    })
}

fn require_filename_safe(id: &str, line_num: usize) -> Result<()> {
    if id.is_empty() {
        return Err(MycoNoteError::QuantSheet(format!(
            "line {line_num}: sample_id is empty"
        )));
    }
    // Allow ASCII letters, digits, underscore, dash. Anything else
    // risks breaking across shells, salmon output-dir semantics, or
    // downstream R scripts that use column names as symbols.
    if !id
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
    {
        return Err(MycoNoteError::QuantSheet(format!(
            "line {line_num}: sample_id '{id}' must match [A-Za-z0-9_-]+ \
             (no spaces, dots, slashes, or unicode)"
        )));
    }
    Ok(())
}

fn resolve_path(p: &str, root_dir: &Path) -> PathBuf {
    let pb = PathBuf::from(p);
    if pb.is_absolute() {
        pb
    } else {
        root_dir.join(pb)
    }
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn root() -> &'static Path {
        Path::new("/tmp")
    }

    #[test]
    fn parse_minimal_paired() {
        let sheet = "sample_id\tfastq_r1\tfastq_r2\n\
                     WT_rep1\tr1.fq.gz\tr2.fq.gz\n";
        let s = parse_sheet_from_str(sheet, root()).unwrap();
        assert_eq!(s.samples.len(), 1);
        let sample = &s.samples[0];
        assert_eq!(sample.sample_id, "WT_rep1");
        assert_eq!(sample.fastq_r1, root().join("r1.fq.gz"));
        assert_eq!(sample.fastq_r2, Some(root().join("r2.fq.gz")));
        assert!(sample.is_paired());
        assert_eq!(sample.strandedness, Strandedness::Auto);
    }

    #[test]
    fn parse_single_end_missing_r2_column() {
        let sheet = "sample_id\tfastq_r1\n\
                     run1\tsingle.fq.gz\n";
        let s = parse_sheet_from_str(sheet, root()).unwrap();
        assert_eq!(s.samples.len(), 1);
        assert!(!s.samples[0].is_paired());
        assert!(s.samples[0].fastq_r2.is_none());
    }

    #[test]
    fn parse_single_end_empty_r2_cell() {
        // An r2 column that exists but is blank should also mean SE.
        let sheet = "sample_id\tfastq_r1\tfastq_r2\n\
                     run1\tsingle.fq.gz\t\n";
        let s = parse_sheet_from_str(sheet, root()).unwrap();
        assert!(!s.samples[0].is_paired());
    }

    #[test]
    fn parse_absolute_path_passes_through() {
        let sheet = "sample_id\tfastq_r1\n\
                     run1\t/abs/path/to/read.fq.gz\n";
        let s = parse_sheet_from_str(sheet, root()).unwrap();
        assert_eq!(
            s.samples[0].fastq_r1,
            PathBuf::from("/abs/path/to/read.fq.gz")
        );
    }

    #[test]
    fn parse_strandedness_variants() {
        let sheet = "sample_id\tfastq_r1\tstrandedness\n\
                     a\tA.fq\tforward\n\
                     b\tB.fq\treverse\n\
                     c\tC.fq\tunstranded\n\
                     d\tD.fq\tauto\n\
                     e\tE.fq\tISR\n";
        let s = parse_sheet_from_str(sheet, root()).unwrap();
        assert_eq!(s.samples[0].strandedness, Strandedness::Forward);
        assert_eq!(s.samples[1].strandedness, Strandedness::Reverse);
        assert_eq!(s.samples[2].strandedness, Strandedness::Unstranded);
        assert_eq!(s.samples[3].strandedness, Strandedness::Auto);
        assert_eq!(s.samples[4].strandedness, Strandedness::Reverse);
    }

    #[test]
    fn reject_duplicate_sample_id() {
        let sheet = "sample_id\tfastq_r1\n\
                     WT\ta.fq\n\
                     WT\tb.fq\n";
        let err = parse_sheet_from_str(sheet, root()).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("duplicate sample_id"), "got: {msg}");
        assert!(msg.contains("line 3"), "got: {msg}");
    }

    #[test]
    fn reject_missing_required_column() {
        let sheet = "sample_id\tother\n\
                     a\tx\n";
        let err = parse_sheet_from_str(sheet, root()).unwrap_err();
        assert!(
            format!("{err}").contains("fastq_r1"),
            "expected missing-fastq_r1 complaint, got: {err}"
        );
    }

    #[test]
    fn reject_unsafe_sample_id() {
        for bad in [
            "has space",
            "has.dot",
            "has/slash",
            "has:colon",
            "naïve",
            "",
        ] {
            let sheet = format!("sample_id\tfastq_r1\n{bad}\tr.fq\n");
            let res = parse_sheet_from_str(&sheet, root());
            assert!(res.is_err(), "sample_id '{bad}' should be rejected");
        }
    }

    #[test]
    fn reject_unknown_strandedness() {
        let sheet = "sample_id\tfastq_r1\tstrandedness\n\
                     a\tr.fq\tsideways\n";
        let err = parse_sheet_from_str(sheet, root()).unwrap_err();
        assert!(format!("{err}").contains("sideways"), "got: {err}");
    }

    #[test]
    fn comments_and_blank_lines_skipped() {
        let sheet = "# a header comment\n\
                     \n\
                     sample_id\tfastq_r1\n\
                     # inline comment before data\n\
                     \n\
                     a\tr.fq\n";
        let s = parse_sheet_from_str(sheet, root()).unwrap();
        assert_eq!(s.samples.len(), 1);
        assert_eq!(s.samples[0].source_line, 6);
    }

    #[test]
    fn unknown_columns_preserved_in_extras() {
        let sheet = "sample_id\tfastq_r1\tbarcode\tread_length\n\
                     a\tr.fq\tACGT\t150\n";
        let s = parse_sheet_from_str(sheet, root()).unwrap();
        let extras = &s.samples[0].extras;
        assert_eq!(extras.get("barcode").map(String::as_str), Some("ACGT"));
        assert_eq!(extras.get("read_length").map(String::as_str), Some("150"));
    }

    #[test]
    fn empty_sheet_is_error() {
        let err = parse_sheet_from_str("", root()).unwrap_err();
        assert!(format!("{err}").contains("empty"));
    }

    #[test]
    fn header_only_is_error() {
        let err = parse_sheet_from_str("sample_id\tfastq_r1\n", root()).unwrap_err();
        assert!(format!("{err}").contains("no sample rows"));
    }

    #[test]
    fn salmon_libtype_mapping() {
        assert_eq!(Strandedness::Auto.salmon_libtype(true), "A");
        assert_eq!(Strandedness::Auto.salmon_libtype(false), "A");
        assert_eq!(Strandedness::Unstranded.salmon_libtype(true), "IU");
        assert_eq!(Strandedness::Unstranded.salmon_libtype(false), "U");
        assert_eq!(Strandedness::Forward.salmon_libtype(true), "ISF");
        assert_eq!(Strandedness::Forward.salmon_libtype(false), "SF");
        assert_eq!(Strandedness::Reverse.salmon_libtype(true), "ISR");
        assert_eq!(Strandedness::Reverse.salmon_libtype(false), "SR");
    }
}
