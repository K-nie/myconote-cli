use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct GFFRecord {
    pub seqid: String,
    pub source: String,
    pub feature_type: String,
    pub start: u64,
    pub end: u64,
    pub score: Option<f64>,
    pub strand: char,
    pub phase: Option<u8>,
    pub attributes: HashMap<String, String>,
}

impl GFFRecord {
    pub fn from_line(line: &str, line_num: usize) -> Result<Self> {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return Err(MycoNoteError::ParseError {
                line: line_num,
                message: "Skipping comment or empty line".to_string(),
            });
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 {
            return Err(MycoNoteError::ParseError {
                line: line_num,
                message: format!("Expected 9 fields, found {}", fields.len()),
            });
        }

        let seqid = fields[0].to_string();
        let source = fields[1].to_string();
        let feature_type = fields[2].to_string();

        let start = fields[3]
            .parse::<u64>()
            .map_err(|_| MycoNoteError::ParseError {
                line: line_num,
                message: format!("Invalid start coordinate: {}", fields[3]),
            })?;

        let end = fields[4]
            .parse::<u64>()
            .map_err(|_| MycoNoteError::ParseError {
                line: line_num,
                message: format!("Invalid end coordinate: {}", fields[4]),
            })?;

        if start > end {
            return Err(MycoNoteError::ParseError {
                line: line_num,
                message: format!("Start ({}) > end ({})", start, end),
            });
        }

        let score = if fields[5] == "." {
            None
        } else {
            Some(
                fields[5]
                    .parse::<f64>()
                    .map_err(|_| MycoNoteError::ParseError {
                        line: line_num,
                        message: format!("Invalid score: {}", fields[5]),
                    })?,
            )
        };

        let strand = if fields[6].len() == 1 {
            fields[6].chars().next().unwrap()
        } else {
            '.'
        };

        let phase = if fields[7] == "." {
            None
        } else {
            Some(
                fields[7]
                    .parse::<u8>()
                    .map_err(|_| MycoNoteError::ParseError {
                        line: line_num,
                        message: format!("Invalid phase: {}", fields[7]),
                    })?,
            )
        };

        let mut attributes = HashMap::new();
        if fields[8] != "." {
            for pair in fields[8].split(';') {
                if pair.is_empty() {
                    continue;
                }
                let parts: Vec<&str> = pair.splitn(2, '=').collect();
                if parts.len() == 2 {
                    attributes.insert(parts[0].to_string(), parts[1].to_string());
                }
            }
        }

        Ok(GFFRecord {
            seqid,
            source,
            feature_type,
            start,
            end,
            score,
            strand,
            phase,
            attributes,
        })
    }

    pub fn id(&self) -> Option<&String> {
        self.attributes.get("ID")
    }

    pub fn parent(&self) -> Option<&String> {
        self.attributes.get("Parent")
    }

    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }

    /// Returns the score field as a GFF3 string (`"."` when absent).
    pub fn score_str(&self) -> String {
        match self.score {
            Some(s) => format!("{}", s),
            None => ".".to_string(),
        }
    }

    /// Returns the phase field as a GFF3 string (`"."` when absent).
    pub fn phase_str(&self) -> String {
        match self.phase {
            Some(p) => format!("{}", p),
            None => ".".to_string(),
        }
    }

    /// Serialise the record back to a GFF3 tab-separated line.
    pub fn to_gff3_line(&self) -> String {
        let attrs: String = self
            .attributes
            .iter()
            .map(|(k, v)| format!("{}={}", k, v))
            .collect::<Vec<_>>()
            .join(";");
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.seqid,
            self.source,
            self.feature_type,
            self.start,
            self.end,
            self.score_str(),
            self.strand,
            self.phase_str(),
            if attrs.is_empty() {
                ".".to_string()
            } else {
                attrs
            }
        )
    }
}

pub struct GFFReader {
    reader: BufReader<File>,
    line_num: usize,
}

impl GFFReader {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(GFFReader {
            reader: BufReader::new(file),
            line_num: 0,
        })
    }

    pub fn next_record(&mut self) -> Result<Option<GFFRecord>> {
        let mut line = String::new();
        loop {
            line.clear();
            let bytes_read = self.reader.read_line(&mut line)?;
            if bytes_read == 0 {
                return Ok(None);
            }
            self.line_num += 1;

            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            match GFFRecord::from_line(line, self.line_num) {
                Ok(record) => return Ok(Some(record)),
                Err(e) => return Err(e),
            }
        }
    }

    pub fn line_num(&self) -> usize {
        self.line_num
    }
}

impl Iterator for GFFReader {
    type Item = Result<GFFRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Fold a transcript's stop codon into its terminal CDS (and coincident exon).
///
/// Ab initio predictors (Augustus, SNAP, GeneMark) emit the stop codon as its
/// own 3-bp feature sitting just *outside* the CDS — e.g. on `+` strand a CDS
/// ending at 2950 with a `stop_codon` at 2951-2953. NCBI / GenBank / RefSeq
/// instead run the terminal CDS *through* the stop codon, so a RefSeq CDS for
/// that gene would end at 2953. A comparator that scores exact CDS structure
/// against RefSeq therefore marks every one of our genes wrong by 3 bp at its
/// 3' end. This folds those 3 bp back in so our emitted CDS matches the
/// convention the reference uses.
///
/// `records` must be the feature rows of a single transcript (its CDS, exon,
/// stop_codon, …). The call is idempotent: if the CDS already runs through the
/// stop codon (the adjacency test below fails) nothing changes, so it is safe
/// to run on output we have already corrected. A model with no `stop_codon`
/// row — a partial / edge gene — is left untouched; we never invent coding
/// bases or push past where a stop was actually called.
pub fn include_stop_codon_in_cds(records: &mut [GFFRecord]) {
    // Must have both a CDS to extend and at least one stop codon to fold in.
    // A model with no stop (partial / edge gene) is left untouched — we never
    // invent coding bases.
    if !records.iter().any(|r| r.feature_type == "CDS") {
        return;
    }
    if !records.iter().any(|r| r.feature_type == "stop_codon") {
        return;
    }

    // Strand drives which end is the 3' (translation) end. Prefer the CDS's
    // strand; fall back to the stop codon, then to any record, so a '.' on the
    // CDS row (seen from some converters) can't silently flip us onto the wrong
    // branch.
    let strand = records
        .iter()
        .find(|r| r.feature_type == "CDS")
        .map(|r| r.strand)
        .filter(|&s| s == '+' || s == '-')
        .or_else(|| {
            records
                .iter()
                .find(|r| r.feature_type == "stop_codon")
                .map(|r| r.strand)
                .filter(|&s| s == '+' || s == '-')
        })
        .unwrap_or('+');

    if strand == '-' {
        // Translation runs high→low coordinate, so the terminal (stop-bearing)
        // CDS segment is the one with the smallest start, and the stop sits
        // immediately below it: stop_end == cds_start - 1. Scan ALL stop rows
        // for the one that is actually flush against that boundary rather than
        // trusting the first stop_codon encountered (multi-segment / multi-
        // isoform models can carry more than one).
        let cds_start = records
            .iter()
            .filter(|r| r.feature_type == "CDS")
            .map(|r| r.start)
            .min()
            .unwrap();
        let fold_to = records
            .iter()
            .filter(|r| r.feature_type == "stop_codon")
            .find(|s| s.start < cds_start && s.end + 1 == cds_start)
            .map(|s| s.start);
        if let Some(new_start) = fold_to {
            for r in records.iter_mut() {
                // Extend the terminal CDS segment and the exon coincident with
                // it (UTR-off predictors give exon.start == cds.start at that
                // boundary); internal segments are left alone.
                if (r.feature_type == "CDS" || r.feature_type == "exon") && r.start == cds_start {
                    r.start = new_start;
                }
            }
        }
    } else {
        // '+' (and unknown '.') — translation runs low→high coordinate, so the
        // terminal CDS segment has the largest end and the stop sits just above
        // it: stop_start == cds_end + 1.
        let cds_end = records
            .iter()
            .filter(|r| r.feature_type == "CDS")
            .map(|r| r.end)
            .max()
            .unwrap();
        let fold_to = records
            .iter()
            .filter(|r| r.feature_type == "stop_codon")
            .find(|s| s.start == cds_end + 1 && s.end > cds_end)
            .map(|s| s.end);
        if let Some(new_end) = fold_to {
            for r in records.iter_mut() {
                if (r.feature_type == "CDS" || r.feature_type == "exon") && r.end == cds_end {
                    r.end = new_end;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(feature_type: &str, start: u64, end: u64, strand: char) -> GFFRecord {
        GFFRecord {
            seqid: "chr1".into(),
            source: "test".into(),
            feature_type: feature_type.into(),
            start,
            end,
            score: None,
            strand,
            phase: Some(0),
            attributes: HashMap::new(),
        }
    }

    fn find(records: &[GFFRecord], feature_type: &str) -> (u64, u64) {
        let r = records
            .iter()
            .find(|r| r.feature_type == feature_type)
            .expect("feature present");
        (r.start, r.end)
    }

    #[test]
    fn plus_strand_folds_stop_into_terminal_cds_and_exon() {
        // The real-run gene: CDS 1802-2950, stop_codon 2951-2953.
        let mut records = vec![
            rec("CDS", 1802, 2950, '+'),
            rec("exon", 1802, 2950, '+'),
            rec("stop_codon", 2951, 2953, '+'),
        ];
        include_stop_codon_in_cds(&mut records);
        assert_eq!(find(&records, "CDS"), (1802, 2953), "CDS end == stop end");
        assert_eq!(find(&records, "exon"), (1802, 2953), "exon follows CDS");
        assert_eq!(
            find(&records, "stop_codon"),
            (2951, 2953),
            "stop kept as-is"
        );
    }

    #[test]
    fn minus_strand_folds_stop_into_terminal_cds_and_exon() {
        // Mirror image: translation 3' end is the low coordinate, stop below it.
        let mut records = vec![
            rec("CDS", 1802, 2950, '-'),
            rec("exon", 1802, 2950, '-'),
            rec("stop_codon", 1799, 1801, '-'),
        ];
        include_stop_codon_in_cds(&mut records);
        assert_eq!(
            find(&records, "CDS"),
            (1799, 2950),
            "CDS start == stop start"
        );
        assert_eq!(find(&records, "exon"), (1799, 2950), "exon follows CDS");
    }

    #[test]
    fn multi_exon_plus_extends_only_terminal_segment() {
        let mut records = vec![
            rec("CDS", 100, 200, '+'),
            rec("exon", 100, 200, '+'),
            rec("CDS", 300, 450, '+'),
            rec("exon", 300, 450, '+'),
            rec("stop_codon", 451, 453, '+'),
        ];
        include_stop_codon_in_cds(&mut records);
        // First segment untouched; last (max-end) segment extended by the codon.
        assert_eq!(records[0].end, 200);
        assert_eq!(records[1].end, 200);
        assert_eq!(records[2].end, 453);
        assert_eq!(records[3].end, 453);
    }

    #[test]
    fn idempotent_when_cds_already_includes_stop() {
        // Second pass over already-corrected output must be a no-op: the stop
        // now lies inside the CDS, so the adjacency test fails.
        let mut records = vec![
            rec("CDS", 1802, 2953, '+'),
            rec("exon", 1802, 2953, '+'),
            rec("stop_codon", 2951, 2953, '+'),
        ];
        include_stop_codon_in_cds(&mut records);
        assert_eq!(find(&records, "CDS"), (1802, 2953));
        assert_eq!(find(&records, "exon"), (1802, 2953));
    }

    #[test]
    fn no_stop_codon_leaves_partial_gene_untouched() {
        let mut records = vec![rec("CDS", 100, 200, '+'), rec("exon", 100, 200, '+')];
        include_stop_codon_in_cds(&mut records);
        assert_eq!(find(&records, "CDS"), (100, 200));
    }

    #[test]
    fn non_adjacent_stop_is_not_folded_in() {
        // A stop codon that is not flush against the terminal CDS (e.g. a split
        // stop across an intron, or a malformed model) is left alone rather
        // than swallowing the intervening bases.
        let mut records = vec![
            rec("CDS", 100, 200, '+'),
            rec("exon", 100, 200, '+'),
            rec("stop_codon", 260, 262, '+'),
        ];
        include_stop_codon_in_cds(&mut records);
        assert_eq!(find(&records, "CDS"), (100, 200));
    }
}
