/// FASTA parser
///
/// Reads multi-FASTA files into `FastaRecord` structs.
/// Provides both sequential iteration and random-access via a
/// `HashMap<String, FastaRecord>` keyed on the bare sequence ID
/// (first whitespace-delimited token of the header line).

use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// Data types
// ─────────────────────────────────────────────────────────────────────────────

/// A single FASTA record.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Full header line (excluding the leading `>`).
    pub header: String,
    /// Bare sequence identifier — first whitespace-delimited token of `header`.
    pub id: String,
    /// Upper-case DNA/RNA/protein sequence, all whitespace stripped.
    pub sequence: String,
}

impl FastaRecord {
    /// Length of the sequence in bases / residues.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Return a sub-sequence using **1-based inclusive** coordinates (GFF3 style).
    /// Clamps to sequence bounds; returns an empty string for out-of-range slices.
    pub fn subsequence(&self, start: u64, end: u64) -> &str {
        if start == 0 || start > end {
            return "";
        }
        let s = (start as usize).saturating_sub(1);
        let e = (end as usize).min(self.sequence.len());
        if s >= e {
            ""
        } else {
            &self.sequence[s..e]
        }
    }

    /// Reverse-complement of the whole sequence (DNA only).
    pub fn reverse_complement(&self) -> String {
        self.sequence
            .chars()
            .rev()
            .map(complement_base)
            .collect()
    }
}

fn complement_base(c: char) -> char {
    match c {
        'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G',
        'a' => 't', 't' => 'a', 'g' => 'c', 'c' => 'g',
        'N' => 'N', 'n' => 'n',
        'U' => 'A', 'u' => 'a', // RNA
        other => other,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Iterator-based reader
// ─────────────────────────────────────────────────────────────────────────────

/// Streaming FASTA reader — yields one `FastaRecord` at a time without
/// loading the whole file into memory.
pub struct FastaReader {
    lines: std::io::Lines<BufReader<File>>,
    /// Pending header from the previous `>` line (header for the *next* record).
    pending_header: Option<String>,
}

impl FastaReader {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref()).map_err(|e| {
            MycoNoteError::Io(e)
        })?;
        Ok(FastaReader {
            lines: BufReader::new(file).lines(),
            pending_header: None,
        })
    }
}

impl Iterator for FastaReader {
    type Item = Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut seq_lines: Vec<String> = Vec::new();
        let mut current_header: Option<String> = self.pending_header.take();

        loop {
            match self.lines.next() {
                None => {
                    // EOF
                    return if let Some(header) = current_header {
                        Some(Ok(build_record(header, seq_lines)))
                    } else {
                        None
                    };
                }
                Some(Err(e)) => return Some(Err(MycoNoteError::Io(e))),
                Some(Ok(line)) => {
                    let trimmed = line.trim_end();
                    if trimmed.is_empty() {
                        continue;
                    }
                    if trimmed.starts_with('>') {
                        let new_header = trimmed[1..].to_string();
                        if let Some(header) = current_header {
                            // Save the new header for next call
                            self.pending_header = Some(new_header);
                            return Some(Ok(build_record(header, seq_lines)));
                        } else {
                            current_header = Some(new_header);
                        }
                    } else if current_header.is_some() {
                        // Sequence line — strip all whitespace and uppercase
                        seq_lines.push(
                            trimmed.chars()
                                .filter(|c| !c.is_whitespace())
                                .flat_map(|c| c.to_uppercase())
                                .collect()
                        );
                    }
                    // Lines before any header are silently ignored.
                }
            }
        }
    }
}

fn build_record(header: String, seq_lines: Vec<String>) -> FastaRecord {
    let id = header
        .split_whitespace()
        .next()
        .unwrap_or(&header)
        .to_string();
    let sequence = seq_lines.concat();
    FastaRecord { header, id, sequence }
}

// ─────────────────────────────────────────────────────────────────────────────
// Convenience load functions
// ─────────────────────────────────────────────────────────────────────────────

/// Load all records from a FASTA file into a `Vec`.
pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<FastaRecord>> {
    FastaReader::from_path(path)?.collect()
}

/// Load all records into a `HashMap` keyed by bare sequence ID.
///
/// If two sequences share the same bare ID, the later one wins
/// (this mirrors `biopython` behaviour).
pub fn read_fasta_index<P: AsRef<Path>>(
    path: P,
) -> Result<HashMap<String, FastaRecord>> {
    let records = read_fasta(path)?;
    let mut map = HashMap::with_capacity(records.len());
    for rec in records {
        map.insert(rec.id.clone(), rec);
    }
    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// Free-function utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Reverse-complement a DNA string (free function version).
pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement_base).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = ">seq1 some description\nATGCATGC\nATGC\n>seq2\nNNNNNN\n";

    fn write_temp(content: &str) -> tempfile::NamedTempFile {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn test_basic_parse() {
        let f = write_temp(SAMPLE);
        let records = read_fasta(f.path()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].header, "seq1 some description");
        assert_eq!(records[0].sequence, "ATGCATGCATGC");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].sequence, "NNNNNN");
    }

    #[test]
    fn test_index() {
        let f = write_temp(SAMPLE);
        let idx = read_fasta_index(f.path()).unwrap();
        assert!(idx.contains_key("seq1"));
        assert!(idx.contains_key("seq2"));
    }

    #[test]
    fn test_subsequence() {
        let f = write_temp(SAMPLE);
        let records = read_fasta(f.path()).unwrap();
        // 1-based: positions 1–4 → "ATGC"
        assert_eq!(records[0].subsequence(1, 4), "ATGC");
        // positions 5–8
        assert_eq!(records[0].subsequence(5, 8), "ATGC");
    }

    #[test]
    fn test_reverse_complement() {
        let f = write_temp(">s\nATGC\n");
        let records = read_fasta(f.path()).unwrap();
        assert_eq!(records[0].reverse_complement(), "GCAT");
    }
}
