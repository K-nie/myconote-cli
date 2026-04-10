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
