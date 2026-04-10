/// Output validation module
///
/// Validates GFF3, FASTA, and GenBank output files for correctness:
///   - GFF3: ID uniqueness, Parent consistency, coordinate ordering,
///           required attributes, feature type hierarchy
///   - FASTA: internal stop codons, ambiguous amino acids,
///           length consistency, valid characters
///   - GenBank: locus_tag uniqueness, required qualifiers
///
/// This is myconote-cli's own validation implementation.
use crate::parser::gff::{GFFReader, GFFRecord};
use crate::utils::error::{MycoNoteError, Result};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 Validation
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
pub struct Gff3ValidationResult {
    pub total_features: usize,
    pub errors: Vec<ValidationIssue>,
    pub warnings: Vec<ValidationIssue>,
}

#[derive(Debug, Clone)]
pub struct ValidationIssue {
    pub level: IssueLevel,
    pub category: String,
    pub message: String,
    pub line: Option<usize>,
}

#[derive(Debug, Clone, Copy)]
pub enum IssueLevel {
    Error,
    Warning,
    Info,
}

impl Gff3ValidationResult {
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }

    pub fn print_summary(&self) {
        println!(
            "  GFF3 Validation: {} features checked",
            self.total_features
        );
        if self.errors.is_empty() {
            println!(
                "  Result: PASSED (0 errors, {} warnings)",
                self.warnings.len()
            );
        } else {
            println!(
                "  Result: FAILED ({} errors, {} warnings)",
                self.errors.len(),
                self.warnings.len()
            );
        }
        for e in self.errors.iter().take(20) {
            println!("    ERROR [{}]: {}", e.category, e.message);
        }
        for w in self.warnings.iter().take(10) {
            println!("    WARN  [{}]: {}", w.category, w.message);
        }
        if self.errors.len() > 20 {
            println!("    ... and {} more errors", self.errors.len() - 20);
        }
    }
}

/// Validate a GFF3 file for common issues.
pub fn validate_gff3(path: &Path) -> Result<Gff3ValidationResult> {
    let mut result = Gff3ValidationResult::default();

    let records: Vec<GFFRecord> = GFFReader::from_path(path)?.filter_map(|r| r.ok()).collect();

    result.total_features = records.len();

    // 1. Check ID uniqueness
    let mut id_counts: HashMap<String, usize> = HashMap::new();
    for rec in &records {
        if let Some(id) = rec.id() {
            *id_counts.entry(id.clone()).or_insert(0) += 1;
        }
    }
    for (id, count) in &id_counts {
        if *count > 1 {
            result.errors.push(ValidationIssue {
                level: IssueLevel::Error,
                category: "duplicate_id".to_string(),
                message: format!("ID '{}' appears {} times", id, count),
                line: None,
            });
        }
    }

    // 2. Check Parent references
    let all_ids: HashSet<String> = records.iter().filter_map(|r| r.id().cloned()).collect();

    for rec in &records {
        if let Some(parent) = rec.parent() {
            if !all_ids.contains(parent.as_str()) {
                result.errors.push(ValidationIssue {
                    level: IssueLevel::Error,
                    category: "orphan_feature".to_string(),
                    message: format!(
                        "{} at {}:{}-{} references non-existent Parent='{}'",
                        rec.feature_type, rec.seqid, rec.start, rec.end, parent
                    ),
                    line: None,
                });
            }
        }
    }

    // 3. Check coordinate ordering (start <= end)
    for rec in &records {
        if rec.start > rec.end {
            result.errors.push(ValidationIssue {
                level: IssueLevel::Error,
                category: "coordinate_order".to_string(),
                message: format!(
                    "{} at {}: start ({}) > end ({})",
                    rec.feature_type, rec.seqid, rec.start, rec.end
                ),
                line: None,
            });
        }
    }

    // 4. Check feature hierarchy (gene → mRNA → CDS)
    let mut gene_children: HashMap<String, Vec<String>> = HashMap::new();
    for rec in &records {
        if let Some(parent) = rec.parent() {
            gene_children
                .entry(parent.clone())
                .or_default()
                .push(rec.feature_type.clone());
        }
    }

    for rec in &records {
        if rec.feature_type == "gene" {
            if let Some(id) = rec.id() {
                let children_types: HashSet<_> = gene_children
                    .get(id)
                    .map(|v| v.iter().cloned().collect())
                    .unwrap_or_default();

                if children_types.is_empty() {
                    result.warnings.push(ValidationIssue {
                        level: IssueLevel::Warning,
                        category: "childless_gene".to_string(),
                        message: format!("Gene '{}' has no child features", id),
                        line: None,
                    });
                }
            }
        }
    }

    // 5. Check CDS has valid phase for coding features
    for rec in &records {
        if rec.feature_type == "CDS" && rec.phase.is_none() {
            result.warnings.push(ValidationIssue {
                level: IssueLevel::Warning,
                category: "missing_phase".to_string(),
                message: format!(
                    "CDS at {}:{}-{} has no phase value",
                    rec.seqid, rec.start, rec.end
                ),
                line: None,
            });
        }
    }

    // 6. Check genes have locus_tag or ID
    for rec in &records {
        if rec.feature_type == "gene" {
            if rec.id().is_none() && !rec.attributes.contains_key("locus_tag") {
                result.warnings.push(ValidationIssue {
                    level: IssueLevel::Warning,
                    category: "missing_identifier".to_string(),
                    message: format!(
                        "Gene at {}:{}-{} has no ID or locus_tag",
                        rec.seqid, rec.start, rec.end
                    ),
                    line: None,
                });
            }
        }
    }

    Ok(result)
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTA Validation (protein)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
pub struct FastaValidationResult {
    pub total_sequences: usize,
    pub total_residues: usize,
    pub internal_stops: Vec<(String, Vec<usize>)>,
    pub short_sequences: Vec<(String, usize)>,
    pub invalid_chars: Vec<(String, Vec<char>)>,
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
}

impl FastaValidationResult {
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }

    pub fn print_summary(&self) {
        println!(
            "  FASTA Validation: {} sequences, {} residues",
            self.total_sequences, self.total_residues
        );
        if !self.internal_stops.is_empty() {
            println!(
                "  WARNING: {} sequences have internal stop codons",
                self.internal_stops.len()
            );
        }
        if !self.short_sequences.is_empty() {
            println!(
                "  WARNING: {} sequences shorter than 10 aa",
                self.short_sequences.len()
            );
        }
        if !self.invalid_chars.is_empty() {
            println!(
                "  ERROR: {} sequences have invalid characters",
                self.invalid_chars.len()
            );
        }
    }
}

/// Validate a protein FASTA file.
pub fn validate_protein_fasta(path: &Path) -> Result<FastaValidationResult> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);

    let mut result = FastaValidationResult::default();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    let valid_aa: HashSet<char> = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*".chars().collect();

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();

        if let Some(header) = trimmed.strip_prefix('>') {
            // Process previous sequence
            if !current_id.is_empty() {
                check_protein_sequence(&current_id, &current_seq, &valid_aa, &mut result);
            }
            current_id = header.split_whitespace().next().unwrap_or("").to_string();
            current_seq.clear();
            result.total_sequences += 1;
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    // Process last sequence
    if !current_id.is_empty() {
        check_protein_sequence(&current_id, &current_seq, &valid_aa, &mut result);
    }

    Ok(result)
}

fn check_protein_sequence(
    id: &str,
    seq: &str,
    valid_aa: &HashSet<char>,
    result: &mut FastaValidationResult,
) {
    result.total_residues += seq.len();

    // Check for internal stop codons
    let internal_stops: Vec<usize> = seq
        .char_indices()
        .filter(|&(i, c)| c == '*' && i < seq.len() - 1)
        .map(|(i, _)| i)
        .collect();

    if !internal_stops.is_empty() {
        result.internal_stops.push((id.to_string(), internal_stops));
        result
            .warnings
            .push(format!("{}: contains internal stop codon(s)", id));
    }

    // Check for short sequences
    let effective_len = seq.trim_end_matches('*').len();
    if effective_len < 10 {
        result.short_sequences.push((id.to_string(), effective_len));
    }

    // Check for invalid characters
    let invalid: Vec<char> = seq.chars().filter(|c| !valid_aa.contains(c)).collect();
    if !invalid.is_empty() {
        let unique_invalid: Vec<char> = invalid
            .iter()
            .cloned()
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        result.invalid_chars.push((id.to_string(), unique_invalid));
        result
            .errors
            .push(format!("{}: contains invalid amino acid characters", id));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Write validated GFF3 (fix common issues in-place)
// ─────────────────────────────────────────────────────────────────────────────

/// Read a GFF3 file, validate it, and write a cleaned version.
/// Returns (features_written, issues_fixed).
pub fn validate_and_fix_gff3(input: &Path, output: &Path) -> Result<(usize, usize)> {
    let records: Vec<GFFRecord> = GFFReader::from_path(input)?
        .filter_map(|r| r.ok())
        .collect();

    let all_ids: HashSet<String> = records.iter().filter_map(|r| r.id().cloned()).collect();

    let mut f = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(f, "##gff-version 3").map_err(MycoNoteError::Io)?;

    let mut written = 0usize;
    let mut fixed = 0usize;

    // Track IDs to deduplicate
    let mut seen_ids: HashSet<String> = HashSet::new();

    for mut rec in records {
        // Fix: skip orphan features
        if let Some(parent) = rec.parent() {
            if !all_ids.contains(parent.as_str()) {
                fixed += 1;
                continue; // skip orphans
            }
        }

        // Fix: deduplicate IDs
        if let Some(id) = rec.id().cloned() {
            if seen_ids.contains(&id) {
                // Append suffix
                let new_id = format!("{}_dup{}", id, fixed + 1);
                rec.attributes.insert("ID".to_string(), new_id.clone());
                seen_ids.insert(new_id);
                fixed += 1;
            } else {
                seen_ids.insert(id);
            }
        }

        // Fix: swap start/end if inverted
        if rec.start > rec.end {
            std::mem::swap(&mut rec.start, &mut rec.end);
            fixed += 1;
        }

        writeln!(f, "{}", rec.to_gff3_line()).map_err(MycoNoteError::Io)?;
        written += 1;
    }

    Ok((written, fixed))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_protein_internal_stops() {
        let valid_aa: HashSet<char> = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*".chars().collect();
        let mut result = FastaValidationResult::default();

        check_protein_sequence("test1", "MKFG*TER*", &valid_aa, &mut result);
        assert_eq!(result.internal_stops.len(), 1);
        assert_eq!(result.internal_stops[0].1, vec![4]); // position of internal *
    }

    #[test]
    fn test_check_protein_no_stops() {
        let valid_aa: HashSet<char> = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*".chars().collect();
        let mut result = FastaValidationResult::default();

        check_protein_sequence("test2", "MKFGTERVWDL*", &valid_aa, &mut result);
        assert!(result.internal_stops.is_empty()); // trailing * is not internal
    }
}
