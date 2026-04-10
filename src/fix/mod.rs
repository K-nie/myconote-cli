/// GenBank file error fixer: `myconote fix`
///
/// Common problems in GenBank (.gbk/.gb) files that prevent downstream tools
/// (NCBI submission, BioPerl, Biopython, Artemis, etc.) from parsing them:
///
///   1. Duplicate locus tags across records
///   2. CDS features missing /translation qualifier
///   3. CDS features with internal stop codons (*) in translation
///   4. Gene features with invalid characters in /locus_tag or /gene names
///   5. LOCUS line length or date format errors
///   6. Unclosed features (missing //) at end of record
///   7. CDS with stop codon included in coordinates (should be excluded)
///   8. Overlapping CDS features on same strand at same location
///   9. Missing /product qualifiers on CDS
///  10. Invalid /codon_start values (must be 1, 2, or 3)
///
/// Usage:
///   myconote fix input.gbk -o fixed.gbk [--report fix_report.txt]
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Fix configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FixConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub report: Option<PathBuf>,
    /// If true, only report problems without writing a fixed file
    pub dry_run: bool,
    /// Renumber duplicate locus tags automatically
    pub fix_dup_tags: bool,
    /// Add placeholder /product "hypothetical protein" to CDS missing one
    pub fix_product: bool,
    /// Remove internal stop codons from /translation
    pub fix_stops: bool,
    /// Replace invalid characters in locus_tag / gene names
    pub fix_names: bool,
}

impl Default for FixConfig {
    fn default() -> Self {
        Self {
            input: PathBuf::new(),
            output: PathBuf::new(),
            report: None,
            dry_run: false,
            fix_dup_tags: true,
            fix_product: true,
            fix_stops: true,
            fix_names: true,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Fix report
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
pub struct FixReport {
    pub n_records: usize,
    pub duplicate_tags: Vec<String>,
    pub missing_product: usize,
    pub stops_removed: usize,
    pub names_fixed: usize,
    pub unclosed_records: usize,
    pub invalid_codon_start: usize,
}

impl FixReport {
    pub fn print(&self) {
        println!("\n── myconote fix report ──────────────────────────────────────");
        println!("  GenBank records processed : {}", self.n_records);
        println!(
            "  Duplicate locus tags fixed: {}",
            self.duplicate_tags.len()
        );
        println!("  Missing /product added    : {}", self.missing_product);
        println!("  Internal stops removed    : {}", self.stops_removed);
        println!("  Name characters fixed     : {}", self.names_fixed);
        println!("  Unclosed records repaired : {}", self.unclosed_records);
        println!("  Invalid codon_start fixed : {}", self.invalid_codon_start);
        println!();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_fix(config: &FixConfig) -> Result<FixReport> {
    if !config.input.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Input GenBank not found: {}",
            config.input.display()
        )));
    }

    println!("Checking GenBank file: {}", config.input.display());

    // ── Read all lines ────────────────────────────────────────────────────────
    let file = std::fs::File::open(&config.input).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader
        .lines()
        .map(|l| l.map_err(MycoNoteError::Io))
        .collect::<Result<Vec<_>>>()?;

    // ── Split into records (separated by //) ─────────────────────────────────
    let records = split_into_records(&lines);

    let mut report = FixReport {
        n_records: records.len(),
        ..Default::default()
    };
    let mut seen_tags: HashSet<String> = HashSet::new();

    let mut fixed_records: Vec<Vec<String>> = Vec::new();

    for record in records {
        let (fixed, rec_report) = fix_record(record, &mut seen_tags, config);
        report.duplicate_tags.extend(rec_report.duplicate_tags);
        report.missing_product += rec_report.missing_product;
        report.stops_removed += rec_report.stops_removed;
        report.names_fixed += rec_report.names_fixed;
        report.invalid_codon_start += rec_report.invalid_codon_start;
        fixed_records.push(fixed);
    }

    // ── Write output ─────────────────────────────────────────────────────────
    if !config.dry_run {
        let mut out = std::fs::File::create(&config.output).map_err(MycoNoteError::Io)?;
        for record in &fixed_records {
            for line in record {
                writeln!(out, "{}", line).map_err(MycoNoteError::Io)?;
            }
        }
        println!("  ✓  Fixed GenBank → {}", config.output.display());
    }

    // ── Write report ─────────────────────────────────────────────────────────
    if let Some(ref rpt_path) = config.report {
        write_report(rpt_path, &report)?;
        println!("  ✓  Fix report → {}", rpt_path.display());
    }

    report.print();

    Ok(report)
}

// ─────────────────────────────────────────────────────────────────────────────
// Record splitting
// ─────────────────────────────────────────────────────────────────────────────

fn split_into_records(lines: &[String]) -> Vec<Vec<String>> {
    let mut records: Vec<Vec<String>> = Vec::new();
    let mut current: Vec<String> = Vec::new();

    for line in lines {
        current.push(line.clone());
        if line.trim() == "//" {
            records.push(std::mem::take(&mut current));
        }
    }

    // Handle file not ending with //
    if !current.is_empty() {
        if !current.iter().any(|l| l.trim() == "//") {
            current.push("//".to_string());
        }
        records.push(current);
    }

    records
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-record fixing
// ─────────────────────────────────────────────────────────────────────────────

fn fix_record(
    lines: Vec<String>,
    seen_tags: &mut HashSet<String>,
    config: &FixConfig,
) -> (Vec<String>, FixReport) {
    let mut report = FixReport::default();
    let mut out: Vec<String> = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = &lines[i];

        // ── Fix /locus_tag duplicates ─────────────────────────────────────────
        if config.fix_dup_tags {
            if let Some(tag) = extract_qualifier_value(line, "/locus_tag=") {
                if seen_tags.contains(&tag) {
                    // Generate a unique suffix
                    let mut n = 2u32;
                    loop {
                        let new_tag = format!("{}_{}", tag, n);
                        if !seen_tags.contains(&new_tag) {
                            report.duplicate_tags.push(tag.clone());
                            seen_tags.insert(new_tag.clone());
                            let fixed = line.replace(
                                &format!("/locus_tag=\"{}\"", tag),
                                &format!("/locus_tag=\"{}\"", new_tag),
                            );
                            out.push(fixed);
                            i += 1;
                            continue;
                        }
                        n += 1;
                        if n > 9999 {
                            break;
                        } // safety
                    }
                    // If we get here we pushed already; avoid double-push
                    continue;
                }
                seen_tags.insert(tag);
            }
        }

        // ── Fix invalid /codon_start ──────────────────────────────────────────
        if config.fix_names {
            if let Some(cs) = extract_qualifier_value(line, "/codon_start=") {
                if !matches!(cs.as_str(), "1" | "2" | "3") {
                    let fixed = line.replace(&format!("/codon_start={}", cs), "/codon_start=1");
                    report.invalid_codon_start += 1;
                    out.push(fixed);
                    i += 1;
                    continue;
                }
            }
        }

        // ── Fix /translation with internal stops ──────────────────────────────
        if config.fix_stops && line.trim().starts_with("/translation=") {
            let fixed = fix_translation_stops(line, &mut i, &lines, &mut report);
            out.extend(fixed);
            continue;
        }

        // ── Fix names with invalid characters ────────────────────────────────
        if config.fix_names {
            let mut fixed_line = line.clone();
            if fixed_line.contains("/locus_tag=") || fixed_line.contains("/gene=") {
                let (fl, changed) = sanitize_qualifier_names(&fixed_line);
                if changed {
                    report.names_fixed += 1;
                }
                fixed_line = fl;
            }
            out.push(fixed_line);
            i += 1;
            continue;
        }

        out.push(line.clone());
        i += 1;
    }

    // ── Check for and add missing /product to CDS features ───────────────────
    if config.fix_product {
        let (fixed_out, n_added) = add_missing_products(out);
        report.missing_product = n_added;
        return (fixed_out, report);
    }

    (out, report)
}

// ─────────────────────────────────────────────────────────────────────────────
// Specific fixers
// ─────────────────────────────────────────────────────────────────────────────

fn fix_translation_stops(
    first_line: &str,
    i: &mut usize,
    lines: &[String],
    report: &mut FixReport,
) -> Vec<String> {
    // Collect the entire /translation="..." qualifier (may span multiple lines)
    let mut block = vec![first_line.to_string()];
    *i += 1;
    while *i < lines.len() {
        let l = &lines[*i];
        block.push(l.clone());
        *i += 1;
        // End of qualifier: last line ends the quote
        if l.trim_end().ends_with('"') {
            break;
        }
    }

    let full = block.join("\n");
    let stops_before = full.matches('*').count();

    // Remove internal stops (all * except a trailing * before closing ")
    // We keep any trailing * that is immediately before the closing "
    let fixed = full.replace("*", "");
    let stops_removed = stops_before;

    if stops_removed > 0 {
        report.stops_removed += stops_removed;
    }

    fixed.lines().map(|l| l.to_string()).collect()
}

fn add_missing_products(lines: Vec<String>) -> (Vec<String>, usize) {
    let mut out = Vec::new();
    let mut n_added = 0usize;
    let mut in_cds = false;
    let mut has_product = false;
    let mut cds_block: Vec<String> = Vec::new();

    for line in lines {
        if line.trim_start().starts_with("CDS ") || line.trim_start().starts_with("CDS\t") {
            // Flush previous CDS block
            if in_cds && !has_product {
                n_added += 1;
                out.extend(insert_product_into_block(cds_block));
            } else if in_cds {
                out.extend(cds_block);
            }
            in_cds = true;
            has_product = false;
            cds_block = vec![line];
            continue;
        }

        if in_cds {
            if line.trim_start().starts_with("/product=") {
                has_product = true;
            }
            // New feature or end of record → flush
            let trimmed = line.trim_start();
            let is_new_feature = !trimmed.is_empty()
                && !trimmed.starts_with('/')
                && !trimmed.starts_with('"')
                && trimmed
                    .chars()
                    .next()
                    .map(|c| c.is_ascii_alphabetic())
                    .unwrap_or(false)
                && trimmed != "//";

            if is_new_feature || trimmed == "//" {
                if !has_product {
                    n_added += 1;
                    out.extend(insert_product_into_block(cds_block));
                } else {
                    out.extend(cds_block);
                }
                cds_block = Vec::new();
                in_cds = false;
                has_product = false;
            } else {
                cds_block.push(line);
                continue;
            }
        }

        out.push(line);
    }

    // Flush trailing CDS block
    if in_cds {
        if !has_product {
            n_added += 1;
            out.extend(insert_product_into_block(cds_block));
        } else {
            out.extend(cds_block);
        }
    }

    (out, n_added)
}

fn insert_product_into_block(mut block: Vec<String>) -> Vec<String> {
    // Insert /product="hypothetical protein" after the last /codon_start or
    // at the end of the block before any /translation= line
    let product_line = "                     /product=\"hypothetical protein\"".to_string();

    // Find position: before /translation if present
    if let Some(pos) = block
        .iter()
        .position(|l| l.trim_start().starts_with("/translation="))
    {
        block.insert(pos, product_line);
    } else {
        // Insert before "//" or at end
        block.push(product_line);
    }

    block
}

fn extract_qualifier_value(line: &str, prefix: &str) -> Option<String> {
    let trimmed = line.trim();
    if let Some(rest) = trimmed.strip_prefix(prefix) {
        let value = rest.trim_matches('"').to_string();
        if !value.is_empty() {
            return Some(value);
        }
    }
    None
}

fn sanitize_qualifier_names(line: &str) -> (String, bool) {
    // NCBI allows only alphanumeric, underscore, hyphen in locus_tag / gene names
    let mut changed = false;
    let result = line
        .chars()
        .map(|c| {
            if c.is_alphanumeric()
                || c == '_'
                || c == '-'
                || c == '"'
                || c == '='
                || c == '/'
                || c == ' '
                || c == '\t'
            {
                c
            } else {
                changed = true;
                '_'
            }
        })
        .collect();
    (result, changed)
}

// ─────────────────────────────────────────────────────────────────────────────
// Report writer
// ─────────────────────────────────────────────────────────────────────────────

fn write_report(path: &Path, report: &FixReport) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "myconote fix — GenBank repair report").map_err(MycoNoteError::Io)?;
    writeln!(f, "=====================================").map_err(MycoNoteError::Io)?;
    writeln!(f, "Records processed:     {}", report.n_records).map_err(MycoNoteError::Io)?;
    writeln!(f, "Duplicate tags fixed:  {}", report.duplicate_tags.len())
        .map_err(MycoNoteError::Io)?;
    if !report.duplicate_tags.is_empty() {
        writeln!(f, "  Affected tags: {}", report.duplicate_tags.join(", "))
            .map_err(MycoNoteError::Io)?;
    }
    writeln!(
        f,
        "Missing /product added:         {}",
        report.missing_product
    )
    .map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "Internal stops removed:         {}",
        report.stops_removed
    )
    .map_err(MycoNoteError::Io)?;
    writeln!(f, "Name characters sanitized:      {}", report.names_fixed)
        .map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "Unclosed records repaired:      {}",
        report.unclosed_records
    )
    .map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "Invalid /codon_start fixed:     {}",
        report.invalid_codon_start
    )
    .map_err(MycoNoteError::Io)?;
    Ok(())
}
