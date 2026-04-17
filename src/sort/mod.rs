/// Genome FASTA pre-processing: sort + rename contigs
///
/// Equivalent to `funannotate sort`:
///   - Sorts contigs by length (longest first)
///   - Removes contigs shorter than --min-length
///   - Renames headers to clean, sequential IDs (scaffold_1, scaffold_2, …)
///   - Optionally strips any text after the first whitespace in the header
///   - Reports a rename table so downstream GFF3 files can be lifted over
///
/// This should be run before `myconote mask` to ensure consistent,
/// clean sequence identifiers throughout the whole pipeline.
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct SortConfig {
    /// Input FASTA (raw assembly)
    pub input: PathBuf,
    /// Output FASTA (sorted, renamed)
    pub output: PathBuf,
    /// Discard contigs shorter than this (bp). Default: 0 (keep all)
    pub min_length: usize,
    /// Prefix for renamed headers. Default: "scaffold"
    pub prefix: String,
    /// Strip everything after first whitespace in the original header
    pub strip_desc: bool,
    /// Write old→new rename table to this path (TSV)
    pub rename_table: Option<PathBuf>,
    /// If true, sort by name instead of by length
    pub sort_by_name: bool,
}

impl Default for SortConfig {
    fn default() -> Self {
        Self {
            input: PathBuf::new(),
            output: PathBuf::new(),
            min_length: 0,
            prefix: "scaffold".to_string(),
            strip_desc: true,
            rename_table: None,
            sort_by_name: false,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal: read FASTA into (header, sequence) pairs
// ─────────────────────────────────────────────────────────────────────────────

struct FastaRecord {
    /// Full original header line (without '>')
    original_header: String,
    /// Bare ID (first token of header)
    original_id: String,
    sequence: String,
}

fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);

    let mut records: Vec<FastaRecord> = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        if line.starts_with('>') {
            if !current_header.is_empty() {
                let id = current_header
                    .split_whitespace()
                    .next()
                    .unwrap_or(&current_header)
                    .to_string();
                records.push(FastaRecord {
                    original_header: current_header.clone(),
                    original_id: id,
                    sequence: current_seq.clone(),
                });
                current_seq.clear();
            }
            current_header = line[1..].to_string();
        } else {
            current_seq.push_str(line.trim());
        }
    }

    // Push final record
    if !current_header.is_empty() {
        let id = current_header
            .split_whitespace()
            .next()
            .unwrap_or(&current_header)
            .to_string();
        records.push(FastaRecord {
            original_header: current_header,
            original_id: id,
            sequence: current_seq,
        });
    }

    Ok(records)
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Run the sort/rename pipeline and return the rename map (old_id → new_id).
pub fn run_sort(config: &SortConfig) -> Result<HashMap<String, String>> {
    // ── Read ──────────────────────────────────────────────────────────────
    let mut records = read_fasta(&config.input)?;

    // ── Filter ────────────────────────────────────────────────────────────
    let before = records.len();
    records.retain(|r| r.sequence.len() >= config.min_length);
    let removed = before - records.len();

    // ── Sort ──────────────────────────────────────────────────────────────
    if config.sort_by_name {
        records.sort_by(|a, b| a.original_id.cmp(&b.original_id));
    } else {
        // Longest first
        records.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));
    }

    // ── Build rename map ──────────────────────────────────────────────────
    let pad = records.len().to_string().len().max(3);
    let mut rename_map: HashMap<String, String> = HashMap::new();
    for (i, rec) in records.iter().enumerate() {
        let new_id = format!("{}_{:0>width$}", config.prefix, i + 1, width = pad);
        rename_map.insert(rec.original_id.clone(), new_id);
    }

    // ── Write renamed FASTA ───────────────────────────────────────────────
    // Create parent directory implicitly — users routinely pass paths like
    // `my_out_dir/sorted.fa` without running `mkdir -p my_out_dir` first,
    // and the resulting "IO error: No such file or directory" is cryptic.
    if let Some(parent) = config.output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(MycoNoteError::Io)?;
        }
    }
    let mut out = std::fs::File::create(&config.output).map_err(MycoNoteError::Io)?;
    for rec in &records {
        let new_id = &rename_map[&rec.original_id];
        if config.strip_desc {
            writeln!(out, ">{}", new_id).map_err(MycoNoteError::Io)?;
        } else {
            // Keep description but update the ID portion
            let desc_rest: String = rec
                .original_header
                .splitn(2, char::is_whitespace)
                .nth(1)
                .map(|s| format!(" {}", s))
                .unwrap_or_default();
            writeln!(out, ">{}{}", new_id, desc_rest).map_err(MycoNoteError::Io)?;
        }
        // Write sequence in 60-char lines
        for chunk in rec.sequence.as_bytes().chunks(60) {
            out.write_all(chunk).map_err(MycoNoteError::Io)?;
            writeln!(out).map_err(MycoNoteError::Io)?;
        }
    }

    // ── Optionally write rename table ─────────────────────────────────────
    if let Some(ref table_path) = config.rename_table {
        if let Some(parent) = table_path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent).map_err(MycoNoteError::Io)?;
            }
        }
        let mut tbl = std::fs::File::create(table_path).map_err(MycoNoteError::Io)?;
        writeln!(tbl, "original_id\tnew_id\tlength").map_err(MycoNoteError::Io)?;
        for rec in &records {
            let new_id = &rename_map[&rec.original_id];
            writeln!(
                tbl,
                "{}\t{}\t{}",
                rec.original_id,
                new_id,
                rec.sequence.len()
            )
            .map_err(MycoNoteError::Io)?;
        }
    }

    // ── Report ────────────────────────────────────────────────────────────
    let total_bp: usize = records.iter().map(|r| r.sequence.len()).sum();
    let longest = records.first().map(|r| r.sequence.len()).unwrap_or(0);
    let shortest = records.last().map(|r| r.sequence.len()).unwrap_or(0);

    println!("  Contigs kept:   {}", records.len());
    if removed > 0 {
        println!(
            "  Contigs removed (< {} bp): {}",
            config.min_length, removed
        );
    }
    println!("  Total assembly: {} bp", total_bp);
    println!("  Longest contig: {} bp", longest);
    println!("  Shortest kept:  {} bp", shortest);
    println!("  Output:         {}", config.output.display());
    if let Some(ref t) = config.rename_table {
        println!("  Rename table:   {}", t.display());
    }

    Ok(rename_map)
}

/// Apply a rename map to an existing GFF3 file, updating seqid (column 1).
/// Useful when you need to lift over annotations after sorting the assembly.
pub fn liftover_gff3(
    gff_input: &Path,
    gff_output: &Path,
    rename_map: &HashMap<String, String>,
) -> Result<usize> {
    let inp = std::fs::File::open(gff_input).map_err(MycoNoteError::Io)?;
    if let Some(parent) = gff_output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(MycoNoteError::Io)?;
        }
    }
    let mut out = std::fs::File::create(gff_output).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(inp);
    let mut updated = 0usize;

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() {
            writeln!(out, "{}", line).map_err(MycoNoteError::Io)?;
            continue;
        }
        let mut fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 9 {
            writeln!(out, "{}", line).map_err(MycoNoteError::Io)?;
            continue;
        }
        if let Some(new_id) = rename_map.get(fields[0]) {
            fields[0] = new_id;
            updated += 1;
        }
        writeln!(out, "{}", fields.join("\t")).map_err(MycoNoteError::Io)?;
    }

    Ok(updated)
}
