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
    /// Rename contigs to clean sequential IDs (`<prefix>_N`). Default `false`:
    /// the original FASTA seqids are preserved end-to-end so downstream
    /// predict/annotate/GFF3 carry the user's own accessions and no lift-back
    /// table is needed. Opt in with `--ncbi-clean` / `--rename-contigs` when
    /// short, submission-safe names are required (e.g. before `submit`).
    pub rename_contigs: bool,
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
            rename_contigs: false,
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
    // Default (`rename_contigs == false`): identity map — every contig keeps
    // its original seqid, so predict/annotate output carries the user's own
    // accessions and no lift-back is required. Opt-in renaming reproduces the
    // historical `<prefix>_N` behaviour (and is what `submit` wants for NCBI's
    // short-name requirement).
    let pad = records.len().to_string().len().max(3);
    let mut rename_map: HashMap<String, String> = HashMap::new();
    for (i, rec) in records.iter().enumerate() {
        let new_id = if config.rename_contigs {
            format!("{}_{:0>width$}", config.prefix, i + 1, width = pad)
        } else {
            rec.original_id.clone()
        };
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
    if config.rename_contigs {
        println!(
            "  Seqids:         renamed to {}_N (--ncbi-clean)",
            config.prefix
        );
    } else {
        println!("  Seqids:         preserved (original accessions kept)");
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn tmp(label: &str, ext: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!(
            "myconote_sort_test_{}_{}.{}",
            std::process::id(),
            label,
            ext
        ));
        p
    }

    /// Two contigs with real-looking accessions; the second is longer so a
    /// by-length sort would reorder them (letting us tell preserve from rename).
    fn write_input(label: &str) -> PathBuf {
        let p = tmp(label, "fa");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, ">NODE_7_length_30 some description").unwrap();
        writeln!(f, "ACGTACGTAC").unwrap();
        writeln!(f, ">CP012345.1 chromosome 1").unwrap();
        writeln!(f, "ACGTACGTACGTACGTACGT").unwrap();
        p
    }

    fn output_headers(path: &Path) -> Vec<String> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|l| l.starts_with('>'))
            .map(|l| l[1..].split_whitespace().next().unwrap().to_string())
            .collect()
    }

    #[test]
    fn default_preserves_original_seqids_through_sort() {
        let input = write_input("preserve");
        let out = tmp("preserve", "out.fa");
        let config = SortConfig {
            input: input.clone(),
            output: out.clone(),
            ..SortConfig::default()
        };
        let map = run_sort(&config).unwrap();

        // Identity map: every original id maps to itself.
        assert_eq!(map.get("CP012345.1"), Some(&"CP012345.1".to_string()));
        assert_eq!(
            map.get("NODE_7_length_30"),
            Some(&"NODE_7_length_30".to_string())
        );

        // Output headers are the original accessions, longest first (CP… is
        // longer than NODE_7…), and contain no scaffold_ renaming.
        let headers = output_headers(&out);
        assert_eq!(headers, vec!["CP012345.1", "NODE_7_length_30"]);
        assert!(!headers.iter().any(|h| h.starts_with("scaffold_")));

        let _ = std::fs::remove_file(&input);
        let _ = std::fs::remove_file(&out);
    }

    #[test]
    fn ncbi_clean_renames_to_scaffold_n_and_emits_table() {
        let input = write_input("rename");
        let out = tmp("rename", "out.fa");
        let table = tmp("rename", "tsv");
        let config = SortConfig {
            input: input.clone(),
            output: out.clone(),
            rename_contigs: true,
            rename_table: Some(table.clone()),
            ..SortConfig::default()
        };
        let map = run_sort(&config).unwrap();

        // Longest contig (CP012345.1) becomes scaffold_001.
        assert_eq!(map.get("CP012345.1"), Some(&"scaffold_001".to_string()));
        assert_eq!(
            map.get("NODE_7_length_30"),
            Some(&"scaffold_002".to_string())
        );

        let headers = output_headers(&out);
        assert_eq!(headers, vec!["scaffold_001", "scaffold_002"]);

        // Rename table maps original → new.
        let tbl = std::fs::read_to_string(&table).unwrap();
        assert!(tbl.contains("CP012345.1\tscaffold_001\t20"));
        assert!(tbl.contains("NODE_7_length_30\tscaffold_002\t10"));

        let _ = std::fs::remove_file(&input);
        let _ = std::fs::remove_file(&out);
        let _ = std::fs::remove_file(&table);
    }

    #[test]
    fn liftover_roundtrips_under_rename_flag() {
        // With --ncbi-clean the rename map lifts a GFF3's seqids over to the
        // new names, and nothing else on the line changes.
        let input = write_input("lift");
        let out = tmp("lift", "out.fa");
        let config = SortConfig {
            input: input.clone(),
            output: out.clone(),
            rename_contigs: true,
            ..SortConfig::default()
        };
        let map = run_sort(&config).unwrap();

        let gff_in = tmp("lift", "in.gff3");
        {
            let mut f = std::fs::File::create(&gff_in).unwrap();
            writeln!(f, "##gff-version 3").unwrap();
            writeln!(f, "CP012345.1\tmyco\tgene\t5\t15\t.\t+\t.\tID=g1").unwrap();
            writeln!(f, "NODE_7_length_30\tmyco\tgene\t1\t8\t.\t-\t.\tID=g2").unwrap();
        }
        let gff_out = tmp("lift", "out.gff3");
        let updated = liftover_gff3(&gff_in, &gff_out, &map).unwrap();
        assert_eq!(updated, 2);

        let lifted = std::fs::read_to_string(&gff_out).unwrap();
        assert!(lifted.contains("scaffold_001\tmyco\tgene\t5\t15\t.\t+\t.\tID=g1"));
        assert!(lifted.contains("scaffold_002\tmyco\tgene\t1\t8\t.\t-\t.\tID=g2"));
        // Original accessions are gone from the lifted GFF.
        assert!(!lifted.contains("CP012345.1\tmyco"));

        for p in [&input, &out, &gff_in, &gff_out] {
            let _ = std::fs::remove_file(p);
        }
    }

    #[test]
    fn preserve_mode_liftover_is_identity() {
        // Under the default (preserve) map, lifting a GFF3 leaves seqids intact.
        let input = write_input("idlift");
        let out = tmp("idlift", "out.fa");
        let config = SortConfig {
            input: input.clone(),
            output: out.clone(),
            ..SortConfig::default()
        };
        let map = run_sort(&config).unwrap();

        let gff_in = tmp("idlift", "in.gff3");
        {
            let mut f = std::fs::File::create(&gff_in).unwrap();
            writeln!(f, "CP012345.1\tmyco\tgene\t5\t15\t.\t+\t.\tID=g1").unwrap();
        }
        let gff_out = tmp("idlift", "out.gff3");
        liftover_gff3(&gff_in, &gff_out, &map).unwrap();
        let lifted = std::fs::read_to_string(&gff_out).unwrap();
        assert!(lifted.contains("CP012345.1\tmyco\tgene\t5\t15"));

        for p in [&input, &out, &gff_in, &gff_out] {
            let _ = std::fs::remove_file(p);
        }
    }
}
