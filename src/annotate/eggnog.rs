/// EggNog-mapper integration
///
/// EggNog-mapper assigns COG/NOG functional categories and deeper
/// functional descriptions by searching against the EggNog 5 database.
///
/// Two modes are supported:
///   1. **Local** — run `emapper.py` if it is in PATH (fastest, needs ~50 GB DB)
///   2. **Pre-computed** — parse an existing `*.emapper.annotations` file
///      (user ran emapper externally, or downloaded results)
///
/// Either way the annotations are merged into the GFF3 attribute field
/// and the TSV annotation table under the keys:
///   eggnog_cog  — COG single-letter category (e.g. "C", "K,T")
///   eggnog_desc — functional description
///   eggnog_og   — best orthologous group (e.g. "COG0001@1|root")
///   eggnog_kegg — KEGG pathway IDs (comma-separated)
///   eggnog_go   — GO terms from EggNog (supplement to UniProt GO)

use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Data type
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct EggNogHit {
    /// Query gene / protein ID
    pub query_id:    String,
    /// Best orthologous group  (field 1 in emapper output)
    pub best_og:     String,
    /// COG functional category letter(s) (field 20)
    pub cog_cat:     String,
    /// Functional description (field 21)
    pub description: String,
    /// Preferred gene name (field 22), e.g. "tubulinA"
    pub gene_name:   String,
    /// GO terms (field 9)
    pub go_terms:    Vec<String>,
    /// KEGG pathway IDs (field 11)
    pub kegg_paths:  Vec<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Parse emapper.annotations file
// ─────────────────────────────────────────────────────────────────────────────

/// Parse an `*.emapper.annotations` file and return a map of
/// query_id → EggNogHit.
pub fn parse_emapper_results(path: &Path) -> Result<HashMap<String, EggNogHit>> {
    let file   = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.is_empty() { continue; }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        // emapper v2 output has at least 22 columns
        if fields.len() < 22 { continue; }

        let query_id    = fields[0].to_string();
        let best_og     = fields[1].to_string();
        // fields[8] = GOs, fields[10] = KEGG pathways
        let go_terms: Vec<String> = if fields[9] == "-" {
            vec![]
        } else {
            fields[9].split(',').map(|s| s.trim().to_string()).collect()
        };
        let kegg_paths: Vec<String> = if fields[11] == "-" {
            vec![]
        } else {
            fields[11].split(',').map(|s| s.trim().to_string()).collect()
        };
        let cog_cat     = fields[20].to_string();
        let description = fields[21].to_string();
        let gene_name   = if fields.len() > 22 { fields[22].to_string() } else { String::new() };

        map.insert(query_id.clone(), EggNogHit {
            query_id,
            best_og,
            cog_cat,
            description,
            gene_name,
            go_terms,
            kegg_paths,
        });
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// Run emapper.py locally
// ─────────────────────────────────────────────────────────────────────────────

/// Run `emapper.py` on a protein FASTA and return the path to the
/// generated `.emapper.annotations` file.
pub fn run_emapper(
    protein_fasta: &Path,
    out_dir:       &Path,
    db_dir:        Option<&Path>,
    threads:       usize,
) -> Result<PathBuf> {
    let emapper = which_emapper()?;

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let out_prefix = out_dir.join("eggnog");

    let mut cmd = Command::new(&emapper);
    cmd.arg("-i").arg(protein_fasta)
       .arg("-o").arg(&out_prefix)
       .arg("--cpu").arg(threads.to_string())
       .arg("--override");

    if let Some(db) = db_dir {
        cmd.arg("--data_dir").arg(db);
    }

    let status = cmd.status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("emapper.py: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "emapper.py exited with non-zero status".to_string()
        ));
    }

    Ok(out_dir.join("eggnog.emapper.annotations"))
}

/// Detect emapper.py in PATH.
pub fn which_emapper() -> Result<String> {
    let candidates = ["emapper.py", "emapper"];
    for c in &candidates {
        if Command::new("which").arg(c).output().map(|o| o.status.success()).unwrap_or(false) {
            return Ok(c.to_string());
        }
    }
    Err(MycoNoteError::ExternalTool(
        "emapper.py not found in PATH. Install with: conda install -c bioconda eggnog-mapper".to_string()
    ))
}

pub fn emapper_available() -> bool {
    which_emapper().is_ok()
}

// ─────────────────────────────────────────────────────────────────────────────
// Write COG summary table
// ─────────────────────────────────────────────────────────────────────────────

/// Write a two-column TSV: gene_id, cog_category, description, gene_name, kegg
pub fn write_eggnog_table(
    hits:   &HashMap<String, EggNogHit>,
    output: &Path,
) -> Result<usize> {
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(out, "gene_id\tbest_og\tcog_category\tdescription\tgene_name\tgo_terms\tkegg_pathways")
        .map_err(MycoNoteError::Io)?;

    let mut sorted: Vec<&EggNogHit> = hits.values().collect();
    sorted.sort_by(|a, b| a.query_id.cmp(&b.query_id));

    for hit in &sorted {
        writeln!(out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            hit.query_id,
            hit.best_og,
            hit.cog_cat,
            hit.description,
            hit.gene_name,
            hit.go_terms.join(","),
            hit.kegg_paths.join(","),
        ).map_err(MycoNoteError::Io)?;
    }

    Ok(sorted.len())
}

// ─────────────────────────────────────────────────────────────────────────────
// COG category summary
// ─────────────────────────────────────────────────────────────────────────────

/// Human-readable COG category descriptions
pub fn cog_category_name(cat: char) -> &'static str {
    match cat {
        'J' => "Translation, ribosomal structure and biogenesis",
        'A' => "RNA processing and modification",
        'K' => "Transcription",
        'L' => "Replication, recombination and repair",
        'B' => "Chromatin structure and dynamics",
        'D' => "Cell cycle control, division, chromosome partitioning",
        'Y' => "Nuclear structure",
        'V' => "Defense mechanisms",
        'T' => "Signal transduction mechanisms",
        'M' => "Cell wall/membrane/envelope biogenesis",
        'N' => "Cell motility",
        'Z' => "Cytoskeleton",
        'W' => "Extracellular structures",
        'U' => "Intracellular trafficking, secretion, vesicular transport",
        'O' => "Posttranslational modification, protein turnover, chaperones",
        'X' => "Mobilome: prophages, transposons",
        'C' => "Energy production and conversion",
        'G' => "Carbohydrate transport and metabolism",
        'E' => "Amino acid transport and metabolism",
        'F' => "Nucleotide transport and metabolism",
        'H' => "Coenzyme transport and metabolism",
        'I' => "Lipid transport and metabolism",
        'P' => "Inorganic ion transport and metabolism",
        'Q' => "Secondary metabolites biosynthesis, transport and catabolism",
        'R' => "General function prediction only",
        'S' => "Function unknown",
        _   => "Unknown category",
    }
}

/// Print a COG category breakdown to stdout.
pub fn print_cog_summary(hits: &HashMap<String, EggNogHit>) {
    let mut cat_counts: HashMap<char, usize> = HashMap::new();
    for hit in hits.values() {
        for c in hit.cog_cat.chars() {
            if c != '-' && c != ' ' {
                *cat_counts.entry(c).or_insert(0) += 1;
            }
        }
    }

    if cat_counts.is_empty() { return; }

    let mut sorted: Vec<(char, usize)> = cat_counts.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));

    println!("  COG category breakdown:");
    for (cat, count) in &sorted {
        println!("    [{}] {:3}  {}", cat, count, cog_category_name(*cat));
    }
}
