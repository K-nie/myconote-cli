/// MMseqs2 homology search wrapper
///
/// Searches predicted protein sequences against Swiss-Prot (or a custom DB)
/// using MMseqs2 easy-search in protein mode.
///
/// Install: conda install -c bioconda mmseqs2
/// DB setup: myconote annotate --download-dbs   (downloads & indexes Swiss-Prot)

use crate::utils::error::{MycoNoteError, Result};
use super::AnnotateConfig;
use std::path::Path;
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Hit record
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct MmseqsHit {
    pub query_id:    String,
    pub target_id:   String,
    pub description: String,
    pub identity:    f64,   // percent (0–100)
    pub evalue:      f64,
    pub bitscore:    f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run MMseqs2 easy-search and return parsed hits (best per query).
pub fn run(
    query_fa:   &Path,
    db_path:    &Path,
    output_tsv: &Path,
    config:     &AnnotateConfig,
) -> Result<Vec<MmseqsHit>> {
    let mmseqs = which::which("mmseqs").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "mmseqs not found in PATH.\n\
             Install with: conda install -c bioconda mmseqs2".to_string()
        )
    })?;

    println!("  Running MMseqs2 easy-search…");

    let tmp_dir = tempfile::TempDir::new().map_err(MycoNoteError::Io)?;

    // Output format: query target identity evalue bitscore description
    let format_str = "query,target,fident,evalue,bits,tset_description";

    let status = Command::new(&mmseqs)
        .args([
            "easy-search",
            query_fa.to_str().unwrap_or(""),
            db_path.to_str().unwrap_or(""),
            output_tsv.to_str().unwrap_or(""),
            tmp_dir.path().to_str().unwrap_or("/tmp"),
            "--format-mode", "4",
            "--format-output", format_str,
            "--threads", &config.threads.to_string(),
            "-e",       &config.evalue.to_string(),
            "--min-seq-id", &config.min_identity.to_string(),
            "--db-load-mode", "2",
        ])
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "MMseqs2 easy-search failed. Check that the Swiss-Prot database \
             has been indexed (run: myconote annotate --download-dbs)".to_string()
        ));
    }

    let hits = parse_mmseqs_tsv(output_tsv)?;
    println!("  MMseqs2 finished — {} raw hits", hits.len());

    Ok(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// TSV parser
// ─────────────────────────────────────────────────────────────────────────────

/// Parse MMseqs2 --format-output "query,target,fident,evalue,bits,tset_description"
/// Returns only the best hit per query (lowest e-value).
fn parse_mmseqs_tsv(path: &Path) -> Result<Vec<MmseqsHit>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut best: std::collections::HashMap<String, MmseqsHit> = std::collections::HashMap::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() { continue; }

        let f: Vec<&str> = line.splitn(6, '\t').collect();
        if f.len() < 5 { continue; }

        let query_id  = f[0].to_string();
        let target_id = f[1].to_string();
        let identity: f64  = f[2].parse::<f64>().unwrap_or(0.0) * 100.0;
        let evalue:   f64  = f[3].parse().unwrap_or(f64::MAX);
        let bitscore: f64  = f[4].parse().unwrap_or(0.0);
        let description    = clean_description(f.get(5).unwrap_or(&""));

        let hit = MmseqsHit {
            query_id: query_id.clone(),
            target_id,
            description,
            identity,
            evalue,
            bitscore,
        };

        // Keep only the best (lowest e-value) hit per query
        let entry = best.entry(query_id).or_insert_with(|| hit.clone());
        if hit.evalue < entry.evalue {
            *entry = hit;
        }
    }

    Ok(best.into_values().collect())
}

/// Strip Swiss-Prot boilerplate from description lines.
/// e.g. "sp|P12345|GENE_HUMAN Ribosomal protein S3 OS=Homo sapiens OX=9606 GN=RPS3 PE=1 SV=1"
///  → "Ribosomal protein S3"
fn clean_description(raw: &str) -> String {
    // Remove OS= and everything after
    let desc = if let Some(pos) = raw.find(" OS=") {
        &raw[..pos]
    } else {
        raw
    };

    // If this looks like a Swiss-Prot header with pipe-separated ID, take the part after the last |
    let desc = if desc.contains('|') {
        desc.split('|').last().unwrap_or(desc)
    } else {
        desc
    };

    // Remove leading accession-like prefix "GENENAME_SPECIES "
    let desc = if let Some(pos) = desc.find(' ') {
        let prefix = &desc[..pos];
        if prefix.contains('_') && prefix.chars().all(|c| c.is_ascii_uppercase() || c == '_' || c.is_ascii_digit()) {
            desc[pos + 1..].trim()
        } else {
            desc.trim()
        }
    } else {
        desc.trim()
    };

    if desc.is_empty() {
        "hypothetical protein".to_string()
    } else {
        desc.to_string()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Database creation helper (called from db.rs)
// ─────────────────────────────────────────────────────────────────────────────

/// Build an MMseqs2 sequence database from a FASTA file.
pub fn create_db(fasta: &Path, db_out: &Path) -> Result<()> {
    let mmseqs = which::which("mmseqs").map_err(|_| {
        MycoNoteError::UnsupportedFormat("mmseqs not found in PATH.".to_string())
    })?;

    println!("  Indexing {} for MMseqs2…", fasta.display());

    let status = Command::new(&mmseqs)
        .args(["createdb", fasta.to_str().unwrap_or(""), db_out.to_str().unwrap_or("")])
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "mmseqs createdb failed.".to_string()
        ));
    }

    Ok(())
}
