/// Pfam domain annotation via hmmscan
///
/// Searches predicted proteins against the Pfam-A HMM database.
/// Results add Pfam domain IDs to gene records.
///
/// Install: conda install -c bioconda hmmer
/// DB setup: myconote annotate --download-dbs   (downloads Pfam-A.hmm)

use crate::utils::error::{MycoNoteError, Result};
use super::AnnotateConfig;
use std::path::Path;
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Hit record
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct PfamHit {
    pub query_id:    String,
    pub domain_id:   String,   // e.g. "PF00001"
    pub domain_name: String,   // e.g. "7tm_1"
    pub evalue:      f64,
    pub bitscore:    f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run(
    query_fa:    &Path,
    pfam_hmm:    &Path,
    output_tsv:  &Path,
    config:      &AnnotateConfig,
) -> Result<Vec<PfamHit>> {
    let hmmscan = which::which("hmmscan").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "hmmscan not found in PATH.\n\
             Install with: conda install -c bioconda hmmer".to_string()
        )
    })?;

    println!("  Running hmmscan (Pfam-A)…");

    // Press the database if .h3i index is missing
    let h3i = pfam_hmm.with_extension("hmm.h3i");
    if !h3i.exists() {
        press_pfam(pfam_hmm)?;
    }

    let domtblout = output_tsv.with_extension("domtblout");

    let status = Command::new(&hmmscan)
        .args([
            "--domtblout", domtblout.to_str().unwrap_or(""),
            "--cpu",       &config.threads.to_string(),
            "-E",          &config.evalue.to_string(),
            "--domE",      &config.evalue.to_string(),
            "--noali",
            pfam_hmm.to_str().unwrap_or(""),
            query_fa.to_str().unwrap_or(""),
        ])
        .stdout(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "hmmscan failed. Check that Pfam-A.hmm is properly pressed.".to_string()
        ));
    }

    let hits = parse_domtblout(&domtblout)?;
    println!("  hmmscan finished — {} domain hits", hits.len());

    // Write a clean TSV summary
    write_pfam_tsv(output_tsv, &hits)?;

    Ok(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// HMM pressing helper
// ─────────────────────────────────────────────────────────────────────────────

fn press_pfam(pfam_hmm: &Path) -> Result<()> {
    let hmmpress = which::which("hmmpress").map_err(|_| {
        MycoNoteError::UnsupportedFormat("hmmpress not found. Install hmmer.".to_string())
    })?;

    println!("  Pressing Pfam-A.hmm (one-time setup)…");
    let status = Command::new(&hmmpress)
        .args(["-f", pfam_hmm.to_str().unwrap_or("")])
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat("hmmpress failed.".to_string()));
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Parser for hmmscan --domtblout format
// ─────────────────────────────────────────────────────────────────────────────

fn parse_domtblout(path: &Path) -> Result<Vec<PfamHit>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut hits: Vec<PfamHit> = Vec::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() { continue; }

        // domtblout columns (space-separated):
        // 0: target(domain) name, 1: accession, 2: target length,
        // 3: query name, 4: accession, 5: query length,
        // 6: full-seq E-value, 7: full-seq score, 8: full-seq bias,
        // 9: domain #, 10: total domains,
        // 11: domain c-evalue, 12: domain i-evalue, 13: domain score, ...
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 14 { continue; }

        let domain_name = f[0].to_string();
        let domain_id   = f[1].trim_start_matches("PF").to_string();
        let domain_id   = if domain_id.is_empty() { domain_name.clone() }
                          else { format!("PF{}", domain_id) };
        let query_id    = f[3].to_string();
        let evalue:    f64 = f[12].parse().unwrap_or(f64::MAX);
        let bitscore:  f64 = f[13].parse().unwrap_or(0.0);

        hits.push(PfamHit { query_id, domain_id, domain_name, evalue, bitscore });
    }

    // Deduplicate: keep best hit per (query, domain) pair
    hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap());
    let mut seen = std::collections::HashSet::new();
    hits.retain(|h| seen.insert((h.query_id.clone(), h.domain_id.clone())));

    Ok(hits)
}

fn write_pfam_tsv(path: &Path, hits: &[PfamHit]) -> Result<()> {
    use std::io::Write;

    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "query_id\tdomain_id\tdomain_name\tevalue\tbitscore")
        .map_err(MycoNoteError::Io)?;
    for h in hits {
        writeln!(f, "{}\t{}\t{}\t{:.2e}\t{:.1}",
            h.query_id, h.domain_id, h.domain_name, h.evalue, h.bitscore)
            .map_err(MycoNoteError::Io)?;
    }
    Ok(())
}
