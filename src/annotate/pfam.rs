use super::AnnotateConfig;
/// Pfam domain annotation via HMMER
///
/// Searches predicted proteins against the Pfam-A HMM database using either:
///   - hmmsearch (preferred): searches each HMM profile against all proteins
///     at once. Much faster for large protein sets (5,000+) because HMMER's
///     internal SSV/MSV acceleration works on the full sequence database.
///   - hmmscan (fallback): searches each protein against all HMM profiles.
///     Slower but available when hmmsearch can't be used.
///
/// Performance note (HMMER docs):
///   "hmmsearch is generally faster than hmmscan for searching a large
///    sequence database with a profile database, because the overhead of
///    loading each profile is amortized across all sequences."
///
/// For a typical fungal genome (5,000 proteins × 20,795 Pfam-A profiles):
///   hmmscan:   ~45 min (searches 5K seqs × 20K profiles one-by-one)
///   hmmsearch: ~8 min  (loads 20K profiles × searches all 5K seqs at once)
///
/// Install: conda install -c bioconda hmmer
/// DB setup: myconote annotate --download-dbs   (downloads Pfam-A.hmm)
use crate::utils::error::{MycoNoteError, Result};
use std::io::Write;
use std::path::Path;
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Hit record
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct PfamHit {
    pub query_id: String,
    pub domain_id: String,   // e.g. "PF00001.21"
    pub domain_name: String, // e.g. "7tm_1"
    pub evalue: f64,
    pub bitscore: f64,
    pub ali_from: u64,
    pub ali_to: u64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run(
    query_fa: &Path,
    pfam_hmm: &Path,
    output_tsv: &Path,
    config: &AnnotateConfig,
) -> Result<Vec<PfamHit>> {
    // Press the database if .h3i index is missing
    let h3i = pfam_hmm.with_extension("hmm.h3i");
    if !h3i.exists() {
        press_pfam(pfam_hmm)?;
    }

    // Try hmmsearch first (much faster), fall back to hmmscan
    let hmmsearch_available = which::which("hmmsearch").is_ok();
    let hmmscan_available = which::which("hmmscan").is_ok();

    if !hmmsearch_available && !hmmscan_available {
        return Err(MycoNoteError::UnsupportedFormat(
            "Neither hmmsearch nor hmmscan found in PATH.\n\
             Install with: conda install -c bioconda hmmer"
                .to_string(),
        ));
    }

    let domtblout = output_tsv.with_extension("domtblout");

    let hits = if hmmsearch_available {
        run_hmmsearch(query_fa, pfam_hmm, &domtblout, config)?
    } else {
        run_hmmscan(query_fa, pfam_hmm, &domtblout, config)?
    };

    // Write a clean TSV summary
    write_pfam_tsv(output_tsv, &hits)?;

    Ok(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// hmmsearch (preferred — faster for large protein sets)
// ─────────────────────────────────────────────────────────────────────────────

fn run_hmmsearch(
    query_fa: &Path,
    pfam_hmm: &Path,
    domtblout: &Path,
    config: &AnnotateConfig,
) -> Result<Vec<PfamHit>> {
    let hmmsearch = which::which("hmmsearch")
        .map_err(|_| MycoNoteError::UnsupportedFormat("hmmsearch not found".to_string()))?;

    println!("  Running hmmsearch (Pfam-A) — faster than hmmscan for large proteomes");

    let pfam_s = pfam_hmm.to_string_lossy();
    let query_s = query_fa.to_string_lossy();
    let dom_s = domtblout.to_string_lossy();

    let status = Command::new(&hmmsearch)
        .args([
            "--domtblout",
            dom_s.as_ref(),
            "--cpu",
            &config.threads.to_string(),
            "-E",
            &config.evalue.to_string(),
            "--domE",
            &config.evalue.to_string(),
            "--noali", // skip alignment output (saves time + disk)
            pfam_s.as_ref(),
            query_s.as_ref(),
        ])
        .stdout(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "hmmsearch failed. Check that Pfam-A.hmm is properly pressed (hmmpress).".to_string(),
        ));
    }

    // hmmsearch --domtblout has the same format as hmmscan --domtblout
    // but columns 0-2 and 3-5 are swapped:
    //   hmmscan:   target=domain, query=protein
    //   hmmsearch: target=protein, query=domain
    let hits = parse_domtblout_hmmsearch(domtblout)?;
    println!("  hmmsearch finished — {} domain hits", hits.len());

    Ok(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// hmmscan (fallback — slower but works everywhere)
// ─────────────────────────────────────────────────────────────────────────────

fn run_hmmscan(
    query_fa: &Path,
    pfam_hmm: &Path,
    domtblout: &Path,
    config: &AnnotateConfig,
) -> Result<Vec<PfamHit>> {
    let hmmscan = which::which("hmmscan")
        .map_err(|_| MycoNoteError::UnsupportedFormat("hmmscan not found".to_string()))?;

    println!("  Running hmmscan (Pfam-A) — fallback mode");

    let pfam_s = pfam_hmm.to_string_lossy();
    let query_s = query_fa.to_string_lossy();
    let dom_s = domtblout.to_string_lossy();

    let status = Command::new(&hmmscan)
        .args([
            "--domtblout",
            dom_s.as_ref(),
            "--cpu",
            &config.threads.to_string(),
            "-E",
            &config.evalue.to_string(),
            "--domE",
            &config.evalue.to_string(),
            "--noali",
            pfam_s.as_ref(),
            query_s.as_ref(),
        ])
        .stdout(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "hmmscan failed. Check that Pfam-A.hmm is properly pressed (hmmpress).".to_string(),
        ));
    }

    let hits = parse_domtblout_hmmscan(domtblout)?;
    println!("  hmmscan finished — {} domain hits", hits.len());

    Ok(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// HMM pressing helper
// ─────────────────────────────────────────────────────────────────────────────

fn press_pfam(pfam_hmm: &Path) -> Result<()> {
    let hmmpress = which::which("hmmpress").map_err(|_| {
        MycoNoteError::UnsupportedFormat("hmmpress not found. Install hmmer.".to_string())
    })?;

    let pfam_s = pfam_hmm.to_string_lossy();
    println!("  Pressing Pfam-A.hmm (one-time setup)...");
    let status = Command::new(&hmmpress)
        .args(["-f", pfam_s.as_ref()])
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat("hmmpress failed.".to_string()));
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Parser for hmmsearch --domtblout format
// ─────────────────────────────────────────────────────────────────────────────
// hmmsearch domtblout columns (space-delimited):
//   0: target name (PROTEIN)    1: target accession  2: target length
//   3: query name (DOMAIN)      4: query accession   5: query length
//   6: full E-value   7: full score   8: full bias
//   9: dom#  10: ndom
//  11: dom c-Evalue  12: dom i-Evalue  13: dom score  14: dom bias
//  15: hmm from  16: hmm to  17: ali from  18: ali to  19: env from  20: env to
//  21: acc  22+: target description

fn parse_domtblout_hmmsearch(path: &Path) -> Result<Vec<PfamHit>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut hits: Vec<PfamHit> = Vec::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 22 {
            continue;
        }

        // In hmmsearch: col 0 = protein (target), col 3 = domain (query)
        let query_id = f[0].to_string(); // protein name
        let domain_name = f[3].to_string(); // domain name
        let domain_acc = f[4].to_string(); // domain accession (PFxxxxx.xx)
        let evalue: f64 = f[12].parse().unwrap_or(f64::MAX); // domain i-Evalue
        let bitscore: f64 = f[13].parse().unwrap_or(0.0); // domain score
        let ali_from: u64 = f[17].parse().unwrap_or(0);
        let ali_to: u64 = f[18].parse().unwrap_or(0);

        // Normalize domain ID to PFxxxxx.xx format
        let domain_id = normalize_pfam_acc(&domain_acc, &domain_name);

        hits.push(PfamHit {
            query_id,
            domain_id,
            domain_name,
            evalue,
            bitscore,
            ali_from,
            ali_to,
        });
    }

    dedup_hits(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// Parser for hmmscan --domtblout format
// ─────────────────────────────────────────────────────────────────────────────
// hmmscan domtblout columns (space-delimited):
//   0: target name (DOMAIN)     1: target accession  2: target length
//   3: query name (PROTEIN)     4: query accession   5: query length
//   ...same structure after that

fn parse_domtblout_hmmscan(path: &Path) -> Result<Vec<PfamHit>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut hits: Vec<PfamHit> = Vec::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 22 {
            continue;
        }

        // In hmmscan: col 0 = domain (target), col 3 = protein (query)
        let domain_name = f[0].to_string();
        let domain_acc = f[1].to_string();
        let query_id = f[3].to_string();
        let evalue: f64 = f[12].parse().unwrap_or(f64::MAX);
        let bitscore: f64 = f[13].parse().unwrap_or(0.0);
        let ali_from: u64 = f[17].parse().unwrap_or(0);
        let ali_to: u64 = f[18].parse().unwrap_or(0);

        let domain_id = normalize_pfam_acc(&domain_acc, &domain_name);

        hits.push(PfamHit {
            query_id,
            domain_id,
            domain_name,
            evalue,
            bitscore,
            ali_from,
            ali_to,
        });
    }

    dedup_hits(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Normalize a Pfam accession to "PFxxxxx.xx" format.
fn normalize_pfam_acc(acc: &str, name: &str) -> String {
    if acc.starts_with("PF") {
        acc.to_string()
    } else if name.starts_with("PF") {
        name.to_string()
    } else {
        // Use the accession as-is (might be a Pfam-B or clan ID)
        acc.to_string()
    }
}

/// Deduplicate: keep best hit per (query, domain) pair by lowest e-value.
fn dedup_hits(mut hits: Vec<PfamHit>) -> Result<Vec<PfamHit>> {
    hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut seen = std::collections::HashSet::new();
    hits.retain(|h| seen.insert((h.query_id.clone(), h.domain_id.clone())));
    Ok(hits)
}

fn write_pfam_tsv(path: &Path, hits: &[PfamHit]) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "query_id\tdomain_id\tdomain_name\tevalue\tbitscore\tali_from\tali_to"
    )
    .map_err(MycoNoteError::Io)?;
    for h in hits {
        writeln!(
            f,
            "{}\t{}\t{}\t{:.2e}\t{:.1}\t{}\t{}",
            h.query_id, h.domain_id, h.domain_name, h.evalue, h.bitscore, h.ali_from, h.ali_to
        )
        .map_err(MycoNoteError::Io)?;
    }
    Ok(())
}
