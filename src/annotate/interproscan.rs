/// InterProScan 5 REST API integration
///
/// Searches predicted proteins against 13+ signature databases simultaneously:
/// Pfam, TIGRFAM, SUPERFAMILY, Gene3D, PRINTS, ProSitePatterns, ProSiteProfiles,
/// SMART, CDD, HAMAP, PIRSF, PANTHER, and MobiDBLite.
///
/// Uses the EBI InterProScan5 public REST API — no local installation required.
/// Sequences are submitted in batches; results are polled until complete.
///
/// API base: https://www.ebi.ac.uk/Tools/services/rest/iprscan5
/// Rate limit: 30 jobs / minute; max 30 sequences / job; 30K AA / job.
///
/// Results are cached in ~/.myconote/iprscan_cache.tsv so re-runs are instant.

use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Duration;

// ─────────────────────────────────────────────────────────────────────────────
// API constants
// ─────────────────────────────────────────────────────────────────────────────

const API_BASE:     &str = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5";
const BATCH_SIZE:   usize = 30;    // max sequences per job (API limit)
const POLL_DELAY:   u64   = 15;    // seconds between status polls
const MAX_POLLS:    usize = 120;   // 30 min max wait per batch
const SUBMIT_DELAY: u64   = 3;     // seconds between job submissions (rate limit)

// ─────────────────────────────────────────────────────────────────────────────
// Result record
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct IprScanResult {
    pub protein_id:    String,
    /// InterPro accession (IPR...) — may be empty if hit has no InterPro entry
    pub ipr_accession: Option<String>,
    /// InterPro description
    pub ipr_desc:      Option<String>,
    /// Source database (e.g. "Pfam", "TIGRFAM", "Gene3D")
    pub database:      String,
    /// Database accession (e.g. "PF00001", "TIGR00001")
    pub db_accession:  String,
    /// Domain description from source database
    pub db_desc:       String,
    /// GO terms associated with this InterPro entry
    pub go_terms:      Vec<String>,
    /// Pathway terms (KEGG, MetaCyc, Reactome)
    pub pathways:      Vec<String>,
    /// E-value from the source database (not all databases report this)
    pub evalue:        Option<f64>,
    /// Match coordinates (1-based, in the protein)
    pub start:         u32,
    pub end:           u32,
}

// ─────────────────────────────────────────────────────────────────────────────
// Cache management
// ─────────────────────────────────────────────────────────────────────────────

fn cache_path() -> PathBuf {
    let home = std::env::var("HOME").unwrap_or_else(|_| "/tmp".to_string());
    PathBuf::from(home).join(".myconote").join("iprscan_cache.tsv")
}

/// Load cached results: protein_id → Vec<IprScanResult>
fn load_cache() -> HashMap<String, Vec<IprScanResult>> {
    let mut map: HashMap<String, Vec<IprScanResult>> = HashMap::new();
    let path = cache_path();
    if let Ok(content) = std::fs::read_to_string(&path) {
        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() { continue; }
            if let Some(r) = parse_cache_line(line) {
                map.entry(r.protein_id.clone()).or_default().push(r);
            }
        }
    }
    map
}

fn parse_cache_line(line: &str) -> Option<IprScanResult> {
    let f: Vec<&str> = line.splitn(12, '\t').collect();
    if f.len() < 10 { return None; }
    Some(IprScanResult {
        protein_id:    f[0].to_string(),
        ipr_accession: if f[1].is_empty() { None } else { Some(f[1].to_string()) },
        ipr_desc:      if f[2].is_empty() { None } else { Some(f[2].to_string()) },
        database:      f[3].to_string(),
        db_accession:  f[4].to_string(),
        db_desc:       f[5].to_string(),
        go_terms:      f[6].split('|').filter(|s| !s.is_empty()).map(String::from).collect(),
        pathways:      f[7].split('|').filter(|s| !s.is_empty()).map(String::from).collect(),
        evalue:        f[8].parse().ok(),
        start:         f[9].parse().unwrap_or(0),
        end:           f.get(10).and_then(|s| s.parse().ok()).unwrap_or(0),
    })
}

fn save_cache(results: &[IprScanResult]) {
    let path = cache_path();
    if let Some(parent) = path.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    // Append to cache
    if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&path) {
        for r in results {
            let _ = writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.protein_id,
                r.ipr_accession.as_deref().unwrap_or(""),
                r.ipr_desc.as_deref().unwrap_or(""),
                r.database,
                r.db_accession,
                r.db_desc,
                r.go_terms.join("|"),
                r.pathways.join("|"),
                r.evalue.map(|v| format!("{:.2e}", v)).unwrap_or_default(),
                r.start,
                r.end,
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run InterProScan for all proteins in `proteins_fa`.
/// Returns a map of protein_id → Vec<IprScanResult>.
pub fn run(
    proteins_fa: &Path,
    output_tsv:  &Path,
    email:       &str,
) -> Result<HashMap<String, Vec<IprScanResult>>> {
    // Load protein sequences
    let proteins = load_protein_fasta(proteins_fa)?;
    println!("  InterProScan: {} proteins to search", proteins.len());

    // Load cache
    let mut cache = load_cache();
    let mut all_results: HashMap<String, Vec<IprScanResult>> = HashMap::new();
    let mut missing: Vec<(String, String)> = Vec::new();  // (id, sequence)

    for (id, seq) in &proteins {
        if let Some(cached) = cache.get(id) {
            all_results.insert(id.clone(), cached.clone());
        } else {
            missing.push((id.clone(), seq.clone()));
        }
    }

    println!("  {} cached, {} to submit to EBI API", all_results.len(), missing.len());

    // Submit in batches
    if !missing.is_empty() {
        let total_batches = (missing.len() + BATCH_SIZE - 1) / BATCH_SIZE;
        for (batch_idx, chunk) in missing.chunks(BATCH_SIZE).enumerate() {
            println!("  Batch {}/{}: submitting {} sequences…",
                batch_idx + 1, total_batches, chunk.len());

            match submit_and_poll(chunk, email) {
                Ok(batch_results) => {
                    // Add to cache and results
                    save_cache(&batch_results);
                    for r in batch_results {
                        all_results.entry(r.protein_id.clone()).or_default().push(r.clone());
                        cache.entry(r.protein_id.clone()).or_default().push(r);
                    }
                }
                Err(e) => {
                    eprintln!("  ⚠  InterProScan batch {} failed: {}", batch_idx + 1, e);
                }
            }

            // Rate limit: 30 jobs/min = 1 job every 2 seconds
            if batch_idx + 1 < total_batches {
                std::thread::sleep(Duration::from_secs(SUBMIT_DELAY));
            }
        }
    }

    // Write TSV output
    write_iprscan_tsv(output_tsv, &all_results)?;

    println!("  InterProScan: {} proteins with at least one hit",
        all_results.values().filter(|v| !v.is_empty()).count());

    Ok(all_results)
}

// ─────────────────────────────────────────────────────────────────────────────
// API: submit job
// ─────────────────────────────────────────────────────────────────────────────

fn submit_and_poll(
    sequences: &[(String, String)],
    email:     &str,
) -> Result<Vec<IprScanResult>> {
    let job_id = submit_job(sequences, email)?;
    poll_job(&job_id)?;
    let tsv = fetch_result_tsv(&job_id)?;
    let results = parse_iprscan_tsv(&tsv);
    Ok(results)
}

fn submit_job(sequences: &[(String, String)], email: &str) -> Result<String> {
    // Build FASTA string for submission
    let mut fasta_str = String::new();
    for (id, seq) in sequences {
        fasta_str.push('>');
        fasta_str.push_str(id);
        fasta_str.push('\n');
        fasta_str.push_str(seq);
        fasta_str.push('\n');
    }

    let url = format!("{}/run", API_BASE);

    // Build form data
    let form_body = format!(
        "email={}&title=myconote&goterms=true&pathways=true&sequence={}",
        url_encode(email),
        url_encode(&fasta_str)
    );

    let resp = ureq::post(&url)
        .set("Content-Type", "application/x-www-form-urlencoded")
        .send_string(&form_body)
        .map_err(|e| MycoNoteError::InvalidFormat(format!("InterProScan submit failed: {}", e)))?;

    let job_id = resp.into_string()
        .map_err(|e| MycoNoteError::InvalidFormat(format!("InterProScan response read error: {}", e)))?
        .trim()
        .to_string();

    if job_id.is_empty() || job_id.contains("error") {
        return Err(MycoNoteError::InvalidFormat(
            format!("InterProScan job submission returned invalid job ID: {}", job_id)
        ));
    }

    Ok(job_id)
}

fn poll_job(job_id: &str) -> Result<()> {
    let status_url = format!("{}/status/{}", API_BASE, job_id);

    for poll in 0..MAX_POLLS {
        std::thread::sleep(Duration::from_secs(POLL_DELAY));

        let status = match ureq::get(&status_url).call() {
            Ok(r)  => r.into_string().unwrap_or_default().trim().to_string(),
            Err(e) => {
                eprintln!("  ⚠  Poll {}: status check failed: {}", poll + 1, e);
                continue;
            }
        };

        match status.as_str() {
            "FINISHED"  => return Ok(()),
            "RUNNING" | "QUEUED" => {
                if poll % 4 == 0 {
                    print!("    [{} min] status: {}…\r",
                        (poll as u64 * POLL_DELAY) / 60, status);
                    let _ = std::io::stdout().flush();
                }
            }
            "FAILED" => return Err(MycoNoteError::InvalidFormat(
                format!("InterProScan job {} failed on EBI servers.", job_id)
            )),
            other => {
                eprintln!("  ⚠  Unexpected job status: {}", other);
            }
        }
    }

    Err(MycoNoteError::InvalidFormat(format!(
        "InterProScan job {} timed out after {} polls.", job_id, MAX_POLLS
    )))
}

fn fetch_result_tsv(job_id: &str) -> Result<String> {
    let url = format!("{}/result/{}/tsv", API_BASE, job_id);
    ureq::get(&url)
        .call()
        .map_err(|e| MycoNoteError::InvalidFormat(format!("InterProScan result fetch failed: {}", e)))?
        .into_string()
        .map_err(|e| MycoNoteError::InvalidFormat(format!("InterProScan result read error: {}", e)))
}

// ─────────────────────────────────────────────────────────────────────────────
// TSV result parser
// ─────────────────────────────────────────────────────────────────────────────

/// Parse InterProScan TSV output format.
///
/// Columns (15 total):
/// 0  Protein accession
/// 1  Sequence MD5 digest
/// 2  Sequence length
/// 3  Analysis (database)
/// 4  Signature accession
/// 5  Signature description
/// 6  Start location
/// 7  Stop location
/// 8  Score
/// 9  Status (T=true)
/// 10 Date
/// 11 InterPro accession
/// 12 InterPro description
/// 13 GO annotations (pipe-separated)
/// 14 Pathways (pipe-separated)
fn parse_iprscan_tsv(tsv: &str) -> Vec<IprScanResult> {
    let mut results = Vec::new();

    for line in tsv.lines() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 { continue; }

        // Skip hits that are not "True" matches
        if f.get(9) == Some(&"FALSE") { continue; }

        let protein_id   = f[0].to_string();
        let database     = f[3].to_string();
        let db_accession = f[4].to_string();
        let db_desc      = f[5].to_string();
        let start:   u32 = f[6].parse().unwrap_or(0);
        let end:     u32 = f[7].parse().unwrap_or(0);
        let evalue       = f[8].parse::<f64>().ok();
        let ipr_acc = f.get(11).and_then(|s| {
            if s.is_empty() || *s == "-" { None } else { Some(s.to_string()) }
        });
        let ipr_desc = f.get(12).and_then(|s| {
            if s.is_empty() || *s == "-" { None } else { Some(s.to_string()) }
        });
        let go_terms: Vec<String> = f.get(13)
            .map(|s| s.split('|').filter(|t| t.starts_with("GO:")).map(String::from).collect())
            .unwrap_or_default();
        let pathways: Vec<String> = f.get(14)
            .map(|s| s.split('|').filter(|t| !t.trim().is_empty()).map(String::from).collect())
            .unwrap_or_default();

        results.push(IprScanResult {
            protein_id,
            ipr_accession: ipr_acc,
            ipr_desc,
            database,
            db_accession,
            db_desc,
            go_terms,
            pathways,
            evalue,
            start,
            end,
        });
    }

    results
}

// ─────────────────────────────────────────────────────────────────────────────
// Output
// ─────────────────────────────────────────────────────────────────────────────

fn write_iprscan_tsv(
    path:    &Path,
    results: &HashMap<String, Vec<IprScanResult>>,
) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "protein_id\tipr_accession\tipr_desc\tdatabase\tdb_accession\tdb_desc\tgo_terms\tevalue\tstart\tend")
        .map_err(MycoNoteError::Io)?;

    let mut protein_ids: Vec<&String> = results.keys().collect();
    protein_ids.sort();

    for pid in protein_ids {
        for r in &results[pid] {
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.protein_id,
                r.ipr_accession.as_deref().unwrap_or(""),
                r.ipr_desc.as_deref().unwrap_or(""),
                r.database,
                r.db_accession,
                r.db_desc,
                r.go_terms.join("|"),
                r.evalue.map(|v| format!("{:.2e}", v)).unwrap_or_default(),
                r.start,
                r.end,
            ).map_err(MycoNoteError::Io)?;
        }
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Protein FASTA loader
// ─────────────────────────────────────────────────────────────────────────────

fn load_protein_fasta(path: &Path) -> Result<Vec<(String, String)>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut proteins: Vec<(String, String)> = Vec::new();
    let mut current_id  = String::new();
    let mut current_seq = String::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim();
        if line.starts_with('>') {
            if !current_id.is_empty() && !current_seq.is_empty() {
                proteins.push((current_id.clone(), current_seq.clone()));
            }
            current_id  = line[1..].split_whitespace().next().unwrap_or("").to_string();
            current_seq = String::new();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_id.is_empty() && !current_seq.is_empty() {
        proteins.push((current_id, current_seq));
    }

    Ok(proteins)
}

// ─────────────────────────────────────────────────────────────────────────────
// URL encoding helper
// ─────────────────────────────────────────────────────────────────────────────

fn url_encode(s: &str) -> String {
    let mut encoded = String::with_capacity(s.len() * 2);
    for byte in s.bytes() {
        match byte {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9'
            | b'-' | b'_' | b'.' | b'~' => {
                encoded.push(byte as char);
            }
            b' ' => encoded.push('+'),
            other => {
                encoded.push('%');
                encoded.push_str(&format!("{:02X}", other));
            }
        }
    }
    encoded
}
