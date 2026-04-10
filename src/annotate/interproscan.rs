/// InterProScan 5 REST API integration
///
/// API base: https://www.ebi.ac.uk/Tools/services/rest/iprscan5
/// Rate limit: 30 jobs / minute; max 30 sequences / job; 30K AA / job.
/// Results cached in ~/.myconote/iprscan_cache.tsv.
///
/// Pipeline:
///   Step 0 — Pre-validate & auto-fix every sequence before any submission.
///   Phase 1 — Rate-limited sequential submission; HTTP 400 batches are
///              binary-searched to isolate any remaining bad sequences.
///   Phase 2 — All submitted jobs polled concurrently; persistent-ERROR
///              batches are binary-searched the same way.
///   Report  — Final table of every rejected protein and the exact reason.
use crate::utils::error::{MycoNoteError, Result};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc,
};
use std::time::Duration;

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

const API_BASE: &str = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5";
const BATCH_SIZE: usize = 30;
const POLL_DELAY: u64 = 15;
const MAX_POLLS: usize = 120;
const SUBMIT_DELAY: u64 = 3;
const MAX_RETRIES: usize = 3;
const MAX_CONSEC_ERRORS: usize = 8;
const DNS_RETRY_DELAY: u64 = 20;
const MAX_DNS_RETRIES: usize = 5;
const MAX_SEQ_LEN: usize = 40_000;
const MIN_SEQ_LEN: usize = 10;
const MAX_X_FRACTION: f64 = 0.80; // reject if >80% unknown residues
const NO_HITS_TAG: &str = "__NO_HITS__";

// ─────────────────────────────────────────────────────────────────────────────
// Public result type
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct IprScanResult {
    pub protein_id: String,
    pub ipr_accession: Option<String>,
    pub ipr_desc: Option<String>,
    pub database: String,
    pub db_accession: String,
    pub db_desc: String,
    pub go_terms: Vec<String>,
    pub pathways: Vec<String>,
    pub evalue: Option<f64>,
    pub start: u32,
    pub end: u32,
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal types
// ─────────────────────────────────────────────────────────────────────────────

/// Outcome of pre-processing one protein sequence.
enum PrepOutcome {
    /// Sequence is ready to submit (possibly with minor fixes applied).
    Ok {
        cleaned: String,
        /// Non-empty if any auto-fix was applied (for user reporting).
        fix_note: Option<String>,
    },
    /// Sequence is unfixable and will be skipped.
    Skip { reason: String },
}

/// A protein that failed and was not submitted to EBI.
#[derive(Debug, Clone)]
struct BadProtein {
    protein_id: String,
    reason: String,
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 0 — Sequence pre-validation and auto-fixing
// ─────────────────────────────────────────────────────────────────────────────

/// Characters EBI InterProScan accepts without complaint.
/// Standard 20 AA + common IUPAC ambiguity codes.
const VALID_AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwyBJOUXZbjoxuz";

/// Try to make a sequence submittable; returns PrepOutcome.
///
/// Auto-fixes applied (in order):
///   1. Strip trailing stop codon '*'.
///   2. Truncate at the first *internal* stop codon (with a note).
///   3. Replace non-standard characters with 'X' (unknown residue).
///   4. Upper-case the sequence.
///
/// Rejection criteria (after fixes):
///   • Length < MIN_SEQ_LEN or > MAX_SEQ_LEN.
///   • Fraction of X residues exceeds MAX_X_FRACTION.
///   • Looks like a nucleotide sequence (>80% ATCG).
fn preprocess_sequence(protein_id: &str, seq: &str) -> PrepOutcome {
    let mut notes: Vec<String> = Vec::new();

    // 1. Strip trailing '*'
    let after_trailing = seq.trim_end_matches('*');

    // 2. Truncate at internal '*'
    let truncated = if let Some(pos) = after_trailing.find('*') {
        notes.push(format!(
            "truncated at internal stop codon at position {} (was {} AA)",
            pos + 1,
            after_trailing.len()
        ));
        &after_trailing[..pos]
    } else {
        after_trailing
    };

    if truncated.len() < MIN_SEQ_LEN {
        return PrepOutcome::Skip {
            reason: format!(
                "sequence too short after stop-codon removal ({} AA; minimum {})",
                truncated.len(),
                MIN_SEQ_LEN
            ),
        };
    }
    if truncated.len() > MAX_SEQ_LEN {
        return PrepOutcome::Skip {
            reason: format!(
                "sequence too long ({} AA; EBI limit {})",
                truncated.len(),
                MAX_SEQ_LEN
            ),
        };
    }

    // 3. Replace non-standard characters; upper-case everything
    let mut cleaned = String::with_capacity(truncated.len());
    let mut replaced_count = 0usize;
    let mut replaced_chars: std::collections::HashSet<char> = std::collections::HashSet::new();

    for ch in truncated.chars() {
        let upper = ch.to_ascii_uppercase();
        if VALID_AA.contains(&(upper as u8)) {
            cleaned.push(upper);
        } else {
            cleaned.push('X');
            replaced_count += 1;
            replaced_chars.insert(ch);
        }
    }

    if replaced_count > 0 {
        let mut chars: Vec<char> = replaced_chars.into_iter().collect();
        chars.sort_unstable();
        notes.push(format!(
            "replaced {} non-standard character(s) {:?} with 'X'",
            replaced_count, chars
        ));
    }

    // 4. Reject if sequence is >80% X (effectively unknown)
    let x_count = cleaned.bytes().filter(|&b| b == b'X').count();
    let x_frac = x_count as f64 / cleaned.len() as f64;
    if x_frac > MAX_X_FRACTION {
        return PrepOutcome::Skip {
            reason: format!(
                "{:.0}% of residues are unknown (X) after cleaning — too ambiguous for InterProScan",
                x_frac * 100.0
            ),
        };
    }

    // 5. Reject if it looks like a nucleotide sequence (>80% A/T/C/G)
    let atcg = cleaned
        .bytes()
        .filter(|&b| matches!(b, b'A' | b'T' | b'C' | b'G'))
        .count();
    if atcg as f64 / cleaned.len() as f64 > 0.80 {
        return PrepOutcome::Skip {
            reason: format!(
                "sequence looks like DNA/RNA (>80% ATCG) — protein expected. \
                 Gene {} may not have been correctly translated.",
                protein_id
            ),
        };
    }

    let fix_note = if notes.is_empty() {
        None
    } else {
        Some(notes.join("; "))
    };
    PrepOutcome::Ok { cleaned, fix_note }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cache
// ─────────────────────────────────────────────────────────────────────────────

fn cache_path() -> PathBuf {
    let home = std::env::var("HOME").unwrap_or_else(|_| "/tmp".to_string());
    PathBuf::from(home)
        .join(".myconote")
        .join("iprscan_cache.tsv")
}

fn load_cache() -> HashMap<String, Vec<IprScanResult>> {
    let mut map: HashMap<String, Vec<IprScanResult>> = HashMap::new();
    if let Ok(content) = std::fs::read_to_string(cache_path()) {
        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.splitn(3, '\t').collect();
            if fields.len() >= 2 && fields[1] == NO_HITS_TAG {
                map.entry(fields[0].to_string()).or_default();
                continue;
            }
            if let Some(r) = parse_cache_line(line) {
                map.entry(r.protein_id.clone()).or_default().push(r);
            }
        }
    }
    map
}

fn parse_cache_line(line: &str) -> Option<IprScanResult> {
    let f: Vec<&str> = line.splitn(12, '\t').collect();
    if f.len() < 10 {
        return None;
    }
    Some(IprScanResult {
        protein_id: f[0].to_string(),
        ipr_accession: if f[1].is_empty() {
            None
        } else {
            Some(f[1].to_string())
        },
        ipr_desc: if f[2].is_empty() {
            None
        } else {
            Some(f[2].to_string())
        },
        database: f[3].to_string(),
        db_accession: f[4].to_string(),
        db_desc: f[5].to_string(),
        go_terms: f[6]
            .split('|')
            .filter(|s| !s.is_empty())
            .map(String::from)
            .collect(),
        pathways: f[7]
            .split('|')
            .filter(|s| !s.is_empty())
            .map(String::from)
            .collect(),
        evalue: f[8].parse().ok(),
        start: f[9].parse().unwrap_or(0),
        end: f.get(10).and_then(|s| s.parse().ok()).unwrap_or(0),
    })
}

fn save_cache(results: &[IprScanResult], submitted_ids: &[String]) {
    let path = cache_path();
    if let Some(p) = path.parent() {
        let _ = std::fs::create_dir_all(p);
    }

    let ids_with_hits: std::collections::HashSet<&str> =
        results.iter().map(|r| r.protein_id.as_str()).collect();

    if let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
    {
        for r in results {
            let _ = writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
        for id in submitted_ids {
            if !ids_with_hits.contains(id.as_str()) {
                let _ = writeln!(f, "{}\t{}", id, NO_HITS_TAG);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Binary-search isolation for batches that still fail after pre-processing
// ─────────────────────────────────────────────────────────────────────────────

/// Submit a batch. On HTTP 400, split in half and recurse until the bad
/// sequence(s) are isolated (base case: single sequence → diagnose and skip).
/// Returns (list of (ids, job_id) for good sub-batches, list of bad proteins).
fn isolate_and_submit(
    sequences: &[(String, String)],
    email: &str,
) -> (Vec<(Vec<String>, String)>, Vec<BadProtein>) {
    if sequences.is_empty() {
        return (vec![], vec![]);
    }

    std::thread::sleep(Duration::from_secs(2)); // brief rate-limit pause

    match submit_job(sequences, email) {
        Ok(job_id) => {
            let ids = sequences.iter().map(|(id, _)| id.clone()).collect();
            (vec![(ids, job_id)], vec![])
        }
        Err(e) if e.to_string().contains("status code 400") => {
            if sequences.len() == 1 {
                // Isolate the culprit — re-run pre-processing for final diagnosis
                let reason = match preprocess_sequence(&sequences[0].0, &sequences[0].1) {
                    PrepOutcome::Skip { reason } => reason,
                    PrepOutcome::Ok { .. } => {
                        "EBI rejected the sequence (HTTP 400) despite passing local checks; \
                         may contain very unusual residues"
                            .to_string()
                    }
                };
                return (
                    vec![],
                    vec![BadProtein {
                        protein_id: sequences[0].0.clone(),
                        reason,
                    }],
                );
            }
            let mid = sequences.len() / 2;
            let (mut jobs, mut bad) = isolate_and_submit(&sequences[..mid], email);
            let (j2, b2) = isolate_and_submit(&sequences[mid..], email);
            jobs.extend(j2);
            bad.extend(b2);
            (jobs, bad)
        }
        Err(e) => {
            // Non-400 error — report all as failed (transient; retry next run)
            let reason = format!("submit error: {}", e);
            let bad = sequences
                .iter()
                .map(|(id, _)| BadProtein {
                    protein_id: id.clone(),
                    reason: reason.clone(),
                })
                .collect();
            (vec![], bad)
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run(
    proteins_fa: &Path,
    output_tsv: &Path,
    email: &str,
) -> Result<HashMap<String, Vec<IprScanResult>>> {
    let raw_proteins = load_protein_fasta(proteins_fa)?;
    println!("  InterProScan: {} proteins to search", raw_proteins.len());

    // ── Step 0: Pre-validate and auto-fix sequences ──────────────────────────
    let mut proteins: Vec<(String, String)> = Vec::new(); // ready to submit
    let mut pre_bad: Vec<BadProtein> = Vec::new(); // unfixable
    let mut fixed_count = 0usize;

    for (id, seq) in &raw_proteins {
        match preprocess_sequence(id, seq) {
            PrepOutcome::Ok { cleaned, fix_note } => {
                if let Some(note) = fix_note {
                    println!("  ⚡ Auto-fixed {}: {}", id, note);
                    fixed_count += 1;
                }
                proteins.push((id.clone(), cleaned));
            }
            PrepOutcome::Skip { reason } => {
                eprintln!("  ✗  Skipping {} — {}", id, reason);
                pre_bad.push(BadProtein {
                    protein_id: id.clone(),
                    reason,
                });
            }
        }
    }

    if fixed_count > 0 || !pre_bad.is_empty() {
        println!(
            "  Pre-processing: {} auto-fixed, {} skipped (unfixable)",
            fixed_count,
            pre_bad.len()
        );
    }

    // ── Cache lookup ─────────────────────────────────────────────────────────
    let cache = load_cache();
    let mut all_results: HashMap<String, Vec<IprScanResult>> = HashMap::new();
    let mut missing: Vec<(String, String)> = Vec::new();

    for (id, seq) in &proteins {
        if let Some(cached) = cache.get(id) {
            all_results.insert(id.clone(), cached.clone());
        } else {
            missing.push((id.clone(), seq.clone()));
        }
    }
    println!(
        "  {} cached, {} to submit to EBI API",
        all_results.len(),
        missing.len()
    );

    let mut api_bad: Vec<BadProtein> = Vec::new(); // failures during API calls

    if !missing.is_empty() {
        let total_batches = (missing.len() + BATCH_SIZE - 1) / BATCH_SIZE;
        println!(
            "  Submitting {} batches ({}s apart), then polling all concurrently…",
            total_batches, SUBMIT_DELAY
        );

        let chunks: Vec<Vec<(String, String)>> =
            missing.chunks(BATCH_SIZE).map(|c| c.to_vec()).collect();

        // ── Phase 1: Rate-limited sequential submission ──────────────────────
        let mut good_jobs: Vec<(Vec<String>, String)> = Vec::new();

        for (idx, chunk) in chunks.iter().enumerate() {
            if idx > 0 {
                std::thread::sleep(Duration::from_secs(SUBMIT_DELAY));
            }
            let ids: Vec<String> = chunk.iter().map(|(id, _)| id.clone()).collect();
            println!(
                "  Batch {}/{}: submitting {} sequences…",
                idx + 1,
                total_batches,
                chunk.len()
            );

            match submit_with_retry(chunk, email) {
                Ok(job_id) => {
                    good_jobs.push((ids, job_id));
                }
                Err(ref e) if e.contains("status code 400") => {
                    eprintln!("  ⚠  Batch {} HTTP 400 — running binary search to isolate bad sequence(s)…",
                        idx + 1);
                    let (sub_jobs, bad) = isolate_and_submit(chunk, email);
                    let n_rescued = ids.len().saturating_sub(bad.len());
                    eprintln!(
                        "  ⚠  Batch {}: {} problem sequence(s) found, {} rescued",
                        idx + 1,
                        bad.len(),
                        n_rescued
                    );
                    for bp in &bad {
                        eprintln!("      ✗  {} — {}", bp.protein_id, bp.reason);
                    }
                    good_jobs.extend(sub_jobs);
                    api_bad.extend(bad);
                }
                Err(e) => {
                    eprintln!(
                        "  ✗  Batch {} failed ({}); will retry on next run.",
                        idx + 1,
                        e
                    );
                }
            }
        }

        let total_submitted = good_jobs.len();
        println!(
            "  {} job(s) queued on EBI. Polling all concurrently…",
            total_submitted
        );

        // ── Phase 2: Concurrent polling with live progress counter ───────────
        let done_count = Arc::new(AtomicUsize::new(0));
        let fail_count = Arc::new(AtomicUsize::new(0));
        let total_arc = total_submitted;

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(total_submitted.min(64).max(1))
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

        // Poll all jobs in parallel; collect (ids, result)
        let batch_results: Vec<(Vec<String>, std::result::Result<Vec<IprScanResult>, String>)> =
            pool.install(|| {
                good_jobs
                    .par_iter()
                    .map(|(ids, job_id)| {
                        let result = poll_and_fetch(job_id).map_err(|e| e.to_string());
                        (ids.clone(), result)
                    })
                    .collect()
            });

        // Collect results; isolate bad sequences from ERROR jobs
        for (ids, result) in batch_results {
            match result {
                Ok(hits) => {
                    save_cache(&hits, &ids);
                    let d = done_count.fetch_add(1, Ordering::Relaxed) + 1;
                    let f = fail_count.load(Ordering::Relaxed);
                    println!(
                        "  ✓  [{}/{} done | {} failed] {} hits",
                        d,
                        total_arc,
                        f,
                        hits.len()
                    );
                    for r in hits {
                        all_results.entry(r.protein_id.clone()).or_default().push(r);
                    }
                }
                Err(ref e) => {
                    let f = fail_count.fetch_add(1, Ordering::Relaxed) + 1;
                    let d = done_count.load(Ordering::Relaxed);

                    let is_error_status = e.contains("returned ERROR") && e.contains("consecutive");

                    if is_error_status {
                        eprintln!("  ⚠  [{}/{} done | {} failed] Persistent ERROR — binary-searching for bad sequences…",
                            d, total_arc, f);

                        // Look up the original cleaned sequences for this batch
                        let chunk: Vec<(String, String)> = ids
                            .iter()
                            .filter_map(|id| missing.iter().find(|(mid, _)| mid == id).cloned())
                            .collect();

                        let (sub_jobs, bad) = isolate_and_submit(&chunk, email);
                        let n_rescued = ids.len().saturating_sub(bad.len());
                        if !bad.is_empty() {
                            eprintln!(
                                "      {} bad sequence(s) found, {} rescued",
                                bad.len(),
                                n_rescued
                            );
                        }
                        for bp in &bad {
                            eprintln!("      ✗  {} — {}", bp.protein_id, bp.reason);
                        }
                        api_bad.extend(bad);

                        // Poll rescued sub-jobs synchronously (they're small)
                        for (sub_ids, sub_job_id) in sub_jobs {
                            match poll_and_fetch(&sub_job_id) {
                                Ok(hits) => {
                                    save_cache(&hits, &sub_ids);
                                    done_count.fetch_add(1, Ordering::Relaxed);
                                    // Adjust failure counter: one failure replaced by success
                                    fail_count.fetch_sub(
                                        fail_count.load(Ordering::Relaxed).min(1),
                                        Ordering::Relaxed,
                                    );
                                    for r in hits {
                                        all_results
                                            .entry(r.protein_id.clone())
                                            .or_default()
                                            .push(r);
                                    }
                                }
                                Err(e2) => {
                                    eprintln!("      ✗  Rescued sub-batch still failed: {}", e2);
                                }
                            }
                        }
                    } else {
                        eprintln!("  ✗  [{}/{} done | {} failed] {}", d, total_arc, f, e);
                    }
                }
            }
        }
    }

    // ── Final rejection report ───────────────────────────────────────────────
    let all_bad: Vec<&BadProtein> = pre_bad.iter().chain(api_bad.iter()).collect();

    if !all_bad.is_empty() {
        println!();
        println!(
            "  ┌─ Rejected proteins ({} total) ────────────────────────────",
            all_bad.len()
        );

        // Group by reason category
        let mut by_reason: HashMap<&str, Vec<&str>> = HashMap::new();
        for bp in &all_bad {
            by_reason
                .entry(bp.reason.as_str())
                .or_default()
                .push(bp.protein_id.as_str());
        }
        let mut reasons: Vec<(&str, Vec<&str>)> = by_reason.into_iter().collect();
        reasons.sort_by_key(|(r, _)| *r);

        for (reason, ids) in &reasons {
            println!("  │");
            println!("  │  Reason: {}", reason);
            println!("  │  Affected ({}):", ids.len());
            for id in ids.iter().take(10) {
                println!("  │    • {}", id);
            }
            if ids.len() > 10 {
                println!("  │    … and {} more", ids.len() - 10);
            }
        }
        println!("  │");
        println!("  └────────────────────────────────────────────────────────────");
    }

    write_iprscan_tsv(output_tsv, &all_results)?;
    println!(
        "  InterProScan: {}/{} proteins with at least one hit",
        all_results.values().filter(|v| !v.is_empty()).count(),
        raw_proteins.len()
    );
    Ok(all_results)
}

// ─────────────────────────────────────────────────────────────────────────────
// Submit helpers
// ─────────────────────────────────────────────────────────────────────────────

fn submit_with_retry(
    sequences: &[(String, String)],
    email: &str,
) -> std::result::Result<String, String> {
    let mut attempt = 0usize;
    loop {
        match submit_job(sequences, email) {
            Ok(id) => return Ok(id),
            Err(e) => {
                let msg = e.to_string();
                if msg.contains("status code 400") {
                    return Err(msg);
                }
                let (max_try, delay) = if msg.contains("Dns Failed") || msg.contains("dns") {
                    (
                        MAX_DNS_RETRIES,
                        DNS_RETRY_DELAY * (1 + attempt as u64).min(3),
                    )
                } else {
                    (MAX_RETRIES, SUBMIT_DELAY * (1u64 << attempt.min(4)))
                };
                attempt += 1;
                if attempt >= max_try {
                    return Err(msg);
                }
                eprintln!(
                    "  ⚠  Submit attempt {} failed: {} — retrying in {}s…",
                    attempt, msg, delay
                );
                std::thread::sleep(Duration::from_secs(delay));
            }
        }
    }
}

fn submit_job(sequences: &[(String, String)], email: &str) -> Result<String> {
    let mut fasta = String::new();
    for (id, seq) in sequences {
        if seq.is_empty() {
            continue;
        }
        fasta.push('>');
        fasta.push_str(id);
        fasta.push('\n');
        fasta.push_str(seq); // already cleaned by preprocess_sequence
        fasta.push('\n');
    }

    let body = format!(
        "email={}&title=myconote&goterms=true&pathways=true&sequence={}",
        url_encode(email),
        url_encode(&fasta)
    );

    let resp = ureq::post(&format!("{}/run", API_BASE))
        .set("Content-Type", "application/x-www-form-urlencoded")
        .timeout(Duration::from_secs(30))
        .send_string(&body)
        .map_err(|e| MycoNoteError::InvalidFormat(format!("InterProScan submit failed: {}", e)))?;

    let job_id = resp
        .into_string()
        .map_err(|e| MycoNoteError::InvalidFormat(format!("Response read error: {}", e)))?
        .trim()
        .to_string();

    if job_id.is_empty() || job_id.to_lowercase().contains("error") {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Bad job ID from EBI: '{}'",
            job_id
        )));
    }
    Ok(job_id)
}

// ─────────────────────────────────────────────────────────────────────────────
// Poll + fetch
// ─────────────────────────────────────────────────────────────────────────────

fn poll_and_fetch(job_id: &str) -> Result<Vec<IprScanResult>> {
    poll_job(job_id)?;
    let tsv = fetch_result_tsv(job_id)?;
    Ok(parse_iprscan_tsv(&tsv))
}

fn poll_job(job_id: &str) -> Result<()> {
    let url = format!("{}/status/{}", API_BASE, job_id);
    let mut consec_errors = 0usize;
    let mut dns_retries = 0usize;

    for poll in 0..MAX_POLLS {
        let status = match ureq::get(&url).timeout(Duration::from_secs(30)).call() {
            Ok(r) => r.into_string().unwrap_or_default().trim().to_string(),
            Err(e) => {
                let msg = e.to_string();
                if msg.contains("Dns Failed") || msg.contains("dns") {
                    dns_retries += 1;
                    if dns_retries >= MAX_DNS_RETRIES {
                        return Err(MycoNoteError::InvalidFormat(format!(
                            "DNS failed {} times for {}: {}",
                            dns_retries, job_id, msg
                        )));
                    }
                    let wait = (DNS_RETRY_DELAY * (dns_retries as u64).min(3)).min(60);
                    eprintln!("  ⚠  DNS error ({}×) — retrying in {}s…", dns_retries, wait);
                    std::thread::sleep(Duration::from_secs(wait));
                    continue;
                }
                eprintln!("  ⚠  Poll {} for {} failed: {}", poll + 1, job_id, msg);
                std::thread::sleep(Duration::from_secs(POLL_DELAY));
                continue;
            }
        };

        dns_retries = 0;
        match status.as_str() {
            "FINISHED" => return Ok(()),
            "RUNNING" | "QUEUED" => {
                consec_errors = 0;
                if poll % 4 == 0 && poll > 0 {
                    let min = (poll as u64 * POLL_DELAY) / 60;
                    println!("  … {} {} ({} min)", job_id, status, min);
                }
            }
            "FAILED" => {
                return Err(MycoNoteError::InvalidFormat(format!(
                    "InterProScan job {} FAILED on EBI servers.",
                    job_id
                )))
            }
            "NOT_FOUND" => {
                return Err(MycoNoteError::InvalidFormat(format!(
                    "InterProScan job {} not found (expired?).",
                    job_id
                )))
            }
            "ERROR" => {
                consec_errors += 1;
                if consec_errors % 3 == 1 {
                    eprintln!("  ⚠  Job {} status: ERROR ({}×)…", job_id, consec_errors);
                }
                if consec_errors >= MAX_CONSEC_ERRORS {
                    return Err(MycoNoteError::InvalidFormat(format!(
                        "InterProScan job {} returned ERROR {} consecutive times. \
                         Likely caused by invalid sequence characters. Skipping.",
                        job_id, consec_errors
                    )));
                }
            }
            other => eprintln!("  ⚠  Job {} unexpected status: {}", job_id, other),
        }
        std::thread::sleep(Duration::from_secs(POLL_DELAY));
    }

    Err(MycoNoteError::InvalidFormat(format!(
        "InterProScan job {} timed out after {} polls.",
        job_id, MAX_POLLS
    )))
}

fn fetch_result_tsv(job_id: &str) -> Result<String> {
    ureq::get(&format!("{}/result/{}/tsv", API_BASE, job_id))
        .timeout(Duration::from_secs(60))
        .call()
        .map_err(|e| MycoNoteError::InvalidFormat(format!("Result fetch failed: {}", e)))?
        .into_string()
        .map_err(|e| MycoNoteError::InvalidFormat(format!("Result read error: {}", e)))
}

// ─────────────────────────────────────────────────────────────────────────────
// TSV parser
// ─────────────────────────────────────────────────────────────────────────────

fn parse_iprscan_tsv(tsv: &str) -> Vec<IprScanResult> {
    let mut out = Vec::new();
    for line in tsv.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 || f.get(9) == Some(&"FALSE") {
            continue;
        }
        out.push(IprScanResult {
            protein_id: f[0].to_string(),
            database: f[3].to_string(),
            db_accession: f[4].to_string(),
            db_desc: f[5].to_string(),
            start: f[6].parse().unwrap_or(0),
            end: f[7].parse().unwrap_or(0),
            evalue: f[8].parse().ok(),
            ipr_accession: f.get(11).and_then(|s| {
                if s.is_empty() || *s == "-" {
                    None
                } else {
                    Some(s.to_string())
                }
            }),
            ipr_desc: f.get(12).and_then(|s| {
                if s.is_empty() || *s == "-" {
                    None
                } else {
                    Some(s.to_string())
                }
            }),
            go_terms: f
                .get(13)
                .map(|s| {
                    s.split('|')
                        .filter(|t| t.starts_with("GO:"))
                        .map(String::from)
                        .collect()
                })
                .unwrap_or_default(),
            pathways: f
                .get(14)
                .map(|s| {
                    s.split('|')
                        .filter(|t| !t.trim().is_empty())
                        .map(String::from)
                        .collect()
                })
                .unwrap_or_default(),
        });
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// TSV writer
// ─────────────────────────────────────────────────────────────────────────────

fn write_iprscan_tsv(path: &Path, results: &HashMap<String, Vec<IprScanResult>>) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "protein_id\tipr_accession\tipr_desc\tdatabase\tdb_accession\tdb_desc\tgo_terms\tpathways\tevalue\tstart\tend")
        .map_err(MycoNoteError::Io)?;
    let mut ids: Vec<&String> = results.keys().collect();
    ids.sort();
    for id in ids {
        for r in &results[id] {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
            )
            .map_err(MycoNoteError::Io)?;
        }
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTA loader
// ─────────────────────────────────────────────────────────────────────────────

fn load_protein_fasta(path: &Path) -> Result<Vec<(String, String)>> {
    use std::io::BufRead;
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut out = Vec::new();
    let mut cur_id = String::new();
    let mut cur_seq = String::new();
    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim();
        if line.starts_with('>') {
            if !cur_id.is_empty() && !cur_seq.is_empty() {
                out.push((cur_id.clone(), cur_seq.clone()));
            }
            cur_id = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            cur_seq = String::new();
        } else {
            cur_seq.push_str(line);
        }
    }
    if !cur_id.is_empty() && !cur_seq.is_empty() {
        out.push((cur_id, cur_seq));
    }
    Ok(out)
}

// ─────────────────────────────────────────────────────────────────────────────
// URL encoding
// ─────────────────────────────────────────────────────────────────────────────

fn url_encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len() * 2);
    for byte in s.bytes() {
        match byte {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => {
                out.push(byte as char)
            }
            b' ' => out.push('+'),
            other => {
                out.push('%');
                out.push_str(&format!("{:02X}", other));
            }
        }
    }
    out
}
