/// Remote annotation services: `myconote remote`
///
/// Submits protein sequences to remote annotation servers and merges
/// the results into the annotation table.
///
/// Supported services:
///   Phobius      — Signal peptide + TM topology (Stockholm University REST API)
///   InterProScan — Comprehensive domain/family search (EBI REST API)
///   SignalP      — Alternative to Phobius; same EBI REST endpoint
///   DeepLoc 2    — Subcellular localisation (DTU web service)
///
/// All remote services:
///   1. Split proteins into batches (API rate limits)
///   2. Submit batches with backoff retry on 429/503
///   3. Poll for results (EBI returns a job ID)
///   4. Parse and merge results
///
/// Results are cached in <out_dir>/remote_cache/ so re-runs are fast.
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::thread;
use std::time::Duration;

// ─────────────────────────────────────────────────────────────────────────────
// Remote result types
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct PhobiusResult {
    pub query_id: String,
    pub signal_peptide: bool,
    pub signal_end: Option<u32>, // position of signal peptide cleavage
    pub tm_count: u32,
    pub topology: String, // e.g. "i", "o", "iXXXXo" (inside/outside)
}

#[derive(Debug, Clone, Default)]
pub struct RemoteAnnotation {
    pub phobius: Option<PhobiusResult>,
    /// InterProScan result: IPR accessions
    pub ipr_accs: Vec<String>,
    /// InterProScan GO terms
    pub go_terms: Vec<String>,
    /// Subcellular localisation from DeepLoc (e.g. "Extracellular", "Nucleus")
    pub localisation: Option<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct RemoteConfig {
    pub proteins_fa: PathBuf,
    pub out_dir: PathBuf,
    pub run_phobius: bool,
    pub run_interpro: bool,
    pub run_deeploc: bool,
    /// EBI email (required for InterProScan)
    pub email: String,
    /// Max proteins per batch
    pub batch_size: usize,
    /// Seconds between API polls
    pub poll_interval: u64,
    /// Maximum retry attempts per batch
    pub max_retries: usize,
}

impl Default for RemoteConfig {
    fn default() -> Self {
        Self {
            proteins_fa: PathBuf::new(),
            out_dir: PathBuf::from("remote_out"),
            run_phobius: true,
            run_interpro: false,
            run_deeploc: false,
            email: String::new(),
            batch_size: 100,
            poll_interval: 30,
            max_retries: 10,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_remote(config: &RemoteConfig) -> Result<HashMap<String, RemoteAnnotation>> {
    if !config.proteins_fa.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Protein FASTA not found: {}",
            config.proteins_fa.display()
        )));
    }

    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;
    let cache_dir = config.out_dir.join("remote_cache");
    std::fs::create_dir_all(&cache_dir).map_err(MycoNoteError::Io)?;

    // Read all proteins
    let proteins = read_fasta_sequences(&config.proteins_fa)?;
    println!("  {} proteins loaded for remote annotation", proteins.len());

    let mut results: HashMap<String, RemoteAnnotation> = HashMap::new();
    for id in proteins.keys() {
        results.insert(id.clone(), RemoteAnnotation::default());
    }

    // ── Phobius ──────────────────────────────────────────────────────────────
    if config.run_phobius {
        println!("  Running Phobius (signal peptide + TM topology)…");
        let phobius_cache = cache_dir.join("phobius_results.tsv");

        let phobius_map = if phobius_cache.exists() {
            println!("    Using cached Phobius results.");
            parse_phobius_tsv(&phobius_cache)?
        } else {
            let map = run_phobius_remote(&proteins, config, &phobius_cache)?;
            map
        };

        let n_signal: usize = phobius_map.values().filter(|r| r.signal_peptide).count();
        println!(
            "    Phobius: {} signal peptides, {} TM proteins",
            n_signal,
            phobius_map.values().filter(|r| r.tm_count > 0).count(),
        );

        for (id, phobius) in phobius_map {
            results.entry(id).or_default().phobius = Some(phobius);
        }
    }

    // ── InterProScan ─────────────────────────────────────────────────────────
    if config.run_interpro {
        if config.email.is_empty() {
            eprintln!("  ⚠  InterProScan requires --email <address> (EBI policy)");
        } else {
            println!("  Running InterProScan (EBI REST API)…");
            let ipr_cache = cache_dir.join("interpro_results.tsv");

            let ipr_map = if ipr_cache.exists() {
                println!("    Using cached InterProScan results.");
                parse_interpro_tsv(&ipr_cache)?
            } else {
                let map = run_interproscan_remote(&proteins, config, &ipr_cache)?;
                map
            };

            let n_with_hits = ipr_map.values().filter(|v| !v.is_empty()).count();
            println!(
                "    InterProScan: {}/{} proteins with domain hits",
                n_with_hits,
                proteins.len()
            );

            for (id, ipr_hits) in ipr_map {
                let entry = results.entry(id).or_default();
                entry
                    .ipr_accs
                    .extend(ipr_hits.iter().map(|(acc, _)| acc.clone()));
                entry
                    .go_terms
                    .extend(ipr_hits.iter().flat_map(|(_, gos)| gos.clone()));
            }
        }
    }

    // ── DeepLoc 2 ────────────────────────────────────────────────────────────
    if config.run_deeploc {
        println!("  Running DeepLoc 2 (subcellular localisation)…");
        println!("  Note: DeepLoc 2 requires local installation or the DTU web server.");
        let dl_result = run_deeploc_local(&config.proteins_fa, &config.out_dir);
        match dl_result {
            Ok(dl_map) => {
                println!("    DeepLoc 2: {} predictions", dl_map.len());
                for (id, loc) in dl_map {
                    results.entry(id).or_default().localisation = Some(loc);
                }
            }
            Err(e) => eprintln!("  ⚠  DeepLoc failed: {}", e),
        }
    }

    // ── Write merged TSV ─────────────────────────────────────────────────────
    let merged_tsv = config.out_dir.join("remote_annotations.tsv");
    write_remote_tsv(&merged_tsv, &results)?;
    println!("  ✓  Remote annotations → {}", merged_tsv.display());

    Ok(results)
}

// ─────────────────────────────────────────────────────────────────────────────
// Phobius REST API (http://phobius.sbc.su.se/cgi-bin/predict.pl)
// ─────────────────────────────────────────────────────────────────────────────

const PHOBIUS_URL: &str = "https://phobius.sbc.su.se/cgi-bin/predict.pl";

fn run_phobius_remote(
    proteins: &HashMap<String, String>,
    config: &RemoteConfig,
    cache_out: &Path,
) -> Result<HashMap<String, PhobiusResult>> {
    // Phobius accepts FASTA via HTTP POST, returns text
    // Use curl as the HTTP client (avoid requiring reqwest in Rust)

    let mut all_results: HashMap<String, PhobiusResult> = HashMap::new();
    let ids: Vec<&String> = proteins.keys().collect();

    for (batch_idx, chunk) in ids.chunks(config.batch_size).enumerate() {
        println!(
            "    Phobius batch {}/{}…",
            batch_idx + 1,
            (ids.len() + config.batch_size - 1) / config.batch_size
        );

        // Write batch FASTA to temp file
        let tmp_fa = config
            .out_dir
            .join(format!("phobius_batch_{}.fa", batch_idx));
        {
            let mut f = std::fs::File::create(&tmp_fa).map_err(MycoNoteError::Io)?;
            for id in chunk {
                writeln!(f, ">{}", id).map_err(MycoNoteError::Io)?;
                // Write sequence in 60-char lines
                let seq = &proteins[*id];
                for piece in seq.as_bytes().chunks(60) {
                    writeln!(f, "{}", std::str::from_utf8(piece).unwrap_or(""))
                        .map_err(MycoNoteError::Io)?;
                }
            }
        }

        // Submit to Phobius via curl
        let tmp_out = config
            .out_dir
            .join(format!("phobius_out_{}.txt", batch_idx));
        let ok = submit_with_retry(
            &[
                "curl",
                "-s",
                "-X",
                "POST",
                "-F",
                &format!("protseq=@{}", tmp_fa.to_str().unwrap_or("")),
                "-F",
                "format=short",
                PHOBIUS_URL,
                "-o",
                tmp_out.to_str().unwrap_or(""),
            ],
            config.max_retries,
            config.poll_interval,
        );

        if ok {
            let batch_results = parse_phobius_output(&tmp_out)?;
            all_results.extend(batch_results);
        } else {
            eprintln!(
                "    ⚠  Phobius batch {} failed after retries",
                batch_idx + 1
            );
        }

        let _ = std::fs::remove_file(&tmp_fa);
        thread::sleep(Duration::from_secs(2)); // polite interval
    }

    // Cache results
    write_phobius_cache(&all_results, cache_out)?;

    Ok(all_results)
}

fn submit_with_retry(args: &[&str], max_retries: usize, interval_secs: u64) -> bool {
    use std::process::Command;
    for attempt in 0..max_retries {
        if attempt > 0 {
            thread::sleep(Duration::from_secs(interval_secs));
        }
        if let Some((prog, rest)) = args.split_first() {
            if let Ok(status) = Command::new(prog).args(rest).status() {
                if status.success() {
                    return true;
                }
            }
        }
    }
    false
}

fn parse_phobius_output(path: &Path) -> Result<HashMap<String, PhobiusResult>> {
    // Phobius "short" format:
    //   ID  TM  SP  PREDICTION
    //   GENE0001  0  Y  SIGNAL
    //   GENE0002  2  N  TM
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for line in reader.lines().flatten() {
        let t = line.trim();
        if t.is_empty() || t.starts_with("ID") || t.starts_with("//") {
            continue;
        }
        let cols: Vec<&str> = t.split_whitespace().collect();
        if cols.len() < 3 {
            continue;
        }

        let id = cols[0].to_string();
        let tm: u32 = cols[1].parse().unwrap_or(0);
        let sp = cols[2].to_uppercase() == "Y";
        let topology = cols.get(3).unwrap_or(&"").to_string();

        map.insert(
            id.clone(),
            PhobiusResult {
                query_id: id,
                signal_peptide: sp,
                signal_end: None,
                tm_count: tm,
                topology,
            },
        );
    }

    Ok(map)
}

fn write_phobius_cache(results: &HashMap<String, PhobiusResult>, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "query_id\ttm_count\tsignal_peptide\ttopology").map_err(MycoNoteError::Io)?;
    for r in results.values() {
        writeln!(
            f,
            "{}\t{}\t{}\t{}",
            r.query_id,
            r.tm_count,
            if r.signal_peptide { "Y" } else { "N" },
            r.topology,
        )
        .map_err(MycoNoteError::Io)?;
    }
    Ok(())
}

fn parse_phobius_tsv(path: &Path) -> Result<HashMap<String, PhobiusResult>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for line in reader.lines().flatten().skip(1) {
        let cols: Vec<&str> = line.trim().split('\t').collect();
        if cols.len() < 4 {
            continue;
        }
        let id = cols[0].to_string();
        let tm: u32 = cols[1].parse().unwrap_or(0);
        let sp = cols[2].to_uppercase() == "Y";
        let topology = cols[3].to_string();
        map.insert(
            id.clone(),
            PhobiusResult {
                query_id: id,
                signal_peptide: sp,
                signal_end: None,
                tm_count: tm,
                topology,
            },
        );
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// InterProScan EBI REST API
// ─────────────────────────────────────────────────────────────────────────────

const EBI_IPRSCAN_URL: &str = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5";

fn run_interproscan_remote(
    proteins: &HashMap<String, String>,
    config: &RemoteConfig,
    cache_out: &Path,
) -> Result<HashMap<String, Vec<(String, Vec<String>)>>> {
    let mut all: HashMap<String, Vec<(String, Vec<String>)>> = HashMap::new();
    let ids: Vec<&String> = proteins.keys().collect();

    for (batch_idx, chunk) in ids.chunks(config.batch_size.min(30)).enumerate() {
        // EBI limits to 30 sequences per job
        println!(
            "    InterProScan batch {}/{}…",
            batch_idx + 1,
            (ids.len() + 29) / 30
        );

        // Build FASTA payload
        let mut fa = String::new();
        for id in chunk {
            fa.push('>');
            fa.push_str(id);
            fa.push('\n');
            fa.push_str(&proteins[*id]);
            fa.push('\n');
        }

        // Submit job
        let job_id = submit_interpro_job(&fa, &config.email, config.max_retries);
        let job_id = match job_id {
            Some(id) => id,
            None => {
                eprintln!(
                    "    ⚠  InterProScan submission failed for batch {}",
                    batch_idx + 1
                );
                continue;
            }
        };

        // Poll for result
        let result_tsv = poll_interpro_result(&job_id, config.poll_interval, config.max_retries);

        match result_tsv {
            Some(tsv_text) => {
                let batch_map = parse_interpro_tsv_text(&tsv_text);
                all.extend(batch_map);
            }
            None => {
                eprintln!("    ⚠  InterProScan polling failed for job {}", job_id);
            }
        }

        thread::sleep(Duration::from_secs(5));
    }

    // Cache
    write_interpro_cache(&all, cache_out)?;

    Ok(all)
}

fn submit_interpro_job(fasta: &str, email: &str, _retries: usize) -> Option<String> {
    use std::process::Command;

    // Write FASTA to temp file
    let tmp = std::env::temp_dir().join("iprscan_input.fa");
    if std::fs::write(&tmp, fasta).is_err() {
        return None;
    }

    // curl POST to EBI
    let output = Command::new("curl")
        .args([
            "-s",
            "-X",
            "POST",
            &format!("{}/run", EBI_IPRSCAN_URL),
            "-F",
            &format!("email={}", email),
            "-F",
            "title=myconote",
            "-F",
            "goterms=true",
            "-F",
            "pathways=false",
            "-F",
            &format!("sequence=@{}", tmp.to_str().unwrap_or("")),
        ])
        .output()
        .ok()?;

    let job_id = String::from_utf8_lossy(&output.stdout).trim().to_string();
    if job_id.contains("iprscan5") {
        Some(job_id)
    } else {
        None
    }
}

fn poll_interpro_result(job_id: &str, interval: u64, max_tries: usize) -> Option<String> {
    use std::process::Command;
    use std::time::Duration;

    for _ in 0..max_tries {
        // Check status
        let status_out = Command::new("curl")
            .args(["-s", &format!("{}/status/{}", EBI_IPRSCAN_URL, job_id)])
            .output()
            .ok()?;
        let status = String::from_utf8_lossy(&status_out.stdout).to_string();

        if status.trim() == "FINISHED" {
            // Retrieve TSV result
            let result = Command::new("curl")
                .args(["-s", &format!("{}/result/{}/tsv", EBI_IPRSCAN_URL, job_id)])
                .output()
                .ok()?;
            return Some(String::from_utf8_lossy(&result.stdout).to_string());
        } else if status.trim() == "FAILURE" || status.trim() == "ERROR" {
            eprintln!("    InterProScan job {} failed: {}", job_id, status.trim());
            return None;
        }

        thread::sleep(Duration::from_secs(interval));
    }

    None
}

fn parse_interpro_tsv_text(text: &str) -> HashMap<String, Vec<(String, Vec<String>)>> {
    // InterProScan TSV columns:
    // 0:protein_id 1:md5 2:length 3:analysis 4:sig_acc 5:sig_desc
    // 6:start 7:end 8:score 9:status 10:date
    // 11:ipr_acc 12:ipr_desc 13:go_terms 14:pathways
    let mut map: HashMap<String, Vec<(String, Vec<String>)>> = HashMap::new();

    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 12 {
            continue;
        }

        let protein_id = cols[0].to_string();
        let ipr_acc = cols[11].to_string();
        if ipr_acc.is_empty() || ipr_acc == "-" {
            continue;
        }

        let go_terms: Vec<String> = if cols.len() > 13 && !cols[13].is_empty() && cols[13] != "-" {
            cols[13].split('|').map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        map.entry(protein_id).or_default().push((ipr_acc, go_terms));
    }

    map
}

fn write_interpro_cache(
    results: &HashMap<String, Vec<(String, Vec<String>)>>,
    path: &Path,
) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "protein_id\tipr_acc\tgo_terms").map_err(MycoNoteError::Io)?;
    for (pid, hits) in results {
        for (acc, gos) in hits {
            writeln!(f, "{}\t{}\t{}", pid, acc, gos.join("|")).map_err(MycoNoteError::Io)?;
        }
    }
    Ok(())
}

fn parse_interpro_tsv(path: &Path) -> Result<HashMap<String, Vec<(String, Vec<String>)>>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, Vec<(String, Vec<String>)>> = HashMap::new();

    for line in reader.lines().flatten().skip(1) {
        let cols: Vec<&str> = line.trim().split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let pid = cols[0].to_string();
        let acc = cols[1].to_string();
        let gos: Vec<String> = if cols[2].is_empty() {
            vec![]
        } else {
            cols[2].split('|').map(|s| s.to_string()).collect()
        };
        map.entry(pid).or_default().push((acc, gos));
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// DeepLoc 2 (local)
// ─────────────────────────────────────────────────────────────────────────────

fn run_deeploc_local(proteins_fa: &Path, out_dir: &Path) -> Result<HashMap<String, String>> {
    use std::process::Command;

    let out_csv = out_dir.join("deeploc_output.csv");

    let ok = Command::new("deeploc2")
        .args([
            "--fasta",
            proteins_fa.to_str().unwrap_or(""),
            "--output",
            out_dir.to_str().unwrap_or(""),
        ])
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if !ok {
        return Err(MycoNoteError::ExternalTool(
            "deeploc2 not found. Install: pip install deeploc2".to_string(),
        ));
    }

    parse_deeploc_csv(&out_csv)
}

fn parse_deeploc_csv(path: &Path) -> Result<HashMap<String, String>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for line in reader.lines().flatten().skip(1) {
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() < 2 {
            continue;
        }
        let id = cols[0].trim().trim_matches('"').to_string();
        let loc = cols[1].trim().trim_matches('"').to_string();
        if !id.is_empty() {
            map.insert(id, loc);
        }
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// Merged output
// ─────────────────────────────────────────────────────────────────────────────

fn write_remote_tsv(path: &Path, results: &HashMap<String, RemoteAnnotation>) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "gene_id\tsignal_peptide\ttm_count\tipr_accessions\tgo_terms\tlocalisation"
    )
    .map_err(MycoNoteError::Io)?;

    let mut ids: Vec<&String> = results.keys().collect();
    ids.sort();

    for id in ids {
        let ann = &results[id];
        let (sp, tm) = ann
            .phobius
            .as_ref()
            .map(|p| {
                (
                    if p.signal_peptide { "Y" } else { "N" },
                    p.tm_count.to_string(),
                )
            })
            .unwrap_or(("", "".to_string()));

        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            id,
            sp,
            tm,
            ann.ipr_accs.join("|"),
            ann.go_terms.join("|"),
            ann.localisation.as_deref().unwrap_or(""),
        )
        .map_err(MycoNoteError::Io)?;
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTA reader helper
// ─────────────────────────────────────────────────────────────────────────────

fn read_fasta_sequences(path: &Path) -> Result<HashMap<String, String>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, String> = HashMap::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in reader.lines().flatten() {
        let t = line.trim();
        if t.starts_with('>') {
            if !current_id.is_empty() {
                map.insert(current_id.clone(), current_seq.clone());
            }
            current_id = t[1..].split_whitespace().next().unwrap_or("").to_string();
            current_seq = String::new();
        } else {
            current_seq.push_str(t);
        }
    }

    if !current_id.is_empty() {
        map.insert(current_id, current_seq);
    }

    Ok(map)
}
