/// GO term assignment from UniProt accessions
///
/// Fetches Gene Ontology annotations for UniProt accessions using the
/// UniProt REST API.  Results are cached locally so the API is only
/// queried once per accession.
///
/// Cache: ~/.myconote/go_cache.tsv  (tab: accession → GO:XXXXXXX,...)

use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Cache management
// ─────────────────────────────────────────────────────────────────────────────

fn cache_path() -> PathBuf {
    let home = std::env::var("HOME").unwrap_or_else(|_| "/tmp".to_string());
    PathBuf::from(home).join(".myconote").join("go_cache.tsv")
}

fn load_cache() -> HashMap<String, Vec<String>> {
    let path = cache_path();
    let mut map = HashMap::new();
    if let Ok(content) = std::fs::read_to_string(&path) {
        for line in content.lines() {
            let mut parts = line.splitn(2, '\t');
            if let (Some(acc), Some(terms)) = (parts.next(), parts.next()) {
                let go_terms: Vec<String> = terms.split(',')
                    .filter(|t| !t.is_empty())
                    .map(String::from)
                    .collect();
                map.insert(acc.to_string(), go_terms);
            }
        }
    }
    map
}

fn save_cache(map: &HashMap<String, Vec<String>>) {
    let path = cache_path();
    if let Some(parent) = path.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    if let Ok(mut f) = std::fs::File::create(&path) {
        for (acc, terms) in map {
            let _ = writeln!(f, "{}\t{}", acc, terms.join(","));
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Given a list of UniProt accessions, return a map of accession → GO term list.
/// Writes results to `output_tsv` and updates the local cache.
pub fn assign_go_terms(
    accs:       &[String],
    output_tsv: &Path,
) -> Result<HashMap<String, Vec<String>>> {
    if accs.is_empty() {
        return Ok(HashMap::new());
    }

    println!("  Assigning GO terms for {} UniProt accessions…", accs.len());

    let mut cache = load_cache();
    let mut result: HashMap<String, Vec<String>> = HashMap::new();
    let mut missing: Vec<String> = Vec::new();

    for acc in accs {
        if let Some(terms) = cache.get(acc) {
            result.insert(acc.clone(), terms.clone());
        } else {
            missing.push(acc.clone());
        }
    }

    if !missing.is_empty() {
        println!("  Fetching GO terms for {} uncached accessions…", missing.len());
        let fetched = fetch_go_terms_batch(&missing)?;
        for (acc, terms) in &fetched {
            cache.insert(acc.clone(), terms.clone());
            result.insert(acc.clone(), terms.clone());
        }
        // Fill empty entries so we don't re-query known misses
        for acc in &missing {
            cache.entry(acc.clone()).or_default();
        }
        save_cache(&cache);
    }

    write_go_tsv(output_tsv, &result)?;

    Ok(result)
}

// ─────────────────────────────────────────────────────────────────────────────
// UniProt REST API
// ─────────────────────────────────────────────────────────────────────────────

const BATCH_SIZE: usize = 50;
const DELAY_MS:   u64   = 400;

fn fetch_go_terms_batch(accs: &[String]) -> Result<HashMap<String, Vec<String>>> {
    let mut all: HashMap<String, Vec<String>> = HashMap::new();

    for chunk in accs.chunks(BATCH_SIZE) {
        let ids = chunk.join(",");
        let url = format!(
            "https://rest.uniprot.org/uniprotkb/accessions?accessions={}&fields=accession,go_id&format=tsv",
            ids
        );

        match ureq::get(&url).call() {
            Ok(resp) => {
                match resp.into_string() {
                    Ok(body) => parse_uniprot_go_tsv(&body, &mut all),
                    Err(e) => eprintln!("  ⚠  UniProt GO read error: {}", e),
                }
            }
            Err(e) => {
                eprintln!("  ⚠  UniProt GO fetch failed: {}", e);
            }
        }

        std::thread::sleep(std::time::Duration::from_millis(DELAY_MS));
    }

    Ok(all)
}

/// Parse UniProt TSV response (columns: Accession, Gene Ontology IDs)
fn parse_uniprot_go_tsv(body: &str, out: &mut HashMap<String, Vec<String>>) {
    for line in body.lines().skip(1) {  // skip header
        let mut parts = line.splitn(2, '\t');
        let acc = match parts.next() { Some(a) if !a.is_empty() => a, _ => continue };
        let go_raw = parts.next().unwrap_or("");

        let go_terms: Vec<String> = go_raw.split(';')
            .map(|t| t.trim().to_string())
            .filter(|t| t.starts_with("GO:"))
            .collect();

        out.insert(acc.to_string(), go_terms);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Output
// ─────────────────────────────────────────────────────────────────────────────

fn write_go_tsv(path: &Path, go_map: &HashMap<String, Vec<String>>) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "uniprot_acc\tgo_terms").map_err(MycoNoteError::Io)?;

    let mut entries: Vec<_> = go_map.iter().collect();
    entries.sort_by_key(|(k, _)| k.as_str());

    for (acc, terms) in entries {
        writeln!(f, "{}\t{}", acc, terms.join(";"))
            .map_err(MycoNoteError::Io)?;
    }

    Ok(())
}
