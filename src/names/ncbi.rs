/// NCBI Entrez Gene API backend
///
/// Uses the Entrez E-utilities (esearch + esummary) to resolve gene IDs
/// or gene symbols to official gene names.
///
/// Rate limit: 3 requests/second without API key, 10/sec with one.
/// We batch IDs and sleep between requests to stay compliant.
use std::collections::HashMap;

const ESEARCH_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
const ESUMMARY_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";
const BATCH_SIZE: usize = 50;

/// Fetch official gene names from NCBI Gene for a list of IDs/symbols.
///
/// `taxon_id` narrows the search to a specific organism (e.g. 5207 for
/// *Cryptococcus neoformans*).  Pass `None` to search all organisms.
///
/// Returns a map of `input_id → official_name`.
/// IDs that cannot be resolved are absent from the map.
pub fn fetch_gene_names(ids: &[String], taxon_id: Option<u32>) -> HashMap<String, String> {
    let mut result = HashMap::new();

    for chunk in ids.chunks(BATCH_SIZE) {
        match fetch_chunk(chunk, taxon_id) {
            Ok(partial) => result.extend(partial),
            Err(e) => eprintln!("    NCBI fetch error: {}", e),
        }
        // Polite delay — stay within 3 req/sec limit
        std::thread::sleep(std::time::Duration::from_millis(350));
    }

    result
}

fn fetch_chunk(
    ids: &[String],
    taxon_id: Option<u32>,
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let client = reqwest::blocking::Client::builder()
        .user_agent("myconote-cli/0.1 (genome annotation tool; contact via GitHub)")
        .timeout(std::time::Duration::from_secs(30))
        .build()?;

    let mut results = HashMap::new();

    for id in ids {
        // Build esearch query: gene symbol + optional taxon filter
        let query = match taxon_id {
            Some(txid) => format!("{}[Gene Name] AND {}[Taxonomy ID]", id, txid),
            None => format!("{}[Gene Name]", id),
        };

        let search_url = format!(
            "{}?db=gene&term={}&retmode=json&retmax=1",
            ESEARCH_URL,
            urlencoding(&query)
        );

        let search_resp: serde_json::Value = client.get(&search_url).send()?.json()?;

        let uid_list = search_resp
            .pointer("/esearchresult/idlist")
            .and_then(|v| v.as_array())
            .cloned()
            .unwrap_or_default();

        if uid_list.is_empty() {
            continue;
        }

        let uid = uid_list[0].as_str().unwrap_or("").to_string();
        if uid.is_empty() {
            continue;
        }

        // Fetch gene summary for this UID
        let summary_url = format!("{}?db=gene&id={}&retmode=json", ESUMMARY_URL, uid);

        let summary_resp: serde_json::Value = client.get(&summary_url).send()?.json()?;

        // Official name is under result.<uid>.description
        // Common name / symbol is under result.<uid>.name
        let name = summary_resp
            .pointer(&format!("/result/{}/description", uid))
            .and_then(|v| v.as_str())
            .or_else(|| {
                summary_resp
                    .pointer(&format!("/result/{}/name", uid))
                    .and_then(|v| v.as_str())
            })
            .unwrap_or("")
            .to_string();

        if !name.is_empty() {
            results.insert(id.clone(), name);
        }

        std::thread::sleep(std::time::Duration::from_millis(350));
    }

    Ok(results)
}

fn urlencoding(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            'A'..='Z' | 'a'..='z' | '0'..='9' | '-' | '_' | '.' | '~' => c.to_string(),
            ' ' => "+".to_string(),
            _ => format!("%{:02X}", c as u32),
        })
        .collect()
}
