/// UniProt REST API backend
///
/// Searches UniProtKB for gene names using the gene ID as a query term.
/// Works well for model organism genes and widely-studied fungal species.
///
/// API docs: https://rest.uniprot.org/docs/

use std::collections::HashMap;

const UNIPROT_SEARCH: &str = "https://rest.uniprot.org/uniprotkb/search";
const BATCH_SIZE: usize = 20; // UniProt prefers smaller batches

/// Resolve gene IDs to gene names via UniProt.
pub fn fetch_gene_names(ids: &[String]) -> HashMap<String, String> {
    let mut result = HashMap::new();

    for chunk in ids.chunks(BATCH_SIZE) {
        match fetch_chunk(chunk) {
            Ok(partial) => result.extend(partial),
            Err(e) => eprintln!("    UniProt fetch error: {}", e),
        }
        std::thread::sleep(std::time::Duration::from_millis(500));
    }

    result
}

fn fetch_chunk(ids: &[String]) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let client = reqwest::blocking::Client::builder()
        .user_agent("myconote-cli/0.1")
        .timeout(std::time::Duration::from_secs(30))
        .build()?;

    let mut results = HashMap::new();

    for id in ids {
        // Query UniProt for entries where the gene name matches the ID
        let url = format!(
            "{}?query=gene_exact:{}+AND+reviewed:true&fields=gene_names,protein_name&format=tsv&size=1",
            UNIPROT_SEARCH,
            urlencoding(id)
        );

        let text = client.get(&url).send()?.text()?;

        // TSV response: header line + data lines
        // Columns: Gene Names, Protein names
        for line in text.lines().skip(1) {
            let cols: Vec<&str> = line.splitn(3, '\t').collect();
            if cols.len() >= 2 {
                // Gene Names column may have multiple space-separated names;
                // take the first (primary) one
                let gene_name = cols[0].split_whitespace().next().unwrap_or("").to_string();
                if !gene_name.is_empty() {
                    results.insert(id.clone(), gene_name);
                    break;
                }
                // Fall back to protein name if no gene name found
                if !cols[1].is_empty() {
                    // Strip isoform/variant suffixes like " (Fragment)"
                    let prot_name = cols[1]
                        .split('(')
                        .next()
                        .unwrap_or(cols[1])
                        .trim()
                        .to_string();
                    if !prot_name.is_empty() {
                        results.insert(id.clone(), prot_name);
                    }
                    break;
                }
            }
        }

        std::thread::sleep(std::time::Duration::from_millis(300));
    }

    Ok(results)
}

fn urlencoding(s: &str) -> String {
    s.chars().map(|c| match c {
        'A'..='Z' | 'a'..='z' | '0'..='9' | '-' | '_' | '.' | '~' => c.to_string(),
        ' ' => "+".to_string(),
        _ => format!("%{:02X}", c as u32),
    }).collect()
}
