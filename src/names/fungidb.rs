/// FungiDB / VEuPathDB REST API backend
///
/// VEuPathDB (which hosts FungiDB, ToxoDB, PlasmoDB, etc.) provides a
/// REST API for querying gene information by stable ID.
///
/// This backend is most useful for species whose IDs have the pattern:
///   CNAG_00001  (Cryptococcus neoformans — FungiDB)
///   AFUB_000010 (Aspergillus fumigatus — FungiDB)
///   Fvpg_G000010 (Fusarium — FungiDB)
///   AT1G01010   (Arabidopsis — PhytozomeDB / EnsemblPlants)
///   ENSMUSG...  (Mouse — Ensembl)
///
/// The VEuPathDB endpoint also accepts PromBase-style IDs for supported
/// yeast species.
///
/// API: https://veupathdb.org/veupathdb/service/record-types/gene/records

use std::collections::HashMap;

const VEUPATHDB_RECORDS: &str =
    "https://veupathdb.org/veupathdb/service/record-types/gene/records";

const ENSEMBL_LOOKUP: &str =
    "https://rest.ensembl.org/lookup/id";

const BATCH_SIZE: usize = 20;

/// Resolve gene IDs via FungiDB / VEuPathDB, then fall back to Ensembl
/// for animal/plant IDs that look like Ensembl stable IDs.
pub fn fetch_gene_names(ids: &[String]) -> HashMap<String, String> {
    let mut result = HashMap::new();

    // Split IDs: Ensembl-style vs VEuPathDB-style
    let (ensembl_ids, veupathdb_ids): (Vec<_>, Vec<_>) = ids.iter()
        .partition(|id| is_ensembl_id(id));

    // VEuPathDB batch
    for chunk in veupathdb_ids.chunks(BATCH_SIZE) {
        match fetch_veupathdb(chunk) {
            Ok(partial) => result.extend(partial),
            Err(e) => eprintln!("    FungiDB fetch error: {}", e),
        }
        std::thread::sleep(std::time::Duration::from_millis(500));
    }

    // Ensembl batch (animals/plants)
    for chunk in ensembl_ids.chunks(BATCH_SIZE) {
        match fetch_ensembl(chunk) {
            Ok(partial) => result.extend(partial),
            Err(e) => eprintln!("    Ensembl fetch error: {}", e),
        }
        std::thread::sleep(std::time::Duration::from_millis(500));
    }

    result
}

fn is_ensembl_id(id: &str) -> bool {
    // Ensembl stable IDs start with ENS (e.g. ENSG, ENSMUSG, ENSDARG)
    id.starts_with("ENS")
}

fn fetch_veupathdb(
    ids: &[&String],
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    use serde_json::json;

    let client = reqwest::blocking::Client::builder()
        .user_agent("myconote-cli/0.1")
        .timeout(std::time::Duration::from_secs(30))
        .build()?;

    let mut results = HashMap::new();

    for id in ids {
        // POST request with the gene's primary key
        let body = json!({
            "primaryKey": [{"name": "source_id", "value": id.as_str()}],
            "tables": []
        });

        let resp = client
            .post(VEUPATHDB_RECORDS)
            .header("Content-Type", "application/json")
            .json(&body)
            .send();

        let Ok(resp) = resp else { continue };
        let Ok(json): Result<serde_json::Value, _> = resp.json() else { continue };

        // Gene name is in attributes.gene_name or attributes.display_name
        let name = json.pointer("/attributes/gene_name/value")
            .or_else(|| json.pointer("/attributes/display_name/value"))
            .or_else(|| json.pointer("/displayName"))
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty() && *s != id.as_str())
            .map(str::to_string);

        if let Some(n) = name {
            results.insert(id.to_string(), n);
        }

        std::thread::sleep(std::time::Duration::from_millis(300));
    }

    Ok(results)
}

fn fetch_ensembl(
    ids: &[&String],
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    use serde_json::json;

    let client = reqwest::blocking::Client::builder()
        .user_agent("myconote-cli/0.1")
        .timeout(std::time::Duration::from_secs(30))
        .build()?;

    // Ensembl POST lookup for multiple IDs at once
    let id_list: Vec<&str> = ids.iter().map(|s| s.as_str()).collect();
    let body = json!({ "ids": id_list });

    let resp = client
        .post(format!("{}", ENSEMBL_LOOKUP))
        .header("Content-Type", "application/json")
        .header("Accept", "application/json")
        .json(&body)
        .send()?;

    let json: serde_json::Value = resp.json()?;
    let mut results = HashMap::new();

    for id in ids {
        if let Some(entry) = json.get(id.as_str()) {
            // Ensembl returns display_name and description
            let name = entry.get("display_name")
                .and_then(|v| v.as_str())
                .filter(|s| !s.is_empty())
                .or_else(|| entry.get("description")
                    .and_then(|v| v.as_str())
                    .map(|s| s.split('[').next().unwrap_or(s).trim())
                    .filter(|s| !s.is_empty()))
                .map(str::to_string);

            if let Some(n) = name {
                results.insert(id.to_string(), n);
            }
        }
    }

    Ok(results)
}
