//! Phenotypic lookup for Y1000+ species.
//!
//! Right now this module surfaces a single phenotype — growth at 37 °C —
//! from `y1000p_growth_at_37.xlsx` in the bundled `phenotypes` subset.
//! Values are canonical: Y (yes), N (no), W (weak), V (variable), S (slow).
//! Callers merge per-species lookups into a Jaccard-weighted vote in
//! `place.rs` to predict the user's thermotolerance.
//!
//! The `GrowthRates_Yeasts.xlsx` file (continuous per-substrate growth
//! rates across 25 conditions) is deliberately not consumed yet; it'll
//! power per-substrate predictions in a later pass.
//!
//! Citation: Opulente DA et al. (2024). Science 384(6694): eadj4503.

use crate::utils::error::{MycoNoteError, Result};
use crate::y1000plus::metabolism::normalise_species;
use calamine::{open_workbook_auto, Data, Reader};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Growth-at-37 entry for one species.
#[derive(Debug, Clone)]
pub struct GrowthAt37 {
    pub species_pretty: String,
    /// Canonical label: "Y" | "N" | "W" | "V" | "S".
    pub label: String,
    pub source: String,
}

#[derive(Debug, Clone, Default)]
pub struct PhenotypeIndex {
    pub growth_at_37: HashMap<String, GrowthAt37>,
}

impl PhenotypeIndex {
    pub fn lookup_37(&self, kegg_species_stem: &str) -> Option<&GrowthAt37> {
        self.growth_at_37.get(&normalise_species(kegg_species_stem))
    }
    pub fn growth_at_37_count(&self) -> usize {
        self.growth_at_37.len()
    }
}

pub fn load_index(cache_root: &Path) -> Result<PhenotypeIndex> {
    let root = cache_root.join("phenotypes");
    let mut idx = PhenotypeIndex::default();
    if let Some(p) = find_file(&root, "y1000p_growth_at_37") {
        idx.growth_at_37 = parse_growth_at_37(&p)?;
    }
    if idx.growth_at_37.is_empty() {
        return Err(MycoNoteError::UnsupportedFormat(format!(
            "No parsable growth-at-37 table under {}. Install with \
             `setup --y1000plus --include phenotypes`.",
            root.display()
        )));
    }
    Ok(idx)
}

fn parse_growth_at_37(path: &Path) -> Result<HashMap<String, GrowthAt37>> {
    let mut wb = open_workbook_auto(path)
        .map_err(|e| MycoNoteError::ExternalTool(format!("open {}: {e}", path.display())))?;
    let sheet_name =
        wb.sheet_names().first().cloned().ok_or_else(|| {
            MycoNoteError::InvalidFormat(format!("Empty xlsx: {}", path.display()))
        })?;
    let range = wb
        .worksheet_range(&sheet_name)
        .map_err(|e| MycoNoteError::ExternalTool(format!("sheet {sheet_name}: {e}")))?;

    let mut rows = range.rows();
    let header = rows.next().ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!("No header row in {}", path.display()))
    })?;
    let col = |want: &str| {
        header.iter().position(|c| match c {
            Data::String(s) => s.trim().eq_ignore_ascii_case(want),
            _ => false,
        })
    };
    // Robust to column re-orderings.
    let species_col = col("assembly_fullID_updated")
        .or_else(|| col("assembly_fullid_updated"))
        .or_else(|| col("Species"))
        .unwrap_or(0);
    let label_col = col("at_37_C").unwrap_or(1);
    let source_col = col("Temp_Source").unwrap_or(2);

    let mut out: HashMap<String, GrowthAt37> = HashMap::new();
    for row in rows {
        let Some(raw_species) = row.get(species_col).and_then(cell_string) else {
            continue;
        };
        if raw_species.is_empty() {
            continue;
        }
        let label_raw = row.get(label_col).and_then(cell_string).unwrap_or_default();
        let label = canonical_label(&label_raw);
        if label.is_empty() {
            continue;
        }
        let source = row
            .get(source_col)
            .and_then(cell_string)
            .unwrap_or_default();
        let key = normalise_species(&raw_species);
        if key.is_empty() {
            continue;
        }
        out.insert(
            key,
            GrowthAt37 {
                species_pretty: raw_species,
                label,
                source,
            },
        );
    }
    Ok(out)
}

/// Collapse any variant spelling of the growth label to a canonical single
/// letter the rest of the code can match on.
fn canonical_label(raw: &str) -> String {
    let t = raw.trim();
    if t.is_empty() {
        return String::new();
    }
    let first = t.chars().next().unwrap().to_ascii_uppercase();
    match first {
        'Y' | 'N' | 'W' | 'V' | 'S' => first.to_string(),
        _ => String::new(),
    }
}

fn cell_string(c: &Data) -> Option<String> {
    match c {
        Data::String(s) => Some(s.trim().to_string()),
        Data::Float(f) => Some(format!("{f}")),
        Data::Int(i) => Some(format!("{i}")),
        _ => None,
    }
}

fn find_file(root: &Path, stem_hint: &str) -> Option<PathBuf> {
    fn walk(p: &Path, stem_hint: &str, out: &mut Option<PathBuf>) {
        if out.is_some() {
            return;
        }
        let Ok(entries) = std::fs::read_dir(p) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                walk(&path, stem_hint, out);
                if out.is_some() {
                    return;
                }
            } else if let Some(name) = path.file_name().and_then(|s| s.to_str()) {
                if name.to_lowercase().contains(&stem_hint.to_lowercase()) {
                    *out = Some(path);
                    return;
                }
            }
        }
    }
    let mut found = None;
    walk(root, stem_hint, &mut found);
    found
}

/// Human-readable tier from a Y/N/W/V/S label.
pub fn label_description(label: &str) -> &'static str {
    match label {
        "Y" => "grows at 37 °C",
        "N" => "does NOT grow at 37 °C",
        "W" => "weak growth at 37 °C",
        "V" => "variable growth at 37 °C",
        "S" => "slow growth at 37 °C",
        _ => "unknown",
    }
}
