//! Carbon / nitrogen metabolic classification lookup for Y1000+ species.
//!
//! Parses the `CarbonNitrogen_Classifications` xlsx bundled with the
//! `metabolism` subset into an in-memory index keyed by a normalised species
//! binomial. Callers look up an entry per top-N Y1000+ hit from `place`
//! and compute a Jaccard-weighted majority vote to predict the user's own
//! carbon / nitrogen lifestyle.
//!
//! Per Opulente et al. 2024, classifications are ternary:
//!   * **Specialist** — grows on very few substrates of that element
//!   * **Standard**   — typical breadth
//!   * **Generalist** — grows on many substrates
//!
//! Citation: Opulente DA et al. (2024). Science 384(6694): eadj4503.

use crate::utils::error::{MycoNoteError, Result};
use calamine::{open_workbook_auto, Data, Reader};
use std::collections::HashMap;
use std::path::Path;

/// One species' metabolic classification row from the bundled xlsx.
#[derive(Debug, Clone)]
pub struct Classification {
    pub species_pretty: String,
    pub carbon_breadth: Option<f64>,
    pub nitrogen_breadth: Option<f64>,
    pub carbon_class: String,
    pub nitrogen_class: String,
}

/// In-memory lookup table built once from the xlsx. Key is the normalised
/// binomial (lowercased "genus species") — callers normalise their own
/// species names before querying.
#[derive(Debug, Clone, Default)]
pub struct MetabolismIndex {
    pub by_species: HashMap<String, Classification>,
}

impl MetabolismIndex {
    pub fn lookup(&self, kegg_species_stem: &str) -> Option<&Classification> {
        self.by_species.get(&normalise_species(kegg_species_stem))
    }

    pub fn len(&self) -> usize {
        self.by_species.len()
    }
}

/// Locate and parse the classifications xlsx. Searches the metabolism/
/// subset tree for `Yeast Breadth Classifications with Growth Data*.xlsx`.
pub fn load_index(cache_root: &Path) -> Result<MetabolismIndex> {
    let metabolism_root = cache_root.join("metabolism");
    let xlsx = find_xlsx(&metabolism_root).ok_or_else(|| {
        MycoNoteError::UnsupportedFormat(format!(
            "No classifications xlsx under {}. Install with `setup --y1000plus --include metabolism`.",
            metabolism_root.display()
        ))
    })?;

    let mut wb = open_workbook_auto(&xlsx)
        .map_err(|e| MycoNoteError::ExternalTool(format!("open {}: {e}", xlsx.display())))?;
    let sheet_name =
        wb.sheet_names().first().cloned().ok_or_else(|| {
            MycoNoteError::InvalidFormat(format!("Empty xlsx: {}", xlsx.display()))
        })?;
    let range = wb
        .worksheet_range(&sheet_name)
        .map_err(|e| MycoNoteError::ExternalTool(format!("sheet {}: {e}", sheet_name)))?;

    // Row 0 is the header — find the column indices we care about by name so
    // we're resilient to minor column re-orderings.
    let mut rows = range.rows();
    let header = rows.next().ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!("No header row in {}", xlsx.display()))
    })?;
    let col = |want: &str| {
        header.iter().position(|c| match c {
            Data::String(s) => s.trim().eq_ignore_ascii_case(want),
            _ => false,
        })
    };
    let species_col = col("Species").ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!("No 'Species' column in {}", xlsx.display()))
    })?;
    let c_breadth_col = col("Carbon Breadth");
    let n_breadth_col = col("Nitrogen Breadth");
    let c_class_col = col("Carbon Class");
    let n_class_col = col("Nitrogen Class");

    let mut index = MetabolismIndex::default();
    for row in rows {
        let Some(species_pretty) = row.get(species_col).and_then(cell_string) else {
            continue;
        };
        if species_pretty.is_empty() {
            continue;
        }
        let key = normalise_species(&species_pretty);
        if key.is_empty() {
            continue;
        }
        let carbon_breadth = c_breadth_col.and_then(|c| row.get(c).and_then(cell_float));
        let nitrogen_breadth = n_breadth_col.and_then(|c| row.get(c).and_then(cell_float));
        let carbon_class = c_class_col
            .and_then(|c| row.get(c).and_then(cell_string))
            .unwrap_or_else(|| "Unknown".to_string());
        let nitrogen_class = n_class_col
            .and_then(|c| row.get(c).and_then(cell_string))
            .unwrap_or_else(|| "Unknown".to_string());
        index.by_species.insert(
            key,
            Classification {
                species_pretty,
                carbon_breadth,
                nitrogen_breadth,
                carbon_class,
                nitrogen_class,
            },
        );
    }
    Ok(index)
}

fn find_xlsx(root: &Path) -> Option<std::path::PathBuf> {
    fn walk(p: &Path, out: &mut Option<std::path::PathBuf>) {
        if out.is_some() {
            return;
        }
        let Ok(entries) = std::fs::read_dir(p) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                walk(&path, out);
                if out.is_some() {
                    return;
                }
            } else if let Some(name) = path.file_name().and_then(|s| s.to_str()) {
                // Prefer the classifications file specifically over the big
                // Y1000_KEGG_Annotations.xlsx in the same subset.
                if name.to_lowercase().contains("breadth classifications") {
                    *out = Some(path);
                    return;
                }
            }
        }
    }
    let mut found = None;
    walk(root, &mut found);
    found
}

fn cell_string(c: &Data) -> Option<String> {
    match c {
        Data::String(s) => Some(s.trim().to_string()),
        Data::Float(f) => Some(format!("{f}")),
        Data::Int(i) => Some(format!("{i}")),
        _ => None,
    }
}

fn cell_float(c: &Data) -> Option<f64> {
    match c {
        Data::Float(f) => Some(*f),
        Data::Int(i) => Some(*i as f64),
        Data::String(s) => s.trim().parse::<f64>().ok(),
        _ => None,
    }
}

/// Normalise a species name to a lookup key: lowercase, drop Y1000+
/// accession prefixes (yHMPuXXXXXXX_), strip strain/date/assembler
/// suffixes (_170307, .haplomerger2, _MASURCA, etc.) and keep only the
/// first two alphabetic word tokens (genus + species).
///
/// Examples:
///   `Candida tropicalis`                            → "candida tropicalis"
///   `candida_tropicalis`                            → "candida tropicalis"
///   `yHMPu5000035693_hanseniaspora_nectarophila_160613.haplomerger2`
///                                                   → "hanseniaspora nectarophila"
pub fn normalise_species(raw: &str) -> String {
    // Expand separators to spaces, lowercase, split into tokens.
    let cleaned: String = raw
        .chars()
        .map(|c| match c {
            '_' | '-' | '.' | '/' => ' ',
            c => c.to_ascii_lowercase(),
        })
        .collect();
    let tokens: Vec<&str> = cleaned
        .split_whitespace()
        .filter(|t| !is_accession(t))
        .collect();
    // Find the first two consecutive purely-alphabetic tokens ≥ 3 chars
    // and return them as "genus species". Falls back to the first two
    // alphabetic tokens if no 3-char pair is found.
    for w in tokens.windows(2) {
        let [a, b] = [w[0], w[1]];
        if is_taxon_token(a) && is_taxon_token(b) {
            return format!("{a} {b}");
        }
    }
    let alpha: Vec<&str> = tokens.into_iter().filter(|t| is_taxon_token(t)).collect();
    if alpha.len() >= 2 {
        format!("{} {}", alpha[0], alpha[1])
    } else if !alpha.is_empty() {
        alpha[0].to_string()
    } else {
        String::new()
    }
}

fn is_taxon_token(t: &str) -> bool {
    t.len() >= 3 && t.chars().all(|c| c.is_ascii_alphabetic())
}

/// True if `t` looks like an assembly/strain accession (mostly digits, or
/// mixed with a short alpha prefix like `yhmpu5000035693`).
fn is_accession(t: &str) -> bool {
    if t.chars().all(|c| c.is_ascii_digit()) {
        return true;
    }
    let digits = t.chars().filter(|c| c.is_ascii_digit()).count();
    digits >= 3 && digits * 2 >= t.len()
}
