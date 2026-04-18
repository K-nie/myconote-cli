//! Species placement against the Y1000+ bundle.
//!
//! v1 (`mode = functional`): Jaccard similarity of KEGG-KO sets between the
//! user's annotated genome and each of the 1,154 Y1000+ yeasts. Fast, no
//! external tools, uses only the already-installed `kegg` subset.
//!
//! A future `mode = phylogenetic` will layer EPA-ng onto the `phylogeny-place`
//! subset's 1,403 marker MSAs + reference tree.
//!
//! Citation: Opulente et al. 2024, Science 384(6694): eadj4503.

use crate::utils::error::{MycoNoteError, Result};
use crate::y1000plus::benchmark::count_user_kos;
use crate::y1000plus::environment::load_index as load_environment_index;
use crate::y1000plus::manifest::{cache_root, load as load_manifest};
use crate::y1000plus::metabolism::{load_index as load_metabolism_index, Classification};
use crate::y1000plus::phenotypes::{
    label_description as growth_label_desc, load_index as load_phenotype_index,
};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// One per-species placement result.
#[derive(Debug, Clone)]
pub struct Placement {
    pub species: String,
    pub jaccard: f64,
    pub shared: usize,
    pub ref_kos: usize,
}

/// Top-N closest species by KEGG-KO Jaccard similarity.
pub struct PlacementReport {
    pub user_kos: usize,
    pub user_gene_rows: usize,
    pub species_total: usize,
    pub top: Vec<Placement>,
    /// Inferred genetic code for the top hit, if the `codontable` subset is
    /// installed. Present only for the #1 closest species.
    pub top_codon_table: Option<CodonTableAdvice>,
    /// Carbon / nitrogen lifestyle prediction when the `metabolism` subset
    /// is installed. Derived from a Jaccard-weighted majority vote over
    /// those members of `top` that have classification entries.
    pub metabolic_prediction: Option<MetabolicPrediction>,
    /// Growth-at-37 °C thermotolerance prediction when the `phenotypes`
    /// subset is installed. Same Jaccard-weighted vote approach as the
    /// metabolism prediction, over Y/N/W/V/S labels.
    pub thermotolerance: Option<ThermoPrediction>,
    /// Ecological niche / isolation-source prediction when the
    /// `environment` subset is installed. Jaccard-weighted majority over
    /// the friendly labels of the top-N neighbours' isolation ontology.
    pub niche: Option<NichePrediction>,
}

#[derive(Debug, Clone)]
pub struct NichePrediction {
    pub vote: Vote,
    pub neighbours: Vec<NicheNeighbour>,
}

#[derive(Debug, Clone)]
pub struct NicheNeighbour {
    pub species: String,
    pub jaccard: f64,
    pub niche_label: String,
}

#[derive(Debug, Clone)]
pub struct ThermoPrediction {
    pub vote: Vote,
    pub neighbours: Vec<ThermoNeighbour>,
}

#[derive(Debug, Clone)]
pub struct ThermoNeighbour {
    pub species: String,
    pub jaccard: f64,
    pub label: String,
}

#[derive(Debug, Clone)]
pub struct MetabolicPrediction {
    pub carbon_vote: Vote,
    pub nitrogen_vote: Vote,
    /// Per-neighbour rows (top-N that had classifications), in the same
    /// order as `PlacementReport::top`.
    pub neighbours: Vec<NeighbourClassification>,
}

#[derive(Debug, Clone, Default)]
pub struct Vote {
    /// Winning label.
    pub label: String,
    /// Weight share the winner received (0.0–1.0).
    pub confidence: f64,
    /// All label → weight pairs that participated in the vote.
    pub tallies: Vec<(String, f64)>,
}

#[derive(Debug, Clone)]
pub struct NeighbourClassification {
    pub species: String,
    pub jaccard: f64,
    pub carbon_class: String,
    pub nitrogen_class: String,
}

#[derive(Debug, Clone)]
pub struct CodonTableAdvice {
    pub top_species: String,
    /// The raw 64-character code string from codetta (amino acids in standard
    /// codon order TTT, TTC, TTA, TTG, TCT, …).
    pub code: String,
    /// Deviations from NCBI standard table: (codon, standard_aa, inferred_aa).
    pub deviations: Vec<(String, char, char)>,
    /// Convenience: whether CTG is reassigned to Ser (CTG-Ser clade marker).
    pub ctg_is_ser: bool,
}

pub struct PlaceOptions<'a> {
    /// Path to myconote `annotate` output TSV — source of the user's KO set.
    pub annotated_tsv: &'a Path,
    /// How many nearest species to return.
    pub top_n: usize,
}

pub fn place_functional(opts: &PlaceOptions) -> Result<PlacementReport> {
    let root = cache_root();
    if !root.exists() {
        return Err(MycoNoteError::UnsupportedFormat(format!(
            "No Y1000+ cache at {}. Run `myconote-cli setup --y1000plus --include kegg` first.",
            root.display()
        )));
    }
    let manifest = load_manifest(&root).unwrap_or_default();
    if !manifest.installed.contains_key("kegg") {
        return Err(MycoNoteError::UnsupportedFormat(
            "The `kegg` subset isn't installed. \
             Run `myconote-cli setup --y1000plus --include kegg` first."
                .to_string(),
        ));
    }

    let (user_ko_count, user_gene_rows) = count_user_kos(opts.annotated_tsv)?;
    let user_kos = load_user_ko_set(opts.annotated_tsv)?;
    if user_kos.is_empty() {
        return Err(MycoNoteError::InvalidFormat(
            "User's annotated TSV had no KEGG KOs — can't place a profile \
             without any functional annotations."
                .to_string(),
        ));
    }

    let kegg_root = root.join("kegg");
    let mut results: Vec<Placement> = Vec::new();
    for (species_file, species_name) in walk_species_files(&kegg_root)? {
        let ref_kos = load_species_ko_set(&species_file)?;
        if ref_kos.is_empty() {
            continue;
        }
        let intersection = user_kos.intersection(&ref_kos).count();
        let union = user_kos.union(&ref_kos).count();
        let jaccard = if union == 0 {
            0.0
        } else {
            intersection as f64 / union as f64
        };
        results.push(Placement {
            species: species_name,
            jaccard,
            shared: intersection,
            ref_kos: ref_kos.len(),
        });
    }

    let species_total = results.len();
    results.sort_by(|a, b| {
        b.jaccard
            .partial_cmp(&a.jaccard)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    results.truncate(opts.top_n);

    // If the codontable subset is installed, derive the inferred genetic
    // code for the #1 hit so we can warn about CTG-Ser clades at the end.
    let top_codon_table = if manifest.installed.contains_key("codontable") {
        results
            .first()
            .and_then(|hit| load_codon_table(&root, &hit.species).ok())
    } else {
        None
    };

    // If the metabolism subset is installed, derive a carbon/nitrogen
    // lifestyle prediction from a Jaccard-weighted majority vote over the
    // top-N neighbours that we can classify.
    let metabolic_prediction = if manifest.installed.contains_key("metabolism") {
        load_metabolism_index(&root)
            .ok()
            .map(|idx| predict_metabolism(&results, &idx))
            .and_then(|p| {
                if p.neighbours.is_empty() {
                    None
                } else {
                    Some(p)
                }
            })
    } else {
        None
    };

    // If the phenotypes subset is installed, predict thermotolerance from
    // the nearest neighbours' growth-at-37 labels.
    let thermotolerance = if manifest.installed.contains_key("phenotypes") {
        load_phenotype_index(&root)
            .ok()
            .map(|idx| predict_thermotolerance(&results, &idx))
            .and_then(|p| {
                if p.neighbours.is_empty() {
                    None
                } else {
                    Some(p)
                }
            })
    } else {
        None
    };

    // If the environment subset is installed, predict ecological niche /
    // isolation source from the nearest neighbours' ontology labels.
    let niche = if manifest.installed.contains_key("environment") {
        load_environment_index(&root)
            .ok()
            .map(|idx| predict_niche(&results, &idx))
            .and_then(|p| {
                if p.neighbours.is_empty() {
                    None
                } else {
                    Some(p)
                }
            })
    } else {
        None
    };

    Ok(PlacementReport {
        user_kos: user_ko_count,
        user_gene_rows,
        species_total,
        top: results,
        top_codon_table,
        metabolic_prediction,
        thermotolerance,
        niche,
    })
}

fn predict_niche(
    top: &[Placement],
    idx: &crate::y1000plus::environment::EnvironmentIndex,
) -> NichePrediction {
    let mut neighbours = Vec::new();
    let mut weights: HashMap<String, f64> = HashMap::new();
    for hit in top {
        let Some(n) = idx.lookup(&hit.species) else {
            continue;
        };
        let weight = hit.jaccard.max(1e-4);
        *weights.entry(n.niche_label.clone()).or_insert(0.0) += weight;
        neighbours.push(NicheNeighbour {
            species: n.species_pretty.clone(),
            jaccard: hit.jaccard,
            niche_label: n.niche_label.clone(),
        });
    }
    NichePrediction {
        vote: vote_from(weights),
        neighbours,
    }
}

/// Jaccard-weighted vote over nearest-neighbours' Y/N/W/V/S growth-at-37
/// labels. Unknown species are silently skipped.
fn predict_thermotolerance(
    top: &[Placement],
    idx: &crate::y1000plus::phenotypes::PhenotypeIndex,
) -> ThermoPrediction {
    let mut neighbours = Vec::new();
    let mut weights: HashMap<String, f64> = HashMap::new();
    for hit in top {
        let Some(g) = idx.lookup_37(&hit.species) else {
            continue;
        };
        let weight = hit.jaccard.max(1e-4);
        *weights.entry(g.label.clone()).or_insert(0.0) += weight;
        neighbours.push(ThermoNeighbour {
            species: g.species_pretty.clone(),
            jaccard: hit.jaccard,
            label: g.label.clone(),
        });
    }
    ThermoPrediction {
        vote: vote_from(weights),
        neighbours,
    }
}

/// Friendly name for a Y/N/W/V/S label so callers can surface it in
/// `place` output without knowing the legend by heart.
pub fn thermo_label_description(label: &str) -> &'static str {
    growth_label_desc(label)
}

/// Jaccard-weighted majority vote over the carbon/nitrogen classes of the
/// top-N neighbours that appear in the metabolism index. Unknown species
/// are silently skipped — the prediction degrades gracefully as overlap
/// with the index shrinks.
fn predict_metabolism(
    top: &[Placement],
    idx: &crate::y1000plus::metabolism::MetabolismIndex,
) -> MetabolicPrediction {
    let mut neighbours = Vec::new();
    let mut carbon_weights: HashMap<String, f64> = HashMap::new();
    let mut nitrogen_weights: HashMap<String, f64> = HashMap::new();

    for hit in top {
        let Some(cls) = idx.lookup(&hit.species) else {
            continue;
        };
        // Use max(epsilon, jaccard) so a 0-Jaccard neighbour still
        // contributes trivially (the earlier sort already put the best
        // ones first, and dropping them entirely would silently lose
        // neighbours the user expects to see in the report).
        let weight = hit.jaccard.max(1e-4);
        *carbon_weights
            .entry(classify_label(&cls.carbon_class))
            .or_insert(0.0) += weight;
        *nitrogen_weights
            .entry(classify_label(&cls.nitrogen_class))
            .or_insert(0.0) += weight;
        neighbours.push(NeighbourClassification {
            species: cls.species_pretty.clone(),
            jaccard: hit.jaccard,
            carbon_class: cls.carbon_class.clone(),
            nitrogen_class: cls.nitrogen_class.clone(),
        });
    }

    MetabolicPrediction {
        carbon_vote: vote_from(carbon_weights),
        nitrogen_vote: vote_from(nitrogen_weights),
        neighbours,
    }
}

fn classify_label(raw: &str) -> String {
    let t = raw.trim();
    if t.is_empty() {
        "Unknown".to_string()
    } else {
        t.to_string()
    }
}

fn vote_from(mut weights: HashMap<String, f64>) -> Vote {
    let total: f64 = weights.values().sum();
    if total == 0.0 {
        return Vote::default();
    }
    let mut tallies: Vec<(String, f64)> = weights.drain().collect();
    tallies.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let (label, top_weight) = tallies[0].clone();
    Vote {
        label,
        confidence: top_weight / total,
        tallies,
    }
}

/// NCBI standard genetic code — 64 single-letter amino acids in the codon
/// order TTT, TTC, TTA, TTG, TCT, … (matches the codetta `.code` layout).
const STANDARD_CODE: &[u8; 64] =
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const CODON_INDEX_CTG: usize = 19;

fn codon_at(idx: usize) -> String {
    const BASES: [char; 4] = ['T', 'C', 'A', 'G'];
    let b1 = BASES[(idx / 16) % 4];
    let b2 = BASES[(idx / 4) % 4];
    let b3 = BASES[idx % 4];
    format!("{b1}{b2}{b3}")
}

/// Find and parse the `.code` file for the given species under the
/// `codetta/` cache. Codetta may record unresolved codons as `?`; those are
/// ignored rather than flagged as deviations.
fn load_codon_table(cache_root: &Path, species_stem: &str) -> Result<CodonTableAdvice> {
    let codetta_root = cache_root.join("codetta");
    let mut found: Option<PathBuf> = None;
    fn walk(p: &Path, stem: &str, found: &mut Option<PathBuf>) -> std::io::Result<()> {
        if found.is_some() {
            return Ok(());
        }
        for entry in std::fs::read_dir(p)?.flatten() {
            let path = entry.path();
            if path.is_dir() {
                walk(&path, stem, found)?;
                if found.is_some() {
                    return Ok(());
                }
            } else if let Some(name) = path.file_name().and_then(|s| s.to_str()) {
                if name.ends_with(".code") && name.contains(stem) {
                    *found = Some(path);
                    return Ok(());
                }
            }
        }
        Ok(())
    }
    walk(&codetta_root, species_stem, &mut found).map_err(MycoNoteError::Io)?;
    let path = found.ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!(
            "No codetta .code file for species '{species_stem}'"
        ))
    })?;
    let raw = std::fs::read_to_string(&path).map_err(MycoNoteError::Io)?;
    let code: String = raw.trim().chars().take(64).collect();
    if code.len() < 64 {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Truncated codetta code in {}: only {} chars",
            path.display(),
            code.len()
        )));
    }
    let code_bytes = code.as_bytes();

    let mut deviations = Vec::new();
    for (i, (inferred, standard)) in code_bytes.iter().zip(STANDARD_CODE.iter()).enumerate() {
        if *inferred == b'?' {
            continue;
        }
        if inferred != standard {
            deviations.push((codon_at(i), *standard as char, *inferred as char));
        }
    }
    let ctg_is_ser = code_bytes[CODON_INDEX_CTG] == b'S';

    Ok(CodonTableAdvice {
        top_species: species_stem.to_string(),
        code,
        deviations,
        ctg_is_ser,
    })
}

fn load_user_ko_set(path: &Path) -> Result<HashSet<String>> {
    let f = File::open(path).map_err(MycoNoteError::Io)?;
    let mut lines = BufReader::new(f).lines();
    let header = lines
        .next()
        .transpose()
        .map_err(MycoNoteError::Io)?
        .ok_or_else(|| {
            MycoNoteError::InvalidFormat(format!("Empty annotated TSV: {}", path.display()))
        })?;
    let headers: Vec<&str> = header.split('\t').collect();
    let kegg_col = headers
        .iter()
        .position(|h| {
            let h = h.trim().to_lowercase();
            h == "kegg_pathways" || h == "kegg_ko" || h == "kegg"
        })
        .ok_or_else(|| {
            MycoNoteError::InvalidFormat(format!(
                "No kegg_pathways/kegg_ko column in {}. Is this a myconote `annotate` output?",
                path.display()
            ))
        })?;

    let mut kos: HashSet<String> = HashSet::new();
    for line in lines.map_while(|l| l.ok()) {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= kegg_col {
            continue;
        }
        for tok in cols[kegg_col].split(',') {
            let tok = tok.trim();
            if looks_like_ko(tok) {
                kos.insert(tok.to_string());
            }
        }
    }
    Ok(kos)
}

fn load_species_ko_set(path: &Path) -> Result<HashSet<String>> {
    let mut kos: HashSet<String> = HashSet::new();
    let f = File::open(path).map_err(MycoNoteError::Io)?;
    for line in BufReader::new(f).lines().map_while(|l| l.ok()) {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            continue;
        }
        let ko = cols[1].trim();
        if looks_like_ko(ko) {
            kos.insert(ko.to_string());
        }
    }
    Ok(kos)
}

fn looks_like_ko(tok: &str) -> bool {
    tok.starts_with('K')
        && (2..=8).contains(&tok.len())
        && tok.chars().skip(1).all(|c| c.is_ascii_digit())
}

/// Yield (path, pretty_species_name) for every per-species KEGG TSV.
fn walk_species_files(kegg_root: &Path) -> Result<Vec<(PathBuf, String)>> {
    let mut out = Vec::new();
    fn walk(p: &Path, out: &mut Vec<(PathBuf, String)>) -> std::io::Result<()> {
        for entry in std::fs::read_dir(p)?.flatten() {
            let path = entry.path();
            if path.is_dir() {
                walk(&path, out)?;
                continue;
            }
            let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
                continue;
            };
            if name.eq_ignore_ascii_case("README.txt")
                || name.eq_ignore_ascii_case("README.md")
                || name.starts_with('.')
            {
                continue;
            }
            if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
                if ext == "tsv" || ext == "txt" {
                    // Strip the last dot-suffix only — preserves species names
                    // that contain underscores or dots naturally.
                    let stem = path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or(name)
                        .to_string();
                    out.push((path, stem));
                }
            }
        }
        Ok(())
    }
    walk(kegg_root, &mut out).map_err(MycoNoteError::Io)?;
    Ok(out)
}
