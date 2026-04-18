//! Reference distributions derived from the Y1000+ bundle, used by
//! `stats --benchmark y1000plus`.
//!
//! The two distributions we derive right now:
//!
//! * **BUSCO completeness** (`busco/Y1000p_BUCO_fulltable/*.full_table.tsv`)
//!   — per-species % of the `saccharomycetes_odb10` BUSCOs recovered as
//!   Complete + Duplicated.
//! * **KEGG-KO count** (`kegg/y1000plus_annotations_pep_kegg/*.txt`) —
//!   number of distinct KOs assigned to each species' proteome.
//!
//! We deliberately skip per-species gene counts here (the raw annotations
//! aren't in the starter bundle). They'll join this module the day the
//! user installs the `annotations` subset.
//!
//! Citation: Opulente et al. 2024, Science 384(6694): eadj4503.

use crate::utils::error::{MycoNoteError, Result};
use crate::y1000plus::manifest::{cache_root, load as load_manifest};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Summary statistics over a numeric distribution.
#[derive(Debug, Clone)]
pub struct DistStats {
    pub n: usize,
    pub mean: f64,
    pub median: f64,
    pub p25: f64,
    pub p75: f64,
    pub min: f64,
    pub max: f64,
    /// Every value, sorted ascending — used to compute percentile ranks for
    /// the user's own genome in `percentile_of`.
    pub sorted: Vec<f64>,
}

impl DistStats {
    fn from(mut values: Vec<f64>) -> Option<Self> {
        if values.is_empty() {
            return None;
        }
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = values.len();
        let mean = values.iter().sum::<f64>() / n as f64;
        let pct = |p: f64| {
            let idx = ((n as f64 - 1.0) * p).round() as usize;
            values[idx.min(n - 1)]
        };
        Some(DistStats {
            n,
            mean,
            median: pct(0.50),
            p25: pct(0.25),
            p75: pct(0.75),
            min: values[0],
            max: values[n - 1],
            sorted: values,
        })
    }

    /// Percentile rank of `x` in this distribution, 0–100. Uses the
    /// "fraction of reference values ≤ x" convention — matches how most
    /// users interpret "my value is at the Nth percentile".
    pub fn percentile_of(&self, x: f64) -> f64 {
        if self.sorted.is_empty() {
            return f64::NAN;
        }
        let below = self.sorted.iter().filter(|v| **v <= x).count();
        100.0 * below as f64 / self.sorted.len() as f64
    }
}

/// Full reference snapshot loaded from an installed Y1000+ cache.
#[derive(Debug, Clone, Default)]
pub struct Reference {
    pub busco_completeness: Option<DistStats>,
    pub kegg_ko_count: Option<DistStats>,
    pub trna_count: Option<DistStats>,
    pub source_root: PathBuf,
    pub species_count: usize,
}

/// Load whatever reference distributions are present in the cache. Missing
/// subsets are reported as `None` instead of errors — benchmarking degrades
/// gracefully based on what the user has installed.
pub fn load_reference() -> Result<Reference> {
    let root = cache_root();
    if !root.exists() {
        return Err(MycoNoteError::UnsupportedFormat(format!(
            "No Y1000+ cache at {}. Run `myconote-cli setup --y1000plus --preset starter` first.",
            root.display()
        )));
    }
    let manifest = load_manifest(&root).unwrap_or_default();

    let busco = if manifest.installed.contains_key("busco") {
        load_busco_distribution(&root.join("busco")).ok()
    } else {
        None
    };

    let kegg = if manifest.installed.contains_key("kegg") {
        load_kegg_distribution(&root.join("kegg")).ok()
    } else {
        None
    };

    let trna = if manifest.installed.contains_key("trna") {
        load_trna_distribution(&root.join("trna")).ok()
    } else {
        None
    };

    let species_count = [
        busco.as_ref().map(|d| d.n),
        kegg.as_ref().map(|d| d.n),
        trna.as_ref().map(|d| d.n),
    ]
    .into_iter()
    .flatten()
    .max()
    .unwrap_or(0);

    Ok(Reference {
        busco_completeness: busco,
        kegg_ko_count: kegg,
        trna_count: trna,
        source_root: root,
        species_count,
    })
}

/// Walk `<root>/busco/**/full_table.tsv` files; for each species compute
/// completeness = (unique Complete + unique Duplicated) / total expected.
/// The per-file header tells us the BUSCO lineage size.
fn load_busco_distribution(busco_root: &Path) -> Result<DistStats> {
    let mut values = Vec::new();
    for entry in walk_tsv_files(busco_root)? {
        if let Ok(pct) = per_species_busco_completeness(&entry) {
            values.push(pct);
        }
    }
    DistStats::from(values).ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!(
            "No parseable BUSCO full_tables under {}",
            busco_root.display()
        ))
    })
}

fn per_species_busco_completeness(path: &Path) -> Result<f64> {
    let mut total_buscos: Option<usize> = None;
    let mut complete: HashSet<String> = HashSet::new();
    let mut duplicated: HashSet<String> = HashSet::new();

    let f = File::open(path).map_err(MycoNoteError::Io)?;
    for line in BufReader::new(f).lines().map_while(|l| l.ok()) {
        if let Some(stripped) = line.strip_prefix('#') {
            // Header example:
            //   # The lineage dataset is: saccharomycetes_odb10 …
            //     number of BUSCOs: 2137
            if let Some(idx) = stripped.find("number of BUSCOs:") {
                let tail = &stripped[idx + "number of BUSCOs:".len()..];
                // Value may be followed by more text — grab the first int.
                let digits: String = tail
                    .chars()
                    .skip_while(|c| !c.is_ascii_digit())
                    .take_while(|c| c.is_ascii_digit())
                    .collect();
                if let Ok(n) = digits.parse::<usize>() {
                    total_buscos = Some(n);
                }
            }
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            continue;
        }
        let busco_id = cols[0].to_string();
        match cols[1] {
            "Complete" => {
                complete.insert(busco_id);
            }
            "Duplicated" => {
                duplicated.insert(busco_id);
            }
            _ => {}
        }
    }
    let total = total_buscos.unwrap_or(0);
    if total == 0 {
        return Err(MycoNoteError::InvalidFormat(format!(
            "BUSCO total not recorded in header of {}",
            path.display()
        )));
    }
    Ok(100.0 * (complete.len() + duplicated.len()) as f64 / total as f64)
}

/// Walk `<root>/trna/**/*.tRNA.gff`; for each species count the number of
/// features with `type == "tRNA"`. Each GFF file represents one species
/// (the `y1000p_tRNA_scan` subset ships one tRNAscan run per genome).
fn load_trna_distribution(trna_root: &Path) -> Result<DistStats> {
    let mut values = Vec::new();
    fn walk(p: &Path, out: &mut Vec<f64>) {
        let Ok(entries) = std::fs::read_dir(p) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                walk(&path, out);
                continue;
            }
            let Some(name) = path.file_name().and_then(|s| s.to_str()) else {
                continue;
            };
            if !name.ends_with(".tRNA.gff") {
                continue;
            }
            if let Ok(f) = File::open(&path) {
                let mut count = 0u64;
                for line in BufReader::new(f).lines().map_while(|l| l.ok()) {
                    if line.starts_with('#') || line.is_empty() {
                        continue;
                    }
                    let cols: Vec<&str> = line.split('\t').collect();
                    if cols.len() >= 3 && cols[2] == "tRNA" {
                        count += 1;
                    }
                }
                if count > 0 {
                    out.push(count as f64);
                }
            }
        }
    }
    walk(trna_root, &mut values);
    DistStats::from(values).ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!(
            "No parseable tRNAscan GFFs under {}",
            trna_root.display()
        ))
    })
}

/// Count tRNA features in a user's GFF3. Used to place the user on the
/// Y1000+ tRNA-count distribution. Returns the raw count; 0 if the GFF
/// has no tRNA features (which is itself meaningful — probably annotate
/// wasn't run with tRNAscan enabled).
pub fn count_user_trnas(gff_path: &Path) -> Result<usize> {
    let f = File::open(gff_path).map_err(MycoNoteError::Io)?;
    let mut count = 0usize;
    for line in BufReader::new(f).lines().map_while(|l| l.ok()) {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() >= 3 && cols[2] == "tRNA" {
            count += 1;
        }
    }
    Ok(count)
}

/// Per-species KEGG-KO files: two-column TSVs `gene_id<TAB>KO`. We count
/// distinct KOs per species — the canonical "functional repertoire size".
fn load_kegg_distribution(kegg_root: &Path) -> Result<DistStats> {
    let mut values = Vec::new();
    for entry in walk_tsv_files(kegg_root)? {
        if let Ok(count) = per_species_ko_count(&entry) {
            values.push(count as f64);
        }
    }
    DistStats::from(values).ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!(
            "No parseable KEGG TSVs under {}",
            kegg_root.display()
        ))
    })
}

fn per_species_ko_count(path: &Path) -> Result<usize> {
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
        if ko.starts_with('K') && ko.len() > 1 {
            kos.insert(ko.to_string());
        }
    }
    Ok(kos.len())
}

/// Count distinct KEGG-KO IDs in a myconote-annotated TSV. Expects the
/// EggNog-derived table written by `write_eggnog_table` (columns include
/// `kegg_pathways`, which despite its name carries KEGG_ko entries straight
/// from EggNog-mapper's column 12, i.e. `K00001`-style KOs).
///
/// Returns (distinct_ko_count, gene_rows_scanned) so callers can report
/// how many genes carried any KEGG annotation at all.
pub fn count_user_kos(annotated_tsv: &Path) -> Result<(usize, usize)> {
    let f = File::open(annotated_tsv).map_err(MycoNoteError::Io)?;
    let mut lines = BufReader::new(f).lines();

    let header_line = lines
        .next()
        .transpose()
        .map_err(MycoNoteError::Io)?
        .ok_or_else(|| {
            MycoNoteError::InvalidFormat(format!(
                "Empty annotated TSV: {}",
                annotated_tsv.display()
            ))
        })?;
    let headers: Vec<&str> = header_line.split('\t').collect();
    let kegg_col = headers
        .iter()
        .position(|h| {
            let h = h.trim().to_lowercase();
            h == "kegg_pathways" || h == "kegg_ko" || h == "kegg"
        })
        .ok_or_else(|| {
            MycoNoteError::InvalidFormat(format!(
                "No kegg_pathways/kegg_ko column in {}. Is this a myconote `annotate` output?",
                annotated_tsv.display()
            ))
        })?;

    let mut kos: HashSet<String> = HashSet::new();
    let mut rows = 0usize;
    for line in lines.map_while(|l| l.ok()) {
        if line.is_empty() {
            continue;
        }
        rows += 1;
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= kegg_col {
            continue;
        }
        for tok in cols[kegg_col].split(',') {
            let tok = tok.trim();
            if tok.starts_with('K')
                && tok.len() >= 2
                && tok.len() <= 8
                && tok.chars().skip(1).all(|c| c.is_ascii_digit())
            {
                kos.insert(tok.to_string());
            }
        }
    }
    Ok((kos.len(), rows))
}

/// Recursive walk for `.tsv` / `.txt` files — handles the fact that each
/// tarball extracts into a versioned subdir like
/// `busco/Y1000p_BUCO_fulltable/<species>.full_table.tsv`.
fn walk_tsv_files(root: &Path) -> Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    fn walk(p: &Path, out: &mut Vec<PathBuf>) -> std::io::Result<()> {
        for entry in std::fs::read_dir(p)?.flatten() {
            let path = entry.path();
            if path.is_dir() {
                walk(&path, out)?;
                continue;
            }
            let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
                continue;
            };
            // Skip bundled READMEs and hidden files — they are not per-species data.
            if name.eq_ignore_ascii_case("README.txt")
                || name.eq_ignore_ascii_case("README.md")
                || name.starts_with('.')
            {
                continue;
            }
            if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
                if ext == "tsv" || ext == "txt" {
                    out.push(path);
                }
            }
        }
        Ok(())
    }
    walk(root, &mut out).map_err(MycoNoteError::Io)?;
    Ok(out)
}
