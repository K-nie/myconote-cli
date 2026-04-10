/// Evidence Modeler-style consensus gene caller
///
/// Takes GFF3 predictions from multiple sources (Augustus, SNAP, etc.)
/// and produces a single high-confidence gene set by:
///
///   1. Parsing all predictions into a unified model
///   2. Grouping overlapping predictions on the same strand
///   3. Scoring each model: weight × number of supporting predictors
///   4. For each overlap group, emitting the highest-scoring model
///   5. Renaming genes with a clean sequential locus tag
///
/// Weights (configurable):
///   Augustus:       10  (most accurate ab initio tool)
///   SNAP:            3  (fast but noisier)
///   Protein align:  20  (homology evidence is strongest)
///   EST/transcript:  8  (good exon support)
use crate::parser::gff::{GFFReader, GFFRecord};
use crate::utils::error::Result;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// Evidence weights
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
pub struct EvidenceWeights {
    pub augustus: f64,
    pub snap: f64,
    pub protein: f64,
    pub est: f64,
    #[serde(default = "default_glimmer_weight")]
    pub glimmerhmm: f64,
    #[serde(default = "default_genemark_weight")]
    pub genemark: f64,
    #[serde(default = "default_miniprot_weight")]
    pub miniprot: f64,
    #[serde(default = "default_trna_weight")]
    pub trnascan: f64,
}

fn default_glimmer_weight() -> f64 {
    2.0
}
fn default_genemark_weight() -> f64 {
    5.0
}
fn default_miniprot_weight() -> f64 {
    20.0
}
fn default_trna_weight() -> f64 {
    15.0
}

impl Default for EvidenceWeights {
    fn default() -> Self {
        Self {
            augustus: 10.0,
            snap: 3.0,
            protein: 20.0,
            est: 8.0,
            glimmerhmm: 2.0,
            genemark: 5.0,
            miniprot: 20.0,
            trnascan: 15.0,
        }
    }
}

impl EvidenceWeights {
    /// Load weights from a TOML configuration file.
    ///
    /// Example TOML:
    /// ```toml
    /// augustus = 10.0
    /// snap = 3.0
    /// protein = 20.0
    /// est = 8.0
    /// glimmerhmm = 2.0
    /// genemark = 5.0
    /// miniprot = 20.0
    /// trnascan = 15.0
    /// ```
    pub fn from_toml(path: &std::path::Path) -> crate::utils::error::Result<Self> {
        let contents =
            std::fs::read_to_string(path).map_err(crate::utils::error::MycoNoteError::Io)?;
        let weights: EvidenceWeights = toml::from_str(&contents).map_err(|e| {
            crate::utils::error::MycoNoteError::InvalidFormat(format!(
                "Failed to parse weights TOML: {}",
                e
            ))
        })?;
        Ok(weights)
    }

    /// Write current weights to a TOML file (for reproducibility).
    pub fn write_toml(&self, path: &std::path::Path) -> crate::utils::error::Result<()> {
        let toml_str = toml::to_string_pretty(self).map_err(|e| {
            crate::utils::error::MycoNoteError::InvalidFormat(format!(
                "Failed to serialize weights: {}",
                e
            ))
        })?;
        std::fs::write(path, toml_str).map_err(crate::utils::error::MycoNoteError::Io)?;
        Ok(())
    }

    /// Get weight by source name string.
    pub fn weight_for(&self, source: &str) -> f64 {
        match source.to_lowercase().as_str() {
            "augustus" => self.augustus,
            "snap" => self.snap,
            "protein" | "protein2genome" => self.protein,
            "est" | "transcript" => self.est,
            "glimmerhmm" | "glimmer" => self.glimmerhmm,
            "genemark" | "genemark-es" | "genemark-et" => self.genemark,
            "miniprot" | "protein_evidence" => self.miniprot,
            "trnascan" | "trnascan-se" => self.trnascan,
            _ => 1.0, // default weight for unknown sources
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Gene model
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
#[allow(dead_code)]
struct GeneModel {
    seqid: String,
    start: u64,
    end: u64,
    strand: char,
    score: f64,
    source: String,
    /// All GFF3 records belonging to this gene (gene + mRNA + CDS + exon)
    records: Vec<GFFRecord>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Parsing helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Load a GFF3 file and extract top-level gene models with their children.
fn load_genes(path: &Path, source_weight: f64, source_name: &str) -> Result<Vec<GeneModel>> {
    let all: Vec<GFFRecord> = GFFReader::from_path(path)?.filter_map(|r| r.ok()).collect();

    // Map ID → record for child lookup
    let mut id_map: HashMap<String, usize> = HashMap::new();
    for (i, rec) in all.iter().enumerate() {
        if let Some(id) = rec.id() {
            id_map.insert(id.clone(), i);
        }
    }

    // Build gene models: find all gene-level features
    let mut models: Vec<GeneModel> = Vec::new();

    for rec in all.iter().filter(|r| r.feature_type == "gene") {
        let gene_id = match rec.id() {
            Some(id) => id.clone(),
            None => continue,
        };

        // Collect all children recursively
        let mut members = vec![rec.clone()];
        collect_children(&gene_id, &all, &mut members);

        models.push(GeneModel {
            seqid: rec.seqid.clone(),
            start: rec.start,
            end: rec.end,
            strand: rec.strand,
            score: source_weight,
            source: source_name.to_string(),
            records: members,
        });
    }

    Ok(models)
}

fn collect_children(parent_id: &str, all: &[GFFRecord], out: &mut Vec<GFFRecord>) {
    for rec in all {
        if rec.parent().map(|p| p.as_str()) == Some(parent_id) {
            let cid = rec.id().cloned();
            out.push(rec.clone());
            if let Some(id) = cid {
                collect_children(&id, all, out);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Overlap resolution
// ─────────────────────────────────────────────────────────────────────────────

/// For a set of models on the same sequence, group overlapping models and
/// keep the highest-scoring one per overlap group.
fn resolve_overlaps(mut models: Vec<GeneModel>) -> Vec<GeneModel> {
    // Sort by start coordinate
    models.sort_by(|a, b| a.start.cmp(&b.start).then(a.seqid.cmp(&b.seqid)));

    let mut kept: Vec<GeneModel> = Vec::new();
    let mut used = vec![false; models.len()];

    for i in 0..models.len() {
        if used[i] {
            continue;
        }

        // Find all models overlapping models[i]
        let mut group: Vec<usize> = vec![i];
        let mut group_end = models[i].end;

        for j in (i + 1)..models.len() {
            if models[j].seqid != models[i].seqid {
                continue;
            }
            if models[j].start > group_end {
                break;
            }
            if models[j].strand != models[i].strand {
                continue;
            }
            // Overlapping on same strand
            group.push(j);
            group_end = group_end.max(models[j].end);
        }

        // Mark all in group as used
        for &idx in &group {
            used[idx] = true;
        }

        // Keep highest-scoring model in the group
        if let Some(best) = group.into_iter().max_by(|&a, &b| {
            models[a]
                .score
                .partial_cmp(&models[b].score)
                .unwrap_or(std::cmp::Ordering::Equal)
        }) {
            kept.push(models[best].clone());
        }
    }

    kept
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Merge predictions from multiple GFF3 files into a single consensus GFF3.
///
/// `inputs` is a list of `(gff3_path, source_name, weight)`.
/// `locus_prefix` sets the gene ID prefix (e.g. "MYCO" → "MYCO_000001").
pub fn merge_predictions(
    inputs: &[(&Path, &str, f64)],
    output_gff: &Path,
    locus_prefix: &str,
) -> Result<usize> {
    let mut all_models: Vec<GeneModel> = Vec::new();

    for &(path, source, weight) in inputs {
        match load_genes(path, weight, source) {
            Ok(models) => {
                println!("  {} predictions from {}", models.len(), source);
                all_models.extend(models);
            }
            Err(e) => eprintln!("  ⚠  Skipping {} ({})", source, e),
        }
    }

    if all_models.is_empty() {
        return Err(crate::utils::error::MycoNoteError::InvalidFormat(
            "No gene predictions from any source.".to_string(),
        ));
    }

    // Resolve overlapping predictions (per chromosome)
    let mut by_chr: HashMap<String, Vec<GeneModel>> = HashMap::new();
    for m in all_models {
        by_chr.entry(m.seqid.clone()).or_default().push(m);
    }

    let mut final_models: Vec<GeneModel> = Vec::new();
    for (_, models) in by_chr {
        final_models.extend(resolve_overlaps(models));
    }

    // Sort by chromosome then position
    final_models.sort_by(|a, b| a.seqid.cmp(&b.seqid).then(a.start.cmp(&b.start)));

    println!(
        "  Consensus: {} genes after overlap resolution",
        final_models.len()
    );

    // Write GFF3 with clean sequential locus tags
    let mut out =
        std::fs::File::create(output_gff).map_err(crate::utils::error::MycoNoteError::Io)?;
    writeln!(out, "##gff-version 3").map_err(crate::utils::error::MycoNoteError::Io)?;

    let total = final_models.len();
    let pad = total.to_string().len().max(6);

    for (gene_num, model) in final_models.iter().enumerate() {
        let locus_tag = format!("{}_{:0>pad$}", locus_prefix, gene_num + 1, pad = pad);
        let gene_id = locus_tag.clone();
        let mrna_id = format!("{}-mRNA1", locus_tag);

        // Rewrite records with unique IDs per feature
        let mut cds_counter = 0usize;
        let mut exon_counter = 0usize;
        let mut other_counter = 0usize;

        for rec in &model.records {
            let mut r = rec.clone();

            match rec.feature_type.as_str() {
                "gene" => {
                    r.attributes.insert("ID".into(), gene_id.clone());
                    r.attributes.insert("locus_tag".into(), locus_tag.clone());
                }
                "mRNA" | "transcript" => {
                    r.attributes.insert("ID".into(), mrna_id.clone());
                    r.attributes.insert("Parent".into(), gene_id.clone());
                }
                "CDS" => {
                    cds_counter += 1;
                    r.attributes.insert("Parent".into(), mrna_id.clone());
                    r.attributes
                        .insert("ID".into(), format!("{}-CDS-{}", locus_tag, cds_counter));
                }
                "exon" => {
                    exon_counter += 1;
                    r.attributes.insert("Parent".into(), mrna_id.clone());
                    r.attributes
                        .insert("ID".into(), format!("{}-exon-{}", locus_tag, exon_counter));
                }
                "five_prime_UTR" | "three_prime_UTR" | "UTR" => {
                    r.attributes.insert("Parent".into(), mrna_id.clone());
                    other_counter += 1;
                    r.attributes.insert(
                        "ID".into(),
                        format!("{}-{}-{}", locus_tag, rec.feature_type, other_counter),
                    );
                }
                _ => {
                    // For any other child feature, reparent to the mRNA
                    other_counter += 1;
                    if rec.parent().is_some() {
                        r.attributes.insert("Parent".into(), mrna_id.clone());
                    }
                    r.attributes.insert(
                        "ID".into(),
                        format!("{}-{}-{}", locus_tag, rec.feature_type, other_counter),
                    );
                }
            }

            writeln!(out, "{}", r.to_gff3_line())
                .map_err(crate::utils::error::MycoNoteError::Io)?;
        }
    }

    Ok(total)
}
