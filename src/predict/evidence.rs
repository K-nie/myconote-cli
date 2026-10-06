/// Evidence Modeler-style consensus gene caller
///
/// Takes GFF3 predictions from multiple sources (Augustus, SNAP, etc.)
/// and produces a single high-confidence gene set by:
///
///   1. Parsing all predictions into a unified model
///   2. Grouping overlapping predictions on the same strand
///   3. Scoring each model by agreement-weighted support: its own source
///      weight plus, for every other model in the group, that model's
///      weight scaled by how much their coding bases agree (CDS Jaccard).
///      A locus two high-weight predictors call identically therefore
///      outranks one high-weight predictor calling it alone, and the
///      winner is the model whose boundaries the evidence most agrees on.
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

/// Extract this model's CDS intervals (1-based inclusive), sorted by start.
/// Falls back to the gene span when a model carries no CDS (e.g. a tRNA
/// gene from tRNAscan), so coding-base agreement is always defined.
fn cds_intervals(m: &GeneModel) -> Vec<(u64, u64)> {
    let mut iv: Vec<(u64, u64)> = m
        .records
        .iter()
        .filter(|r| r.feature_type == "CDS")
        .map(|r| (r.start, r.end))
        .collect();
    if iv.is_empty() {
        iv.push((m.start, m.end));
    }
    iv.sort_by_key(|&(s, _)| s);
    iv
}

/// Overlapping coding bases between two start-sorted interval lists.
fn interval_overlap_bp(a: &[(u64, u64)], b: &[(u64, u64)]) -> u64 {
    let (mut i, mut j, mut total) = (0usize, 0usize, 0u64);
    while i < a.len() && j < b.len() {
        let lo = a[i].0.max(b[j].0);
        let hi = a[i].1.min(b[j].1);
        if lo <= hi {
            total += hi - lo + 1;
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    total
}

/// Jaccard of covered coding nucleotides: |A∩B| / |A∪B|. 1.0 means the two
/// models call exactly the same coding bases (identical boundaries); 0.0
/// means they share none. This is the agreement term in the consensus score.
fn cds_jaccard(a: &[(u64, u64)], b: &[(u64, u64)]) -> f64 {
    let inter = interval_overlap_bp(a, b);
    if inter == 0 {
        return 0.0;
    }
    let len: fn(&[(u64, u64)]) -> u64 = |iv| iv.iter().map(|&(s, e)| e.saturating_sub(s) + 1).sum();
    let union = len(a) + len(b) - inter;
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

/// For a set of models on the same sequence, group overlapping models and
/// keep the single best-supported one per overlap group.
///
/// The winner is chosen by an agreement-weighted consensus score rather than
/// by raw source weight. For a candidate model `x` in a group, the score is
///
///   score(x) = weight(x) + Σ_{y≠x in group} weight(y) · cds_jaccard(x, y)
///
/// so a model that several predictors (each by their own weight) agree with
/// at the coding-base level outranks a lone high-weight call. Because the
/// Jaccard term peaks when boundaries coincide, the emitted model tends to
/// carry the consensus exon boundaries — the lever for exact-match (strict)
/// gene accuracy. Ties fall back to raw source weight, then to the longer
/// total CDS, so the ordering is deterministic.
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

        // Precompute CDS intervals once per group member, then score each by
        // agreement-weighted consensus and keep the best.
        let cds: Vec<Vec<(u64, u64)>> = group
            .iter()
            .map(|&idx| cds_intervals(&models[idx]))
            .collect();

        let consensus_score = |x: usize| -> f64 {
            let mut s = models[group[x]].score;
            for y in 0..group.len() {
                if y == x {
                    continue;
                }
                s += models[group[y]].score * cds_jaccard(&cds[x], &cds[y]);
            }
            s
        };

        let cds_len =
            |x: usize| -> u64 { cds[x].iter().map(|&(s, e)| e.saturating_sub(s) + 1).sum() };

        if let Some(best) = (0..group.len()).max_by(|&a, &b| {
            consensus_score(a)
                .partial_cmp(&consensus_score(b))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(
                    models[group[a]]
                        .score
                        .partial_cmp(&models[group[b]].score)
                        .unwrap_or(std::cmp::Ordering::Equal),
                )
                .then(cds_len(a).cmp(&cds_len(b)))
        }) {
            kept.push(models[group[best]].clone());
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

        // First pass: rewrite records + collect CDS spans and codon presence
        // so we can synthesize the missing exon / start_codon / stop_codon
        // records after the model's own records are emitted. Rationale:
        //
        //   * F3a — canonical GFF3 expects `exon` features; SNAP / Augustus
        //     with `--UTR=off` emit CDS only, so downstream tools (IGV,
        //     table2asn, JBrowse2) that read exons see empty transcripts.
        //     We derive one exon per CDS span when no explicit exon
        //     records exist.
        //   * F3b — different predictors label the transcript feature
        //     "mRNA" vs. "transcript"; normalize the emitted line to
        //     "mRNA" so every consensus gene is SOFA-compliant.
        //   * F3c — codon markers are optional in some predictor output
        //     even though the underlying CDS is complete. Synthesize a
        //     coordinate-only start_codon / stop_codon from the outermost
        //     CDS span when the model does not carry one. The 3-bp span
        //     is derived from the strand-aware 5' / 3' end of the CDS
        //     union; we do not verify the FASTA (evidence.rs has no
        //     genome handle), so the synthesized record is a placeholder
        //     that fixes tool interop and does not assert biology.
        let mut cds_counter = 0usize;
        let mut exon_counter = 0usize;
        let mut other_counter = 0usize;

        let mut cds_spans: Vec<(u64, u64, char)> = Vec::new(); // (start, end, strand)
        let mut has_exon_records = false;
        let mut has_start_codon = false;
        let mut has_stop_codon = false;
        let mut mrna_source: Option<String> = None;

        for rec in &model.records {
            let mut r = rec.clone();

            match rec.feature_type.as_str() {
                "gene" => {
                    r.attributes.insert("ID".into(), gene_id.clone());
                    r.attributes.insert("locus_tag".into(), locus_tag.clone());
                }
                "mRNA" | "transcript" => {
                    // F3b: normalize whichever label the predictor used.
                    r.feature_type = "mRNA".to_string();
                    r.attributes.insert("ID".into(), mrna_id.clone());
                    r.attributes.insert("Parent".into(), gene_id.clone());
                    mrna_source = Some(r.source.clone());
                }
                "CDS" => {
                    cds_counter += 1;
                    r.attributes.insert("Parent".into(), mrna_id.clone());
                    r.attributes
                        .insert("ID".into(), format!("{}-CDS-{}", locus_tag, cds_counter));
                    cds_spans.push((rec.start, rec.end, rec.strand));
                }
                "exon" => {
                    has_exon_records = true;
                    exon_counter += 1;
                    r.attributes.insert("Parent".into(), mrna_id.clone());
                    r.attributes
                        .insert("ID".into(), format!("{}-exon-{}", locus_tag, exon_counter));
                }
                "start_codon" => {
                    has_start_codon = true;
                    other_counter += 1;
                    r.attributes.insert("Parent".into(), mrna_id.clone());
                    r.attributes.insert(
                        "ID".into(),
                        format!("{}-start_codon-{}", locus_tag, other_counter),
                    );
                }
                "stop_codon" => {
                    has_stop_codon = true;
                    other_counter += 1;
                    r.attributes.insert("Parent".into(), mrna_id.clone());
                    r.attributes.insert(
                        "ID".into(),
                        format!("{}-stop_codon-{}", locus_tag, other_counter),
                    );
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

        // Synthesize missing exon / codon records.
        if !cds_spans.is_empty() {
            let source = mrna_source.as_deref().unwrap_or(&model.source);

            if !has_exon_records {
                // F3a: one exon per CDS span; UTRs are absent, so exon
                // coordinates equal CDS coordinates.
                for (i, (s, e, strand)) in cds_spans.iter().enumerate() {
                    exon_counter += 1;
                    writeln!(
                        out,
                        "{}\t{}\texon\t{}\t{}\t.\t{}\t.\tID={}-exon-{};Parent={}",
                        model.seqid,
                        source,
                        s,
                        e,
                        strand,
                        locus_tag,
                        i + 1,
                        mrna_id
                    )
                    .map_err(crate::utils::error::MycoNoteError::Io)?;
                }
            }

            // F3c: derive missing codon markers from the 5' / 3' extremes
            // of the CDS union (strand-aware). Placeholder only — not FASTA-
            // verified — so downstream tools that expect the record type
            // are satisfied without falsely asserting the biology of the
            // ATG / STOP nucleotides.
            let cds_min = cds_spans.iter().map(|c| c.0).min().unwrap();
            let cds_max = cds_spans.iter().map(|c| c.1).max().unwrap();
            let strand = cds_spans[0].2;

            if !has_start_codon {
                let (s, e) = if strand == '-' {
                    (cds_max.saturating_sub(2), cds_max)
                } else {
                    (cds_min, (cds_min + 2).max(cds_min))
                };
                other_counter += 1;
                writeln!(
                    out,
                    "{}\t{}\tstart_codon\t{}\t{}\t.\t{}\t0\tID={}-start_codon-{};Parent={};Note=synthesized_from_cds_bounds",
                    model.seqid, source, s, e, strand, locus_tag, other_counter, mrna_id
                )
                .map_err(crate::utils::error::MycoNoteError::Io)?;
            }

            if !has_stop_codon {
                let (s, e) = if strand == '-' {
                    (cds_min, (cds_min + 2).max(cds_min))
                } else {
                    (cds_max.saturating_sub(2), cds_max)
                };
                other_counter += 1;
                writeln!(
                    out,
                    "{}\t{}\tstop_codon\t{}\t{}\t.\t{}\t0\tID={}-stop_codon-{};Parent={};Note=synthesized_from_cds_bounds",
                    model.seqid, source, s, e, strand, locus_tag, other_counter, mrna_id
                )
                .map_err(crate::utils::error::MycoNoteError::Io)?;
            }
        }
    }

    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gene_model(start: u64, end: u64, score: f64, cds: &[(u64, u64)]) -> GeneModel {
        let mut records = Vec::new();
        for &(s, e) in cds {
            records.push(GFFRecord {
                seqid: "chr1".into(),
                source: "test".into(),
                feature_type: "CDS".into(),
                start: s,
                end: e,
                score: None,
                strand: '+',
                phase: Some(0),
                attributes: HashMap::new(),
            });
        }
        GeneModel {
            seqid: "chr1".into(),
            start,
            end,
            strand: '+',
            score,
            source: "test".into(),
            records,
        }
    }

    #[test]
    fn cds_jaccard_identical_is_one() {
        let a = [(100u64, 200u64)];
        assert!((cds_jaccard(&a, &a) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn cds_jaccard_disjoint_is_zero() {
        assert_eq!(cds_jaccard(&[(100, 200)], &[(300, 400)]), 0.0);
    }

    // Two low-weight predictors that agree on identical boundaries must beat a
    // single higher-weight predictor with different boundaries. Under the old
    // raw-weight rule the lone score-10 model won; the agreement-weighted
    // consensus flips the winner to the score-6 pair's boundaries.
    #[test]
    fn agreement_outranks_lone_higher_weight() {
        let a = gene_model(100, 200, 6.0, &[(100, 200)]);
        let b = gene_model(100, 200, 6.0, &[(100, 200)]);
        let c = gene_model(150, 260, 10.0, &[(150, 260)]);

        let kept = resolve_overlaps(vec![a, b, c]);
        assert_eq!(kept.len(), 1, "one locus → one gene");
        let w = &kept[0];
        assert_eq!(
            (w.start, w.end),
            (100, 200),
            "consensus should keep the agreed-upon boundaries, not the lone higher-weight model"
        );
    }

    // Sanity: raw max-weight would have kept (150, 260) here.
    #[test]
    fn lone_model_passes_through_unchanged() {
        let c = gene_model(150, 260, 10.0, &[(150, 260)]);
        let kept = resolve_overlaps(vec![c]);
        assert_eq!(kept.len(), 1);
        assert_eq!((kept[0].start, kept[0].end), (150, 260));
    }
}
