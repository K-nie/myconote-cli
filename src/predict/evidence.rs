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
// Consensus false-positive filter (v0.7.7)
// ─────────────────────────────────────────────────────────────────────────────

/// Post-consensus false-positive filter.
///
/// `resolve_overlaps` keeps exactly one model per overlap group, unconditionally.
/// On the 6-genome HTCondor benchmark that admitted false-positive genes: turning
/// on a weak secondary predictor (SNAP + GeneMark) raised the raw gene count
/// 5465 → 5781 but LOWERED loose gene-F1 0.906 → 0.892, because a locus called by
/// a single low-weight predictor is emitted whether or not anything corroborates
/// it. This filter drops those singletons.
///
/// A group's winning model is emitted only if EITHER
///   * it is backed by `>= min_predictor_support` distinct predictor sources
///     whose CDS overlaps the winner — the winner's own source counts as one,
///     plus every other source in the group with `cds_jaccard > 0` to it; OR
///   * (rescue) its agreement-weighted consensus score is `>= min_consensus_score`.
///
/// The defaults (`min_predictor_support = 1`, `min_consensus_score = None`)
/// reproduce the pre-filter behaviour byte-for-byte: support is always >= 1 (the
/// winner counts itself), so every winner is kept. The recommended fungal
/// setting is `min_predictor_support = 2`, which requires a second predictor to
/// corroborate each locus; `min_consensus_score` then acts as a rescue so a lone
/// but very-high-scoring call is not discarded.
#[derive(Debug, Clone)]
pub struct ConsensusFilter {
    /// Minimum number of distinct predictor sources that must support a locus.
    /// 1 = keep everything (default, pre-filter behaviour).
    pub min_predictor_support: usize,
    /// Optional rescue threshold on the agreement-weighted consensus score. A
    /// winner that fails the support test is still kept when its score reaches
    /// this value. `None` disables the rescue.
    pub min_consensus_score: Option<f64>,
}

impl Default for ConsensusFilter {
    fn default() -> Self {
        Self {
            min_predictor_support: 1,
            min_consensus_score: None,
        }
    }
}

impl ConsensusFilter {
    /// Decide whether a winning model is kept given its support count (distinct
    /// corroborating sources, incl. itself) and its consensus score.
    fn keeps(&self, support: usize, consensus_score: f64) -> bool {
        if support >= self.min_predictor_support {
            return true;
        }
        match self.min_consensus_score {
            Some(min) => consensus_score >= min,
            None => false,
        }
    }

    /// True when this filter cannot drop anything — lets callers skip the
    /// per-locus bookkeeping and logging on the default (opt-out) path.
    fn is_noop(&self) -> bool {
        self.min_predictor_support <= 1
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
/// Back-compat shim used by the existing overlap-resolution tests, which assert
/// pre-filter behaviour. Production merges go through `merge_predictions_filtered`
/// → `resolve_overlaps_filtered` directly, so this is only compiled for tests.
#[cfg(test)]
fn resolve_overlaps(models: Vec<GeneModel>) -> Vec<GeneModel> {
    resolve_overlaps_filtered(models, &ConsensusFilter::default()).0
}

/// Like `resolve_overlaps` but applies a `ConsensusFilter` to each group's
/// winner. Returns `(kept, dropped)` where `dropped` counts loci discarded by
/// the filter. With a no-op filter this is identical to `resolve_overlaps` and
/// `dropped == 0`.
fn resolve_overlaps_filtered(
    mut models: Vec<GeneModel>,
    filter: &ConsensusFilter,
) -> (Vec<GeneModel>, usize) {
    // Sort by start coordinate
    models.sort_by(|a, b| a.start.cmp(&b.start).then(a.seqid.cmp(&b.seqid)));

    let mut kept: Vec<GeneModel> = Vec::new();
    let mut dropped = 0usize;
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
            // Fast path: the default no-op filter keeps every winner, so skip
            // the support count entirely and preserve byte-identical output.
            if filter.is_noop() {
                kept.push(models[group[best]].clone());
                continue;
            }

            let winner_score = consensus_score(best);

            // Support = distinct predictor sources whose CDS overlaps the
            // winner, counting the winner's own source. Two models from the
            // *same* source do not inflate support (we want corroboration from
            // an independent predictor, which is what the benchmark FP source
            // lacked).
            let mut sources: std::collections::HashSet<&str> =
                std::collections::HashSet::new();
            sources.insert(models[group[best]].source.as_str());
            for y in 0..group.len() {
                if y == best {
                    continue;
                }
                if cds_jaccard(&cds[best], &cds[y]) > 0.0 {
                    sources.insert(models[group[y]].source.as_str());
                }
            }
            let support = sources.len();

            if filter.keeps(support, winner_score) {
                kept.push(models[group[best]].clone());
            } else {
                dropped += 1;
                // No silent data loss: name every locus the filter removes so
                // non-TTY cluster logs record what the gene set lost and why.
                eprintln!(
                    "  ⚠  consensus filter dropped {}:{}-{} ({}) — support={} source(s), \
                     score={:.2} (need support ≥ {}{})",
                    models[group[best]].seqid,
                    models[group[best]].start,
                    models[group[best]].end,
                    models[group[best]].source,
                    support,
                    winner_score,
                    filter.min_predictor_support,
                    match filter.min_consensus_score {
                        Some(m) => format!(" or score ≥ {:.2}", m),
                        None => String::new(),
                    }
                );
            }
        }
    }

    (kept, dropped)
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Merge predictions from multiple GFF3 files into a single consensus GFF3.
///
/// `inputs` is a list of `(gff3_path, source_name, weight)`.
/// `locus_prefix` sets the gene ID prefix (e.g. "MYCO" → "MYCO_000001").
///
/// This keeps the pre-filter signature; it delegates to
/// `merge_predictions_filtered` with a no-op `ConsensusFilter`, so the output is
/// byte-identical to earlier releases.
pub fn merge_predictions(
    inputs: &[(&Path, &str, f64)],
    output_gff: &Path,
    locus_prefix: &str,
) -> Result<usize> {
    merge_predictions_filtered(inputs, output_gff, locus_prefix, &ConsensusFilter::default())
}

/// Merge predictions into a consensus GFF3, dropping false-positive loci that
/// do not meet `filter`. With `ConsensusFilter::default()` this is identical to
/// `merge_predictions`.
pub fn merge_predictions_filtered(
    inputs: &[(&Path, &str, f64)],
    output_gff: &Path,
    locus_prefix: &str,
    filter: &ConsensusFilter,
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
    let mut total_dropped = 0usize;
    // Iterate chromosomes in sorted order so the dropped-locus log (stderr) is
    // deterministic regardless of HashMap iteration order.
    let mut chr_keys: Vec<String> = by_chr.keys().cloned().collect();
    chr_keys.sort();
    for chr in chr_keys {
        let models = by_chr.remove(&chr).unwrap();
        let (kept, dropped) = resolve_overlaps_filtered(models, filter);
        total_dropped += dropped;
        final_models.extend(kept);
    }

    // Sort by chromosome then position
    final_models.sort_by(|a, b| a.seqid.cmp(&b.seqid).then(a.start.cmp(&b.start)));

    if total_dropped > 0 {
        println!(
            "  Consensus filter: dropped {} false-positive locus/loci \
             (min-predictor-support={}{})",
            total_dropped,
            filter.min_predictor_support,
            match filter.min_consensus_score {
                Some(m) => format!(", min-consensus-score={:.2}", m),
                None => String::new(),
            }
        );
    }

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

        // Fold the stop codon into the terminal CDS (and its exon) before we
        // emit anything, matching NCBI / RefSeq convention. Doing it here means
        // the CDS spans collected below — and therefore any synthesized exon —
        // already carry the stop-codon 3 bp. A model is one transcript, so the
        // whole record set belongs to a single mRNA.
        let mut norm_records = model.records.clone();
        crate::parser::gff::include_stop_codon_in_cds(&mut norm_records);

        for rec in &norm_records {
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

    /// Like `gene_model` but lets a test set the predictor source name, which
    /// is what the consensus filter counts for support.
    fn gene_model_src(
        start: u64,
        end: u64,
        score: f64,
        cds: &[(u64, u64)],
        source: &str,
    ) -> GeneModel {
        let mut m = gene_model(start, end, score, cds);
        m.source = source.to_string();
        m
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

    // Pull the (start, end) of the single CDS row for a given locus tag out of
    // an emitted GFF3, so the stop-codon fold can be asserted on real output.
    fn cds_bounds(gff: &str, locus_tag: &str) -> (u64, u64) {
        for line in gff.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() == 9 && f[2] == "CDS" && f[8].contains(locus_tag) {
                return (f[3].parse().unwrap(), f[4].parse().unwrap());
            }
        }
        panic!("no CDS row for {}", locus_tag);
    }

    fn count_feature(gff: &str, feature: &str) -> usize {
        gff.lines()
            .filter(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                f.len() == 9 && f[2] == feature
            })
            .count()
    }

    // End-to-end: a + strand and a - strand gene, each with the stop codon as a
    // separate 3-bp feature outside the CDS (exactly what Augustus emits), run
    // through the real consensus writer. The emitted terminal CDS must now
    // swallow the stop codon to match NCBI / RefSeq.
    #[test]
    fn merge_predictions_folds_stop_codon_into_cds_both_strands() {
        let dir = tempfile::tempdir().unwrap();

        // + strand gene: CDS 1802-2950, stop_codon 2951-2953.
        // - strand gene: CDS 4000-5000, stop_codon 3997-3999.
        let input = dir.path().join("in.gff3");
        std::fs::write(
            &input,
            "##gff-version 3\n\
chrA\tAugustus\tgene\t1802\t2953\t.\t+\t.\tID=gP\n\
chrA\tAugustus\tmRNA\t1802\t2953\t.\t+\t.\tID=gP.t1;Parent=gP\n\
chrA\tAugustus\tCDS\t1802\t2950\t.\t+\t0\tID=gP.cds;Parent=gP.t1\n\
chrA\tAugustus\texon\t1802\t2950\t.\t+\t.\tID=gP.exon;Parent=gP.t1\n\
chrA\tAugustus\tstop_codon\t2951\t2953\t.\t+\t0\tID=gP.stop;Parent=gP.t1\n\
chrA\tAugustus\tgene\t3997\t5000\t.\t-\t.\tID=gM\n\
chrA\tAugustus\tmRNA\t3997\t5000\t.\t-\t.\tID=gM.t1;Parent=gM\n\
chrA\tAugustus\tCDS\t4000\t5000\t.\t-\t0\tID=gM.cds;Parent=gM.t1\n\
chrA\tAugustus\texon\t4000\t5000\t.\t-\t.\tID=gM.exon;Parent=gM.t1\n\
chrA\tAugustus\tstop_codon\t3997\t3999\t.\t-\t0\tID=gM.stop;Parent=gM.t1\n",
        )
        .unwrap();

        let out = dir.path().join("consensus.gff3");
        let inputs: Vec<(&Path, &str, f64)> = vec![(input.as_path(), "augustus", 10.0)];
        let n = merge_predictions(&inputs, &out, "TEST").unwrap();
        assert_eq!(n, 2, "two genes in, two genes out");

        let gff = std::fs::read_to_string(&out).unwrap();

        // Locus tags are assigned by genomic order: TEST_000001 (+), TEST_000002 (-).
        let plus = cds_bounds(&gff, "TEST_000001");
        assert_eq!(plus.1, 2953, "+ strand CDS end extends through stop codon");
        assert_eq!(plus, (1802, 2953));

        let minus = cds_bounds(&gff, "TEST_000002");
        assert_eq!(
            minus.0, 3997,
            "- strand CDS start extends through stop codon"
        );
        assert_eq!(minus, (3997, 5000));

        // The informational stop_codon rows are preserved as-is (one per gene).
        assert_eq!(count_feature(&gff, "stop_codon"), 2, "stop_codon rows kept");
        assert!(
            gff.lines()
                .any(|l| l.contains("\tstop_codon\t2951\t2953\t")),
            "+ strand stop_codon row preserved unchanged"
        );
        assert!(
            gff.lines()
                .any(|l| l.contains("\tstop_codon\t3997\t3999\t")),
            "- strand stop_codon row preserved unchanged"
        );

        // The exon coincident with each terminal CDS must track the extension,
        // keeping exon ⊇ CDS (what RefSeq single-exon genes look like).
        assert!(
            gff.lines()
                .any(|l| l.contains("\texon\t1802\t2953\t") && l.contains("TEST_000001")),
            "+ strand exon extended through stop codon"
        );
        assert!(
            gff.lines()
                .any(|l| l.contains("\texon\t3997\t5000\t") && l.contains("TEST_000002")),
            "- strand exon extended through stop codon"
        );
    }

    // Faithful reproduction of the REAL cluster run: Augustus emits the
    // transcript as `transcript` (not `mRNA`), the start_codon/CDS/stop_codon
    // children carry only a Parent (no ID), ids follow the c0_g1 / c0_g1.t1
    // scheme, AND a second predictor (GeneMark) calls the same locus so
    // resolve_overlaps actually runs and picks a winner. This is the shape the
    // synthetic test above failed to capture.
    #[test]
    fn merge_predictions_folds_stop_codon_real_augustus_shape() {
        let dir = tempfile::tempdir().unwrap();

        let aug = dir.path().join("augustus.gff3");
        std::fs::write(
            &aug,
            "##gff-version 3\n\
scaffold_001\tAUGUSTUS\tgene\t1802\t2953\t1\t+\t.\tID=c0_g1\n\
scaffold_001\tAUGUSTUS\ttranscript\t1802\t2953\t1\t+\t.\tID=c0_g1.t1;Parent=c0_g1\n\
scaffold_001\tAUGUSTUS\tstart_codon\t1802\t1804\t.\t+\t0\tParent=c0_g1.t1\n\
scaffold_001\tAUGUSTUS\tCDS\t1802\t2950\t1\t+\t0\tID=c0_g1.t1.cds;Parent=c0_g1.t1\n\
scaffold_001\tAUGUSTUS\tstop_codon\t2951\t2953\t.\t+\t0\tParent=c0_g1.t1\n",
        )
        .unwrap();

        // GeneMark calls the same locus with slightly different boundaries.
        let gm = dir.path().join("genemark.gff3");
        std::fs::write(
            &gm,
            "##gff-version 3\n\
scaffold_001\tGeneMark.hmm\tgene\t1805\t2953\t.\t+\t.\tID=gm_1\n\
scaffold_001\tGeneMark.hmm\tmRNA\t1805\t2953\t.\t+\t.\tID=gm_1.t1;Parent=gm_1\n\
scaffold_001\tGeneMark.hmm\tCDS\t1805\t2950\t.\t+\t0\tID=gm_1.cds;Parent=gm_1.t1\n\
scaffold_001\tGeneMark.hmm\tstop_codon\t2951\t2953\t.\t+\t0\tParent=gm_1.t1\n",
        )
        .unwrap();

        let out = dir.path().join("consensus.gff3");
        let inputs: Vec<(&Path, &str, f64)> = vec![
            (aug.as_path(), "augustus", 10.0),
            (gm.as_path(), "genemark", 5.0),
        ];
        let n = merge_predictions(&inputs, &out, "SCE").unwrap();
        assert_eq!(n, 1, "one consensus gene for the shared locus");

        let gff = std::fs::read_to_string(&out).unwrap();
        let cds = cds_bounds(&gff, "SCE_000001");
        assert_eq!(
            cds.1, 2953,
            "real-shape consensus CDS must fold the stop codon (got {:?})\n{}",
            cds, gff
        );
    }

    // ── Feature 1: consensus false-positive filter ────────────────────────────

    // A locus called by one low-weight predictor alone is kept under the default
    // (support ≥ 1) but dropped once support ≥ 2 is required.
    #[test]
    fn singleton_dropped_at_support_2_kept_at_support_1() {
        let lone = gene_model_src(100, 200, 3.0, &[(100, 200)], "snap");

        let (kept1, dropped1) =
            resolve_overlaps_filtered(vec![lone.clone()], &ConsensusFilter::default());
        assert_eq!(kept1.len(), 1, "default keeps the lone locus");
        assert_eq!(dropped1, 0);

        let f2 = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: None,
        };
        let (kept2, dropped2) = resolve_overlaps_filtered(vec![lone], &f2);
        assert_eq!(kept2.len(), 0, "lone low-weight predictor dropped at support ≥ 2");
        assert_eq!(dropped2, 1);
    }

    // Two independent predictors calling the same locus → support 2, kept at
    // both thresholds.
    #[test]
    fn multi_predictor_locus_kept_at_both_thresholds() {
        let a = gene_model_src(100, 200, 10.0, &[(100, 200)], "augustus");
        let b = gene_model_src(100, 200, 5.0, &[(100, 200)], "genemark");

        let (k1, d1) = resolve_overlaps_filtered(
            vec![a.clone(), b.clone()],
            &ConsensusFilter::default(),
        );
        assert_eq!(k1.len(), 1);
        assert_eq!(d1, 0);

        let f2 = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: None,
        };
        let (k2, d2) = resolve_overlaps_filtered(vec![a, b], &f2);
        assert_eq!(k2.len(), 1, "two independent predictors corroborate → kept at support ≥ 2");
        assert_eq!(d2, 0);
    }

    // Two models from the *same* predictor are not independent corroboration, so
    // support stays 1 and the locus is dropped at support ≥ 2.
    #[test]
    fn same_source_twice_does_not_inflate_support() {
        let a = gene_model_src(100, 200, 3.0, &[(100, 200)], "snap");
        let b = gene_model_src(100, 200, 3.0, &[(100, 200)], "snap");
        let f2 = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: None,
        };
        let (k, d) = resolve_overlaps_filtered(vec![a, b], &f2);
        assert_eq!(k.len(), 0, "same-source duplicates are not independent support");
        assert_eq!(d, 1);
    }

    // A lone but very-high-scoring call is rescued by --min-consensus-score even
    // when it fails the support threshold, and dropped once the bar exceeds it.
    #[test]
    fn high_score_singleton_rescued_by_min_consensus_score() {
        let lone = gene_model_src(100, 200, 20.0, &[(100, 200)], "protein");

        let rescue = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: Some(15.0),
        };
        let (k, d) = resolve_overlaps_filtered(vec![lone.clone()], &rescue);
        assert_eq!(k.len(), 1, "score 20 ≥ 15 rescues the lone high-confidence locus");
        assert_eq!(d, 0);

        let strict = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: Some(25.0),
        };
        let (k2, d2) = resolve_overlaps_filtered(vec![lone], &strict);
        assert_eq!(k2.len(), 0, "score 20 < 25 → not rescued, dropped");
        assert_eq!(d2, 1);
    }

    // End-to-end emitted-gene totals through the real consensus writer: one
    // Augustus+SNAP-corroborated locus, one Augustus-only locus, one SNAP-only
    // locus. Default keeps all three; support ≥ 2 keeps only the corroborated one.
    #[test]
    fn merge_filter_drops_single_predictor_loci_by_count() {
        let dir = tempfile::tempdir().unwrap();

        let aug = dir.path().join("augustus.gff3");
        std::fs::write(
            &aug,
            "##gff-version 3\n\
chr1\tAUGUSTUS\tgene\t100\t200\t.\t+\t.\tID=a1\n\
chr1\tAUGUSTUS\tmRNA\t100\t200\t.\t+\t.\tID=a1.t1;Parent=a1\n\
chr1\tAUGUSTUS\tCDS\t100\t200\t.\t+\t0\tID=a1.cds;Parent=a1.t1\n\
chr1\tAUGUSTUS\tgene\t1000\t1100\t.\t+\t.\tID=a2\n\
chr1\tAUGUSTUS\tmRNA\t1000\t1100\t.\t+\t.\tID=a2.t1;Parent=a2\n\
chr1\tAUGUSTUS\tCDS\t1000\t1100\t.\t+\t0\tID=a2.cds;Parent=a2.t1\n",
        )
        .unwrap();

        let snap = dir.path().join("snap.gff3");
        std::fs::write(
            &snap,
            "##gff-version 3\n\
chr1\tSNAP\tgene\t100\t200\t.\t+\t.\tID=s1\n\
chr1\tSNAP\tmRNA\t100\t200\t.\t+\t.\tID=s1.t1;Parent=s1\n\
chr1\tSNAP\tCDS\t100\t200\t.\t+\t0\tID=s1.cds;Parent=s1.t1\n\
chr1\tSNAP\tgene\t5000\t5100\t.\t+\t.\tID=s2\n\
chr1\tSNAP\tmRNA\t5000\t5100\t.\t+\t.\tID=s2.t1;Parent=s2\n\
chr1\tSNAP\tCDS\t5000\t5100\t.\t+\t0\tID=s2.cds;Parent=s2.t1\n",
        )
        .unwrap();

        let inputs: Vec<(&Path, &str, f64)> = vec![
            (aug.as_path(), "augustus", 10.0),
            (snap.as_path(), "snap", 3.0),
        ];

        // Default: a1/s1 merge, a2 and s2 pass through → 3 genes.
        let out_def = dir.path().join("def.gff3");
        let n_def =
            merge_predictions_filtered(&inputs, &out_def, "D", &ConsensusFilter::default())
                .unwrap();
        assert_eq!(n_def, 3, "default keeps all three loci");
        assert_eq!(
            count_feature(&std::fs::read_to_string(&out_def).unwrap(), "gene"),
            3
        );

        // support ≥ 2: only the Augustus+SNAP locus (100-200) survives.
        let out_f = dir.path().join("filt.gff3");
        let filter = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: None,
        };
        let n_f = merge_predictions_filtered(&inputs, &out_f, "F", &filter).unwrap();
        assert_eq!(n_f, 1, "single-predictor loci dropped at support ≥ 2");

        let gff = std::fs::read_to_string(&out_f).unwrap();
        assert_eq!(count_feature(&gff, "gene"), 1);
        assert!(
            gff.lines().any(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                f.len() == 9 && f[2] == "gene" && f[3] == "100" && f[4] == "200"
            }),
            "surviving gene must be the corroborated 100-200 locus\n{}",
            gff
        );
    }
}
