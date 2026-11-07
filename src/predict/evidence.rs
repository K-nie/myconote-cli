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
use crate::parser::fasta::FastaRecord;
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

// Serde fallbacks for a weights TOML that omits a field. They match the
// current (corroborated-combiner) defaults below.
fn default_glimmer_weight() -> f64 {
    2.0
}
fn default_genemark_weight() -> f64 {
    3.0
}
fn default_miniprot_weight() -> f64 {
    8.0
}
fn default_trna_weight() -> f64 {
    6.0
}

/// Default weights (v0.9.1+), tuned for the corroborated combiner.
///
/// The old 10 / 5 / 3 (Augustus / GeneMark / SNAP) split made Augustus
/// structurally decisive: with Augustus wrong and SNAP + GeneMark agreeing on
/// the correct model, the agreement score is `10 + 8j` for Augustus vs
/// `8 + 10j` for the pair, so Augustus wins for every `j < 1`. The Nov 2026
/// accuracy diagnostic (docs/paper/htcondor_benchmark/diagnostics) found that
/// 60–91% of the loci funannotate gets exactly right and MycoNote does not were
/// already called correctly by SNAP and/or GeneMark and then outvoted this way.
///
/// The new ratios keep the ordering Augustus > GeneMark > SNAP (so Augustus
/// still breaks three-way disagreements) but satisfy `snap + genemark >
/// augustus`, so two independent predictors that agree exactly can overrule a
/// lone Augustus call. The largest ab-initio ratio is 2:1, the same ceiling
/// funannotate's EVM configuration uses. Extrinsic tracks are scaled by the
/// same 0.4 factor as Augustus so their standing relative to Augustus is
/// unchanged.
impl Default for EvidenceWeights {
    fn default() -> Self {
        Self {
            augustus: 4.0,
            snap: 2.0,
            protein: 8.0,
            est: 3.0,
            glimmerhmm: 2.0,
            genemark: 3.0,
            miniprot: 8.0,
            trnascan: 6.0,
        }
    }
}

impl EvidenceWeights {
    /// The pre-v0.9.1 weights (Augustus 10, GeneMark 5, SNAP 3, ...). Used by
    /// `--combiner legacy` so old benchmark numbers can be reproduced exactly.
    pub fn legacy() -> Self {
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

    /// Load weights from a TOML configuration file.
    ///
    /// Example TOML:
    /// ```toml
    /// augustus = 4.0
    /// snap = 2.0
    /// protein = 8.0
    /// est = 3.0
    /// glimmerhmm = 2.0
    /// genemark = 3.0
    /// miniprot = 8.0
    /// trnascan = 6.0
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
    /// Minimum distinct-predictor support required specifically for a
    /// **single-exon, purely-ab-initio** coding locus (BRAKER3's "remove
    /// single-exon genes predicted purely ab initio" rule). Single-exon models
    /// called by one ab-initio predictor alone are the dominant false-positive
    /// source; multi-exon models are trustworthy and stay governed by
    /// `min_predictor_support`. `1` (default) applies no extra requirement, so
    /// behaviour is unchanged. The effective bar for such a locus is
    /// `max(min_predictor_support, min_support_single_exon)`, so this is never
    /// weaker than the general threshold. A locus with any extrinsic
    /// (protein / EST / miniprot) support, or any non-coding locus (no CDS,
    /// e.g. a tRNA), is exempt and uses `min_predictor_support`.
    pub min_support_single_exon: usize,
    /// When true, drop coding models that fail a genome-free ORF-completeness
    /// check (concatenated CDS length not a multiple of 3). `false` (default)
    /// keeps everything. This is a predict-stage structural gate only; it does
    /// NOT verify start/stop nucleotides or internal stop codons, which need
    /// the genome sequence and are handled at the annotate stage.
    pub require_complete_orf: bool,
}

impl Default for ConsensusFilter {
    fn default() -> Self {
        Self {
            min_predictor_support: 1,
            min_consensus_score: None,
            min_support_single_exon: 1,
            require_complete_orf: false,
        }
    }
}

impl ConsensusFilter {
    /// Decide whether a winning model is kept given its support count (distinct
    /// corroborating sources, incl. itself), its consensus score, and whether
    /// it is a single-exon purely-ab-initio coding locus (the only case to
    /// which `min_support_single_exon` applies).
    ///
    /// The ORF-completeness gate is applied separately by the caller (it needs
    /// the model's CDS records, not just these scalars) and is independent of
    /// the score rescue: an incomplete ORF is dropped regardless of score when
    /// `require_complete_orf` is set.
    fn keeps(&self, support: usize, consensus_score: f64, single_exon_ab_initio: bool) -> bool {
        let needed = if single_exon_ab_initio {
            self.min_predictor_support.max(self.min_support_single_exon)
        } else {
            self.min_predictor_support
        };
        if support >= needed {
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
            && self.min_support_single_exon <= 1
            && !self.require_complete_orf
    }
}

/// Ab-initio protein-coding predictors, for the single-exon rule. A locus whose
/// support comes only from these (no protein / EST / miniprot evidence) is
/// "purely ab initio". tRNAscan is deliberately excluded: its genes are
/// single-exon by nature and must never be caught by the single-exon drop.
fn is_ab_initio_source(source: &str) -> bool {
    matches!(
        source.to_lowercase().as_str(),
        "augustus" | "snap" | "genemark" | "genemark-es" | "genemark-et" | "glimmerhmm" | "glimmer"
    )
}

/// Count a model's CDS records. 0 means a non-coding model (e.g. tRNA), which
/// the single-exon and ORF checks treat as exempt.
fn cds_record_count(m: &GeneModel) -> usize {
    m.records.iter().filter(|r| r.feature_type == "CDS").count()
}

/// Genome-free ORF-completeness heuristic: a complete coding sequence spans
/// whole codons, so its concatenated CDS length (1-based inclusive) is a
/// multiple of 3. A model with no CDS is non-coding and exempt (returns true).
/// This does not check the ATG / stop nucleotides or internal stops — those
/// need the FASTA and are verified at the annotate stage.
fn orf_length_complete(m: &GeneModel) -> bool {
    let coding: Vec<u64> = m
        .records
        .iter()
        .filter(|r| r.feature_type == "CDS")
        .map(|r| r.end.saturating_sub(r.start) + 1)
        .collect();
    if coding.is_empty() {
        return true;
    }
    coding.iter().sum::<u64>() % 3 == 0
}

// ─────────────────────────────────────────────────────────────────────────────
// Verified ORF-completeness check (v0.9.0, Feature #4)
//
// Extends the genome-free frame gate (`orf_length_complete`) with the three
// nucleotide checks it could not make without the FASTA: the start codon, the
// terminal stop codon, and the absence of an internal stop. When a genome index
// is supplied the consensus filter runs the full check; without one it degrades
// to the length-only gate (so the flag still does something in a genome-less
// unit test). The check is opt-in via `--require-complete-orf` and never alters
// the default output.
// ─────────────────────────────────────────────────────────────────────────────

/// Why a coding model failed (or passed) the ORF-completeness check. Used to
/// give the dropped-locus log a precise, actionable reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OrfStatus {
    Complete,
    /// Concatenated CDS length is not a multiple of 3 — broken frame.
    BadLength,
    /// First codon is not a start codon (ATG).
    NoStart,
    /// Final codon is not a stop codon for the active translation table.
    NoStop,
    /// A stop codon appears before the terminal codon (truncating protein).
    InternalStop,
    /// The genome lacks this model's sequence, so the nucleotide checks could
    /// not run. Treated as "keep" — we never drop a model over a tool-side gap.
    Unverifiable,
}

impl OrfStatus {
    fn reason(self) -> &'static str {
        match self {
            OrfStatus::Complete => "complete ORF",
            OrfStatus::BadLength => "incomplete ORF (CDS length not a multiple of 3)",
            OrfStatus::NoStart => "incomplete ORF (no ATG start codon)",
            OrfStatus::NoStop => "incomplete ORF (no terminal stop codon)",
            OrfStatus::InternalStop => "incomplete ORF (internal stop codon)",
            OrfStatus::Unverifiable => "ORF unverifiable (sequence not found in genome)",
        }
    }
}

/// Stop codons for an NCBI translation table. Tables 1 (standard) and 12
/// (the *Candida* CTG clade, the common fungal exceptions) share the standard
/// stop set; only internal-codon reassignment differs there, not the stops.
/// Table 6 (some ciliates — rare in fungi) reassigns TAA/TAG to Gln, leaving
/// TGA as the sole stop. Any other table falls back to the standard set.
fn stop_codons_for(code: u8) -> &'static [&'static str] {
    match code {
        6 => &["TGA"],
        _ => &["TAA", "TAG", "TGA"],
    }
}

/// Reconstruct a model's spliced coding sequence (5'→3') from the genome.
/// CDS intervals are concatenated in ascending genomic order; a minus-strand
/// model is then reverse-complemented, which yields the exons in coding order
/// each complemented — exactly the mRNA coding sequence. Returns `None` when
/// the model has no CDS, or the genome has no record for its sequence.
fn spliced_cds(m: &GeneModel, genome: &HashMap<String, FastaRecord>) -> Option<String> {
    let rec = genome.get(&m.seqid)?;
    let mut intervals: Vec<(u64, u64)> = m
        .records
        .iter()
        .filter(|r| r.feature_type == "CDS")
        .map(|r| (r.start, r.end))
        .collect();
    if intervals.is_empty() {
        return None;
    }
    intervals.sort_by_key(|&(s, _)| s);
    let mut seq = String::new();
    for (s, e) in &intervals {
        seq.push_str(rec.subsequence(*s, *e));
    }
    if m.strand == '-' {
        Some(crate::parser::fasta::reverse_complement(&seq))
    } else {
        Some(seq)
    }
}

/// Classify a model's ORF completeness. A non-coding model (0 CDS) is exempt
/// (`Complete`). With no genome this reduces to the length-only frame gate;
/// with a genome it adds the start / stop / internal-stop nucleotide checks.
fn orf_status(
    m: &GeneModel,
    genome: Option<&HashMap<String, FastaRecord>>,
    genetic_code: u8,
) -> OrfStatus {
    if cds_record_count(m) == 0 {
        return OrfStatus::Complete; // non-coding (e.g. tRNA) — not an ORF
    }
    if !orf_length_complete(m) {
        return OrfStatus::BadLength;
    }
    let genome = match genome {
        Some(g) => g,
        None => return OrfStatus::Complete, // length-only mode (frame already OK)
    };
    let cds = match spliced_cds(m, genome) {
        Some(s) => s,
        None => return OrfStatus::Unverifiable,
    };
    if cds.len() < 3 || cds.len() % 3 != 0 {
        return OrfStatus::BadLength;
    }
    if &cds[0..3] != "ATG" {
        return OrfStatus::NoStart;
    }
    let stops = stop_codons_for(genetic_code);
    let last = &cds[cds.len() - 3..];
    if !stops.contains(&last) {
        return OrfStatus::NoStop;
    }
    let n_codons = cds.len() / 3;
    for i in 0..(n_codons - 1) {
        let codon = &cds[i * 3..i * 3 + 3];
        if stops.contains(&codon) {
            return OrfStatus::InternalStop;
        }
    }
    OrfStatus::Complete
}

// ─────────────────────────────────────────────────────────────────────────────
// Post-consensus passes: fused-gene splitting (Feature #2) and repeat filtering
// (Feature #3). Both run AFTER overlap resolution, on the winning models, and
// both are opt-in — a default run never constructs them.
// ─────────────────────────────────────────────────────────────────────────────

/// Extra, opt-in knobs applied to the consensus gene set after overlap
/// resolution. All fields default to "off / unchanged", so the default merge
/// path is byte-identical to earlier releases.
#[derive(Debug, Clone)]
pub struct ConsensusExtras {
    /// Split fused gene models at implausibly long inter-CDS gaps
    /// (funannotate's `--jaccard_clip` analog; see `split_fused_model`).
    pub jaccard_clip: bool,
    /// Maximum plausible intron (bp). An inter-CDS gap larger than this cannot
    /// be a real intron, so under `--jaccard-clip` the model is split there.
    pub jaccard_clip_max_intron: usize,
    /// Minimum coding length (bp) each half must keep for a split to proceed.
    /// A model that would fragment into a sub-`min` piece is left intact.
    pub jaccard_clip_min_cds: u64,
    /// Drop gene models whose CDS lies mostly in soft-masked repeat.
    pub repeat_filter: bool,
    /// Per-gene CDS-in-repeat fraction above which a model is dropped. `None`
    /// means pick it adaptively from the genome-wide soft-mask fraction.
    pub repeat_filter_threshold: Option<f64>,
    /// NCBI translation table for the verified-ORF nucleotide check.
    pub genetic_code: u8,
    /// How each overlap group is resolved (see [`Combiner`]).
    pub combiner: Combiner,
}

/// Overlap-group resolution strategy for the consensus merge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Combiner {
    /// v0.9.1+ default. Two changes over `Legacy`:
    ///
    /// 1. **Per-source agreement.** Each other source contributes at most one
    ///    vote to a candidate, `weight · max_jaccard` over that source's models,
    ///    instead of one vote per overlapping model. Under `Legacy` a fused
    ///    model spanning two genes collected partial credit from *every* piece
    ///    the other predictors called, so it outscored the correctly split
    ///    calls; this is a driver of the "merged ≫ split" error.
    /// 2. **Several winners per group.** After the best model is kept, the
    ///    next-best model that does not overlap any kept model is kept too, and
    ///    so on. `Legacy` keeps exactly one model per transitively-overlapping
    ///    group, so when a correct non-fused call wins, its neighbour (linked
    ///    to the group only through the fused model) is silently dropped.
    #[default]
    Corroborated,
    /// Pre-v0.9.1 behaviour: sum agreement over every overlapping model and
    /// keep exactly one winner per overlap group.
    Legacy,
}

impl Combiner {
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "corroborated" | "default" | "v2" => Some(Self::Corroborated),
            "legacy" | "v1" => Some(Self::Legacy),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Corroborated => "corroborated",
            Self::Legacy => "legacy",
        }
    }
}

impl Default for ConsensusExtras {
    fn default() -> Self {
        Self {
            jaccard_clip: false,
            jaccard_clip_max_intron: 10_000,
            jaccard_clip_min_cds: 150,
            repeat_filter: false,
            repeat_filter_threshold: None,
            genetic_code: 1,
            combiner: Combiner::default(),
        }
    }
}

/// Everything the consensus merge needs beyond the GFF inputs and the support
/// `ConsensusFilter`: an optional genome index (for the verified-ORF check),
/// an optional soft-mask map (for repeat filtering), and the extra knobs.
/// `Default` is all-off, which keeps the merge byte-identical to prior output.
pub struct ConsensusContext<'a> {
    pub genome: Option<&'a HashMap<String, FastaRecord>>,
    pub repeat_mask: Option<&'a RepeatMask>,
    pub extras: ConsensusExtras,
}

impl Default for ConsensusContext<'_> {
    fn default() -> Self {
        Self {
            genome: None,
            repeat_mask: None,
            extras: ConsensusExtras::default(),
        }
    }
}

/// Split a fused gene model wherever two consecutive CDS exons are separated by
/// a gap longer than the maximum plausible intron. A gap that large cannot be a
/// real intron, so the model almost certainly fuses two adjacent genes — the
/// dominant over-merge error funannotate's `--jaccard_clip` targets in
/// gene-dense fungal genomes.
///
/// MycoNote differs from funannotate here, deliberately and loudly: funannotate
/// clips with PASA *transcript-alignment* Jaccard, which needs assembled RNA
/// evidence. The consensus stage has no transcripts, so MycoNote uses the
/// structural signal it does have — an inter-CDS gap beyond the biological
/// intron ceiling — and confirms the two halves are disjoint (coding-base
/// Jaccard 0 by construction) and each still carries `min_cds` coding bases.
/// The split is conservative: if *any* resulting piece would fall below
/// `min_cds`, the model is left intact rather than fragmented.
fn split_fused_model(m: &GeneModel, max_intron: usize, min_cds: u64) -> Vec<GeneModel> {
    let mut cds: Vec<&GFFRecord> = m
        .records
        .iter()
        .filter(|r| r.feature_type == "CDS")
        .collect();
    if cds.len() < 2 {
        return vec![m.clone()]; // single-exon: no inter-CDS gap to split on
    }
    cds.sort_by_key(|r| r.start);

    // Partition CDS indices into contiguous runs, breaking at over-long gaps.
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut cur: Vec<usize> = vec![0];
    for i in 1..cds.len() {
        let gap = cds[i]
            .start
            .saturating_sub(cds[i - 1].end)
            .saturating_sub(1);
        if gap > max_intron as u64 {
            groups.push(std::mem::take(&mut cur));
        }
        cur.push(i);
    }
    groups.push(cur);

    if groups.len() < 2 {
        return vec![m.clone()]; // no over-long gap → nothing to split
    }

    let group_cds_len =
        |g: &[usize]| -> u64 { g.iter().map(|&i| cds[i].end - cds[i].start + 1).sum() };

    // Conservative guard: never fragment a model into a sub-`min_cds` piece.
    if groups.iter().any(|g| group_cds_len(g) < min_cds) {
        return vec![m.clone()];
    }

    let gene_rec = m.records.iter().find(|r| r.feature_type == "gene").cloned();
    let tx_rec = m
        .records
        .iter()
        .find(|r| matches!(r.feature_type.as_str(), "mRNA" | "transcript"))
        .cloned();

    let mut out: Vec<GeneModel> = Vec::with_capacity(groups.len());
    for g in &groups {
        let gstart = g.iter().map(|&i| cds[i].start).min().unwrap();
        let gend = g.iter().map(|&i| cds[i].end).max().unwrap();

        let mut recs: Vec<GFFRecord> = Vec::new();
        if let Some(mut gr) = gene_rec.clone() {
            gr.start = gstart;
            gr.end = gend;
            recs.push(gr);
        }
        if let Some(mut tr) = tx_rec.clone() {
            tr.start = gstart;
            tr.end = gend;
            recs.push(tr);
        }
        for &i in g {
            recs.push(cds[i].clone());
        }

        out.push(GeneModel {
            seqid: m.seqid.clone(),
            start: gstart,
            end: gend,
            strand: m.strand,
            score: m.score,
            source: m.source.clone(),
            records: recs,
        });
    }
    out
}

/// Apply `split_fused_model` across a set of winners. Returns the (possibly
/// larger) model list and the number of *extra* loci created by splitting.
fn apply_jaccard_clip(
    models: Vec<GeneModel>,
    max_intron: usize,
    min_cds: u64,
) -> (Vec<GeneModel>, usize) {
    let mut out: Vec<GeneModel> = Vec::with_capacity(models.len());
    let mut extra = 0usize;
    for m in models {
        let parts = split_fused_model(&m, max_intron, min_cds);
        if parts.len() > 1 {
            extra += parts.len() - 1;
            eprintln!(
                "  ⚠  jaccard-clip split fused model {}:{}-{} ({}) into {} genes \
                 (inter-CDS gap > {} bp)",
                m.seqid,
                m.start,
                m.end,
                m.source,
                parts.len(),
                max_intron
            );
        }
        out.extend(parts);
    }
    (out, extra)
}

/// Soft-mask (lowercase) intervals read from a masked FASTA with case
/// **preserved**. The predict input is soft-masked (repeats lowercased), so
/// lowercase runs are the repeat annotation. Built only when `--repeat-filter`
/// is on. The streaming reader keeps memory to the per-sequence run list, not
/// the whole genome.
pub struct RepeatMask {
    /// seqid → sorted, 1-based inclusive lowercase runs.
    intervals: HashMap<String, Vec<(u64, u64)>>,
    /// Genome-wide fraction of bases that are soft-masked (lowercase).
    pub masked_fraction: f64,
}

impl RepeatMask {
    /// Scan a masked FASTA, preserving case, and record every lowercase run.
    pub fn from_fasta(path: &Path) -> Result<Self> {
        use std::io::{BufRead, BufReader};
        let f = std::fs::File::open(path).map_err(crate::utils::error::MycoNoteError::Io)?;
        let reader = BufReader::new(f);

        let mut intervals: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut cur_id: Option<String> = None;
        let mut pos: u64 = 0;
        let mut run_start: Option<u64> = None;
        let mut total: u64 = 0;
        let mut masked: u64 = 0;

        for line in reader.lines() {
            let line = line.map_err(crate::utils::error::MycoNoteError::Io)?;
            let t = line.trim_end();
            if t.is_empty() {
                continue;
            }
            if let Some(stripped) = t.strip_prefix('>') {
                // Close any run open at the end of the previous sequence.
                if let (Some(id), Some(rs)) = (&cur_id, run_start.take()) {
                    intervals.entry(id.clone()).or_default().push((rs, pos));
                }
                let id = stripped.split_whitespace().next().unwrap_or("").to_string();
                cur_id = Some(id);
                pos = 0;
                run_start = None;
                continue;
            }
            for ch in t.chars() {
                if ch.is_whitespace() {
                    continue;
                }
                pos += 1;
                total += 1;
                if ch.is_ascii_lowercase() {
                    masked += 1;
                    if run_start.is_none() {
                        run_start = Some(pos);
                    }
                } else if let Some(rs) = run_start.take() {
                    if let Some(id) = &cur_id {
                        intervals.entry(id.clone()).or_default().push((rs, pos - 1));
                    }
                }
            }
        }
        // Close a run that reaches the end of the final sequence.
        if let (Some(id), Some(rs)) = (&cur_id, run_start.take()) {
            intervals.entry(id.clone()).or_default().push((rs, pos));
        }

        let masked_fraction = if total == 0 {
            0.0
        } else {
            masked as f64 / total as f64
        };
        Ok(RepeatMask {
            intervals,
            masked_fraction,
        })
    }

    /// Fraction of a model's CDS bases that fall within soft-masked repeat.
    fn cds_repeat_fraction(&self, m: &GeneModel) -> f64 {
        let runs = match self.intervals.get(&m.seqid) {
            Some(r) => r,
            None => return 0.0,
        };
        let mut cds: Vec<(u64, u64)> = m
            .records
            .iter()
            .filter(|r| r.feature_type == "CDS")
            .map(|r| (r.start, r.end))
            .collect();
        if cds.is_empty() {
            return 0.0;
        }
        cds.sort_by_key(|&(s, _)| s);
        let total: u64 = cds.iter().map(|&(s, e)| e - s + 1).sum();
        if total == 0 {
            return 0.0;
        }
        let overlap = interval_overlap_bp(&cds, runs);
        overlap as f64 / total as f64
    }

    /// Adaptive per-gene drop threshold, chosen from the genome-wide masked
    /// fraction when the user did not pin one with `--repeat-filter-threshold`.
    ///
    /// Compact, repeat-poor genomes are treated leniently — only a model whose
    /// CDS is almost entirely repeat (≥ 90 %) is dropped. Repeat-rich genomes
    /// trigger aggressive filtering down to a 50 % CDS-in-repeat threshold. The
    /// threshold ramps linearly between 10 % and 40 % genome-wide soft-mask
    /// fraction (funannotate flips `--repeat_filter` on by hand; MycoNote scales
    /// it so a compact genome is not over-filtered while a transposon-rich one
    /// is).
    pub fn adaptive_threshold(&self) -> f64 {
        let (lo, hi) = (0.10_f64, 0.40_f64);
        let f = self.masked_fraction.clamp(lo, hi);
        let t = (f - lo) / (hi - lo); // 0 at lo, 1 at hi
        0.90 - t * (0.90 - 0.50) // 0.90 (lenient) → 0.50 (aggressive)
    }
}

/// Drop gene models whose CDS-in-repeat fraction exceeds `threshold`. Returns
/// the kept models and the number dropped. Non-coding models (0 CDS) are never
/// dropped. Every drop is logged by name and reason (no silent data loss).
fn apply_repeat_filter(
    models: Vec<GeneModel>,
    mask: &RepeatMask,
    threshold: f64,
) -> (Vec<GeneModel>, usize) {
    let mut kept: Vec<GeneModel> = Vec::with_capacity(models.len());
    let mut dropped = 0usize;
    for m in models {
        if cds_record_count(&m) == 0 {
            kept.push(m);
            continue;
        }
        let frac = mask.cds_repeat_fraction(&m);
        if frac > threshold {
            dropped += 1;
            eprintln!(
                "  ⚠  repeat filter dropped {}:{}-{} ({}) — {:.0}% of CDS in soft-masked \
                 repeat (> {:.0}% threshold)",
                m.seqid,
                m.start,
                m.end,
                m.source,
                frac * 100.0,
                threshold * 100.0,
            );
        } else {
            kept.push(m);
        }
    }
    (kept, dropped)
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

/// CDS intervals with the stop codon folded into the terminal segment, for
/// agreement scoring. Predictors disagree on the convention: Augustus and
/// GeneMark emit CDS without the stop (plus a separate `stop_codon` row), SNAP
/// emits CDS that already includes it. Compared raw, an Augustus and a SNAP
/// model with identical structure never reach Jaccard 1 (on the bundled
/// Brettanomyces run 0 of 5218 Augustus models matched SNAP exactly raw vs
/// 3439 after folding). The write path already folds stops, so this only
/// aligns scoring with what is emitted.
fn cds_intervals_stop_folded(m: &GeneModel) -> Vec<(u64, u64)> {
    if !m.records.iter().any(|r| r.feature_type == "stop_codon") {
        return cds_intervals(m);
    }
    let mut recs: Vec<GFFRecord> = m
        .records
        .iter()
        .filter(|r| r.feature_type == "CDS" || r.feature_type == "stop_codon")
        .cloned()
        .collect();
    crate::parser::gff::include_stop_codon_in_cds(&mut recs);
    let folded = GeneModel {
        seqid: String::new(),
        start: m.start,
        end: m.end,
        strand: m.strand,
        score: m.score,
        source: String::new(),
        records: recs,
    };
    cds_intervals(&folded)
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
    resolve_overlaps_filtered(
        models,
        &ConsensusFilter::default(),
        None,
        1,
        Combiner::default(),
    )
    .0
}

/// Like `resolve_overlaps` but applies a `ConsensusFilter` to each group's
/// winner. Returns `(kept, dropped)` where `dropped` counts loci discarded by
/// the filter. With a no-op filter this is identical to `resolve_overlaps` and
/// `dropped == 0`.
fn resolve_overlaps_filtered(
    mut models: Vec<GeneModel>,
    filter: &ConsensusFilter,
    genome: Option<&HashMap<String, FastaRecord>>,
    genetic_code: u8,
    combiner: Combiner,
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
        // agreement-weighted consensus.
        let cds: Vec<Vec<(u64, u64)>> = group
            .iter()
            .map(|&idx| match combiner {
                Combiner::Legacy => cds_intervals(&models[idx]),
                Combiner::Corroborated => cds_intervals_stop_folded(&models[idx]),
            })
            .collect();

        let scores: Vec<f64> = (0..group.len())
            .map(|x| match combiner {
                Combiner::Legacy => {
                    let mut s = models[group[x]].score;
                    for y in 0..group.len() {
                        if y == x {
                            continue;
                        }
                        s += models[group[y]].score * cds_jaccard(&cds[x], &cds[y]);
                    }
                    s
                }
                Combiner::Corroborated => {
                    // One vote per source: that source's best-agreeing model.
                    let mut best_by_source: HashMap<&str, f64> = HashMap::new();
                    for y in 0..group.len() {
                        if y == x {
                            continue;
                        }
                        let v = models[group[y]].score * cds_jaccard(&cds[x], &cds[y]);
                        let e = best_by_source
                            .entry(models[group[y]].source.as_str())
                            .or_insert(0.0);
                        if v > *e {
                            *e = v;
                        }
                    }
                    models[group[x]].score + best_by_source.values().sum::<f64>()
                }
            })
            .collect();

        let cds_len =
            |x: usize| -> u64 { cds[x].iter().map(|&(s, e)| e.saturating_sub(s) + 1).sum() };

        // Rank candidates best-first: consensus score, then raw source weight,
        // then longer CDS, then position (fully deterministic).
        let mut order: Vec<usize> = (0..group.len()).collect();
        order.sort_by(|&a, &b| {
            scores[b]
                .partial_cmp(&scores[a])
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(
                    models[group[b]]
                        .score
                        .partial_cmp(&models[group[a]].score)
                        .unwrap_or(std::cmp::Ordering::Equal),
                )
                .then(cds_len(b).cmp(&cds_len(a)))
                .then(models[group[a]].start.cmp(&models[group[b]].start))
                .then(models[group[a]].end.cmp(&models[group[b]].end))
        });

        // Legacy keeps one model per group. Corroborated keeps, greedily in
        // rank order, every model whose span overlaps no model already picked.
        let mut picked: Vec<usize> = Vec::new();
        for &cand in &order {
            let m = &models[group[cand]];
            if picked
                .iter()
                .any(|&p| models[group[p]].start <= m.end && m.start <= models[group[p]].end)
            {
                continue;
            }
            picked.push(cand);
            if combiner == Combiner::Legacy {
                break;
            }
        }

        for best in picked {
            // Fast path: the default no-op filter keeps every winner, so skip
            // the support count entirely.
            if filter.is_noop() {
                kept.push(models[group[best]].clone());
                continue;
            }

            let winner_score = scores[best];

            // Support = distinct predictor sources whose CDS overlaps the
            // winner, counting the winner's own source. Two models from the
            // *same* source do not inflate support (we want corroboration from
            // an independent predictor, which is what the benchmark FP source
            // lacked).
            let mut sources: std::collections::HashSet<&str> = std::collections::HashSet::new();
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

            // Single-exon purely-ab-initio determination (BRAKER3 rule). A
            // locus is subject to the stricter single-exon threshold only when
            // (a) it is protein-coding with exactly one CDS, and (b) every
            // corroborating source is an ab-initio predictor — any extrinsic
            // (protein / EST / miniprot) support exempts it, as does a
            // non-coding model (0 CDS, e.g. tRNA → not single-exon here).
            let cds_count = cds_record_count(&models[group[best]]);
            let purely_ab_initio = sources.iter().all(|s| is_ab_initio_source(s));
            let single_exon_ab_initio = cds_count == 1 && purely_ab_initio;

            // Two independent gates, both must pass to keep. The ORF gate is
            // independent of the score rescue (an incomplete ORF is structural,
            // not a confidence question).
            let keeps_support = filter.keeps(support, winner_score, single_exon_ab_initio);
            let orf = if filter.require_complete_orf {
                orf_status(&models[group[best]], genome, genetic_code)
            } else {
                OrfStatus::Complete
            };
            let orf_ok = orf == OrfStatus::Complete;

            if keeps_support && orf_ok {
                kept.push(models[group[best]].clone());
            } else {
                dropped += 1;
                // No silent data loss: name every locus the filter removes and
                // why, so non-TTY cluster logs record what the gene set lost.
                let reason = if !orf_ok {
                    orf.reason().to_string()
                } else {
                    let need = if single_exon_ab_initio {
                        filter
                            .min_predictor_support
                            .max(filter.min_support_single_exon)
                    } else {
                        filter.min_predictor_support
                    };
                    format!(
                        "support={} source(s){}, need ≥ {}{}",
                        support,
                        if single_exon_ab_initio {
                            " [single-exon ab-initio]"
                        } else {
                            ""
                        },
                        need,
                        match filter.min_consensus_score {
                            Some(m) => format!(" or score ≥ {:.2}", m),
                            None => String::new(),
                        }
                    )
                };
                eprintln!(
                    "  ⚠  consensus filter dropped {}:{}-{} ({}) — score={:.2}; {}",
                    models[group[best]].seqid,
                    models[group[best]].start,
                    models[group[best]].end,
                    models[group[best]].source,
                    winner_score,
                    reason,
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
    merge_predictions_filtered(
        inputs,
        output_gff,
        locus_prefix,
        &ConsensusFilter::default(),
        &ConsensusContext::default(),
    )
}

/// Merge predictions into a consensus GFF3, dropping false-positive loci that
/// do not meet `filter`. With `ConsensusFilter::default()` this is identical to
/// `merge_predictions`.
pub fn merge_predictions_filtered(
    inputs: &[(&Path, &str, f64)],
    output_gff: &Path,
    locus_prefix: &str,
    filter: &ConsensusFilter,
    ctx: &ConsensusContext,
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
        let (kept, dropped) = resolve_overlaps_filtered(
            models,
            filter,
            ctx.genome,
            ctx.extras.genetic_code,
            ctx.extras.combiner,
        );
        total_dropped += dropped;
        final_models.extend(kept);
    }

    // ── Post-consensus pass A: fused-gene splitting (--jaccard-clip) ──────────
    // Runs after overlap resolution so it operates on the final winners. Off by
    // default; when on it only ever ADDS loci (splitting one into several).
    if ctx.extras.jaccard_clip {
        let (split_models, extra) = apply_jaccard_clip(
            final_models,
            ctx.extras.jaccard_clip_max_intron,
            ctx.extras.jaccard_clip_min_cds,
        );
        final_models = split_models;
        if extra > 0 {
            println!(
                "  Jaccard-clip: split fused models into {} extra locus/loci \
                 (inter-CDS gap > {} bp)",
                extra, ctx.extras.jaccard_clip_max_intron
            );
        }
    }

    // ── Post-consensus pass B: adaptive repeat filtering (--repeat-filter) ────
    // Needs the soft-mask map. The per-gene threshold is the user's pin, or the
    // adaptive value derived from the genome-wide soft-mask fraction.
    let mut repeat_dropped = 0usize;
    if ctx.extras.repeat_filter {
        if let Some(mask) = ctx.repeat_mask {
            let threshold = ctx
                .extras
                .repeat_filter_threshold
                .unwrap_or_else(|| mask.adaptive_threshold());
            println!(
                "  Repeat filter: genome soft-mask fraction {:.1}% → per-gene drop \
                 threshold {:.0}% CDS-in-repeat{}",
                mask.masked_fraction * 100.0,
                threshold * 100.0,
                if ctx.extras.repeat_filter_threshold.is_some() {
                    " (user-pinned)"
                } else {
                    " (adaptive)"
                },
            );
            let (kept, dropped) = apply_repeat_filter(final_models, mask, threshold);
            final_models = kept;
            repeat_dropped = dropped;
            if dropped > 0 {
                println!(
                    "  Repeat filter: dropped {} repeat-overlapping locus/loci",
                    dropped
                );
            }
        } else {
            eprintln!(
                "  ⚠  --repeat-filter requested but no soft-mask map was built \
                 (needs the masked FASTA) — skipping repeat filtering."
            );
        }
    }
    let _ = repeat_dropped;

    // Sort by chromosome then position
    final_models.sort_by(|a, b| a.seqid.cmp(&b.seqid).then(a.start.cmp(&b.start)));

    if total_dropped > 0 {
        println!(
            "  Consensus filter: dropped {} false-positive locus/loci \
             (min-predictor-support={}, min-single-exon-support={}{}{})",
            total_dropped,
            filter.min_predictor_support,
            filter.min_support_single_exon,
            match filter.min_consensus_score {
                Some(m) => format!(", min-consensus-score={:.2}", m),
                None => String::new(),
            },
            if filter.require_complete_orf {
                ", require-complete-orf"
            } else {
                ""
            },
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

        // A GFF3 parent feature MUST span all its children. Some predictors
        // (and the consensus merge, which can combine CDS from more than one
        // predictor) emit a `gene` / `mRNA` whose own start/end cover only a
        // subset of the model's CDS — e.g. a single exon — leaving the other
        // CDS "outside mRNA bounds". That is invalid GFF and makes table2asn /
        // NCBI reject the record. Recompute the true span (min child start,
        // max child end) over every coordinate-bearing child and clamp the
        // parent features to it when we emit them below. Codons are already
        // folded into the CDS, so they never extend this span.
        let (child_min, child_max) = norm_records
            .iter()
            .filter(|r| !matches!(r.feature_type.as_str(), "gene" | "mRNA" | "transcript"))
            .fold((u64::MAX, 0u64), |(mn, mx), r| {
                (mn.min(r.start), mx.max(r.end))
            });
        let have_children = child_min <= child_max;

        for rec in &norm_records {
            let mut r = rec.clone();

            match rec.feature_type.as_str() {
                "gene" => {
                    r.attributes.insert("ID".into(), gene_id.clone());
                    r.attributes.insert("locus_tag".into(), locus_tag.clone());
                    // Clamp the gene to span all its children (see above).
                    if have_children {
                        r.start = child_min;
                        r.end = child_max;
                    }
                }
                "mRNA" | "transcript" => {
                    // F3b: normalize whichever label the predictor used.
                    r.feature_type = "mRNA".to_string();
                    r.attributes.insert("ID".into(), mrna_id.clone());
                    r.attributes.insert("Parent".into(), gene_id.clone());
                    mrna_source = Some(r.source.clone());
                    // Clamp the mRNA to span all its children (see above).
                    if have_children {
                        r.start = child_min;
                        r.end = child_max;
                    }
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

    // Diagnostic case (Nov 2026): Augustus calls a wrong structure, SNAP and
    // GeneMark agree exactly on the correct one. Under the legacy 10/5/3
    // weights Augustus always wins (`10 + 8j > 8 + 10j` for j < 1); under the
    // default weights the corroborated pair wins.
    fn augustus_vs_agreeing_pair(w: &EvidenceWeights, c: Combiner) -> Vec<GeneModel> {
        // Wrong: alternative 3' exon boundary. Correct: two exons.
        let aug = gene_model_src(100, 900, w.augustus, &[(100, 300), (400, 900)], "Augustus");
        let snap = gene_model_src(100, 800, w.snap, &[(100, 300), (400, 800)], "SNAP");
        let gm = gene_model_src(100, 800, w.genemark, &[(100, 300), (400, 800)], "GeneMark");
        resolve_overlaps_filtered(vec![aug, snap, gm], &ConsensusFilter::default(), None, 1, c).0
    }

    #[test]
    fn agreeing_pair_outvotes_lone_augustus_by_default() {
        let kept = augustus_vs_agreeing_pair(&EvidenceWeights::default(), Combiner::default());
        assert_eq!(kept.len(), 1);
        assert_eq!((kept[0].start, kept[0].end), (100, 800));
    }

    #[test]
    fn legacy_weights_and_combiner_reproduce_augustus_dominance() {
        let kept = augustus_vs_agreeing_pair(&EvidenceWeights::legacy(), Combiner::Legacy);
        assert_eq!(kept.len(), 1);
        assert_eq!((kept[0].start, kept[0].end), (100, 900));
        assert_eq!(kept[0].source, "Augustus");
    }

    // With no corroboration either way, Augustus still wins a two-way
    // disagreement against SNAP or GeneMark alone.
    #[test]
    fn default_weights_keep_augustus_on_two_way_disagreement() {
        let w = EvidenceWeights::default();
        for (other, ow) in [("SNAP", w.snap), ("GeneMark", w.genemark)] {
            let aug = gene_model_src(100, 900, w.augustus, &[(100, 300), (400, 900)], "Augustus");
            let o = gene_model_src(100, 800, ow, &[(100, 300), (400, 800)], other);
            let kept = resolve_overlaps(vec![aug, o]);
            assert_eq!(kept.len(), 1);
            assert_eq!(kept[0].source, "Augustus", "vs {}", other);
        }
    }

    #[test]
    fn default_weights_snap_plus_genemark_exceed_augustus() {
        let w = EvidenceWeights::default();
        assert!(w.snap + w.genemark > w.augustus);
        assert!(w.augustus > w.genemark && w.genemark >= w.snap);
        assert!(
            w.augustus / w.snap <= 2.0 + 1e-9,
            "ab-initio ratio capped at 2:1"
        );
    }

    // Augustus fuses two adjacent genes; SNAP and GeneMark call them apart.
    // Legacy: the fused model collects partial agreement from all four pieces
    // and is the only model emitted for the group. Corroborated: both real
    // genes are emitted.
    fn fused_case(w: &EvidenceWeights) -> Vec<GeneModel> {
        let fused = gene_model_src(
            100,
            2000,
            w.augustus,
            &[(100, 900), (1200, 2000)],
            "Augustus",
        );
        vec![
            fused,
            gene_model_src(100, 900, w.snap, &[(100, 900)], "SNAP"),
            gene_model_src(1200, 2000, w.snap, &[(1200, 2000)], "SNAP"),
            gene_model_src(100, 900, w.genemark, &[(100, 900)], "GeneMark"),
            gene_model_src(1200, 2000, w.genemark, &[(1200, 2000)], "GeneMark"),
        ]
    }

    #[test]
    fn corroborated_splits_fused_augustus_model() {
        let kept = resolve_overlaps(fused_case(&EvidenceWeights::default()));
        let mut spans: Vec<(u64, u64)> = kept.iter().map(|m| (m.start, m.end)).collect();
        spans.sort();
        assert_eq!(spans, vec![(100, 900), (1200, 2000)]);
    }

    #[test]
    fn legacy_keeps_fused_augustus_model() {
        let w = EvidenceWeights::legacy();
        let (kept, _) = resolve_overlaps_filtered(
            fused_case(&w),
            &ConsensusFilter::default(),
            None,
            1,
            Combiner::Legacy,
        );
        assert_eq!(kept.len(), 1);
        assert_eq!((kept[0].start, kept[0].end), (100, 2000));
    }

    // Two genes linked into one overlap group only through a long model that
    // loses: the neighbour must not be silently dropped.
    #[test]
    fn corroborated_keeps_non_overlapping_neighbour_in_group() {
        let a1 = gene_model_src(100, 500, 4.0, &[(100, 500)], "Augustus");
        let s1 = gene_model_src(100, 500, 2.0, &[(100, 500)], "SNAP");
        let bridge = gene_model_src(450, 1100, 2.0, &[(450, 1100)], "GeneMark");
        let a2 = gene_model_src(1000, 1500, 4.0, &[(1000, 1500)], "Augustus");
        let kept = resolve_overlaps(vec![a1, s1, bridge, a2]);
        let mut spans: Vec<(u64, u64)> = kept.iter().map(|m| (m.start, m.end)).collect();
        spans.sort();
        assert_eq!(spans, vec![(100, 500), (1000, 1500)]);
    }

    // Same structure, different stop conventions (Augustus: CDS + separate
    // stop_codon row; SNAP: stop inside CDS) must score as exact agreement.
    #[test]
    fn stop_codon_convention_does_not_break_exact_agreement() {
        let mut aug = gene_model_src(100, 899, 4.0, &[(100, 300), (400, 896)], "Augustus");
        let mut stop = aug.records[0].clone();
        stop.feature_type = "stop_codon".into();
        stop.start = 897;
        stop.end = 899;
        aug.records.push(stop);
        let snap = gene_model_src(100, 899, 2.0, &[(100, 300), (400, 899)], "SNAP");
        let a = cds_intervals_stop_folded(&aug);
        let b = cds_intervals_stop_folded(&snap);
        assert_eq!(a, b);
        assert!((cds_jaccard(&a, &b) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn combiner_parse_round_trips() {
        for c in [Combiner::Corroborated, Combiner::Legacy] {
            assert_eq!(Combiner::parse(c.name()), Some(c));
        }
        assert_eq!(Combiner::parse("nope"), None);
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

    // Regression: a parent gene/mRNA whose own bounds cover only a SUBSET of
    // its CDS (the SCE_000005 shape — parent 15601-15609 but a second CDS at
    // 15025-15318) must be clamped on emit to span ALL children. Otherwise the
    // other CDS is "outside mRNA bounds" → invalid GFF, table2asn rejects it.
    #[test]
    fn parent_bounds_span_all_cds_children() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.gff3");
        std::fs::write(
            &input,
            "##gff-version 3\n\
scaffold_001\ttest\tgene\t15601\t15609\t.\t-\t.\tID=g1\n\
scaffold_001\ttest\tmRNA\t15601\t15609\t.\t-\t.\tID=g1.t1;Parent=g1\n\
scaffold_001\ttest\tCDS\t15601\t15609\t.\t-\t0\tID=g1.cds1;Parent=g1.t1\n\
scaffold_001\ttest\tCDS\t15025\t15318\t.\t-\t2\tID=g1.cds2;Parent=g1.t1\n",
        )
        .unwrap();

        let out = dir.path().join("consensus.gff3");
        let inputs: Vec<(&Path, &str, f64)> = vec![(input.as_path(), "test", 5.0)];
        let n = merge_predictions(&inputs, &out, "SCE").unwrap();
        assert_eq!(n, 1);

        let gff = std::fs::read_to_string(&out).unwrap();
        // Every emitted gene and mRNA must span the full 15025..15609 CDS union.
        for want_type in ["gene", "mRNA"] {
            let line = gff
                .lines()
                .find(|l| l.split('\t').nth(2) == Some(want_type))
                .unwrap_or_else(|| panic!("no {want_type} line\n{gff}"));
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f[3], "15025", "{want_type} start must be min CDS\n{gff}");
            assert_eq!(f[4], "15609", "{want_type} end must be max CDS\n{gff}");
        }
        // And no CDS may fall outside the parent span.
        for l in gff.lines().filter(|l| l.split('\t').nth(2) == Some("CDS")) {
            let f: Vec<&str> = l.split('\t').collect();
            let (s, e): (u64, u64) = (f[3].parse().unwrap(), f[4].parse().unwrap());
            assert!(
                s >= 15025 && e <= 15609,
                "CDS {s}-{e} outside parent\n{gff}"
            );
        }
    }

    // ── Feature 1: consensus false-positive filter ────────────────────────────

    // A locus called by one low-weight predictor alone is kept under the default
    // (support ≥ 1) but dropped once support ≥ 2 is required.
    #[test]
    fn singleton_dropped_at_support_2_kept_at_support_1() {
        let lone = gene_model_src(100, 200, 3.0, &[(100, 200)], "snap");

        let (kept1, dropped1) = resolve_overlaps_filtered(
            vec![lone.clone()],
            &ConsensusFilter::default(),
            None,
            1,
            Combiner::default(),
        );
        assert_eq!(kept1.len(), 1, "default keeps the lone locus");
        assert_eq!(dropped1, 0);

        let f2 = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: None,
            ..ConsensusFilter::default()
        };
        let (kept2, dropped2) =
            resolve_overlaps_filtered(vec![lone], &f2, None, 1, Combiner::default());
        assert_eq!(
            kept2.len(),
            0,
            "lone low-weight predictor dropped at support ≥ 2"
        );
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
            None,
            1,
            Combiner::default(),
        );
        assert_eq!(k1.len(), 1);
        assert_eq!(d1, 0);

        let f2 = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: None,
            ..ConsensusFilter::default()
        };
        let (k2, d2) = resolve_overlaps_filtered(vec![a, b], &f2, None, 1, Combiner::default());
        assert_eq!(
            k2.len(),
            1,
            "two independent predictors corroborate → kept at support ≥ 2"
        );
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
            ..ConsensusFilter::default()
        };
        let (k, d) = resolve_overlaps_filtered(vec![a, b], &f2, None, 1, Combiner::default());
        assert_eq!(
            k.len(),
            0,
            "same-source duplicates are not independent support"
        );
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
            ..ConsensusFilter::default()
        };
        let (k, d) =
            resolve_overlaps_filtered(vec![lone.clone()], &rescue, None, 1, Combiner::default());
        assert_eq!(
            k.len(),
            1,
            "score 20 ≥ 15 rescues the lone high-confidence locus"
        );
        assert_eq!(d, 0);

        let strict = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: Some(25.0),
            ..ConsensusFilter::default()
        };
        let (k2, d2) = resolve_overlaps_filtered(vec![lone], &strict, None, 1, Combiner::default());
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
        let n_def = merge_predictions_filtered(
            &inputs,
            &out_def,
            "D",
            &ConsensusFilter::default(),
            &ConsensusContext::default(),
        )
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
            ..ConsensusFilter::default()
        };
        let n_f =
            merge_predictions_filtered(&inputs, &out_f, "F", &filter, &ConsensusContext::default())
                .unwrap();
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

    // ── Feature 1b: single-exon-aware filter (BRAKER3 rule) ───────────────────

    // A single-exon, purely-ab-initio locus is dropped by
    // `min_support_single_exon` even when the general `min_predictor_support`
    // would keep it (default 1). This is the targeted BRAKER3 behaviour: cull
    // lone single-exon ab-initio calls, keep everything else.
    #[test]
    fn single_exon_ab_initio_dropped_while_multiexon_kept() {
        // Single-exon lone Augustus locus (one CDS).
        let single = gene_model_src(100, 200, 10.0, &[(100, 200)], "augustus");
        // Multi-exon lone Augustus locus (two CDS) — same support of 1.
        let multi = gene_model_src(1000, 1400, 10.0, &[(1000, 1100), (1300, 1400)], "augustus");

        let f = ConsensusFilter {
            min_predictor_support: 1,   // general: keep everything
            min_support_single_exon: 2, // single-exon: require corroboration
            ..ConsensusFilter::default()
        };
        let (kept, dropped) =
            resolve_overlaps_filtered(vec![single, multi], &f, None, 1, Combiner::default());
        assert_eq!(dropped, 1, "the single-exon lone locus is dropped");
        assert_eq!(kept.len(), 1, "the multi-exon lone locus is kept");
        assert_eq!(
            (kept[0].start, kept[0].end),
            (1000, 1400),
            "survivor is the multi-exon locus"
        );
    }

    // A single-exon locus with EXTRINSIC (protein) support is exempt from the
    // single-exon rule — it falls back to the general threshold.
    #[test]
    fn single_exon_with_protein_support_is_exempt() {
        let lone_protein = gene_model_src(100, 200, 20.0, &[(100, 200)], "protein");
        let f = ConsensusFilter {
            min_predictor_support: 1,
            min_support_single_exon: 2,
            ..ConsensusFilter::default()
        };
        let (kept, dropped) =
            resolve_overlaps_filtered(vec![lone_protein], &f, None, 1, Combiner::default());
        assert_eq!(dropped, 0, "evidence-backed single-exon locus is exempt");
        assert_eq!(kept.len(), 1);
    }

    // A non-coding model (no CDS, e.g. a tRNA) is never treated as single-exon
    // for this rule, so a lone tRNA survives even with a high single-exon bar.
    #[test]
    fn noncoding_locus_exempt_from_single_exon_rule() {
        // No CDS records → cds_record_count == 0.
        let trna = gene_model_src(100, 200, 15.0, &[], "trnascan");
        let f = ConsensusFilter {
            min_predictor_support: 1,
            min_support_single_exon: 5,
            ..ConsensusFilter::default()
        };
        let (kept, dropped) =
            resolve_overlaps_filtered(vec![trna], &f, None, 1, Combiner::default());
        assert_eq!(
            dropped, 0,
            "non-coding locus is exempt from the single-exon drop"
        );
        assert_eq!(kept.len(), 1);
    }

    // The single-exon bar is never weaker than the general bar: a multi-exon
    // locus still obeys `min_predictor_support` even when the single-exon bar
    // is lower.
    #[test]
    fn single_exon_bar_never_below_general_bar() {
        let single = gene_model_src(100, 200, 10.0, &[(100, 200)], "augustus");
        let f = ConsensusFilter {
            min_predictor_support: 2,   // general bar is the stricter one
            min_support_single_exon: 1, // explicitly "off" for single-exon
            ..ConsensusFilter::default()
        };
        // effective bar for the single-exon locus = max(2, 1) = 2 > support 1.
        let (kept, dropped) =
            resolve_overlaps_filtered(vec![single], &f, None, 1, Combiner::default());
        assert_eq!(dropped, 1, "general bar still applies to single-exon loci");
        assert_eq!(kept.len(), 0);
    }

    // ── Feature 2: ORF-completeness gate ──────────────────────────────────────

    #[test]
    fn orf_length_complete_detects_frame() {
        // One CDS of length 101 (200-100+1) → 101 % 3 == 2 → incomplete.
        let bad = gene_model_src(100, 200, 10.0, &[(100, 200)], "augustus");
        assert!(!orf_length_complete(&bad));
        // One CDS of length 102 (201-100+1) → 102 % 3 == 0 → complete.
        let good = gene_model_src(100, 201, 10.0, &[(100, 201)], "augustus");
        assert!(orf_length_complete(&good));
        // Multi-CDS summing to a multiple of 3: 51 + 51 = 102.
        let good_multi = gene_model_src(100, 400, 10.0, &[(100, 150), (300, 350)], "augustus");
        assert!(orf_length_complete(&good_multi));
        // No CDS → non-coding → exempt (true).
        let noncoding = gene_model_src(100, 200, 10.0, &[], "trnascan");
        assert!(orf_length_complete(&noncoding));
    }

    #[test]
    fn require_complete_orf_drops_frame_broken_model() {
        let bad = gene_model_src(100, 200, 10.0, &[(100, 200)], "augustus"); // len 101
        let good = gene_model_src(1000, 1101, 10.0, &[(1000, 1101)], "augustus"); // len 102

        // Off (default): both kept.
        let (kept_off, dropped_off) = resolve_overlaps_filtered(
            vec![bad.clone(), good.clone()],
            &ConsensusFilter::default(),
            None,
            1,
            Combiner::default(),
        );
        assert_eq!(dropped_off, 0);
        assert_eq!(kept_off.len(), 2);

        // On: the frame-broken model is dropped, the clean one kept.
        let f = ConsensusFilter {
            require_complete_orf: true,
            ..ConsensusFilter::default()
        };
        let (kept_on, dropped_on) =
            resolve_overlaps_filtered(vec![bad, good], &f, None, 1, Combiner::default());
        assert_eq!(dropped_on, 1, "frame-broken ORF dropped");
        assert_eq!(kept_on.len(), 1);
        assert_eq!((kept_on[0].start, kept_on[0].end), (1000, 1101));
    }

    // The ORF gate is independent of the score rescue: an incomplete ORF is
    // dropped even when its score would otherwise rescue it.
    #[test]
    fn require_complete_orf_overrides_score_rescue() {
        let bad = gene_model_src(100, 200, 100.0, &[(100, 200)], "augustus"); // len 101
        let f = ConsensusFilter {
            min_predictor_support: 2,
            min_consensus_score: Some(1.0), // would rescue on score alone
            require_complete_orf: true,
            ..ConsensusFilter::default()
        };
        let (kept, dropped) =
            resolve_overlaps_filtered(vec![bad], &f, None, 1, Combiner::default());
        assert_eq!(dropped, 1, "incomplete ORF dropped despite score rescue");
        assert_eq!(kept.len(), 0);
    }

    // Default filter remains a strict no-op even now that it carries the two new
    // fields: nothing is dropped and the fast path is taken.
    #[test]
    fn default_filter_is_still_noop() {
        let f = ConsensusFilter::default();
        assert!(f.is_noop());
        let single = gene_model_src(100, 200, 3.0, &[(100, 200)], "snap"); // len 101, lone
        let (kept, dropped) =
            resolve_overlaps_filtered(vec![single], &f, None, 1, Combiner::default());
        assert_eq!(dropped, 0, "default drops nothing");
        assert_eq!(kept.len(), 1);
    }

    // ── Feature #4: verified-ORF nucleotide check ─────────────────────────────

    fn genome1(seq: &str) -> HashMap<String, FastaRecord> {
        let mut m = HashMap::new();
        m.insert(
            "chr1".to_string(),
            FastaRecord {
                header: "chr1".into(),
                id: "chr1".into(),
                sequence: seq.to_string(),
            },
        );
        m
    }

    fn model_cds(cds: &[(u64, u64)], strand: char) -> GeneModel {
        let records: Vec<GFFRecord> = cds
            .iter()
            .map(|&(s, e)| GFFRecord {
                seqid: "chr1".into(),
                source: "augustus".into(),
                feature_type: "CDS".into(),
                start: s,
                end: e,
                score: None,
                strand,
                phase: Some(0),
                attributes: HashMap::new(),
            })
            .collect();
        let start = cds.iter().map(|c| c.0).min().unwrap();
        let end = cds.iter().map(|c| c.1).max().unwrap();
        GeneModel {
            seqid: "chr1".into(),
            start,
            end,
            strand,
            score: 10.0,
            source: "augustus".into(),
            records,
        }
    }

    #[test]
    fn orf_status_complete_plus_strand() {
        // ATG AAA TAA — start, body, stop; frame OK, no internal stop.
        let g = genome1("ATGAAATAA");
        let m = model_cds(&[(1, 9)], '+');
        assert_eq!(orf_status(&m, Some(&g), 1), OrfStatus::Complete);
    }

    #[test]
    fn orf_status_flags_missing_start_stop_and_internal_stop() {
        assert_eq!(
            orf_status(&model_cds(&[(1, 9)], '+'), Some(&genome1("CCCAAATAA")), 1),
            OrfStatus::NoStart
        );
        assert_eq!(
            orf_status(&model_cds(&[(1, 9)], '+'), Some(&genome1("ATGAAACCC")), 1),
            OrfStatus::NoStop
        );
        // ATG TAA TAA — internal TAA at codon 1, terminal TAA is a valid stop.
        assert_eq!(
            orf_status(&model_cds(&[(1, 9)], '+'), Some(&genome1("ATGTAATAA")), 1),
            OrfStatus::InternalStop
        );
    }

    #[test]
    fn orf_status_bad_length_before_genome_lookup() {
        // 8 bp CDS — frame broken; flagged even with a genome present.
        let g = genome1("ATGAAATA");
        assert_eq!(
            orf_status(&model_cds(&[(1, 8)], '+'), Some(&g), 1),
            OrfStatus::BadLength
        );
    }

    #[test]
    fn orf_status_length_only_without_genome() {
        // No genome → degrade to the frame-only gate (3-divisible = Complete).
        assert_eq!(
            orf_status(&model_cds(&[(1, 9)], '+'), None, 1),
            OrfStatus::Complete
        );
        assert_eq!(
            orf_status(&model_cds(&[(1, 8)], '+'), None, 1),
            OrfStatus::BadLength
        );
    }

    #[test]
    fn orf_status_noncoding_is_exempt() {
        // A model with no CDS (e.g. tRNA) is never an ORF.
        let trna = GeneModel {
            seqid: "chr1".into(),
            start: 1,
            end: 70,
            strand: '+',
            score: 1.0,
            source: "trnascan".into(),
            records: vec![GFFRecord {
                seqid: "chr1".into(),
                source: "trnascan".into(),
                feature_type: "tRNA".into(),
                start: 1,
                end: 70,
                score: None,
                strand: '+',
                phase: None,
                attributes: HashMap::new(),
            }],
        };
        assert_eq!(
            orf_status(&trna, Some(&genome1("ACGT")), 1),
            OrfStatus::Complete
        );
    }

    #[test]
    fn orf_status_minus_strand_complete() {
        // Place the reverse-complement of a complete ORF on the plus strand, so a
        // minus-strand model recovers ATG…TAA after reverse-complementing.
        let plus = crate::parser::fasta::reverse_complement("ATGAAATAA");
        let g = genome1(&plus);
        let m = model_cds(&[(1, plus.len() as u64)], '-');
        assert_eq!(orf_status(&m, Some(&g), 1), OrfStatus::Complete);
    }

    #[test]
    fn orf_status_table6_only_tga_is_stop() {
        // Table 6 reassigns TAA/TAG to Gln; a terminal TAA is then NOT a stop.
        let g = genome1("ATGAAATAA");
        assert_eq!(
            orf_status(&model_cds(&[(1, 9)], '+'), Some(&g), 6),
            OrfStatus::NoStop
        );
        // …but a terminal TGA still is.
        let g2 = genome1("ATGAAATGA");
        assert_eq!(
            orf_status(&model_cds(&[(1, 9)], '+'), Some(&g2), 6),
            OrfStatus::Complete
        );
    }

    // ── Feature #2: fused-gene splitting (jaccard-clip) ───────────────────────

    #[test]
    fn jaccard_clip_splits_on_overlong_gap() {
        // Two CDS clusters 19 600 bp apart — far beyond any real intron.
        let m = model_cds(&[(100, 400), (20_000, 20_300)], '+');
        let parts = split_fused_model(&m, 10_000, 150);
        assert_eq!(parts.len(), 2, "fused model split into two genes");
        assert_eq!(parts[0].start, 100);
        assert_eq!(parts[1].start, 20_000);
    }

    #[test]
    fn jaccard_clip_keeps_normal_introns_intact() {
        // 59 bp gap — an ordinary intron; the model stays one gene.
        let m = model_cds(&[(100, 400), (460, 760)], '+');
        let parts = split_fused_model(&m, 10_000, 150);
        assert_eq!(parts.len(), 1);
    }

    #[test]
    fn jaccard_clip_is_conservative_about_tiny_fragments() {
        // Over-long gap, but the second piece is only 101 bp (< 150): don't split.
        let m = model_cds(&[(100, 400), (20_000, 20_100)], '+');
        let parts = split_fused_model(&m, 10_000, 150);
        assert_eq!(parts.len(), 1, "never fragment below the min coding length");
    }

    #[test]
    fn jaccard_clip_single_exon_never_splits() {
        let m = model_cds(&[(100, 5_000)], '+');
        assert_eq!(split_fused_model(&m, 10_000, 150).len(), 1);
    }

    // ── Feature #3: adaptive repeat filtering ─────────────────────────────────

    #[test]
    fn repeat_mask_reads_lowercase_runs_and_fraction() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("masked.fa");
        // 12 bases, positions 5-8 lowercase (soft-masked repeat).
        std::fs::write(&fa, ">chr1\nACGTacgtACGT\n").unwrap();
        let mask = RepeatMask::from_fasta(&fa).unwrap();
        assert!((mask.masked_fraction - 4.0 / 12.0).abs() < 1e-9);

        // CDS fully inside the repeat run → 100%.
        assert!((mask.cds_repeat_fraction(&model_cds(&[(5, 8)], '+')) - 1.0).abs() < 1e-9);
        // CDS outside any repeat → 0%.
        assert!(mask.cds_repeat_fraction(&model_cds(&[(1, 4)], '+')).abs() < 1e-9);
        // CDS half in repeat → 50%.
        assert!((mask.cds_repeat_fraction(&model_cds(&[(1, 8)], '+')) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn repeat_mask_adaptive_threshold_scales_with_repeat_content() {
        // Compact genome (≤10% masked) → lenient 0.90.
        let dir = tempfile::tempdir().unwrap();
        let fa_lo = dir.path().join("lo.fa");
        std::fs::write(&fa_lo, ">c\nACGTACGTAC\n").unwrap(); // 0% masked
        assert!((RepeatMask::from_fasta(&fa_lo).unwrap().adaptive_threshold() - 0.90).abs() < 1e-9);

        // Repeat-rich genome (≥40% masked) → aggressive 0.50.
        let fa_hi = dir.path().join("hi.fa");
        std::fs::write(&fa_hi, ">c\nacgtacgtac\n").unwrap(); // 100% masked, clamps to 0.40
        assert!((RepeatMask::from_fasta(&fa_hi).unwrap().adaptive_threshold() - 0.50).abs() < 1e-9);
    }

    #[test]
    fn repeat_filter_drops_repeat_heavy_models_only() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("m.fa");
        // chr1: 10 bp, all lowercase (full repeat).
        std::fs::write(&fa, ">chr1\nacgtacgtac\n").unwrap();
        let mask = RepeatMask::from_fasta(&fa).unwrap();
        let in_repeat = model_cds(&[(1, 10)], '+');
        let (kept, dropped) = apply_repeat_filter(vec![in_repeat], &mask, 0.5);
        assert_eq!(dropped, 1);
        assert_eq!(kept.len(), 0);
    }
}
