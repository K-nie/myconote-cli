//! Personalized-transcriptome construction for ASE.
//!
//! Given a reference CDS FASTA (as emitted by `convert --to cds`), a
//! GFF3 (to map genome coordinates to CDS coordinates), and a list of
//! phased variants, produce one personalized CDS FASTA per haplotype.
//!
//! The CDS FASTA from `convert --to cds` is already in transcript
//! orientation — reverse-complemented for '-' strand transcripts,
//! spliced across exons. Here we apply phased variants to those
//! spliced sequences, taking strand into account (REF/ALT are in
//! genome-forward-strand orientation in the VCF and need to be
//! reverse-complemented before applying to '-' strand transcripts).
//!
//! Correctness rules (per scratch/ase_spec.md):
//! - Apply variants in descending CDS-position order on each transcript
//!   so an indel at a later position doesn't shift earlier positions.
//! - Variants that overlap in cis on the same transcript's same
//!   haplotype are rejected loudly (the sequential-apply model
//!   can't produce a well-defined output for them).
//! - Variants that span an exon boundary (REF extends past the end of
//!   an exon) are skipped with a recorded reason — we can't safely
//!   apply them without the intron content.
//! - The REF allele is verified against the reference CDS at the
//!   computed position; a mismatch indicates a VCF ↔ assembly
//!   disagreement and is recorded as a skip, not a fatal error
//!   (the user can investigate in variants_skipped.tsv).

use crate::ase::vcf::{PhasedVariant, VariantKind};
use crate::parser::gff::{GFFReader, GFFRecord};
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// Transcript / exon indexing (genome ↔ CDS coordinate mapping)
// ─────────────────────────────────────────────────────────────────────────────

/// One CDS exon of a transcript, as recorded in the GFF3.
#[derive(Debug, Clone, PartialEq)]
pub struct CdsExon {
    pub seqid: String,
    /// 1-based inclusive.
    pub start: u64,
    /// 1-based inclusive.
    pub end: u64,
    /// `+` or `-`. `.` is treated as `+` (unstranded features shouldn't
    /// appear in CDS rows, so this is a defensive default).
    pub strand: char,
    /// GFF3 `phase` column (0, 1, 2). Not used by coordinate math here;
    /// carried through so downstream consumers don't have to re-parse
    /// the GFF3.
    pub phase: u8,
}

/// transcript_id → its CDS exons in genome order (ascending by `start`).
#[derive(Debug, Clone, Default)]
pub struct TranscriptIndex {
    by_tid: HashMap<String, Vec<CdsExon>>,
}

impl TranscriptIndex {
    /// Build the index by walking a GFF3 and collecting CDS rows per
    /// transcript. Transcripts with zero CDS rows are omitted; the
    /// caller already has them in the CDS FASTA if they're protein
    /// coding, so silence isn't a real hazard.
    pub fn from_gff3(path: &Path) -> Result<Self> {
        let mut reader = GFFReader::from_path(path)?;
        let mut by_tid: HashMap<String, Vec<CdsExon>> = HashMap::new();

        while let Some(rec) = reader.next_record()? {
            if rec.feature_type != "CDS" {
                continue;
            }
            let Some(parent) = cds_parent(&rec) else {
                // CDS rows must carry a Parent linking back to the mRNA.
                // Rows without one can't be matched to a transcript ID,
                // so we skip — the summary can mention the count later.
                continue;
            };
            by_tid.entry(parent).or_default().push(CdsExon {
                seqid: rec.seqid.clone(),
                start: rec.start,
                end: rec.end,
                strand: rec.strand,
                phase: rec.phase.unwrap_or(0),
            });
        }

        // Normalize: sort each transcript's exons by ascending genome
        // start. The genome-to-CDS map honors strand when iterating.
        for exons in by_tid.values_mut() {
            exons.sort_by_key(|e| e.start);
        }

        Ok(Self { by_tid })
    }

    pub fn get(&self, tid: &str) -> Option<&Vec<CdsExon>> {
        self.by_tid.get(tid)
    }

    pub fn transcripts(&self) -> impl Iterator<Item = &String> {
        self.by_tid.keys()
    }

    /// Every transcript whose CDS covers `(chrom, pos)` (inclusive on
    /// both endpoints). Returns (transcript_id, exon_index) pairs.
    /// Most variants will hit 0–2 transcripts; overlapping CDS on
    /// opposite strands is the common case for 2.
    pub fn transcripts_at(&self, chrom: &str, pos: u64) -> Vec<(String, usize)> {
        let mut out = Vec::new();
        for (tid, exons) in &self.by_tid {
            for (i, e) in exons.iter().enumerate() {
                if e.seqid == chrom && pos >= e.start && pos <= e.end {
                    out.push((tid.clone(), i));
                }
            }
        }
        out
    }
}

/// Extract the `Parent=` attribute from a CDS row. GFF3 allows
/// multiple parents (comma-separated) but for CDS that is extremely
/// rare and usually indicates a malformed file. We take the first.
fn cds_parent(rec: &GFFRecord) -> Option<String> {
    rec.attributes
        .get("Parent")
        .map(|s| s.split(',').next().unwrap_or("").trim().to_string())
}

/// Map a genome position to a 0-based CDS position for a given
/// transcript, honoring strand. Returns `None` if the genome position
/// doesn't land inside any CDS exon of the transcript.
///
/// For `+` strand: the CDS is read low-to-high along the genome,
/// offset within an exon = pos - exon.start.
/// For `-` strand: the CDS is read high-to-low along the genome
/// after reverse-complementing, so offset within an exon = exon.end
/// - pos.
pub fn genome_to_cds_pos(exons: &[CdsExon], pos: u64) -> Option<u64> {
    if exons.is_empty() {
        return None;
    }
    let strand = exons[0].strand;
    // Validate all exons agree on strand. A mixed-strand transcript
    // is malformed; silently drop it.
    if exons.iter().any(|e| e.strand != strand) {
        return None;
    }
    // Iterate in transcript order.
    let mut iter: Vec<&CdsExon> = exons.iter().collect();
    if strand == '+' {
        iter.sort_by_key(|e| e.start);
    } else {
        iter.sort_by_key(|e| std::cmp::Reverse(e.start));
    }
    let mut cumulative: u64 = 0;
    for exon in iter {
        if pos >= exon.start && pos <= exon.end {
            let offset = if strand == '+' {
                pos - exon.start
            } else {
                exon.end - pos
            };
            return Some(cumulative + offset);
        }
        cumulative += exon.end - exon.start + 1;
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Variant application
// ─────────────────────────────────────────────────────────────────────────────

/// One variant that got applied to one transcript, for the audit
/// trail (variants_applied.tsv).
#[derive(Debug, Clone)]
pub struct VariantApplication {
    pub transcript_id: String,
    pub haplotype: String,
    pub chrom: String,
    pub genome_pos: u64,
    pub cds_pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub kind: VariantKind,
    pub strand: char,
}

/// One variant that could not be applied to a transcript, with reason
/// (variants_skipped.tsv).
#[derive(Debug, Clone)]
pub struct VariantSkip {
    pub transcript_id: String,
    pub haplotype: String,
    pub chrom: String,
    pub genome_pos: u64,
    pub reason: ApplicationSkipReason,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ApplicationSkipReason {
    /// Variant position falls outside any CDS of this transcript.
    OutsideCds,
    /// Variant REF extends past the end of an exon (spans an intron
    /// boundary). We can't safely apply it without the intron content.
    SpansExonBoundary,
    /// Observed reference base(s) at the computed CDS position don't
    /// match the VCF REF — indicates a VCF ↔ assembly disagreement.
    RefMismatch { observed: String, expected: String },
    /// Another variant on the same haplotype already modified an
    /// overlapping CDS range. Rejected per scratch/ase_spec.md Q7
    /// compound-het handling: the sequential-apply model can't
    /// produce well-defined output for in-cis overlapping variants.
    InCisOverlap,
}

/// Per-haplotype personalized CDS: transcript_id → sequence bytes.
#[derive(Debug, Clone)]
pub struct HaplotypeCds {
    pub name: String,
    pub sequences: HashMap<String, Vec<u8>>,
}

/// Full output of the personalization step.
#[derive(Debug)]
pub struct PersonalizationOutput {
    pub haplotypes: Vec<HaplotypeCds>,
    pub applied: Vec<VariantApplication>,
    pub skipped: Vec<VariantSkip>,
}

/// Build personalized CDS sequences for two haplotypes given the
/// reference CDS (already read into memory), a transcript → exon
/// index, the phased variants, and the two haplotype names.
///
/// The input `reference_cds` is consumed and cloned once per
/// haplotype; callers can reload from the on-disk FASTA if they need
/// the reference again.
pub fn personalize(
    reference_cds: &HashMap<String, Vec<u8>>,
    tx_index: &TranscriptIndex,
    variants: &[PhasedVariant],
    haplotype_names: &[String; 2],
) -> Result<PersonalizationOutput> {
    let mut applied: Vec<VariantApplication> = Vec::new();
    let mut skipped: Vec<VariantSkip> = Vec::new();
    let mut haplotypes: Vec<HaplotypeCds> = Vec::with_capacity(2);

    for (hap_idx, hap_name) in haplotype_names.iter().enumerate() {
        let mut sequences: HashMap<String, Vec<u8>> = reference_cds.clone();
        let mut per_tx_applied: HashMap<String, Vec<ApplicationRecord>> = HashMap::new();

        // Phase 1: map each variant to every transcript it hits and
        // compute the CDS position / range. Don't mutate yet — we
        // need all hits before we can sort by CDS position.
        for v in variants {
            let hits = tx_index.transcripts_at(&v.chrom, v.pos);
            if hits.is_empty() {
                // Variant doesn't hit any annotated CDS. Not a
                // per-transcript skip — record once per variant in
                // the VCF-level skipped list rather than fanning
                // out. Here, we simply don't emit anything.
                continue;
            }
            for (tid, _exon_idx) in hits {
                let Some(exons) = tx_index.get(&tid) else {
                    continue;
                };
                let alt_for_hap = &v.hap_alleles[hap_idx];
                // No-op when this haplotype carries the REF allele —
                // no substitution to apply, no audit record to emit.
                // Keeps variants_applied.tsv focused on actual edits.
                if alt_for_hap == &v.ref_allele {
                    continue;
                }
                let strand = exons[0].strand;
                let ref_end = v.pos + v.ref_allele.len() as u64 - 1;

                // Span-an-exon-boundary check: for the variant's REF
                // span [pos, ref_end], every base must map to a CDS
                // position inside the same transcript.
                let cds_start = match genome_to_cds_pos(exons, v.pos) {
                    Some(p) => p,
                    None => {
                        skipped.push(VariantSkip {
                            transcript_id: tid.clone(),
                            haplotype: hap_name.clone(),
                            chrom: v.chrom.clone(),
                            genome_pos: v.pos,
                            reason: ApplicationSkipReason::OutsideCds,
                        });
                        continue;
                    }
                };
                let cds_end_opt = genome_to_cds_pos(exons, ref_end);
                let cds_end = match cds_end_opt {
                    Some(p) => p,
                    None => {
                        skipped.push(VariantSkip {
                            transcript_id: tid.clone(),
                            haplotype: hap_name.clone(),
                            chrom: v.chrom.clone(),
                            genome_pos: v.pos,
                            reason: ApplicationSkipReason::SpansExonBoundary,
                        });
                        continue;
                    }
                };

                // cds_start / cds_end are the 0-based CDS positions
                // of the FIRST and LAST bases of REF. For '+' strand
                // cds_start <= cds_end; for '-' strand the reverse.
                let (lo, hi) = if cds_start <= cds_end {
                    (cds_start, cds_end)
                } else {
                    (cds_end, cds_start)
                };

                per_tx_applied
                    .entry(tid.clone())
                    .or_default()
                    .push(ApplicationRecord {
                        variant: v.clone(),
                        cds_lo: lo as usize,
                        cds_hi: hi as usize,
                        strand,
                        alt_for_hap: alt_for_hap.clone(),
                    });
            }
        }

        // Phase 2: per transcript, sort variants by descending CDS
        // position (high-first). Apply in that order so indels that
        // shift downstream positions don't invalidate not-yet-applied
        // CDS coordinates for other variants on the same transcript.
        for (tid, mut recs) in per_tx_applied {
            recs.sort_by_key(|r| std::cmp::Reverse(r.cds_lo));

            // Check for in-cis overlap between any two variants on
            // the same transcript for this haplotype. Overlap means
            // the [cds_lo, cds_hi] ranges intersect. Compute the
            // overlap decision in a separate pass so the later
            // application loop can mutate `recs` without aliasing
            // the check's immutable borrow.
            let mut overlap_found = false;
            'outer: for (i, a) in recs.iter().enumerate() {
                for b in recs.iter().skip(i + 1) {
                    if ranges_overlap(a.cds_lo, a.cds_hi, b.cds_lo, b.cds_hi) {
                        overlap_found = true;
                        // Record ALL records on this transcript as
                        // skipped — user wants to know every variant
                        // they lost, not just the pair that collided.
                        for r in &recs {
                            skipped.push(VariantSkip {
                                transcript_id: tid.clone(),
                                haplotype: hap_name.clone(),
                                chrom: r.variant.chrom.clone(),
                                genome_pos: r.variant.pos,
                                reason: ApplicationSkipReason::InCisOverlap,
                            });
                        }
                        break 'outer;
                    }
                }
            }
            if overlap_found {
                // Abort all edits on this transcript for this
                // haplotype — leaving it as reference is safer
                // than applying partial edits.
                recs.clear();
            }

            let Some(seq) = sequences.get_mut(&tid) else {
                // No reference CDS entry for this transcript — the
                // GFF3 had a CDS row but convert --to cds didn't
                // emit the record. Skip silently; the bundle will
                // reflect the count via transcripts-in-gff-not-in-cds.
                continue;
            };

            for rec in recs {
                // REF / ALT orientation. For '-' strand, the CDS
                // sequence stored is the reverse-complement of the
                // genome; both REF and ALT must be RC'd before
                // applying.
                let ref_bytes: Vec<u8> = if rec.strand == '+' {
                    rec.variant.ref_allele.bytes().collect()
                } else {
                    reverse_complement(rec.variant.ref_allele.as_bytes())
                };
                let alt_bytes: Vec<u8> = if rec.strand == '+' {
                    rec.alt_for_hap.bytes().collect()
                } else {
                    reverse_complement(rec.alt_for_hap.as_bytes())
                };

                // The [rec.cds_lo, rec.cds_hi] range is inclusive
                // and has length rec.cds_hi - rec.cds_lo + 1. On '+'
                // strand that equals REF.len(); on '-' strand it
                // also equals REF.len() because both strand maps
                // preserve length.
                let span_len = rec.cds_hi - rec.cds_lo + 1;
                if span_len != ref_bytes.len() {
                    // Defensive: should never hit this — means the
                    // coordinate math disagreed with the length of
                    // REF. Record as a skip and move on.
                    skipped.push(VariantSkip {
                        transcript_id: tid.clone(),
                        haplotype: hap_name.clone(),
                        chrom: rec.variant.chrom.clone(),
                        genome_pos: rec.variant.pos,
                        reason: ApplicationSkipReason::RefMismatch {
                            observed: format!("span {span_len} bases"),
                            expected: format!("REF {} bases", ref_bytes.len()),
                        },
                    });
                    continue;
                }

                // Verify the observed reference bases match VCF REF.
                if seq.len() < rec.cds_hi + 1 {
                    // Transcript shorter than the computed position —
                    // also defensive, shouldn't happen if convert --to
                    // cds and the GFF3 agree.
                    skipped.push(VariantSkip {
                        transcript_id: tid.clone(),
                        haplotype: hap_name.clone(),
                        chrom: rec.variant.chrom.clone(),
                        genome_pos: rec.variant.pos,
                        reason: ApplicationSkipReason::OutsideCds,
                    });
                    continue;
                }
                let observed: Vec<u8> = seq[rec.cds_lo..=rec.cds_hi]
                    .iter()
                    .map(|b| b.to_ascii_uppercase())
                    .collect();
                let expected: Vec<u8> = ref_bytes.iter().map(|b| b.to_ascii_uppercase()).collect();
                if observed != expected {
                    skipped.push(VariantSkip {
                        transcript_id: tid.clone(),
                        haplotype: hap_name.clone(),
                        chrom: rec.variant.chrom.clone(),
                        genome_pos: rec.variant.pos,
                        reason: ApplicationSkipReason::RefMismatch {
                            observed: String::from_utf8_lossy(&observed).into_owned(),
                            expected: String::from_utf8_lossy(&expected).into_owned(),
                        },
                    });
                    continue;
                }

                // Apply the edit: replace the REF bytes with ALT
                // bytes. For SNVs and MNPs, the replacement is the
                // same length so this is an in-place overwrite. For
                // indels, the sequence length changes.
                seq.splice(rec.cds_lo..=rec.cds_hi, alt_bytes.iter().copied());

                applied.push(VariantApplication {
                    transcript_id: tid.clone(),
                    haplotype: hap_name.clone(),
                    chrom: rec.variant.chrom.clone(),
                    genome_pos: rec.variant.pos,
                    cds_pos: rec.cds_lo as u64,
                    ref_allele: rec.variant.ref_allele.clone(),
                    alt_allele: rec.alt_for_hap.clone(),
                    kind: rec.variant.kind,
                    strand: rec.strand,
                });
            }
        }

        haplotypes.push(HaplotypeCds {
            name: hap_name.clone(),
            sequences,
        });
    }

    Ok(PersonalizationOutput {
        haplotypes,
        applied,
        skipped,
    })
}

/// Internal per-hit record used during sorting before application.
struct ApplicationRecord {
    variant: PhasedVariant,
    cds_lo: usize,
    cds_hi: usize,
    strand: char,
    alt_for_hap: String,
}

fn ranges_overlap(a_lo: usize, a_hi: usize, b_lo: usize, b_hi: usize) -> bool {
    a_lo <= b_hi && b_lo <= a_hi
}

/// Reverse-complement a DNA sequence. Uppercase in, uppercase out.
/// Non-ACGTN characters pass through unchanged (IUPAC ambiguity codes
/// aren't emitted by most VCF pipelines, so we don't handle them).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            other => other,
        })
        .collect()
}

/// Read a FASTA file into a transcript_id → sequence map. Uses the
/// first whitespace-delimited token of each `>` header as the ID, to
/// match how `convert --to cds` writes records.
pub fn read_cds_fasta(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let reader = BufReader::new(File::open(path).map_err(MycoNoteError::Io)?);
    let mut out: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_id: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                out.insert(id, std::mem::take(&mut current_seq));
            }
            current_id = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else if current_id.is_some() {
            for b in line.bytes() {
                if !b.is_ascii_whitespace() {
                    current_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    if let Some(id) = current_id {
        out.insert(id, current_seq);
    }
    Ok(out)
}

/// Write a transcript_id → sequence map back to a FASTA file with a
/// fixed 60-character line length. Records emitted in sorted order
/// for determinism.
pub fn write_cds_fasta(path: &Path, sequences: &HashMap<String, Vec<u8>>) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};
    let mut w = BufWriter::new(File::create(path).map_err(MycoNoteError::Io)?);
    let mut ids: Vec<&String> = sequences.keys().collect();
    ids.sort();
    for id in ids {
        let seq = &sequences[id];
        writeln!(w, ">{id}").map_err(MycoNoteError::Io)?;
        for chunk in seq.chunks(60) {
            w.write_all(chunk).map_err(MycoNoteError::Io)?;
            writeln!(w).map_err(MycoNoteError::Io)?;
        }
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ase::vcf::{PhasedVariant, VariantKind};

    fn exon(start: u64, end: u64, strand: char) -> CdsExon {
        CdsExon {
            seqid: "chrI".to_string(),
            start,
            end,
            strand,
            phase: 0,
        }
    }

    // ── genome_to_cds_pos: '+' strand ───────────────────────────────────────

    #[test]
    fn plus_strand_single_exon() {
        let exons = vec![exon(100, 200, '+')];
        // Exon [100,200] = 101 bp at CDS positions 0..=100.
        assert_eq!(genome_to_cds_pos(&exons, 100), Some(0));
        assert_eq!(genome_to_cds_pos(&exons, 150), Some(50));
        assert_eq!(genome_to_cds_pos(&exons, 200), Some(100));
        assert_eq!(genome_to_cds_pos(&exons, 99), None);
        assert_eq!(genome_to_cds_pos(&exons, 201), None);
    }

    #[test]
    fn plus_strand_multiple_exons() {
        // Exon1 [100,200] length 101 → CDS 0..=100
        // Exon2 [300,400] length 101 → CDS 101..=201
        // Exon3 [500,600] length 101 → CDS 202..=302
        let exons = vec![
            exon(100, 200, '+'),
            exon(300, 400, '+'),
            exon(500, 600, '+'),
        ];
        assert_eq!(genome_to_cds_pos(&exons, 100), Some(0));
        assert_eq!(genome_to_cds_pos(&exons, 200), Some(100));
        assert_eq!(genome_to_cds_pos(&exons, 300), Some(101));
        assert_eq!(genome_to_cds_pos(&exons, 350), Some(151));
        assert_eq!(genome_to_cds_pos(&exons, 500), Some(202));
        assert_eq!(genome_to_cds_pos(&exons, 600), Some(302));
        // Intron positions map to None.
        assert_eq!(genome_to_cds_pos(&exons, 250), None);
        assert_eq!(genome_to_cds_pos(&exons, 450), None);
    }

    // ── genome_to_cds_pos: '-' strand ───────────────────────────────────────

    #[test]
    fn minus_strand_single_exon() {
        // '-' strand, exon [100,200]. CDS reads 200 → 100 (RC'd).
        // Genome 200 → CDS 0, genome 100 → CDS 100.
        let exons = vec![exon(100, 200, '-')];
        assert_eq!(genome_to_cds_pos(&exons, 200), Some(0));
        assert_eq!(genome_to_cds_pos(&exons, 150), Some(50));
        assert_eq!(genome_to_cds_pos(&exons, 100), Some(100));
    }

    #[test]
    fn minus_strand_multiple_exons() {
        // '-' strand, exons [100,200] and [300,400].
        // Transcript order: exon2 first (higher start), then exon1.
        // CDS 0..=100 = RC(genome[300..=400]).
        // CDS 101..=201 = RC(genome[100..=200]).
        // Genome 400 → CDS 0, genome 300 → CDS 100.
        // Genome 200 → CDS 101, genome 100 → CDS 201.
        let exons = vec![exon(100, 200, '-'), exon(300, 400, '-')];
        assert_eq!(genome_to_cds_pos(&exons, 400), Some(0));
        assert_eq!(genome_to_cds_pos(&exons, 350), Some(50));
        assert_eq!(genome_to_cds_pos(&exons, 300), Some(100));
        assert_eq!(genome_to_cds_pos(&exons, 200), Some(101));
        assert_eq!(genome_to_cds_pos(&exons, 100), Some(201));
    }

    // ── reverse_complement ──────────────────────────────────────────────────

    #[test]
    fn reverse_complement_basic() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT".to_vec());
        assert_eq!(reverse_complement(b"AATTCCGG"), b"CCGGAATT".to_vec());
        assert_eq!(reverse_complement(b"A"), b"T".to_vec());
    }

    #[test]
    fn reverse_complement_preserves_n() {
        assert_eq!(reverse_complement(b"ACNGT"), b"ACNGT".to_vec());
    }

    // ── personalize: SNVs on '+' strand ─────────────────────────────────────

    fn snv(chrom: &str, pos: u64, r: &str, h0: &str, h1: &str) -> PhasedVariant {
        PhasedVariant {
            chrom: chrom.to_string(),
            pos,
            ref_allele: r.to_string(),
            hap_alleles: [h0.to_string(), h1.to_string()],
            informative: h0 != h1,
            kind: VariantKind::Snv,
            source_line: 1,
        }
    }

    fn build_tx_index(entries: &[(&str, Vec<CdsExon>)]) -> TranscriptIndex {
        let mut tx = TranscriptIndex::default();
        for (tid, exons) in entries {
            tx.by_tid.insert(tid.to_string(), exons.clone());
        }
        tx
    }

    #[test]
    fn snv_plus_strand_applies_to_one_transcript() {
        // tx1 on '+' strand: exon [100,200]. Ref CDS = 101 'A's.
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '+')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'A'; 101]);

        // SNV at genome 150: REF=A, HAP0=A, HAP1=G. Hap1 differs.
        let variants = vec![snv("chrI", 150, "A", "A", "G")];
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &variants, &haps).unwrap();
        assert_eq!(out.haplotypes.len(), 2);
        // hap0 = all A (unchanged)
        assert_eq!(out.haplotypes[0].sequences["tx1"], vec![b'A'; 101]);
        // hap1 = A...A G A...A with G at CDS 50
        let mut expected = vec![b'A'; 101];
        expected[50] = b'G';
        assert_eq!(out.haplotypes[1].sequences["tx1"], expected);
        assert_eq!(out.applied.len(), 1, "one application on hap1 tx1");
    }

    #[test]
    fn snv_minus_strand_reverse_complements_alt() {
        // tx1 on '-' strand: exon [100,200]. Reference CDS = 101 'T's
        // (which is RC of 101 'A's in the genome). VCF REF is
        // genome-strand — we pass REF='A', HAP1='C'. Applying to the
        // CDS should insert 'G' (RC of 'C') at the right CDS position.
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '-')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'T'; 101]);

        // Variant at genome pos 150 → CDS pos 50. REF='A' (genome),
        // which in CDS is RC('A')='T' at position 50. After applying
        // hap1 ALT='C', CDS position 50 becomes 'G' (RC of 'C').
        let variants = vec![snv("chrI", 150, "A", "A", "C")];
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &variants, &haps).unwrap();
        let mut expected = vec![b'T'; 101];
        expected[50] = b'G';
        assert_eq!(out.haplotypes[1].sequences["tx1"], expected);
        // hap0 unchanged.
        assert_eq!(out.haplotypes[0].sequences["tx1"], vec![b'T'; 101]);
    }

    #[test]
    fn homozygous_snv_same_on_both_haps_but_not_reference() {
        // Hom-ALT: GT=1/1 → both haplotypes carry G. CDS at both haps
        // should have the SNV; reference CDS is unchanged.
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '+')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'A'; 101]);

        let variants = vec![snv("chrI", 150, "A", "G", "G")];
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &variants, &haps).unwrap();
        let mut expected = vec![b'A'; 101];
        expected[50] = b'G';
        assert_eq!(out.haplotypes[0].sequences["tx1"], expected);
        assert_eq!(out.haplotypes[1].sequences["tx1"], expected);
        assert_eq!(
            out.applied.len(),
            2,
            "both haps record the hom-alt application"
        );
    }

    #[test]
    fn variant_outside_cds_is_skipped() {
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '+')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'A'; 101]);

        // Variant at 500 — outside the CDS. transcripts_at returns
        // empty; no per-transcript skip emitted. No applications either.
        let variants = vec![snv("chrI", 500, "A", "A", "G")];
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &variants, &haps).unwrap();
        assert_eq!(out.applied.len(), 0);
        assert_eq!(out.haplotypes[0].sequences["tx1"], vec![b'A'; 101]);
        assert_eq!(out.haplotypes[1].sequences["tx1"], vec![b'A'; 101]);
    }

    #[test]
    fn ref_mismatch_produces_skip_not_panic() {
        // CDS is all 'A'; variant claims REF='T' at position 150.
        // Tool should record a RefMismatch skip, not panic or silently
        // apply with garbage.
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '+')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'A'; 101]);

        let variants = vec![snv("chrI", 150, "T", "T", "C")];
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &variants, &haps).unwrap();
        assert_eq!(out.applied.len(), 0);
        let ref_mismatches = out
            .skipped
            .iter()
            .filter(|s| matches!(s.reason, ApplicationSkipReason::RefMismatch { .. }))
            .count();
        assert!(
            ref_mismatches >= 1,
            "should record ref mismatch: {:?}",
            out.skipped
        );
    }

    // ── indels + in-cis overlap detection ───────────────────────────────────

    #[test]
    fn in_cis_overlapping_variants_rejected() {
        // Two SNVs at nearby positions on the same transcript + same
        // haplotype. They're at distinct positions (non-overlapping
        // ranges), so they should BOTH apply successfully.
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '+')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'A'; 101]);

        let v1 = snv("chrI", 150, "A", "A", "G"); // hap1 G
        let v2 = snv("chrI", 160, "A", "A", "C"); // hap1 C
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &[v1, v2], &haps).unwrap();
        let h1_seq = &out.haplotypes[1].sequences["tx1"];
        assert_eq!(h1_seq[50], b'G');
        assert_eq!(h1_seq[60], b'C');
    }

    #[test]
    fn overlapping_indel_and_snv_flagged() {
        // Indel at CDS 50 spanning 3 bases; SNV at CDS 51. Their
        // ranges [50,52] and [51,51] overlap → both get flagged as
        // InCisOverlap and neither is applied.
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '+')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'A'; 101]);

        let indel = PhasedVariant {
            chrom: "chrI".to_string(),
            pos: 150,
            ref_allele: "AAA".to_string(),
            hap_alleles: ["AAA".to_string(), "A".to_string()],
            informative: true,
            kind: VariantKind::Indel,
            source_line: 1,
        };
        let s = snv("chrI", 151, "A", "A", "G");
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &[indel, s], &haps).unwrap();
        // No applications on hap1 tx1 (all rejected).
        let h1_applied = out.applied.iter().filter(|a| a.haplotype == "hap1").count();
        assert_eq!(
            h1_applied, 0,
            "no applications when in-cis overlap detected"
        );
        let overlap_skips = out
            .skipped
            .iter()
            .filter(|s| matches!(s.reason, ApplicationSkipReason::InCisOverlap))
            .count();
        assert!(overlap_skips >= 2, "both variants logged as overlap skips");
    }

    #[test]
    fn insertion_applies_correctly_on_plus_strand() {
        // REF=A at genome 150 → CDS pos 50 on a '+' strand transcript
        // of 101 bp all 'A'. hap1 ALT=ATCG. After apply, hap1 CDS
        // length = 101 - 1 + 4 = 104. Position 50 becomes A; positions
        // 51-53 become T,C,G; positions 54-103 are the remaining
        // reference A's.
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '+')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'A'; 101]);

        let ins = PhasedVariant {
            chrom: "chrI".to_string(),
            pos: 150,
            ref_allele: "A".to_string(),
            hap_alleles: ["A".to_string(), "ATCG".to_string()],
            informative: true,
            kind: VariantKind::Indel,
            source_line: 1,
        };
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &[ins], &haps).unwrap();
        let h1 = &out.haplotypes[1].sequences["tx1"];
        assert_eq!(h1.len(), 104, "insertion extends CDS by 3 bases");
        assert_eq!(h1[50], b'A');
        assert_eq!(h1[51], b'T');
        assert_eq!(h1[52], b'C');
        assert_eq!(h1[53], b'G');
        assert_eq!(h1[54], b'A');
    }

    #[test]
    fn deletion_applies_correctly_on_plus_strand() {
        // REF=AAA at genome 150 (spans 150,151,152 → CDS 50,51,52)
        // on '+' strand. hap1 ALT=A (2-base deletion). New length
        // 101 - 3 + 1 = 99.
        let tx_index = build_tx_index(&[("tx1", vec![exon(100, 200, '+')])]);
        let mut cds: HashMap<String, Vec<u8>> = HashMap::new();
        cds.insert("tx1".to_string(), vec![b'A'; 101]);

        let del = PhasedVariant {
            chrom: "chrI".to_string(),
            pos: 150,
            ref_allele: "AAA".to_string(),
            hap_alleles: ["AAA".to_string(), "A".to_string()],
            informative: true,
            kind: VariantKind::Indel,
            source_line: 1,
        };
        let haps = ["hap0".to_string(), "hap1".to_string()];

        let out = personalize(&cds, &tx_index, &[del], &haps).unwrap();
        let h1 = &out.haplotypes[1].sequences["tx1"];
        assert_eq!(h1.len(), 99, "deletion removes 2 bases");
        assert_eq!(h1.iter().all(|&b| b == b'A'), true);
    }

    // ── CDS FASTA round-trip ────────────────────────────────────────────────

    #[test]
    fn fasta_round_trip_preserves_content() {
        let tmp = tempfile::TempDir::new().unwrap();
        let p = tmp.path().join("cds.fa");

        let mut in_map: HashMap<String, Vec<u8>> = HashMap::new();
        in_map.insert("tx1".to_string(), b"ACGTACGTACGT".to_vec());
        in_map.insert("tx2".to_string(), b"AAAAATTTTT".to_vec());

        write_cds_fasta(&p, &in_map).unwrap();
        let back = read_cds_fasta(&p).unwrap();
        assert_eq!(back["tx1"], b"ACGTACGTACGT".to_vec());
        assert_eq!(back["tx2"], b"AAAAATTTTT".to_vec());
    }

    #[test]
    fn fasta_reader_uppercases_and_strips_whitespace() {
        let tmp = tempfile::TempDir::new().unwrap();
        let p = tmp.path().join("cds.fa");
        std::fs::write(
            &p,
            ">tx1 some extra header stuff\nacgt\nACGT\n\n>tx2\nAAAA\n",
        )
        .unwrap();
        let m = read_cds_fasta(&p).unwrap();
        assert_eq!(m["tx1"], b"ACGTACGT".to_vec());
        assert_eq!(m["tx2"], b"AAAA".to_vec());
    }
}
