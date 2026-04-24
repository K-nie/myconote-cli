//! Phased-VCF parser for allele-specific expression (`ase`).
//!
//! Scope (per Q3 in `scratch/ase_spec.md`): zero new crates, ~200 LoC
//! parser that reads gzipped or plain VCF 4.2+, extracts phased SNVs,
//! MNPs (as sequential SNVs), and small indels (ALT ≤ 50 bp) from the
//! first sample column. Rejects unphased heterozygous sites loudly
//! with line numbers.
//!
//! What we deliberately do NOT do:
//! - BCF parsing (use noodles-vcf or bcftools view first).
//! - Multi-sample extraction (we only need one genotype for ASE).
//! - INFO / FILTER filtering (caller can pre-filter with bcftools).
//! - Structural variants (ALT starting with `<` — skipped with a
//!   reason).
//! - Indels longer than the max-indel-size cap (default 50 bp —
//!   skipped; override with `--max-indel-size`).
//! - Statistical phasing — we consume phased output of SHAPEIT /
//!   Beagle / WhatsHap; we do not phase.
//!
//! The parser emits a `Vec<PhasedVariant>` keyed on (chrom, pos). The
//! caller (src/ase/personalize.rs) maps each variant to CDS
//! coordinates.

use crate::utils::error::{MycoNoteError, Result};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// One variant with its phased diploid genotype extracted. We store
/// the REF and per-haplotype ALT allele explicitly — this is enough
/// to apply the variant to a personalized transcriptome later.
#[derive(Debug, Clone, PartialEq)]
pub struct PhasedVariant {
    pub chrom: String,
    /// 1-based VCF POS (the reference coordinate of the first base of REF).
    pub pos: u64,
    pub ref_allele: String,
    /// Per-haplotype ALT. Length is always 2 for diploid input. Each
    /// entry is either a copy of `ref_allele` (when that haplotype
    /// carries REF) or an alternate allele string.
    pub hap_alleles: [String; 2],
    /// True when the two haplotypes differ — i.e. this site is
    /// informative for ASE. Homozygous sites are retained (so reads
    /// from them don't decoy-map) but marked uninformative.
    pub informative: bool,
    /// Best-guess variant category, for the skipped-variants audit
    /// trail and per-haplotype summary counts.
    pub kind: VariantKind,
    /// Source VCF line number (1-based; for error messages).
    pub source_line: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantKind {
    /// REF and ALT are both length 1.
    Snv,
    /// REF and ALT are equal length > 1 (decomposes to multiple SNVs).
    Mnp,
    /// Small indel (max(|REF|, |ALT|) ≤ max_indel_size, not both length 1).
    Indel,
}

/// Why a VCF row was not emitted as a `PhasedVariant`.
#[derive(Debug, Clone, PartialEq)]
pub struct SkippedVariant {
    pub chrom: String,
    pub pos: u64,
    pub reason: SkipReason,
    pub source_line: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SkipReason {
    /// ALT starts with `<` (e.g. `<DEL>`, `<INV>`) or a bracket
    /// (breakend). SVs are out of scope; use an SV-aware tool.
    StructuralVariant,
    /// max(|REF|, |ALT|) > max_indel_size.
    IndelTooLarge { max_allowed: usize, observed: usize },
    /// Multi-allelic site (comma-separated ALT). Callers should
    /// normalize + split these with `bcftools norm -m -any` first.
    MultiAllelic,
    /// Genotype field missing or `./.`.
    MissingGenotype,
    /// Hemizygous (e.g. genotype `0`) — we require diploid for ASE.
    NotDiploid,
}

/// A fatal parse error carries the line number it saw.
#[derive(Debug)]
pub struct VcfParseError {
    pub line: usize,
    pub message: String,
}

/// Parser configuration.
pub struct VcfOptions {
    /// Maximum length, in nucleotides, of an indel to retain. Larger
    /// variants are recorded in the skipped list. Default 50 (same
    /// convention as GATK / bcftools).
    pub max_indel_size: usize,
}

impl Default for VcfOptions {
    fn default() -> Self {
        Self { max_indel_size: 50 }
    }
}

/// Parsed VCF result: the kept variants plus the skipped-variant audit.
#[derive(Debug)]
pub struct VcfParsed {
    pub variants: Vec<PhasedVariant>,
    pub skipped: Vec<SkippedVariant>,
}

/// Parse a phased VCF from disk. Opens `.vcf` or `.vcf.gz`; sniffs by
/// magic bytes rather than file extension so a `.gz` file with no
/// extension is still handled. Phased heterozygous sites are required;
/// the first unphased het fails the parse with a concrete line
/// number. See `parse_phased_vcf_reader` for the reader-based variant
/// used by tests.
pub fn parse_phased_vcf(path: &Path, opts: &VcfOptions) -> Result<VcfParsed> {
    let f = File::open(path).map_err(MycoNoteError::Io)?;
    let mut head = [0u8; 2];
    let mut buf_reader = BufReader::new(f);
    use std::io::BufRead as _;
    let peek = buf_reader.fill_buf()?;
    if peek.len() >= 2 {
        head.copy_from_slice(&peek[..2]);
    }
    // Gzip magic: 1f 8b.
    if head == [0x1f, 0x8b] {
        let gz = GzDecoder::new(buf_reader);
        parse_phased_vcf_reader(BufReader::new(gz), opts).map_err(to_err)
    } else {
        parse_phased_vcf_reader(buf_reader, opts).map_err(to_err)
    }
}

fn to_err(e: VcfParseError) -> MycoNoteError {
    MycoNoteError::QuantSheet(format!("phased VCF line {}: {}", e.line, e.message))
}

/// Parse a phased VCF from any `BufRead` source. Returns a
/// `VcfParseError` on the first unrecoverable issue (malformed
/// header, missing sample column, unphased heterozygous site, etc.).
pub fn parse_phased_vcf_reader<R: BufRead>(
    mut reader: R,
    opts: &VcfOptions,
) -> std::result::Result<VcfParsed, VcfParseError> {
    let mut variants = Vec::new();
    let mut skipped = Vec::new();
    let mut line = String::new();
    let mut line_num = 0usize;
    let mut saw_header = false;
    let mut sample_col_idx: Option<usize> = None;

    loop {
        line.clear();
        let n = reader.read_line(&mut line).map_err(|e| VcfParseError {
            line: line_num,
            message: format!("io error: {e}"),
        })?;
        if n == 0 {
            break;
        }
        line_num += 1;
        let trimmed = line.trim_end_matches('\n').trim_end_matches('\r');
        if trimmed.is_empty() {
            continue;
        }

        // `##` = meta-information, ignore.
        if trimmed.starts_with("##") {
            continue;
        }

        // `#CHROM` = column header. Must appear exactly once and
        // before any data row. First sample column is at index 9.
        if trimmed.starts_with("#CHROM") {
            let cols: Vec<&str> = trimmed.split('\t').collect();
            if cols.len() < 10 {
                return Err(VcfParseError {
                    line: line_num,
                    message: format!(
                        "column header has {} fields; need at least 10 (8 mandatory + FORMAT + at least 1 sample)",
                        cols.len()
                    ),
                });
            }
            sample_col_idx = Some(9);
            saw_header = true;
            continue;
        }

        if !saw_header {
            return Err(VcfParseError {
                line: line_num,
                message: "data row encountered before #CHROM header".into(),
            });
        }

        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 10 {
            return Err(VcfParseError {
                line: line_num,
                message: format!("data row has {} columns; expected ≥10", cols.len()),
            });
        }

        let chrom = cols[0].to_string();
        let pos: u64 = cols[1].parse().map_err(|_| VcfParseError {
            line: line_num,
            message: format!("POS '{}' is not an integer", cols[1]),
        })?;
        let ref_allele = cols[3].to_ascii_uppercase();
        let alt_field = cols[4];
        let format = cols[8];
        let sample = cols[sample_col_idx.unwrap()];

        // Multi-allelic rows — rejected with a skip record. Users
        // should pre-normalize with `bcftools norm -m -any`.
        if alt_field.contains(',') {
            skipped.push(SkippedVariant {
                chrom,
                pos,
                reason: SkipReason::MultiAllelic,
                source_line: line_num,
            });
            continue;
        }

        // Structural variants / breakends — ALT begins with `<` or `[` or `]`.
        let first_alt_char = alt_field.chars().next().unwrap_or('N');
        if matches!(first_alt_char, '<' | '[' | ']') {
            skipped.push(SkippedVariant {
                chrom,
                pos,
                reason: SkipReason::StructuralVariant,
                source_line: line_num,
            });
            continue;
        }

        let alt_allele = alt_field.to_ascii_uppercase();

        // Size cap — max(|REF|, |ALT|) ≤ max_indel_size.
        let observed_size = ref_allele.len().max(alt_allele.len());
        let is_snv = ref_allele.len() == 1 && alt_allele.len() == 1;
        let is_mnp = ref_allele.len() == alt_allele.len() && ref_allele.len() > 1;
        if !is_snv && !is_mnp && observed_size > opts.max_indel_size {
            skipped.push(SkippedVariant {
                chrom,
                pos,
                reason: SkipReason::IndelTooLarge {
                    max_allowed: opts.max_indel_size,
                    observed: observed_size,
                },
                source_line: line_num,
            });
            continue;
        }

        let kind = if is_snv {
            VariantKind::Snv
        } else if is_mnp {
            VariantKind::Mnp
        } else {
            VariantKind::Indel
        };

        // Locate GT in FORMAT. Required field per VCF 4.2.
        let gt_idx = format
            .split(':')
            .position(|k| k == "GT")
            .ok_or_else(|| VcfParseError {
                line: line_num,
                message: format!("FORMAT '{format}' has no GT field"),
            })?;
        let gt_value = sample.split(':').nth(gt_idx).ok_or_else(|| VcfParseError {
            line: line_num,
            message: "sample column missing GT value".into(),
        })?;

        // Missing genotype → skip with reason.
        if gt_value == "." || gt_value == "./." {
            skipped.push(SkippedVariant {
                chrom,
                pos,
                reason: SkipReason::MissingGenotype,
                source_line: line_num,
            });
            continue;
        }

        // Unphased heterozygote → fatal. The whole point of this
        // parser is to reject unphased input loudly.
        let is_phased = gt_value.contains('|');
        let is_het = gt_value.contains('/') || is_phased;
        let sep = if is_phased { '|' } else { '/' };

        let parts: Vec<&str> = gt_value.split(sep).collect();
        if parts.len() != 2 {
            skipped.push(SkippedVariant {
                chrom,
                pos,
                reason: SkipReason::NotDiploid,
                source_line: line_num,
            });
            continue;
        }

        let a: u8 = parts[0].parse().map_err(|_| VcfParseError {
            line: line_num,
            message: format!("GT allele index '{}' is not a number", parts[0]),
        })?;
        let b: u8 = parts[1].parse().map_err(|_| VcfParseError {
            line: line_num,
            message: format!("GT allele index '{}' is not a number", parts[1]),
        })?;

        let is_het_site = a != b;
        if is_het_site && !is_phased {
            return Err(VcfParseError {
                line: line_num,
                message: format!(
                    "unphased heterozygous site {}:{} (GT='{}'); `ase` requires a fully \
                     phased VCF. Run SHAPEIT / Beagle / WhatsHap first, or fix this site \
                     manually.",
                    chrom, pos, gt_value
                ),
            });
        }
        let _ = is_het; // silence unused for homs; kept for clarity.

        let hap = |idx: u8| -> std::result::Result<String, VcfParseError> {
            match idx {
                0 => Ok(ref_allele.clone()),
                1 => Ok(alt_allele.clone()),
                other => Err(VcfParseError {
                    line: line_num,
                    message: format!(
                        "GT allele index {other} out of range (only ALT 1 supported; \
                         split multi-allelics with `bcftools norm -m -any` first)"
                    ),
                }),
            }
        };

        let h0 = hap(a)?;
        let h1 = hap(b)?;
        let informative = h0 != h1;

        variants.push(PhasedVariant {
            chrom,
            pos,
            ref_allele,
            hap_alleles: [h0, h1],
            informative,
            kind,
            source_line: line_num,
        });
    }

    if !saw_header {
        return Err(VcfParseError {
            line: line_num,
            message: "no #CHROM header line found — is this a VCF file?".into(),
        });
    }

    Ok(VcfParsed { variants, skipped })
}

/// Helper used by tests + the main parser to probe gzip by magic bytes
/// without re-reading the file twice on disk.
pub fn looks_like_gzip(buf: &[u8]) -> bool {
    buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
}

/// Count variants per kind — used for the bundle manifest summary.
pub fn summarize_variants(variants: &[PhasedVariant]) -> VariantSummary {
    let mut s = VariantSummary::default();
    for v in variants {
        match v.kind {
            VariantKind::Snv => s.snv += 1,
            VariantKind::Mnp => s.mnp += 1,
            VariantKind::Indel => s.indel += 1,
        }
        if v.informative {
            s.informative += 1;
        } else {
            s.homozygous += 1;
        }
    }
    s.total = variants.len();
    s
}

#[derive(Debug, Default, Clone)]
pub struct VariantSummary {
    pub total: usize,
    pub snv: usize,
    pub mnp: usize,
    pub indel: usize,
    pub informative: usize,
    pub homozygous: usize,
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn parse(body: &str) -> std::result::Result<VcfParsed, VcfParseError> {
        parse_phased_vcf_reader(Cursor::new(body), &VcfOptions::default())
    }

    const MIN_HEADER: &str = "##fileformat=VCFv4.2\n\
                              ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n\
                              ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

    #[test]
    fn parses_single_phased_snv() {
        let vcf = format!("{MIN_HEADER}chrI\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n");
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 1);
        assert_eq!(parsed.skipped.len(), 0);
        let v = &parsed.variants[0];
        assert_eq!(v.chrom, "chrI");
        assert_eq!(v.pos, 100);
        assert_eq!(v.ref_allele, "A");
        assert_eq!(v.hap_alleles, ["A".to_string(), "G".to_string()]);
        assert!(v.informative);
        assert_eq!(v.kind, VariantKind::Snv);
    }

    #[test]
    fn parses_phased_indel_insertion() {
        let vcf = format!("{MIN_HEADER}chrI\t50\t.\tA\tAGGG\t.\tPASS\t.\tGT\t1|0\n");
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 1);
        let v = &parsed.variants[0];
        assert_eq!(v.kind, VariantKind::Indel);
        // Phase: HAP0 gets ALT (first index=1), HAP1 gets REF (second=0).
        assert_eq!(v.hap_alleles[0], "AGGG");
        assert_eq!(v.hap_alleles[1], "A");
        assert!(v.informative);
    }

    #[test]
    fn parses_mnp_as_mnp_kind() {
        let vcf = format!("{MIN_HEADER}chrI\t10\t.\tAC\tTG\t.\tPASS\t.\tGT\t0|1\n");
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 1);
        assert_eq!(parsed.variants[0].kind, VariantKind::Mnp);
    }

    #[test]
    fn rejects_unphased_heterozygote_with_line_number() {
        let vcf = format!(
            "{MIN_HEADER}chrI\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n\
             chrI\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
        );
        let err = parse(&vcf).unwrap_err();
        // Header is 4 lines (3 ##, 1 #CHROM) + 1 OK data + 1 bad = bad on line 6.
        assert_eq!(err.line, 6);
        assert!(err.message.contains("unphased"));
        assert!(err.message.contains("20"));
    }

    #[test]
    fn homozygous_unphased_is_ok() {
        // `0/0` and `1/1` are unambiguous regardless of phasing.
        let vcf = format!(
            "{MIN_HEADER}chrI\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\n\
             chrI\t20\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n"
        );
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 2);
        assert!(!parsed.variants[0].informative);
        assert!(!parsed.variants[1].informative);
        // Hom-alt: both haplotypes carry ALT.
        assert_eq!(
            parsed.variants[1].hap_alleles,
            ["T".to_string(), "T".to_string()]
        );
    }

    #[test]
    fn skips_structural_variants() {
        let vcf = format!("{MIN_HEADER}chrI\t100\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0|1\n");
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 0);
        assert_eq!(parsed.skipped.len(), 1);
        assert!(matches!(
            parsed.skipped[0].reason,
            SkipReason::StructuralVariant
        ));
    }

    #[test]
    fn skips_indels_over_size_cap() {
        let long_alt = "A".to_string() + &"G".repeat(60);
        let vcf = format!("{MIN_HEADER}chrI\t1\t.\tA\t{long_alt}\t.\tPASS\t.\tGT\t0|1\n");
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 0);
        assert_eq!(parsed.skipped.len(), 1);
        assert!(matches!(
            parsed.skipped[0].reason,
            SkipReason::IndelTooLarge {
                max_allowed: 50,
                observed: 61
            }
        ));
    }

    #[test]
    fn skips_multiallelic_rows() {
        let vcf = format!("{MIN_HEADER}chrI\t10\t.\tA\tT,G\t.\tPASS\t.\tGT\t0|1\n");
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 0);
        assert!(matches!(parsed.skipped[0].reason, SkipReason::MultiAllelic));
    }

    #[test]
    fn skips_missing_genotype() {
        let vcf = format!(
            "{MIN_HEADER}chrI\t10\t.\tA\tT\t.\tPASS\t.\tGT\t.\n\
             chrI\t20\t.\tA\tT\t.\tPASS\t.\tGT\t./.\n"
        );
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 0);
        assert_eq!(parsed.skipped.len(), 2);
        assert!(parsed
            .skipped
            .iter()
            .all(|s| matches!(s.reason, SkipReason::MissingGenotype)));
    }

    #[test]
    fn rejects_missing_chrom_header() {
        let vcf = "##fileformat=VCFv4.2\nchrI\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n";
        let err = parse(vcf).unwrap_err();
        assert!(err.message.contains("#CHROM") || err.message.contains("header"));
    }

    #[test]
    fn rejects_truncated_data_row() {
        let vcf = format!("{MIN_HEADER}chrI\t10\t.\tA\n");
        let err = parse(&vcf).unwrap_err();
        assert!(err.message.contains("columns"));
    }

    #[test]
    fn gzip_magic_detected_correctly() {
        // 1f 8b are the gzip magic bytes.
        assert!(looks_like_gzip(&[0x1f, 0x8b, 0x08, 0x00]));
        assert!(!looks_like_gzip(b"##"));
        assert!(!looks_like_gzip(&[]));
    }

    #[test]
    fn summarize_counts_by_kind_and_informative() {
        let vcf = format!(
            "{MIN_HEADER}\
             chrI\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n\
             chrI\t20\t.\tAC\tTG\t.\tPASS\t.\tGT\t0|1\n\
             chrI\t30\t.\tA\tAG\t.\tPASS\t.\tGT\t0|1\n\
             chrI\t40\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\n"
        );
        let parsed = parse(&vcf).unwrap();
        let s = summarize_variants(&parsed.variants);
        assert_eq!(s.total, 4);
        assert_eq!(s.snv, 2);
        assert_eq!(s.mnp, 1);
        assert_eq!(s.indel, 1);
        assert_eq!(s.informative, 3);
        assert_eq!(s.homozygous, 1);
    }

    #[test]
    fn handles_additional_format_fields_around_gt() {
        // GT isn't always first in FORMAT — we search for it by name.
        let vcf = format!("{MIN_HEADER}chrI\t10\t.\tA\tT\t.\tPASS\t.\tDP:GT:GQ\t50:0|1:99\n");
        let parsed = parse(&vcf).unwrap();
        assert_eq!(parsed.variants.len(), 1);
        assert_eq!(
            parsed.variants[0].hap_alleles,
            ["A".to_string(), "T".to_string()]
        );
    }

    #[test]
    fn preserves_source_line_numbers() {
        let vcf = format!(
            "{MIN_HEADER}\
             chrI\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n\
             chrI\t20\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n"
        );
        let parsed = parse(&vcf).unwrap();
        // Header is 4 lines (##fileformat, ##INFO, ##FORMAT, #CHROM).
        // First data variant = line 5; second = line 6.
        assert_eq!(parsed.variants[0].source_line, 5);
        assert_eq!(parsed.variants[1].source_line, 6);
    }
}
