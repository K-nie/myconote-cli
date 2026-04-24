//! Reproducibility manifest for an ASE run (`ase_bundle.json`).
//!
//! Extends the schema documented in `scratch/ase_spec.md` — same
//! pattern as `src/quant/bundle.rs` but carries haplotype-specific
//! fields: personalized-transcriptome SHA256s, per-haplotype
//! variant-application counts, per-sample per-haplotype mapping
//! rates, and the VCF SHA256 alongside the CDS + genome hashes.

use crate::quant::bundle::{ExternalTools, FastpSummary, HashedPath};
use crate::utils::error::Result;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// Top-level bundle
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Serialize, Deserialize)]
pub struct AseBundle {
    pub version: String,
    pub run_timestamp: String,
    pub command: String,
    pub external_tools: ExternalTools,
    pub inputs: AseBundleInputs,
    pub haplotypes: Vec<String>,
    pub personalized_transcriptomes: Vec<HaplotypeTranscriptome>,
    pub variants: VariantApplicationSummary,
    pub samples: Vec<AseSampleSummary>,
    pub seed: u64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct AseBundleInputs {
    pub cds_fa: HashedPath,
    pub genome_fa: HashedPath,
    pub vcf: HashedPath,
    pub sample_sheet: HashedPath,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct HaplotypeTranscriptome {
    pub name: String,
    pub path: String,
    pub sha256: String,
    pub n_variants_applied: usize,
    pub n_transcripts_modified: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct VariantApplicationSummary {
    pub total_in_vcf: usize,
    pub applied_hap0: usize,
    pub applied_hap1: usize,
    pub skipped_unphased: usize,
    pub skipped_structural: usize,
    pub skipped_multiallelic: usize,
    pub skipped_missing_genotype: usize,
    pub skipped_indel_too_large: usize,
    pub skipped_outside_cds: usize,
    pub skipped_spans_exon_boundary: usize,
    pub skipped_ref_mismatch: usize,
    pub skipped_in_cis_overlap: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct AseSampleSummary {
    pub sample_id: String,
    pub mapping_rate_hap0: f64,
    pub mapping_rate_hap1: f64,
    pub library_size_hap0: u64,
    pub library_size_hap1: u64,
    /// `true` when |mapping_rate_hap0 - mapping_rate_hap1| exceeds
    /// the configured threshold. See `scratch/ase_spec.md` §
    /// "Mapping rate warning".
    pub asymmetry_flag: bool,
    pub fastp: FastpSummary,
}

impl AseBundle {
    pub fn write(&self, path: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(self)?;
        let mut f = File::create(path)?;
        f.write_all(json.as_bytes())?;
        f.write_all(b"\n")?;
        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn make_hashed(p: &str) -> HashedPath {
        HashedPath {
            path: p.to_string(),
            sha256: "a".repeat(64),
        }
    }

    fn make_fastp() -> FastpSummary {
        FastpSummary {
            reads_before_filtering: 10_000_000,
            reads_after_filtering: 9_800_000,
            reads_passed_pct: 98.0,
            q30_rate_before: 0.95,
            q30_rate_after: 0.97,
            adapter_trimmed_reads: Some(1_000),
            adapter_trimmed_bases: Some(20_000),
            duplication_rate: 0.05,
            insert_size_peak: Some(200),
        }
    }

    #[test]
    fn bundle_roundtrips_through_serde() {
        let b = AseBundle {
            version: "0.5.0".to_string(),
            run_timestamp: "2026-04-24T17:30:00Z".to_string(),
            command: "myconote-cli ase cds.fa --vcf v.vcf.gz --samples s.tsv --genome g.fa"
                .to_string(),
            external_tools: ExternalTools {
                fastp: "1.3.2".to_string(),
                salmon: "1.11.4".to_string(),
            },
            inputs: AseBundleInputs {
                cds_fa: make_hashed("cds.fa"),
                genome_fa: make_hashed("genome.fa"),
                vcf: make_hashed("v.vcf.gz"),
                sample_sheet: make_hashed("s.tsv"),
            },
            haplotypes: vec!["hap0".to_string(), "hap1".to_string()],
            personalized_transcriptomes: vec![
                HaplotypeTranscriptome {
                    name: "hap0".to_string(),
                    path: "ase_out/cds_hap0.fa".to_string(),
                    sha256: "b".repeat(64),
                    n_variants_applied: 5123,
                    n_transcripts_modified: 2100,
                },
                HaplotypeTranscriptome {
                    name: "hap1".to_string(),
                    path: "ase_out/cds_hap1.fa".to_string(),
                    sha256: "c".repeat(64),
                    n_variants_applied: 5198,
                    n_transcripts_modified: 2111,
                },
            ],
            variants: VariantApplicationSummary {
                total_in_vcf: 14502,
                applied_hap0: 5123,
                applied_hap1: 5198,
                skipped_unphased: 0,
                skipped_structural: 83,
                skipped_multiallelic: 12,
                skipped_missing_genotype: 0,
                skipped_indel_too_large: 44,
                skipped_outside_cds: 2341,
                skipped_spans_exon_boundary: 13,
                skipped_ref_mismatch: 2,
                skipped_in_cis_overlap: 4,
            },
            samples: vec![AseSampleSummary {
                sample_id: "WT_rep1".to_string(),
                mapping_rate_hap0: 0.89,
                mapping_rate_hap1: 0.88,
                library_size_hap0: 12_345_678,
                library_size_hap1: 12_234_567,
                asymmetry_flag: false,
                fastp: make_fastp(),
            }],
            seed: 42,
        };

        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("ase_bundle.json");
        b.write(&p).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let parsed: AseBundle = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed.version, "0.5.0");
        assert_eq!(parsed.haplotypes, vec!["hap0", "hap1"]);
        assert_eq!(parsed.personalized_transcriptomes.len(), 2);
        assert_eq!(parsed.variants.total_in_vcf, 14502);
        assert_eq!(parsed.variants.applied_hap0, 5123);
        assert_eq!(parsed.samples[0].sample_id, "WT_rep1");
        assert!(!parsed.samples[0].asymmetry_flag);
    }

    #[test]
    fn bundle_json_shape_lists_all_skip_categories() {
        // Schema check: the JSON carries every enumerated skip
        // category so downstream consumers (e.g. a future explain
        // stage for ASE) can count without re-deriving the list.
        let b = AseBundle {
            version: "0.5.0".to_string(),
            run_timestamp: "t".to_string(),
            command: "cmd".to_string(),
            external_tools: ExternalTools {
                fastp: "x".to_string(),
                salmon: "y".to_string(),
            },
            inputs: AseBundleInputs {
                cds_fa: make_hashed("c"),
                genome_fa: make_hashed("g"),
                vcf: make_hashed("v"),
                sample_sheet: make_hashed("s"),
            },
            haplotypes: vec!["hap0".to_string(), "hap1".to_string()],
            personalized_transcriptomes: vec![],
            variants: VariantApplicationSummary {
                total_in_vcf: 0,
                applied_hap0: 0,
                applied_hap1: 0,
                skipped_unphased: 0,
                skipped_structural: 0,
                skipped_multiallelic: 0,
                skipped_missing_genotype: 0,
                skipped_indel_too_large: 0,
                skipped_outside_cds: 0,
                skipped_spans_exon_boundary: 0,
                skipped_ref_mismatch: 0,
                skipped_in_cis_overlap: 0,
            },
            samples: vec![],
            seed: 0,
        };
        let json = serde_json::to_string_pretty(&b).unwrap();
        for key in [
            "skipped_unphased",
            "skipped_structural",
            "skipped_multiallelic",
            "skipped_missing_genotype",
            "skipped_indel_too_large",
            "skipped_outside_cds",
            "skipped_spans_exon_boundary",
            "skipped_ref_mismatch",
            "skipped_in_cis_overlap",
        ] {
            assert!(json.contains(key), "bundle JSON missing field `{key}`");
        }
    }
}
