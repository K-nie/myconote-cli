/// Ploidy awareness for gene prediction
///
/// Detects and handles polyploid genomes by:
///   1. Estimating ploidy from k-mer frequency distributions
///   2. Adjusting gene prediction parameters for allelic variation
///   3. Filtering/flagging allelic duplicates in the final gene set
///   4. Supporting haplotype-aware assembly inputs
///
/// This is myconote-cli's own implementation — completely independent
/// of any other annotation pipeline.
use crate::utils::error::{MycoNoteError, Result};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Ploidy configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Ploidy {
    Haploid,
    Diploid,
    Triploid,
    Tetraploid,
    /// User-specified or auto-detected higher ploidy
    Polyploid(u8),
}

impl Ploidy {
    pub fn from_n(n: u8) -> Self {
        match n {
            0 | 1 => Self::Haploid,
            2 => Self::Diploid,
            3 => Self::Triploid,
            4 => Self::Tetraploid,
            _ => Self::Polyploid(n),
        }
    }

    pub fn level(&self) -> u8 {
        match self {
            Self::Haploid => 1,
            Self::Diploid => 2,
            Self::Triploid => 3,
            Self::Tetraploid => 4,
            Self::Polyploid(n) => *n,
        }
    }

    pub fn display_name(&self) -> String {
        match self {
            Self::Haploid => "haploid (1n)".to_string(),
            Self::Diploid => "diploid (2n)".to_string(),
            Self::Triploid => "triploid (3n)".to_string(),
            Self::Tetraploid => "tetraploid (4n)".to_string(),
            Self::Polyploid(n) => format!("polyploid ({}n)", n),
        }
    }

    pub fn is_polyploid(&self) -> bool {
        self.level() > 1
    }
}

#[derive(Debug, Clone)]
pub struct PloidyConfig {
    /// User-specified ploidy (None = auto-detect)
    pub ploidy: Option<Ploidy>,
    /// Minimum sequence identity to consider allelic (0-100)
    pub allelic_identity_threshold: f64,
    /// Minimum alignment coverage to consider allelic (0-1)
    pub allelic_coverage_threshold: f64,
    /// Whether to collapse allelic duplicates
    pub collapse_alleles: bool,
    /// Report allelic pairs without collapsing
    pub report_alleles: bool,
}

impl Default for PloidyConfig {
    fn default() -> Self {
        Self {
            ploidy: None,
            allelic_identity_threshold: 95.0,
            allelic_coverage_threshold: 0.80,
            collapse_alleles: false,
            report_alleles: true,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ploidy estimation from assembly metrics
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug)]
pub struct PloidyEstimate {
    pub estimated_ploidy: Ploidy,
    pub genome_size_bp: u64,
    pub expected_haploid_size: u64,
    pub ratio: f64,
    pub confidence: PloidyConfidence,
    pub method: String,
}

#[derive(Debug, Clone, Copy)]
pub enum PloidyConfidence {
    High,
    Medium,
    Low,
}

/// Estimate ploidy from assembly size vs. expected haploid genome size.
///
/// This is a simple heuristic: if the assembly is ~2x the expected size,
/// it's likely diploid. Works well for fungi where haploid sizes are
/// well-characterized.
pub fn estimate_ploidy_from_size(
    genome_fasta: &Path,
    expected_haploid_mbp: f64,
) -> Result<PloidyEstimate> {
    // Count total bases in the FASTA
    let file = std::fs::File::open(genome_fasta).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut total_bp: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();
        if !trimmed.starts_with('>') && !trimmed.is_empty() {
            total_bp += trimmed.len() as u64;
        }
    }

    let expected_haploid = (expected_haploid_mbp * 1_000_000.0) as u64;
    let ratio = total_bp as f64 / expected_haploid as f64;

    let (ploidy, confidence) = if ratio < 1.3 {
        (Ploidy::Haploid, PloidyConfidence::High)
    } else if ratio < 1.7 {
        (Ploidy::Haploid, PloidyConfidence::Low) // borderline
    } else if ratio < 2.3 {
        (Ploidy::Diploid, PloidyConfidence::High)
    } else if ratio < 2.7 {
        (Ploidy::Diploid, PloidyConfidence::Medium)
    } else if ratio < 3.3 {
        (Ploidy::Triploid, PloidyConfidence::Medium)
    } else if ratio < 4.5 {
        (Ploidy::Tetraploid, PloidyConfidence::Medium)
    } else {
        (
            Ploidy::Polyploid((ratio.round() as u8).min(8)),
            PloidyConfidence::Low,
        )
    };

    Ok(PloidyEstimate {
        estimated_ploidy: ploidy,
        genome_size_bp: total_bp,
        expected_haploid_size: expected_haploid,
        ratio,
        confidence,
        method: "assembly_size_ratio".to_string(),
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Allelic duplicate detection
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct AllelicPair {
    pub gene_a: String,
    pub gene_b: String,
    pub identity: f64,
    pub coverage: f64,
}

/// Detect allelic duplicates using self-alignment of predicted proteins.
/// Runs MMseqs2 or DIAMOND easy-search of proteins against themselves,
/// then identifies high-identity pairs as putative alleles.
pub fn detect_allelic_duplicates(
    proteins_fa: &Path,
    out_dir: &Path,
    config: &PloidyConfig,
    threads: usize,
) -> Result<Vec<AllelicPair>> {
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let self_hits = out_dir.join("self_alignment.tsv");

    // Use MMseqs2 if available, else DIAMOND
    let mmseqs_ok = Command::new("mmseqs")
        .arg("version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);

    if mmseqs_ok {
        let tmp = out_dir.join("tmp_mmseqs");
        std::fs::create_dir_all(&tmp).map_err(MycoNoteError::Io)?;

        let status = Command::new("mmseqs")
            .args([
                "easy-search",
                proteins_fa.to_str().unwrap_or(""),
                proteins_fa.to_str().unwrap_or(""),
                self_hits.to_str().unwrap_or(""),
                tmp.to_str().unwrap_or(""),
                "--min-seq-id",
                &format!("{}", config.allelic_identity_threshold / 100.0),
                "-c",
                &format!("{}", config.allelic_coverage_threshold),
                "--threads",
                &threads.to_string(),
                "--format-output",
                "query,target,pident,qcov",
            ])
            .status()
            .map_err(|e| MycoNoteError::ExternalTool(format!("mmseqs self-search: {}", e)))?;

        if !status.success() {
            return Err(MycoNoteError::ExternalTool(
                "MMseqs2 self-alignment failed".to_string(),
            ));
        }
    } else {
        // Fallback: DIAMOND
        let diamond_ok = Command::new("diamond")
            .arg("version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);

        if !diamond_ok {
            return Err(MycoNoteError::ExternalTool(
                "Neither mmseqs nor diamond found for allelic detection".to_string(),
            ));
        }

        let db = out_dir.join("self.dmnd");
        Command::new("diamond")
            .args([
                "makedb",
                "--in",
                proteins_fa.to_str().unwrap_or(""),
                "--db",
                db.to_str().unwrap_or(""),
            ])
            .status()
            .map_err(|e| MycoNoteError::ExternalTool(format!("diamond makedb: {}", e)))?;

        Command::new("diamond")
            .args([
                "blastp",
                "--query",
                proteins_fa.to_str().unwrap_or(""),
                "--db",
                db.to_str().unwrap_or(""),
                "--out",
                self_hits.to_str().unwrap_or(""),
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "qcovhsp",
                "--id",
                &format!("{}", config.allelic_identity_threshold),
                "--query-cover",
                &format!("{}", config.allelic_coverage_threshold * 100.0),
                "--threads",
                &threads.to_string(),
            ])
            .status()
            .map_err(|e| MycoNoteError::ExternalTool(format!("diamond blastp: {}", e)))?;
    }

    // Parse self-alignment hits
    let pairs = parse_allelic_hits(&self_hits, config)?;

    // Write allelic pairs report
    if config.report_alleles && !pairs.is_empty() {
        let report = out_dir.join("allelic_pairs.tsv");
        write_allelic_report(&report, &pairs)?;
    }

    Ok(pairs)
}

fn parse_allelic_hits(hits_file: &Path, config: &PloidyConfig) -> Result<Vec<AllelicPair>> {
    let file = std::fs::File::open(hits_file).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);

    let mut pairs: Vec<AllelicPair> = Vec::new();
    let mut seen: std::collections::HashSet<(String, String)> = std::collections::HashSet::new();

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let cols: Vec<&str> = line.trim().split('\t').collect();
        if cols.len() < 4 {
            continue;
        }

        let query = cols[0];
        let target = cols[1];

        // Skip self-hits
        if query == target {
            continue;
        }

        let identity: f64 = cols[2].parse().unwrap_or(0.0);
        let coverage: f64 = cols[3].parse().unwrap_or(0.0);

        if identity < config.allelic_identity_threshold {
            continue;
        }
        if coverage < config.allelic_coverage_threshold * 100.0 {
            continue;
        }

        // Deduplicate (A,B) and (B,A)
        let pair_key = if query < target {
            (query.to_string(), target.to_string())
        } else {
            (target.to_string(), query.to_string())
        };

        if seen.contains(&pair_key) {
            continue;
        }
        seen.insert(pair_key);

        pairs.push(AllelicPair {
            gene_a: query.to_string(),
            gene_b: target.to_string(),
            identity,
            coverage,
        });
    }

    Ok(pairs)
}

fn write_allelic_report(path: &Path, pairs: &[AllelicPair]) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "gene_a\tgene_b\tidentity\tcoverage").map_err(MycoNoteError::Io)?;
    for pair in pairs {
        writeln!(
            f,
            "{}\t{}\t{:.1}\t{:.1}",
            pair.gene_a, pair.gene_b, pair.identity, pair.coverage
        )
        .map_err(MycoNoteError::Io)?;
    }
    Ok(())
}

/// Adjust prediction parameters for polyploid genomes.
pub fn adjust_for_ploidy(ploidy: &Ploidy) -> PloidyAdjustments {
    match ploidy {
        Ploidy::Haploid => PloidyAdjustments {
            overlap_tolerance: 0.1,
            min_cds_identity_collapse: 99.0,
            expect_allelic_pairs: false,
            report_note: None,
        },
        Ploidy::Diploid => PloidyAdjustments {
            overlap_tolerance: 0.3,
            min_cds_identity_collapse: 95.0,
            expect_allelic_pairs: true,
            report_note: Some(
                "Diploid assembly: allelic duplicates may inflate gene count".to_string(),
            ),
        },
        _ => PloidyAdjustments {
            overlap_tolerance: 0.5,
            min_cds_identity_collapse: 90.0,
            expect_allelic_pairs: true,
            report_note: Some(format!(
                "Polyploid assembly ({}): expect {}-fold gene duplication from alleles",
                ploidy.display_name(),
                ploidy.level()
            )),
        },
    }
}

#[derive(Debug)]
pub struct PloidyAdjustments {
    /// How much overlap to tolerate between gene models (fraction)
    pub overlap_tolerance: f64,
    /// Identity threshold for collapsing allelic CDS
    pub min_cds_identity_collapse: f64,
    /// Whether to expect allelic pairs
    pub expect_allelic_pairs: bool,
    /// Note for reproducibility report
    pub report_note: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ploidy_from_n() {
        assert_eq!(Ploidy::from_n(1).level(), 1);
        assert_eq!(Ploidy::from_n(2).level(), 2);
        assert_eq!(Ploidy::from_n(4).level(), 4);
        assert_eq!(Ploidy::from_n(6).level(), 6);
    }

    #[test]
    fn test_ploidy_display() {
        assert_eq!(Ploidy::Haploid.display_name(), "haploid (1n)");
        assert_eq!(Ploidy::Diploid.display_name(), "diploid (2n)");
    }

    #[test]
    fn test_is_polyploid() {
        assert!(!Ploidy::Haploid.is_polyploid());
        assert!(Ploidy::Diploid.is_polyploid());
        assert!(Ploidy::Polyploid(6).is_polyploid());
    }
}
