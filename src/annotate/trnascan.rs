/// tRNAscan-SE integration for tRNA gene prediction
///
/// Runs tRNAscan-SE on a genome FASTA and produces GFF3 features
/// for all predicted tRNA genes. Supports both eukaryotic and
/// organellar (mitochondrial) scanning modes.
///
/// This is myconote-cli's own independent implementation wrapping
/// the tRNAscan-SE binary. It does NOT share any code with other
/// annotation pipelines.
///
/// Install: conda install -c bioconda trnascan-se
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct TrnaScanConfig {
    /// Search mode: "eukaryotic", "mitochondrial", "general"
    pub mode: String,
    /// Number of threads
    pub threads: usize,
    /// Minimum score to report (default: 20.0)
    pub min_score: f64,
    /// Also predict pseudogenes
    pub pseudogenes: bool,
}

impl Default for TrnaScanConfig {
    fn default() -> Self {
        Self {
            mode: "eukaryotic".to_string(),
            threads: 4,
            min_score: 20.0,
            pseudogenes: false,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// tRNA prediction result
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct TrnaGene {
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub amino_acid: String,
    pub anticodon: String,
    pub score: f64,
    pub is_pseudo: bool,
    pub intron_start: Option<u64>,
    pub intron_end: Option<u64>,
}

#[derive(Debug, Default)]
pub struct TrnaScanResult {
    pub trnas: Vec<TrnaGene>,
    pub total: usize,
    pub by_amino_acid: HashMap<String, usize>,
    pub pseudogenes: usize,
    pub with_introns: usize,
}

impl TrnaScanResult {
    pub fn summarize(&self) {
        println!("  tRNA genes found: {}", self.total);
        if self.pseudogenes > 0 {
            println!("  tRNA pseudogenes: {}", self.pseudogenes);
        }
        if self.with_introns > 0 {
            println!("  tRNAs with introns: {}", self.with_introns);
        }

        // Print amino acid distribution
        let mut aa_counts: Vec<_> = self.by_amino_acid.iter().collect();
        aa_counts.sort_by(|a, b| b.1.cmp(a.1));
        if !aa_counts.is_empty() {
            let summary: Vec<String> = aa_counts
                .iter()
                .map(|(aa, n)| format!("{}:{}", aa, n))
                .collect();
            println!("  Anticodon distribution: {}", summary.join(", "));
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Check availability
// ─────────────────────────────────────────────────────────────────────────────

pub fn trnascan_available() -> bool {
    Command::new("tRNAscan-SE")
        .arg("--version")
        .output()
        .map(|o| o.status.success() || !o.stderr.is_empty()) // tRNAscan prints version to stderr
        .unwrap_or(false)
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run tRNAscan-SE on a genome FASTA and return predicted tRNA genes.
pub fn run_trnascan(
    genome_fasta: &Path,
    out_dir: &Path,
    config: &TrnaScanConfig,
) -> Result<TrnaScanResult> {
    if !trnascan_available() {
        return Err(MycoNoteError::ExternalTool(
            "tRNAscan-SE not found. Install: conda install -c bioconda trnascan-se".to_string(),
        ));
    }

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let raw_out = out_dir.join("trnascan_raw.txt");
    let struct_out = out_dir.join("trnascan_structures.ss");
    let gff_out = out_dir.join("trnascan.gff3");
    let stats_out = out_dir.join("trnascan_stats.txt");

    // Build command
    let mut cmd = Command::new("tRNAscan-SE");

    // Mode selection
    match config.mode.as_str() {
        "mitochondrial" | "mito" => {
            cmd.arg("-M");
        }
        "general" | "bacterial" => {
            cmd.arg("-G");
        }
        _ => {
            cmd.arg("-E");
        } // eukaryotic (default)
    }

    cmd.arg("--thread").arg(config.threads.to_string());
    cmd.arg("-o").arg(&raw_out);
    cmd.arg("-f").arg(&struct_out);
    cmd.arg("-m").arg(&stats_out);

    if config.pseudogenes {
        cmd.arg("--pseudo");
    }

    // GFF3 output
    cmd.arg("--gff").arg(&gff_out);

    // Minimum score
    cmd.arg("-X").arg(config.min_score.to_string());

    cmd.arg(genome_fasta);

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("tRNAscan-SE: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "tRNAscan-SE exited with non-zero status".to_string(),
        ));
    }

    // Parse results
    let result = parse_trnascan_output(&raw_out)?;

    // Write summary
    write_trna_summary(out_dir, &result)?;

    Ok(result)
}

/// Write tRNA predictions as GFF3 features (our own format, independent of
/// tRNAscan-SE's --gff output which may not include all fields we want).
pub fn write_trna_gff3(trnas: &[TrnaGene], output: &Path, locus_prefix: &str) -> Result<usize> {
    let mut f = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(f, "##gff-version 3").map_err(MycoNoteError::Io)?;

    let mut count = 0usize;
    for (i, trna) in trnas.iter().enumerate() {
        let gene_id = format!("{}_tRNA_{:04}", locus_prefix, i + 1);
        let product = format!("tRNA-{}", trna.amino_acid);
        let note = if trna.is_pseudo {
            format!("anticodon:{};pseudo", trna.anticodon)
        } else {
            format!("anticodon:{}", trna.anticodon)
        };

        // Gene feature
        writeln!(
            f,
            "{}\ttRNAscan-SE\tgene\t{}\t{}\t{:.1}\t{}\t.\tID={};Name={};product={};Note={}",
            trna.seqid,
            trna.start,
            trna.end,
            trna.score,
            trna.strand,
            gene_id,
            product,
            product,
            note,
        )
        .map_err(MycoNoteError::Io)?;

        // tRNA feature
        writeln!(
            f,
            "{}\ttRNAscan-SE\ttRNA\t{}\t{}\t{:.1}\t{}\t.\tID={}-tRNA;Parent={};product={}",
            trna.seqid, trna.start, trna.end, trna.score, trna.strand, gene_id, gene_id, product,
        )
        .map_err(MycoNoteError::Io)?;

        // Exon features (handling introns if present)
        if let (Some(int_s), Some(int_e)) = (trna.intron_start, trna.intron_end) {
            // Exon 1: before intron
            if trna.start < int_s {
                writeln!(
                    f,
                    "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tParent={}-tRNA",
                    trna.seqid,
                    trna.start,
                    int_s - 1,
                    trna.strand,
                    gene_id,
                )
                .map_err(MycoNoteError::Io)?;
            }
            // Exon 2: after intron
            if int_e < trna.end {
                writeln!(
                    f,
                    "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tParent={}-tRNA",
                    trna.seqid,
                    int_e + 1,
                    trna.end,
                    trna.strand,
                    gene_id,
                )
                .map_err(MycoNoteError::Io)?;
            }
        } else {
            // Single exon
            writeln!(
                f,
                "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tParent={}-tRNA",
                trna.seqid, trna.start, trna.end, trna.strand, gene_id,
            )
            .map_err(MycoNoteError::Io)?;
        }

        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// Parser for tRNAscan-SE tabular output
// ─────────────────────────────────────────────────────────────────────────────

fn parse_trnascan_output(path: &Path) -> Result<TrnaScanResult> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);

    let mut result = TrnaScanResult::default();

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();

        // Skip header lines (tRNAscan-SE output has 3 header lines)
        if trimmed.is_empty()
            || trimmed.starts_with("Sequence")
            || trimmed.starts_with("Name")
            || trimmed.starts_with('-')
        {
            continue;
        }

        let cols: Vec<&str> = trimmed.split_whitespace().collect();
        if cols.len() < 9 {
            continue;
        }

        // Columns: SeqName  tRNA#  Begin  End  Type  Anti  IntronBegin  IntronEnd  Score  Note
        let seqid = cols[0].to_string();
        let begin: u64 = cols[2].parse().unwrap_or(0);
        let end: u64 = cols[3].parse().unwrap_or(0);
        let amino_acid = cols[4].to_string();
        let anticodon = cols[5].to_string();

        let intron_begin: u64 = cols[6].parse().unwrap_or(0);
        let intron_end: u64 = cols[7].parse().unwrap_or(0);
        let score: f64 = cols[8].parse().unwrap_or(0.0);

        let (start, stop, strand) = if begin <= end {
            (begin, end, '+')
        } else {
            (end, begin, '-')
        };

        let is_pseudo = amino_acid == "Pseudo"
            || amino_acid == "Undet"
            || (cols.len() > 9 && cols[9..].join(" ").contains("pseudo"));

        let (int_s, int_e) = if intron_begin > 0 && intron_end > 0 {
            let (a, b) = if intron_begin <= intron_end {
                (intron_begin, intron_end)
            } else {
                (intron_end, intron_begin)
            };
            (Some(a), Some(b))
        } else {
            (None, None)
        };

        let trna = TrnaGene {
            seqid,
            start: start,
            end: stop,
            strand,
            amino_acid: amino_acid.clone(),
            anticodon,
            score,
            is_pseudo,
            intron_start: int_s,
            intron_end: int_e,
        };

        if is_pseudo {
            result.pseudogenes += 1;
        }
        if int_s.is_some() {
            result.with_introns += 1;
        }

        *result.by_amino_acid.entry(amino_acid).or_insert(0) += 1;
        result.trnas.push(trna);
    }

    result.total = result.trnas.len();
    Ok(result)
}

fn write_trna_summary(out_dir: &Path, result: &TrnaScanResult) -> Result<()> {
    let summary_path = out_dir.join("trna_summary.txt");
    let mut f = std::fs::File::create(&summary_path).map_err(MycoNoteError::Io)?;

    writeln!(f, "myconote tRNAscan-SE Summary").map_err(MycoNoteError::Io)?;
    writeln!(f, "============================").map_err(MycoNoteError::Io)?;
    writeln!(f, "Total tRNA genes: {}", result.total).map_err(MycoNoteError::Io)?;
    writeln!(f, "Pseudogenes:      {}", result.pseudogenes).map_err(MycoNoteError::Io)?;
    writeln!(f, "With introns:     {}", result.with_introns).map_err(MycoNoteError::Io)?;
    writeln!(f).map_err(MycoNoteError::Io)?;
    writeln!(f, "Amino acid distribution:").map_err(MycoNoteError::Io)?;

    let mut aa_counts: Vec<_> = result.by_amino_acid.iter().collect();
    aa_counts.sort_by(|a, b| b.1.cmp(a.1));
    for (aa, count) in &aa_counts {
        writeln!(f, "  {:<6} {}", aa, count).map_err(MycoNoteError::Io)?;
    }

    Ok(())
}
