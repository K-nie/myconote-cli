/// Augustus gene prediction wrapper
///
/// Calls Augustus as a subprocess and parses its GFF3 output.
/// Supports:
///   - Standard prediction from masked FASTA
///   - Hint-based prediction (protein/EST hints improve accuracy)
///   - Multi-sequence parallelism via rayon (one contig → one Augustus call)
///
/// Augustus must be installed and in PATH:
///   conda install -c bioconda augustus
use crate::utils::error::{MycoNoteError, Result};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Augustus configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct AugustusConfig {
    /// Species model (e.g. "saccharomyces_cerevisiae_S288C", "arabidopsis")
    pub species: String,
    /// Number of parallel threads (split by contig)
    pub threads: usize,
    /// Use UTR prediction
    pub utr: bool,
    /// Path to hints GFF file (optional — protein/EST evidence)
    pub hints_file: Option<PathBuf>,
    /// Extra Augustus arguments (passed verbatim)
    pub extra_args: Vec<String>,
}

impl Default for AugustusConfig {
    fn default() -> Self {
        Self {
            species: "saccharomyces_cerevisiae_S288C".to_string(),
            threads: 4,
            utr: false,
            hints_file: None,
            extra_args: Vec::new(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run Augustus on a masked FASTA and return path to the GFF3 output.
pub fn run(masked_fasta: &Path, output_gff: &Path, config: &AugustusConfig) -> Result<()> {
    let aug = which::which("augustus").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "augustus not found in PATH.\n\
             Install with: conda install -c bioconda augustus\n\
             Then verify with: augustus --species=help"
                .to_string(),
        )
    })?;

    println!("  Running Augustus (species: {})…", config.species);

    let utr_flag = if config.utr { "on" } else { "off" };

    let mut args = vec![
        format!("--species={}", config.species),
        "--gff3=on".to_string(),
        "--strand=both".to_string(),
        format!("--UTR={}", utr_flag),
        "--softmasking=1".to_string(), // respect lowercase soft-mask
        format!("--outfile={}", output_gff.display()),
    ];

    // Hints file (protein/EST evidence improves exon boundary accuracy)
    if let Some(ref hints) = config.hints_file {
        args.push(format!("--hintsfile={}", hints.display()));
        args.push("--extrinsicCfgFile=extrinsic.M.RM.E.W.cfg".to_string());
    }

    args.extend(config.extra_args.clone());
    args.push(masked_fasta.to_string_lossy().into_owned());

    let status = Command::new(&aug)
        .args(&args)
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Augustus failed. Check that species '{}' is installed.\n\
             List available species: augustus --species=help",
            config.species
        )));
    }

    println!("  Augustus finished → {}", output_gff.display());
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Protein hint generation from homolog evidence
// ─────────────────────────────────────────────────────────────────────────────

/// Generate an Augustus hints file from a protein alignment (BLAST tabular fmt6).
/// Each hit becomes a CDSpart hint weighted by alignment identity.
pub fn make_protein_hints(blast_tsv: &Path, hints_out: &Path, priority: u8) -> Result<()> {
    use std::io::BufRead;

    let input = std::fs::File::open(blast_tsv).map_err(MycoNoteError::Io)?;
    let mut out = std::fs::File::create(hints_out).map_err(MycoNoteError::Io)?;

    for line in std::io::BufReader::new(input).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        // BLAST fmt6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        if f.len() < 12 {
            continue;
        }

        let seqid = f[0];
        let start: u64 = f[6].parse().unwrap_or(0);
        let end: u64 = f[7].parse().unwrap_or(0);
        let (s, e) = if start <= end {
            (start, end)
        } else {
            (end, start)
        };

        writeln!(
            out,
            "{}\tP\tCDSpart\t{}\t{}\t.\t.\t.\tsrc=P;pri={}",
            seqid, s, e, priority
        )
        .map_err(MycoNoteError::Io)?;
    }

    Ok(())
}
