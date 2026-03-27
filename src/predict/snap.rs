/// SNAP gene predictor wrapper
///
/// SNAP (Semi-HMM-based Nucleic Acid Parser) is a fast ab initio predictor.
/// Used as a secondary predictor whose results are fed into the Evidence
/// Modeler consensus alongside Augustus.
///
/// Install: conda install -c bioconda snap

use crate::utils::error::{MycoNoteError, Result};
use std::path::Path;
use std::process::Command;

#[derive(Debug, Clone)]
pub struct SnapConfig {
    /// Path to SNAP HMM parameter file
    /// Common files: fungal, pombe, worm, fly, human
    pub hmm: String,
    pub threads: usize,
}

impl Default for SnapConfig {
    fn default() -> Self {
        Self {
            hmm:     "fungal".to_string(),
            threads: 4,
        }
    }
}

/// Run SNAP and convert output to GFF3.
/// Returns path to the GFF3 output file.
pub fn run(
    masked_fasta: &Path,
    output_gff:   &Path,
    config:       &SnapConfig,
) -> Result<()> {
    let snap = which::which("snap").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "snap not found in PATH.\n\
             Install with: conda install -c bioconda snap".to_string()
        )
    })?;

    println!("  Running SNAP (HMM: {})…", config.hmm);

    // SNAP writes ZFF format; we capture stdout then convert
    let _snap_out = output_gff.with_extension("zff");

    let output = Command::new(&snap)
        .args([
            &config.hmm,
            masked_fasta.to_str().unwrap_or(""),
            "-gff",   // output GFF format directly
        ])
        .output()
        .map_err(MycoNoteError::Io)?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(MycoNoteError::InvalidFormat(format!(
            "SNAP failed: {}", stderr
        )));
    }

    // Convert SNAP GFF to GFF3 (SNAP outputs GFF2 style)
    let snap_gff2 = String::from_utf8_lossy(&output.stdout);
    let gff3 = convert_snap_to_gff3(&snap_gff2);

    std::fs::write(output_gff, gff3).map_err(MycoNoteError::Io)?;

    println!("  SNAP finished → {}", output_gff.display());
    Ok(())
}

/// Convert SNAP's GFF2-style output to GFF3.
/// SNAP emits lines like:
///   contig1  SNAP  Esngl  100  300  .  +  0  gene.1
fn convert_snap_to_gff3(snap_output: &str) -> String {
    let mut lines = Vec::new();
    lines.push("##gff-version 3".to_string());

    let mut gene_counter = 0u64;
    let mut current_gene: Option<String> = None;
    let mut mrna_counter = 0u64;

    for line in snap_output.lines() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 { continue; }

        let seqid   = f[0];
        let ftype   = f[2];
        let start   = f[3];
        let end     = f[4];
        let score   = f[5];
        let strand  = f[6];
        let phase   = f[7];
        let gene_id = f[8].trim();

        // Map SNAP feature types to GFF3
        let (gff3_type, need_gene) = match ftype {
            "Esngl" | "Einit" | "Eterm" | "Exon" => ("CDS", true),
            "Intron" => ("intron", false),
            _ => (ftype, false),
        };

        if need_gene {
            // Emit a gene/mRNA wrapper on first encounter
            if current_gene.as_deref() != Some(gene_id) {
                gene_counter += 1;
                mrna_counter += 1;
                let gid = format!("SNAP_gene_{}", gene_counter);
                let mid = format!("SNAP_mRNA_{}", mrna_counter);
                lines.push(format!(
                    "{}\tSNAP\tgene\t{}\t{}\t{}\t{}\t.\tID={}",
                    seqid, start, end, score, strand, gid
                ));
                lines.push(format!(
                    "{}\tSNAP\tmRNA\t{}\t{}\t{}\t{}\t.\tID={};Parent={}",
                    seqid, start, end, score, strand, mid, gid
                ));
                current_gene = Some(gene_id.to_string());
            }

            let mid = format!("SNAP_mRNA_{}", mrna_counter);
            lines.push(format!(
                "{}\tSNAP\t{}\t{}\t{}\t{}\t{}\t{}\tParent={}",
                seqid, gff3_type, start, end, score, strand, phase, mid
            ));
        }
    }

    lines.join("\n")
}
