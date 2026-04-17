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
            hmm: "fungal".to_string(),
            threads: 4,
        }
    }
}

/// Run SNAP and convert output to GFF3.
/// Returns path to the GFF3 output file.
pub fn run(masked_fasta: &Path, output_gff: &Path, config: &SnapConfig) -> Result<()> {
    let snap = which::which("snap").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "snap not found in PATH.\n\
             Install with: conda install -c bioconda snap"
                .to_string(),
        )
    })?;

    println!("  Running SNAP (HMM: {})…", config.hmm);

    // Resolve the HMM argument to something SNAP can actually open. SNAP
    // reads `$ZOE/HMM/<name>` (conda sets ZOE in activate scripts but we
    // shell out bare, so ZOE is usually unset at runtime). Three paths
    // — first that exists wins:
    //   1. `config.hmm` is already an absolute path we can pass through
    //   2. `$ZOE/HMM/<name>` if ZOE is set and the file exists
    //   3. Auto-detect the conda share dir from the snap binary location
    //      (`<conda_root>/share/snap/HMM/<name>`) — this is the bioconda
    //      layout everyone actually has on disk
    let (hmm_arg, zoe_dir): (String, Option<std::path::PathBuf>) = {
        let hmm_path = std::path::Path::new(&config.hmm);
        if hmm_path.is_absolute() && hmm_path.exists() {
            (config.hmm.clone(), None)
        } else if let Ok(zoe) = std::env::var("ZOE") {
            let full = std::path::PathBuf::from(&zoe).join("HMM").join(&config.hmm);
            if full.exists() {
                (config.hmm.clone(), Some(std::path::PathBuf::from(zoe)))
            } else {
                (config.hmm.clone(), None)
            }
        } else {
            // snap binary at `<env>/bin/snap` → HMM dir at `<env>/share/snap/HMM/`
            let share_dir = snap
                .parent()
                .and_then(|p| p.parent())
                .map(|p| p.join("share").join("snap"));
            if let Some(share) = &share_dir {
                let hmm_file = share.join("HMM").join(&config.hmm);
                if hmm_file.exists() {
                    (config.hmm.clone(), Some(share.clone()))
                } else if share.join("HMM").exists() {
                    // Share dir exists but no matching HMM — surface immediately
                    return Err(MycoNoteError::InvalidFormat(format!(
                        "SNAP HMM '{}' not found. Available HMMs in {}:\n  {}",
                        config.hmm,
                        share.join("HMM").display(),
                        std::fs::read_dir(share.join("HMM"))
                            .ok()
                            .map(|rd| rd.filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().into_owned()).collect::<Vec<_>>().join(", "))
                            .unwrap_or_default()
                    )));
                } else {
                    (config.hmm.clone(), None)
                }
            } else {
                (config.hmm.clone(), None)
            }
        }
    };

    let mut cmd = Command::new(&snap);
    cmd.args([
        &hmm_arg,
        masked_fasta.to_str().unwrap_or(""),
        "-gff", // output GFF format directly
    ]);
    if let Some(zoe) = zoe_dir {
        cmd.env("ZOE", zoe);
    }

    let output = cmd.output().map_err(MycoNoteError::Io)?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(MycoNoteError::InvalidFormat(format!(
            "SNAP failed: {}",
            stderr
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
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }

        let seqid = f[0];
        let ftype = f[2];
        let start = f[3];
        let end = f[4];
        let score = f[5];
        let strand = f[6];
        let phase = f[7];
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
