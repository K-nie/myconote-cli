/// Trinity transcript assembler wrapper

use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;
use super::TrainConfig;

pub fn trinity_available() -> bool {
    Command::new("which").arg("Trinity")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

pub fn check_trinity() -> Result<()> {
    if !trinity_available() {
        return Err(MycoNoteError::ExternalTool(
            "Trinity not found in PATH.\n  \
             Install: conda install -c bioconda trinity\n  \
             Or download from: https://github.com/trinityrnaseq/trinityrnaseq".to_string()
        ));
    }
    Ok(())
}

pub fn run_trinity(config: &TrainConfig, output: &Path) -> Result<PathBuf> {
    let trinity_dir = config.out_dir.join("trinity_out");

    let mut cmd = Command::new("Trinity");

    // Input reads
    if !config.left_reads.is_empty() {
        let left: Vec<String> = config.left_reads.iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
        let right: Vec<String> = config.right_reads.iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
        cmd.arg("--left").arg(left.join(","));
        if !right.is_empty() {
            cmd.arg("--right").arg(right.join(","));
            cmd.arg("--seqType").arg("fq");
        }
    } else if !config.single_reads.is_empty() {
        let singles: Vec<String> = config.single_reads.iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
        cmd.arg("--single").arg(singles.join(","));
        cmd.arg("--seqType").arg("fq");
    } else {
        return Err(MycoNoteError::ExternalTool(
            "No RNA-seq reads provided. Use --left/--right or --single.".to_string()
        ));
    }

    // Genome-guided mode
    cmd.arg("--genome_guided_bam").arg(config.out_dir.join("trinity_aligned.bam"));
    cmd.arg("--genome_guided_max_intron").arg(config.max_intron.to_string());

    cmd.arg("--output").arg(&trinity_dir)
       .arg("--CPU").arg(config.threads.to_string())
       .arg("--max_memory").arg(&config.trinity_memory)
       .arg("--full_cleanup");

    if !config.strand.is_empty() {
        cmd.arg("--SS_lib_type").arg(&config.strand);
    }

    let status = cmd.status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("Trinity: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool("Trinity failed".to_string()));
    }

    // Trinity writes to GG.fasta in genome-guided mode
    let gg_fasta = trinity_dir.join("Trinity-GG.fasta");
    if gg_fasta.exists() {
        std::fs::copy(&gg_fasta, output).map_err(MycoNoteError::Io)?;
    } else {
        // de novo output
        let dn_fasta = trinity_dir.join("Trinity.fasta");
        if dn_fasta.exists() {
            std::fs::copy(&dn_fasta, output).map_err(MycoNoteError::Io)?;
        } else {
            return Err(MycoNoteError::ExternalTool(
                "Trinity output FASTA not found".to_string()
            ));
        }
    }

    Ok(output.to_path_buf())
}
