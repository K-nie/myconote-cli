/// SNAP HMM training wrapper
use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

pub fn snap_available() -> bool {
    Command::new("which")
        .arg("snap")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

pub fn train_snap(training_gff3: &Path, genome_fasta: &Path, out_dir: &Path) -> Result<PathBuf> {
    if !snap_available() {
        return Err(MycoNoteError::ExternalTool(
            "SNAP not found. Install: conda install -c bioconda snap".to_string(),
        ));
    }

    let snap_dir = out_dir.join("snap_training");
    std::fs::create_dir_all(&snap_dir).map_err(MycoNoteError::Io)?;

    // 1. Convert GFF3 → ZFF format (SNAP's native format)
    let zff_path = snap_dir.join("training.ann");
    let dna_path = snap_dir.join("training.dna");

    let zff_f = std::fs::File::create(&zff_path).map_err(MycoNoteError::Io)?;
    let dna_f = std::fs::File::create(&dna_path).map_err(MycoNoteError::Io)?;

    let status = Command::new("gff3_to_zff.pl")
        .arg(training_gff3)
        .arg(genome_fasta)
        .stdout(std::process::Stdio::from(zff_f))
        .status();

    // If the perl script isn't available, use our built-in converter
    if status.map(|s| !s.success()).unwrap_or(true) {
        convert_gff3_to_zff(training_gff3, genome_fasta, &zff_path, &dna_path)?;
    } else {
        // Copy genome as DNA file
        std::fs::copy(genome_fasta, &dna_path).map_err(MycoNoteError::Io)?;
        drop(dna_f);
    }

    // 2. Run SNAP training (forge)
    let status = Command::new("snap")
        .arg("-train")
        .arg(&zff_path)
        .arg(&dna_path)
        .current_dir(&snap_dir)
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("snap -train: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "SNAP training failed".to_string(),
        ));
    }

    // SNAP writes HMM to the current directory
    let hmm_out = snap_dir.join("HMM");
    if !hmm_out.exists() {
        return Err(MycoNoteError::ExternalTool(
            "SNAP training completed but HMM file not found".to_string(),
        ));
    }

    // Rename to something descriptive
    let final_hmm = out_dir.join("snap_trained.hmm");
    std::fs::copy(&hmm_out, &final_hmm).map_err(MycoNoteError::Io)?;

    Ok(final_hmm)
}

/// Minimal GFF3 → ZFF conversion (fallback if gff3_to_zff.pl not available).
/// ZFF format: lines starting with '>' are sequence names;
/// feature lines are tab-separated: type start end strand [optional]
fn convert_gff3_to_zff(gff3: &Path, _genome: &Path, zff_out: &Path, _dna: &Path) -> Result<()> {
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, Write};

    let file = std::fs::File::open(gff3).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut out = std::fs::File::create(zff_out).map_err(MycoNoteError::Io)?;

    // Group exons/CDS by mRNA parent, keyed by seqid
    let mut by_seq: HashMap<String, Vec<(String, u64, u64, char)>> = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let t = line.trim();
        if t.starts_with('#') || t.is_empty() {
            continue;
        }

        let fields: Vec<&str> = t.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        if fields[2] != "CDS" {
            continue;
        }

        let seqid = fields[0].to_string();
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand: char = fields[6].chars().next().unwrap_or('+');
        let feat_type = "Exon".to_string();

        by_seq
            .entry(seqid)
            .or_default()
            .push((feat_type, start, end, strand));
    }

    for (seqid, features) in &by_seq {
        writeln!(out, ">{}", seqid).map_err(MycoNoteError::Io)?;
        for (ftype, start, end, strand) in features {
            writeln!(out, "{}\t{}\t{}\t{}", ftype, start, end, strand)
                .map_err(MycoNoteError::Io)?;
        }
    }

    Ok(())
}
