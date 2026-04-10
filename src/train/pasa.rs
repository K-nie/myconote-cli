use super::TrainConfig;
/// PASA (Program to Assemble Spliced Alignments) wrapper
///
/// PASA aligns Trinity transcripts to the genome, assembles overlapping
/// alignments, and builds a database of transcript assemblies.
/// High-confidence complete models (with start + stop codon) are then
/// used as the gold-standard training set for Augustus/SNAP.
use crate::utils::error::{MycoNoteError, Result};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

pub fn pasa_available() -> bool {
    Command::new("which")
        .arg("Launch_PASA_pipeline.pl")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

pub fn run_pasa(
    config: &TrainConfig,
    transcripts: &Path,
    db_path: &Path,
    gff3_output: &Path,
) -> Result<usize> {
    if !pasa_available() {
        return Err(MycoNoteError::ExternalTool(
            "PASA (Launch_PASA_pipeline.pl) not found.\n  \
             Install: conda install -c bioconda pasa"
                .to_string(),
        ));
    }

    let pasa_dir = config.out_dir.join("pasa_run");
    std::fs::create_dir_all(&pasa_dir).map_err(MycoNoteError::Io)?;

    // Write PASA config file
    let conf_path = pasa_dir.join("pasa.config");
    write_pasa_config(&conf_path, db_path, config.max_intron)?;

    // Run PASA pipeline
    let status = Command::new("Launch_PASA_pipeline.pl")
        .arg("--config")
        .arg(&conf_path)
        .arg("--genome")
        .arg(&config.masked_fasta)
        .arg("--transcripts")
        .arg(transcripts)
        .arg("--CPU")
        .arg(config.threads.to_string())
        .arg("--ALIGNERS")
        .arg("minimap2")
        .arg("--stringent_alignment_overlap")
        .arg("30.0")
        .arg("--transcript_db")
        .arg(db_path)
        .current_dir(&pasa_dir)
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("PASA: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "PASA pipeline exited with non-zero status".to_string(),
        ));
    }

    // Find and copy the GFF3 output
    let gff3_src = find_pasa_gff3(&pasa_dir)?;
    std::fs::copy(&gff3_src, gff3_output).map_err(MycoNoteError::Io)?;

    // Count assemblies
    let count = count_gff3_genes(gff3_output);
    Ok(count)
}

fn write_pasa_config(path: &Path, db_path: &Path, max_intron: usize) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "# PASA configuration file").map_err(MycoNoteError::Io)?;
    writeln!(f, "DATABASE={}", db_path.display()).map_err(MycoNoteError::Io)?;
    writeln!(f, "MAX_INTRON_LENGTH={}", max_intron).map_err(MycoNoteError::Io)?;
    writeln!(f, "MYSQLDB=sqlite3").map_err(MycoNoteError::Io)?;
    Ok(())
}

fn find_pasa_gff3(pasa_dir: &Path) -> Result<PathBuf> {
    if let Ok(entries) = std::fs::read_dir(pasa_dir) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.extension().and_then(|e| e.to_str()) == Some("gff3") {
                return Ok(p);
            }
        }
    }
    Err(MycoNoteError::ExternalTool(
        "PASA GFF3 output file not found".to_string(),
    ))
}

fn count_gff3_genes(gff3: &Path) -> usize {
    use std::io::{BufRead, BufReader};
    let Ok(f) = std::fs::File::open(gff3) else {
        return 0;
    };
    BufReader::new(f)
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            fields.get(2).map(|t| *t == "gene").unwrap_or(false)
        })
        .count()
}
