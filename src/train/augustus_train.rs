/// Augustus training wrapper

use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

pub fn train_augustus(
    training_gff3:  &Path,
    genome_fasta:   &Path,
    species_name:   &str,
    out_dir:        &Path,
    threads:        usize,
) -> Result<PathBuf> {
    // Check for augustus_species_dir script
    let has_augustus = Command::new("which").arg("augustus")
        .output().map(|o| o.status.success()).unwrap_or(false);
    let has_etraining = Command::new("which").arg("etraining")
        .output().map(|o| o.status.success()).unwrap_or(false);

    if !has_augustus || !has_etraining {
        return Err(MycoNoteError::ExternalTool(
            "Augustus training tools (augustus, etraining) not found.\n  \
             Install: conda install -c bioconda augustus".to_string()
        ));
    }

    let species_dir = out_dir.join("augustus_training").join(species_name);
    std::fs::create_dir_all(&species_dir).map_err(MycoNoteError::Io)?;

    // 1. Create new species in Augustus config
    let _ = Command::new("new_species.pl")
        .arg(&format!("--species={}", species_name))
        .status();

    // 2. Convert GFF3 to Augustus training format (genbank-style)
    let gb_train = out_dir.join("training.gb");
    let status = Command::new("gff2gbSmallDNA.pl")
        .arg(training_gff3)
        .arg(genome_fasta)
        .arg("1000")   // flanking region
        .arg(&gb_train)
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("gff2gbSmallDNA.pl: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "gff2gbSmallDNA.pl failed (needed to convert training data for Augustus)".to_string()
        ));
    }

    // 3. Split into train/test sets (90/10 split)
    let gb_test  = out_dir.join("test.gb");
    let _ = Command::new("randomSplit.pl")
        .arg(&gb_train)
        .arg("100")  // 100 genes for test set
        .status();

    // 4. Run etraining
    let status = Command::new("etraining")
        .arg(&format!("--species={}", species_name))
        .arg(&gb_train)
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("etraining: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool("etraining failed".to_string()));
    }

    // 5. Evaluate accuracy on test set
    if gb_test.exists() {
        let eval_out = out_dir.join("augustus_eval.txt");
        let eval_f = std::fs::File::create(&eval_out).map_err(MycoNoteError::Io)?;
        let _ = Command::new("augustus")
            .arg(&format!("--species={}", species_name))
            .arg(&gb_test)
            .stdout(std::process::Stdio::from(eval_f))
            .status();

        // Parse and print gene-level sensitivity/specificity
        if let Ok(content) = std::fs::read_to_string(&eval_out) {
            for line in content.lines() {
                if line.contains("gene level") || line.contains("sensitivity") || line.contains("specificity") {
                    println!("  Augustus eval: {}", line.trim());
                }
            }
        }
    }

    // 6. Optimize (optional — can be very slow, skip by default)
    // optimize_augustus.pl is omitted here; users can run it manually

    let _ = threads; // threads not directly used by etraining (single-threaded)
    println!("  Augustus species '{}' trained successfully", species_name);

    Ok(species_dir)
}
