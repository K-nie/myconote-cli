/// RepeatMasker wrapper
///
/// Calls RepeatMasker as a subprocess and parses its .out file to extract
/// masked regions.  RepeatMasker must be installed and in PATH.
///
/// Install: conda install -c bioconda repeatmasker
use super::{MaskConfig, MaskedRegion};
use crate::utils::error::{MycoNoteError, Result};
use std::io::Write;
use std::process::Command;

/// Run RepeatMasker on the given sequences and return masked regions.
pub fn run(sequences: &[(String, String)], config: &MaskConfig) -> Result<Vec<MaskedRegion>> {
    // Check RepeatMasker is available
    let rm_path = which::which("RepeatMasker").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "RepeatMasker not found in PATH.\n\
             Install with: conda install -c bioconda repeatmasker\n\
             Or use --engine minimap2 for the built-in repeat finder."
                .to_string(),
        )
    })?;

    // Write sequences to a temp FASTA
    let tmp_dir = tempfile::TempDir::new().map_err(MycoNoteError::Io)?;
    let fasta_path = tmp_dir.path().join("input.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).map_err(MycoNoteError::Io)?;
        for (id, seq) in sequences {
            writeln!(f, ">{}", id).map_err(MycoNoteError::Io)?;
            for chunk in seq.as_bytes().chunks(60) {
                writeln!(f, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                    .map_err(MycoNoteError::Io)?;
            }
        }
    }

    // Build RepeatMasker arguments
    let mut args: Vec<String> = vec![
        "-pa".into(),
        config.threads.to_string(),
        "-xsmall".into(), // soft-mask (lowercase)
        "-nolow".into(),  // skip low-complexity (handled separately)
        "-dir".into(),
        tmp_dir.path().to_string_lossy().into_owned(),
    ];

    if config.hard_mask {
        // Replace -xsmall with default (hard masking = N)
        args.retain(|a| a != "-xsmall");
    }

    if let Some(ref lib) = config.repeat_lib {
        args.push("-lib".into());
        args.push(lib.to_string_lossy().into_owned());
    } else if let Some(ref species) = config.species {
        args.push("-species".into());
        args.push(species.clone());
    } else {
        // Generic masking without a species database
        args.push("-species".into());
        args.push("fungi".into());
    }

    args.push(fasta_path.to_string_lossy().into_owned());

    println!("  Running RepeatMasker…");
    let status = Command::new(&rm_path)
        .args(&args)
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "RepeatMasker exited with non-zero status. \
             Check that the species database is installed: \
             `RepeatMasker -species fungi -help`"
                .to_string(),
        ));
    }

    // Parse the .out file produced by RepeatMasker
    let out_file = tmp_dir.path().join("input.fa.out");
    if !out_file.exists() {
        return Err(MycoNoteError::InvalidFormat(
            "RepeatMasker did not produce an output file.".to_string(),
        ));
    }

    parse_rm_out(&out_file, config.min_length)
}

/// Parse a RepeatMasker .out file and return masked regions.
///
/// Format (space-delimited, skip first 3 header lines):
///   score  div  del  ins  query  qStart  qEnd  qLeft  strand  repeat  class  rStart  rEnd  rLeft  id
fn parse_rm_out(path: &std::path::Path, min_length: usize) -> Result<Vec<MaskedRegion>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = std::io::BufReader::new(file);
    let mut regions = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        if i < 3 {
            continue;
        } // skip header lines
        let line = line.map_err(MycoNoteError::Io)?;
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 9 {
            continue;
        }

        let seqid = cols[4].to_string();
        let start: usize = cols[5].parse().unwrap_or(0);
        let end: usize = cols[6].parse().unwrap_or(0);

        if end.saturating_sub(start) < min_length {
            continue;
        }

        // RepeatMasker uses 1-based inclusive; convert to 0-based half-open
        regions.push(MaskedRegion {
            seqid,
            start: start.saturating_sub(1),
            end,
        });
    }

    println!(
        "  RepeatMasker: {} repeat regions identified",
        regions.len()
    );
    Ok(regions)
}
