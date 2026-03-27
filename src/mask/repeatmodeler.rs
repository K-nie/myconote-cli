/// RepeatModeler de novo repeat library builder
///
/// Discovers organism-specific repeat families directly from the genome
/// without relying on RepBase. This is critical for non-model organisms
/// where RepBase coverage is poor.
///
/// Workflow:
///   1. `BuildDatabase` — index genome for RepeatModeler
///   2. `RepeatModeler -engine ncbi -pa <threads> -database <db>` — discover repeats
///   3. Output: `<prefix>-families.fa` — de novo repeat library
///   4. Feed library to RepeatMasker for masking
///
/// Install: conda install -c bioconda repeatmodeler

use crate::utils::error::{MycoNoteError, Result};
use super::{MaskConfig, MaskedRegion};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point: build library + mask
// ─────────────────────────────────────────────────────────────────────────────

/// Build a de novo repeat library with RepeatModeler, then mask the genome
/// with RepeatMasker using that library.
/// Returns masked regions.
pub fn run(
    sequences: &[(String, String)],
    config:    &MaskConfig,
) -> Result<Vec<MaskedRegion>> {
    let tmp_dir = tempfile::TempDir::new().map_err(MycoNoteError::Io)?;
    let tmp_path = tmp_dir.path();

    // Write genome FASTA to temp dir
    let genome_fa = tmp_path.join("genome.fa");
    write_fasta_temp(&genome_fa, sequences)?;

    // ── 1. Build RepeatModeler database ───────────────────────────────────────
    let db_name = tmp_path.join("rmdb");
    build_database(&genome_fa, &db_name)?;

    // ── 2. Run RepeatModeler ──────────────────────────────────────────────────
    let repeat_lib = run_repeatmodeler(&db_name, tmp_path, config.threads)?;
    println!("  RepeatModeler complete: {}", repeat_lib.display());

    // ── 3. Run RepeatMasker with the de novo library ──────────────────────────
    let rm_out_dir = tmp_path.join("rm_out");
    std::fs::create_dir_all(&rm_out_dir).map_err(MycoNoteError::Io)?;

    run_repeatmasker_with_lib(&genome_fa, &repeat_lib, &rm_out_dir, config)?;

    // ── 4. Parse RepeatMasker output ──────────────────────────────────────────
    let rm_out_file = rm_out_dir.join("genome.fa.out");
    let regions = if rm_out_file.exists() {
        parse_rm_out(&rm_out_file, config.min_length)?
    } else {
        eprintln!("  ⚠  RepeatMasker .out file not found — no regions masked");
        Vec::new()
    };

    // Optionally copy the repeat library to out_dir for reuse
    if let Some(parent) = config.output.parent() {
        let lib_dest = parent.join(format!(
            "{}-families.fa",
            config.output.file_stem().unwrap_or_default().to_string_lossy()
        ));
        let _ = std::fs::copy(&repeat_lib, &lib_dest);
        println!("  De novo repeat library saved: {}", lib_dest.display());
    }

    println!("  {} masked regions from de novo library", regions.len());
    Ok(regions)
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 1: BuildDatabase
// ─────────────────────────────────────────────────────────────────────────────

fn build_database(genome_fa: &Path, db_name: &Path) -> Result<()> {
    let build_db = which::which("BuildDatabase").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "BuildDatabase not found in PATH.\n\
             Install RepeatModeler: conda install -c bioconda repeatmodeler".to_string()
        )
    })?;

    println!("  Building RepeatModeler database…");

    let status = Command::new(&build_db)
        .args([
            "-name",   db_name.to_str().unwrap_or(""),
            "-engine", "ncbi",
            genome_fa.to_str().unwrap_or(""),
        ])
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "BuildDatabase failed. Check that RepeatModeler is correctly installed.".to_string()
        ));
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 2: RepeatModeler
// ─────────────────────────────────────────────────────────────────────────────

fn run_repeatmodeler(db_name: &Path, work_dir: &Path, threads: usize) -> Result<PathBuf> {
    let repeatmodeler = which::which("RepeatModeler").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "RepeatModeler not found in PATH.\n\
             Install with: conda install -c bioconda repeatmodeler".to_string()
        )
    })?;

    println!("  Running RepeatModeler (this may take 30–120 min for a typical genome)…");

    let status = Command::new(&repeatmodeler)
        .current_dir(work_dir)
        .args([
            "-engine", "ncbi",
            "-pa",     &threads.to_string(),
            "-database", db_name.to_str().unwrap_or(""),
        ])
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "RepeatModeler failed. Check that RECON and RepeatScout are installed.".to_string()
        ));
    }

    // Find the output families FASTA — RepeatModeler creates it as
    // <db_basename>-families.fa in the working directory
    let db_stem = db_name.file_name()
        .unwrap_or_default()
        .to_string_lossy();

    // RepeatModeler v2 puts output in a RM_* subdirectory
    let families_v2 = find_families_v2(work_dir);
    if let Some(p) = families_v2 {
        return Ok(p);
    }

    let families = work_dir.join(format!("{}-families.fa", db_stem));
    if families.exists() {
        return Ok(families);
    }

    Err(MycoNoteError::InvalidFormat(
        "RepeatModeler ran but output families FASTA not found.".to_string()
    ))
}

/// RepeatModeler v2 puts families in a RM_<PID>.<date>/families-classified.fa
fn find_families_v2(work_dir: &Path) -> Option<PathBuf> {
    let entries = std::fs::read_dir(work_dir).ok()?;
    for entry in entries.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with("RM_") {
            let classified = entry.path().join("families-classified.fa");
            if classified.exists() { return Some(classified); }
            let families = entry.path().join("families.fa");
            if families.exists() { return Some(families); }
        }
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 3: RepeatMasker with custom library
// ─────────────────────────────────────────────────────────────────────────────

fn run_repeatmasker_with_lib(
    genome_fa:  &Path,
    repeat_lib: &Path,
    out_dir:    &Path,
    config:     &MaskConfig,
) -> Result<()> {
    let rm = which::which("RepeatMasker").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "RepeatMasker not found in PATH.\n\
             Install with: conda install -c bioconda repeatmasker".to_string()
        )
    })?;

    println!("  Running RepeatMasker with de novo library…");

    let cmd_args = [
        "-lib",  repeat_lib.to_str().unwrap_or(""),
        "-dir",  out_dir.to_str().unwrap_or(""),
        "-pa",   // threads arg comes next (in extra_args)
    ];

    let threads_str = config.threads.to_string();
    // RepeatMasker: -xsmall gives soft-masked (lowercase) output;
    // no flag gives the default hard-masked (N) output.
    let soft_flag = if !config.hard_mask {
        Some("-xsmall".to_string())
    } else {
        None
    };
    let mut extra_args: Vec<&str> = vec![&threads_str];
    if let Some(ref flag) = soft_flag {
        extra_args.push(flag.as_str());
    }
    extra_args.push(genome_fa.to_str().unwrap_or(""));

    let status = Command::new(&rm)
        .args(&cmd_args)
        .args(&extra_args)
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        eprintln!("  ⚠  RepeatMasker returned non-zero — partial results may be available");
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// RepeatMasker .out parser (same as repeatmasker.rs)
// ─────────────────────────────────────────────────────────────────────────────

fn parse_rm_out(path: &Path, min_length: usize) -> Result<Vec<MaskedRegion>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut regions = Vec::new();

    for (line_num, line) in std::io::BufReader::new(file).lines().enumerate() {
        let line = line.map_err(MycoNoteError::Io)?;
        // Skip 3 header lines
        if line_num < 3 { continue; }
        let line = line.trim();
        if line.is_empty() { continue; }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 7 { continue; }

        let seqid = fields[4].to_string();
        let start: usize = match fields[5].parse::<usize>() {
            Ok(n) => n.saturating_sub(1),  // 1-based → 0-based
            Err(_) => continue,
        };
        let end: usize = match fields[6].parse::<usize>() {
            Ok(n) => n,
            Err(_) => continue,
        };

        if end.saturating_sub(start) < min_length { continue; }

        regions.push(MaskedRegion { seqid, start, end });
    }

    Ok(regions)
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn write_fasta_temp(path: &Path, sequences: &[(String, String)]) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    for (id, seq) in sequences {
        writeln!(f, ">{}", id).map_err(MycoNoteError::Io)?;
        for chunk in seq.as_bytes().chunks(60) {
            writeln!(f, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                .map_err(MycoNoteError::Io)?;
        }
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Availability check
// ─────────────────────────────────────────────────────────────────────────────

pub fn is_available() -> bool {
    which::which("RepeatModeler").is_ok() && which::which("BuildDatabase").is_ok()
}
