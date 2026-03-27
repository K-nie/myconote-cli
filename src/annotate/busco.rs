/// BUSCO completeness assessment wrapper
///
/// Assesses gene set completeness using BUSCO against kingdom-appropriate
/// OrthoDB lineage datasets.
///
/// Install: conda install -c conda-forge -c bioconda busco=5
/// Lineage datasets are auto-downloaded by BUSCO on first run.

use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Summary record
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct BuscoSummary {
    pub lineage:     String,
    pub single_copy: usize,
    pub duplicated:  usize,
    pub fragmented:  usize,
    pub missing:     usize,
    pub total:       usize,
}

impl BuscoSummary {
    pub fn complete(&self) -> usize {
        self.single_copy + self.duplicated
    }

    pub fn percent_complete(&self) -> f64 {
        if self.total == 0 { return 0.0; }
        self.complete() as f64 / self.total as f64 * 100.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run BUSCO on a protein FASTA and return the summary.
pub fn run(
    proteins_fa: &Path,
    lineage:     &str,
    out_dir:     &Path,
    threads:     usize,
) -> Result<BuscoSummary> {
    let busco = which::which("busco").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "busco not found in PATH.\n\
             Install with: conda install -c conda-forge -c bioconda busco=5".to_string()
        )
    })?;

    println!("  Running BUSCO (lineage: {})…", lineage);

    // BUSCO needs to write to a directory; use a subdirectory
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let run_name = "myconote_busco";

    let status = Command::new(&busco)
        .current_dir(out_dir)
        .args([
            "-m",  "proteins",
            "-i",  proteins_fa.to_str().unwrap_or(""),
            "-o",  run_name,
            "-l",  lineage,
            "--cpu", &threads.to_string(),
            "--offline",   // don't try to download datasets at runtime
            "--force",     // overwrite existing run
        ])
        .status()
        .map_err(MycoNoteError::Io)?;

    // BUSCO may exit non-zero even when partially succeeding — try to parse
    if !status.success() {
        eprintln!("  ⚠  BUSCO exited with non-zero status — attempting to parse results anyway…");
    }

    // Find the short summary file
    let summary_path = find_summary(out_dir, run_name);
    match summary_path {
        Some(p) => parse_busco_summary(&p, lineage),
        None => {
            // Re-try without --offline (dataset may need downloading)
            run_with_download(proteins_fa, lineage, out_dir, threads, run_name, &busco)
        }
    }
}

fn run_with_download(
    proteins_fa: &Path,
    lineage:     &str,
    out_dir:     &Path,
    threads:     usize,
    run_name:    &str,
    busco:       &Path,
) -> Result<BuscoSummary> {
    println!("  Downloading BUSCO dataset {} (first time only)…", lineage);

    let status = Command::new(busco)
        .current_dir(out_dir)
        .args([
            "-m",  "proteins",
            "-i",  proteins_fa.to_str().unwrap_or(""),
            "-o",  run_name,
            "-l",  lineage,
            "--cpu", &threads.to_string(),
            "--force",
        ])
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "BUSCO failed. Check that the lineage dataset can be downloaded.".to_string()
        ));
    }

    let summary_path = find_summary(out_dir, run_name)
        .ok_or_else(|| MycoNoteError::InvalidFormat(
            "BUSCO finished but no short_summary file found.".to_string()
        ))?;

    parse_busco_summary(&summary_path, lineage)
}

// ─────────────────────────────────────────────────────────────────────────────
// File helpers
// ─────────────────────────────────────────────────────────────────────────────

fn find_summary(out_dir: &Path, run_name: &str) -> Option<PathBuf> {
    // BUSCO 5 puts summary in: {out_dir}/{run_name}/short_summary.specific.*.txt
    let run_dir = out_dir.join(run_name);
    if let Ok(entries) = std::fs::read_dir(&run_dir) {
        for entry in entries.flatten() {
            let name = entry.file_name();
            let name = name.to_string_lossy();
            if name.starts_with("short_summary") && name.ends_with(".txt") {
                return Some(entry.path());
            }
        }
    }
    // BUSCO 4 path
    let v4 = out_dir.join(run_name).join(format!("short_summary.{}.txt", run_name));
    if v4.exists() { return Some(v4); }

    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Summary parser
// ─────────────────────────────────────────────────────────────────────────────

/// Parse BUSCO short_summary.txt file.
/// Handles BUSCO v4 and v5 formats.
fn parse_busco_summary(path: &Path, lineage: &str) -> Result<BuscoSummary> {
    let content = std::fs::read_to_string(path).map_err(MycoNoteError::Io)?;

    let mut s = 0usize;
    let mut d = 0usize;
    let mut f = 0usize;
    let mut m = 0usize;
    let mut total = 0usize;

    for line in content.lines() {
        let line = line.trim();
        // e.g. "C:98.0%[S:97.2%,D:0.8%],F:0.5%,M:1.5%,n:758"
        if line.starts_with("C:") {
            // Parse the compact format
            parse_compact_line(line, &mut s, &mut d, &mut f, &mut m, &mut total);
        }
        // Verbose format: "  XXX Complete BUSCOs (C)"
        if let Some(n) = extract_busco_num(line, "Complete BUSCOs") {
            if s + d == 0 { s = n; }
        }
        if let Some(n) = extract_busco_num(line, "Complete and single-copy BUSCOs") {
            s = n;
        }
        if let Some(n) = extract_busco_num(line, "Complete and duplicated BUSCOs") {
            d = n;
        }
        if let Some(n) = extract_busco_num(line, "Fragmented BUSCOs") {
            f = n;
        }
        if let Some(n) = extract_busco_num(line, "Missing BUSCOs") {
            m = n;
        }
        if let Some(n) = extract_busco_num(line, "Total BUSCO groups searched") {
            total = n;
        }
    }

    if total == 0 { total = s + d + f + m; }

    Ok(BuscoSummary {
        lineage: lineage.to_string(),
        single_copy: s,
        duplicated:  d,
        fragmented:  f,
        missing:     m,
        total,
    })
}

/// Parse "C:98.0%[S:97.2%,D:0.8%],F:0.5%,M:1.5%,n:758" style line.
fn parse_compact_line(line: &str, s: &mut usize, d: &mut usize,
                       f: &mut usize, m: &mut usize, n: &mut usize) {
    // Extract n: (total count)
    if let Some(pos) = line.find("n:") {
        let rest = &line[pos + 2..];
        let end = rest.find(|c: char| !c.is_ascii_digit()).unwrap_or(rest.len());
        *n = rest[..end].parse().unwrap_or(0);
    }
    // We don't parse percentages here — counts come from verbose lines
    let _ = (s, d, f, m);
}

fn extract_busco_num(line: &str, label: &str) -> Option<usize> {
    if line.contains(label) {
        let num_str: String = line.chars().take_while(|c| c.is_ascii_digit()).collect();
        num_str.parse().ok()
    } else {
        None
    }
}
