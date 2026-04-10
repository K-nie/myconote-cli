/// minimap2 self-alignment repeat finder
///
/// Aligns the genome against itself using minimap2 in asm-to-asm mode.
/// Regions that align to multiple locations are repeats.  No external
/// database required — works purely from the input sequence.
///
/// Also detects tandem repeats natively in Rust using a sliding-window
/// period-detection algorithm (covers satellite DNA and microsatellites
/// that minimap2 misses due to its seed length).
use super::{MaskConfig, MaskedRegion};
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// minimap2 self-alignment
// ─────────────────────────────────────────────────────────────────────────────

pub fn run(sequences: &[(String, String)], config: &MaskConfig) -> Result<Vec<MaskedRegion>> {
    let mut all_regions: Vec<MaskedRegion> = Vec::new();

    // Tandem repeat detection (pure Rust — always runs)
    println!("  Detecting tandem repeats…");
    let tr_regions = find_tandem_repeats(sequences, config.min_length);
    println!("  {} tandem repeat regions found", tr_regions.len());
    all_regions.extend(tr_regions);

    // minimap2 self-alignment (interspersed repeats)
    match which::which("minimap2") {
        Ok(mm2) => {
            println!("  Running minimap2 self-alignment…");
            let mm2_regions = run_minimap2_self(sequences, config, &mm2)?;
            println!("  {} interspersed repeat regions found", mm2_regions.len());
            all_regions.extend(mm2_regions);
        }
        Err(_) => {
            eprintln!("  ℹ  minimap2 not in PATH — skipping interspersed repeat detection.");
            eprintln!("     Install with: conda install -c bioconda minimap2");
            eprintln!("     Tandem repeats will still be masked.");
        }
    }

    Ok(all_regions)
}

fn run_minimap2_self(
    sequences: &[(String, String)],
    config: &MaskConfig,
    mm2: &std::path::Path,
) -> Result<Vec<MaskedRegion>> {
    let tmp_dir = tempfile::TempDir::new().map_err(MycoNoteError::Io)?;
    let fasta_path = tmp_dir.path().join("genome.fa");
    let paf_path = tmp_dir.path().join("self.paf");

    // Write temp FASTA
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

    // Run minimap2: genome vs. itself, asm5 preset
    let status = Command::new(mm2)
        .args([
            "-x",
            "asm5",
            "-t",
            &config.threads.to_string(),
            "-D", // skip identical self-hits (diagonal)
            fasta_path.to_str().unwrap_or(""),
            fasta_path.to_str().unwrap_or(""),
        ])
        .stdout(std::fs::File::create(&paf_path).map_err(MycoNoteError::Io)?)
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "minimap2 self-alignment failed.".to_string(),
        ));
    }

    parse_paf_to_regions(&paf_path, config.min_length)
}

/// Parse PAF file → masked regions (both query and target coordinates).
fn parse_paf_to_regions(path: &std::path::Path, min_len: usize) -> Result<Vec<MaskedRegion>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut regions = Vec::new();
    // Track which query positions are multiply mapped
    let mut hit_counts: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }

        let qname = f[0];
        let qstart: usize = f[2].parse().unwrap_or(0);
        let qend: usize = f[3].parse().unwrap_or(0);
        let tname = f[5];
        let tstart: usize = f[7].parse().unwrap_or(0);
        let tend: usize = f[8].parse().unwrap_or(0);

        // Skip the trivial self-alignment (same seq, same coords)
        if qname == tname && qstart == tstart && qend == tend {
            continue;
        }

        if qend.saturating_sub(qstart) >= min_len {
            hit_counts
                .entry(qname.to_string())
                .or_default()
                .push((qstart, qend));
        }
        if tend.saturating_sub(tstart) >= min_len {
            hit_counts
                .entry(tname.to_string())
                .or_default()
                .push((tstart, tend));
        }
    }

    for (seqid, spans) in hit_counts {
        for (start, end) in spans {
            regions.push(MaskedRegion {
                seqid: seqid.clone(),
                start,
                end,
            });
        }
    }

    Ok(regions)
}

// ─────────────────────────────────────────────────────────────────────────────
// Native tandem repeat finder
// ─────────────────────────────────────────────────────────────────────────────

/// Detect tandem repeats using a simple sliding-window autocorrelation.
///
/// For each window of `window` bases, checks whether the sequence has
/// a periodic structure with period in range [2, max_period].
/// If the fraction of bases that match the period pattern exceeds
/// `threshold`, the window is flagged as a tandem repeat.
fn find_tandem_repeats(sequences: &[(String, String)], min_length: usize) -> Vec<MaskedRegion> {
    let window = 64usize; // detection window
    let step = 16usize; // step between windows
    let max_period = 50usize; // max repeat unit length (catches satellites)
    let threshold = 0.75f64; // fraction of matching bases to call as repeat

    let mut regions = Vec::new();

    for (id, seq) in sequences {
        let bytes = seq.as_bytes();
        let len = bytes.len();
        if len < window {
            continue;
        }

        // Collect windows that look like tandem repeats
        let mut repeat_windows: Vec<bool> = vec![false; len];

        let mut w_start = 0;
        while w_start + window <= len {
            let w = &bytes[w_start..w_start + window];
            'period: for period in 2..=max_period {
                let mut matches = 0usize;
                let check_len = window - period;
                for i in 0..check_len {
                    if w[i] == w[i + period] {
                        matches += 1;
                    }
                }
                if matches as f64 / check_len as f64 >= threshold {
                    for i in w_start..w_start + window {
                        if i < len {
                            repeat_windows[i] = true;
                        }
                    }
                    break 'period;
                }
            }
            w_start += step;
        }

        // Merge adjacent repeat windows into regions
        let mut in_repeat = false;
        let mut region_start = 0usize;
        for i in 0..len {
            match (in_repeat, repeat_windows[i]) {
                (false, true) => {
                    in_repeat = true;
                    region_start = i;
                }
                (true, false) => {
                    let region_len = i - region_start;
                    if region_len >= min_length {
                        regions.push(MaskedRegion {
                            seqid: id.clone(),
                            start: region_start,
                            end: i,
                        });
                    }
                    in_repeat = false;
                }
                _ => {}
            }
        }
        // Close open region at end of sequence
        if in_repeat && len - region_start >= min_length {
            regions.push(MaskedRegion {
                seqid: id.clone(),
                start: region_start,
                end: len,
            });
        }
    }

    regions
}
