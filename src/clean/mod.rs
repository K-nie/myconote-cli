//! Duplicate-contig cleanup for FASTA assemblies.
//!
//! `myconote-cli clean --mode contigs <genome.fa>` runs minimap2 self-alignment
//! and drops contigs that are largely redundant — i.e. covered above a
//! coverage threshold by a longer contig at high identity. Common in
//! draft assemblies where alternative haplotigs sneak in alongside the
//! primary contigs (purge_dups / purge_haplotigs territory).
//!
//! GFF3-mode cleanup lives in `main.rs::handle_clean`; this module is
//! reachable only via `--mode contigs`.

use crate::utils::error::{MycoNoteError, Result};
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct ContigCleanConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    /// TSV report listing dropped contigs. Defaults to
    /// `<input_stem>_dropped.tsv` next to the output FASTA.
    pub report: PathBuf,
    /// Minimum query-coverage of the shorter contig that must be matched
    /// by a longer one before we consider it redundant.
    pub coverage: f64,
    /// Minimum gap-compressed identity over the matched span.
    pub identity: f64,
    pub threads: usize,
    /// Path to minimap2; defaults to "minimap2" from PATH.
    pub minimap2: PathBuf,
}

impl Default for ContigCleanConfig {
    fn default() -> Self {
        Self {
            input: PathBuf::new(),
            output: PathBuf::new(),
            report: PathBuf::new(),
            coverage: 0.95,
            identity: 0.95,
            threads: 4,
            minimap2: PathBuf::from("minimap2"),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_contig_clean(cfg: &ContigCleanConfig) -> Result<ContigCleanReport> {
    if !cfg.input.is_file() {
        return Err(MycoNoteError::ExternalTool(format!(
            "clean --mode contigs: input FASTA not found at {}",
            cfg.input.display()
        )));
    }
    if !(0.0..=1.0).contains(&cfg.coverage) {
        return Err(MycoNoteError::ExternalTool(format!(
            "clean --mode contigs: --coverage must be in [0,1], got {}",
            cfg.coverage
        )));
    }
    if !(0.0..=1.0).contains(&cfg.identity) {
        return Err(MycoNoteError::ExternalTool(format!(
            "clean --mode contigs: --identity must be in [0,1], got {}",
            cfg.identity
        )));
    }

    // Pre-flight: minimap2 must be on PATH (or wherever cfg.minimap2 points).
    let minimap2_ok = Command::new(&cfg.minimap2)
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !minimap2_ok {
        return Err(MycoNoteError::ExternalTool(format!(
            "minimap2 not found (looked for `{}`). \
             Install: conda install -c bioconda minimap2",
            cfg.minimap2.display()
        )));
    }

    println!(
        "Cleaning contigs in {} (coverage>={}, identity>={})",
        cfg.input.display(),
        cfg.coverage,
        cfg.identity
    );

    // Pass 1: read FASTA to learn each contig's length.
    let lengths = read_contig_lengths(&cfg.input)?;
    if lengths.is_empty() {
        eprintln!("  ⚠ no contigs in {}", cfg.input.display());
    }

    // Pass 2: minimap2 self-alignment → PAF on stdout.
    let paf_path = cfg.output.with_extension("paf");
    run_minimap2_self(cfg, &paf_path)?;

    // Pass 3: parse PAF, decide which contigs to drop.
    let drops = decide_drops(&paf_path, &lengths, cfg.coverage, cfg.identity)?;

    // Pass 4: write cleaned FASTA + drop report.
    let kept = write_cleaned_fasta(&cfg.input, &cfg.output, &drops)?;
    write_drop_report(&cfg.report, &drops, &lengths)?;

    // Tidy up the intermediate PAF — keep on failure for debugging,
    // remove on success.
    let _ = fs::remove_file(&paf_path);

    println!("\n  Contigs kept:    {}", kept);
    println!("  Contigs dropped: {}", drops.len());
    println!("  Cleaned FASTA:   {}", cfg.output.display());
    println!("  Drop report:     {}", cfg.report.display());

    Ok(ContigCleanReport {
        contigs_input: lengths.len(),
        contigs_kept: kept,
        contigs_dropped: drops.len(),
    })
}

#[derive(Debug, Clone, Copy)]
pub struct ContigCleanReport {
    pub contigs_input: usize,
    pub contigs_kept: usize,
    pub contigs_dropped: usize,
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTA helpers — minimal hand-rolled streaming reader
// ─────────────────────────────────────────────────────────────────────────────

/// Read sequence lengths only. Streams the file; never holds the
/// genome in RAM.
pub fn read_contig_lengths(path: &Path) -> Result<BTreeMap<String, usize>> {
    let f = File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(f);
    let mut lengths: BTreeMap<String, usize> = BTreeMap::new();
    let mut current: Option<String> = None;
    let mut acc: usize = 0;

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        if let Some(stripped) = line.strip_prefix('>') {
            if let Some(name) = current.take() {
                lengths.insert(name, acc);
            }
            // Header is up to the first whitespace — matches samtools/seqkit.
            let id = stripped.split_whitespace().next().unwrap_or("").to_string();
            current = Some(id);
            acc = 0;
        } else {
            acc += line.trim().len();
        }
    }
    if let Some(name) = current.take() {
        lengths.insert(name, acc);
    }
    Ok(lengths)
}

/// Write a new FASTA dropping every record whose ID is in `drop_set`.
/// Returns the number of records kept.
fn write_cleaned_fasta(
    input: &Path,
    output: &Path,
    drops: &BTreeMap<String, DropEntry>,
) -> Result<usize> {
    let in_f = File::open(input).map_err(MycoNoteError::Io)?;
    let out_f = File::create(output).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(in_f);
    let mut writer = std::io::BufWriter::new(out_f);

    let drop_set: HashSet<&String> = drops.keys().collect();
    let mut kept_records = 0usize;
    let mut keep_current = true;

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        if let Some(stripped) = line.strip_prefix('>') {
            let id = stripped.split_whitespace().next().unwrap_or("").to_string();
            keep_current = !drop_set.contains(&id);
            if keep_current {
                kept_records += 1;
                writeln!(writer, "{}", line).map_err(MycoNoteError::Io)?;
            }
        } else if keep_current {
            writeln!(writer, "{}", line).map_err(MycoNoteError::Io)?;
        }
    }
    Ok(kept_records)
}

// ─────────────────────────────────────────────────────────────────────────────
// minimap2 invocation
// ─────────────────────────────────────────────────────────────────────────────

fn run_minimap2_self(cfg: &ContigCleanConfig, paf_out: &Path) -> Result<()> {
    // -X drops self-mapping, asm5 is appropriate for highly similar
    // contigs from the same assembly. PAF goes to stdout; we redirect
    // to a file so the parser can stream it back.
    let paf_file = File::create(paf_out).map_err(MycoNoteError::Io)?;
    let status = Command::new(&cfg.minimap2)
        .args([
            "-x",
            "asm5",
            "-X",
            "-t",
            &cfg.threads.to_string(),
            cfg.input.to_str().unwrap_or(""),
            cfg.input.to_str().unwrap_or(""),
        ])
        .stdout(paf_file)
        .stderr(std::process::Stdio::inherit())
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("minimap2 failed to launch: {}", e)))?;
    if !status.success() {
        return Err(MycoNoteError::ExternalTool(format!(
            "minimap2 self-alignment failed (exit {})",
            status.code().unwrap_or(-1)
        )));
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// PAF parsing + drop decision
// ─────────────────────────────────────────────────────────────────────────────

/// One row recording why a contig is being dropped.
#[derive(Debug, Clone)]
pub struct DropEntry {
    pub subsumed_by: String,
    pub coverage: f64,
    pub identity: f64,
}

/// Per-pair best alignment summary used to make a drop call. Public for
/// unit tests; not exposed in the CLI.
#[derive(Debug, Clone, Copy, Default)]
struct PairStats {
    /// Sum of matched bases on the shorter contig.
    matched: usize,
    /// Sum of aligned bases on the shorter contig (incl. mismatches & gaps).
    aligned: usize,
}

/// Walk the PAF and decide which contigs to drop.
/// Returns a map from dropped-contig-id → DropEntry with the contig
/// that subsumed it and the (coverage, identity) of the call.
pub fn decide_drops(
    paf_path: &Path,
    lengths: &BTreeMap<String, usize>,
    min_coverage: f64,
    min_identity: f64,
) -> Result<BTreeMap<String, DropEntry>> {
    let f = File::open(paf_path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(f);

    // Aggregate per (shorter, longer) ordered pair: alignments may
    // appear in many slices, and we want total coverage of the shorter
    // contig by the longer.
    let mut pair_stats: BTreeMap<(String, String), PairStats> = BTreeMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let row = match parse_paf_line(&line) {
            Some(r) => r,
            None => continue,
        };

        // Decide which is shorter (the candidate to drop) and which
        // could subsume it. Skip self-pairs (minimap2 -X already does
        // this, but we double-check).
        if row.qname == row.tname {
            continue;
        }

        let qlen = *lengths.get(&row.qname).unwrap_or(&row.qlen);
        let tlen = *lengths.get(&row.tname).unwrap_or(&row.tlen);

        // Tie-break by lexicographic order so the call is
        // deterministic when two contigs are exactly the same length.
        // On ties the *lex-larger* name is treated as shorter (the
        // candidate to drop) so the lex-smaller name always survives.
        let (shorter, longer, matched, aligned) =
            if qlen < tlen || (qlen == tlen && row.qname > row.tname) {
                (
                    row.qname.clone(),
                    row.tname.clone(),
                    row.matches,
                    row.aln_block_len,
                )
            } else {
                // The PAF row is q-on-t; matched bases are the same in both
                // orientations, but the aligned span on the *target* side
                // is (tend-tstart). Use that for the "shorter is t" case.
                (
                    row.tname.clone(),
                    row.qname.clone(),
                    row.matches,
                    (row.tend - row.tstart),
                )
            };

        let entry = pair_stats.entry((shorter, longer)).or_default();
        entry.matched += matched;
        entry.aligned += aligned;
    }

    // For each candidate-shorter, pick the best subsuming partner that
    // clears both thresholds. Lengths are looked up from the FASTA pass
    // so coverage is computed against the true contig length.
    let mut best: BTreeMap<String, DropEntry> = BTreeMap::new();
    for ((shorter, longer), stats) in &pair_stats {
        let s_len = match lengths.get(shorter) {
            Some(&n) if n > 0 => n,
            _ => continue,
        };
        let coverage = stats.aligned as f64 / s_len as f64;
        let identity = if stats.aligned > 0 {
            stats.matched as f64 / stats.aligned as f64
        } else {
            0.0
        };

        if coverage < min_coverage || identity < min_identity {
            continue;
        }

        // Keep the best (highest coverage × identity) subsumer.
        let score = coverage * identity;
        let prev = best
            .get(shorter)
            .map(|d| d.coverage * d.identity)
            .unwrap_or(-1.0);
        if score > prev {
            best.insert(
                shorter.clone(),
                DropEntry {
                    subsumed_by: longer.clone(),
                    coverage,
                    identity,
                },
            );
        }
    }

    // Resolve mutual subsumption: if A is going to be dropped because
    // B subsumes it, but B is also going to be dropped because of A,
    // keep the longer one (or, on a tie, the lexicographically smaller
    // one — same rule as above).
    let drop_keys: Vec<String> = best.keys().cloned().collect();
    let mut resolved: BTreeMap<String, DropEntry> = BTreeMap::new();
    for k in &drop_keys {
        let by = &best[k].subsumed_by;
        if best.contains_key(by) {
            // Both want to drop each other. Keep the longer.
            let kl = lengths.get(k).copied().unwrap_or(0);
            let bl = lengths.get(by).copied().unwrap_or(0);
            // Keep the longer; on ties keep the lex-smaller name.
            if kl > bl || (kl == bl && k < by) {
                // k is the survivor; skip dropping it.
                continue;
            }
        }
        resolved.insert(k.clone(), best[k].clone());
    }

    Ok(resolved)
}

#[derive(Debug)]
struct PafRow {
    qname: String,
    qlen: usize,
    #[allow(dead_code)]
    qstart: usize,
    #[allow(dead_code)]
    qend: usize,
    tname: String,
    tlen: usize,
    tstart: usize,
    tend: usize,
    matches: usize,
    aln_block_len: usize,
}

fn parse_paf_line(line: &str) -> Option<PafRow> {
    if line.trim().is_empty() {
        return None;
    }
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 12 {
        return None;
    }
    Some(PafRow {
        qname: f[0].to_string(),
        qlen: f[1].parse().ok()?,
        qstart: f[2].parse().ok()?,
        qend: f[3].parse().ok()?,
        // strand at f[4]
        tname: f[5].to_string(),
        tlen: f[6].parse().ok()?,
        tstart: f[7].parse().ok()?,
        tend: f[8].parse().ok()?,
        matches: f[9].parse().ok()?,
        aln_block_len: f[10].parse().ok()?,
        // mapq at f[11], optional tags after
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Drop-report writer
// ─────────────────────────────────────────────────────────────────────────────

fn write_drop_report(
    report_path: &Path,
    drops: &BTreeMap<String, DropEntry>,
    lengths: &BTreeMap<String, usize>,
) -> Result<()> {
    let f = File::create(report_path).map_err(MycoNoteError::Io)?;
    let mut w = std::io::BufWriter::new(f);
    writeln!(
        w,
        "dropped_contig\tdropped_length\tsubsumed_by\tsubsumer_length\tcoverage\tidentity"
    )
    .map_err(MycoNoteError::Io)?;
    for (id, entry) in drops {
        let dl = lengths.get(id).copied().unwrap_or(0);
        let sl = lengths.get(&entry.subsumed_by).copied().unwrap_or(0);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
            id, dl, entry.subsumed_by, sl, entry.coverage, entry.identity
        )
        .map_err(MycoNoteError::Io)?;
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;

    fn write_fasta(dir: &Path, name: &str, records: &[(&str, &str)]) -> PathBuf {
        let p = dir.join(name);
        let mut f = File::create(&p).unwrap();
        for (id, seq) in records {
            writeln!(f, ">{}", id).unwrap();
            // Wrap at 60 like everyone else.
            for chunk in seq.as_bytes().chunks(60) {
                f.write_all(chunk).unwrap();
                writeln!(f).unwrap();
            }
        }
        p
    }

    #[test]
    fn read_lengths_handles_multi_line_seqs() {
        let dir = tempfile::tempdir().unwrap();
        let p = write_fasta(
            &dir.path(),
            "in.fa",
            &[("ctg1", &"A".repeat(150)), ("ctg2", "ACGT")],
        );
        let lens = read_contig_lengths(&p).unwrap();
        assert_eq!(lens.get("ctg1"), Some(&150));
        assert_eq!(lens.get("ctg2"), Some(&4));
    }

    #[test]
    fn read_lengths_strips_after_first_whitespace() {
        // Same convention samtools and seqkit use.
        let dir = tempfile::tempdir().unwrap();
        let p = write_fasta(
            &dir.path(),
            "in.fa",
            &[("ctg1 description here", "ACGTACGT")],
        );
        let lens = read_contig_lengths(&p).unwrap();
        assert_eq!(lens.get("ctg1"), Some(&8));
    }

    #[test]
    fn parse_paf_line_basic() {
        let line = "ctg2\t100\t0\t100\t+\tctg1\t1000\t0\t100\t98\t100\t60";
        let row = parse_paf_line(line).unwrap();
        assert_eq!(row.qname, "ctg2");
        assert_eq!(row.qlen, 100);
        assert_eq!(row.tname, "ctg1");
        assert_eq!(row.tlen, 1000);
        assert_eq!(row.matches, 98);
        assert_eq!(row.aln_block_len, 100);
    }

    #[test]
    fn parse_paf_line_rejects_short_rows() {
        assert!(parse_paf_line("not enough fields").is_none());
        assert!(parse_paf_line("").is_none());
    }

    #[test]
    fn decide_drops_picks_shorter_when_subsumed() {
        let dir = tempfile::tempdir().unwrap();
        let paf = dir.path().join("self.paf");
        // ctg2 (100 bp) is fully covered (100/100) at 98% identity by
        // ctg1 (1000 bp). Should be dropped.
        let mut f = File::create(&paf).unwrap();
        writeln!(f, "ctg2\t100\t0\t100\t+\tctg1\t1000\t0\t100\t98\t100\t60").unwrap();
        // ctg1 also has a tiny match against ctg2 in the reverse pair —
        // we should not drop ctg1 from that.
        writeln!(f, "ctg1\t1000\t0\t100\t+\tctg2\t100\t0\t100\t98\t100\t60").unwrap();
        drop(f);

        let mut lengths = BTreeMap::new();
        lengths.insert("ctg1".to_string(), 1000);
        lengths.insert("ctg2".to_string(), 100);

        let drops = decide_drops(&paf, &lengths, 0.95, 0.95).unwrap();
        assert!(drops.contains_key("ctg2"));
        assert!(!drops.contains_key("ctg1"));
        assert_eq!(drops["ctg2"].subsumed_by, "ctg1");
        assert!(drops["ctg2"].coverage >= 0.95);
        assert!(drops["ctg2"].identity >= 0.95);
    }

    #[test]
    fn decide_drops_skips_partial_coverage() {
        // ctg2 only has 50/100 bp covered — under 0.95 threshold.
        let dir = tempfile::tempdir().unwrap();
        let paf = dir.path().join("self.paf");
        let mut f = File::create(&paf).unwrap();
        writeln!(f, "ctg2\t100\t0\t50\t+\tctg1\t1000\t0\t50\t49\t50\t60").unwrap();
        drop(f);

        let mut lengths = BTreeMap::new();
        lengths.insert("ctg1".to_string(), 1000);
        lengths.insert("ctg2".to_string(), 100);

        let drops = decide_drops(&paf, &lengths, 0.95, 0.95).unwrap();
        assert!(drops.is_empty());
    }

    #[test]
    fn decide_drops_skips_low_identity() {
        // ctg2 is fully covered but only 70% identity — below 0.95.
        let dir = tempfile::tempdir().unwrap();
        let paf = dir.path().join("self.paf");
        let mut f = File::create(&paf).unwrap();
        writeln!(f, "ctg2\t100\t0\t100\t+\tctg1\t1000\t0\t100\t70\t100\t60").unwrap();
        drop(f);

        let mut lengths = BTreeMap::new();
        lengths.insert("ctg1".to_string(), 1000);
        lengths.insert("ctg2".to_string(), 100);

        let drops = decide_drops(&paf, &lengths, 0.95, 0.95).unwrap();
        assert!(drops.is_empty());
    }

    #[test]
    fn write_cleaned_fasta_drops_named_records() {
        let dir = tempfile::tempdir().unwrap();
        let in_p = write_fasta(
            &dir.path(),
            "in.fa",
            &[("ctg1", "ACGT"), ("ctg2", "TTTT"), ("ctg3", "GGGG")],
        );
        let out_p = dir.path().join("out.fa");

        let mut drops: BTreeMap<String, DropEntry> = BTreeMap::new();
        drops.insert(
            "ctg2".to_string(),
            DropEntry {
                subsumed_by: "ctg1".to_string(),
                coverage: 1.0,
                identity: 1.0,
            },
        );

        let kept = write_cleaned_fasta(&in_p, &out_p, &drops).unwrap();
        assert_eq!(kept, 2);
        let body = fs::read_to_string(&out_p).unwrap();
        assert!(body.contains(">ctg1"));
        assert!(!body.contains(">ctg2"));
        assert!(body.contains(">ctg3"));
    }

    #[test]
    fn write_drop_report_has_header_and_rows() {
        let dir = tempfile::tempdir().unwrap();
        let report = dir.path().join("dropped.tsv");
        let mut drops: BTreeMap<String, DropEntry> = BTreeMap::new();
        drops.insert(
            "ctg2".to_string(),
            DropEntry {
                subsumed_by: "ctg1".to_string(),
                coverage: 0.99,
                identity: 0.98,
            },
        );
        let mut lengths = BTreeMap::new();
        lengths.insert("ctg1".to_string(), 1000);
        lengths.insert("ctg2".to_string(), 100);

        write_drop_report(&report, &drops, &lengths).unwrap();
        let body = fs::read_to_string(&report).unwrap();
        assert!(body.starts_with("dropped_contig\t"));
        let lines: Vec<&str> = body.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[1].starts_with("ctg2\t100\tctg1\t1000\t"));
    }

    #[test]
    fn decide_drops_breaks_mutual_subsumption_by_length() {
        // Pathological case: ctg1 == ctg2 in length and content, both
        // fully cover each other. We must keep exactly one. The rule
        // is "longer wins, tiebreak lex".
        let dir = tempfile::tempdir().unwrap();
        let paf = dir.path().join("self.paf");
        let mut f = File::create(&paf).unwrap();
        writeln!(f, "ctg1\t100\t0\t100\t+\tctg2\t100\t0\t100\t100\t100\t60").unwrap();
        writeln!(f, "ctg2\t100\t0\t100\t+\tctg1\t100\t0\t100\t100\t100\t60").unwrap();
        drop(f);

        let mut lengths = BTreeMap::new();
        lengths.insert("ctg1".to_string(), 100);
        lengths.insert("ctg2".to_string(), 100);

        let drops = decide_drops(&paf, &lengths, 0.95, 0.95).unwrap();
        // Exactly one of them is dropped; lex-smaller "ctg1" survives.
        assert_eq!(drops.len(), 1);
        assert!(drops.contains_key("ctg2"));
    }
}
