/// Repeat masking module
///
/// Soft-masks repetitive elements in a genome FASTA so that gene predictors
/// (Augustus, SNAP, etc.) ignore them.  Lowercase = masked in soft-mask mode.
///
/// Engine options:
///   1. SelfAlign     — minimap2 self-alignment + native TRF; no database needed
///   2. RepeatMasker  — RepeatMasker with RepBase species database
///   3. Both          — RepeatMasker first, SelfAlign fills gaps
///   4. RepeatModeler — build de novo repeat library first, then RepeatMasker
///   5. Full          — RepeatModeler (de novo) + RepeatMasker + SelfAlign
///
/// Output: soft-masked FASTA (`genome_masked.fa`) + stats report

pub mod repeatmasker;
pub mod self_align;
pub mod repeatmodeler;

use crate::utils::error::{MycoNoteError, Result};
use std::path::PathBuf;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
pub enum MaskEngine {
    SelfAlign,      // minimap2 self-alignment (no external DB needed)
    RepeatMasker,   // RepeatMasker with RepBase
    Both,           // RepeatMasker + SelfAlign merged
    RepeatModeler,  // De novo library (RepeatModeler) → RepeatMasker
    Full,           // RepeatModeler → RepeatMasker → SelfAlign (most thorough)
}

impl MaskEngine {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "repeatmasker" | "rm"           => MaskEngine::RepeatMasker,
            "repeatmodeler" | "denovo"      => MaskEngine::RepeatModeler,
            "full" | "all"                  => MaskEngine::Full,
            "both"                          => MaskEngine::Both,
            _                               => MaskEngine::SelfAlign,
        }
    }
}

#[derive(Debug, Clone)]
pub struct MaskConfig {
    /// Input genome FASTA (unmasked)
    pub input: PathBuf,
    /// Output soft-masked FASTA (default: <input>_masked.fa)
    pub output: PathBuf,
    /// Which masking engine to use
    pub engine: MaskEngine,
    /// Organism species name for RepeatMasker database lookup
    /// e.g. "fungi", "arabidopsis", "human"
    pub species: Option<String>,
    /// Custom repeat library FASTA (overrides species database)
    pub repeat_lib: Option<PathBuf>,
    /// Number of threads for parallel execution
    pub threads: usize,
    /// Minimum repeat length to mask (bp)
    pub min_length: usize,
    /// If true, hard-mask (N) instead of soft-mask (lowercase)
    pub hard_mask: bool,
}

impl Default for MaskConfig {
    fn default() -> Self {
        Self {
            input:      PathBuf::new(),
            output:     PathBuf::new(),
            engine:     MaskEngine::SelfAlign,
            species:    None,
            repeat_lib: None,
            threads:    4,
            min_length: 200,
            hard_mask:  false,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Masking stats
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
pub struct MaskStats {
    pub total_bases:    u64,
    pub masked_bases:   u64,
    pub repeat_regions: usize,
}

impl MaskStats {
    pub fn percent_masked(&self) -> f64 {
        if self.total_bases == 0 { 0.0 }
        else { self.masked_bases as f64 / self.total_bases as f64 * 100.0 }
    }

    pub fn print_summary(&self) {
        println!("\n  Total bases:    {}", self.total_bases);
        println!("  Masked bases:   {} ({:.1}%)",
            self.masked_bases, self.percent_masked());
        println!("  Repeat regions: {}", self.repeat_regions);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Masked region: a span [start, end) on a given sequence (0-based half-open)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct MaskedRegion {
    pub seqid: String,
    pub start: usize,   // 0-based
    pub end:   usize,   // exclusive
}

// ─────────────────────────────────────────────────────────────────────────────
// Soft/hard masking of a FASTA given a list of regions
// ─────────────────────────────────────────────────────────────────────────────

/// Apply a set of masked regions to sequences and write the masked FASTA.
/// Regions outside any sequence are silently ignored.
pub fn apply_masking(
    sequences: &[(String, String)],  // (id, seq)
    regions:   &[MaskedRegion],
    config:    &MaskConfig,
) -> (Vec<(String, String)>, MaskStats) {
    use std::collections::HashMap;

    // Group regions by seqid
    let mut by_seq: HashMap<&str, Vec<(usize, usize)>> = HashMap::new();
    for r in regions {
        by_seq.entry(&r.seqid).or_default().push((r.start, r.end));
    }
    // Sort and merge overlapping regions per sequence
    for spans in by_seq.values_mut() {
        spans.sort_by_key(|s| s.0);
        let mut merged: Vec<(usize, usize)> = Vec::new();
        for &(s, e) in spans.iter() {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 { last.1 = last.1.max(e); continue; }
            }
            merged.push((s, e));
        }
        *spans = merged;
    }

    let mut stats = MaskStats::default();
    let mut masked_seqs: Vec<(String, String)> = Vec::new();

    for (id, seq) in sequences {
        stats.total_bases += seq.len() as u64;
        let mut chars: Vec<char> = seq.chars().collect();

        if let Some(spans) = by_seq.get(id.as_str()) {
            for &(s, e) in spans {
                let e = e.min(chars.len());
                if s >= e { continue; }
                stats.repeat_regions += 1;
                for c in &mut chars[s..e] {
                    stats.masked_bases += 1;
                    *c = if config.hard_mask {
                        'N'
                    } else {
                        c.to_lowercase().next().unwrap_or('n')
                    };
                }
            }
        }

        masked_seqs.push((id.clone(), chars.iter().collect()));
    }

    (masked_seqs, stats)
}

/// Write sequences as FASTA, wrapping at 60 chars per line.
pub fn write_fasta<W: std::io::Write>(
    w: &mut W,
    sequences: &[(String, String)],
) -> std::io::Result<()> {
    for (id, seq) in sequences {
        writeln!(w, ">{}", id)?;
        for chunk in seq.as_bytes().chunks(60) {
            writeln!(w, "{}", std::str::from_utf8(chunk).unwrap_or(""))?;
        }
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_masking(config: &MaskConfig) -> Result<MaskStats> {
    use crate::parser::read_fasta;
    use crate::progress;

    println!("── Repeat masking ───────────────────────────────────────────");
    println!("  Engine : {:?}", config.engine);
    if let Some(ref sp) = config.species {
        println!("  Species: {}", sp);
    }

    // Load genome
    let pb = progress::spinner(format!("Loading {}", config.input.display()));
    let records = read_fasta(&config.input)?;
    let sequences: Vec<(String, String)> = records.iter()
        .map(|r| (r.id.clone(), r.sequence.clone()))
        .collect();
    let total_bases: u64 = sequences.iter().map(|(_, s)| s.len() as u64).sum();
    progress::finish_spinner(&pb, format!(
        "{} sequences ({:.1} Mbp)", sequences.len(), total_bases as f64 / 1e6
    ));

    // Gather masked regions
    let regions: Vec<MaskedRegion> = match &config.engine {
        MaskEngine::SelfAlign => {
            let pb2 = progress::spinner("Running self-alignment masking…");
            let r = self_align::run(&sequences, config)?;
            progress::finish_spinner(&pb2, format!("{} regions found", r.len()));
            r
        }
        MaskEngine::RepeatMasker => {
            let pb2 = progress::spinner("Running RepeatMasker…");
            let r = repeatmasker::run(&sequences, config)?;
            progress::finish_spinner(&pb2, format!("{} regions found", r.len()));
            r
        }
        MaskEngine::Both => {
            let pb2 = progress::spinner("Running RepeatMasker…");
            let mut r = repeatmasker::run(&sequences, config)
                .unwrap_or_else(|e| {
                    progress::warn_spinner(&pb2, format!("RepeatMasker failed: {}", e));
                    Vec::new()
                });
            progress::finish_spinner(&pb2, format!("{} RepeatMasker regions", r.len()));
            let pb3 = progress::spinner("Running self-alignment masking…");
            let self_r = self_align::run(&sequences, config)?;
            progress::finish_spinner(&pb3, format!("{} self-align regions", self_r.len()));
            r.extend(self_r);
            r
        }
        MaskEngine::RepeatModeler => {
            let pb2 = progress::spinner(
                "Building de novo repeat library with RepeatModeler…"
            );
            let r = repeatmodeler::run(&sequences, config)?;
            progress::finish_spinner(&pb2, format!("{} de novo repeat regions", r.len()));
            r
        }
        MaskEngine::Full => {
            // Step 1: RepeatModeler de novo library → RepeatMasker
            let pb2 = progress::spinner("Step 1/3: RepeatModeler (de novo library)…");
            let mut r = repeatmodeler::run(&sequences, config)
                .unwrap_or_else(|e| {
                    progress::warn_spinner(&pb2, format!("RepeatModeler failed: {}", e));
                    Vec::new()
                });
            progress::finish_spinner(&pb2, format!("{} de novo regions", r.len()));

            // Step 2: RepeatMasker with species database (catches known families)
            let pb3 = progress::spinner("Step 2/3: RepeatMasker (species database)…");
            let rm_r = repeatmasker::run(&sequences, config)
                .unwrap_or_else(|e| {
                    progress::warn_spinner(&pb3, format!("RepeatMasker failed: {}", e));
                    Vec::new()
                });
            progress::finish_spinner(&pb3, format!("{} RepeatMasker regions", rm_r.len()));
            r.extend(rm_r);

            // Step 3: self-alignment (catches tandem repeats and novel elements)
            let pb4 = progress::spinner("Step 3/3: Self-alignment masking…");
            let sa_r = self_align::run(&sequences, config)?;
            progress::finish_spinner(&pb4, format!("{} self-align regions", sa_r.len()));
            r.extend(sa_r);
            r
        }
    };

    // Apply masking
    let (masked_seqs, stats) = apply_masking(&sequences, &regions, config);

    // Write output
    let out_path = if config.output.as_os_str().is_empty() {
        let stem = config.input.file_stem().unwrap_or_default().to_string_lossy();
        config.input.with_file_name(format!("{}_masked.fa", stem))
    } else {
        config.output.clone()
    };

    let mut out = std::fs::File::create(&out_path).map_err(MycoNoteError::Io)?;
    write_fasta(&mut out, &masked_seqs).map_err(MycoNoteError::Io)?;

    println!("✓ Masked genome: {}", out_path.display());
    stats.print_summary();

    Ok(stats)
}
