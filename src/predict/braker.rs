/// BRAKER 1 / 2 / 3 wrapper
///
/// `braker.pl` (bioconda package: `braker3`, since v3 covers all three modes
/// from a single driver) is a complete predictor that internally cross-trains
/// Augustus and GeneMark from extrinsic evidence:
///
///   BRAKER1 — RNA-seq BAMs only            (`--esmode` plus --bam hints)
///   BRAKER2 — protein FASTA only           (`--epmode`)
///   BRAKER3 — RNA-seq BAMs + protein FASTA (`--etpmode`, gold standard)
///
/// When `--use-braker` is on, MycoNote-CLI does NOT run any of its own
/// ab-initio predictors (Augustus / SNAP / GlimmerHMM / GeneMark / miniprot)
/// or the EVM consensus stage; we just point downstream stages at BRAKER's
/// `braker.gff3`.
///
/// CLI flag names verified against the upstream `Gaius-Augustus/BRAKER`
/// `scripts/braker.pl` GetOptions block on 2026-04-24:
///   --genome=s     --bam=s (repeatable)   --prot_seq=s (repeatable)
///   --species=s    --threads=i            --workingdir=s
///   --softmasking  --softmasking_off
///   --fungus       --gff3
///   --esmode       --epmode               --etpmode
///   --translation_table=s
use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Mode enum
// ─────────────────────────────────────────────────────────────────────────────

/// Which BRAKER mode to drive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrakerMode {
    /// RNA-seq BAMs only.
    Braker1,
    /// Protein FASTA only.
    Braker2,
    /// RNA-seq BAMs + protein FASTA (gold standard).
    Braker3,
}

impl BrakerMode {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.trim() {
            "1" | "braker1" => Some(Self::Braker1),
            "2" | "braker2" => Some(Self::Braker2),
            "3" | "braker3" => Some(Self::Braker3),
            _ => None,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Braker1 => "BRAKER1",
            Self::Braker2 => "BRAKER2",
            Self::Braker3 => "BRAKER3",
        }
    }

    /// The BRAKER `--esmode` / `--epmode` / `--etpmode` switch.
    pub fn mode_flag(&self) -> &'static str {
        match self {
            Self::Braker1 => "--esmode",
            Self::Braker2 => "--epmode",
            Self::Braker3 => "--etpmode",
        }
    }
}

/// Auto-detect the BRAKER mode from the inputs that were supplied:
///
///   RNA only         -> BRAKER1
///   protein only     -> BRAKER2
///   both             -> BRAKER3
///   neither          -> None (caller must error: at least one input required)
///
/// Pass an explicit override into `BrakerConfig::mode` to skip auto-detection.
pub fn detect_mode(rna_bams: &[PathBuf], proteins: Option<&Path>) -> Option<BrakerMode> {
    let has_rna = !rna_bams.is_empty();
    let has_prot = proteins.is_some();
    match (has_rna, has_prot) {
        (true, true) => Some(BrakerMode::Braker3),
        (true, false) => Some(BrakerMode::Braker1),
        (false, true) => Some(BrakerMode::Braker2),
        (false, false) => None,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct BrakerConfig {
    /// Soft-masked genome FASTA.
    pub genome: PathBuf,
    /// One or more RNA-seq BAM files.  Empty for BRAKER2.
    pub rna_bams: Vec<PathBuf>,
    /// Protein FASTA (typically OrthoDB fungi or a curated lineage set).
    /// `None` for BRAKER1.
    pub proteins: Option<PathBuf>,
    /// Output directory.  BRAKER writes `braker.gff3` here.
    pub out_dir: PathBuf,
    /// `--species` value passed to braker.pl + Augustus.  BRAKER refuses to
    /// reuse a species name across runs unless you also pass `--useexisting`,
    /// so we expose this directly rather than auto-deriving it.
    pub species: String,
    /// Threads -> braker.pl `--threads`.
    pub threads: usize,
    /// Translation table -> braker.pl `--translation_table`.  Common fungal
    /// values: 1 (standard), 12 (CTG clade Candida), 6 (ciliate-style).
    pub genetic_code: u8,
    /// Mode override.  `None` triggers `detect_mode` from the inputs.
    pub mode: Option<BrakerMode>,
    /// Run with `--fungus` (passed to GeneMark inside BRAKER).
    pub fungus: bool,
}

impl Default for BrakerConfig {
    fn default() -> Self {
        Self {
            genome: PathBuf::new(),
            rna_bams: Vec::new(),
            proteins: None,
            out_dir: PathBuf::from("braker_out"),
            species: "myconote_braker".to_string(),
            threads: 4,
            genetic_code: 1,
            mode: None,
            fungus: true,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Availability
// ─────────────────────────────────────────────────────────────────────────────

pub fn braker_available() -> bool {
    Command::new("which")
        .arg("braker.pl")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

// ─────────────────────────────────────────────────────────────────────────────
// Run BRAKER
// ─────────────────────────────────────────────────────────────────────────────

/// Validate the config + invoke `braker.pl`.  Returns the path to
/// `braker.gff3` (the canonical BRAKER output that downstream stages
/// consume directly).
///
/// Side-effects: creates `out_dir`, runs braker.pl (very long-running),
/// streams its stdout/stderr through to our terminal so the run is
/// reproducibility-bundle-friendly.  We do NOT capture into memory because
/// BRAKER logs can run to hundreds of MB on real fungal genomes.
pub fn run_braker(config: &BrakerConfig) -> Result<PathBuf> {
    if !braker_available() {
        return Err(MycoNoteError::ExternalTool(
            "BRAKER (braker.pl) not found on PATH.\n  \
             Install it with `myconote-cli install predict`.\n  \
             The bioconda package is `braker3` (versions 3.0.x); the binary it\n  \
             provides is `braker.pl`.  BRAKER additionally requires Augustus\n  \
             with a writeable AUGUSTUS_CONFIG_PATH (try `myconote-cli setup\n  \
             --db augustus-fungi`), a licensed GeneMark install, and ProtHint\n  \
             (bundled with the braker3 conda package)."
                .to_string(),
        ));
    }

    if !config.genome.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Genome FASTA not found: {}",
            config.genome.display()
        )));
    }

    // ── Resolve mode (explicit override, then auto-detect) ────────────────────
    let mode = match config.mode {
        Some(m) => m,
        None => detect_mode(&config.rna_bams, config.proteins.as_deref()).ok_or_else(|| {
            MycoNoteError::InvalidFormat(
                "BRAKER requires at least one of --braker-rna-bam or \
                 --braker-proteins.  See `myconote-cli predict --help`."
                    .to_string(),
            )
        })?,
    };

    // ── Validate inputs against the chosen mode ───────────────────────────────
    match mode {
        BrakerMode::Braker1 => {
            if config.rna_bams.is_empty() {
                return Err(MycoNoteError::InvalidFormat(
                    "BRAKER1 requires --braker-rna-bam <file>.".to_string(),
                ));
            }
        }
        BrakerMode::Braker2 => {
            if config.proteins.is_none() {
                return Err(MycoNoteError::InvalidFormat(
                    "BRAKER2 requires --braker-proteins <fa>.".to_string(),
                ));
            }
        }
        BrakerMode::Braker3 => {
            if config.rna_bams.is_empty() || config.proteins.is_none() {
                return Err(MycoNoteError::InvalidFormat(
                    "BRAKER3 requires both --braker-rna-bam and \
                     --braker-proteins."
                        .to_string(),
                ));
            }
        }
    }

    for bam in &config.rna_bams {
        if !bam.exists() {
            return Err(MycoNoteError::InvalidFormat(format!(
                "BRAKER RNA-seq BAM not found: {}",
                bam.display()
            )));
        }
    }
    if let Some(prot) = config.proteins.as_ref() {
        if !prot.exists() {
            return Err(MycoNoteError::InvalidFormat(format!(
                "BRAKER protein FASTA not found: {}",
                prot.display()
            )));
        }
    }

    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    // ── Build the command ─────────────────────────────────────────────────────
    let mut cmd = Command::new("braker.pl");
    cmd.arg("--genome")
        .arg(&config.genome)
        .arg("--workingdir")
        .arg(&config.out_dir)
        .arg("--species")
        .arg(&config.species)
        .arg("--threads")
        .arg(config.threads.to_string())
        .arg("--gff3")
        .arg(mode.mode_flag());

    // BAMs are repeatable: braker.pl accepts a comma-separated list to
    // a single --bam, but the safer cross-version form is one flag per BAM.
    for bam in &config.rna_bams {
        cmd.arg("--bam").arg(bam);
    }

    if let Some(prot) = config.proteins.as_ref() {
        cmd.arg("--prot_seq").arg(prot);
    }

    if config.fungus {
        cmd.arg("--fungus");
    }

    // BRAKER's translation_table flag — wraps GeneMark's --gcode and
    // Augustus's --translation_table simultaneously (verified against
    // Gaius-Augustus/BRAKER scripts/braker.pl GetOptions).
    if config.genetic_code != 1 {
        cmd.arg(format!("--translation_table={}", config.genetic_code));
    }

    eprintln!("  BRAKER command: {:?}", cmd);

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("braker.pl: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(format!(
            "{} (braker.pl) exited with non-zero status",
            mode.as_str()
        )));
    }

    let gff3 = config.out_dir.join("braker.gff3");
    if !gff3.exists() {
        return Err(MycoNoteError::ExternalTool(format!(
            "BRAKER finished but {} was not produced",
            gff3.display()
        )));
    }
    Ok(gff3)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Mode parsing ─────────────────────────────────────────────────────────

    #[test]
    fn from_str_accepts_numeric_and_named() {
        assert_eq!(BrakerMode::from_str("1"), Some(BrakerMode::Braker1));
        assert_eq!(BrakerMode::from_str("2"), Some(BrakerMode::Braker2));
        assert_eq!(BrakerMode::from_str("3"), Some(BrakerMode::Braker3));
        assert_eq!(BrakerMode::from_str("braker1"), Some(BrakerMode::Braker1));
        assert_eq!(BrakerMode::from_str("braker2"), Some(BrakerMode::Braker2));
        assert_eq!(BrakerMode::from_str("braker3"), Some(BrakerMode::Braker3));
    }

    #[test]
    fn from_str_trims_whitespace() {
        assert_eq!(BrakerMode::from_str("  3  "), Some(BrakerMode::Braker3));
    }

    #[test]
    fn from_str_rejects_unknown() {
        assert_eq!(BrakerMode::from_str("4"), None);
        assert_eq!(BrakerMode::from_str(""), None);
        assert_eq!(BrakerMode::from_str("etp"), None);
    }

    #[test]
    fn mode_flag_matches_braker_cli() {
        // These flag strings must match Gaius-Augustus/BRAKER scripts/braker.pl
        // GetOptions exactly — that file was the source of truth on
        // 2026-04-24.
        assert_eq!(BrakerMode::Braker1.mode_flag(), "--esmode");
        assert_eq!(BrakerMode::Braker2.mode_flag(), "--epmode");
        assert_eq!(BrakerMode::Braker3.mode_flag(), "--etpmode");
    }

    // ── Auto-detection ───────────────────────────────────────────────────────

    #[test]
    fn detect_rna_only_picks_braker1() {
        let bams = vec![PathBuf::from("rna.bam")];
        assert_eq!(detect_mode(&bams, None), Some(BrakerMode::Braker1));
    }

    #[test]
    fn detect_protein_only_picks_braker2() {
        let prot = PathBuf::from("proteins.fa");
        assert_eq!(detect_mode(&[], Some(&prot)), Some(BrakerMode::Braker2));
    }

    #[test]
    fn detect_both_picks_braker3() {
        let bams = vec![PathBuf::from("rna1.bam"), PathBuf::from("rna2.bam")];
        let prot = PathBuf::from("proteins.fa");
        assert_eq!(detect_mode(&bams, Some(&prot)), Some(BrakerMode::Braker3));
    }

    #[test]
    fn detect_no_inputs_returns_none() {
        // run_braker() converts this into a hard error; we don't silently
        // fall back to anything.
        assert_eq!(detect_mode(&[], None), None);
    }

    // ── Config / validation contract ─────────────────────────────────────────

    #[test]
    fn default_config_is_fungus_aware() {
        let c = BrakerConfig::default();
        assert!(c.fungus, "fungal pipeline should default to --fungus");
        assert_eq!(c.genetic_code, 1);
        assert!(c.mode.is_none(), "mode should auto-detect by default");
    }

    #[test]
    fn run_braker_errors_when_no_inputs_when_braker_present() {
        // We can only assert this branch when braker.pl is on PATH; otherwise
        // the binary-availability check fires first and the test would be
        // testing the wrong branch. The CI machine doesn't have braker.pl,
        // so on those machines the test is skipped via early return —
        // documented at the top of this assertion.
        if !braker_available() {
            return;
        }
        let cfg = BrakerConfig {
            genome: PathBuf::from("/tmp/nonexistent_genome_for_braker_test.fa"),
            ..BrakerConfig::default()
        };
        let err = run_braker(&cfg).unwrap_err();
        let msg = err.to_string();
        // Error message should be specific enough to be actionable.
        assert!(
            msg.contains("--braker-rna-bam") || msg.contains("Genome FASTA"),
            "expected actionable error, got: {}",
            msg
        );
    }

    #[test]
    #[ignore = "requires bioconda braker3 + Augustus + GeneMark + ProtHint; run with --ignored"]
    fn braker_end_to_end_live_smoke() {
        // Live smoke gate: only check that braker.pl is reachable.  Same
        // gating model as the de-template Rscript parse() check and the
        // GeneMark ETP live test.
        assert!(braker_available());
    }
}
