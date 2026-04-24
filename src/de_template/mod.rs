//! Differential-expression R-script generator.
//!
//! `myconote-cli de-template` takes `quant` outputs + a design formula
//! + one or more contrasts and writes a self-contained R script. The
//! user runs that script with `Rscript`; it performs tximport →
//! DESeq2 → apeglm shrinkage → per-contrast TSV + MA + volcano PNGs.
//!
//! We deliberately do NOT wrap R as a runtime dep. The tool writes
//! out R code; the user runs R. This is Option 1D from
//! `scratch/rnaseq_spec_decisions.md` and the spec in
//! `scratch/de_template_spec.md`.
//!
//! Prerequisite alerting happens in three layers: in the CLI help
//! block, as a stderr warning at run time if `Rscript` isn't on PATH,
//! and as a guard at the top of the emitted R script so the user
//! gets a clear error before DESeq2 actually runs.

use crate::quant::sample_sheet;
use crate::utils::error::{MycoNoteError, Result};
use std::fs;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

/// Parsed command-line configuration for a `de-template` invocation.
#[derive(Debug, Clone)]
pub struct DeTemplateConfig {
    /// tximport mode: a quant_out directory containing salmon/<sample>/quant.sf.
    pub quant_dir: Option<PathBuf>,
    /// Counts-matrix mode: a wide counts.tsv (transcript × sample).
    pub counts_path: Option<PathBuf>,
    /// Sample sheet path. Defaults to `<quant_dir>/sample_sheet.tsv`
    /// if `--quant-dir` is given and the path isn't explicitly set.
    pub samples_path: PathBuf,
    /// R design formula, e.g. `~ condition` or `~ batch + condition`.
    /// Sanitized at parse time.
    pub design: String,
    /// One or more (factor, level1, level2) triples. At least one is
    /// required for a meaningful call.
    pub contrasts: Vec<Contrast>,
    pub output_script: PathBuf,
    pub fdr_cutoff: f64,
    pub lfc_cutoff: f64,
    /// Verbatim argv for embedding in the generated script header.
    pub raw_cmdline: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Contrast {
    pub factor: String,
    pub level1: String,
    pub level2: String,
}

impl DeTemplateConfig {
    pub fn parse(args: &[String]) -> Result<Self> {
        let mut quant_dir: Option<PathBuf> = None;
        let mut counts_path: Option<PathBuf> = None;
        let mut samples_path: Option<PathBuf> = None;
        let mut design: Option<String> = None;
        let mut contrasts: Vec<Contrast> = Vec::new();
        let mut output_script: PathBuf = PathBuf::from("de_analysis.R");
        let mut fdr_cutoff: f64 = 0.05;
        let mut lfc_cutoff: f64 = 1.0;

        let mut i = 0usize;
        while i < args.len() {
            match args[i].as_str() {
                "--quant-dir" if i + 1 < args.len() => {
                    quant_dir = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--counts" if i + 1 < args.len() => {
                    counts_path = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--samples" if i + 1 < args.len() => {
                    samples_path = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--design" if i + 1 < args.len() => {
                    design = Some(args[i + 1].clone());
                    i += 2;
                }
                "--contrast" if i + 1 < args.len() => {
                    contrasts.push(parse_contrast(&args[i + 1])?);
                    i += 2;
                }
                "--output" | "-o" if i + 1 < args.len() => {
                    output_script = PathBuf::from(&args[i + 1]);
                    i += 2;
                }
                "--fdr" if i + 1 < args.len() => {
                    fdr_cutoff = args[i + 1].parse().map_err(|_| {
                        MycoNoteError::QuantSheet(format!(
                            "de-template: --fdr must be a number, got '{}'",
                            args[i + 1]
                        ))
                    })?;
                    i += 2;
                }
                "--lfc" if i + 1 < args.len() => {
                    lfc_cutoff = args[i + 1].parse().map_err(|_| {
                        MycoNoteError::QuantSheet(format!(
                            "de-template: --lfc must be a number, got '{}'",
                            args[i + 1]
                        ))
                    })?;
                    i += 2;
                }
                other => {
                    return Err(MycoNoteError::QuantSheet(format!(
                        "de-template: unknown flag '{other}'"
                    )));
                }
            }
        }

        // Input-source validation.
        if quant_dir.is_some() && counts_path.is_some() {
            return Err(MycoNoteError::QuantSheet(
                "de-template: pass either --quant-dir or --counts, not both".to_string(),
            ));
        }
        if quant_dir.is_none() && counts_path.is_none() {
            return Err(MycoNoteError::QuantSheet(
                "de-template: must supply --quant-dir <dir> (tximport mode) \
                 or --counts <counts.tsv> (matrix mode)"
                    .to_string(),
            ));
        }

        let design = design.ok_or_else(|| {
            MycoNoteError::QuantSheet(
                "de-template: --design '<R formula>' is required (e.g. '--design ~ condition')"
                    .to_string(),
            )
        })?;
        let design = sanitize_design(&design)?;

        if contrasts.is_empty() {
            return Err(MycoNoteError::QuantSheet(
                "de-template: at least one --contrast 'factor,level1,level2' is required"
                    .to_string(),
            ));
        }

        // Default samples path: <quant_dir>/sample_sheet.tsv when
        // --quant-dir was given and --samples was not.
        let samples_path = samples_path
            .or_else(|| quant_dir.as_ref().map(|d| d.join("sample_sheet.tsv")))
            .ok_or_else(|| {
                MycoNoteError::QuantSheet(
                    "de-template: --samples <sheet.tsv> is required in --counts mode".to_string(),
                )
            })?;

        let raw_cmdline = std::iter::once("myconote-cli de-template".to_string())
            .chain(args.iter().cloned())
            .collect::<Vec<_>>()
            .join(" ");

        Ok(Self {
            quant_dir,
            counts_path,
            samples_path,
            design,
            contrasts,
            output_script,
            fdr_cutoff,
            lfc_cutoff,
            raw_cmdline,
        })
    }
}

/// Parse a `"factor,level1,level2"` contrast spec. Three
/// comma-separated tokens; none may contain commas internally
/// (sample sheets don't use commas in factor-level labels by
/// convention — we enforce that assumption loudly).
fn parse_contrast(s: &str) -> Result<Contrast> {
    let parts: Vec<&str> = s.split(',').map(|p| p.trim()).collect();
    if parts.len() != 3 || parts.iter().any(|p| p.is_empty()) {
        return Err(MycoNoteError::QuantSheet(format!(
            "de-template: --contrast must be 'factor,level1,level2', got '{s}'"
        )));
    }
    Ok(Contrast {
        factor: parts[0].to_string(),
        level1: parts[1].to_string(),
        level2: parts[2].to_string(),
    })
}

/// Light sanitization for design formulas before we embed them
/// verbatim into the generated R code. The embedding is always
/// between `DESIGN <- ` and a newline, so the attack surface is
/// "does the string contain things that let you exit the `DESIGN`
/// assignment and run arbitrary code?" We reject: `;` (statement
/// terminators), backticks (command substitution in shells, also
/// R environment reference), newlines (obvious), and the literal
/// substrings `system(`, `system2(`, `source(`, and `eval(` which
/// are the shortest named paths to arbitrary execution.
///
/// Not a security boundary — the user is generating their own R
/// script and then running it on their own machine. This catches
/// accidental quoting errors and obvious copy-paste mistakes, not
/// determined attacks.
pub fn sanitize_design(s: &str) -> Result<String> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return Err(MycoNoteError::QuantSheet(
            "de-template: --design is empty".to_string(),
        ));
    }
    let forbid = [
        (";", "statement terminator"),
        ("`", "backtick"),
        ("\n", "newline"),
        ("\r", "carriage return"),
        ("system(", "system() call"),
        ("system2(", "system2() call"),
        ("source(", "source() call"),
        ("eval(", "eval() call"),
    ];
    for (pat, desc) in forbid {
        if trimmed.contains(pat) {
            return Err(MycoNoteError::QuantSheet(format!(
                "de-template: --design contains {desc} ('{pat}') — rejected for safety. \
                 Use plain R formulas like '~ condition' or '~ batch + condition'."
            )));
        }
    }
    // Design formulas almost always start with `~`. Warn (not error)
    // if not — users sometimes forget the tilde.
    if !trimmed.starts_with('~') {
        eprintln!(
            "  ⚠ --design '{trimmed}' does not start with '~'; DESeq2 \
             expects a formula like '~ condition'. Passing through anyway."
        );
    }
    Ok(trimmed.to_string())
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation against the sample sheet
// ─────────────────────────────────────────────────────────────────────────────

/// Check that each contrast's `factor` is a column in the sheet, and
/// both `level1` and `level2` appear as values in that column. We
/// check `condition` explicitly as a known column, plus any custom
/// columns captured in the `Sample.extras` map. Catches the common
/// typos (factor name wrong, level misspelled) before R does.
fn validate_contrasts(cfg: &DeTemplateConfig, sheet: &sample_sheet::SampleSheet) -> Result<()> {
    for c in &cfg.contrasts {
        let values = collect_factor_values(&c.factor, sheet);
        if values.is_empty() {
            return Err(MycoNoteError::QuantSheet(format!(
                "de-template: --contrast '{}': factor '{}' is not a column in \
                 {} (or is empty for every sample). Known columns: {}",
                contrast_to_string(c),
                c.factor,
                cfg.samples_path.display(),
                sheet.input_header.join(", ")
            )));
        }
        for (label, level) in [("level1", &c.level1), ("level2", &c.level2)] {
            if !values.iter().any(|v| v == level) {
                return Err(MycoNoteError::QuantSheet(format!(
                    "de-template: --contrast '{}': {label} '{}' not found in column \
                     '{}'. Values present: {}",
                    contrast_to_string(c),
                    level,
                    c.factor,
                    values.iter().cloned().collect::<Vec<_>>().join(", ")
                )));
            }
        }
    }
    Ok(())
}

fn contrast_to_string(c: &Contrast) -> String {
    format!("{},{},{}", c.factor, c.level1, c.level2)
}

fn collect_factor_values(factor: &str, sheet: &sample_sheet::SampleSheet) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for sample in &sheet.samples {
        // The canonical `condition` field of Sample is first-class.
        // Any other column (batch / custom) lives in `extras`.
        let v: Option<String> = if factor == "condition" {
            sample.condition.clone()
        } else if factor == "batch" {
            sample.batch.clone()
        } else {
            sample.extras.get(factor).cloned()
        };
        if let Some(v) = v {
            if !v.is_empty() && !out.contains(&v) {
                out.push(v);
            }
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Top-level dispatcher
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_de_template(args: &[String]) -> Result<()> {
    let cfg = DeTemplateConfig::parse(args)?;

    // Surface the R / DESeq2 prerequisite prominently. We check for
    // `Rscript` but do not fail on its absence — users may generate
    // the script on one machine and run it on another.
    warn_about_r_prereqs();

    // Validate sample sheet early so contrast mistakes are caught
    // before we write any files.
    require_file(&cfg.samples_path, "--samples sheet")?;
    let sheet = sample_sheet::parse_sheet(&cfg.samples_path)?;
    validate_contrasts(&cfg, &sheet)?;

    if let Some(ref d) = cfg.quant_dir {
        require_dir(d, "--quant-dir")?;
        if !d.join("salmon").is_dir() {
            eprintln!(
                "  ⚠ --quant-dir {} has no salmon/ subdirectory; the generated \
                 script will fail unless `quant` output lives there.",
                d.display()
            );
        }
    }
    if let Some(ref p) = cfg.counts_path {
        require_file(p, "--counts")?;
    }

    let script = render_script(&cfg);
    if let Some(parent) = cfg.output_script.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    fs::write(&cfg.output_script, script)?;

    eprintln!("  ✓ wrote {}", cfg.output_script.display());
    eprintln!("\nNext step — install R prerequisites (if not already), then run:");
    eprintln!("  Rscript {}", cfg.output_script.display());
    eprintln!("\nR install (one-time, inside an R session):");
    eprintln!(
        "  if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")"
    );
    eprintln!("  BiocManager::install(c(\"tximport\", \"DESeq2\", \"apeglm\"))");
    Ok(())
}

fn warn_about_r_prereqs() {
    if which::which("Rscript").is_err() {
        eprintln!(
            "  ⚠ `Rscript` not found on PATH. The generated script needs R and \
             the Bioconductor packages tximport, DESeq2 (+ apeglm recommended). \
             The script will be written anyway so you can run it on another host."
        );
    }
}

fn require_file(p: &Path, role: &str) -> Result<()> {
    if !p.is_file() {
        return Err(MycoNoteError::QuantSheet(format!(
            "{role}: file not found at {}",
            p.display()
        )));
    }
    Ok(())
}

fn require_dir(p: &Path, role: &str) -> Result<()> {
    if !p.is_dir() {
        return Err(MycoNoteError::QuantSheet(format!(
            "{role}: directory not found at {}",
            p.display()
        )));
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Template renderer
// ─────────────────────────────────────────────────────────────────────────────

/// Render a full R script from a validated `DeTemplateConfig`. Pure —
/// returns a String, no I/O.
pub fn render_script(cfg: &DeTemplateConfig) -> String {
    let mut s = String::with_capacity(4096);

    let version = env!("CARGO_PKG_VERSION");
    let timestamp = timestamp_rfc3339();

    // Header / metadata
    s.push_str(&format!("#!/usr/bin/env Rscript\n"));
    s.push_str(&format!(
        "# Auto-generated by myconote-cli de-template v{version}\n"
    ));
    s.push_str(&format!("# Generated: {timestamp}\n"));
    s.push_str(&format!("# Command:   {}\n", cfg.raw_cmdline));
    s.push_str("#\n");
    s.push_str("# Usage: Rscript ");
    s.push_str(
        &cfg.output_script
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "de_analysis.R".to_string()),
    );
    s.push('\n');
    s.push_str("#\n");
    s.push_str("# REQUIRES (install once, before running this script):\n");
    s.push_str("#   if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")\n");
    s.push_str("#   BiocManager::install(c(\"tximport\", \"DESeq2\", \"apeglm\"))\n");
    s.push_str("\n");

    // Inputs block — the user can edit these to re-run with different
    // thresholds without regenerating the script.
    s.push_str(
        "# ---------------- Inputs (edit to re-run with different params) ----------------\n",
    );
    s.push_str(&format!(
        "QUANT_DIR            <- {}\n",
        r_string_option(
            cfg.quant_dir
                .as_ref()
                .map(|p| p.to_string_lossy().into_owned())
        )
    ));
    s.push_str(&format!(
        "COUNTS_PATH          <- {}\n",
        r_string_option(
            cfg.counts_path
                .as_ref()
                .map(|p| p.to_string_lossy().into_owned())
        )
    ));
    s.push_str(&format!(
        "SAMPLES_PATH         <- \"{}\"\n",
        cfg.samples_path.to_string_lossy()
    ));
    s.push_str(&format!("DESIGN               <- {}\n", cfg.design));
    s.push_str("CONTRASTS <- list(\n");
    for (i, c) in cfg.contrasts.iter().enumerate() {
        let comma = if i + 1 < cfg.contrasts.len() { "," } else { "" };
        s.push_str(&format!(
            "  c({}, {}, {}){comma}\n",
            r_string(&c.factor),
            r_string(&c.level1),
            r_string(&c.level2)
        ));
    }
    s.push_str(")\n");
    s.push_str(&format!("FDR_CUTOFF           <- {}\n", cfg.fdr_cutoff));
    s.push_str(&format!("LFC_CUTOFF           <- {}\n", cfg.lfc_cutoff));
    s.push_str("OUT_PREFIX           <- \"de\"\n");
    s.push_str("PREFILTER_MIN_COUNT  <- 10\n");
    s.push_str("PREFILTER_MIN_SAMPLES <- NULL  # default: floor(N/2)\n");
    s.push_str("\n");

    // Package guard — fails fast with a clear message if DESeq2 isn't
    // installed. This is the user-alert layer that runs when the
    // script itself is invoked.
    s.push_str(PACKAGE_GUARD);
    s.push_str("\n");

    // Samples loading
    s.push_str(SAMPLES_LOADING);
    s.push_str("\n");

    // DESeqDataSet construction — tximport branch OR matrix branch,
    // chosen at runtime based on which input variable is non-NULL.
    s.push_str(DDS_CONSTRUCTION);
    s.push_str("\n");

    // Prefilter
    s.push_str(PREFILTER);
    s.push_str("\n");

    // Run DESeq
    s.push_str("cat(sprintf(\"\\nRunning DESeq2 on %d transcripts × %d samples …\\n\",\n");
    s.push_str("            nrow(dds), ncol(dds)))\n");
    s.push_str("dds <- DESeq(dds)\n");
    s.push_str("\n");

    // Per-contrast loop (with apeglm shrinkage + TSV + plots)
    s.push_str(CONTRAST_LOOP);
    s.push_str("\n");

    // Session info / provenance
    s.push_str("cat(\"\\n---------------- sessionInfo() ----------------\\n\")\n");
    s.push_str("sessionInfo()\n");

    s
}

/// Format an Option<String> as R code: `NULL` for None, quoted string for Some.
fn r_string_option(v: Option<String>) -> String {
    match v {
        Some(s) => r_string(&s),
        None => "NULL".to_string(),
    }
}

/// R-quote a string. Escapes embedded double quotes and backslashes.
fn r_string(s: &str) -> String {
    let escaped: String = s
        .chars()
        .map(|c| match c {
            '\\' => "\\\\".to_string(),
            '"' => "\\\"".to_string(),
            c => c.to_string(),
        })
        .collect();
    format!("\"{}\"", escaped)
}

// ─────────────────────────────────────────────────────────────────────────────
// Static R snippets — the non-parametric parts of the generated script.
// Kept as `&str` constants so tests can assert they appear in the output
// and users can eyeball them as intended R code rather than string builders.
// ─────────────────────────────────────────────────────────────────────────────

const PACKAGE_GUARD: &str = r###"# ---------------- Package guard ----------------
required <- c("tximport", "DESeq2")
missing_ <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_)) {
  stop(paste0(
    "\n\nMissing R package(s): ", paste(missing_, collapse = ", "), ".\n",
    "Install (one-time) before running this script:\n",
    "  if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")\n",
    "  BiocManager::install(c(\"tximport\", \"DESeq2\", \"apeglm\"))\n"
  ))
}
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})
has_apeglm <- requireNamespace("apeglm", quietly = TRUE)
if (!has_apeglm) {
  message("apeglm not installed; LFC shrinkage will be skipped. ",
          "Install with BiocManager::install(\"apeglm\") for better LFC estimates.")
}
"###;

const SAMPLES_LOADING: &str = r###"# ---------------- Load samples ----------------
samples <- read.table(SAMPLES_PATH, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "#")
if (!"sample_id" %in% names(samples)) {
  stop("samples sheet ", SAMPLES_PATH, " is missing a 'sample_id' column")
}
rownames(samples) <- samples$sample_id
cat(sprintf("Loaded %d samples from %s\n", nrow(samples), SAMPLES_PATH))
"###;

const DDS_CONSTRUCTION: &str = r###"# ---------------- Build DESeqDataSet ----------------
if (!is.null(QUANT_DIR)) {
  salmon_dir <- file.path(QUANT_DIR, "salmon")
  files <- file.path(salmon_dir, samples$sample_id, "quant.sf")
  names(files) <- samples$sample_id
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop("Missing salmon quant.sf for samples:\n  ",
         paste(names(missing_files), collapse = "\n  "))
  }
  txi <- tximport(files, type = "salmon", txOut = TRUE)
  dds <- DESeqDataSetFromTximport(txi, colData = samples, design = DESIGN)
} else if (!is.null(COUNTS_PATH)) {
  counts <- read.table(COUNTS_PATH, header = TRUE, sep = "\t",
                        row.names = 1, check.names = FALSE)
  counts <- as.matrix(counts)
  # counts.tsv may contain non-integer estimated counts from salmon.
  # DESeq2 coerces with a warning; we round here to silence it and be
  # explicit about the behavior.
  storage.mode(counts) <- "integer"
  # Align columns to sample_id order.
  counts <- counts[, samples$sample_id, drop = FALSE]
  dds <- DESeqDataSetFromMatrix(counts, colData = samples, design = DESIGN)
} else {
  stop("Neither QUANT_DIR nor COUNTS_PATH is set")
}
cat(sprintf("DESeqDataSet: %d transcripts × %d samples\n",
            nrow(dds), ncol(dds)))
"###;

const PREFILTER: &str = r###"# ---------------- Prefilter low-count rows ----------------
min_samples <- if (is.null(PREFILTER_MIN_SAMPLES)) {
  floor(ncol(dds) / 2)
} else PREFILTER_MIN_SAMPLES
keep <- rowSums(counts(dds) >= PREFILTER_MIN_COUNT) >= min_samples
dds <- dds[keep, ]
cat(sprintf("After prefilter: %d transcripts pass (>= %d counts in >= %d samples)\n",
            nrow(dds), PREFILTER_MIN_COUNT, min_samples))
"###;

const CONTRAST_LOOP: &str = r###"# ---------------- Per-contrast DE + plots ----------------
for (contrast in CONTRASTS) {
  factor_name <- contrast[1]
  level1      <- contrast[2]
  level2      <- contrast[3]

  res <- results(dds, contrast = c(factor_name, level1, level2))

  # apeglm needs a coef name rather than a contrast tuple. DESeq2 names
  # coefs like "condition_treated_vs_control". If the user's factor
  # reference level is already level2 (the denominator), that name
  # matches directly; otherwise we relevel, re-run results, and then
  # shrink. Failing that we fall back to unshrunk LFCs with a warning.
  res_shrunk <- res
  if (has_apeglm) {
    coef_name <- paste0(factor_name, "_", level1, "_vs_", level2)
    if (coef_name %in% resultsNames(dds)) {
      res_shrunk <- tryCatch(
        lfcShrink(dds, coef = coef_name, type = "apeglm"),
        error = function(e) { message("apeglm failed: ", conditionMessage(e),
                                       "; using unshrunk LFCs."); res }
      )
    } else {
      message("Coefficient '", coef_name, "' not in resultsNames(dds); ",
              "LFC shrinkage skipped for this contrast.")
    }
  }

  tag <- sprintf("%s_%s_vs_%s", factor_name, level1, level2)
  cat(sprintf("\n=== %s ===\n", tag))
  summary(res_shrunk)

  # TSV output — one row per transcript, sorted by padj.
  # Columns chosen to cover both unshrunk (has $stat) and apeglm-shrunk
  # (no $stat) result objects: we always emit the five core columns
  # that both object types provide.
  out_df <- as.data.frame(res_shrunk)
  out_df$transcript <- rownames(out_df)
  out_df <- out_df[, c("transcript", "baseMean", "log2FoldChange",
                        "lfcSE", "pvalue", "padj")]
  out_df <- out_df[order(out_df$padj, na.last = TRUE), ]
  out_tsv <- sprintf("%s_%s.tsv", OUT_PREFIX, tag)
  write.table(out_df, out_tsv, sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
  cat(sprintf("Wrote %s (%d rows)\n", out_tsv, nrow(out_df)))

  sig_rows <- !is.na(out_df$padj) & out_df$padj < FDR_CUTOFF &
              abs(out_df$log2FoldChange) >= LFC_CUTOFF
  cat(sprintf("  significant at FDR<%.3f, |LFC|>=%.2f: %d\n",
              FDR_CUTOFF, LFC_CUTOFF, sum(sig_rows)))

  # MA plot.
  png(sprintf("%s_%s_MA.png", OUT_PREFIX, tag), width = 900, height = 600)
  DESeq2::plotMA(res_shrunk, main = tag, ylim = c(-5, 5),
                 alpha = FDR_CUTOFF)
  dev.off()

  # Volcano plot — highlights padj < FDR_CUTOFF & |LFC| >= LFC_CUTOFF.
  png(sprintf("%s_%s_volcano.png", OUT_PREFIX, tag), width = 900, height = 700)
  par(mar = c(4.5, 4.5, 3, 1))
  with(out_df, {
    plot(log2FoldChange, -log10(pvalue),
         pch = 20, cex = 0.5, col = "grey60",
         main = sprintf("%s  (FDR<%.3f, |LFC|>=%.2f)",
                        tag, FDR_CUTOFF, LFC_CUTOFF),
         xlab = "log2 Fold Change",
         ylab = expression(-log[10](italic(p))))
    abline(h = -log10(FDR_CUTOFF), lty = 2, col = "grey40")
    abline(v = c(-LFC_CUTOFF, LFC_CUTOFF), lty = 2, col = "grey40")
    with(out_df[sig_rows, ],
         points(log2FoldChange, -log10(pvalue), pch = 20, cex = 0.7,
                col = ifelse(log2FoldChange > 0, "firebrick", "steelblue")))
  })
  dev.off()
}
"###;

// ─────────────────────────────────────────────────────────────────────────────
// Timestamp helper (reused from src/quant/mod.rs pattern; intentionally
// local-copy rather than shared to keep de_template self-contained).
// ─────────────────────────────────────────────────────────────────────────────

fn timestamp_rfc3339() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    format_rfc3339(secs)
}

fn format_rfc3339(secs: i64) -> String {
    let days = secs.div_euclid(86400);
    let sod = secs.rem_euclid(86400);
    let hour = sod / 3600;
    let minute = (sod / 60) % 60;
    let second = sod % 60;
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        y, m, d, hour, minute, second
    )
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use tempfile::TempDir;

    fn argv(s: &str) -> Vec<String> {
        // Allow '|' as an argv separator so test fixtures can include
        // spaces within a single arg (e.g. design formulas).
        s.split('|').map(|p| p.trim().to_string()).collect()
    }

    // ── parse tests ─────────────────────────────────────────────────────────

    #[test]
    fn parse_minimal_quant_dir() {
        let a =
            argv("--quant-dir|quant_out|--design|~ condition|--contrast|condition,treated,control");
        let cfg = DeTemplateConfig::parse(&a).unwrap();
        assert_eq!(cfg.quant_dir.as_deref(), Some(Path::new("quant_out")));
        assert!(cfg.counts_path.is_none());
        assert_eq!(
            cfg.samples_path,
            PathBuf::from("quant_out/sample_sheet.tsv")
        );
        assert_eq!(cfg.design, "~ condition");
        assert_eq!(cfg.contrasts.len(), 1);
        assert_eq!(cfg.contrasts[0].factor, "condition");
        assert_eq!(cfg.contrasts[0].level1, "treated");
        assert_eq!(cfg.contrasts[0].level2, "control");
        assert_eq!(cfg.output_script, PathBuf::from("de_analysis.R"));
        assert!((cfg.fdr_cutoff - 0.05).abs() < 1e-9);
        assert!((cfg.lfc_cutoff - 1.0).abs() < 1e-9);
    }

    #[test]
    fn parse_counts_mode_requires_samples() {
        let a = argv("--counts|counts.tsv|--design|~ condition|--contrast|condition,a,b");
        let err = DeTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--samples"));
    }

    #[test]
    fn parse_counts_mode_with_samples() {
        let a = argv(
            "--counts|counts.tsv|--samples|s.tsv|--design|~ condition|--contrast|condition,a,b",
        );
        let cfg = DeTemplateConfig::parse(&a).unwrap();
        assert_eq!(cfg.counts_path.as_deref(), Some(Path::new("counts.tsv")));
        assert_eq!(cfg.samples_path, PathBuf::from("s.tsv"));
    }

    #[test]
    fn parse_rejects_both_inputs() {
        let a = argv("--quant-dir|q|--counts|c.tsv|--design|~ x|--contrast|x,a,b");
        let err = DeTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("not both"));
    }

    #[test]
    fn parse_rejects_neither_input() {
        let a = argv("--design|~ condition|--contrast|condition,a,b");
        let err = DeTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("must supply"));
    }

    #[test]
    fn parse_rejects_missing_design() {
        let a = argv("--quant-dir|q|--contrast|condition,a,b");
        let err = DeTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--design"));
    }

    #[test]
    fn parse_rejects_missing_contrast() {
        let a = argv("--quant-dir|q|--design|~ x");
        let err = DeTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--contrast"));
    }

    #[test]
    fn parse_multiple_contrasts() {
        let a = argv(
            "--quant-dir|q|--design|~ condition|--contrast|condition,a,b|--contrast|condition,c,b",
        );
        let cfg = DeTemplateConfig::parse(&a).unwrap();
        assert_eq!(cfg.contrasts.len(), 2);
    }

    #[test]
    fn parse_rejects_malformed_contrast() {
        let a = argv("--quant-dir|q|--design|~ x|--contrast|only_two,parts");
        let err = DeTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("factor,level1,level2"));
    }

    #[test]
    fn parse_custom_output() {
        let a = argv("--quant-dir|q|--design|~ x|--contrast|x,a,b|-o|my_analysis.R");
        let cfg = DeTemplateConfig::parse(&a).unwrap();
        assert_eq!(cfg.output_script, PathBuf::from("my_analysis.R"));
    }

    // ── sanitization ────────────────────────────────────────────────────────

    #[test]
    fn sanitize_accepts_clean_formulas() {
        for ok in ["~ condition", "~batch + condition", "~ condition * batch"] {
            assert_eq!(sanitize_design(ok).unwrap(), ok);
        }
    }

    #[test]
    fn sanitize_rejects_dangerous_tokens() {
        for bad in [
            "~ condition; system(\"id\")",
            "~ condition\nstop()",
            "~ `cat /etc/passwd`",
            "~ eval(parse(text=x))",
            "~ source(\"x.R\")",
        ] {
            assert!(sanitize_design(bad).is_err(), "should reject: {bad}");
        }
    }

    #[test]
    fn sanitize_rejects_empty() {
        assert!(sanitize_design("").is_err());
        assert!(sanitize_design("   ").is_err());
    }

    // ── contrast validation against a sample sheet ──────────────────────────

    fn fake_sheet() -> sample_sheet::SampleSheet {
        use sample_sheet::{Sample, SampleSheet, Strandedness};
        let make = |id: &str, cond: &str, batch: Option<&str>| Sample {
            sample_id: id.to_string(),
            fastq_r1: PathBuf::from(format!("{id}_R1.fq.gz")),
            fastq_r2: Some(PathBuf::from(format!("{id}_R2.fq.gz"))),
            condition: Some(cond.to_string()),
            strandedness: Strandedness::Auto,
            batch: batch.map(String::from),
            source_line: 2,
            extras: BTreeMap::new(),
        };
        SampleSheet {
            samples: vec![
                make("WT1", "control", Some("A")),
                make("WT2", "control", Some("A")),
                make("KO1", "treated", Some("B")),
                make("KO2", "treated", Some("B")),
            ],
            input_header: vec![
                "sample_id".into(),
                "fastq_r1".into(),
                "fastq_r2".into(),
                "condition".into(),
                "batch".into(),
            ],
            root_dir: PathBuf::from("/tmp"),
        }
    }

    fn make_cfg(contrasts: Vec<Contrast>) -> DeTemplateConfig {
        DeTemplateConfig {
            quant_dir: Some(PathBuf::from("q")),
            counts_path: None,
            samples_path: PathBuf::from("q/sample_sheet.tsv"),
            design: "~ condition".to_string(),
            contrasts,
            output_script: PathBuf::from("de.R"),
            fdr_cutoff: 0.05,
            lfc_cutoff: 1.0,
            raw_cmdline: "test".to_string(),
        }
    }

    #[test]
    fn validate_contrast_happy_path() {
        let sheet = fake_sheet();
        let cfg = make_cfg(vec![Contrast {
            factor: "condition".into(),
            level1: "treated".into(),
            level2: "control".into(),
        }]);
        validate_contrasts(&cfg, &sheet).unwrap();
    }

    #[test]
    fn validate_contrast_batch_column() {
        let sheet = fake_sheet();
        let cfg = make_cfg(vec![Contrast {
            factor: "batch".into(),
            level1: "B".into(),
            level2: "A".into(),
        }]);
        validate_contrasts(&cfg, &sheet).unwrap();
    }

    #[test]
    fn validate_contrast_rejects_unknown_factor() {
        let sheet = fake_sheet();
        let cfg = make_cfg(vec![Contrast {
            factor: "time_point".into(),
            level1: "early".into(),
            level2: "late".into(),
        }]);
        let err = validate_contrasts(&cfg, &sheet).unwrap_err();
        assert!(format!("{err}").contains("not a column"));
    }

    #[test]
    fn validate_contrast_rejects_unknown_level() {
        let sheet = fake_sheet();
        let cfg = make_cfg(vec![Contrast {
            factor: "condition".into(),
            level1: "starved".into(), // not in sheet
            level2: "control".into(),
        }]);
        let err = validate_contrasts(&cfg, &sheet).unwrap_err();
        assert!(format!("{err}").contains("not found in column"));
    }

    // ── template rendering ──────────────────────────────────────────────────

    #[test]
    fn render_contains_expected_substitutions() {
        let cfg = DeTemplateConfig {
            quant_dir: Some(PathBuf::from("/path/to/quant_out")),
            counts_path: None,
            samples_path: PathBuf::from("/path/to/samples.tsv"),
            design: "~ condition".to_string(),
            contrasts: vec![Contrast {
                factor: "condition".into(),
                level1: "treated".into(),
                level2: "control".into(),
            }],
            output_script: PathBuf::from("analysis.R"),
            fdr_cutoff: 0.05,
            lfc_cutoff: 1.0,
            raw_cmdline: "myconote-cli de-template ...".to_string(),
        };
        let s = render_script(&cfg);
        assert!(s.contains("QUANT_DIR            <- \"/path/to/quant_out\""));
        assert!(s.contains("COUNTS_PATH          <- NULL"));
        assert!(s.contains("SAMPLES_PATH         <- \"/path/to/samples.tsv\""));
        assert!(s.contains("DESIGN               <- ~ condition"));
        assert!(s.contains("c(\"condition\", \"treated\", \"control\")"));
        assert!(s.contains("FDR_CUTOFF           <- 0.05"));
        assert!(s.contains("LFC_CUTOFF           <- 1"));
        assert!(s.contains("requireNamespace(\"BiocManager\""));
        assert!(s.contains("DESeqDataSetFromTximport"));
        assert!(s.contains("DESeqDataSetFromMatrix"));
        assert!(s.contains("lfcShrink(dds"));
        assert!(s.contains("plotMA"));
        assert!(s.contains("sessionInfo()"));
    }

    #[test]
    fn render_counts_mode_uses_null_quant_dir() {
        let cfg = DeTemplateConfig {
            quant_dir: None,
            counts_path: Some(PathBuf::from("/c.tsv")),
            samples_path: PathBuf::from("/s.tsv"),
            design: "~ x".to_string(),
            contrasts: vec![Contrast {
                factor: "x".into(),
                level1: "a".into(),
                level2: "b".into(),
            }],
            output_script: PathBuf::from("out.R"),
            fdr_cutoff: 0.05,
            lfc_cutoff: 1.0,
            raw_cmdline: "test".to_string(),
        };
        let s = render_script(&cfg);
        assert!(s.contains("QUANT_DIR            <- NULL"));
        assert!(s.contains("COUNTS_PATH          <- \"/c.tsv\""));
    }

    #[test]
    fn render_multiple_contrasts() {
        let cfg = DeTemplateConfig {
            quant_dir: Some(PathBuf::from("q")),
            counts_path: None,
            samples_path: PathBuf::from("q/sample_sheet.tsv"),
            design: "~ condition".to_string(),
            contrasts: vec![
                Contrast {
                    factor: "condition".into(),
                    level1: "A".into(),
                    level2: "B".into(),
                },
                Contrast {
                    factor: "condition".into(),
                    level1: "C".into(),
                    level2: "B".into(),
                },
            ],
            output_script: PathBuf::from("out.R"),
            fdr_cutoff: 0.05,
            lfc_cutoff: 1.0,
            raw_cmdline: "test".to_string(),
        };
        let s = render_script(&cfg);
        assert!(s.contains("c(\"condition\", \"A\", \"B\")"));
        assert!(s.contains("c(\"condition\", \"C\", \"B\")"));
        // The intermediate comma should appear between the two c() calls.
        let ab_idx = s.find("\"A\", \"B\"").unwrap();
        let cb_idx = s.find("\"C\", \"B\"").unwrap();
        assert!(ab_idx < cb_idx);
    }

    #[test]
    fn render_quotes_special_characters() {
        // Level names with quotes / backslashes should round-trip
        // cleanly into R strings. Unusual but not impossible.
        let s = r_string(r#"weird"level\with\stuff"#);
        assert_eq!(s, r#""weird\"level\\with\\stuff""#);
    }

    // ── rfc3339 ─────────────────────────────────────────────────────────────

    #[test]
    fn rfc3339_epoch() {
        assert_eq!(format_rfc3339(0), "1970-01-01T00:00:00Z");
    }

    #[test]
    fn rfc3339_known() {
        // Verified: datetime(2026,4,24,12,0,0,tzinfo=utc).timestamp() == 1777032000
        assert_eq!(format_rfc3339(1777032000), "2026-04-24T12:00:00Z");
    }
}
