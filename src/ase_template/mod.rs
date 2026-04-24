//! Allele-specific expression R-script generator.
//!
//! `myconote-cli ase-template` takes the outputs of
//! `myconote-cli ase` (`ase_counts.tsv` + `ase_summary.tsv`) and emits
//! a self-contained R script. The user runs that script with
//! `Rscript`; it performs a per-transcript **binomial exact test** on
//! `(count_hap0, count_hap1)` with a null proportion derived from the
//! per-sample total haplotype library-size ratio, BH-adjusts across
//! transcripts, and plots the allelic-imbalance distribution.
//!
//! Deliberately a sibling of `de-template` — same Option 1D pattern
//! (write R, user runs R, no Rust-side R runtime). The statistics are
//! different enough (binomial exact, not NB-GLM) that mixing them
//! under one subcommand would hurt clarity. See
//! `scratch/ase_spec.md` Q6 for the decision.
//!
//! Prerequisite alerting mirrors `de-template`'s three-layer surface:
//! CLI help, runtime stderr if `Rscript` is missing, and an in-script
//! package guard so the R process fails fast with a clear message.

use crate::utils::error::{MycoNoteError, Result};
use std::fs;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

/// Parsed command-line configuration for an `ase-template` invocation.
#[derive(Debug, Clone)]
pub struct AseTemplateConfig {
    /// Path to `ase_counts.tsv` produced by `myconote-cli ase`. First
    /// column is `transcript`; subsequent columns are `<sample>.<hap>`.
    pub counts_path: PathBuf,
    /// Path to `ase_summary.tsv` produced by `myconote-cli ase`. Used
    /// to filter to informative transcripts (those where at least one
    /// variant differs between haplotypes). Optional: if absent, the
    /// generated script tests every transcript.
    pub summary_path: Option<PathBuf>,
    pub output_script: PathBuf,
    /// Two haplotype names in the column suffixes. Default
    /// `["hap0", "hap1"]` matches the `ase` dispatcher default.
    pub haplotype_names: [String; 2],
    /// Benjamini–Hochberg significance threshold.
    pub fdr_cutoff: f64,
    /// Minimum total read count per transcript per sample to be
    /// included in the test. Below this the binomial test is
    /// uninformative (huge null variance).
    pub min_reads_per_sample: u32,
    /// Whether to skip uninformative transcripts (informative=false
    /// in `ase_summary.tsv`). Default true; `--include-uninformative`
    /// disables.
    pub filter_uninformative: bool,
    /// Verbatim argv for embedding in the generated script header.
    pub raw_cmdline: String,
}

impl AseTemplateConfig {
    pub fn parse(args: &[String]) -> Result<Self> {
        let mut counts_path: Option<PathBuf> = None;
        let mut summary_path: Option<PathBuf> = None;
        let mut output_script: PathBuf = PathBuf::from("ase_analysis.R");
        let mut haplotype_names = ["hap0".to_string(), "hap1".to_string()];
        let mut fdr_cutoff: f64 = 0.05;
        let mut min_reads_per_sample: u32 = 20;
        let mut filter_uninformative = true;

        let mut i = 0usize;
        while i < args.len() {
            match args[i].as_str() {
                "--counts" if i + 1 < args.len() => {
                    counts_path = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--summary" if i + 1 < args.len() => {
                    summary_path = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--ase-dir" if i + 1 < args.len() => {
                    // Shortcut: `--ase-dir ase_out` fills both counts
                    // + summary from their canonical positions.
                    let d = PathBuf::from(&args[i + 1]);
                    if counts_path.is_none() {
                        counts_path = Some(d.join("ase_counts.tsv"));
                    }
                    if summary_path.is_none() {
                        summary_path = Some(d.join("ase_summary.tsv"));
                    }
                    i += 2;
                }
                "--output" | "-o" if i + 1 < args.len() => {
                    output_script = PathBuf::from(&args[i + 1]);
                    i += 2;
                }
                "--haplotype-names" if i + 1 < args.len() => {
                    let parts: Vec<&str> = args[i + 1].split(',').collect();
                    if parts.len() != 2 {
                        return Err(MycoNoteError::QuantSheet(format!(
                            "ase-template: --haplotype-names expects 'NAME1,NAME2', got '{}'",
                            args[i + 1]
                        )));
                    }
                    haplotype_names = [parts[0].trim().to_string(), parts[1].trim().to_string()];
                    if haplotype_names[0] == haplotype_names[1] {
                        return Err(MycoNoteError::QuantSheet(
                            "ase-template: --haplotype-names: the two names must differ"
                                .to_string(),
                        ));
                    }
                    i += 2;
                }
                "--fdr" if i + 1 < args.len() => {
                    fdr_cutoff = args[i + 1].parse().map_err(|_| {
                        MycoNoteError::QuantSheet(format!(
                            "ase-template: --fdr must be a number, got '{}'",
                            args[i + 1]
                        ))
                    })?;
                    i += 2;
                }
                "--min-reads" if i + 1 < args.len() => {
                    min_reads_per_sample = args[i + 1].parse().map_err(|_| {
                        MycoNoteError::QuantSheet(format!(
                            "ase-template: --min-reads must be a non-negative integer, got '{}'",
                            args[i + 1]
                        ))
                    })?;
                    i += 2;
                }
                "--include-uninformative" => {
                    filter_uninformative = false;
                    i += 1;
                }
                other => {
                    return Err(MycoNoteError::QuantSheet(format!(
                        "ase-template: unknown flag '{other}'"
                    )));
                }
            }
        }

        let counts_path = counts_path.ok_or_else(|| {
            MycoNoteError::QuantSheet(
                "ase-template: --counts <ase_counts.tsv> is required (or \
                 use --ase-dir <dir> to pick up both counts + summary)."
                    .to_string(),
            )
        })?;

        let raw_cmdline = std::iter::once("myconote-cli ase-template".to_string())
            .chain(args.iter().cloned())
            .collect::<Vec<_>>()
            .join(" ");

        Ok(Self {
            counts_path,
            summary_path,
            output_script,
            haplotype_names,
            fdr_cutoff,
            min_reads_per_sample,
            filter_uninformative,
            raw_cmdline,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Top-level dispatcher
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_ase_template(args: &[String]) -> Result<()> {
    let cfg = AseTemplateConfig::parse(args)?;

    warn_about_r_prereqs();

    require_file(&cfg.counts_path, "--counts")?;
    if let Some(ref s) = cfg.summary_path {
        require_file(s, "--summary")?;
    } else if cfg.filter_uninformative {
        eprintln!(
            "  ⚠ no --summary given; the generated script cannot filter \
             uninformative transcripts. Pass --summary <ase_summary.tsv> \
             or --include-uninformative to acknowledge."
        );
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
    eprintln!("  # base R >= 4.0 is enough; no Bioconductor dependency for binomial ASE.");
    Ok(())
}

fn warn_about_r_prereqs() {
    if which::which("Rscript").is_err() {
        eprintln!(
            "  ⚠ `Rscript` not found on PATH. The generated script uses \
             only base R (no DESeq2 / Bioconductor needed). The script \
             will be written anyway so you can run it on another host."
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

// ─────────────────────────────────────────────────────────────────────────────
// Template renderer
// ─────────────────────────────────────────────────────────────────────────────

/// Render a full R script from a validated `AseTemplateConfig`. Pure —
/// returns a String, no I/O.
pub fn render_script(cfg: &AseTemplateConfig) -> String {
    let mut s = String::with_capacity(4096);

    let version = env!("CARGO_PKG_VERSION");
    let timestamp = timestamp_rfc3339();

    // Header
    s.push_str("#!/usr/bin/env Rscript\n");
    s.push_str(&format!(
        "# Auto-generated by myconote-cli ase-template v{version}\n"
    ));
    s.push_str(&format!("# Generated: {timestamp}\n"));
    s.push_str(&format!("# Command:   {}\n", cfg.raw_cmdline));
    s.push_str("#\n");
    s.push_str("# Usage: Rscript ");
    s.push_str(
        &cfg.output_script
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "ase_analysis.R".to_string()),
    );
    s.push('\n');
    s.push_str("#\n");
    s.push_str("# What this does:\n");
    s.push_str("#   - Reads ASE_COUNTS (transcript × sample.hap counts matrix)\n");
    s.push_str("#   - Optionally filters to informative transcripts via ASE_SUMMARY\n");
    s.push_str("#   - Per-sample, per-transcript binomial exact test of hap0 vs hap1\n");
    s.push_str("#     reads, with null proportion set by that sample's total\n");
    s.push_str("#     hap0 : hap1 library-size ratio (corrects for unequal mapping).\n");
    s.push_str("#   - BH-adjusts p-values across transcripts within each sample\n");
    s.push_str("#   - Writes a long-format ase_results.tsv and a histogram PDF.\n");
    s.push_str("#\n");
    s.push_str("# Requirements: base R >= 4.0. No Bioconductor packages needed.\n");
    s.push_str("\n");

    // Inputs block
    s.push_str(
        "# ---------------- Inputs (edit to re-run with different params) ----------------\n",
    );
    s.push_str(&format!(
        "ASE_COUNTS         <- \"{}\"\n",
        cfg.counts_path.to_string_lossy()
    ));
    s.push_str(&format!(
        "ASE_SUMMARY        <- {}\n",
        r_string_option(
            cfg.summary_path
                .as_ref()
                .map(|p| p.to_string_lossy().into_owned())
        )
    ));
    s.push_str(&format!(
        "HAP_NAMES          <- c({}, {})\n",
        r_string(&cfg.haplotype_names[0]),
        r_string(&cfg.haplotype_names[1])
    ));
    s.push_str(&format!("FDR_CUTOFF         <- {}\n", cfg.fdr_cutoff));
    s.push_str(&format!(
        "MIN_READS_PER_SAMPLE <- {}\n",
        cfg.min_reads_per_sample
    ));
    s.push_str(&format!(
        "FILTER_UNINFORMATIVE <- {}\n",
        if cfg.filter_uninformative {
            "TRUE"
        } else {
            "FALSE"
        }
    ));
    s.push_str("OUT_PREFIX         <- \"ase\"\n");
    s.push_str("\n");

    s.push_str(LOAD_INPUTS);
    s.push_str("\n");
    s.push_str(FILTER_INFORMATIVE);
    s.push_str("\n");
    s.push_str(COLUMN_PARSER);
    s.push_str("\n");
    s.push_str(PER_SAMPLE_NULL);
    s.push_str("\n");
    s.push_str(BINOMIAL_LOOP);
    s.push_str("\n");
    s.push_str(WRITE_RESULTS);
    s.push_str("\n");
    s.push_str(PLOTS);
    s.push_str("\n");

    // Session info
    s.push_str("cat(\"\\n---------------- sessionInfo() ----------------\\n\")\n");
    s.push_str("sessionInfo()\n");

    s
}

fn r_string_option(v: Option<String>) -> String {
    match v {
        Some(s) => r_string(&s),
        None => "NULL".to_string(),
    }
}

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
// Static R snippets
// ─────────────────────────────────────────────────────────────────────────────

const LOAD_INPUTS: &str = r###"# ---------------- Load counts matrix ----------------
counts <- read.table(ASE_COUNTS, header = TRUE, sep = "\t",
                      row.names = 1, check.names = FALSE,
                      comment.char = "")
cat(sprintf("Loaded %d transcripts × %d columns from %s\n",
            nrow(counts), ncol(counts), ASE_COUNTS))

# Columns must be of the form "<sample>.<hap>" — check before going further.
if (ncol(counts) %% 2 != 0) {
  stop(sprintf(
    "ASE counts matrix has an odd number of columns (%d). Expected pairs of <sample>.<hap>.",
    ncol(counts)))
}
"###;

const FILTER_INFORMATIVE: &str = r###"# ---------------- Filter to informative transcripts ----------------
if (FILTER_UNINFORMATIVE && !is.null(ASE_SUMMARY)) {
  summary_df <- read.table(ASE_SUMMARY, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, comment.char = "")
  if (!all(c("transcript", "informative") %in% names(summary_df))) {
    stop("ASE_SUMMARY ", ASE_SUMMARY,
         " is missing required columns 'transcript' and 'informative'")
  }
  informative_ids <- summary_df$transcript[as.logical(summary_df$informative)]
  before <- nrow(counts)
  counts <- counts[rownames(counts) %in% informative_ids, , drop = FALSE]
  cat(sprintf("Filtered to informative transcripts: %d → %d\n",
              before, nrow(counts)))
  if (nrow(counts) == 0) {
    stop("No informative transcripts remain after filtering. ",
         "Use --include-uninformative on the Rust side if you want to test all transcripts.")
  }
}
"###;

const COLUMN_PARSER: &str = r###"# ---------------- Parse <sample>.<hap> columns ----------------
# Columns are named "<sample>.<hap_name>" by the Rust side. We split
# on the LAST period so sample IDs may contain dots without confusing
# the parser.
split_col <- function(col) {
  last_dot <- max(gregexpr("\\.", col)[[1]])
  if (last_dot <= 0) {
    stop(sprintf("Column '%s' does not contain a '.<hap>' suffix", col))
  }
  list(sample = substr(col, 1, last_dot - 1),
       hap    = substr(col, last_dot + 1, nchar(col)))
}
parsed <- lapply(colnames(counts), split_col)
col_samples <- vapply(parsed, function(x) x$sample, character(1))
col_haps    <- vapply(parsed, function(x) x$hap,    character(1))

unique_haps <- unique(col_haps)
if (!setequal(unique_haps, HAP_NAMES)) {
  stop(sprintf(
    "Haplotype labels in %s (%s) don't match --haplotype-names (%s).",
    ASE_COUNTS,
    paste(sort(unique_haps), collapse = ","),
    paste(sort(HAP_NAMES),   collapse = ",")))
}

unique_samples <- unique(col_samples)
for (s in unique_samples) {
  present <- col_haps[col_samples == s]
  if (!setequal(present, HAP_NAMES)) {
    stop(sprintf("Sample '%s' is missing one or both haplotype columns; has: %s",
                 s, paste(sort(present), collapse = ",")))
  }
}
cat(sprintf("Found %d samples × %d haplotypes\n",
            length(unique_samples), length(unique_haps)))
"###;

const PER_SAMPLE_NULL: &str = r###"# ---------------- Per-sample null proportion ----------------
# For each sample, compute hap0_total / (hap0_total + hap1_total) across
# all transcripts. This is the sample's expected hap0 fraction under a
# no-ASE null and corrects for global mapping-rate differences between
# the two personalized transcriptomes.
null_p0 <- numeric(length(unique_samples))
names(null_p0) <- unique_samples
for (s in unique_samples) {
  hap0_col <- which(col_samples == s & col_haps == HAP_NAMES[1])
  hap1_col <- which(col_samples == s & col_haps == HAP_NAMES[2])
  hap0_total <- sum(counts[, hap0_col])
  hap1_total <- sum(counts[, hap1_col])
  if (hap0_total + hap1_total == 0) {
    null_p0[s] <- 0.5
  } else {
    null_p0[s] <- hap0_total / (hap0_total + hap1_total)
  }
}
cat("Per-sample null p(hap0):\n")
print(round(null_p0, 4))
"###;

const BINOMIAL_LOOP: &str = r###"# ---------------- Binomial ASE test ----------------
# For each (transcript, sample): binom.test(hap0, hap0+hap1, p = null_p0[sample]).
# Long-format output: one row per transcript × sample.
results_list <- list()
row_i <- 1L
for (tx in rownames(counts)) {
  for (s in unique_samples) {
    h0 <- counts[tx, which(col_samples == s & col_haps == HAP_NAMES[1])]
    h1 <- counts[tx, which(col_samples == s & col_haps == HAP_NAMES[2])]
    total <- h0 + h1
    if (total < MIN_READS_PER_SAMPLE) {
      results_list[[row_i]] <- data.frame(
        transcript = tx, sample = s,
        hap0_count = h0, hap1_count = h1,
        total = total,
        hap0_frac = NA_real_,
        imbalance = NA_real_,
        null_p    = null_p0[s],
        pvalue    = NA_real_,
        padj      = NA_real_,
        reason    = "below_min_reads",
        stringsAsFactors = FALSE)
    } else {
      pv <- tryCatch(
        binom.test(round(h0), round(total), p = null_p0[s])$p.value,
        error = function(e) NA_real_)
      frac <- h0 / total
      imbal <- frac - null_p0[s]
      results_list[[row_i]] <- data.frame(
        transcript = tx, sample = s,
        hap0_count = h0, hap1_count = h1,
        total = total,
        hap0_frac = frac,
        imbalance = imbal,
        null_p    = null_p0[s],
        pvalue    = pv,
        padj      = NA_real_,
        reason    = "tested",
        stringsAsFactors = FALSE)
    }
    row_i <- row_i + 1L
  }
}
results <- do.call(rbind, results_list)
cat(sprintf("Built %d (transcript × sample) rows, %d below min-reads\n",
            nrow(results), sum(results$reason == "below_min_reads")))

# BH adjustment within each sample (multiple-test correction is per
# sample, not pooled — each sample is its own independent experiment).
for (s in unique_samples) {
  idx <- which(results$sample == s & !is.na(results$pvalue))
  results$padj[idx] <- p.adjust(results$pvalue[idx], method = "BH")
}

results$significant <- !is.na(results$padj) & results$padj < FDR_CUTOFF
cat(sprintf("Significant at FDR<%.3f: %d rows (%d unique transcripts)\n",
            FDR_CUTOFF,
            sum(results$significant),
            length(unique(results$transcript[results$significant]))))
"###;

const WRITE_RESULTS: &str = r###"# ---------------- Write results ----------------
out_tsv <- sprintf("%s_results.tsv", OUT_PREFIX)
write.table(results, out_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("Wrote %s (%d rows)\n", out_tsv, nrow(results)))
"###;

const PLOTS: &str = r###"# ---------------- Allelic-imbalance histogram ----------------
# One panel per sample. X = hap0_frac, vertical line at null_p0.
pdf_path <- sprintf("%s_imbalance.pdf", OUT_PREFIX)
pdf(pdf_path, width = 4 * max(1, length(unique_samples)), height = 3.5)
op <- par(mfrow = c(1, length(unique_samples)), mar = c(4, 4, 2.5, 1))
for (s in unique_samples) {
  vals <- results$hap0_frac[results$sample == s & !is.na(results$hap0_frac)]
  if (length(vals) == 0) {
    plot.new()
    title(main = sprintf("%s  (no tested transcripts)", s))
    next
  }
  hist(vals, breaks = 40,
       main = sprintf("%s  (n=%d)", s, length(vals)),
       xlab = sprintf("hap0 fraction  (null = %.3f)", null_p0[s]),
       col = "grey80", border = "grey40")
  abline(v = null_p0[s], col = "firebrick", lty = 2, lwd = 2)
  abline(v = 0.5,         col = "grey40",  lty = 3)
}
par(op)
dev.off()
cat(sprintf("Wrote %s\n", pdf_path))
"###;

// ─────────────────────────────────────────────────────────────────────────────
// Timestamp helpers — local copy, same pattern as de_template.rs.
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

    fn argv(s: &str) -> Vec<String> {
        s.split('|').map(|p| p.trim().to_string()).collect()
    }

    // ── parse tests ─────────────────────────────────────────────────────────

    #[test]
    fn parse_minimal_counts_only() {
        let a = argv("--counts|ase_counts.tsv");
        let cfg = AseTemplateConfig::parse(&a).unwrap();
        assert_eq!(cfg.counts_path, PathBuf::from("ase_counts.tsv"));
        assert!(cfg.summary_path.is_none());
        assert_eq!(cfg.output_script, PathBuf::from("ase_analysis.R"));
        assert_eq!(
            cfg.haplotype_names,
            ["hap0".to_string(), "hap1".to_string()]
        );
        assert!((cfg.fdr_cutoff - 0.05).abs() < 1e-9);
        assert_eq!(cfg.min_reads_per_sample, 20);
        assert!(cfg.filter_uninformative);
    }

    #[test]
    fn parse_ase_dir_expands_to_counts_and_summary() {
        let a = argv("--ase-dir|ase_out");
        let cfg = AseTemplateConfig::parse(&a).unwrap();
        assert_eq!(cfg.counts_path, PathBuf::from("ase_out/ase_counts.tsv"));
        assert_eq!(
            cfg.summary_path.as_deref(),
            Some(Path::new("ase_out/ase_summary.tsv"))
        );
    }

    #[test]
    fn parse_explicit_flags_beat_ase_dir() {
        let a = argv("--counts|x.tsv|--summary|y.tsv|--ase-dir|ase_out");
        let cfg = AseTemplateConfig::parse(&a).unwrap();
        assert_eq!(cfg.counts_path, PathBuf::from("x.tsv"));
        assert_eq!(cfg.summary_path.as_deref(), Some(Path::new("y.tsv")));
    }

    #[test]
    fn parse_rejects_missing_counts() {
        let a = argv("--fdr|0.1");
        let err = AseTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--counts"));
    }

    #[test]
    fn parse_custom_fdr_and_min_reads() {
        let a = argv("--counts|c.tsv|--fdr|0.01|--min-reads|50");
        let cfg = AseTemplateConfig::parse(&a).unwrap();
        assert!((cfg.fdr_cutoff - 0.01).abs() < 1e-9);
        assert_eq!(cfg.min_reads_per_sample, 50);
    }

    #[test]
    fn parse_custom_haplotype_names() {
        let a = argv("--counts|c.tsv|--haplotype-names|paternal,maternal");
        let cfg = AseTemplateConfig::parse(&a).unwrap();
        assert_eq!(
            cfg.haplotype_names,
            ["paternal".to_string(), "maternal".to_string()]
        );
    }

    #[test]
    fn parse_rejects_duplicate_haplotype_names() {
        let a = argv("--counts|c.tsv|--haplotype-names|hap0,hap0");
        let err = AseTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("must differ"));
    }

    #[test]
    fn parse_rejects_bad_haplotype_names_format() {
        let a = argv("--counts|c.tsv|--haplotype-names|only_one");
        let err = AseTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("NAME1,NAME2"));
    }

    #[test]
    fn parse_include_uninformative_flips_filter() {
        let a = argv("--counts|c.tsv|--include-uninformative");
        let cfg = AseTemplateConfig::parse(&a).unwrap();
        assert!(!cfg.filter_uninformative);
    }

    #[test]
    fn parse_rejects_bad_fdr_value() {
        let a = argv("--counts|c.tsv|--fdr|not_a_number");
        let err = AseTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--fdr"));
    }

    #[test]
    fn parse_rejects_unknown_flag() {
        let a = argv("--counts|c.tsv|--bogus|42");
        let err = AseTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("unknown flag"));
    }

    // ── render tests ────────────────────────────────────────────────────────

    fn sample_cfg() -> AseTemplateConfig {
        AseTemplateConfig {
            counts_path: PathBuf::from("ase_out/ase_counts.tsv"),
            summary_path: Some(PathBuf::from("ase_out/ase_summary.tsv")),
            output_script: PathBuf::from("ase_analysis.R"),
            haplotype_names: ["hap0".to_string(), "hap1".to_string()],
            fdr_cutoff: 0.05,
            min_reads_per_sample: 20,
            filter_uninformative: true,
            raw_cmdline: "myconote-cli ase-template --ase-dir ase_out".to_string(),
        }
    }

    #[test]
    fn render_embeds_input_paths() {
        let s = render_script(&sample_cfg());
        assert!(s.contains("ase_out/ase_counts.tsv"));
        assert!(s.contains("ase_out/ase_summary.tsv"));
    }

    #[test]
    fn render_carries_parameters_as_constants() {
        let s = render_script(&sample_cfg());
        assert!(s.contains("FDR_CUTOFF         <- 0.05"));
        assert!(s.contains("MIN_READS_PER_SAMPLE <- 20"));
        assert!(s.contains("FILTER_UNINFORMATIVE <- TRUE"));
        assert!(s.contains("HAP_NAMES          <- c(\"hap0\", \"hap1\")"));
    }

    #[test]
    fn render_includes_binomial_test_call() {
        let s = render_script(&sample_cfg());
        assert!(s.contains("binom.test("));
        assert!(s.contains("p.adjust("));
        assert!(s.contains("method = \"BH\""));
    }

    #[test]
    fn render_handles_missing_summary_as_null() {
        let mut cfg = sample_cfg();
        cfg.summary_path = None;
        let s = render_script(&cfg);
        assert!(s.contains("ASE_SUMMARY        <- NULL"));
    }

    #[test]
    fn render_handles_include_uninformative() {
        let mut cfg = sample_cfg();
        cfg.filter_uninformative = false;
        let s = render_script(&cfg);
        assert!(s.contains("FILTER_UNINFORMATIVE <- FALSE"));
    }

    #[test]
    fn render_embeds_command_line() {
        let s = render_script(&sample_cfg());
        assert!(s.contains("# Command:   myconote-cli ase-template --ase-dir ase_out"));
    }

    #[test]
    fn render_has_shebang_and_sessioninfo() {
        let s = render_script(&sample_cfg());
        assert!(s.starts_with("#!/usr/bin/env Rscript\n"));
        assert!(s.contains("sessionInfo()"));
    }

    #[test]
    fn r_string_escapes_backslashes_and_quotes() {
        assert_eq!(r_string("foo"), "\"foo\"");
        assert_eq!(r_string("a\\b"), "\"a\\\\b\"");
        assert_eq!(r_string("a\"b"), "\"a\\\"b\"");
    }

    #[test]
    fn r_string_option_none_is_NULL() {
        assert_eq!(r_string_option(None), "NULL");
        assert_eq!(r_string_option(Some("x".to_string())), "\"x\"");
    }

    #[test]
    fn format_rfc3339_epoch_and_known() {
        assert_eq!(format_rfc3339(0), "1970-01-01T00:00:00Z");
        // Value computed from the same algorithm; guards against
        // accidental drift in the Zeller-style day-of-year math.
        assert_eq!(format_rfc3339(1_777_217_400), "2026-04-26T15:30:00Z");
    }
}
