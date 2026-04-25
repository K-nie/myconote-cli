//! GO enrichment R-script generator.
//!
//! `myconote-cli go-template` takes a `de-template` results TSV +
//! a `myconote-cli annotate` `annotations.tsv` (carrying a `go_terms`
//! column) and emits a self-contained R script. The user runs that
//! script with `Rscript`; it builds the gene → GO mapping, marks
//! significant genes from the DE results, and runs `topGO`'s classic
//! Fisher's exact test across BP / MF / CC ontologies.
//!
//! Same Option 1D pattern as `de-template` and `ase-template`. We do
//! not wrap R as a runtime dep — the tool writes R, the user runs R.
//! Prerequisite alerting in three layers: CLI help block, runtime
//! stderr if `Rscript` is missing, and a `requireNamespace("topGO")`
//! guard at the top of the emitted script with a clear
//! `BiocManager::install("topGO")` hint.

use crate::utils::error::{MycoNoteError, Result};
use std::fs;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

/// Parsed command-line configuration for a `go-template` invocation.
#[derive(Debug, Clone)]
pub struct GoTemplateConfig {
    /// DE results TSV produced by `de-template` (e.g.
    /// `de_condition_treated_vs_control.tsv`). Columns: transcript,
    /// baseMean, log2FoldChange, lfcSE, pvalue, padj.
    pub de_results: PathBuf,
    /// Functional-annotation TSV produced by `annotate`. Columns
    /// include `locus_tag` and `go_terms` (pipe-separated by default).
    pub annotations: PathBuf,
    /// Output R script path.
    pub output_script: PathBuf,
    /// FDR threshold on `padj` for marking a gene "significant" in the
    /// classic Fisher's test against the universe of tested genes.
    pub fdr_cutoff: f64,
    /// Which ontology branch(es) to test. `OntologyChoice::All` runs
    /// BP + MF + CC and writes one TSV per ontology; the others run
    /// just the named branch.
    pub ontology: OntologyChoice,
    /// Top-N enriched terms shown in the dot plot.
    pub top_terms: u32,
    /// Column name in the DE results TSV that carries the gene/transcript
    /// identifier. Default `transcript` (matches `de-template` output).
    pub de_id_col: String,
    /// Column name in the annotations TSV that carries the gene
    /// identifier. Default `locus_tag` (matches `annotate` output).
    pub ann_id_col: String,
    /// Column name in the annotations TSV that carries the GO terms.
    /// Default `go_terms`.
    pub go_col: String,
    /// Separator used inside the GO column to delimit multiple terms.
    /// Default `|` (matches `annotate` output).
    pub go_separator: String,
    /// Verbatim argv for embedding in the generated script header.
    pub raw_cmdline: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OntologyChoice {
    Bp,
    Mf,
    Cc,
    All,
}

impl OntologyChoice {
    fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_uppercase().as_str() {
            "BP" => Some(Self::Bp),
            "MF" => Some(Self::Mf),
            "CC" => Some(Self::Cc),
            "ALL" => Some(Self::All),
            _ => None,
        }
    }

    fn as_r_vector(&self) -> &'static str {
        match self {
            Self::Bp => "c(\"BP\")",
            Self::Mf => "c(\"MF\")",
            Self::Cc => "c(\"CC\")",
            Self::All => "c(\"BP\", \"MF\", \"CC\")",
        }
    }
}

impl GoTemplateConfig {
    pub fn parse(args: &[String]) -> Result<Self> {
        let mut de_results: Option<PathBuf> = None;
        let mut annotations: Option<PathBuf> = None;
        let mut output_script: PathBuf = PathBuf::from("go_enrichment.R");
        let mut fdr_cutoff: f64 = 0.05;
        let mut ontology = OntologyChoice::All;
        let mut top_terms: u32 = 30;
        let mut de_id_col = "transcript".to_string();
        let mut ann_id_col = "locus_tag".to_string();
        let mut go_col = "go_terms".to_string();
        let mut go_separator = "|".to_string();

        let mut i = 0usize;
        while i < args.len() {
            match args[i].as_str() {
                "--de-results" if i + 1 < args.len() => {
                    de_results = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--annotations" if i + 1 < args.len() => {
                    annotations = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--output" | "-o" if i + 1 < args.len() => {
                    output_script = PathBuf::from(&args[i + 1]);
                    i += 2;
                }
                "--fdr" if i + 1 < args.len() => {
                    fdr_cutoff = args[i + 1].parse().map_err(|_| {
                        MycoNoteError::QuantSheet(format!(
                            "go-template: --fdr must be a number, got '{}'",
                            args[i + 1]
                        ))
                    })?;
                    i += 2;
                }
                "--ontology" if i + 1 < args.len() => {
                    ontology = OntologyChoice::parse(&args[i + 1]).ok_or_else(|| {
                        MycoNoteError::QuantSheet(format!(
                            "go-template: --ontology must be BP|MF|CC|all, got '{}'",
                            args[i + 1]
                        ))
                    })?;
                    i += 2;
                }
                "--top" if i + 1 < args.len() => {
                    top_terms = args[i + 1].parse().map_err(|_| {
                        MycoNoteError::QuantSheet(format!(
                            "go-template: --top must be a positive integer, got '{}'",
                            args[i + 1]
                        ))
                    })?;
                    if top_terms == 0 {
                        return Err(MycoNoteError::QuantSheet(
                            "go-template: --top must be > 0".to_string(),
                        ));
                    }
                    i += 2;
                }
                "--de-id-col" if i + 1 < args.len() => {
                    de_id_col = args[i + 1].clone();
                    i += 2;
                }
                "--ann-id-col" if i + 1 < args.len() => {
                    ann_id_col = args[i + 1].clone();
                    i += 2;
                }
                "--go-col" if i + 1 < args.len() => {
                    go_col = args[i + 1].clone();
                    i += 2;
                }
                "--go-separator" if i + 1 < args.len() => {
                    if args[i + 1].is_empty() {
                        return Err(MycoNoteError::QuantSheet(
                            "go-template: --go-separator must be a non-empty string".to_string(),
                        ));
                    }
                    go_separator = args[i + 1].clone();
                    i += 2;
                }
                other => {
                    return Err(MycoNoteError::QuantSheet(format!(
                        "go-template: unknown flag '{other}'"
                    )));
                }
            }
        }

        let de_results = de_results.ok_or_else(|| {
            MycoNoteError::QuantSheet(
                "go-template: --de-results <de_*.tsv> is required".to_string(),
            )
        })?;
        let annotations = annotations.ok_or_else(|| {
            MycoNoteError::QuantSheet(
                "go-template: --annotations <annotations.tsv> is required".to_string(),
            )
        })?;

        let raw_cmdline = std::iter::once("myconote-cli go-template".to_string())
            .chain(args.iter().cloned())
            .collect::<Vec<_>>()
            .join(" ");

        Ok(Self {
            de_results,
            annotations,
            output_script,
            fdr_cutoff,
            ontology,
            top_terms,
            de_id_col,
            ann_id_col,
            go_col,
            go_separator,
            raw_cmdline,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Top-level dispatcher
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_go_template(args: &[String]) -> Result<()> {
    let cfg = GoTemplateConfig::parse(args)?;

    warn_about_r_prereqs();

    require_file(&cfg.de_results, "--de-results")?;
    require_file(&cfg.annotations, "--annotations")?;

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
    eprintln!("  BiocManager::install(c(\"topGO\"))");
    Ok(())
}

fn warn_about_r_prereqs() {
    if which::which("Rscript").is_err() {
        eprintln!(
            "  ⚠ `Rscript` not found on PATH. The generated script needs R and \
             the Bioconductor package topGO. The script will be written anyway \
             so you can run it on another host."
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

/// Render a full R script from a validated `GoTemplateConfig`. Pure —
/// returns a String, no I/O.
pub fn render_script(cfg: &GoTemplateConfig) -> String {
    let mut s = String::with_capacity(4096);

    let version = env!("CARGO_PKG_VERSION");
    let timestamp = timestamp_rfc3339();

    // Header
    s.push_str("#!/usr/bin/env Rscript\n");
    s.push_str(&format!(
        "# Auto-generated by myconote-cli go-template v{version}\n"
    ));
    s.push_str(&format!("# Generated: {timestamp}\n"));
    s.push_str(&format!("# Command:   {}\n", cfg.raw_cmdline));
    s.push_str("#\n");
    s.push_str("# Usage: Rscript ");
    s.push_str(
        &cfg.output_script
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "go_enrichment.R".to_string()),
    );
    s.push('\n');
    s.push_str("#\n");
    s.push_str("# What this does:\n");
    s.push_str("#   - Reads DE_RESULTS and marks each gene as significant (padj < FDR_CUTOFF)\n");
    s.push_str("#   - Reads ANNOTATIONS and builds a gene -> GO mapping from GO_COL\n");
    s.push_str("#   - Joins on (DE_ID_COL == ANN_ID_COL) to define the universe of tested genes\n");
    s.push_str("#   - Runs topGO classic Fisher's exact test for each ontology in ONTOLOGIES\n");
    s.push_str("#   - Writes one TSV per ontology + a combined dot plot PDF\n");
    s.push_str("#\n");
    s.push_str("# REQUIRES (install once, before running this script):\n");
    s.push_str("#   if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")\n");
    s.push_str("#   BiocManager::install(c(\"topGO\"))\n");
    s.push_str("\n");

    // Inputs block
    s.push_str(
        "# ---------------- Inputs (edit to re-run with different params) ----------------\n",
    );
    s.push_str(&format!(
        "DE_RESULTS    <- \"{}\"\n",
        cfg.de_results.to_string_lossy()
    ));
    s.push_str(&format!(
        "ANNOTATIONS   <- \"{}\"\n",
        cfg.annotations.to_string_lossy()
    ));
    s.push_str(&format!("FDR_CUTOFF    <- {}\n", cfg.fdr_cutoff));
    s.push_str(&format!(
        "ONTOLOGIES    <- {}\n",
        cfg.ontology.as_r_vector()
    ));
    s.push_str(&format!("TOP_TERMS     <- {}\n", cfg.top_terms));
    s.push_str(&format!("DE_ID_COL     <- {}\n", r_string(&cfg.de_id_col)));
    s.push_str(&format!("ANN_ID_COL    <- {}\n", r_string(&cfg.ann_id_col)));
    s.push_str(&format!("GO_COL        <- {}\n", r_string(&cfg.go_col)));
    s.push_str(&format!(
        "GO_SEPARATOR  <- {}\n",
        r_string(&cfg.go_separator)
    ));
    s.push_str("OUT_PREFIX    <- \"go\"\n");
    s.push_str("\n");

    s.push_str(PACKAGE_GUARD);
    s.push_str("\n");
    s.push_str(LOAD_DE);
    s.push_str("\n");
    s.push_str(LOAD_ANNOTATIONS);
    s.push_str("\n");
    s.push_str(BUILD_GENE2GO);
    s.push_str("\n");
    s.push_str(BUILD_FACTOR);
    s.push_str("\n");
    s.push_str(RUN_TOPGO);
    s.push_str("\n");
    s.push_str(WRITE_AND_PLOT);
    s.push_str("\n");

    // Session info
    s.push_str("cat(\"\\n---------------- sessionInfo() ----------------\\n\")\n");
    s.push_str("sessionInfo()\n");

    s
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
// Static R snippets
// ─────────────────────────────────────────────────────────────────────────────

const PACKAGE_GUARD: &str = r###"# ---------------- Package guard ----------------
if (!requireNamespace("topGO", quietly = TRUE)) {
  stop(paste0(
    "\n\nMissing R package: topGO.\n",
    "Install (one-time) before running this script:\n",
    "  if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")\n",
    "  BiocManager::install(c(\"topGO\"))\n"
  ))
}
suppressPackageStartupMessages({
  library(topGO)
})
"###;

const LOAD_DE: &str = r###"# ---------------- Load DE results ----------------
de <- read.table(DE_RESULTS, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, comment.char = "")
if (!(DE_ID_COL %in% names(de))) {
  stop(sprintf("DE id column '%s' not found in %s. Available: %s",
               DE_ID_COL, DE_RESULTS, paste(names(de), collapse = ", ")))
}
if (!("padj" %in% names(de))) {
  stop(sprintf("DE results %s is missing the 'padj' column", DE_RESULTS))
}
de_ids   <- as.character(de[[DE_ID_COL]])
de_padj  <- suppressWarnings(as.numeric(de$padj))
cat(sprintf("Loaded %d rows from %s\n", nrow(de), DE_RESULTS))
"###;

const LOAD_ANNOTATIONS: &str = r###"# ---------------- Load annotations TSV ----------------
ann <- read.table(ANNOTATIONS, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, comment.char = "",
                  quote = "")
if (!(ANN_ID_COL %in% names(ann))) {
  stop(sprintf("Annotation id column '%s' not found in %s. Available: %s",
               ANN_ID_COL, ANNOTATIONS, paste(names(ann), collapse = ", ")))
}
if (!(GO_COL %in% names(ann))) {
  stop(sprintf("GO terms column '%s' not found in %s. Available: %s",
               GO_COL, ANNOTATIONS, paste(names(ann), collapse = ", ")))
}
cat(sprintf("Loaded %d annotation rows from %s\n", nrow(ann), ANNOTATIONS))
"###;

const BUILD_GENE2GO: &str = r###"# ---------------- Build gene -> GO mapping ----------------
# topGO wants a named list: names = gene IDs, values = char vec of GO IDs.
ann_ids <- as.character(ann[[ANN_ID_COL]])
ann_go  <- as.character(ann[[GO_COL]])

split_terms <- function(x) {
  if (is.na(x) || x == "" || x == "-") return(character(0))
  parts <- strsplit(x, GO_SEPARATOR, fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nchar(parts) > 0]
  # Keep only well-formed GO IDs ("GO:" + 7 digits). Dropping malformed
  # rows here keeps the topGO call from failing far downstream.
  parts <- parts[grepl("^GO:[0-9]{7}$", parts)]
  unique(parts)
}

gene2go <- lapply(ann_go, split_terms)
names(gene2go) <- ann_ids
# Drop genes with no GO terms at all from the mapping itself.
gene2go <- gene2go[vapply(gene2go, length, integer(1)) > 0]
cat(sprintf("gene -> GO map: %d genes with >=1 well-formed GO term\n",
            length(gene2go)))
if (length(gene2go) == 0) {
  stop(sprintf("No GO terms recovered from %s. Check --go-col and --go-separator.",
               ANNOTATIONS))
}
"###;

const BUILD_FACTOR: &str = r###"# ---------------- Build the topGO foreground / background factor ----------------
# Universe = genes that appear in BOTH the DE table and the gene -> GO map.
# Foreground = subset of universe with padj < FDR_CUTOFF.
universe <- intersect(de_ids, names(gene2go))
if (length(universe) == 0) {
  stop("No overlap between DE_ID_COL values and ANN_ID_COL values. ",
       "Check that DE results and annotations share an identifier ",
       "(set --de-id-col / --ann-id-col).")
}
de_padj_named <- setNames(de_padj, de_ids)
sig_universe  <- universe[!is.na(de_padj_named[universe]) &
                          de_padj_named[universe] < FDR_CUTOFF]
cat(sprintf("Universe: %d genes; significant at FDR<%.3f: %d\n",
            length(universe), FDR_CUTOFF, length(sig_universe)))
if (length(sig_universe) == 0) {
  warning("No genes pass FDR_CUTOFF — topGO will run but no enrichment ",
          "is meaningful. Lower --fdr or check the DE results.")
}

all_genes <- factor(as.integer(universe %in% sig_universe), levels = c(0, 1))
names(all_genes) <- universe
"###;

const RUN_TOPGO: &str = r###"# ---------------- Run topGO classic Fisher's test ----------------
# We run the "classic" algorithm — straight Fisher's exact test on the
# 2x2 (sig vs not, in-term vs not-in-term) table. topGO has fancier
# algorithms (elim, weight01) that decorrelate parent/child terms; the
# classic test is the most interpretable starting point. Users wanting
# elim/weight01 can edit the algorithm = "classic" below.

results_per_ont <- list()
for (ont in ONTOLOGIES) {
  cat(sprintf("\n=== %s ===\n", ont))

  godata <- new("topGOdata",
                ontology    = ont,
                allGenes    = all_genes,
                annot       = annFUN.gene2GO,
                gene2GO     = gene2go,
                nodeSize    = 5)

  test_stat <- new("classicCount",
                   testStatistic = GOFisherTest,
                   name = "Fisher")
  result <- getSigGroups(godata, test_stat)

  pvals <- score(result)
  res_df <- GenTable(godata,
                     classic = result,
                     orderBy  = "classic",
                     topNodes = length(pvals))
  # GenTable returns p-values as character; ensure numeric for sorting
  # and BH adjustment.
  res_df$classic_p <- suppressWarnings(as.numeric(res_df$classic))
  # GenTable serialises very small p-values as "<1e-30" — preserve that
  # in a separate column for the user.
  res_df$classic_p_raw <- res_df$classic
  # BH-adjust across all tested terms in this ontology.
  res_df$classic_padj <- p.adjust(res_df$classic_p, method = "BH")
  res_df <- res_df[order(res_df$classic_p, na.last = TRUE), ]
  res_df$ontology <- ont
  results_per_ont[[ont]] <- res_df

  cat(sprintf("  %d GO terms tested in %s; %d with classic p < 0.05\n",
              nrow(res_df), ont,
              sum(!is.na(res_df$classic_p) & res_df$classic_p < 0.05)))
}
"###;

const WRITE_AND_PLOT: &str = r###"# ---------------- Write per-ontology TSVs ----------------
for (ont in names(results_per_ont)) {
  res_df <- results_per_ont[[ont]]
  out_tsv <- sprintf("%s_%s_enrichment.tsv", OUT_PREFIX, ont)
  write.table(res_df, out_tsv, sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
  cat(sprintf("Wrote %s (%d rows)\n", out_tsv, nrow(res_df)))
}

# ---------------- Combined dot plot ----------------
# One panel per ontology. Top TOP_TERMS by classic_p, x = -log10(p),
# size = significant gene count, color = classic_padj.
pdf_path <- sprintf("%s_dotplot.pdf", OUT_PREFIX)
n_panels <- length(results_per_ont)
pdf(pdf_path, width = 6 * max(1, n_panels), height = 8)
op <- par(mfrow = c(1, n_panels), mar = c(4, 18, 3, 1), las = 1)
for (ont in names(results_per_ont)) {
  res_df <- results_per_ont[[ont]]
  if (nrow(res_df) == 0) {
    plot.new(); title(main = sprintf("%s  (no terms tested)", ont))
    next
  }
  top <- head(res_df, TOP_TERMS)
  top <- top[!is.na(top$classic_p), ]
  if (nrow(top) == 0) {
    plot.new(); title(main = sprintf("%s  (no significant terms)", ont))
    next
  }
  top <- top[order(top$classic_p, decreasing = TRUE), ]  # bottom = most significant
  labels <- sprintf("%s  %s", top$GO.ID,
                    substr(top$Term, 1, 50))
  x <- -log10(top$classic_p)
  cex <- 0.6 + 1.2 * (top$Significant / max(top$Significant, na.rm = TRUE))
  col_val <- ifelse(top$classic_padj < 0.05, "firebrick", "grey50")

  plot(x, seq_along(x), pch = 19, cex = cex, col = col_val,
       yaxt = "n", xlab = expression(-log[10](italic(p))),
       ylab = "",
       main = sprintf("%s  (top %d by classic Fisher)", ont, nrow(top)))
  axis(2, at = seq_along(x), labels = labels, cex.axis = 0.7)
  abline(v = -log10(0.05), lty = 2, col = "grey40")
}
par(op)
dev.off()
cat(sprintf("Wrote %s\n", pdf_path))
"###;

// ─────────────────────────────────────────────────────────────────────────────
// Timestamp helpers — local copy, same pattern as de_template.rs / ase_template.rs.
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
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(s: &str) -> Vec<String> {
        s.split('|').map(|p| p.trim().to_string()).collect()
    }

    // ── parse tests ─────────────────────────────────────────────────────────

    #[test]
    fn parse_minimal() {
        let a = argv("--de-results|de.tsv|--annotations|ann.tsv");
        let cfg = GoTemplateConfig::parse(&a).unwrap();
        assert_eq!(cfg.de_results, PathBuf::from("de.tsv"));
        assert_eq!(cfg.annotations, PathBuf::from("ann.tsv"));
        assert_eq!(cfg.output_script, PathBuf::from("go_enrichment.R"));
        assert!((cfg.fdr_cutoff - 0.05).abs() < 1e-9);
        assert_eq!(cfg.ontology, OntologyChoice::All);
        assert_eq!(cfg.top_terms, 30);
        assert_eq!(cfg.de_id_col, "transcript");
        assert_eq!(cfg.ann_id_col, "locus_tag");
        assert_eq!(cfg.go_col, "go_terms");
        assert_eq!(cfg.go_separator, "|");
    }

    #[test]
    fn parse_rejects_missing_de_results() {
        let a = argv("--annotations|ann.tsv");
        let err = GoTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--de-results"));
    }

    #[test]
    fn parse_rejects_missing_annotations() {
        let a = argv("--de-results|de.tsv");
        let err = GoTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--annotations"));
    }

    #[test]
    fn parse_custom_thresholds_and_columns() {
        let a = argv(
            "--de-results|de.tsv|--annotations|ann.tsv|--fdr|0.01|--top|50|\
             --de-id-col|gene_id|--ann-id-col|gene_id|--go-col|GO|--go-separator|;",
        );
        let cfg = GoTemplateConfig::parse(&a).unwrap();
        assert!((cfg.fdr_cutoff - 0.01).abs() < 1e-9);
        assert_eq!(cfg.top_terms, 50);
        assert_eq!(cfg.de_id_col, "gene_id");
        assert_eq!(cfg.ann_id_col, "gene_id");
        assert_eq!(cfg.go_col, "GO");
        assert_eq!(cfg.go_separator, ";");
    }

    #[test]
    fn parse_ontology_choices() {
        for (s, expected) in [
            ("BP", OntologyChoice::Bp),
            ("mf", OntologyChoice::Mf),
            ("CC", OntologyChoice::Cc),
            ("all", OntologyChoice::All),
        ] {
            let a = argv(&format!(
                "--de-results|de.tsv|--annotations|ann.tsv|--ontology|{s}"
            ));
            let cfg = GoTemplateConfig::parse(&a).unwrap();
            assert_eq!(cfg.ontology, expected);
        }
    }

    #[test]
    fn parse_rejects_bad_ontology() {
        let a = argv("--de-results|de.tsv|--annotations|ann.tsv|--ontology|XX");
        let err = GoTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--ontology"));
    }

    #[test]
    fn parse_rejects_zero_top() {
        let a = argv("--de-results|de.tsv|--annotations|ann.tsv|--top|0");
        let err = GoTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--top"));
    }

    #[test]
    fn parse_rejects_bad_fdr_value() {
        let a = argv("--de-results|de.tsv|--annotations|ann.tsv|--fdr|nope");
        let err = GoTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--fdr"));
    }

    #[test]
    fn parse_rejects_unknown_flag() {
        let a = argv("--de-results|de.tsv|--annotations|ann.tsv|--bogus|x");
        let err = GoTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("unknown flag"));
    }

    #[test]
    fn parse_rejects_empty_separator() {
        let a = argv("--de-results|de.tsv|--annotations|ann.tsv|--go-separator|");
        let err = GoTemplateConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--go-separator"));
    }

    // ── render tests ────────────────────────────────────────────────────────

    fn sample_cfg() -> GoTemplateConfig {
        GoTemplateConfig {
            de_results: PathBuf::from("/path/de_condition_treated_vs_control.tsv"),
            annotations: PathBuf::from("/path/annotations.tsv"),
            output_script: PathBuf::from("go_enrichment.R"),
            fdr_cutoff: 0.05,
            ontology: OntologyChoice::All,
            top_terms: 30,
            de_id_col: "transcript".to_string(),
            ann_id_col: "locus_tag".to_string(),
            go_col: "go_terms".to_string(),
            go_separator: "|".to_string(),
            raw_cmdline:
                "myconote-cli go-template --de-results de.tsv --annotations annotations.tsv"
                    .to_string(),
        }
    }

    #[test]
    fn render_embeds_input_paths() {
        let s = render_script(&sample_cfg());
        assert!(s.contains("/path/de_condition_treated_vs_control.tsv"));
        assert!(s.contains("/path/annotations.tsv"));
    }

    #[test]
    fn render_carries_thresholds_and_columns() {
        let s = render_script(&sample_cfg());
        assert!(s.contains("FDR_CUTOFF    <- 0.05"));
        assert!(s.contains("TOP_TERMS     <- 30"));
        assert!(s.contains("DE_ID_COL     <- \"transcript\""));
        assert!(s.contains("ANN_ID_COL    <- \"locus_tag\""));
        assert!(s.contains("GO_COL        <- \"go_terms\""));
        assert!(s.contains("GO_SEPARATOR  <- \"|\""));
    }

    #[test]
    fn render_includes_topGO_calls() {
        let s = render_script(&sample_cfg());
        assert!(s.contains("requireNamespace(\"topGO\""));
        assert!(s.contains("BiocManager::install(c(\"topGO\"))"));
        assert!(s.contains("library(topGO)"));
        assert!(s.contains("new(\"topGOdata\""));
        assert!(s.contains("annFUN.gene2GO"));
        assert!(s.contains("GOFisherTest"));
        assert!(s.contains("getSigGroups"));
        assert!(s.contains("p.adjust(res_df$classic_p, method = \"BH\")"));
    }

    #[test]
    fn render_includes_ontology_vector() {
        let mut cfg = sample_cfg();
        cfg.ontology = OntologyChoice::All;
        let s = render_script(&cfg);
        assert!(s.contains("ONTOLOGIES    <- c(\"BP\", \"MF\", \"CC\")"));

        cfg.ontology = OntologyChoice::Bp;
        let s = render_script(&cfg);
        assert!(s.contains("ONTOLOGIES    <- c(\"BP\")"));
    }

    #[test]
    fn render_has_shebang_and_sessioninfo() {
        let s = render_script(&sample_cfg());
        assert!(s.starts_with("#!/usr/bin/env Rscript\n"));
        assert!(s.contains("sessionInfo()"));
    }

    #[test]
    fn render_embeds_command_line() {
        let s = render_script(&sample_cfg());
        assert!(s.contains("# Command:   myconote-cli go-template"));
    }

    #[test]
    fn r_string_escapes_backslashes_and_quotes() {
        assert_eq!(r_string("foo"), "\"foo\"");
        assert_eq!(r_string("a\\b"), "\"a\\\\b\"");
        assert_eq!(r_string("a\"b"), "\"a\\\"b\"");
    }

    // ── rfc3339 ─────────────────────────────────────────────────────────────

    #[test]
    fn format_rfc3339_epoch_and_known() {
        assert_eq!(format_rfc3339(0), "1970-01-01T00:00:00Z");
        // Same anchor used in de_template / ase_template tests for
        // determinism — mid-day on 2026-04-24.
        assert_eq!(format_rfc3339(1777032000), "2026-04-24T12:00:00Z");
    }
}
