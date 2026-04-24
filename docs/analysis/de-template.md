# myconote de-template

Generate a ready-to-run R script that performs **DESeq2** differential-expression analysis on `quant` output. The tool writes the script; you run it with `Rscript`.

This keeps R out of MycoNote-CLI's runtime dependencies while still giving users ~90 % of what a wrapped DE command would provide — you don't have to write tximport + DESeq2 boilerplate yourself, and the script comes configured with apeglm LFC shrinkage, prefiltering, and MA + volcano plots.

## Prerequisites — install BEFORE running the emitted script

You need **R** and these **Bioconductor** packages:

- `tximport`
- `DESeq2`
- `apeglm` (optional but recommended — enables LFC shrinkage)

One-time install from an R session:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("tximport", "DESeq2", "apeglm"))
```

If these aren't installed, the emitted script will fail fast with a clear error message pointing back at this install line. `myconote-cli` itself does not try to install R packages for you.

## Quickstart

After `myconote-cli quant` finishes, you have a `quant_out/` directory with `sample_sheet.tsv` and `salmon/<sample>/quant.sf` files. Generate the DE script:

```bash
myconote-cli de-template \
    --quant-dir quant_out \
    --design '~ condition' \
    --contrast 'condition,treated,control' \
    -o analysis.R
```

Then run it:

```bash
Rscript analysis.R
```

You get three output files per contrast:

- `de_condition_treated_vs_control.tsv` — results table, sorted by `padj`
- `de_condition_treated_vs_control_MA.png` — MA plot
- `de_condition_treated_vs_control_volcano.png` — volcano plot

Plus DESeq2's `summary()` output and an R `sessionInfo()` for reproducibility, both to stdout.

## Inputs

| flag | meaning |
|---|---|
| `--quant-dir <dir>` | `quant_out/` directory. Uses tximport against `quant.sf` files (preferred — retains length offsets). |
| `--counts <counts.tsv>` | Fallback: a wide counts matrix (transcript × sample). Use when you don't have per-sample `quant.sf`. |
| `--samples <sheet.tsv>` | Sample sheet. Defaults to `<quant-dir>/sample_sheet.tsv`; required in `--counts` mode. |
| `--design '<R formula>'` | DESeq2 design, e.g. `~ condition` or `~ batch + condition`. Required. |
| `--contrast '<factor,level1,level2>'` | DE contrast to extract. Repeat the flag for multiple contrasts. |
| `--output <file.R>` / `-o` | Output path (default `de_analysis.R`). |
| `--fdr <n>` | Significance threshold (default `0.05`). |
| `--lfc <n>` | `|LFC|` threshold used to highlight points in the volcano plot (default `1.0`). |

One of `--quant-dir` or `--counts` is required — not both. `--design` and at least one `--contrast` are always required.

## Validation the tool does for you before writing the script

1. **Input exists** — `--quant-dir` or `--counts` file is reachable.
2. **Sample sheet parses** — same TSV schema as `quant` (reused).
3. **Factor column exists in the sheet** — `condition`, `batch`, or any custom column. Rejects typos before DESeq2 fails obscurely.
4. **Both contrast levels appear as values in the factor column** — rejects misspelled level names and lists what values *are* present to help the user fix it.
5. **Design formula sanitization** — rejects obvious injection patterns (`;`, backticks, `system(...)`, `eval(...)`, newlines). Not a security boundary; catches copy-paste mistakes.

## What's in the emitted R script

Hand-maintained template. Structure:

```r
# Header with generation metadata
# Inputs block (edit here to re-run with different thresholds)
# Package guard — stops with install hint if DESeq2 missing
# Sample sheet loading
# DESeqDataSet construction (tximport branch OR matrix branch)
# Prefilter low-count rows
# dds <- DESeq(dds)
# For each contrast:
#   - results()
#   - apeglm lfcShrink() with fallback to unshrunk LFCs
#   - Write de_<contrast>.tsv sorted by padj
#   - MA plot
#   - Volcano plot (dashed lines at FDR + LFC thresholds)
# sessionInfo() for reproducibility
```

Every substitution is a named R variable at the top (`QUANT_DIR`, `DESIGN`, `CONTRASTS`, `FDR_CUTOFF`, `PREFILTER_MIN_COUNT`, …). Users can edit those to re-run without regenerating the script.

## Design formula examples

The tool passes your design string through to R verbatim (after sanitization). Whatever DESeq2 accepts is accepted here.

```bash
# Single factor
--design '~ condition'

# Two-factor with batch adjustment
--design '~ batch + condition'

# Interaction term
--design '~ batch + condition + batch:condition'

# Time-course (continuous)
--design '~ timepoint'
```

If your design needs random effects or complex grouping, use `limma` or `dream` directly — those aren't in scope for this template. Edit the emitted `analysis.R` if you want to swap the `DESeq()` step for something else.

## Multiple contrasts in one call

Repeat `--contrast`:

```bash
myconote-cli de-template \
    --quant-dir quant_out \
    --design '~ condition' \
    --contrast 'condition,24h,0h' \
    --contrast 'condition,48h,0h' \
    --contrast 'condition,48h,24h'
```

The script runs `results()` once per contrast and writes one TSV + MA + volcano per tuple.

## Reading the outputs

### `de_<factor>_<l1>_vs_<l2>.tsv`

Columns: `transcript`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`. Rows sorted by `padj` ascending (`NA`s last).

- `padj < 0.05` = significant at the 5% FDR.
- `log2FoldChange`: **post-apeglm shrinkage** unless apeglm failed or isn't installed. Don't compare shrunk LFCs directly to raw LFCs from other pipelines.
- `NA` in `padj`: DESeq2's independent filtering removed the row. Expected, not a bug.

### MA plot

`baseMean` on x, `log2FoldChange` on y. Good runs: symmetric cloud around y = 0 for low-count rows; high-count rows should cluster toward y = 0 unless there's real biology. Systematic drift (all points offset from zero) flags a normalization problem.

### Volcano plot

`log2FoldChange` on x, `-log10(pvalue)` on y. Dashed lines at the `--fdr` and `±--lfc` thresholds. Points above `-log10(fdr)` AND outside `±lfc` are coloured (red up, blue down).

Tall & narrow volcano = few DEGs / underpowered. Wide & scattered = many DEGs / strong treatment effect.

## What `de-template` does NOT do

- **Doesn't run DE.** The tool writes the script. You run `Rscript analysis.R`.
- **Doesn't wrap edgeR / limma.** DESeq2 only. If there's demand, a `--tool edgeR` variant could land in a follow-up.
- **Doesn't install R packages.** That's your R session's job.
- **Doesn't do pathway / GO enrichment.** Separate concern; use `clusterProfiler`, `topGO`, `enrichplot` downstream.

## Troubleshooting

- **"Missing R package(s): DESeq2"** — you haven't run the `BiocManager::install` line from the top of this page. Do that first.
- **"Coefficient '…' not in resultsNames(dds)"** — apeglm needs the coefficient name that DESeq2 derived, which depends on your design and the reference level. DESeq2 chooses the alphabetically-first level as the reference by default. Edit the generated script to set `dds$condition <- relevel(dds$condition, ref = "control")` before `DESeq(dds)` if you want a specific reference.
- **"factor 'X' is not a column in ..."** — the contrast factor doesn't match any column in the sample sheet. Check spelling and rerun.
- **Volcano plot looks empty / tiny** — your FDR threshold may be too strict for the effect sizes present. Try relaxing `--fdr 0.1`, or lower the `--lfc` threshold, then regenerate.

## Related

- `myconote-cli quant` — produces the input `quant_out/` directory.
- `myconote-cli convert --to cds` — prerequisite to `quant`; generates the transcriptome FASTA.
- DESeq2 tutorial: <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>
- tximport tutorial: <https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html>
