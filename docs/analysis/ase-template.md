# myconote ase-template

Generate a ready-to-run R script that performs a **binomial exact test** for allele-specific expression on `myconote-cli ase` output. The tool writes the script; you run it with `Rscript`.

Companion to `ase` — same design philosophy as `de-template` (write R, user runs R, no R runtime inside MycoNote-CLI). The statistics are simple enough that we don't need Bioconductor: base R's `binom.test()` + `p.adjust()` covers it.

## Prerequisites — install before running the emitted script

Base **R ≥ 4.0**. Nothing else.

No `DESeq2`, no `tximport`, no `BiocManager::install`. The emitted script uses only `read.table`, `binom.test`, `p.adjust`, and `hist`.

Verify you have R:

```bash
Rscript --version
```

If `Rscript` isn't on your PATH, the tool warns but still writes the script — you can run it on another host.

## Quickstart

Run it on the output of `myconote-cli ase`:

```bash
myconote-cli ase-template --ase-dir ase_out -o ase_analysis.R
Rscript ase_analysis.R
```

`--ase-dir` is shorthand for `--counts ase_out/ase_counts.tsv --summary ase_out/ase_summary.tsv`.

You get two output files:

- `ase_results.tsv` — long-format table: one row per `transcript × sample` with `pvalue`, BH-adjusted `padj`, `hap0_frac`, `imbalance`, and a `reason` field telling you why a row was skipped (e.g. below `--min-reads`).
- `ase_imbalance.pdf` — per-sample histogram of hap0-fraction, with a dashed red line at that sample's null proportion.

## Inputs

| flag | meaning |
|---|---|
| `--ase-dir <dir>` | Shortcut: fills `--counts` and `--summary` from their canonical positions inside the `ase` output directory. |
| `--counts <ase_counts.tsv>` | Explicit path to the counts matrix from `ase`. Required if `--ase-dir` isn't given. |
| `--summary <ase_summary.tsv>` | Optional; enables filtering to informative transcripts. Falling through without it still works but tests every transcript. |
| `--output <file.R>` / `-o` | Output R script (default `ase_analysis.R`). |
| `--haplotype-names N1,N2` | Must match what `ase` used (default `hap0,hap1`). |
| `--fdr <n>` | BH significance threshold (default `0.05`). |
| `--min-reads <n>` | Skip `(transcript × sample)` pairs with total read count below this (default `20`). Below this the binomial test is uninformative — too few reads to resolve small departures from the null. |
| `--include-uninformative` | Test every transcript even when `ase_summary.tsv` marks them as uninformative (haplotype sequences identical). Off by default. |

## The binomial test

For each transcript × sample pair, the script runs:

```r
binom.test(hap0_count, hap0_count + hap1_count, p = null_p0[sample])
```

The null proportion `null_p0[sample]` is **sample-specific**, computed as:

```
null_p0[s] = sum(hap0 counts for sample s across all transcripts) /
             sum(total counts for sample s across all transcripts)
```

This is the key correction: if sample 1 maps 52 % of its reads to hap0 globally (e.g. because hap0 has better-quality assembly), the null is 0.52 for that sample, not 0.5. Without this correction, every transcript in sample 1 would appear to show ASE toward hap0.

BH adjustment (`p.adjust(..., method = "BH")`) is applied **per sample**, not pooled. Each sample is its own independent experiment.

## Filtering

Two filters run in order:

1. **Informative filter** (on by default, requires `--summary`): drops transcripts where both haplotypes have identical CDS sequence. Can't test what can't differ.
2. **Min-reads filter** (per-row): if `hap0_count + hap1_count < --min-reads` for a given `(transcript, sample)` pair, the row is kept in the output with `reason = "below_min_reads"` and `pvalue = NA`. It does not contribute to BH adjustment.

## What's in the emitted R script

Hand-maintained template. Structure:

```r
# Header with generation metadata
# Inputs block (edit here to re-run with different params)
# Load ase_counts.tsv
# Filter to informative transcripts if ASE_SUMMARY is set
# Parse <sample>.<hap> column names (splits on the LAST dot, so sample IDs can contain dots)
# Verify both haplotypes are present for every sample
# Compute per-sample null_p0
# Per-transcript-per-sample binom.test loop
# BH adjustment within each sample
# Write ase_results.tsv
# Per-sample imbalance histogram PDF
# sessionInfo() for reproducibility
```

Every substitution is a named R variable at the top (`ASE_COUNTS`, `HAP_NAMES`, `FDR_CUTOFF`, `MIN_READS_PER_SAMPLE`, `FILTER_UNINFORMATIVE`). Edit and re-run without regenerating.

## Reading `ase_results.tsv`

| column | meaning |
|---|---|
| `transcript` | Transcript ID from `ase_counts.tsv`. |
| `sample` | Sample ID. |
| `hap0_count`, `hap1_count` | Raw counts. Fractional (salmon's estimated counts) preserved. |
| `total` | `hap0_count + hap1_count`. |
| `hap0_frac` | `hap0_count / total`. |
| `imbalance` | `hap0_frac − null_p0[sample]`. Signed: positive = excess hap0 expression, negative = excess hap1. |
| `null_p` | This sample's null hap0 proportion. |
| `pvalue` | Two-sided binomial exact test p-value. `NA` when `total < MIN_READS_PER_SAMPLE`. |
| `padj` | BH-adjusted p-value within this sample. `NA` for skipped rows. |
| `reason` | `"tested"` or `"below_min_reads"`. |
| `significant` | `padj < FDR_CUTOFF`. `FALSE` for `NA` rows. |

Sort by `(sample, padj)` in R or shell — the script writes in transcript × sample order so you can join it back to other tables cleanly.

## What `ase-template` does NOT do

- **Doesn't run the test.** The tool writes the script. You run `Rscript ase_analysis.R`.
- **Doesn't do cis/trans regression.** For that, look at [Wittkopp et al. 2004 / Landry et al. 2005]-style per-transcript linear models — the emitted script is the *first* pass; keep the output of `ase_counts.tsv` around and model it yourself if you need cis/trans decomposition.
- **Doesn't handle >2 haplotypes.** Binomial testing on pairs only. Polyploid ASE (e.g. triploid) needs a multinomial / Dirichlet-multinomial formulation — out of scope for 0.5.0.
- **Doesn't do gene-level aggregation.** Operates at the transcript level, matching salmon's output. Summing isoforms by gene changes the null (Simpson's paradox territory); do it deliberately, not by default.

## Troubleshooting

- **"Haplotype labels in ase_counts.tsv (hap0,hap1) don't match --haplotype-names (paternal,maternal)"** — you passed custom `--haplotype-names` to `ase` but the default to `ase-template`, or vice versa. They must match.
- **"Sample 'X' is missing one or both haplotype columns"** — the counts matrix was hand-edited, or an earlier salmon run failed partway. Re-run `ase` to regenerate.
- **All `pvalue` are `NA`** — `--min-reads` is higher than your per-transcript totals. Lower it (e.g. `--min-reads 5`) for a smoke test, but don't over-interpret: the binomial test is weak at small N.
- **`null_p0[sample]` looks extreme (e.g. 0.8)** — one haplotype dominates global mapping. Check the `asymmetry_flag` column in `ase_bundle.json` from the `ase` run; if flagged, fix the phasing or assembly before trusting any per-transcript result.

## Related

- `myconote-cli ase` — produces the `ase_counts.tsv` + `ase_summary.tsv` this script consumes.
- `myconote-cli de-template` — same "emit R, user runs R" pattern for DESeq2 differential expression.
- Binomial exact test background: <https://en.wikipedia.org/wiki/Binomial_test>
- BH (Benjamini–Hochberg) FDR: <https://en.wikipedia.org/wiki/False_discovery_rate>
