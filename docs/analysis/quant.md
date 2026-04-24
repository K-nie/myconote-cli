# myconote quant

Quantify RNA-seq expression against an annotated fungal genome. Reads pass through `fastp` for QC and adapter trimming, then `salmon` with a decoy-aware index produces per-sample `quant.sf` files plus wide count and TPM matrices. Differential expression stays in R (DESeq2 / edgeR / limma) — `quant` stops at the count matrix on purpose.

## Quickstart

```bash
myconote-cli predict genome_masked.fa --kingdom fungi
myconote-cli annotate predict_out/consensus.gff3 --fasta genome.fa --kingdom fungi

# Derive the transcriptome for salmon
myconote-cli convert annotate_out/annotated.gff3 --to cds \
    --fasta genome.fa -o cds.fa

# Quantify
myconote-cli quant cds.fa \
    --samples samples.tsv \
    --genome genome.fa \
    --output quant_out \
    --threads 16
```

Then ~10 lines of R to run DESeq2. The `quant_out/salmon/<sample>/quant.sf` files are consumed directly by `tximport`.

## Sample sheet

TSV, header required:

| column         | required | meaning                                                           |
| -------------- | -------- | ----------------------------------------------------------------- |
| `sample_id`    | yes      | unique identifier, `[A-Za-z0-9_-]+`                               |
| `fastq_r1`     | yes      | R1 FASTQ path (gzipped or plain)                                  |
| `fastq_r2`     | no       | R2 for paired-end; blank or missing → single-end                  |
| `condition`    | no       | DE grouping label; passed through to the bundle                   |
| `strandedness` | no       | `unstranded` / `forward` / `reverse` / `auto` (default `auto`)    |
| `batch`        | no       | technical batch label                                             |

Relative FASTQ paths resolve against the directory containing the sheet. Comment lines starting `#` and blank lines are ignored. Unknown columns are preserved.

Example:

```tsv
sample_id	fastq_r1	fastq_r2	condition	strandedness
WT_rep1	fastq/WT_rep1_R1.fq.gz	fastq/WT_rep1_R2.fq.gz	control	reverse
WT_rep2	fastq/WT_rep2_R1.fq.gz	fastq/WT_rep2_R2.fq.gz	control	reverse
WT_rep3	fastq/WT_rep3_R1.fq.gz	fastq/WT_rep3_R2.fq.gz	control	reverse
C_lim1	fastq/C_lim1_R1.fq.gz	fastq/C_lim1_R2.fq.gz	carbon_limited	reverse
C_lim2	fastq/C_lim2_R1.fq.gz	fastq/C_lim2_R2.fq.gz	carbon_limited	reverse
C_lim3	fastq/C_lim3_R1.fq.gz	fastq/C_lim3_R2.fq.gz	carbon_limited	reverse
```

## Outputs

Under `--output` (default `quant_out/`):

| path                             | purpose                                                        |
| -------------------------------- | -------------------------------------------------------------- |
| `counts.tsv`                     | wide estimated-counts matrix (transcript × sample, sorted rows) |
| `tpm.tsv`                        | wide TPM matrix (same layout)                                  |
| `salmon/<sample>/quant.sf`       | per-sample salmon output — tximport consumes this directly     |
| `salmon/<sample>/logs/...`       | salmon run logs (mapping rate, EM iterations)                  |
| `fastp/<sample>.json`            | full fastp QC report per sample                                |
| `fastp/<sample>.html`            | browsable fastp report (per-base quality, adapter content)     |
| `quant_bundle.json`              | fixed-schema reproducibility manifest                          |
| `sample_sheet.tsv`               | copy of the input sheet                                        |

### `quant_bundle.json`

One-screen summary of the run: tool versions (fastp, salmon), SHA256 of every input file, salmon index key and sizing, and a per-sample block with mapping rate, library size, and nine fastp QC fields (reads before/after, Q30 before/after, pass %, adapter trimmed reads/bases, duplication rate, insert-size peak). The full fastp JSON stays under `fastp/` for deep dives.

## Salmon index cache

The decoy-aware index is expensive to build and cheap to reuse. `quant` caches it by a content hash of `(cds.fa, genome.fa, salmon version, k)`. Cache directory resolution in precedence order:

1. `--index-cache <dir>` (CLI flag)
2. `MYCONOTE_INDEX_CACHE` (env var)
3. `$XDG_CACHE_HOME/myconote/salmon_index/`
4. `~/.cache/myconote/salmon_index/` (Linux default when XDG is unset)
5. `~/.myconote/salmon_index/` (legacy macOS fallback)

No silent fallback to `/tmp` — a silent location shuffle would break the bundle's reproducibility claim. If the resolved directory is read-only, `quant` fails fast and names the path that failed.

Pass `--rebuild-index` to delete the cache entry for this input tuple and rebuild.

## Library strandedness

Salmon auto-detects strandedness from the first few thousand reads when `strandedness = auto` in the sheet. For Illumina TruSeq stranded dUTP protocols (the common fungal RNA-seq prep) the correct value is `reverse` — set it explicitly once and skip the auto-detection overhead. `quant` passes each sample's strandedness to salmon as `--libType`; paired vs single-end is inferred from the presence of `fastq_r2`.

## Parallelism

0.3.0 ships with `--jobs 1` (one salmon process at a time, each using `--threads`). Outer parallelism across samples is not yet benchmarked — the A. niger 6-sample benchmark on a 16-core host is pending. `--jobs N` is accepted but currently logs a warning and runs serially; a real default will land in a follow-up release.

## Required external tools

- `salmon` ≥ 1.10 — `conda install -c bioconda salmon`
- `fastp` ≥ 0.23 — `conda install -c bioconda fastp`

Both are registered in `myconote-cli install` and `myconote-cli check` under `used_by = "quant"`. Run `myconote-cli check` to verify both are on PATH.

## Feeding into DESeq2

```r
library(tximport)
library(DESeq2)

samples <- read.table("quant_out/sample_sheet.tsv", header = TRUE, sep = "\t")
files <- file.path("quant_out", "salmon", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

# For a gene-level summary: supply a tx2gene table. If you only want
# transcript-level results, set txOut = TRUE and skip tx2gene.
txi <- tximport(files, type = "salmon", txOut = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
```

## Troubleshooting

- **Mapping rate < 50 %** — usually adapter contamination or the wrong transcriptome. Check `fastp/<sample>.html`; a right-skewed per-base-quality plot on read 3′ ends flags the former. Verify the `cds.fa` was derived from the same genome build as the reads.
- **`"not writable"` error on cache dir** — `$HOME` is quota-full on an HPC node. Pass `--index-cache <writable dir>` or set `MYCONOTE_INDEX_CACHE` to a scratch path.
- **All samples report `library_size = 0`** — no reads mapped. Almost always a strandedness mismatch or a cds / genome pair that don't correspond. Re-check by setting `strandedness = auto` for one sample and comparing the logged `--libType` detection against your sheet.
- **`salmon index` crashes on k=31 with very small fixtures** — increase k-mer permissible minimum: use `-k 21` for tutorial-scale inputs.

## What this does not do

- **DE analysis** (intentional — 10 lines of R with DESeq2 is the standard, and wrapping it in the tool would drag in R as a runtime dependency).
- **scRNA-seq** (different tool — single-cell needs Cell Ranger / STARsolo / kallisto|bustools and downstream Scanpy / Seurat; bulk and single-cell are different data shapes, not just different scales).
- **Variant-aware / allele-specific quant** (salmon ReefClusters or a phased-alignment pipeline — specialized enough that the right answer is almost always a dedicated tool, not a shoe-horn into bulk quant).
- **HISAT2-backed alignment for BAM output** — not in scope; users who need a BAM for IGV review should run HISAT2 / STAR externally on the same FASTQs.
- **Read ingestion** was in this list until 0.3.1; `myconote-cli fetch-rna` now covers SRA / ENA accession download.
