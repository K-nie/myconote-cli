# myconote go-template

Generate a ready-to-run R script that performs **Gene Ontology enrichment** on a `de-template` results table using `topGO`'s classic Fisher's exact test. The tool writes the script; you run it with `Rscript`.

Same Option 1D pattern as `de-template` and `ase-template` — write R, the user runs R, no R runtime inside MycoNote-CLI.

## Prerequisites — install before running the emitted script

R ≥ 4.0 plus Bioconductor's `topGO`. Install once:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("topGO"))
```

If `Rscript` isn't on your PATH, the tool warns but still writes the script — you can run it on another host.

## Quickstart

After running `de-template` (which produces `de_<factor>_<level1>_vs_<level2>.tsv`) and `annotate` (which produces `annotations.tsv` with a `go_terms` column):

```bash
myconote-cli go-template \
    --de-results de_condition_treated_vs_control.tsv \
    --annotations annotate_out/annotations.tsv \
    -o go_enrichment.R
Rscript go_enrichment.R
```

You get one TSV per ontology and a combined dot plot:

- `go_BP_enrichment.tsv` — biological process terms (one row per term, sorted by classic Fisher *p*).
- `go_MF_enrichment.tsv` — molecular function.
- `go_CC_enrichment.tsv` — cellular component.
- `go_dotplot.pdf` — one panel per ontology, top-N terms by *p*-value, dot size scales with significant gene count, dots in red when BH-adjusted *p* < 0.05.

## Inputs

| flag | meaning |
|---|---|
| `--de-results <tsv>` | DE results TSV from `de-template` (must contain a `padj` column and an identifier column — `transcript` by default). Required. |
| `--annotations <tsv>` | Functional annotation TSV from `annotate` (must contain a gene-ID column — `locus_tag` by default — and a GO terms column — `go_terms` by default). Required. |
| `--output <file.R>` / `-o` | Output R script (default `go_enrichment.R`). |
| `--fdr <n>` | Threshold on DE `padj` for marking a gene "significant" in the foreground set (default `0.05`). |
| `--ontology <BP\|MF\|CC\|all>` | Which ontology branch(es) to test (default `all`). When `all`, three TSVs and one combined PDF are written. |
| `--top <n>` | Top-N enriched terms shown in the dot plot per ontology (default `30`). The TSVs always carry every term that topGO tested. |
| `--de-id-col <name>` | DE identifier column (default `transcript`). |
| `--ann-id-col <name>` | Annotations identifier column (default `locus_tag`). Set this to whatever shared identifier you have between the two TSVs. |
| `--go-col <name>` | GO terms column (default `go_terms`). |
| `--go-separator <str>` | Delimiter inside the GO column (default `|`). Set to `,` if you converted from EggNOG output, etc. |

## What the script does

1. **Reads `DE_RESULTS`**, validates that `padj` and the chosen ID column exist, and pulls them out.
2. **Reads `ANNOTATIONS`**, validates that the chosen ID and GO columns exist.
3. **Builds the `gene -> GO` map** by splitting each cell on `GO_SEPARATOR`, trimming whitespace, and keeping only well-formed GO IDs (`GO:` followed by 7 digits). Malformed entries are silently dropped here so the topGO call doesn't fail downstream — if you suspect your annotations are malformed, inspect the TSV directly.
4. **Defines the universe** as the intersection of DE IDs and annotated genes. Genes in the DE table without any well-formed GO term are not in the universe (they would be uninformative for enrichment).
5. **Defines the foreground** as the subset of the universe with `padj < FDR_CUTOFF`.
6. **Runs `topGO`** with `algorithm = "classic"` and `statistic = "fisher"` for each ontology in `ONTOLOGIES`. `nodeSize = 5` (drops GO terms annotated to fewer than 5 universe genes — the typical default).
7. **BH-adjusts** the per-term *p*-values within each ontology (this is `classic_padj`; the raw classic *p* is also kept).
8. **Writes one TSV per ontology** with `(GO.ID, Term, Annotated, Significant, Expected, classic, classic_p, classic_p_raw, classic_padj, ontology)`, sorted by `classic_p`.
9. **Writes a dot plot PDF** with one panel per ontology showing the top-N terms.

The "classic" Fisher's test is the most interpretable choice — straight 2×2 (significant vs not, in-term vs not-in-term). topGO's `elim` and `weight01` algorithms decorrelate parent/child terms (often what you actually want for publication-grade GO enrichment); to switch, edit `algorithm = "classic"` near the bottom of the emitted script.

## Joining DE results to annotations

The script does an inner join on `DE_ID_COL == ANN_ID_COL`. By default, `de-template` emits `transcript` IDs (mRNA-level from your CDS FASTA) and `annotate` emits `locus_tag` IDs (gene-level). These often **don't match** out of the box — you have a few options:

- **If your transcripts are 1:1 with genes** (most fungi without alternative splicing), set `--de-id-col transcript --ann-id-col transcript`. You may need to add a `transcript` column to your annotations TSV first (e.g. with `awk` from the GFF3).
- **If your locus tags match transcript IDs after stripping a suffix** (e.g. `GENE_001` ↔ `GENE_001.1`), pre-process the DE TSV to drop the suffix before passing it in.
- **If your DE was run at the gene level** with `tximport(..., tx2gene = ...)`, change the DE-template invocation to emit gene IDs and then set `--de-id-col gene_id --ann-id-col locus_tag` here.

The script will warn loudly with `No overlap between DE_ID_COL values and ANN_ID_COL values` if it can't join — at which point check what's actually in your two TSVs (`head -1` on each).

## Reading the enrichment TSV

| column | meaning |
|---|---|
| `GO.ID` | GO term ID (e.g. `GO:0006281`). |
| `Term` | Human-readable description. |
| `Annotated` | Genes in the universe annotated to this term. |
| `Significant` | Genes in the foreground (DE-significant) annotated to this term. |
| `Expected` | Expected count under the null. |
| `classic` | Classic Fisher *p*-value as topGO returns it (string; can be `<1e-30` for very small values). |
| `classic_p` | Numeric form of `classic` (NA when topGO returned a `<1e-30`-style string). |
| `classic_p_raw` | The original string form, preserved. |
| `classic_padj` | BH-adjusted classic *p* across all tested terms in this ontology. |
| `ontology` | `BP`, `MF`, or `CC`. |

## What `go-template` does NOT do

- **Doesn't run the test.** The tool writes the script. You run `Rscript go_enrichment.R`.
- **Doesn't fetch GO obo / goa files.** topGO uses the gene → GO map you provide; it doesn't look up parents or ancestors against a separate ontology file. (This is the standard topGO workflow with `annFUN.gene2GO`.)
- **Doesn't do GSEA / pre-ranked enrichment.** Classic over-representation only. For ranked enrichment use `fgsea` and pre-rank by `log2FoldChange × -log10(padj)`.
- **Doesn't summarise across genes-per-locus.** Operates on the IDs you pass in. If your DE was at the transcript level, your enrichment is at the transcript level — be deliberate about that.

## Troubleshooting

- **`No overlap between DE_ID_COL values and ANN_ID_COL values`** — your DE and annotations TSVs use different identifiers. Check `head -1` on each and either re-do DE with matching IDs or set `--de-id-col` / `--ann-id-col` appropriately.
- **`No GO terms recovered from <annotations.tsv>`** — the `go_terms` column was empty, or the separator is wrong. Inspect the column manually and set `--go-col` / `--go-separator` to match.
- **`No genes pass FDR_CUTOFF`** — your DE results have no significant genes at the threshold. Lower `--fdr` or check the DE first.
- **Many terms have `Annotated < 5`** — `nodeSize = 5` filters them out before testing. If you want them in, edit the script to `nodeSize = 1`.

## Related

- `myconote-cli de-template` — produces the DE results TSV this script consumes.
- `myconote-cli annotate` — produces the `annotations.tsv` this script consumes.
- `topGO` reference: <https://bioconductor.org/packages/release/bioc/html/topGO.html>
- Alexa & Rahnenführer (2009), *Gene set enrichment analysis with topGO*.
- Benjamini–Hochberg FDR: <https://en.wikipedia.org/wiki/False_discovery_rate>
