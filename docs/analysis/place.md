# `place` — Y1000+ functional placement

`place` locates a user-annotated genome within the **Opulente et al. (2024)
Y1000+** reference of 1,154 sequenced yeasts by comparing KEGG-KO
functional profiles. Output: a ranked list of closest yeasts by Jaccard
similarity, plus — when the relevant Y1000+ subsets are installed —
auto-generated predictions of codon table, metabolic lifestyle,
thermotolerance, and ecological niche.

## Prerequisites

`place` consumes the TSV produced by `annotate`:

```bash
myconote-cli annotate genes.gff3 --fasta genome.fa ...
# → annotate_out/annotations.tsv  (contains one KEGG-KO column per row)
```

The Y1000+ reference bundle must be installed. The `kegg` subset is
required; the others are optional and enable extra predictions:

```bash
# Minimum: KEGG-KO reference distributions (~required)
myconote-cli setup --y1000plus --include kegg

# Recommended: pull the bundle subsets that unlock full place output
myconote-cli setup --y1000plus --include kegg,codontable,metabolism,phenotypes,environment

# Preview what each subset contains first
myconote-cli setup --y1000plus --list
```

## Usage

```
myconote-cli place --annotated <path.tsv> [options]

Options:
  --annotated <file.tsv>   myconote `annotate` output (required)
  --top <n>                How many closest species to report (default: 10)
  --format <human|tsv>     Output format (default: human)
```

## What it reports

### Top-N closest species (always shown)

KEGG-KO Jaccard similarity between your genome and each Y1000+ reference.
The top hit is your functional nearest neighbour in the 1,154-yeast
dataset — a useful starting point for inferring likely phenotypes,
genetic code, and lifestyle.

### Inferred genetic code (`codontable` subset)

If the `codontable` subset is installed, `place` reports the NCBI genetic
code table associated with the top hit. For most ascomycetes this is
**1 (standard)**; for Candida CTG-clade species it is typically **12**.

### Carbon / nitrogen lifestyle (`metabolism` subset)

Weighted vote across the top-N hits over the C/N lifestyle labels
annotated in the Y1000+ metabolism table. Useful for predicting whether
a newly sequenced yeast is likely to grow on specific substrates before
running wet-lab assays.

### Thermotolerance at 37 °C (`phenotypes` subset)

Weighted vote over the growth-at-37 labels (**Y / N / W / V / S**) in the
Y1000+ phenotype matrix. Useful for triaging species for fermentation or
pathogenesis studies.

### Ecological niche (`environment` subset)

Jaccard-weighted majority over the EnvO/OWL-derived isolation-source
labels. Predicts the most likely ecological niche (e.g. *insect gut*,
*rotting wood*, *flower nectar*) for your species.

## Examples

```bash
# Basic placement — KEGG-KO only
myconote-cli place --annotated annotate_out/annotations.tsv

# Top 20 hits, machine-readable
myconote-cli place --annotated annotate_out/annotations.tsv --top 20 --format tsv \
    > placement.tsv

# Full placement with all predictions — requires all subsets installed
myconote-cli setup --y1000plus --include kegg,codontable,metabolism,phenotypes,environment
myconote-cli place --annotated annotate_out/annotations.tsv --top 10
```

## Citation

> Opulente DA et al. (2024). *Genomic factors shape carbon and nitrogen metabolic niche breadth across Saccharomycotina yeasts*. **Science** 384(6694): eadj4503. <https://doi.org/10.1126/science.adj4503>

The Y1000+ reference data is redistributed with attribution via the
`setup --y1000plus` opt-in bundle; see the setup `--list` output for
per-subset licences.

## Caveats

- All predictions are **suggestive, not definitive**. They reflect the
  phenotypes of the functional-nearest yeasts in Y1000+, not direct
  measurements of your genome's behaviour.
- Non-yeast genomes (plants, animals, basidiomycete fungi, protists) can
  still be placed, but the returned neighbours and predictions will be
  biologically meaningless — the reference is ascomycete-biased.
- KEGG-KO profiles require that your `annotate` run populated the
  KEGG-KO column (via EggNOG-mapper). Check
  `annotate_out/annotation_report.txt` before running `place`.
