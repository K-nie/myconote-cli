# Quick Start

This guide walks you through annotating a fungal genome from start to finish. Estimated time: 45–90 minutes (excluding database downloads).

---

## 1. Prepare your genome

Place your genome assembly FASTA in a working directory:

```bash
mkdir my_annotation && cd my_annotation
cp /path/to/mygenome.fas .
```

---

## 2. Sort and validate scaffolds

```bash
myconote sort --genome mygenome.fas --out 01_sorted/
```

This renames scaffolds to a consistent format, removes sequences below a minimum length (default 500 bp), and produces a summary report.

---

## 3. Mask repeats

```bash
myconote mask --genome 01_sorted/genome.fas --out 02_masked/ --threads 8
```

Runs RepeatModeler2 to build a de novo repeat library, then RepeatMasker to soft-mask the genome. Soft-masked bases are lowercased; gene predictors can still read through them.

---

## 4. Train gene predictors

```bash
myconote train \
  --genome 02_masked/genome.masked.fas \
  --rna-bam rnaseq_aligned.bam \
  --species my_fungus \
  --out 03_training/
```

!!! tip
    If you don't have RNA-seq data, omit `--rna-bam`. Augustus will use its built-in *Saccharomyces cerevisiae* model as the starting point.

---

## 5. Predict gene models

```bash
myconote predict \
  --genome 02_masked/genome.masked.fas \
  --training 03_training/ \
  --kingdom fungi \
  --out 04_predictions/
```

Runs Augustus, GlimmerHMM, and SNAP in parallel, then merges results using EvidenceModeler (EVM).

---

## 6. Update with RNA evidence

```bash
myconote update \
  --genome 02_masked/genome.masked.fas \
  --predictions 04_predictions/ \
  --rna-bam rnaseq_aligned.bam \
  --out 05_updated/
```

---

## 7. Functional annotation

```bash
myconote annotate \
  --genome 02_masked/genome.masked.fas \
  --gff 05_updated/final.gff3 \
  --kingdom fungi \
  --out 06_annotation/
```

Runs BLAST, eggNOG-mapper, dbCAN (CAZymes), and antiSMASH (secondary metabolite clusters).

---

## 8. View results

```bash
# Summary statistics
myconote stats --gff 06_annotation/final.annotated.gff3

# Open in JBrowse2
myconote view jbrowse --genome 02_masked/genome.masked.fas --gff 06_annotation/final.annotated.gff3
```

---

## Output files

| File | Description |
|------|-------------|
| `final.annotated.gff3` | Full gene models with functional annotations |
| `proteins.faa` | Predicted protein sequences (FASTA) |
| `transcripts.fna` | Predicted transcript sequences (FASTA) |
| `annotation_summary.txt` | Gene count, N50, BUSCO scores |
| `cazymes.tsv` | CAZyme assignments |
| `antismash/` | Secondary metabolite cluster predictions |

---

## Next steps

- Run a phylogenomic analysis: [phylogeny](analysis/phylogeny.md)
- Visualise your genome: [plot](analysis/plot.md)
- Follow the full workshop tutorial: [Workshop Lesson](lesson.md)
