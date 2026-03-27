# Myconote_CLI

**Genome Annotation Pipeline** — a fast, modular command-line tool for annotating fungal (and other eukaryotic) genomes, from raw assembly to functional annotation.

> Developed by **Benjamin Narh-Madey** · Hittinger Lab · University of Wisconsin–Madison

---

## What is Myconote_CLI?

Myconote_CLI takes a genome assembly (FASTA) through a structured, reproducible annotation pipeline:

```
Sort → Mask → Train → Predict → Update → Annotate
```

Each step is a separate subcommand, giving you full control over where to start, pause, or resume. All intermediate files use standard formats (GFF3, FASTA, GenBank) compatible with Geneious, IGV, JBrowse2, and the UCSC Genome Browser.

---

## Key Features

- **Repeat masking** via RepeatModeler2 + RepeatMasker
- **Ab initio gene prediction** with Augustus, GlimmerHMM, and SNAP
- **Evidence-based consensus** via EvidenceModeler (EVM)
- **Functional annotation** — BLAST, eggNOG-mapper, dbCAN (CAZymes), InterProScan, antiSMASH
- **Phylogenomics** via IQ-TREE2 (UFBoot2 + ModelFinder Plus)
- **Visualisation** — circular genome plots, linear tracks, JBrowse2 integration
- **Multi-kingdom support** — fungi, plants, animals, protists, with sensible defaults per kingdom
- 39/39 integration tests passing on *Candida tropicalis* test data

---

## Quick Example

```bash
# Run the full pipeline on a fungal genome
myconote sort   --genome mygenome.fas --out sorted/
myconote mask   --genome sorted/genome.fas --out masked/
myconote train  --genome masked/genome.fas --rna-bam rna.bam --out training/
myconote predict --genome masked/genome.fas --training training/ --out predictions/
myconote update  --genome masked/genome.fas --predictions predictions/ --rna-bam rna.bam --out updated/
myconote annotate --genome masked/genome.fas --gff updated/final.gff3 --out annotation/
```

---

## Getting Started

- [Installation](installation.md) — one-shot installer, conda, or pre-built binary
- [Quick Start](quickstart.md) — annotate a genome in under an hour
- [Workshop Lesson](lesson.md) — full Carpentries-style tutorial (~3 hours)
