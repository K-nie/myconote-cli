# Myconote_CLI

**Blazing-fast genome annotation pipeline** -- as genome sequencing becomes cheaper and long-read assemblies become routine, the bottleneck in genomics has shifted from sequencing to annotation. Existing pipelines are slow, narrowly scoped, and produce outputs that require extensive manual cleanup before submission. myconote-cli addresses this gap: a high-performance Rust CLI that takes a eukaryotic genome assembly from raw contigs to NCBI-ready submission, integrating 15 annotation sources across 5 kingdoms with built-in validation and reproducibility tracking.

> Developed by **Benjamin Narh-Madey** - Hittinger Lab, Laboratory of Genetics - UW-Madison

---

## What is Myconote_CLI?

Myconote_CLI takes a genome assembly (FASTA) through a structured, reproducible annotation pipeline:

```
Sort -> Mask -> Train -> Predict -> Update -> Annotate -> Submit
```

Each step is a separate subcommand, giving you full control over where to start, pause, or resume. All intermediate files use standard formats (GFF3, FASTA, GenBank) compatible with Geneious, IGV, JBrowse2, and UCSC Genome Browser.

---

## Key Features

- **23 CLI commands** covering the full genome annotation lifecycle
- **15 annotation sources** -- MMseqs2, Pfam (hmmsearch), InterProScan, EggNOG, CAZyme, MEROPS, BUSCO, antiSMASH, tRNAscan-SE, secretome, GO terms, and more
- **Multi-tool gene prediction** -- Augustus, SNAP, GlimmerHMM, GeneMark, miniprot protein evidence, configurable EVM consensus
- **AI-powered interpretation** -- `myconote-cli explain` uses a local LLM (Ollama) to interpret results with grounded citations and command recommendations. No data leaves your machine.
- **Multi-genome batch mode** -- `myconote-cli batch` annotates many genomes from a directory or sample sheet with HTCondor support, resume capability, and a live dashboard
- **18 NCBI genetic code tables** -- Candida CTG clade, mitochondrial genomes, etc.
- **NCBI submission prep** -- GFF3 validation, .tbl generation, table2asn integration
- **Ploidy awareness** -- allelic duplicate detection for polyploid genomes
- **Comparative genomics** -- phylogenetics (IQ-TREE), synteny diagrams, multi-genome comparison
- **Visualization** -- circular/linear genome plots, JBrowse2, UCSC browser
- **15+ format conversions** -- GFF3/GTF/BED/GenBank/FASTA/FASTQ/PHYLIP/NEXUS/VCF
- **Interactive tutorial** -- `myconote-cli learn` (8 lessons, swirl-style)
- **5 kingdoms supported** -- fungi, plants, animals, insects, protists
- **Reproducibility** -- workflow reports (JSON), database version tracking, output validation
- **Docker + Singularity** containers for HPC and cloud environments
- **CI/CD** -- GitHub Actions (fmt, clippy, test, build, security audit)
- **332 tests** passing on real genome data (*Brettanomyces bruxellensis*, *Candida tropicalis*)

---

## Quick Example

```bash
# Install tools and databases
myconote-cli install --yes
myconote-cli setup

# Run the full pipeline
myconote-cli sort assembly.fa --min-length 500
myconote-cli mask assembly_sorted.fa --engine repeatmodeler --threads 8
myconote-cli predict assembly_masked.fa --kingdom fungi --locus-prefix MYORG
myconote-cli annotate predict_out/consensus.gff3 --fasta assembly.fa \
    --trnascan --interproscan --email you@email.edu
myconote-cli submit annotate_out/annotated.gff3 --fasta assembly.fa \
    --organism "Genus species"
```

---

## Performance

Tested on *Brettanomyces bruxellensis* (12.9 Mb, 30 contigs):

| Step | Time | Result |
|------|------|--------|
| sort | <1s | 30 contigs sorted |
| predict (Augustus) | ~3 min | 5,218 genes |
| annotate (MMseqs2) | ~1.5 min | 75% genes with product names |
| annotate (Pfam/hmmsearch) | ~7 min | 85.6% genes with domains |
| annotate (InterProScan) | ~5 min (cached) | 98.3% genes annotated |
| NCBI validation | <1s | PASSED (0 errors) |

Pfam uses hmmsearch instead of hmmscan for a ~6x speedup on large proteomes.

---

## Getting Started

- [Installation](installation.md) -- Docker, Singularity, conda, or manual build
- [Quick Start](quickstart.md) -- annotate a genome in under an hour
- [Workshop Lesson](lesson.md) -- full tutorial (~3 hours)
- **Interactive tutorial**: `myconote-cli learn` -- 8 self-paced lessons in your terminal
- [Explain (AI interpreter)](analysis/explain.md) -- local LLM-powered result interpretation
- [Batch annotation](analysis/batch.md) -- multi-genome annotation with HTCondor support
