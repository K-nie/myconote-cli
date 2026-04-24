# Myconote_CLI

**Blazing-fast fungal genome annotation pipeline** -- as fungal genome sequencing becomes cheaper and long-read assemblies routine, the bottleneck has shifted from sequencing to annotation. Existing pipelines are slow, fragile, and leave outputs that take days of manual cleanup before NCBI will accept them. myconote-cli addresses this gap: a high-performance Rust CLI that takes a fungal genome assembly from raw contigs to an NCBI-ready submission, integrating 15 annotation sources with built-in validation and reproducibility tracking.

The tool runs on other eukaryotes, but defaults, benchmarks, and the `learn` tutorial are tuned for fungi.

> Developed by **Benjamin Narh-Madey** - Hittinger Lab, Laboratory of Genetics - UW-Madison

---

## What is Myconote_CLI?

Myconote_CLI takes a genome assembly (FASTA) through a structured, reproducible annotation pipeline:

```
Sort -> Mask -> Train -> Predict -> Update -> Annotate -> Submit
```

Each step is a separate subcommand, giving you full control over where to start, pause, or resume. All intermediate files use standard formats (GFF3, FASTA, GenBank) compatible with Proksee (web), IGV (desktop), and clinker (cross-species gene-cluster synteny) — the tools MycoNote hands off to rather than reimplementing.

RNA-seq is a first-class input at four distinct points: `train` and `update` use reads as *evidence for gene models*; `quant` uses reads for *expression quantification* (salmon + decoy-aware index, tximport-ready outputs); `fetch-rna` pulls reads by SRA/ENA accession so the entire RNA-seq stack runs end-to-end from an accession list.

---

## Key Features

- **23 CLI commands** covering the full fungal genome annotation lifecycle
- **15 annotation sources** -- MMseqs2, Pfam (hmmsearch), InterProScan, EggNOG, CAZyme, MEROPS, BUSCO, antiSMASH, tRNAscan-SE, secretome, GO terms, and more
- **Multi-tool gene prediction** -- Augustus, SNAP, GlimmerHMM, GeneMark, miniprot protein evidence, configurable EVM consensus
- **RNA-seq expression quantification** -- `quant` wraps fastp + salmon with a decoy-aware index, emits tximport-ready `quant.sf` per sample plus wide count and TPM matrices
- **DESeq2 differential expression** -- `de-template` emits a ready-to-run R script (tximport → DESeq2 → apeglm shrinkage → volcano + MA plots); no R runtime dependency in the CLI itself
- **SRA / ENA ingestion** -- `fetch-rna` pulls public RNA-seq by run, study, BioProject, sample, or experiment accession; ENA REST default, sra-toolkit fallback
- **AI-powered interpretation** -- `explain` uses a local LLM (Ollama) to interpret results with grounded citations and command recommendations. No data leaves your machine.
- **Multi-genome batch mode** -- `batch` annotates many genomes from a directory or sample sheet with HTCondor support, resume capability, and a live dashboard
- **25 NCBI genetic code tables** -- tables 1–6, 9–14, 16, 21–31, 33, every table verified against the NCBI reference. Critical for *Candida* CTG clade (Table 12) and other yeast-specific codes; includes caveats for context-dependent tables 27, 28, 31
- **NCBI submission prep** -- GFF3 validation, .tbl generation, table2asn integration
- **Ploidy awareness** -- allelic duplicate detection for polyploid genomes
- **N-genome orthology** -- `compare` wraps OrthoFinder with tiered genome-count caps tuned for fungal clades (5 fungi / 3 small plants / 2 large)
- **15+ format conversions** -- GFF3/GTF/BED/GenBank/FASTA/FASTQ/PHYLIP/NEXUS/VCF, plus `convert --to cds` for spliced transcriptome extraction
- **Interactive tutorial** -- `myconote-cli learn` (swirl-style, terminal-based)
- **Reproducibility** -- workflow reports (JSON), per-run bundles for quant and explain, database version tracking, output validation
- **Docker + Singularity** containers for HPC and cloud environments
- **CI/CD** -- GitHub Actions (fmt, clippy, test, build, security audit)
- **390+ tests** passing on real fungal genome data (*Brettanomyces bruxellensis*, *Candida tropicalis*)

### Multi-kingdom note

`--kingdom plant|animal|insect|protist` continues to work because the non-fungal Augustus species lists, BUSCO lineages, and intron-size ranges are still registered in the code. They are **not** validated on large non-fungal genomes and are not a supported use case at this release. Treat them as experimental; for serious plant / animal annotation, BRAKER or MAKER are the right tools.

---

## Quick Example

```bash
# Install tools and databases
myconote-cli install --yes
myconote-cli setup

# Annotate a fungal genome end-to-end
myconote-cli sort assembly.fa --min-length 500
myconote-cli mask assembly_sorted.fa --engine repeatmodeler --threads 8
myconote-cli predict assembly_masked.fa --kingdom fungi --locus-prefix MYORG
myconote-cli annotate predict_out/consensus.gff3 --fasta assembly.fa \
    --trnascan --interproscan --email you@email.edu
myconote-cli submit annotate_out/annotated.gff3 --fasta assembly.fa \
    --organism "Genus species" --genetic-code 12  # table 12 for Candida CTG clade

# Optional: quantify RNA-seq expression against the annotation
myconote-cli convert annotated.gff3 --to cds --fasta assembly.fa -o cds.fa
myconote-cli fetch-rna SRR12345678 SRR12345679 -o rna/
myconote-cli quant cds.fa --samples rna/samples.tsv --genome assembly.fa
myconote-cli de-template --quant-dir quant_out --design '~ condition' \
    --contrast 'condition,treated,control' -o analysis.R
Rscript analysis.R   # needs R + tximport + DESeq2 + apeglm installed
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

Pfam uses hmmsearch instead of hmmscan for a ~6× speedup on large proteomes.

---

## Getting Started

- [Installation](installation.md) -- Docker, Singularity, conda, or manual build
- [Quick Start](quickstart.md) -- annotate a genome in under an hour
- [Workshop Lesson](lesson.md) -- full tutorial (~3 hours)
- **Interactive tutorial**: `myconote-cli learn` -- self-paced lessons in your terminal
- [Explain (AI interpreter)](analysis/explain.md) -- local LLM-powered result interpretation
- [Batch annotation](analysis/batch.md) -- multi-genome annotation with HTCondor support
- [RNA-seq quantification](analysis/quant.md) -- fastp + salmon + tximport
- [DE analysis template](analysis/de-template.md) -- generates DESeq2 R script
- [SRA/ENA ingestion](analysis/fetch-rna.md) -- download FASTQ by accession
