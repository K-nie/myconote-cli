# Myconote_CLI

**Blazing-fast genome annotation pipeline** — as genome sequencing becomes cheaper and long-read assemblies become routine, the bottleneck in genomics has shifted from sequencing to annotation. Existing pipelines are slow, narrowly scoped, and produce outputs that require extensive manual cleanup before submission. myconote-cli addresses this gap: a high-performance Rust CLI that takes a eukaryotic genome assembly from raw contigs to NCBI-ready submission, integrating 15 annotation sources across 5 kingdoms with built-in validation and reproducibility tracking.

> Developed by **Benjamin Narh-Madey** · Hittinger Lab, Laboratory of Genetics · UW-Madison
> narhmadey@wisc.edu

---

## Overview

Myconote_CLI takes a genome assembly (FASTA) through a structured, reproducible pipeline:

```
Sort → Mask → Train → Predict → Update → Annotate → Submit
```

Each step is a separate subcommand. All intermediate files use standard formats (GFF3, FASTA, GenBank) compatible with Geneious, IGV, JBrowse2, and UCSC Genome Browser.

**New to myconote?** Run `myconote-cli learn` for an interactive, swirl-style tutorial right in your terminal.

---

## Installation

### Prerequisites

- **macOS** (Apple Silicon or Intel) or **Linux** (x86-64)
- [Miniconda or Anaconda](https://docs.conda.io/en/latest/miniconda.html)
- Rust >= 1.74

### Option A: Docker (recommended for reproducibility)

```bash
docker pull ghcr.io/k-nie/myconote-cli:latest
docker run -v $(pwd):/data myconote-cli predict /data/genome.fa --kingdom fungi
```

### Option B: Singularity (HPC clusters)

```bash
singularity pull myconote-cli.sif docker://ghcr.io/k-nie/myconote-cli:latest
singularity run myconote-cli.sif predict genome.fa --kingdom fungi
```

### Option C: One-shot installer

```bash
bash install.sh
```

Installs Miniconda (if needed), creates a conda environment, builds the binary, installs all 30+ external tools, and downloads annotation databases.

### Option D: Manual build

```bash
cargo build --release
sudo cp target/release/myconote-cli /usr/local/bin/
myconote-cli install --yes    # installs external tools via conda
myconote-cli setup            # downloads annotation databases
```

---

## Verifying your installation

```bash
myconote-cli --version      # shows banner + version
myconote-cli check          # reports status of all 30+ external tools
myconote-cli setup --check  # shows annotation database download status
```

---

## Quick start

```bash
# 1. Sort and rename contigs (longest first)
myconote-cli sort genome.fa --min-length 500

# 2. Soft-mask repeats
myconote-cli mask genome_sorted.fa --engine repeatmodeler --threads 8

# 3. Train gene predictors with RNA-seq (optional but recommended)
myconote-cli train genome_masked.fa --left R1.fq --right R2.fq --species myorg

# 4. Predict genes
myconote-cli predict genome_masked.fa --kingdom fungi --locus-prefix MYORG

# 5. Functionally annotate
myconote-cli annotate predict_out/consensus.gff3 --fasta genome.fa \
    --trnascan --interproscan --email you@email.edu

# 6. Validate and prepare NCBI submission
myconote-cli submit annotate_out/annotated.gff3 --fasta genome.fa \
    --organism "Genus species" --locus-prefix MYORG
```

---

## Pipeline commands (run in order)

| Command | Description |
|---------|-------------|
| `sort` | Sort contigs by length, rename headers, filter short scaffolds |
| `mask` | Identify and soft-mask repeats (5 engines: self, repeatmasker, repeatmodeler, both, full) |
| `train` | RNA-seq-guided training of Augustus and SNAP via Trinity + PASA |
| `predict` | Multi-tool gene prediction (Augustus + SNAP + GlimmerHMM + GeneMark + protein evidence + EVM consensus) |
| `update` | Refine gene models with RNA-seq evidence (PASA or lightweight UTR extension) |
| `annotate` | Functional annotation (MMseqs2, Pfam, InterProScan, EggNOG, CAZyme, MEROPS, BUSCO, antiSMASH, tRNAscan-SE, secretome) |
| `submit` | NCBI GenBank submission prep (GFF3 validation + .tbl + table2asn + .sqn) |

---

## Analysis commands

| Command | Description |
|---------|-------------|
| `stats` | Gene counts, lengths, GC content, N50, isoform stats (JSON / CSV / human-readable) |
| `plot` | Genome maps (linear PNG, circular PNG) |
| `phylogeny` | Maximum-likelihood tree with IQ-TREE 2 (ModelFinder + UFBoot) |
| `compare` | Multi-genome comparison (MMseqs2 / BLAST / MUMmer + phylogeny + synteny) |
| `view` | Interactive genome browser (JBrowse2 or UCSC custom track) |
| `synteny` | Ribbon diagram comparing two genomes (minimap2-based) |
| `convert` | Format conversion: GFF3 <-> GTF / BED / GenBank / protein FASTA; FASTA <-> FASTQ / PHYLIP / NEXUS; VCF conversions |
| `clean` | Validate and repair GFF3 annotation files |
| `fix` | Repair errors in GenBank (.gbk) files |
| `blast` | NCBI BLAST searches via remote API |
| `align` | Sequence alignment (BLAST, MMseqs2, MUMmer, minimap2) |

---

## Utility commands

| Command | Description |
|---------|-------------|
| `install` | Install missing tools via conda/mamba (30+ bioinformatics tools) |
| `check` | Check which external tools are installed and their versions |
| `setup` | Download and index annotation databases (Swiss-Prot, Pfam, EggNOG, BUSCO, dbCAN, MEROPS) |
| `remote` | Submit proteins to remote servers (Phobius, InterProScan, DeepLoc) |
| `species` | List all Augustus species models (grouped by kingdom) |
| `learn` | Interactive tutorial system -- 8 lessons, swirl-style, right in your terminal |

---

## Annotation sources

| Source | Tool | What it finds |
|--------|------|---------------|
| Swiss-Prot homology | MMseqs2 | Product names, UniProt accessions |
| Pfam domains | hmmsearch (6x faster than hmmscan) | Protein domain architecture |
| InterProScan | EBI REST API (cached) | InterPro, TIGRFAM, Gene3D, SMART, Superfamily |
| GO terms | UniProt + InterProScan | Gene Ontology functional categories |
| BUSCO | BUSCO 5 | Genome/proteome completeness |
| EggNOG | eggNOG-mapper | COG/NOG categories, KEGG pathways |
| CAZyme | dbCAN (DIAMOND + HMMER) | Carbohydrate-active enzymes |
| Secretome | SignalP/DeepSig + DeepTMHMM | Signal peptides, transmembrane topology |
| BGC clusters | antiSMASH | Secondary metabolite biosynthetic gene clusters |
| Proteases | MEROPS (DIAMOND) | Peptidase families and clans |
| tRNA | tRNAscan-SE | tRNA genes (eukaryotic, mitochondrial, general) |

---

## Genetic code support

For organisms with non-standard genetic codes:

```bash
myconote-cli annotate genes.gff3 --fasta genome.fa --genetic-code 12   # Candida CTG clade
```

| Code | Table | Organisms |
|------|-------|-----------|
| 1 | Standard | Most eukaryotes (default) |
| 12 | Alternative Yeast Nuclear | Candida CTG clade (CTG = Ser) |
| 3 | Yeast Mitochondrial | Yeast mitochondria |
| 4 | Mold Mitochondrial | Mold/protozoan mitochondria |
| 2 | Vertebrate Mitochondrial | Vertebrate mitochondria |

18 NCBI translation tables supported in total.

---

## Evidence weighting

Customize how the Evidence Modeler scores predictions with a TOML file:

```toml
# weights.toml
augustus = 10.0
snap = 3.0
protein = 25.0
est = 8.0
glimmerhmm = 2.0
genemark = 5.0
```

```bash
myconote-cli predict genome.fa --weights weights.toml
```

---

## Supported kingdoms

| Kingdom | Default Augustus model | BUSCO lineage | Intron range |
|---------|----------------------|---------------|--------------|
| `fungi` | saccharomyces_cerevisiae_S288C | fungi_odb10 | 40-2,000 bp |
| `plant` | arabidopsis | viridiplantae_odb10 | 40-50,000 bp |
| `animal` | human | metazoa_odb10 | 40-500,000 bp |
| `insect` | fly | insecta_odb10 | 40-50,000 bp |
| `protist` | toxoplasma | eukaryota_odb10 | 20-1,000 bp |

---

## Databases

| Database | Size | Used by |
|----------|------|---------|
| UniProt/Swiss-Prot (MMseqs2) | ~1 GB | `annotate` |
| Pfam-A HMM profiles | ~300 MB | `annotate` |
| BUSCO lineage datasets | ~500 MB | `annotate` |
| CAZyme (dbCAN DIAMOND) | ~200 MB | `annotate --cazyme` |
| MEROPS (DIAMOND) | ~100 MB | `annotate --merops` |
| **eggNOG-mapper** | **~50 GB** | `annotate --eggnog` |

Database versions are tracked with timestamps for reproducibility. Run `myconote-cli setup --check` to see download dates.

---

## Output formats

- **GFF3** -- gene models (IGV, JBrowse2, UCSC, Geneious compatible)
- **GenBank (.gbk)** -- Geneious, SnapGene, BioPython
- **NCBI .tbl + .sqn** -- GenBank submission-ready files
- **FASTA** -- protein and nucleotide sequences
- **TSV** -- functional annotation tables (R/Python analysis)
- **PNG/SVG** -- linear and circular genome maps
- **HTML** -- interactive JBrowse2 and UCSC browser views
- **Newick** -- IQ-TREE phylogenetic trees
- **JSON** -- reproducibility reports (tool versions, parameters, checksums)

---

## Output validation

myconote-cli automatically validates all outputs:

- **GFF3**: ID uniqueness, Parent reference consistency, coordinate ordering, feature hierarchy
- **Protein FASTA**: internal stop codons, invalid amino acid characters, minimum length
- **NCBI compliance**: duplicate locus_tags, orphan features, missing qualifiers

---

## Interactive tutorial

Learn myconote-cli step by step, like R's swirl:

```bash
myconote-cli learn          # list all 8 lessons
myconote-cli learn 1        # start lesson 1
myconote-cli learn --resume # pick up where you left off
```

8 lessons covering the full pipeline, from basics to NCBI submission. Progress is saved between sessions.

---

## Validated on real genomes

Tested end-to-end on *Brettanomyces bruxellensis* (12.9 Mb, 5,218 genes):

| Metric | Result |
|--------|--------|
| Genes predicted | 5,218 (Augustus) |
| Functionally annotated | 98.3% (with InterProScan) |
| GO terms | 93.1% of genes |
| Pfam domains | 85.6% of genes |
| NCBI validation | PASSED (0 errors) |

---

## Citation

If you use Myconote_CLI in your research, please cite:

> Narh-Madey, B. (2026). *Myconote_CLI: A modular genome annotation pipeline for eukaryotes.*
> Hittinger Lab, Laboratory of Genetics, University of Wisconsin-Madison.

A CITATION.cff file is included for automated citation tools.

Please also cite the underlying tools used in your analysis (Augustus, SNAP, GlimmerHMM, PASA, Trinity, eggNOG-mapper, BUSCO, antiSMASH, RepeatModeler, tRNAscan-SE, IQ-TREE, etc.).

---

## License

MIT License -- see [LICENSE](LICENSE) for details.

---

## Contact

**Benjamin Narh-Madey**
Hittinger Lab, Laboratory of Genetics - UW-Madison
narhmadey@wisc.edu
