# Myconote_CLI

**Genome Annotation Pipeline** — a fast, modular command-line tool for annotating fungal (and other eukaryotic) genomes from raw assembly to functional annotation.

> Developed by **Benjamin Narh-Madey** · Hittinger Lab · University of Wisconsin–Madison
> narhmadey@wisc.edu

---

## Overview

Myconote_CLI takes a genome assembly (FASTA) through a structured, reproducible pipeline:

```
Sort → Mask → Train → Predict → Update → Annotate
```

Each step is a separate subcommand, giving you full control over where to start, pause, or resume. All intermediate files use standard formats (GFF3, FASTA, GenBank) compatible with Geneious, IGV, JBrowse2, and UCSC Genome Browser.

---

## Installation

### Prerequisites

- **macOS** (Apple Silicon or Intel) or **Linux** (x86-64)
- [Miniconda or Anaconda](https://docs.conda.io/en/latest/miniconda.html)
- Rust ≥ 1.74 (installed automatically if missing — see below)
- ~20 GB disk space for core tools; ~75 GB if downloading eggNOG databases

### Option A: One-shot installer (recommended)

```bash
bash install.sh
```

This script handles everything in order:

1. Installs Miniconda if conda is not found
2. Creates a dedicated `myconote` conda environment
3. Installs or updates Rust (requires ≥ 1.74)
4. Builds the `myconote-cli` binary
5. Installs all 25 bioinformatics tools via conda/mamba
6. Downloads annotation databases (~2.5 GB)

I recommend Option A because it installs everything in one step and correctly handles Rust version requirements.

### Option B: Manual build

If you already have conda and Rust ≥ 1.74:

```bash
conda create -n myconote python=3.10
conda activate myconote
cargo build --release
sudo cp target/release/myconote-cli /usr/local/bin/
myconote-cli install --yes
```

### Option C: Build only (no tool install)

If you only want the binary and will install external tools yourself:

```bash
cargo build --release
# Binary is at: target/release/myconote-cli
```

---

## Verifying your installation

```bash
myconote-cli --version      # shows banner + version
myconote-cli check          # reports status of all 27 external tools
myconote-cli setup --list   # shows annotation database download status
```

---

## Quick start

```bash
# 1. Sort and rename contigs (longest first)
myconote-cli sort genome.fa -o sorted.fa

# 2. Soft-mask repeats
myconote-cli mask sorted.fa --engine repeatmodeler --threads 8

# 3. Train gene predictors with RNA-seq (optional but recommended)
myconote-cli train sorted.fa.masked --rna-bam rna.bam --prefix myorg

# 4. Predict genes
myconote-cli predict sorted.fa.masked --kingdom fungi --threads 8

# 5. Refine with RNA-seq evidence
myconote-cli update sorted.fa.masked --gff predictions.gff3 --bam rna.bam

# 6. Functionally annotate
myconote-cli annotate sorted.fa.masked --gff updated.gff3 \
    --eggnog --pfam --cazyme --busco --threads 8
```

---

## Pipeline commands (run in order)

| Command | Description |
|---------|-------------|
| `sort` | Sort contigs by length, rename headers, filter short scaffolds |
| `mask` | Identify and soft-mask repeats (RepeatModeler + RepeatMasker) |
| `train` | RNA-seq–guided training of Augustus and SNAP via Trinity + PASA |
| `predict` | Multi-tool gene prediction (Augustus + SNAP + GlimmerHMM + GeneMark + EVM) |
| `update` | Refine gene models with RNA-seq evidence (PASA UTR extension) |
| `annotate` | Functional annotation (MMseqs2, Pfam, eggNOG, CAZyme, MEROPS, BUSCO, antiSMASH) |

---

## Analysis commands

| Command | Description |
|---------|-------------|
| `stats` | Gene counts, lengths, GC content, isoform stats — JSON, CSV, or human-readable |
| `convert` | GFF3 → GTF / BED6 / BED12 / TSV / protein FASTA / GenBank |
| `clean` | Validate and repair GFF3 annotation files |
| `fix` | Repair errors in GenBank (.gbk) files |
| `plot` | Genome maps — linear PNG or circular PNG |
| `phylogeny` | Maximum-likelihood tree with IQ-TREE 2 |
| `compare` | Multi-genome comparison (MMseqs2 / BLAST / MUMmer) |
| `view` | Interactive genome browser — JBrowse2 or UCSC custom track |
| `synteny` | Ribbon diagram comparing two genomes |

---

## Utility commands

| Command | Description |
|---------|-------------|
| `install` | Install missing tools via conda/mamba |
| `check` | Check which external tools are installed and their versions |
| `setup` | Download annotation databases (Swiss-Prot, Pfam, eggNOG, BUSCO, CAZyme, MEROPS) |
| `remote` | Submit proteins to remote servers (Phobius, InterProScan) |
| `species` | List all Augustus species models |

---

## Example: annotating a fungal genome

```bash
# Full pipeline with RNA-seq support
myconote-cli sort brettanomyces.fa --prefix bret -o bret_sorted.fa

myconote-cli mask bret_sorted.fa \
    --engine repeatmodeler --threads 16

myconote-cli train bret_sorted.fa.masked \
    --rna-bam rna_sorted.bam --prefix bret --threads 16

myconote-cli predict bret_sorted.fa.masked \
    --kingdom fungi \
    --species saccharomyces_cerevisiae_S288C \
    --augustus --snap --glimmerhmm \
    --threads 16 \
    --out bret_predict/

myconote-cli update bret_sorted.fa.masked \
    --gff bret_predict/evm.gff3 \
    --bam rna_sorted.bam \
    --out bret_update/

myconote-cli annotate bret_sorted.fa.masked \
    --gff bret_update/updated.gff3 \
    --eggnog --pfam --cazyme --merops --busco --antismash \
    --threads 16 \
    --out bret_annotation/
```

---

## Databases

Myconote_CLI uses several annotation databases downloaded by `myconote-cli setup`:

| Database | Size | Used by |
|----------|------|---------|
| UniProt/Swiss-Prot (MMseqs2) | ~1 GB | `annotate` |
| Pfam-A HMM profiles | ~300 MB | `annotate --pfam` |
| BUSCO fungi_odb10 | ~50 MB | `annotate --busco` |
| CAZyme (DIAMOND) | ~200 MB | `annotate --cazyme` |
| MEROPS (DIAMOND) | ~100 MB | `annotate --merops` |
| **eggNOG-mapper** | **~50 GB** | `annotate --eggnog` |

> **Note:** The eggNOG database is not downloaded by default due to its size (~50 GB). Download it manually when needed:
> ```bash
> download_eggnog_data.py -y --data_dir ~/.eggnog_mapper/data
> ```

### Tools requiring manual installation

Two tools cannot be installed automatically due to licensing:

| Tool | Reason | Instructions |
|------|--------|-------------|
| **SignalP 6** | Academic license required | [https://services.healthtech.dtu.dk/services/SignalP-6.0/](https://services.healthtech.dtu.dk/services/SignalP-6.0/) |
| **GeneMark-ES** | License agreement required | [http://topaz.gatech.edu/GeneMark/](http://topaz.gatech.edu/GeneMark/) |

---

## Disk space summary

| Component | Space |
|-----------|-------|
| myconote-cli binary | ~5 MB |
| conda environment + tools | ~8–12 GB |
| Core databases (Swiss-Prot, Pfam, BUSCO, CAZyme, MEROPS) | ~2.5 GB |
| eggNOG database (optional) | ~50 GB |
| **Total (without eggNOG)** | **~15 GB** |
| **Total (with eggNOG)** | **~65 GB** |

---

## Output formats

All pipeline outputs use standard bioinformatics formats:

- **GFF3** — gene models compatible with IGV, JBrowse2, UCSC, Geneious
- **GenBank (.gbk)** — compatible with Geneious, SnapGene, BioPython
- **FASTA** — protein sequences for downstream analysis
- **TSV** — functional annotation tables for R/Python analysis
- **PNG** — linear and circular genome maps
- **HTML** — interactive JBrowse2 and UCSC browser views
- **Newick** — IQ-TREE phylogenetic trees

---

## Supported kingdoms

| Kingdom | Default Augustus model | BUSCO lineage |
|---------|----------------------|---------------|
| `fungi` | saccharomyces_cerevisiae_S288C | fungi_odb10 |
| `plant` | arabidopsis | viridiplantae_odb10 |
| `animal` | human | metazoa_odb10 |
| `insect` | fly | insecta_odb10 |
| `protist` | toxoplasma | eukaryota_odb10 |

---

## Citation

If you use Myconote_CLI in your research, please cite:

> Narh-Madey, B. (2026). *Myconote_CLI: A modular genome annotation pipeline for eukaryotes.*
> Hittinger Lab, University of Wisconsin–Madison.

Please also cite the underlying tools used in your analysis (Augustus, SNAP, GlimmerHMM, PASA, Trinity, eggNOG-mapper, BUSCO, antiSMASH, RepeatModeler, etc.).

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

**Benjamin Narh-Madey**
Hittinger Lab · University of Wisconsin–Madison
narhmadey@wisc.edu
