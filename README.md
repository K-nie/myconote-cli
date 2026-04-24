# Myconote_CLI

**Blazing-fast fungal genome annotation pipeline** — as genome sequencing becomes cheaper and long-read assemblies become routine, the bottleneck in fungal genomics has shifted from sequencing to annotation. Existing pipelines are slow, narrowly scoped, and produce outputs that require extensive manual cleanup before submission. myconote-cli addresses this gap: a high-performance Rust CLI that takes a fungal genome assembly from raw contigs to NCBI-ready submission, integrating 11 annotation sources with built-in validation and reproducibility tracking. It runs on other eukaryotes too — the defaults and benchmarks are tuned for fungi.

> Developed by **Benjamin Narh-Madey** · Hittinger Lab, Laboratory of Genetics · UW-Madison
> narhmadey@wisc.edu

---

## Overview

Myconote_CLI takes a genome assembly (FASTA) through a structured, reproducible pipeline:

```
Sort → Mask → Train → Predict → Update → Annotate → Submit
```

Each step is a separate subcommand. All intermediate files use standard formats (GFF3, FASTA, GenBank) that hand off cleanly to **Proksee** (web circular maps), **IGV** (desktop browser), and **clinker** (cross-species synteny).

**New in v0.2.0** — `explain` (local-LLM result interpreter) and `batch` (multi-genome + HTCondor). See [CHANGELOG.md](CHANGELOG.md) for the full list.

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

Installs Miniconda (if needed), creates a conda environment, builds the binary, installs all ~30 external tools, and downloads annotation databases.

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
myconote-cli check          # reports status of all ~30 external tools
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
| `batch` | Annotate multiple genomes in one command (directory or sample sheet, HTCondor support, resume) |

---

## Analysis commands

| Command | Description |
|---------|-------------|
| `stats` | Gene counts, lengths, GC content, N50, isoform stats, with taxon-aware expected-range warnings |
| `quant` | RNA-seq expression quantification: fastp QC → salmon with decoy-aware index → wide count + TPM matrices + tximport-ready per-sample `quant.sf` + reproducibility bundle |
| `compare` | N-genome ortholog inference (OrthoFinder wrapper) — pan-genome summary + rooted species tree; tiered genome-count caps (5 fungi / 3 small plants / 2 large) |
| `convert` | Format conversion: GFF3 ↔ GTF / BED / GenBank / CDS / protein FASTA; FASTA ↔ FASTQ / PHYLIP / NEXUS; VCF conversions |
| `clean` | Validate and repair GFF3 annotation files |
| `fix` | Repair errors in GenBank (.gbk) files |

**Visualization hand-off** — MycoNote-CLI does not ship its own genome browser or synteny renderer. Use the `convert` command to produce standard outputs, then hand off to best-in-class external tools:

- **Proksee** (web, <https://proksee.ca/>) — circular genome maps from GenBank
- **IGV** (desktop) — interactive browsing of GFF3 + FASTA
- **clinker** (`pip install clinker`) — cross-species gene-cluster synteny from GenBank files

---

## AI-powered interpretation

| Command | Description |
|---------|-------------|
| `explain` | LLM-powered interpreter for any pipeline stage -- grounded findings, paper citations, next-command recommendations |

---

## Utility commands

| Command | Description |
|---------|-------------|
| `install` | Install missing tools via conda/mamba (~30 bioinformatics tools; `signalp6`, `table2asn`, and `gmes_petap.pl` require manual licence steps) |
| `check` | Check which external tools are installed and their versions |
| `setup` | Download and index annotation databases (Swiss-Prot, Pfam, EggNOG, BUSCO, dbCAN, MEROPS, Ollama, paper corpus) |
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

- **GFF3** -- gene models (load in IGV, JBrowse2, UCSC, Geneious)
- **GenBank (.gbk)** -- upload to Proksee or open in Geneious / SnapGene / BioPython; feed into clinker for synteny
- **NCBI .tbl + .sqn** -- GenBank submission-ready files
- **FASTA** -- protein and nucleotide sequences
- **TSV** -- functional annotation tables (R / Python analysis)
- **Newick** -- IQ-TREE phylogenetic trees
- **JSON** -- reproducibility reports (tool versions, parameters, checksums); `explain` reproducibility bundles in `explain_<stage>_<timestamp>/`

---

## Output validation

myconote-cli automatically validates all outputs:

- **GFF3**: ID uniqueness, Parent reference consistency, coordinate ordering, feature hierarchy
- **Protein FASTA**: internal stop codons, invalid amino acid characters, minimum length
- **NCBI compliance**: duplicate locus_tags, orphan features, missing qualifiers

---

## Explain: AI-powered result interpreter

Interpret any pipeline stage's output with grounded, citation-backed explanations. Runs locally via Ollama -- no data leaves your machine.

```bash
# Full interpretation (deterministic rules + local LLM)
myconote-cli explain predict

# Rules only, no LLM needed
myconote-cli explain predict --no-llm

# Analyze pasted output (auto-detects GFF3, FASTA, NCBI errors, logs)
echo "ERROR: SEQ_FEAT.NoStop" | myconote-cli explain --paste

# Verbose mode shows all findings with evidence
myconote-cli explain predict --verbose

# Dry-run: see the assembled prompt without calling the LLM
myconote-cli explain predict --dry-run
```

**How it works:**

1. **Deterministic rule engine** evaluates your stage output against TOML-defined rules (gene count ranges, masking thresholds, validation error patterns)
2. **BM25 retrieval** searches a built-in knowledge base and optional Q1 paper corpus for relevant context
3. **Local LLM** (via Ollama) produces a natural-language interpretation grounded in the retrieved context
4. **Citation validator** strips any claim the LLM cannot back with a retrieved source
5. **Command recommender** suggests 1-3 copy-pasteable next commands

**Smart model selection:** The tool auto-detects your system RAM and picks the most capable local model:

| RAM | Model | Quality |
|-----|-------|---------|
| 48+ GB | llama3.3:70b-instruct-q4_K_M | Best |
| 24+ GB | qwen2.5:32b-instruct-q4_K_M | Excellent |
| 16+ GB | mistral-small:22b | Strong |
| 12+ GB | qwen2.5:14b | Good |
| 8+ GB | llama3.1:8b | Baseline |

**Privacy first:** All inference runs locally. No data ever leaves your machine. Override the model with `--model <name>` or `MYCONOTE_CHAT_MODEL` env var.

**Scientific disclaimer:** All interpretations are suggestive, not definitive. Findings must be independently verified in the context of your research.

Setup:

```bash
myconote-cli setup ollama        # install Ollama + pull best model for your hardware
myconote-cli setup chat-corpus   # download Q1 open-access papers for grounded citations
```

---

## Batch: multi-genome annotation

Annotate multiple genomes in one command. Accepts a directory of FASTA files or a TSV sample sheet with per-genome settings.

```bash
# Annotate all FASTAs in a directory
myconote-cli batch genomes/ --kingdom fungi --threads 8

# Use a sample sheet for per-genome settings
myconote-cli batch samples.tsv --parallel 4

# Select specific stages
myconote-cli batch genomes/ --stages sort,mask,predict

# Resume after interruption
myconote-cli batch --resume batch_out/
```

**Sample sheet format** (TSV, header required):

```
name        fasta                kingdom    species          genetic_code    locus_prefix
isolate_A   /data/isolate_A.fa   fungi      saccharomyces    1               ISOA
isolate_B   /data/isolate_B.fa   fungi      auto             12              ISOB
isolate_C   /data/isolate_C.fa   plant      arabidopsis      1               ISOC
```

Only the `fasta` column is required; all others use defaults.

**Dashboard:** Auto-detects your environment. On interactive terminals, shows live progress bars per genome. On HPC batch jobs (no TTY), prints timestamped log lines suitable for `tail -f`.

**Resume:** A `status.json` file tracks per-genome/per-stage progress. If your run is interrupted (server reboot, SSH disconnect), `--resume batch_out/` picks up where it left off.

**HTCondor support:**

```bash
# Generate HTCondor submit files (does not run locally)
myconote-cli batch genomes/ --condor --condor-mem 64G --condor-cpus 16

# Then submit to the cluster
condor_submit batch_out/condor.sub

# Monitor
condor_q
tail -f batch_out/condor_logs/job_0.out
```

Generates `condor.sub`, `run_genome.sh` (per-job wrapper), and `condor_genomes.txt` (argument list). Each genome runs as a separate job. Works with shared filesystems.

**Next steps hint:** When 2+ genomes succeed, prints a concrete `compare` command recipe for downstream ortholog inference.

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
