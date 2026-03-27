---
title: "Genome Annotation with Myconote_CLI"
author: "Benjamin Narh-Madey"
affiliation: "Hittinger Lab, University of Wisconsin–Madison"
date: "2026"
---

# Genome Annotation with Myconote_CLI

## Lesson Overview

**Duration:** ~3 hours (self-paced) or one half-day workshop

**Software required:** Myconote_CLI, conda, a Unix terminal (macOS or Linux)

---

### Learning Objectives

By the end of this lesson you will be able to:

1. Explain the stages of a eukaryotic genome annotation pipeline
2. Install and verify Myconote_CLI and its external tool dependencies
3. Run each pipeline step — sort, mask, predict, and annotate — on a real genome
4. Interpret annotation statistics and quality metrics
5. Visualise gene models in JBrowse2 and export to standard formats

---

### Prerequisites

- Basic command-line navigation (`cd`, `ls`, `mkdir`)
- Familiarity with FASTA and GFF3 file formats
- A conda installation (Miniconda or Anaconda)

---

## Episode 1: Setting Up

**Time:** ~30 minutes

### Background

Genome annotation is the process of identifying the locations and functions of genes in a DNA sequence. A typical eukaryotic annotation pipeline involves four core stages:

1. **Repeat masking** — identifying repetitive elements so they do not confuse gene predictors
2. **Gene prediction** — using statistical models to find gene structures (exons, introns, UTRs)
3. **Evidence integration** — combining predictions from multiple tools using evidence weights
4. **Functional annotation** — assigning biological functions to predicted proteins

Myconote_CLI wraps all of these stages into a single tool with consistent interfaces and standard output formats.

### Installing Myconote_CLI

There are three ways to install Myconote_CLI. **Option C is recommended for most users.**

**Option A: One-shot installer**

Runs everything automatically — Rust, conda environment, tools, and databases:

```bash
git clone https://github.com/K-nie/myconote-cli
cd myconote-cli
bash install.sh
```

**Option B: Manual installation**

For users who want control over each step:

```bash
conda create -n myconote python=3.10
conda activate myconote
cargo build --release
myconote-cli install --yes
```

**Option C: Use a pre-built release (recommended)**

Download the pre-compiled binary directly — no Rust compiler needed:

```bash
# macOS (Apple Silicon)
curl -L https://github.com/K-nie/myconote-cli/releases/latest/download/myconote-cli-aarch64-apple-darwin \
    -o myconote-cli && chmod +x myconote-cli

# Linux (x86-64)
curl -L https://github.com/K-nie/myconote-cli/releases/latest/download/myconote-cli-x86_64-unknown-linux-gnu \
    -o myconote-cli && chmod +x myconote-cli
```

I recommend Option C because it gives you a clean starting point without requiring a Rust toolchain, and it is always the latest tested release.

### Verifying the installation

```bash
myconote-cli --version
```

You should see the banner and version number. Then run the tool checker:

```bash
myconote-cli check
```

This reports the install status of all 27 external tools. A typical new installation will show most tools as installed via conda, with SignalP and GeneMark flagged as requiring manual download.

> **Key point:** You do not need every tool installed to run the pipeline. Myconote will skip steps where the required tool is missing and tell you what to install.

---

### Challenge 1.1

Run `myconote-cli check` and identify:

1. How many tools are installed?
2. Which tools are listed as MISSING?
3. What command does the output suggest for installing missing tools?

---

## Episode 2: Preparing Your Genome

**Time:** ~20 minutes

### Why we sort contigs first

Most genome assemblers output contigs in arbitrary order with long, non-informative headers like `NODE_1_length_2474448_cov_85.2`. The `sort` command:

- Orders contigs from longest to shortest (annotation tools perform best on long contigs)
- Renames headers to clean sequential IDs (`scaffold_1`, `scaffold_2`, ...)
- Filters out very short contigs that add noise without contributing gene models
- Saves a rename table so you can trace back to original IDs

### Running `sort`

```bash
myconote-cli sort genome.fa \
    --prefix scaffold \
    --min-length 1000 \
    --rename-table rename_table.tsv \
    -o sorted.fa
```

The output includes a summary:

```
Sorting genome: genome.fa
  Contigs kept:   30
  Total assembly: 13.2 Mbp
  Longest contig: 2.47 Mbp
  Shortest kept:  1,024 bp
  Output:         sorted.fa
  Rename table:   rename_table.tsv
```

> **Key point:** Always run `sort` before masking. Consistent, clean headers prevent errors in downstream tools that parse sequence IDs from GFF3 and FASTA files together.

---

### Challenge 2.1

Sort the practice genome, then answer:

```bash
myconote-cli sort tests/data/candida_tropicalis.fas \
    --min-length 10000 \
    -o ct_sorted.fa
```

1. How many contigs passed the 10,000 bp filter?
2. What fraction of the total assembly length do they represent?
3. What would the first contig header look like in `ct_sorted.fa`?

---

## Episode 3: Repeat Masking

**Time:** ~30 minutes (runtime varies with genome size)

### Background

Eukaryotic genomes contain 20–80% repetitive elements (transposons, tandem repeats, satellites). Gene predictors trained on unmasked genomes frequently predict false genes inside repeat regions. Soft-masking converts repeat nucleotides to lowercase letters — gene predictors can see the sequence context but are trained to ignore lowercase regions.

Myconote supports two masking engines:

| Engine | Approach | When to use |
|--------|----------|-------------|
| `repeatmasker` | Uses a pre-built repeat library (Dfam) | Fast; works well for well-studied organisms |
| `repeatmodeler` | Builds a custom repeat library *de novo* | Slower but more accurate for novel genomes |

**Option A: RepeatMasker with Dfam library**

```bash
myconote-cli mask sorted.fa \
    --engine repeatmasker \
    --species fungi \
    --threads 8
```

**Option B: RepeatModeler + RepeatMasker (recommended for novel genomes)**

```bash
myconote-cli mask sorted.fa \
    --engine repeatmodeler \
    --threads 16
```

I recommend Option B for any genome not closely related to a reference assembly in Dfam, because it builds a repeat library tailored to your specific organism.

The output is `sorted.fa.masked` — identical to the input but with repeat regions in lowercase.

---

### Challenge 3.1

After masking, calculate the repeat content:

```bash
# Count lowercase (masked) bases
python3 -c "
seq = open('sorted.fa.masked').read()
seq = ''.join(seq.split()[1::2] if '>' in seq else [seq])
total = sum(1 for c in seq if c.isalpha())
masked = sum(1 for c in seq if c.islower())
print(f'Masked: {masked/total*100:.1f}%')
"
```

What percentage of your genome was masked? Is this typical for a fungal genome?

> **Hint:** Most ascomycete fungi have 2–10% repeat content; basidiomycetes can have up to 40%.

---

## Episode 4: Gene Prediction

**Time:** ~45 minutes

### How multi-tool prediction works

Myconote runs up to four ab initio gene predictors simultaneously and combines their predictions using **EvidenceModeler (EVM)**. Each predictor is assigned a weight based on its reliability for the target kingdom:

| Predictor | Approach | Default weight (fungi) |
|-----------|----------|------------------------|
| Augustus | Generalised HMM | 6 |
| SNAP | Hidden Markov Model | 2 |
| GlimmerHMM | Generalised HMM | 2 |
| GeneMark-ES | Self-training HMM | 4 |

EVM produces a consensus gene set that is generally more accurate than any single predictor alone.

### Running `predict`

```bash
myconote-cli predict sorted.fa.masked \
    --kingdom fungi \
    --augustus --snap --glimmerhmm \
    --threads 8 \
    --out predict_out/
```

Key options:

- `--kingdom` — sets the default Augustus species model, BUSCO lineage, and EVM weights
- `--species` — override the Augustus species model (default: `saccharomyces_cerevisiae_S288C`)
- `--hints` — provide a protein or EST hints GFF for hint-guided prediction (improves accuracy)

### Choosing an Augustus species model

The default species model for fungi is `saccharomyces_cerevisiae_S288C`. If your organism is more closely related to a different well-studied fungus, override the default:

```bash
# For a filamentous ascomycete
myconote-cli predict sorted.fa.masked \
    --kingdom fungi \
    --species aspergillus_fumigatus \
    --threads 8

# See all available models
myconote-cli species --grouped
```

---

### Challenge 4.1

Run prediction on the practice genome and examine the output:

```bash
myconote-cli predict tests/data/candida_tropicalis.fas.masked \
    --kingdom fungi --threads 4 --out ct_predict/
```

Then check the gene statistics:

```bash
myconote-cli stats ct_predict/evm.gff3
```

1. How many genes were predicted?
2. What is the mean gene length?
3. How does this compare to the reference annotation (6,290 genes)?

---

## Episode 5: Functional Annotation

**Time:** ~40 minutes

### What functional annotation does

Gene prediction gives you *where* genes are. Functional annotation tells you *what* they do by comparing predicted proteins against curated databases:

| Module | Database | Provides |
|--------|----------|----------|
| `--mmseqs` | Swiss-Prot / UniRef | Gene names, species of best hit |
| `--pfam` | Pfam-A | Protein domain families |
| `--eggnog` | EggNOG / COG | GO terms, KEGG pathways, COG categories |
| `--cazyme` | CAZy | Carbohydrate-active enzyme families |
| `--merops` | MEROPS | Protease families |
| `--busco` | BUSCO fungi_odb10 | Completeness assessment |
| `--antismash` | antiSMASH | Secondary metabolite gene clusters |

### Running `annotate`

**Option A: Core annotation only (fast, ~30 min)**

```bash
myconote-cli annotate sorted.fa.masked \
    --gff evm.gff3 \
    --pfam --busco \
    --threads 8 \
    --out annotation/
```

**Option B: Full annotation (recommended)**

```bash
myconote-cli annotate sorted.fa.masked \
    --gff evm.gff3 \
    --eggnog --pfam --cazyme --merops --busco --antismash \
    --threads 8 \
    --out annotation/
```

I recommend Option B for a publication-quality annotation because it provides GO terms, KEGG pathways, and secondary metabolite clusters that are expected in most fungal genome papers.

> **Note:** Option B requires the eggNOG database (~50 GB). If it is not yet downloaded, Myconote will skip the `--eggnog` step and warn you.

---

### Challenge 5.1

After annotation, examine the functional coverage:

```bash
myconote-cli stats annotation/final.gff3 --format json | python3 -c "
import json, sys
d = json.load(sys.stdin)
print('Total genes:', d['gene_count'])
print('Annotated (any function):', d.get('annotated_count', 'N/A'))
"
```

What percentage of predicted genes received at least one functional annotation? What does a low percentage suggest about the organism's novelty?

---

## Episode 6: Visualising Results

**Time:** ~20 minutes

### Genome plots

```bash
# Linear map of the largest scaffold
myconote-cli plot annotation/final.gff3 \
    --type linear \
    --region scaffold_1:1-500000 \
    --output scaffold1_map.png

# Circular whole-genome map
myconote-cli plot annotation/final.gff3 \
    --type circular \
    --output whole_genome.png
```

### Interactive browser

```bash
myconote-cli view annotation/final.gff3 \
    --fasta sorted.fa \
    --output genome_browser.html
```

Open `genome_browser.html` in any web browser for a full JBrowse2 interactive view — no internet connection required.

### Converting to other formats

```bash
# For R/Python analysis
myconote-cli convert annotation/final.gff3 --to table -o genes.tsv

# For Geneious / SnapGene
myconote-cli convert annotation/final.gff3 --to genbank \
    --fasta sorted.fa -o annotation.gbk

# Protein sequences for BLAST or OrthoFinder
myconote-cli convert annotation/final.gff3 --to protein -o proteins.faa
```

---

## Key Points

- Run pipeline steps **in order**: sort → mask → (train) → predict → update → annotate
- Use `--kingdom fungi` for fungal genomes; this sets correct defaults for all downstream steps
- `myconote-cli check` tells you exactly which tools are installed before you start
- All outputs use **standard formats** (GFF3, FASTA, GenBank) compatible with other tools
- The `stats` command gives you a quality summary at any point in the pipeline
- eggNOG functional annotation is optional but strongly recommended for publication
- Two tools (SignalP 6, GeneMark-ES) require manual download due to licensing

---

## Further Reading

- Augustus: [https://bioinf.uni-greifswald.de/augustus/](https://bioinf.uni-greifswald.de/augustus/)
- EvidenceModeler: [https://evidencemodeler.github.io/](https://evidencemodeler.github.io/)
- BUSCO: [https://busco.ezlab.org/](https://busco.ezlab.org/)
- antiSMASH: [https://antismash.secondarymetabolites.org/](https://antismash.secondarymetabolites.org/)
- eggNOG-mapper: [http://eggnog-mapper.embl.de/](http://eggnog-mapper.embl.de/)
- PASA: [https://github.com/PASApipeline/PASApipeline](https://github.com/PASApipeline/PASApipeline)

---

## Instructor Notes

**Episode 1 (Setup):** Allow extra time — conda environment creation and Rust compilation can take 5–10 minutes on first run. Pre-build the binary ahead of time for a classroom setting.

**Episode 3 (Masking):** RepeatModeler on a 13 Mbp fungal genome takes ~1–2 hours. For a workshop, pre-run masking and provide the `.masked` file so participants start from Episode 4.

**Episode 4 (Prediction):** Full EVM prediction on 13 Mbp takes ~30–60 minutes depending on threads. Same advice — provide pre-computed output for time-constrained workshops.

**Episode 5 (Annotation):** Without eggNOG (~50 GB), annotation runs in ~15 minutes. The `--pfam --busco` flags alone are sufficient for the challenges.
