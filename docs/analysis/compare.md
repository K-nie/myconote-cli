# myconote compare

N-genome ortholog inference and pan-genome analysis. Wraps [OrthoFinder](https://github.com/davidemms/OrthoFinder) (Emms & Kelly 2019, *Genome Biology*) — the research-grade standard — and surfaces its outputs as MycoNote-native files.

```bash
myconote-cli compare <g1.gff3> <g1.fa> <g2.gff3> <g2.fa> [...] [options]
```

## What you get

- **`ortholog_table.tsv`** — one row per orthogroup, columns = genomes. Each cell lists the gene IDs from that genome that fell into the cluster (empty = not present). Category column classifies the row as `core`, `core_single_copy`, `accessory`, or `singleton`.
- **`pangenome_summary.tsv`** — Tettelin-style bins: `core` / `core_single_copy` / `soft_core` / `shell` / `cloud` / `singletons`, with counts and fractions.
- **`species_tree.nwk`** — rooted species tree from STAG/STRIDE (produced by OrthoFinder).
- **`compare_report.txt`** — human-readable summary.
- **Untouched OrthoFinder output** preserved under `compare_out/proteins/OrthoFinder/Results_<date>/` — includes gene trees, duplication events, multiple sequence alignments (with `--msa`), and per-species single-copy ortholog FASTAs.

## Genome-count caps

Scope is capped by genome size because OrthoFinder's all-vs-all similarity search scales with `total_proteins²`. Above a certain total, commodity hardware swaps to disk and the run becomes unreliable.

| Tier | Proteins per genome | Typical examples | Max genomes |
|---|---|---|---|
| **Small** | ≤ 15,000 | fungi, most protists | **5** |
| **Medium** | 15,001 – 30,000 | small/mid plants, many invertebrates | **3** |
| **Large** | > 30,000 | major crops, vertebrates | **2** |

The cap is set by the **largest** input: one plant genome mixed in with four yeasts triggers the plant tier. Compare detects this from the actual protein counts after primary-transcript filtering, so you get an accurate tier rather than a guess from the file name.

Over-cap input is rejected with an actionable error message pointing you at:

1. Reducing the input count
2. Running OrthoFinder directly on HPC
3. `--force-cap <n>` (advanced escape hatch; not shown in `--help`)

## Quick examples

```bash
# Three closely related yeasts
myconote-cli compare \
  sc.gff3 sc.fa \
  sp.gff3 sp.fa \
  cg.gff3 cg.fa

# Two Arabidopsis accessions with all cores
myconote-cli compare \
  col0.gff3 col0.fa \
  cvi.gff3 cvi.fa \
  --threads 16 --msa

# Single-copy orthologs are emitted to compare_out/proteins/OrthoFinder/Results_*/
# Single_Copy_Orthologue_Sequences/ — ready to feed into external ML-tree
# builders (IQ-TREE, RAxML-NG) after alignment with MAFFT.
```

## Options

| Flag | Default | Description |
|---|---|---|
| `--output <dir>` / `-o` | `compare_out` | Output directory |
| `--threads <n>` / `-t` | all cores | Threads for OrthoFinder |
| `--sensitive` | on | Use `diamond_ultra_sens` — higher accuracy, ~2× slower |
| `--fast` | — | Use default `diamond` — faster, less accurate |
| `--msa` | off | Multiple-sequence-alignment tree refinement (2–3× slower) |
| `--genetic-code <n>` | 1 | NCBI translation table (e.g. 12 for *Candida* CTG clade) |
| `--soft-core <frac>` | 0.95 | Soft-core presence threshold |
| `--cloud <frac>` | 0.15 | Cloud upper-bound threshold |

## Primary-transcript filtering

Compare keeps **one protein per gene** — the longest-CDS isoform — before running OrthoFinder. Alternative-splice isoforms cluster into the same orthogroup anyway, so feeding all of them inflates the all-vs-all search without adding information. Plant GFF3s with 5–8 isoforms per gene see the biggest speed-up.

This is always on. Users who truly need all isoforms can pre-extract proteins themselves and run OrthoFinder directly.

Pseudogenes and genes with internal stop codons are also dropped — they add noise to orthology inference and aren't useful tree-building material.

## Categories in `ortholog_table.tsv`

| Category | Meaning |
|---|---|
| `core_single_copy` | Present in every input genome with **exactly one** copy each. These are the rows used for species-tree inference. |
| `core` | Present in every input genome, but at least one has a paralog pair (duplications). |
| `accessory` | Present in some but not all genomes (2 ≤ n < N). |
| `singleton` | Present in only one genome. |

`core_single_copy` is usually the answer to "which genes should I align to build a species tree?" — these are automatically separated by OrthoFinder into the `Single_Copy_Orthologue_Sequences/` directory.

## Requirements

```bash
conda install -c bioconda orthofinder
```

OrthoFinder itself brings in DIAMOND, MAFFT (if `--msa`), FastTree, and the required Python stack. If you already use OrthoFinder in other workflows, no extra installation is needed.

## Typical runtimes

Measured on commodity hardware (8-core M3 laptop, 16 GB RAM):

| Input | Runtime |
|---|---|
| 5 yeasts | ~10 min |
| 3 Arabidopsis-class dicots | ~45 min |
| 2 maize-class monocots | ~3 h |
| 3 fungi with `--msa` | ~30 min |

HPC nodes with 32+ cores scale roughly linearly for the DIAMOND step.

## Related commands

- For visual pairwise or multi-way synteny, convert each genome to GenBank (`convert --to genbank`) and load the `.gbk` files into **clinker** (`pip install clinker`).
- For a species tree from compare's `Single_Copy_Orthologue_Sequences/`, align with MAFFT and run IQ-TREE or RAxML-NG externally.
- [`annotate`](../pipeline/annotate.md) — produces the GFF3 + FASTA pairs that compare consumes.
