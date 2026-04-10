# MycoNote-CLI Reproducibility Package

This document provides the exact commands, software versions, and input files needed to reproduce all results in the MycoNote-CLI manuscript.

**Manuscript:** Narh-Madey B et al. (2026). MycoNote-CLI: an integrated, validated, and reproducible eukaryotic genome annotation pipeline.

---

## 1. Software Versions

The following software versions were used for all reported benchmarks:

| Software | Version | Source |
|----------|---------|--------|
| MycoNote-CLI | v0.1.0 | https://github.com/K-nie/myconote-cli/releases/tag/v0.1.0 |
| Rust toolchain | 1.85.0 | rustup default stable |
| Augustus | 3.5.0 | bioconda channel |
| SNAP | 2013_11_29 | bioconda channel |
| HMMER | 3.4 | bioconda channel |
| MMseqs2 | 13.45111 | bioconda channel |
| BLAST+ | 2.15.0 | bioconda channel (for comparison runs) |
| miniprot | 0.13 | bioconda channel |
| minimap2 | 2.28 | bioconda channel |
| samtools | 1.20 | bioconda channel |
| BUSCO | 5.7.1 | bioconda channel |
| Trinity | 2.15.1 | bioconda channel |
| RepeatModeler2 | 2.0.5 | bioconda channel |
| RepeatMasker | 4.1.5 | bioconda channel |
| tRNAscan-SE | 2.0.12 | bioconda channel |
| EggNOG-mapper | 2.1.12 | bioconda channel |
| dbCAN | 4.1.4 | bioconda channel |

For comparison runs:
| Software | Version | Source |
|----------|---------|--------|
| funannotate | 1.8.17 | bioconda channel (planned) |
| MAKER | 3.01.04 | bioconda channel (planned) |
| BRAKER | 3.0.8 | bioconda channel (planned) |

---

## 2. Database Versions

| Database | Version / Date | Size |
|----------|---------------|------|
| UniProt/Swiss-Prot | Release 2026_04 (2026-04-01) | ~270 MB compressed |
| Pfam-A | Pfam 36.0 (2025-09-15) | ~300 MB compressed |
| BUSCO fungi_odb10 | 2024-01-08 | ~50 MB |
| dbCAN HMMs | V12 | ~50 MB |
| MEROPS | Release 12.5 | ~15 MB |

---

## 3. Hardware Configuration

All benchmarks reported in the manuscript were measured on:

- **CPU:** Apple M3 Pro (12-core: 6 performance + 6 efficiency cores)
- **Memory:** 18 GB unified memory
- **Storage:** Apple SSD (built-in)
- **OS:** macOS 15.2 (Sequoia)
- **Kernel:** Darwin 24.5.0
- **Power state:** AC power, no thermal throttling observed

For thread-count flags, we used `--threads 4` to match the count of performance cores typically allocated to a single process during testing.

---

## 4. Test Datasets

### 4.1 Brettanomyces bruxellensis

- Assembly: 12.9 Mb across 30 contigs
- Source: Hittinger Lab in-house assembly (will be deposited in GenBank prior to manuscript publication)
- File: `brettanomyces_test/01_sorted/genome.fas` in the MycoNote-CLI repository
- Already-sorted version: `brettanomyces_test/01_sorted/genome.fas`
- Already-masked version: `brettanomyces_test/01_sorted/genome_masked.fas`
- Pre-existing predictions: `brettanomyces_test/02_predict/consensus.gff3`

### 4.2 Candida tropicalis

- Assembly: 14.2 Mb across 24 contigs
- Source: NCBI accession (specific accession to be confirmed)
- File: `tests/data/candida_tropicalis.fas` in the MycoNote-CLI repository
- Pre-existing annotation: `tests/data/candida_tropicalis.final.gff3`

---

## 5. Exact Commands

### 5.1 Setup

```bash
# Install MycoNote-CLI v0.1.0
git clone https://github.com/K-nie/myconote-cli.git
cd myconote-cli
git checkout v0.1.0
cargo build --release
sudo cp target/release/myconote-cli /usr/local/bin/

# Install all external tools (one-time)
myconote-cli install --yes

# Download all databases (one-time)
myconote-cli setup
```

### 5.2 Brettanomyces bruxellensis Pipeline

```bash
cd ~/myconote_paper_runs/

# Step 1: Sort
myconote-cli sort brettanomyces.fa \
  --output sorted.fa \
  --min-length 500 \
  --rename-table rename.tsv

# Step 2: Mask
myconote-cli mask sorted.fa \
  --output masked.fa \
  --engine repeatmodeler \
  --threads 4

# Step 3: Predict (without RNA-seq for this benchmark)
myconote-cli predict masked.fa \
  --kingdom fungi \
  --locus-prefix BRET \
  --output predict_out \
  --threads 4

# Step 4: Annotate
myconote-cli annotate predict_out/consensus.gff3 \
  --fasta masked.fa \
  --output annotate_out \
  --kingdom fungi \
  --locus-prefix BRET \
  --no-busco \
  --threads 4

# Step 5: Run hmmsearch-based Pfam separately for timing
time myconote-cli annotate predict_out/consensus.gff3 \
  --fasta masked.fa \
  --output annotate_pfam \
  --kingdom fungi \
  --no-busco \
  --threads 4

# Step 6: InterProScan via EBI API
time myconote-cli annotate predict_out/consensus.gff3 \
  --fasta masked.fa \
  --output annotate_ipr \
  --kingdom fungi \
  --no-busco --no-pfam \
  --interproscan --email YOUR_EMAIL_HERE \
  --threads 4

# Step 7: Validate for NCBI submission
myconote-cli submit annotate_out/annotated.gff3 \
  --fasta masked.fa \
  --organism "Brettanomyces bruxellensis" \
  --validate-only
```

### 5.3 Candida tropicalis Pipeline (CTG Clade)

```bash
# Annotate with Table 12 (Alternative Yeast Nuclear Code)
myconote-cli annotate tests/data/candida_tropicalis.final.gff3 \
  --fasta tests/data/candida_tropicalis.fas \
  --output candida_ctg \
  --kingdom fungi \
  --locus-prefix CTRO \
  --genetic-code 12 \
  --no-busco --no-pfam \
  --threads 4

# Compare with Table 1 (Standard) for validation
myconote-cli annotate tests/data/candida_tropicalis.final.gff3 \
  --fasta tests/data/candida_tropicalis.fas \
  --output candida_std \
  --kingdom fungi \
  --locus-prefix CTRO \
  --genetic-code 1 \
  --no-busco --no-pfam \
  --threads 4

# Compare protein outputs
diff <(grep -v "^>" candida_ctg/proteins.fa | tr -d '\n') \
     <(grep -v "^>" candida_std/proteins.fa | tr -d '\n') \
     | wc -c
```

### 5.4 hmmscan vs hmmsearch Comparison

```bash
PROTEINS="annotate_pfam/proteins.fa"
PFAM_DB="$HOME/.myconote/dbs/pfam/Pfam-A.hmm"

# hmmsearch (MycoNote-CLI default)
time hmmsearch --domtblout hmmsearch_out.domtblout \
  --cpu 4 -E 1e-5 --domE 1e-5 --noali \
  $PFAM_DB $PROTEINS

# hmmscan (funannotate default)
time hmmscan --domtblout hmmscan_out.domtblout \
  --cpu 4 -E 1e-5 --domE 1e-5 --noali \
  $PFAM_DB $PROTEINS

# Compare hit counts
wc -l hmmsearch_out.domtblout hmmscan_out.domtblout
```

---

## 6. Planned Comparative Benchmark (For Final Manuscript)

**This section describes the comprehensive benchmark we plan to run before final manuscript submission. These results are not yet available; the v0.1.0 release contains only the two-genome demonstration.**

### 6.1 Reference Genomes (Planned)

We plan to benchmark on the following 8 reference genomes spanning four eukaryotic kingdoms, all with curated reference annotations:

| Organism | Kingdom | Source | Reference annotation |
|----------|---------|--------|---------------------|
| *Saccharomyces cerevisiae* S288C | Fungi | SGD | sacCer3 / R64-1-1 |
| *Aspergillus nidulans* FGSC A4 | Fungi | AspGD | NCBI RefSeq |
| *Cryptococcus neoformans* JEC21 | Fungi | NCBI RefSeq | NCBI RefSeq |
| *Candida albicans* SC5314 | Fungi (CTG) | CGD | CGD curated |
| *Arabidopsis thaliana* | Plant | TAIR | TAIR10 |
| *Drosophila melanogaster* | Insect | FlyBase | FlyBase r6 |
| *Caenorhabditis elegans* | Animal | WormBase | WormBase WS283 |
| *Plasmodium falciparum* 3D7 | Protist | PlasmoDB | PlasmoDB-66 |

### 6.2 Comparison Tools (Planned)

| Tool | Version | Purpose |
|------|---------|---------|
| MycoNote-CLI | v0.1.0 | This manuscript |
| funannotate | 1.8.17 | Established fungal pipeline |
| MAKER3 | 3.01.04 | Flexible multi-kingdom |
| BRAKER3 | 3.0.8 | RNA-seq based gold standard |

### 6.3 Accuracy Metrics (Planned)

For each tool on each genome:

- **Sensitivity**: True positives / (True positives + False negatives) at exon, transcript, and gene levels
- **Specificity (Precision)**: True positives / (True positives + False positives)
- **F1 score**: Harmonic mean of sensitivity and specificity
- **BUSCO completeness**: Versus appropriate lineage
- **Functional annotation rate**: Fraction of genes with at least one functional annotation
- **Annotation concordance**: Agreement between tool's GO term assignments and reference

### 6.4 Performance Metrics (Planned)

For each tool on each genome:

- **Wall-clock time**: Total pipeline runtime
- **Peak memory usage**: Maximum RSS during execution
- **CPU utilization**: User + system time / wall time
- **Storage footprint**: Total intermediate and output file size
- **Number of replicates**: 3 independent runs per tool per genome

### 6.5 Statistical Analysis (Planned)

- Paired t-tests comparing MycoNote-CLI vs each alternative on each metric
- Effect sizes (Cohen's d) for practical significance
- Multiple testing correction (Benjamini-Hochberg)
- 95% confidence intervals for all reported differences

---

## 7. Reproducibility Reports

Each MycoNote-CLI run automatically generates a JSON workflow report containing:

- MycoNote-CLI version
- Pipeline steps executed
- External tool versions detected at runtime
- Database paths and download dates
- Input file paths and sizes
- Command-line parameters
- Runtime duration
- System information

These JSON files (one per run) are available in the `reproducibility_reports/` directory of the supplementary materials.

---

## 8. Known Limitations of This Reproducibility Package

We are committed to making MycoNote-CLI as reproducible as possible, but several limitations should be acknowledged:

1. **Bit-for-bit reproducibility is not guaranteed.** Augustus, MMseqs2, and several other tools use heuristics that may produce slightly different output on different hardware, with different thread counts, or in subsequent runs. Differences are typically <1% of predictions.

2. **Database snapshots are not archived.** We record the download date but do not archive the database files themselves. Re-running with newer database versions will produce different results. For long-term reproducibility, users should archive their database files.

3. **External tool versions are detected, not pinned.** The reproducibility report records the version at runtime but does not enforce them in future runs. The Docker container provides better version pinning.

4. **Hardware effects on floating-point operations.** SIMD instructions used by HMMER and similar tools may produce slightly different floating-point results on different CPU architectures.

5. **Random number seeds are not currently controlled.** Operations using stochastic methods (e.g., bootstrapping in IQ-TREE) are not seeded by MycoNote-CLI.

For maximum reproducibility, we recommend:
1. Using the Docker container with the v0.1.0 tag
2. Archiving your database files (not just version metadata)
3. Documenting your hardware in publications
4. Running multiple replicates and reporting variability

---

## 9. Contact for Reproducibility Issues

If you encounter difficulties reproducing any result in the manuscript, please:

1. File an issue at https://github.com/K-nie/myconote-cli/issues with the "reproducibility" label
2. Include your reproducibility JSON report
3. Include the output of `myconote-cli check`
4. Email the corresponding author at cthittinger@wisc.edu

We commit to addressing reproducibility issues within 2 weeks.

---

*Reproducibility Package v1.0 prepared April 2026 in support of Narh-Madey et al. (2026).*
