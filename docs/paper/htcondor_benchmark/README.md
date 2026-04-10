# MycoNote-CLI HTCondor Benchmark Suite

A complete benchmarking package for comparing MycoNote-CLI against funannotate, MAKER, and BRAKER3 across 8 reference eukaryotic genomes on the GLBRC HTCondor pool.

This package implements the comprehensive comparative benchmark requested by Reviewer #2 of the MycoNote-CLI manuscript.

---

## Overview

**Tools compared**: 4 (MycoNote-CLI v0.1.0, funannotate v1.8.17, MAKER v3.01.04, BRAKER v3.0.8)

**Genomes**: 8 reference genomes spanning fungi, plants, insects, animals, and protists (Table 1)

**Replicates**: 3 per tool per genome

**Total jobs**: 96 (4 tools x 8 genomes x 3 replicates)

**Estimated total compute**: 500-1,000 CPU-hours

**Estimated wall time on HTCondor**: 12-24 hours (with parallel slot allocation)

**Output**: Comparison tables suitable for inclusion in the manuscript Results section, sensitivity/specificity metrics at gene/exon/nucleotide levels for each tool-genome combination, performance metrics (runtime, memory), and BUSCO completeness scores.

---

## Test Genomes

### Table 1. Reference genomes used in the benchmark.

| ID | Organism | Size (Mb) | Kingdom | Reference annotation source |
|----|----------|-----------|---------|----------------------------|
| 1 | *Saccharomyces cerevisiae* S288C | 12 | Fungi | SGD R64-1-1 |
| 2 | *Candida albicans* SC5314 | 14 | Fungi (CTG) | CGD A22 |
| 3 | *Aspergillus nidulans* FGSC A4 | 30 | Fungi | AspGD / FungiDB |
| 4 | *Cryptococcus neoformans* JEC21 | 19 | Fungi | NCBI RefSeq |
| 5 | *Arabidopsis thaliana* Col-0 | 135 | Plant | TAIR10 / Araport11 |
| 6 | *Drosophila melanogaster* | 144 | Insect | FlyBase r6.55 |
| 7 | *Caenorhabditis elegans* | 100 | Animal | WormBase WS292 |
| 8 | *Plasmodium falciparum* 3D7 | 23 | Protist | PlasmoDB-66 |

All eight have manually curated reference annotations that serve as gold standards for sensitivity/specificity calculations.

---

## Directory Structure

```
htcondor_benchmark/
  README.md                          # this file
  setup.sh                           # one-time setup: download genomes, install tools
  submit_all.sh                      # submit all 96 jobs to HTCondor
  collect_results.sh                 # aggregate results into tables
  scripts/
    compare_annotations.py           # compute sens/spec from GFF3 vs reference
    extract_protein_metrics.py       # functional annotation metrics
    measure_performance.sh           # wrap runs with /usr/bin/time -v
    aggregate_metrics.py             # combine all per-job results
    download_references.sh           # download genomes + annotations
  configs/
    genomes.tsv                      # genome metadata
    tools.tsv                        # tool versions
  jobs/
    myconote.sub                     # MycoNote-CLI HTCondor submit file
    funannotate.sub                  # funannotate submit file
    maker.sub                        # MAKER submit file
    braker.sub                       # BRAKER submit file
    run_myconote.sh                  # MycoNote-CLI wrapper
    run_funannotate.sh               # funannotate wrapper
    run_maker.sh                     # MAKER wrapper
    run_braker.sh                    # BRAKER wrapper
```

---

## Quick Start

### 1. Clone this benchmark to your HTCondor submit node

```bash
git clone https://github.com/K-nie/myconote-cli.git
cd myconote-cli/docs/paper/htcondor_benchmark
```

### 2. Configure for your environment

Edit `configs/htcondor.conf` to set:
- Your username on the cluster
- Your conda installation path
- Output directory (should be on shared storage like `/staging/`)
- Email for job notifications

### 3. One-time setup

```bash
bash setup.sh
```

This downloads all 8 reference genomes and their curated annotations (~5 GB), creates conda environments for the comparison tools, and verifies the setup.

### 4. Submit all 96 jobs

```bash
bash submit_all.sh
```

This submits all jobs to HTCondor. Jobs run in parallel as slots become available.

### 5. Monitor progress

```bash
condor_q $(whoami)
```

### 6. Collect results when complete

```bash
bash collect_results.sh
```

This aggregates per-job metrics into manuscript-ready tables in `results/`.

---

## Metrics Computed

### Per-job (one tool, one genome, one replicate)

**Gene-level accuracy** (vs gold-standard reference):
- True positives, false positives, false negatives
- Sensitivity (recall) at three matching stringencies (loose / moderate / strict)
- Specificity (precision) at three stringencies
- F1 score

**Exon-level accuracy**:
- Exon sensitivity and specificity
- Exact exon match (both boundaries)
- Internal exon match (boundaries from splice sites)

**Nucleotide-level accuracy**:
- Nucleotide sensitivity (coding nucleotides correctly predicted as coding)
- Nucleotide specificity (predicted coding nucleotides that are actually coding)

**Functional annotation quality**:
- Fraction of genes with product description
- Fraction with high-confidence Swiss-Prot hit (>50% identity, e<1e-50)
- Fraction with Pfam domain
- Fraction with GO term assignment
- BUSCO completeness against appropriate lineage (fungi_odb10, viridiplantae_odb10, etc.)

**Performance metrics** (from `/usr/bin/time -v`):
- Wall-clock time (seconds)
- CPU user time + system time
- Peak resident set size (memory)
- Major and minor page faults
- File system inputs and outputs

**Submission readiness**:
- GFF3 validity (no duplicate IDs, no orphan parents, valid coordinates)
- Protein FASTA validity (no internal stops)
- NCBI submit --validate-only pass/fail

### Per-tool (aggregated across all genomes and replicates)

- Mean sensitivity, specificity, F1 with 95% confidence intervals
- Mean runtime per genome size
- Mean peak memory per genome size
- BUSCO completeness distribution

### Statistical analyses

- Paired t-tests for each metric, MycoNote-CLI vs each alternative
- Cohen's d effect sizes
- Multiple testing correction (Benjamini-Hochberg)
- Wilcoxon signed-rank as a non-parametric backup

---

## Hardware Requirements per Job

| Tool | Cores | Memory | Disk | Wall time |
|------|-------|--------|------|-----------|
| MycoNote-CLI | 8 | 16 GB | 5 GB | 30-90 min |
| funannotate | 8 | 16 GB | 30 GB | 1-4 hours |
| MAKER | 8 | 32 GB | 20 GB | 2-12 hours |
| BRAKER | 16 | 64 GB | 30 GB | 2-8 hours |

For the largest genomes (Arabidopsis, Drosophila), runtimes can be substantially longer.

---

## Submission Strategy

The `submit_all.sh` script uses HTCondor's `queue` mechanism to submit jobs in priority order:

1. **First**: All MycoNote-CLI jobs (smallest, fastest)
2. **Second**: All funannotate jobs (well-tested, moderate runtime)
3. **Third**: All BRAKER jobs (longest, most memory)
4. **Fourth**: All MAKER jobs (highly variable runtime)

This ensures that if compute time is limited, the most informative results (MycoNote-CLI accuracy) are obtained first.

Each submit file uses HTCondor's `getenv = true` to inherit the conda environment, and `transfer_input_files` to bring the genome and reference annotation to the execute node.

---

## Expected Outcomes

**For the manuscript**: A comparison table showing each tool's mean ± SD for sensitivity, specificity, F1, runtime, and memory across 8 genomes, with statistical significance markers. This is the central result the reviewer asked for.

**For users**: Confidence in MycoNote-CLI's accuracy relative to established tools. If the benchmarks show MycoNote-CLI achieves comparable accuracy to funannotate/MAKER/BRAKER while running faster, that justifies adoption. If it shows MycoNote-CLI is less accurate, that informs honest scoping of where it should be used.

**For the field**: A reproducible benchmark protocol that other developers can use to evaluate new annotation tools.

---

## Reproducibility

All results are reproducible by:
1. Cloning this repository at a specific commit
2. Running `bash setup.sh` (downloads pinned versions)
3. Running `bash submit_all.sh`
4. The aggregate metrics in `results/` will match those in the manuscript within statistical noise

Each individual job records its environment in a JSON file (`*.reproducibility.json`).

---

## Contact

Questions about this benchmark:
- Benjamin Narh-Madey (narhmadey@wisc.edu)
- Hittinger Lab, Laboratory of Genetics, UW-Madison

Issues with the benchmark scripts:
- Open an issue at https://github.com/K-nie/myconote-cli/issues with the "benchmark" label
