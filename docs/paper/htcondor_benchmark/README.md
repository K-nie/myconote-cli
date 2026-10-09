# MycoNote-CLI HTCondor Benchmark Suite

A complete benchmarking package for comparing MycoNote-CLI against funannotate, MAKER, and BRAKER3 across six fungal reference genomes on the GLBRC HTCondor pool.

This package implements the comprehensive comparative benchmark requested by Reviewer #2 of the MycoNote-CLI manuscript.

---

## Overview

**Tools compared**: 4 (MycoNote-CLI v0.7.5, funannotate v1.8.17, MAKER v3.01.04, BRAKER v3.0.8)

**Genomes**: 6 fungal reference genomes spanning Saccharomycotina, Pezizomycotina, and Basidiomycota (Table 1). Non-fungal taxa are out of scope for this panel; the `--kingdom plant|animal|insect|protist` paths are experimental and not validated at this release.

**Replicates**: 3 per tool per genome

**Total jobs**: 72 (4 tools × 6 genomes × 3 replicates)

**Estimated total compute**: 500-1,000 CPU-hours

**Estimated wall time on HTCondor**: 12-24 hours (with parallel slot allocation)

**Output**: Comparison tables suitable for inclusion in the manuscript Results section, sensitivity/specificity metrics at gene/exon/nucleotide levels for each tool-genome combination, performance metrics (runtime, memory), and BUSCO completeness scores.

---

## Test Genomes

### Table 1. Reference genomes used in the benchmark.

| ID | Organism | Size (Mb) | Clade | NCBI RefSeq assembly |
|----|----------|-----------|-------|----------------------|
| sce | *Saccharomyces cerevisiae* S288C | 12 | Saccharomycotina (budding yeast) | GCF_000146045.2 (R64) |
| cal | *Candida albicans* SC5314 | 14 | Saccharomycotina (CTG clade) | GCF_000182965.3 (ASM18296v3) |
| ylp | *Yarrowia lipolytica* CLIB122 | 21 | Saccharomycotina (dimorphic) | GCF_000002525.2 (ASM252v1) |
| ani | *Aspergillus nidulans* FGSC A4 | 30 | Pezizomycotina (filamentous ascomycete) | GCF_000011425.1 (ASM1142v1) |
| ncr | *Neurospora crassa* OR74A | 40 | Pezizomycotina (filamentous ascomycete) | GCF_000182925.2 (NC12) |
| cne | *Cryptococcus neoformans* JEC21 | 19 | Basidiomycota (human pathogen) | GCF_000091045.1 (ASM9104v1) |

All six assemblies carry manually curated NCBI RefSeq annotations (which mirror the canonical community resources — SGD, CGD, AspGD, FungiDB, JGI MycoCosm — but through a single stable download endpoint). These serve as gold standards for sensitivity/specificity calculations.

The panel spans the three subphyla that contain almost all fungal-annotation targets: Saccharomycotina (budding yeasts), Pezizomycotina (filamentous ascomycetes), and Basidiomycota (mushrooms, pathogens). Non-fungal genomes are deliberately excluded — they're out of scope for the tool's current validated use case.

---

## Directory Structure

```
htcondor_benchmark/
  README.md                          # this file
  setup.sh                           # one-time setup: download genomes, install tools
  submit_all.sh                      # submit all 72 jobs to HTCondor
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

This downloads all six fungal reference genomes and their curated annotations (~5 GB), creates conda environments for the comparison tools, provisions the MycoNote-CLI reference databases (`myconote-cli setup`), and verifies the setup.

The MycoNote-CLI database step downloads Swiss-Prot, Pfam, BUSCO, EggNog, dbCAN (CAZyme), and MEROPS into `MYCONOTE_DB_DIR` (default `~/.myconote/dbs`), plus the Augustus fungal species bundle into `~/.myconote/augustus_config`. These download once and are reused across all 18 MycoNote replicates. Before submitting, export both locations so the jobs (which inherit the shell via `getenv = true`) can find them:

```bash
source configs/htcondor.conf   # sets MYCONOTE_DB_DIR and AUGUSTUS_CONFIG_PATH
```

On a multi-node pool, point `MYCONOTE_DB_DIR` and `AUGUSTUS_CONFIG_PATH` at shared storage (e.g. `/staging/$USER/...`) so every execute node sees one copy rather than re-downloading per node.

### 4. Submit all 72 jobs

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
- BUSCO completeness against the appropriate fungal lineage (saccharomycetes_odb10, eurotiomycetes_odb10, sordariomycetes_odb10, or tremellomycetes_odb10, chosen per genome)

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

## MycoNote-CLI Annotation Workflow

To keep the comparison fair, MycoNote-CLI's `annotate` step runs the full functional-annotation suite that matches funannotate's default workflow rather than a stripped-down configuration:

- **Swiss-Prot** (MMseqs2) and **Pfam** (HMM) — always on; product names and protein domains.
- **`--eggnog`** — EggNog-mapper COG/NOG categories and GO terms.
- **`--cazyme`** — carbohydrate-active enzyme families (dbCAN / DIAMOND).
- **`--merops`** — protease families (MEROPS / DIAMOND).

All four resolve their databases from `--db-dir` (`$MYCONOTE_DB_DIR`). Gene prediction (`predict`) reads Augustus fungal species models from `AUGUSTUS_CONFIG_PATH`.

Two sources that funannotate can also run are **deliberately excluded** here because they require network access from execute nodes, which the GLBRC pool does not grant to jobs: **antiSMASH** (secondary-metabolite BGC clusters) and **InterProScan** (EBI REST API). Excluding them keeps every replicate deterministic and network-free. `--secretome` and `--trnascan` are likewise omitted to avoid depending on external tool binaries (DeepSig/DeepTMHMM, tRNAscan-SE) on the execute nodes.

---

## Hardware Requirements per Job

| Tool | Cores | Memory | Disk | Wall time |
|------|-------|--------|------|-----------|
| MycoNote-CLI | 8 | 16 GB | 5 GB | 30-90 min |
| funannotate | 8 | 16 GB | 30 GB | 1-4 hours |
| MAKER | 8 | 32 GB | 20 GB | 2-12 hours |
| BRAKER | 16 | 64 GB | 30 GB | 2-8 hours |

Runtimes are dominated by `annotate` (functional annotation against MMseqs2 / Pfam / InterProScan), not by gene prediction. The largest genomes in the panel are *Aspergillus nidulans* (30 Mb) and *Neurospora crassa* (40 Mb); the full fungal panel fits comfortably inside the time budgets above.

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

## GLBRC Execute-Node Notes (from the actual run)

These are the environment quirks the MycoNote arm hit on the GLBRC `scarcity`
pool and how `jobs/run_myconote.sh` handles each. They are recorded here because
they are not obvious from the tool docs and a reviewer re-running the benchmark
on a different cluster will need the same fixes.

- **Conda tool-env split.** `predict`/`annotate` shell out to external binaries
  by bare `PATH` lookup (augustus, snap, diamond, mmseqs, hmmsearch/hmmscan,
  minimap2, run_dbcan, emapper.py, busco). On this pool they live in two conda
  envs: the `myconote` env has working snap/diamond/hmmer/mmseqs/minimap2/
  run_dbcan/emapper.py/busco but a **broken augustus** (links the wrong boost);
  the `myconote_augustus` env has a working augustus 3.3.3. The wrapper prepends
  `myconote` env `bin`, then `myconote_augustus` `bin`, so each tool resolves to
  a working copy.

- **`conda activate` silently no-ops on execute nodes.** The base conda's
  `profile.d/conda.sh` lives under `/opt/bifxapps`, which is **not mounted on
  every execute node**. When activation no-ops, most tools vanish from `PATH`,
  `annotate` skips Swiss-Prot/Pfam/EggNog/CAZyme, and `predict` falls back to
  Augustus-only — producing a run that exits 0 but is ~0 % functionally
  annotated. The fix is an **unconditional env-`bin` prepend** after the
  activate attempt: the env `bin` is on shared cephfs and its binaries are
  RPATH-linked to `../lib`, so prepending `bin` directly works whether or not
  activation succeeded. Verified end-to-end: 97.9 % of genes carry a product,
  96.5 % a GO term, 94.7 % a Pfam domain; BUSCO 96.4 % (fungi_odb10) on *S.
  cerevisiae*.

- **Augustus config skeleton.** `predict` reads Augustus species models from a
  fixed `AUGUSTUS_CONFIG_PATH`, but Augustus also needs the standard config
  skeleton (`model/`, `extrinsic/`, `profile/`, `cgp/`) alongside `species/`.
  A `species/`-only directory fails immediately with
  `Could not find config file .../model/states_shadow.cfg`. Copy the skeleton
  from a working augustus env's `config/` into `AUGUSTUS_CONFIG_PATH`.

- **`+RequestRuntime` removed from all `.sub` files.** GLBRC machines advertise
  `Runtime = undefined`; HTCondor expands `+RequestRuntime` to
  `TARGET.Runtime >= RequestRuntime`, which no slot satisfies, so jobs sit idle
  forever. The requirement is dropped from `myconote.sub`, `funannotate.sub`,
  `maker.sub`, and `braker.sub`.

- **EggNog DB path.** `annotate --eggnog` resolves the emapper.py database from
  `--eggnog-db`, **not** from `--db-dir` (the tool never derives it). The
  wrapper points `--eggnog-db` at `$MYCONOTE_DB_DIR/eggnog` (eggNOG 5.0.2,
  emapper 2.1.15).

- **Seqid rename lift-over before scoring.** `sort` renames contigs to
  `scaffold_N` (longest first) for NCBI-clean output, so the predicted GFF3 is
  in the `scaffold_N` namespace while `reference.gff3` keeps original accessions
  (e.g. `NC_001133.9`). `compare_annotations.py` matches features only on
  identical seqids, so without a lift-over **every metric is 0**. The wrapper
  now captures `sort --rename-table` and passes it to
  `compare_annotations.py --rename-table`, which inverts it (`new_id →
  original_id`) to lift predicted seqids back to the reference namespace before
  comparison. On the *S. cerevisiae* smoke this took gene-level F1 from 0.000 to
  0.912 (loose), with nucleotide-level F1 0.985.

- **Per-genome Augustus species model.** `predict --kingdom fungi` defaults the
  Augustus species to `saccharomyces_cerevisiae_S288C` for *every* genome
  (`src/predict/kingdom.rs`). That budding-yeast model is near-intron-less and
  wrecks gene structure on intron-rich fungi — an early `cne` (*Cryptococcus*)
  run under the default scored BUSCO 14.9 % and strict gene F1 ~0 despite
  loose gene F1 0.87 (genes on the right loci, wrong exon/intron structure).
  funannotate/BRAKER train species-specific models, so a fair comparison must
  give MycoNote the matching model. The wrapper maps each genome to an installed
  clade-appropriate Augustus species via `predict --species` (strain-exact where
  available): sce→`saccharomyces_cerevisiae_S288C`, cal→`candida_albicans`,
  ylp→`yarrowia_lipolytica`, ani→`aspergillus_nidulans`, ncr→`neurospora_crassa`,
  cne→`cryptococcus_neoformans_neoformans_JEC21`.

- **EVM weighting is Augustus-dominant.** With SNAP and Augustus both feeding
  EVidenceModeler, the consensus is Augustus-driven (5,465 genes on *S.
  cerevisiae*, deterministic). SNAP contributes little to the final gene set —
  a MycoNote characteristic, not a benchmark artifact.

- **Timing instrumentation (Step 7) fixed.** `run_myconote.sh` records per-stage
  wall/CPU/RSS by parsing `/usr/bin/time -v` logs. The original parser matched
  Elapsed with a regex requiring two colons *and* fractional seconds, so it
  silently failed on both of GNU time's real forms — `m:ss.ss` (sort/mask/
  predict) and `h:mm:ss` with integer seconds (annotate). The result was a
  `total_wall_seconds` that dropped the `annotate` stage entirely (the dominant
  one) and read 0 or mirrored a single stage. The parser now splits the trailing
  token on `:` and sums 1/2/3 fields. `performance.json` also now reports
  `total_cpu_seconds` (user+system) beside `total_wall_seconds`: wall time swings
  2–3× between replicates of the same job on the contended `scarcity` pool (e.g.
  `cne` 50–148 min) even with `request_cpus = 8` giving dedicated cores, because
  the variance is node-level memory-bandwidth / I/O contention and EggNog DB
  load, not CPU oversubscription. CPU-time is steadier but not perfectly
  contention-invariant here (it tracks some real EggNog/hmmer thread-scaling
  work), so the manuscript reports **median wall over the three replicates with
  the range**, and peak RSS (stable at ~6.8–6.9 GB) as the memory figure.

---

## Results (complete four-tool run; MycoNote v0.7.4–0.7.6 era)

All four tools completed 18/18 jobs (6 genomes × 3 replicates). Accuracy follows
the Eilbeck et al. (2009) gene/nucleotide sensitivity–specificity–F1 scheme.
Strict gene F1 uses the **CDS-structure** definition (`genes_match_strict`:
identical ordered CDS set, same strand) — not the old outer-bound match, which
confounded UTR-annotation convention with gene-structure accuracy and is
superseded. **This comparison is apples-to-apples on the predictor set: all
four tools ran with GeneMark** (funannotate's first run skipped GeneMark for an
unset `$GENEMARK_PATH`; it was re-run with GeneMark provisioned so both full
pipelines carry the same ab-initio stack). Each tool used its clade-appropriate
Augustus species (MycoNote via explicit `--species`; funannotate via
`--busco_seed_species`). Replicates are tight (per-rep gene-count drift
~0.1–0.5 %); values below are the per-genome mean of 3 reps.

### Accuracy — gene-level loose F1 / strict F1, nucleotide F1 (myco / funannotate / MAKER / BRAKER)

| genome | loose gene F1 | strict gene F1 | nucleotide F1 |
|---|---|---|---|
| sce | 0.891 / **0.939** / 0.912 / 0.917 | 0.788 / 0.813 / **0.827** / 0.729 | **0.985** / 0.968 / 0.986 / 0.980 |
| cal | 0.925 / **0.956** / 0.952 / 0.958 | 0.803 / 0.810 / **0.839** / 0.709 | 0.986 / 0.978 / **0.989** / 0.990 |
| ylp | 0.908 / **0.947** / 0.913 / 0.918 | 0.726 / 0.672 / **0.734** / 0.728 | 0.972 / 0.959 / 0.974 / **0.977** |
| ani | 0.897 / 0.921 / 0.919 / **0.932** | 0.445 / 0.488 / 0.465 / **0.524** | 0.941 / 0.936 / 0.944 / **0.952** |
| ncr | 0.897 / 0.909 / 0.898 / **0.914** | 0.564 / 0.520 / 0.568 / **0.647** | 0.964 / 0.936 / 0.964 / **0.978** |
| cne | 0.921 / **0.949** / 0.936 / 0.940 | 0.526 / 0.490 / 0.541 / **0.650** | 0.945 / 0.947 / 0.947 / **0.961** |

**Honest reading.** MycoNote is **competitive, not the accuracy leader**:
- *Loose gene F1* — MycoNote is lowest/near-lowest on every genome; funannotate leads.
- *Strict gene F1* — a split with funannotate (each wins 3/6); MAKER leads the
  yeasts (sce/cal/ylp), BRAKER leads the intron-rich (ani/ncr/cne). MycoNote is mid-pack.
- *Nucleotide F1* — MycoNote edges funannotate on 5/6 and is top-tier overall.

### Speed and memory — median over 3 reps, 8 threads (BRAKER 16)

| tool | gene prediction | functional annotation | total | peak RSS |
|---|---|---|---|---|
| **MycoNote** | **18.9 min** | 114.8 min (EggNog)¹ | 146.6 min | 6.6 GB |
| funannotate | 74.5 min² | 8.0 min | 81.8 min² | 2.8 GB |
| MAKER | 30.4 min (whole run; prediction only) | — none | 30.4 min | 0.3 GB |
| BRAKER | ~94 min² (prediction only; log-span) | — none | ~94 min² | — |

**MycoNote has the fastest gene prediction by a wide margin** (18.9 min vs
30–94). Three caveats stated plainly:
1. ¹ MycoNote's 114.8-min functional annotation is the **pre-v0.7.8 EggNog**
   (`--dmnd_iterate yes`). A same-node A/B showed `--dmnd_iterate no` runs
   ~4.4× faster (3h22m→46m) with **identical recall** (9,736/10,655 proteins);
   v0.7.8 makes `no` the default, projecting annotation to ~10–28 min/genome.
   MAKER and BRAKER do **no** functional annotation, so their totals are
   prediction-only and not comparable to the full-pipeline totals of
   MycoNote/funannotate.
2. ² funannotate and BRAKER use Augustus 3.5.0, which is AVX2-compiled and
   **SIGILLs on non-AVX2 execute nodes**; their jobs were gated to AVX2-capable
   (newer, faster) nodes, so their wall-times are **not directly comparable** to
   MycoNote/MAKER, which ran the full pool. MycoNote's Augustus 3.3.3 is not
   AVX2-compiled and runs everywhere — a genuine portability advantage.
3. Peak RSS: MycoNote's 6.6 GB is the EggNog-mapper DB load; its prediction
   stage alone is light.

---

### Function-by-function capability (what each tool actually ran)

| function | MycoNote | funannotate | MAKER | BRAKER |
|---|---|---|---|---|
| repeat masking | ✓ SelfAlign (no DB) | ✓ tantan (simple) | ✗ off | ✗ raw |
| Augustus / SNAP / GeneMark | ✓ / ✓ / ✓(ES) | ✓ / ✓ / ✓(EP) | ✓ / ✗ / ✗ | ✓ / ✗ / ✓(EP) |
| GlimmerHMM | ✗ | ✓ | ✗ | ✗ |
| consensus | own agreement-weighted | EVM | single predictor | TSEBRA |
| protein homology / Pfam | ✓ MMseqs2 / ✓ | ✓ Diamond / ✓ | ✗ | ✗ |
| **EggNog orthology** | ✓ | ✗ (default) | ✗ | ✗ |
| CAZyme / MEROPS / BUSCO | ✓ / ✓ / ✓ | ✓ / ✓ / ✓ | ✗ | ✗ |
| NCBI submission prep | ✓ table2asn | ✓ (.tbl) | ✗ | ✗ |

MycoNote and funannotate are the only full predict+functional-annotation
pipelines; MAKER and BRAKER are gene-prediction only.

## Bottom line (honest positioning)

MycoNote is **not the most accurate** fungal annotator in this panel — it is
competitive (beats funannotate on nucleotide F1, splits strict, trails on
loose; mid-pack against MAKER/BRAKER). Its defensible, durable advantages are:

1. **Fastest gene prediction** (18.9 min median; 1.6–5× faster than the others).
2. **Most complete functional annotation** out of the box — the only tool here
   doing EggNog orthology by default (MAKER/BRAKER do none); after v0.7.8 this
   runs at competitive speed.
3. **Portability / low friction** — single Rust binary, DB-free default masking
   (no RepeatMasker/RepBase/Dfam), Augustus 3.3.3 that runs the full compute
   pool (funannotate/BRAKER need AVX2 nodes), modern `table2asn`.

### funannotate-weakness coverage (summary)

A systematic audit of funannotate's full issue tracker (1,093 issues) + docs
mapped ~22 distinct weakness themes against MycoNote. MycoNote cleanly
**resolves ~4** (Augustus-PPX path avoided, DB-free default masking, no
Python-legacy in core, no CodingQuarry failures), **partially mitigates ~8**
(install footprint, DB/BUSCO provisioning, EggNog speed, SignalP/Phobius
licensing, IPRScan→REST, submission, determinism, Augustus-training class),
**still shares ~6** (GeneMark license, Augustus binary, PASA lock, tRNAscan,
contig-rename, silent-drop design), is **unverified on ~2** (large-assembly
scaling, table2asn validation), and **introduces ~3 of its own** (the
S288C-for-all-fungi default, BUSCO offline-dataset fallback, IPRScan
network-only). Several of these are addressed by later versions: v0.7.9 adds a
clade-aware species default + loud guard (S288C), opt-in lossless contig naming,
and submit structural pre-validation. **MycoNote does not remove funannotate's
hard external-tool dependencies (GeneMark, Augustus, emapper, PASA) — it
repackages and speeds up around them.**

## Caveats for the manuscript

- Node-class asymmetry: funannotate/BRAKER ran AVX2-gated; wall-times are not
  directly comparable across arms (CPU-time is steadier). Report prediction
  speed and total runtime separately.
- MycoNote used explicit per-genome `--species`; the naive default (pre-v0.7.9)
  is S288C for all fungi and collapses on divergent fungi — state the species
  provenance.
- EggNog timing above is the pre-v0.7.8 slow path; cite the v0.7.8 speedup
  explicitly if quoting MycoNote's total runtime.
- BRAKER `performance.json` timing was unusable (all-zero); its runtime is a
  `braker.log` log-span estimate, not the `/usr/bin/time` figure.

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
