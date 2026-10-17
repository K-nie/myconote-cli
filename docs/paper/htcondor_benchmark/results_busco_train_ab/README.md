# A/B: Augustus training modes — stock vs self vs BUSCO-seeded (v0.9.0)

Reproducibility + methods record for the A/B that tests whether MycoNote v0.9.0's
new **BUSCO-seeded Augustus training** (`--augustus-training busco`, Feature #1)
raises gene-prediction accuracy over the stock pretrained model and over the
circular self-training mode — AND validates that the BUSCO-training path actually
runs end-to-end on real fungal genomes (it had only been unit-tested on synthetic
fixtures). Written against commit `2f06d16` (v0.9.0).

## Question

Three questions, per arm, across the 14-genome fungal benchmark:
1. Does BUSCO-seeded training (funannotate's accuracy mechanism: train Augustus on
   the genome's own Complete single-copy BUSCO genes) beat the stock pretrained
   Augustus model on gene-level F1?
2. Does it beat circular self-training (`--augustus-training self`)?
3. **Does the BUSCO path actually train per genome, or silently fall back to
   stock?** The code fails safe to stock on any BUSCO failure, printing a loud
   multi-line banner (`BUSCO-SEEDED TRAINING DID NOT RUN`). A per-genome
   trained-vs-fell-back table is the core deliverable; if it falls back
   everywhere, the busco arm == stock arm and the feature is not actually tested.

Baseline to beat (stock, 14-genome mean, from the matched benchmark):
**gene F1 loose 0.899 / strict 0.609.**

## Design — 3 arms × 14 genomes, predict + score only

Identical inputs across arms; the ONLY variable is `--augustus-training`.
- **stock**: default pretrained species model (control).
- **self**: `--augustus-training self` (trains on MycoNote's own first-pass calls).
- **busco**: `--augustus-training busco` (BUSCO genome-mode → Complete single-copy
  genes → Augustus etraining → predict).

`sort`/`mask` are NOT re-run or re-measured: each job reuses the Stage-1 matched
rep1 `masked.fa` + `rename_table.tsv` verbatim, runs only `predict`, and scores
`predict_out/consensus.gff3` (not the annotated GFF — annotate does not change gene
structure, and skipping it removes the slow Pfam/hmmsearch step). GeneMark-ES and
SNAP are enabled identically in all three arms (same provisioning as
`jobs/run_myconote.sh`), so the consensus differs across arms ONLY in its Augustus
model.

## BUSCO lineage — IMPORTANT finding

There is **no `--busco-lineage` / `--busco-db` CLI flag**. The BUSCO-training path
resolves the lineage from `--kingdom` via `Kingdom::busco_lineage()`
(`src/predict/kingdom.rs`), which hardcodes **`fungi_odb10`** for `--kingdom fungi`
— for ALL 14 genomes. The clade-specific `busco_lineage` column in
`configs/genomes.tsv` (saccharomycetes_odb10, sordariomycetes_odb10, …) is NOT
consumed by this feature and is recorded in the per-job log as informational only.
`fungi_odb10` (758 genes, 2024-01-08) is a valid broad fungal set for every genome
in the panel; funannotate likewise trains from a broad fungal BUSCO set.

The training call (`src/predict/busco_train.rs::run_busco_genome`) runs
`busco -m genome --augustus --augustus_species <clade-model> -l fungi_odb10
--offline` first; on this cluster the offline attempt fails (no cwd-relative
`busco_downloads`) and BUSCO retries WITH download — execute nodes have network to
busco-data.ezlab.org (verified in an existing annotate run.log: fungi_odb10
downloaded in ~10 s). So the path self-provisions the lineage per job.

## Inputs

- Masked FASTA + rename table: `results_timing_stage1_matched/myconote/<g>/rep1/{masked.fa,rename_table.tsv}`
  (all 14 present; verified). Seqids are `scaffold_N`; `rename_table.tsv` maps
  `scaffold_N → original accession`; the scorer lifts predicted seqids back with it.
- Reference annotations: `data/<g>/reference.gff3` (all 14 present; verified).
- 14 genomes: sce cal ylp ani ncr cne klac kpha fgra afum umay ccin spom tmel.

## Commands

Harness: `jobs/run_busco_ab.sh <genome_id> <arm> <benchmark_dir> <results_dir>`.
Per job, predict (the only measured stage):

```
myconote-cli predict <masked.fa> \
    --kingdom <K> --species <AUGUSTUS_SPECIES> \
    --locus-prefix <PREFIX> --genetic-code <CODE> \
    --output <scratch>/predict_out --threads 8 \
    [--augustus-training self | --augustus-training busco]
```

`--kingdom`, `--genetic-code` from `configs/genomes.tsv`; `--locus-prefix` =
uppercased genome id; `<AUGUSTUS_SPECIES>` from the clade-model case block copied
verbatim from `jobs/run_myconote.sh`. `--genetic-code` is passed for spec
compliance but is a verified no-op on this predict path (the CLI forwards it only
to BRAKER or `--require-complete-orf`, neither used), so stock stays identical to
the published baseline.

Scoring (coding-only, gene/exon/nucleotide F1 per Eilbeck et al. 2009):

```
python3 scripts/compare_annotations.py <scratch>/predict_out/consensus.gff3 \
    data/<g>/reference.gff3 --rename-table <rep1>/rename_table.tsv \
    --label <arm>_<g> --output <out>/metrics.json
```

Metrics reported: `gene_level.loose.F1`, `gene_level.strict.F1`,
`predicted_summary.n_genes`. Strict = identical ordered CDS structure; loose =
same-strand locus overlap.

Submit (from the benchmark dir, env exported so getenv=true propagates it):

```
export BENCHMARK_DIR=$PWD RESULTS_DIR=$PWD/results_busco_train_ab/results
mkdir -p "$RESULTS_DIR/logs"
condor_submit jobs/busco_ab_stock_self.sub    # 28 jobs (stock+self)
condor_submit jobs/busco_ab_busco_smoke.sub   # 1 job  (kpha busco — validate path)
# after kpha busco logs "Trained species (BUSCO-seeded)":
condor_submit jobs/busco_ab_busco_rest.sub    # 13 jobs (busco, rest)
```

Resources per job: 8 CPU, 16 GB, 20 GB disk; pinned to `has_avx2 && Microarch
x86_64-v3` for cross-arm wall-time comparability. Each job copies
AUGUSTUS_CONFIG_PATH (48 MB) into its Condor scratch (race-free; shared config
left pristine). Heavy BUSCO/Augustus I/O runs on execute-node scratch; only
`consensus.gff3`, `metrics.json`, `run.log`, `time_predict.log`, and BUSCO
training evidence (`short_summary*`, `busco_full_table.tsv`,
`busco_training.gff3`) are copied to the shared results dir.

## Engagement verdict (per job → `train_verdict.txt`)

- busco TRAINED  : run.log has `Trained species (BUSCO-seeded):` (+ records the
  `Built BUSCO training set: N` training-gene count).
- busco FELL_BACK: run.log has `BUSCO-SEEDED TRAINING DID NOT RUN` (+ captured reason).
- self  TRAINED / FELL_BACK: `Trained species:` vs `SELF-TRAINING DID NOT RUN`.

## Tool / software versions (verified on cluster, 2026-10-16)

- MycoNote-CLI **v0.9.0** (commit `2f06d16`), built by Condor job 146869
  (`cargo install rc=0`, ExitCode 0), synced `~/.cargo/bin` → `~/.local/bin`
  (run scripts prepend `~/.local/bin`).
- BUSCO **5.7.1** (supports `--augustus`, `--augustus_species`); lineage
  fungi_odb10 (eukaryota, 2024-01-08, 758 genes).
- Augustus **3.3.3** (non-AVX2; runs full pool) + companion scripts
  (new_species.pl, etraining, gff2gbSmallDNA.pl) in the `myconote_augustus` env.
- BLAST+ (tblastn, makeblastdb), hmmsearch, SNAP, DIAMOND, minimap2 in the
  `myconote` env; GeneMark-ES (gmes_linux_64_4).
- Scorer: `scripts/compare_annotations.py` (stdlib Python 3).

## Outputs

- Per (arm, genome): `results_busco_train_ab/results/<arm>/<g>/{consensus.gff3,
  metrics.json, run.log, time_predict.log, train_verdict.txt}` and, for busco,
  `short_summary*.txt`, `busco_full_table.tsv`, `busco_training.gff3`.
- Aggregate + per-genome tables: filled in the Results section below after the run.

## Results

Run: Condor cluster 146876 (clean batch, all 30 jobs ExitCode 0), 2026-10-16/17.
A prior batch (146872/4/5) was discarded after an executable-overwrite race (the
wrapper was re-scp'd while jobs were launching) corrupted some jobs mid-read; the
wrapper was then sanitized to pure ASCII and redeployed with no jobs live, and the
whole batch rerun clean. busco was run on 2 genomes (kpha proxy-model, sce
exact-model); the other 12 busco jobs were deliberately NOT run once the fallback
cause was shown to be genome-independent (see below) — running them would only
reproduce the stock numbers and waste the expensive tmel/ccin/fgra busco passes.

### (a) Per-arm mean gene F1 + mean predicted n_genes

| arm   | n  | mean F1 loose | mean F1 strict | mean n_genes | mean predict wall (s) |
|-------|----|---------------|----------------|--------------|-----------------------|
| stock | 14 | **0.8992**    | **0.6095**     | 8413         | 2081                  |
| self  | 14 | 0.7729        | 0.5499         | 5770         | 2745                  |
| busco | 2  | 0.9228†       | 0.6944†        | 5233†        | 508†                  |

stock REPRODUCES the published baseline (0.899 / 0.609) — harness + stock arm
validated. †busco n=2 and both FELL BACK to stock, so these are stock values for
sce+kpha, NOT a BUSCO-trained result.

**Verdict: neither training lever beats stock. self-training HURTS; busco never
trained.**

- **self vs stock: self loses on loose F1 in 14/14 genomes** (mean Δloose
  **−0.1263**, Δstrict **−0.0596**; strict loses in 9/14, marginally improves in 5).
  Driven by under-prediction — self-trained models call far fewer genes (e.g. sce
  5776→2552, umay 6982→2748, tmel 15989→3246). All 14 self jobs genuinely trained
  (verdict TRAINED; verified), so this is a real trained-model regression, not a
  fallback artefact.

  Per-genome Δ(self−stock) loose: sce −0.341, cal −0.215, ylp −0.003, ani −0.042,
  ncr −0.103, cne −0.136, klac −0.123, kpha −0.054, fgra −0.023, afum −0.031,
  umay −0.362, ccin −0.066, spom −0.051, tmel −0.220.

### (b) busco arm — per-genome trained vs fell-back (THE validation result)

| genome | model type | verdict | reason |
|--------|-----------|---------|--------|
| kpha   | proxy (pichia_stipitis) | **FELL_BACK** | BUSCO genome-mode run failed (Augustus protein-extraction) |
| sce    | exact (S288C)           | **FELL_BACK** | BUSCO genome-mode run failed (Augustus protein-extraction) |
| (12 others) | — | NOT RUN (cause shown genome-independent) | — |

**busco training did NOT run on any genome tested.** The failure is genome- and
species-independent (confirmed on an exact-model genome, sce, ruling out the
pichia_stipitis proxy). Two integration blockers, in order of discovery:

1. **FIXED (harness):** BUSCO's BBTools runner (`name="bbtools"`, `cmd="stats.sh"`)
   looks up `stats.sh` on PATH; BBMap lives only in the `myconote_busco` conda env,
   which was not on the harness PATH, so busco genome-mode aborted with *"bbtools
   tool cannot be found"*. Fixed by appending `~/.conda/envs/myconote_busco/bin` at
   the TAIL of PATH (gives `stats.sh` + its java without shadowing Augustus 3.3.3 /
   tblastn / busco from the higher-priority envs).
2. **NOT fixable in the harness (integration/version gap):** after the PATH fix,
   BUSCO 5.7.1 genome-mode runs bbtools → makeblastdb → tblastn → Augustus
   (867 candidate jobs, using the seed species) → then fails at protein extraction:
   *"Augustus did not recognize any genes matching the dataset fungi_odb10 in the
   input file."* The pool's Augustus for the MycoNote training/prediction stack is
   **3.3.3**; **BUSCO 5.x requires Augustus ≥3.4** (the only ≥3.4 build, 3.5.0, is in
   the `myconote_busco` env whose `busco` import is itself broken, and using it for
   the FINAL prediction would confound the busco-vs-stock comparison with an Augustus
   version change). MycoNote's own parsing (`src/predict/busco_train.rs`) never runs —
   BUSCO fails before writing `full_table.tsv`/gff — so this is NOT a MycoNote parsing
   bug; it is an Augustus↔BUSCO version incompatibility in the cluster environment.

Because busco falls back to stock on every genome, **the busco arm == the stock
arm** (sce 0.9335 vs stock 0.9328; kpha 0.9121 vs stock 0.9137 — identical within
one-rep noise).

### (c) busco vs stock + compute cost

- **busco does NOT beat stock** — it cannot, because it never trained a BUSCO-seeded
  model. Δ(busco−stock) ≈ 0 (loose −0.0005, strict +0.0004, n=2 paired) = the stochastic
  noise of a stock-fallback rerun.
- **Compute cost:** a SUCCESSFUL busco-training wall cost could not be measured (it
  never succeeded). The FAILED BUSCO attempt added negligible wall — BUSCO aborts
  fast at the extraction step (~1–2 min), so busco predict wall ≈ stock (sce busco
  468 s vs sce stock 535 s; kpha busco 548 s vs kpha stock 538 s). A successful run
  would additionally incur the full BUSCO genome-mode + `etraining` time, unmeasured here.

### Fix options for the BUSCO path (handoff to tool-builder / env owner)

1. Provision a BUSCO stack built against Augustus ≥3.4 whose trained model is still
   usable by the 3.3.3 prediction path (or upgrade the MycoNote prediction Augustus
   to ≥3.4 and re-baseline stock/self on it).
2. Or change `busco_train.rs` to run the BUSCO subprocess with a BUSCO-compatible
   Augustus (≥3.4) independent of the Augustus used for the final prediction — the
   BUSCO internal Augustus only needs to produce a trainable gene set, so its version
   need not match the prediction Augustus and would not confound the comparison.
3. Interim: add a MycoNote pre-flight that checks `stats.sh` on PATH and the Augustus
   version BUSCO will use, and refuse `--augustus-training busco` with a clear message
   rather than running BUSCO for minutes only to fall back.
