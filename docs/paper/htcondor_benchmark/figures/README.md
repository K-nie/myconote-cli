# Benchmark figures + backing tables — n=5 Stage-1 (production) dataset

Manuscript figure set and figure-backing tables for the four-tool fungal-annotation
benchmark (MycoNote vs funannotate / MAKER / BRAKER). This is manuscript Methods /
Results raw material. Generated 2026-10-16.

## Dataset (what this is)
Stage-1 n=5 IS the production dataset (decision 2026-10-16: no separate Phase-5
matrix). 14 genomes × 4 tools × 5 replicates = 280 headline runs (+70 MycoNote-full
runs for the eggNOG figure = 350 total). All five condor clusters finished 70/70
with ExitCode 0 (verified via condor_history).

Collapse rules (per BENCHMARK_PLAN §1): **accuracy → per-genome MEAN** (near-deterministic),
**timing → per-genome MEDIAN** (HPC-noise variance). All figures plot the collapsed values.

## Data sources (verified before use)
| tool | accuracy + timing source | cluster | notes |
|------|--------------------------|---------|-------|
| MycoNote (HEADLINE) | results_timing_stage1_matched/myconote/ | 146854 | v0.7.16, **no eggNOG** (matched source set) |
| MycoNote (full, REF) | results_timing_stage1/myconote/ | 146848 | v0.7.15, **with eggNOG**; fig12 ONLY; ran at peak load → raw totals contention-confounded |
| funannotate | results_timing_stage1/funannotate/ | 146849 | 1.8.17 |
| MAKER | results_timing_stage1/maker/ | 146850 | 3.01.04; totals only, **no stage split**, RSS ~230 MB |
| BRAKER | results_timing_stage1/braker/ (accuracy) + condor_history 146852 (timing) | 146852 | 3.0.8; its performance.json is all-zero → wall = RemoteWallClockTime, CPU = RemoteUserCpu+RemoteSysCpu, peak RSS = MemoryUsage. All 70 NumJobStarts=1 (clean single-run wall, no restart overcount). |

## Timing fairness (how to read the speed figures)
All timing arms were pinned to one microarch generation, **x86_64-v3**, because strict
same-node pinning is infeasible on this shared, heavily-loaded pool (see
PHASE4_5_RUN_STATUS.md). Wall-clock is the headline; **CPU-seconds (fig5) is the
co-tenant-contention-robust check**. True single-tenant isolation is impossible on a
partitionable-slot shared pool regardless of pinning.

## Figures (each ships .pdf + 300 dpi .png; Nature style: lowercase bold panel letters,
## sans-serif, NO on-figure title, Okabe-Ito tool colours — MycoNote #0072B2,
## funannotate #E69F00, MAKER #009E73, BRAKER #CC79A7 — consistent across all figures)

- **fig1_accuracy_by_genome** — **Per-genome annotation accuracy across four tools.** Grouped bars, 14 genomes × 4 tools, four panels: (a) gene F1 loose, (b) gene F1 strict, (c) exon F1 exact, (d) nucleotide F1. Backing: `accuracy.tsv`.
- **fig2_accuracy_distribution** — **Distribution of each accuracy metric across the 14-genome panel, per tool.** Boxplots (spread across genomes), four panels a–d (one metric each). Backing: `accuracy.tsv`.
- **fig3_runtime_by_genome** — **Total wall-clock time per genome and tool.** Grouped bars, log y. Backing: `timing_wall.tsv`.
- **fig4_runtime_decomposition** — **Where annotation time is spent: predict vs annotate stage wall-clock.** Stacked (predict solid + annotate pale) for MycoNote and funannotate; MAKER and BRAKER shown as hatched totals (no stage instrumentation). Log y; note stacking on a log axis is visually approximate. Backing: `timing_wall.tsv`.
- **fig5_cpu_seconds** — **Total CPU-seconds per genome and tool (contention-robust speed axis).** Grouped bars, log y. Distinct from fig3 (CPU-time is insensitive to co-tenant wall-clock contention). Backing: `timing_cpu.tsv`.
- **fig6_peak_memory** — **Peak memory footprint and its stage origin.** (a) peak RSS per genome × tool, log y, all four tools; (b) predict vs annotate peak RSS (grouped) for the two tools with stage instrumentation — MycoNote is annotate-memory-bound (~2.5 GB flat), funannotate is predict-memory-bound (up to ~8.7 GB). Backing: `rss.tsv`.
- **fig7_myconote_annotate_sources** — **MycoNote annotation cost broken out by evidence source.** Stacked bars per genome (mmseqs, pfam, busco, cazyme, merops); Pfam/hmmsearch dominates. Backing: `annotate_per_source.tsv` (variant=matched).
- **fig8_scaling** — **Runtime vs genome size, with per-tool power-law fit.** Log–log scatter (one point per genome), fit line per tool; tmel (125 Mb) is the large point. Backing: `scaling.tsv`.
- **fig9_accuracy_speed_tradeoff** — **Accuracy–speed trade-off, one marker per tool (Pareto view).** Mean gene F1 (loose) vs median wall-min across genomes. Backing: `accuracy.tsv` + `timing_wall.tsv`.
- **fig10_accuracy_by_clade** — **Accuracy stratified by clade (descriptive stratum, not a tested factor).** Grouped bars over Saccharomycotina / other-Ascomycota / Basidiomycota; (a) gene F1 loose, (b) gene F1 strict. Backing: `clade_strata.tsv`.
- **fig11_gene_count_calibration** — **Predicted vs reference gene count per tool.** Scatter (56 points) with y = x line; above the line = over-prediction. Backing: `gene_counts.tsv`.
- **fig12_eggnog_cost** — **Per-source annotate time, MycoNote matched (headline) vs full (with eggNOG).** Grouped bars per source, log y. eggNOG is a genuine added source (~5800 s median); the small shared-source differences are peak-load contention noise — **raw matched-vs-full totals are contention-confounded and are NOT reported as a speedup.** Backing: `annotate_per_source.tsv`.

## Backing tables (TSV, one header row each)
- `accuracy.tsv` (56) — genome, clade, tool, gene_f1_loose, gene_f1_strict, exon_f1_exact, nt_f1 (per-genome MEANs).
- `timing_wall.tsv` (56) — genome, clade, size_mb, tool, wall_s_median, predict_wall_s_median, annotate_wall_s_median (stage cols blank for MAKER/BRAKER).
- `timing_cpu.tsv` (56) — genome, tool, cpu_s_median.
- `rss.tsv` (56) — genome, tool, peak_rss_mb_median, predict_rss_mb_median, annotate_rss_mb_median.
- `annotate_per_source.tsv` (154) — genome, variant(matched|full), source, seconds_median.
- `scaling.tsv` (56) — genome, size_mb, clade, tool, wall_min_median.
- `gene_counts.tsv` (56) — genome, tool, n_pred_genes_mean, n_ref_genes.
- `clade_strata.tsv` (56) — same columns as accuracy.tsv (clade-stratified view).
- `per_run_raw.tsv` (350) — every run: variant, tool, genome, clade, size_mb, rep, all four F1 metrics, gene counts, wall_s, cpu_s, peak_rss_mb, predict/annotate wall + RSS. **This is the per-run table the Phase-6 statistics consume.**
- `completeness.tsv` (28) — BUSCO-of-predicted-proteins, **MycoNote + funannotate ONLY** (MAKER/BRAKER dirs carry no protein-BUSCO). CAVEAT: different lineages/BUSCO versions (myconote fungi_odb10 n=758 BUSCO v5; funannotate dikarya_odb9 n=1312 BUSCO v3) → NOT a fair cross-tool axis; no figure drawn from it.
- `_braker_history_146852.tsv` — raw condor_history dump backing BRAKER timing (ProcId, ExitCode, RemoteWallClockTime, RemoteUserCpu, RemoteSysCpu, MemoryUsage, ResidentSetSize, NumJobStarts, Args).

## Reproduce
```
# on the cluster (light JSON parsing, submit-node safe):
cd ~/myconote-cli/docs/paper/htcondor_benchmark
condor_history 146852 -af:t ProcId ExitCode RemoteWallClockTime RemoteUserCpu \
   RemoteSysCpu MemoryUsage ResidentSetSize NumJobStarts Args > figures/_braker_history_146852.tsv
python3 scripts/aggregate_bench.py            # writes figures/*.tsv
# locally (needs matplotlib):
python3 scripts/plot_bench.py <tsv_dir> <out_dir>   # writes fig1..fig12 .pdf + .png
```

## Dimensions NOT built (missing data, not fabricated)
- **Cross-tool BUSCO completeness figure** — protein-BUSCO exists only for MycoNote + funannotate, at different lineages/versions. `completeness.tsv` is emitted with the caveat; no figure.
- **BRAKER / MAKER predict-vs-annotate stage split** (fig4, fig6b) — neither tool emits per-stage instrumentation; shown as totals (fig4 hatched) or omitted from the stage panel (fig6b).
