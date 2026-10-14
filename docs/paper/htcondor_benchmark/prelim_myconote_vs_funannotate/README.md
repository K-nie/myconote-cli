# Preliminary MycoNote vs funannotate comparison (1-rep smoke)

**Status: PRELIMINARY.** This is a 1-replicate comparison built from the
expanded-panel smoke runs while the MAKER/BRAKER smokes finish. The production
matrix (n=10 timing reps, 3 accuracy reps, all four tools, proper statistics)
supersedes these numbers — do not cite these as final.

## What this does / why

Side-by-side MycoNote (v0.7.15) vs funannotate (v1.8.17) on the 14-genome fungal
panel: gene / exon / nucleotide accuracy vs the RefSeq reference, plus wall-clock
time. Gives an early read on where MycoNote stands before the full run.

## Inputs / provenance (pinned)

- **MycoNote** v0.7.15 — `results_smoke_expanded_v0715/myconote/<id>/rep1`
  (all 14 genomes, fresh this session; HMMER CAZyme path).
- **funannotate** v1.8.17:
  - 8 new genomes (klac, kpha, fgra, afum, umay, ccin, spom, tmel) — fresh this
    session, `results_smoke_expanded_competitors/funannotate/<id>/rep1`.
  - 6 original genomes (sce, cal, ylp, ani, ncr, cne) — prior fair-run `rep1`
    from `results/funannotate/<id>/rep1` (same v1.8.17, same references, same
    `compare_annotations.py`). CAVEAT: different run session/node than the 8 new
    ones, so cross-genome *timing* carries that provenance difference; accuracy
    is version-/reference-matched and comparable.
- Accuracy metric: `scripts/compare_annotations.py` (Eilbeck-style gene/exon/
  nucleotide TP/FP/FN → F1), identical for both tools, against the same
  `data/<id>/reference.gff3`.
- Timing: each run's `performance.json` (`/usr/bin/time -v` per stage; AVX2-
  pinned nodes both arms).

## Reproduce

```
# pull the per-run JSONs, collate, and render the figure
python3 compare_myco_fun.py          # (uses ~/miniconda3 python with matplotlib)
```

`compare_myco_fun.py` reads the pulled `metrics.json` + `performance.json` for
both tools and writes `comparison_myconote_vs_funannotate.tsv` and
`figure_myconote_vs_funannotate.{pdf,png}`.

## Outputs

- `comparison_myconote_vs_funannotate.tsv` — per-genome, per-tool: predicted
  gene count, gene F1 (loose/strict), exon F1 (exact), nucleotide F1, predict /
  annotate / total wall-minutes, peak RSS (GB).
- `figure_myconote_vs_funannotate.{pdf,png}` — 4 panels: (a) gene F1 loose,
  (b) gene F1 strict, (c) nucleotide F1, (d) total wall-clock. MycoNote (blue) vs
  funannotate (orange), Okabe-Ito; `*` marks the 8 new genomes.

## Results (win tally across 14 genomes)

| metric | better | MycoNote | funannotate | tie |
|---|---|---|---|---|
| gene F1 (loose) | higher | 1 | 13 | 0 |
| gene F1 (strict) | higher | 6 | 8 | 0 |
| exon F1 (exact) | higher | 5 | 9 | 0 |
| nucleotide F1 | higher | 10 | 4 | 0 |
| total wall-clock | lower | 10 | 4 | 0 |

Read: funannotate has a small, consistent edge on gene-level detection (loose
F1, margins ~0.01–0.04); the two are roughly even on exact CDS structure
(strict); MycoNote leads on nucleotide-level accuracy and is faster overall on
most genomes.

## Caveats / things to run down before production

- **Proxy Augustus models underperform (expected).** MycoNote's two genomes with
  no species-exact Augustus model — kpha (`pichia_stipitis`) and tmel
  (`aspergillus_nidulans`) — are its weakest: kpha strict F1 0.568 vs 0.758,
  tmel nucleotide F1 0.561 vs 0.664. This is the flagged nearest-clade-proxy
  penalty, not a general accuracy gap. funannotate trains Augustus per-genome via
  `--busco_seed_species`, which sidesteps it (a design difference the benchmark
  characterizes).
- **spom likely a scoring artifact.** Both tools score gene-loose F1 ~0.58 while
  nucleotide F1 is ~0.98. High nucleotide + low gene-level for *both* tools
  points to a spom reference/gene-ID matching issue in `compare_annotations.py`,
  not a real failure — investigate before the production run.
- **tmel is hard for both** (repeat-rich 125 Mb truffle); low F1 across the board.
- **1 replicate only.** No variance, no significance testing. The production
  matrix adds reps and the Friedman/Wilcoxon/LMM statistics from BENCHMARK_PLAN.md.
