# Benchmark figures — per-test, work-matched (v0.7.13, matched AVX2)

Five independent comparisons from the definitive four-tool HTCondor benchmark
(MycoNote / funannotate / MAKER / BRAKER; 6 fungal genomes × 3 reps; all arms on
the same AVX2-capable node class). Each test is its own figure + data file so
**gene predictors and complete annotation pipelines are never conflated on a
speed axis** — the methodological fix that motivated this redesign.

Regenerate everything: `python3 ../scripts/make_figures.py --results-dir ../results --out-dir .`
(Python 3 + matplotlib + numpy.) Per-run source data is in `../data/`.

## The five tests

| # | figure | data | what it shows |
|---|--------|------|---------------|
| 1 | `figure_accuracy.{pdf,png}` | `../data/accuracy.tsv` | gene-level loose / strict / nucleotide F1, all 4 tools |
| 2 | `figure_prediction_speed.{pdf,png}` | `../data/prediction_speed.tsv` | gene-prediction wall-time, all 4 tools (BRAKER timing unusable) |
| 3 | `figure_annotation_matched.{pdf,png}` | `../data/annotation_matched.tsv` | functional-annotation wall-time at **matched sources** — MycoNote (no EggNog) vs funannotate, with MycoNote-full (+EggNog) for reference |
| 4 | — (table) | `../data/completeness.tsv` | which functional outputs each tool produces |
| 5 | — (table) | `../data/portability.tsv` | AVX2 node-class availability |

## Results (means of 18 runs)

**Accuracy** — competitive, not the leader. Loose: funannotate/BRAKER lead;
strict: MAKER (yeasts) / BRAKER (intron-rich); nucleotide: BRAKER leads, MycoNote
beats funannotate pairwise. BRAKER3 is the strongest accuracy tool overall.

**Speed — work-matched.** MAKER and BRAKER do **no functional annotation**
(prediction only), so they appear only in the prediction and accuracy tests, never
in an end-to-end claim.

| stage | MycoNote | funannotate | MAKER |
|-------|----------|-------------|-------|
| gene prediction | **26.0 min** | 77.0 min | 30.2 min |
| functional annotation, matched sources (no EggNog) | 17.2 min | **9.0 min** | — none — |
| functional annotation, MycoNote full (+EggNog) | 73.2 min | — | — |
| **end-to-end at matched sources** | **~43 min** | ~86 min | n/a |

- **Prediction: MycoNote fastest** (~3× funannotate).
- **End-to-end at matched functional output: MycoNote ~2× faster** (~43 vs ~86 min).
  The earlier "funannotate faster end-to-end" reading was an artifact of MycoNote
  also running EggNog, which funannotate did not.
- **Functional annotation at matched sources: funannotate is faster** (9 vs 17 min).
  Profiled per-source (cne): the gap is **entirely Pfam/hmmsearch** (~11.7 min,
  ~78% of the matched annotate; MMseqs2 ~2.8, BUSCO ~0.6, CAZyme/MEROPS ~0).
  funannotate wins this sub-stage because it uses pyhmmer (same HMMER engine,
  faster bindings). MycoNote deliberately keeps the standard, portable,
  dependency-light `hmmsearch` binary — documented here as an honest cost, not a
  positioning problem (MycoNote still wins prediction and matched end-to-end).
- **EggNog** adds ~56 min for orthology/GO/KEGG that funannotate (as-run) and
  MAKER/BRAKER do not produce — a completeness advantage (Test 4), not a speed
  penalty.

**Portability (Test 5):** MycoNote (Augustus 3.3.3) and MAKER run the full node
pool; funannotate/BRAKER (Augustus 3.5.0, AVX2-compiled) require AVX2 nodes.

## Caveats
- BRAKER wall-time is excluded from the speed figure: its `performance.json`
  records all-zeros (a known BRAKER instrumentation issue).
- antiSMASH and InterProScan are excluded from both pipelines' functional runs
  (they need network access from execute nodes).
- Timing means are right-skewed by the larger intron-rich genomes; see the
  per-run data files for the full distribution.
