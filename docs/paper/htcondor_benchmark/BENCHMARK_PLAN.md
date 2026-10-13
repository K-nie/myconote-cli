# Benchmark plan — MycoNote-CLI vs funannotate / MAKER / BRAKER (publication-grade)

Governing design document and Methods raw material for the definitive four-tool
fungal-annotation benchmark. Synthesises the harness-engineering and
statistical-design reviews (2026-10-13). Execution is phased with the steps in §7.

## 1. Replication unit (the central design decision)
The independent replication unit is the **genome**, not the per-tool jobs.
- **Accuracy** is near-deterministic (within-genome rep SD ≤ 0.001 F1). The 3
  reps are a determinism audit only — collapse to the per-genome **mean** and
  analyse at the genome level. Treating 18 jobs as n=18 is pseudoreplication.
- **Timing** reps carry real HPC-noise variance — keep them, analyse on
  per-genome **medians** / log-ratios; reps are nested within genome.

## 2. Genome panel (expanded + balanced + scaling)
Target ~13–16 genomes, ≥4 per subphylum, ≥3 Basidiomycota, +1 Taphrinomycotina
outgroup, +scaling tests. Curated references required for all.

| clade | genomes (accession) |
|-------|---------------------|
| Saccharomycotina | sce, cal, ylp (existing) + *K. lactis* GCF_000002515.2, *K. phaffii* GCF_001708105.1 |
| Pezizomycotina | ani, ncr (existing) + *F. graminearum* GCF_000240135.3, *A. fumigatus* GCF_000002655.1 |
| Basidiomycota | cne (existing) + *U. maydis* GCF_000328475.2, *C. cinerea* GCF_000182895.1 |
| Taphrinomycotina | *S. pombe* GCF_000002945.2 (PomBase, gold-standard outgroup) |
| scaling | *Tuber melanosporum* GCF_000151645.1 (~125 Mb, repeat-rich) or *Z. tritici* GCF_000219625.1 (~40 Mb); + one controlled fragmentation test (shred an existing genome to draft N50, same reference) |

## 3. Measurement
- **Enforce AVX2** `requirements = (TARGET.has_avx2 =?= true)` on ALL 4 `.sub`
  files (committed) + **same-node pinning** for the timing production run.
- **BRAKER timing:** wrap `braker.pl` in `/usr/bin/time -v` (its own
  `performance.json` is all-zeros); log-span fallback only if wrapping fails.
- **RSS:** report per-stage (predict vs annotate) peak RSS on a log-axis figure.
- **Per-source annotate timing:** instrument `annotate` to emit per-source
  wall-time (swissprot/pfam/busco/cazyme/merops/eggnog) — the Pfam-dominates
  finding needs error bars.
- Pinned-node **wall-time** is the headline; CPU-seconds is the robustness check.

## 4. Statistics (randomized complete block; genome = block)
- Per metric (loose/strict/nt F1): **Friedman** omnibus + **Nemenyi** post-hoc
  (exploratory); **exact Wilcoxon signed-rank** for the pre-registered
  MycoNote-vs-funannotate primary contrast.
- Effect sizes: rank-biserial, Cliff's δ, **exact 64-permutation CI**; Kendall's W.
- Timing: LMM `log(time) ~ tool + (1|genome) + (tool|genome)` (Kenward–Roger),
  with **median-based Friedman as the outlier-robust primary**. Fair arms only.
- Multiplicity: **Holm** for the confirmatory family (myco-vs-fun × 3 accuracy +
  2 speed); **BH-FDR** for the exploratory 4-tool grid; per-metric separate;
  report corrected + uncorrected.
- Clade = **descriptive stratum, not a tested factor** (Basidiomycota n=1 →
  tool×clade non-estimable). Report the pattern descriptively.
- Equivalence claims need **TOST** (margin ±0.02 F1), not NHST non-significance.
- Power: n=6 detects only ΔF1 ≳ 0.05 consistent; ~10–15 for ΔF1 ≈ 0.05; ~25–50
  for ΔF1 ≈ 0.02–0.03 → the panel expansion (§2) is the power lever.

## 5. Deeper risks to pre-empt
- **Reference circularity:** RefSeq/Gnomon references are themselves pipeline
  products → add a reference-independent axis (**BUSCO completeness of each
  tool's protein set**), note per-genome reference provenance (SGD/PomBase manual
  vs Gnomon), frame accuracy as "agreement with the community reference."
- **Reproducibility artifacts:** update `tools.tsv` (0.7.5 → actual tag/SHA);
  create `db_versions.tsv`; implement the per-run `hardware.log`; pin both
  Augustus builds + GeneMark.

## 6. Framing (locked)
MycoNote is **competitive, not the accuracy leader** (BRAKER3 strongest overall);
durable wins = fastest prediction, out-of-box EggNog completeness, x86-64-v2
portability. Every null/loss stays in the text.

## 7. Phased execution (status tracked here as it runs)
1. **Harness fixes (local, no cluster):** AVX2 on all subs + pinning; BRAKER
   time-wrap; `hardware.log`; `db_versions.tsv`; `tools.tsv`; per-genome
   dispersion + RSS + scaling in `make_figures.py`; per-source annotate
   instrumentation (Rust code change).
2. **2-rep smoke** on the existing 6 genomes to validate every change.
3. **Panel download + curation** (new genomes + scaling genomes).
4. **Stage-1 (n=5 timing), same-node pinned,** full expanded panel.
5. **Production** (n=10 headline timing; 3-rep determinism audit for accuracy)
   → `collect_results.sh` → `make_figures.py`.
6. **Statistics** (tests/CIs/effect sizes) → Methods/Results.

READMEs (this file, README.md, figures/README.md, db_versions, hardware) are
updated as each stage runs — they are the manuscript Methods raw material.
