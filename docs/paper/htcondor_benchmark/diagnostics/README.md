# MycoNote accuracy diagnostic — where the gap is, and the one fix worth building

Author: Benjamin Narh-Madey
Date: 2026-11-04
Scope: reviewer-requested diagnostic over the EXISTING 14-genome ab-initio benchmark
(`results_timing_stage1_matched/`, `results_timing_stage1/`). **No fixes implemented** —
this is an evidence-gathering pass that localises the accuracy gap and names the single
highest-value fix before any code is written.

## What this does and why
Reviewer #2 asked where MycoNote loses accuracy to funannotate and which single change is
worth building. This decomposes the gap by genome, by error type, by exon class, and —
decisively — asks whether MycoNote's own predictors already contain the correct gene models
that the consensus combiner then discards.

## Inputs and data versions
- MycoNote installed build on cluster: **v0.9.0** (`~/.local/bin/myconote-cli --version`).
- Benchmark prediction outputs (scored by `figures/per_run_raw.tsv`):
  - MycoNote: `results_timing_stage1_matched/myconote/<g>/rep1/predict_out/`
    (`consensus.gff3`, `augustus.gff3`, `snap.gff3` [yeasts only], `genemark/genemark.gff3`,
    `evidence_weights.toml`, `predict_summary.txt`) + `rep1/rename_table.tsv`.
  - funannotate: `results_timing_stage1/funannotate/<g>/rep1/` (`predict_out/predict_results/*.gff3`,
    `predict_out/annotate_results/*.gff3`, `rename_table.tsv`, `run.log`, `metrics.json`).
  - References: `data/<g>/reference.gff3` (RefSeq, NC_ accessions).
- Scorer: `scripts/compare_annotations.py` (Eilbeck et al. 2009; coding-only;
  strict = identical ordered CDS interval set; loose = same-strand span overlap). Seqids lifted
  to the reference namespace via each tool's `rename_table.tsv`.
- 14 genomes: sce cal ylp ani ncr cne klac kpha fgra afum umay ccin spom tmel.

## How these numbers were produced (reproduce)
All analysis ran LOCALLY on GFF3 copies (no submit-node compute). The local scorer was
verified to reproduce the benchmark exactly for sce (consensus vs ref → loose 0.933,
strict 0.819 = `per_run_raw` 0.9326 / 0.8187).

```
# copy GFF bundle (5 genomes) from cluster
ssh glbrc 'cd <benchmark>; tar czf - <consensus/augustus/snap/genemark/rename/ref/funannotate gff3>' > bundle.tgz
# scripts in this dir (edit the hardcoded SCORER path at the top of each if the repo moves):
python3 q3_taxonomy.py  <myc consensus.gff3> <reference.gff3> <myc rename_table.tsv> "<label>"
python3 q6_exon_split.py <pred gff3>          <reference.gff3> <rename_table.tsv>    "<label>"
python3 q5_lost_votes.py <g> <ref> <myc_rename> <fun_gff> <fun_rename> <consensus> <augustus> <snap|none> <genemark|none>
python3 q5b_oracle.py    <g> <ref> <myc_rename> <consensus> <augustus> <snap|none> <genemark|none>
```

---

## Findings by question

### Q4 — the gap is concentrated, not uniform (from `figures/per_run_raw.tsv`, mean over 5 reps)
MycoNote−funannotate gene-F1 per genome (positive = MycoNote wins):

| genome | clade | myc loose | fun loose | Δloose | myc strict | fun strict | **Δstrict** |
|--------|-------|-----------|-----------|--------|------------|------------|-------------|
| sce  | Sacch | 0.9326 | 0.9481 | −0.0155 | 0.8187 | 0.8644 | −0.0458 |
| cal  | Sacch | 0.9487 | 0.9626 | −0.0140 | 0.8171 | 0.8326 | −0.0155 |
| ylp  | Sacch | 0.9554 | 0.9533 | +0.0021 | 0.7635 | 0.7320 | **+0.0315** |
| ani  | Asco  | 0.8995 | 0.9287 | −0.0293 | 0.4448 | 0.4914 | −0.0466 |
| ncr  | Asco  | 0.9354 | 0.9252 | +0.0102 | 0.5882 | 0.5559 | **+0.0323** |
| cne  | Basidio | 0.9336 | 0.9490 | −0.0155 | 0.5316 | 0.5014 | **+0.0302** |
| klac | Sacch | 0.9614 | 0.9709 | −0.0094 | 0.8454 | 0.8573 | −0.0119 |
| **kpha** | Sacch | 0.9126 | 0.9438 | −0.0312 | 0.5691 | 0.7500 | **−0.1809** |
| fgra | Asco  | 0.9217 | 0.9097 | +0.0121 | 0.5415 | 0.5382 | +0.0033 |
| afum | Asco  | 0.9092 | 0.9396 | −0.0303 | 0.5135 | 0.5386 | −0.0250 |
| umay | Basidio | 0.9291 | 0.9549 | −0.0258 | 0.7172 | 0.7603 | −0.0430 |
| ccin | Basidio | 0.8925 | 0.8985 | −0.0060 | 0.4468 | 0.4390 | +0.0078 |
| spom | Asco  | 0.9346 | 0.9680 | −0.0333 | 0.7768 | 0.8688 | −0.0920 |
| tmel | Asco  | 0.5175 | 0.5476 | −0.0301 | 0.1564 | 0.1771 | −0.0208 |

- **Loose (locus-finding) is close everywhere** (|Δ| ≤ 0.033; MycoNote wins 4/14). MycoNote is not
  missing loci — it finds the genes.
- **The gap lives in STRICT (exact CDS structure)** and is **concentrated**: kpha (−0.181), spom
  (−0.092), ani (−0.047), sce (−0.046), umay (−0.043) carry nearly all of it. MycoNote **wins strict**
  on ylp/ncr/cne (+0.030 to +0.032) and ties fgra/ccin.
- MycoNote over-predicts gene count almost everywhere (e.g. sce 5773 vs ref 6021 but 5488 for fun;
  tmel 16007 vs ref 7420 — tmel is a shared failure for both tools, an outlier).

### Q3 — error taxonomy (q3_taxonomy.py; MycoNote consensus vs reference)
| category | sce | ncr | fgra |
|----------|-----|-----|------|
| exact (correct CDS structure) | 4738 (82.0% of pred) | 5506 (60.8%) | 6844 (57.0%) |
| **correct-locus-wrong-boundaries** | 390 (6.7%) | **2974 (32.9%)** | **4184 (34.9%)** |
| merged (1 myc spans ≥2 ref) | 210 (covers 445 ref) | 208 (431 ref) | 561 (1163 ref) |
| split (≥2 myc in 1 ref) | 101 | 52 | 35 |
| tangled (≥2×≥2) | 7 | 17 | 14 |
| extra (myc, no ref overlap) | 225 (3.9%) | 221 (2.4%) | 314 (2.6%) |
| missed (ref, no myc overlap) | 329 (5.5% of ref) | 758 (7.8%) | 1055 (7.9%) |
| **among 1:1-matched loci: wrong-boundaries** | **7.6%** | **35.1%** | **37.9%** |

- The dominant error is **correct-locus-wrong-boundaries** — the gene is found, the CDS structure is
  wrong. It is ~1/3 of all matched loci in the intron-rich genomes.
- **Merged ≫ split** in every genome (sce 210 vs 101; ncr 208 vs 52; fgra 561 vs 35): MycoNote fuses
  adjacent genes far more than it fragments them.
- Missed (~6–8%) and extra (~2–4%) are small and comparable across genomes — not the driver.
- This refines the earlier sce audit ("14.3% structural"): of the structural bucket, roughly half is
  boundary error, a quarter merged, the rest split/tangled.

### Q6 — single- vs multi-exon strict F1 (q6_exon_split.py)
| genome | class | myc strict F1 | fun strict F1 |
|--------|-------|---------------|---------------|
| sce  | single | 0.8627 | 0.8822 |
| sce  | **multi**  | **0.3116** | 0.5117 |
| ncr  | single | **0.6861** | 0.6396 |
| ncr  | multi  | **0.5593** | 0.5306 |
| fgra | single | **0.6650** | 0.6302 |
| fgra | multi  | 0.5025 | 0.5084 |

**Overturns the expectation that single-exon is the weak spot.** In the intron-rich genomes MycoNote
*beats* funannotate on single-exon strict. The sce gap is driven by **multi-exon over-calling**:
MycoNote emits **608** multi-exon models where the reference has **329** (fun emits 226). Multi-exon
sensitivity is comparable (myc 146/329 vs fun 142/329) but specificity collapses (myc 146/608 = 0.24
vs fun 142/226 = 0.63) — MycoNote invents spurious introns/structures in a genome that is 95%
single-exon.

### Q5 — lost consensus votes (q5_lost_votes.py) — **the crux**
On reference loci where funannotate is strict-correct but MycoNote consensus is NOT, how often did a
MycoNote **raw predictor** (augustus/snap/genemark) already hold the exact-correct CDS structure?

| genome | gap-loci (fun right, myc wrong) | raw predictor HAD the exact model | which predictor held it |
|--------|------|------|------|
| sce  | 418 | **335 (80.1%)** | snap 209, snap+gm 95, gm 31 |
| kpha | 1210 | **1102 (91.1%)** | snap+gm 724, snap 213, gm 165 |
| spom | 418 | **286 (68.4%)** | snap+gm 120, gm 112, snap 54 |
| ncr  | 535 | **321 (60.0%)** | genemark 321 (no SNAP ran) |
| fgra | 773 | **555 (71.8%)** | genemark 555 (no SNAP ran) |

**Augustus never appears as a sole holder of a discarded-correct model** — because Augustus (weight 10)
always wins its locus. Every discarded-correct model was SNAP or GeneMark being overruled.

### Q5b — oracle-combiner headroom (q5b_oracle.py) — size of the prize
Across ALL reference coding genes: strict-correct genes the shipped consensus gets, vs what a perfect
combiner (always pick whichever raw predictor is exactly correct) would get.

| genome | consensus strict | oracle strict | recoverable | funannotate strict (sens) |
|--------|------------------|---------------|-------------|---------------------------|
| sce  | 4832 (80.3%) | 5235 (86.9%) | **+403 (+6.7 pts)** | 4974 (82.6%) |
| kpha | 2773 (55.0%) | 4168 (82.7%) | **+1395 (+27.7 pts)** | 3730 (74.0%) |
| spom | 4143 (80.9%) | 4509 (88.0%) | **+366 (+7.1 pts)** | 4374 (85.4%) |
| ncr  | 5532 (56.7%) | 6169 (63.2%) | **+637 (+6.5 pts)** | 5103 (52.3%) |
| fgra | 6855 (51.5%) | 7739 (58.1%) | **+884 (+6.6 pts)** | 6700 (50.3%) |

**A better combiner alone would lift strict sensitivity by +6.5 to +27.7 points and push MycoNote above
funannotate's strict sensitivity in all five genomes** — including kpha (55%→83% vs fun 74%), the single
worst genome in Q4.

### Mechanism — why correct models are discarded (confirmed in source)
`src/predict/evidence.rs:862`, the winner-take-all consensus score per overlap group:
```
score(x) = weight(x) + Σ_{y≠x in group} weight(y) · cds_jaccard(x, y)
```
Default fungal weights (`evidence_weights.toml`, confirmed in every `predict_summary.txt`):
Augustus **10**, GeneMark **5**, SNAP **3** (hand-picked, never tuned — Q8).

With Augustus wrong (structure A) and SNAP+GeneMark both correct (structure B, jaccard(S,G)=1), let
j = jaccard(A,B):
- score(Augustus) = 10 + (3+5)·j = 10 + 8j
- score(SNAP)     = 3 + 10·j + 5·1 = 8 + 10j
- score(GeneMark) = 5 + 10·j + 3·1 = 8 + 10j

Augustus wins whenever `10 + 8j > 8 + 10j` ⇔ `j < 1`. So **for any disagreement (j<1), two corroborating
and correct lower-weight predictors (combined weight 8) cannot outvote one wrong Augustus (weight 10).**
This is exactly the empirical Q5 result. The 10-vs-3/5 weighting makes Augustus structurally decisive;
when Augustus is wrong (bad proxy species — kpha uses `pichia_stipitis`, tmel `aspergillus_nidulans`),
the combiner has no path to the correct model its own predictors produced.

### Q2 — funannotate's EVM inputs + weights (from `run.log`, identical across genomes)
```
Source         Weight   Count (sce / ncr)
Augustus       1        2788 / 6876
Augustus HiQ   2        1608 / 352     (>90% exon-evidence subset)
GeneMark       1        5415 / 9282
GlimmerHMM     1        5681 / 9879
snap           1        5047 / 8905
```
Plus 562,726 proteins mapped (diamond+exonerate) as exon evidence (this feeds the HiQ subset — a nuance
vs the relayed "no external proteins": funannotate used its bundled protein evidence set, even with no
user proteins). **Contrast with MycoNote:** funannotate runs **5 near-flat tracks** (max ratio 2:1) and
lets EVM optimise each partition, so a non-Augustus source can win a locus. MycoNote runs **2–3 tracks
at 10:5:3** with winner-take-all, so Augustus cannot be outvoted.

### Q9 (relayed) is CONTRADICTED by the artifacts — flag
Relayed: "benchmark MycoNote ran Augustus+SNAP only (GeneMark off by default)." Evidence
(`predict_summary.txt`):
- sce (both matched and non-matched runs): **Augustus + SNAP + GeneMark** (weights 10/3/5), consensus 5780/5785.
- ncr, fgra: **Augustus + GeneMark only — SNAP did NOT run** (not in "Predictors used").
So GeneMark was ON in the scored benchmark, and SNAP was dropped specifically on the two intron-rich
genomes. Both halves of Q9 are wrong. This strengthens, not weakens, the diagnosis: GeneMark's
correct-but-discarded models are a real, scored phenomenon (Q5). **Open sub-question for Benjamin: why
did SNAP silently not run on ncr/fgra?** (worth a look, but orthogonal to the combiner fix.)

### Q10 — annotate time by source + threads (`figures/annotate_per_source.tsv`, median seconds)
All five functional steps are invoked with `--cpu/--threads = config.threads` (=**8** in the matched run;
confirmed `pfam.rs:108`, `mmseqs.rs:67`, `merops.rs:178/263`, `cazyme.rs:83/147/219`, `busco.rs:194/235`;
benchmark thread count `jobs/run_myconote.sh:215/226/249`).

| source | sce | ncr | fgra | ccin | tmel |
|--------|-----|-----|------|------|------|
| **pfam (hmmsearch)** | **1024** | **1526** | **1449** | **1722** | **1854** |
| mmseqs (Swiss-Prot) | 53 | 81 | 111 | 86 | 104 |
| busco | 52 | 73 | 74 | 74 | 73 |
| cazyme (dbCAN) | 31 | 42 | 39 | 45 | 57 |
| merops | 1.6 | 2.6 | 2.5 | 2.2 | 3.6 |

**Pfam/hmmsearch is 85–95% of matched-annotate wall** even at 8 threads (hmmsearch `--cpu` scales
sublinearly). Consistent with the standing decision (`annotate_matched_speed_gap`) to keep hmmsearch for
portability. The "full" variant adds eggnog (3.4k–7.5k s), by far the largest single cost. Annotate cost
is a SPEED issue, not an accuracy one — it does not bear on the gap.

### Q7 — protein-hints A/B (prior evidence + current-build confirmation PENDING)
- **Prior (v0.7.11, commit `5ea962a`, documented in `AB_v0711_selftrain_prothints_README.md`):**
  protein-hints did **not** lift strict F1 and slightly lowered it — ani 0.445→0.432, cne 0.526→0.491;
  loose flat (verified from `results/ab/prothints70/*/metrics.json`). Self-training also no lift. Logged
  as "accuracy levers closed."
- **Current build (v0.9.0) confirmation: SUBMITTED, not yet complete.** `condor_submit` cluster **148567**
  (`~/ab_v077/ab_q7diag.sub`), 4 jobs = {stock, prothints70} × {ncr, cne}. A `stock` (no-flag) arm was
  added to `~/ab_v077/run_ab.sh` (backup `run_ab.sh.bak_q7diag`). Collect when done:
  ```
  for a in stock prothints70; do for g in ncr cne; do
    python3 -c "import json;d=json.load(open('results/ab/$a/$g/metrics.json'));print('$a $g',d['gene_level']['loose']['F1'],d['gene_level']['strict']['F1'])"; done; done
  ```
  **Expectation (mechanistic):** protein-hints adds a weight-20 track but only where miniprot produces a
  confident close-relative alignment; it does not touch the Augustus-dominance problem on the bulk of
  loci, so a large lift is not expected. Jobs were idle behind ~200 of Benjamin's other jobs at submit.

---

## THE ONE FIX WORTH BUILDING: the consensus combiner (vote logic / weights), NOT a new/better predictor

The data converge on a single conclusion:
1. MycoNote **finds** the genes (loose F1 ≈ funannotate everywhere).
2. The gap is **strict/boundary structure**, concentrated in specific genomes (Q4, Q3, Q6).
3. On **60–91%** of the loci where MycoNote loses, **its own predictors already produced the exact-correct
   model** (Q5), which the combiner discarded.
4. A perfect combiner recovers **+6.5 to +27.7 strict points** and **beats funannotate's strict
   sensitivity in all 5 tested genomes** (Q5b) — with zero new predictors and zero new compute.
5. The mechanism is proven, not inferred: the `10 + 8j > 8 + 10j` arithmetic (`evidence.rs:862`) means two
   corroborating correct predictors can never outvote one wrong Augustus. The weights are hand-picked and
   never tuned (Q8).

**Why not the other candidates:** better Augustus training / protein hints were already closed (no lift,
Q7 prior) — and the combiner finding explains *why*: improving one predictor cannot help when the combiner
discards correct models from the others. Enabling GeneMark is moot — it is already on (Q9). Annotate speed
(Q10) is not an accuracy lever.

**Concrete directions to prototype (for the build phase, not done here):**
- Rebalance weights (Augustus is 10× SNAP; funannotate's working system is ≤2:1) and/or add a corroboration
  bonus so ≥2 agreeing predictors can overrule a lone high-weight one; tune against the 5-genome oracle.
- Or replace winner-take-all with EVM-style per-partition optimisation (funannotate's proven approach).
- Validate against the oracle ceilings in Q5b (sce 86.9%, kpha 82.7%, spom 88.0%, ncr 63.2%, fgra 58.1%).

## Needs a re-run vs answerable from existing files
- **Answerable from existing files (done):** Q2, Q3, Q4, Q5, Q5b, Q6, Q9, Q10 — all per-predictor GFFs
  (`augustus.gff3`, `snap.gff3`, `genemark/genemark.gff3`) ARE retained in `predict_out/`, so Q5 needed no
  re-run.
- **Needed a re-run (submitted, pending):** Q7 current-build (v0.9.0) protein-hints confirmation — cluster
  148567. Prior v0.7.11 evidence stands in the interim.

## Caveats
- Q3/Q5/Q5b/Q6 computed on 3–5 representative genomes (sce easy yeast; ncr/fgra intron-rich; kpha/spom the
  largest strict losers), not all 14 — chosen to span the gap, not sample it uniformly.
- Strict = exact CDS-interval match (Eilbeck); it is unforgiving of single-codon boundary shifts, which is
  the intended stringency for a boundary diagnosis.
- The oracle combiner is an upper bound (it assumes perfect selection among raw models); a real combiner
  will recover a fraction of it. The point is the *headroom exists inside MycoNote's own outputs*.
- tmel is a shared outlier (both tools over-predict ~2×); excluded from the fix rationale.
