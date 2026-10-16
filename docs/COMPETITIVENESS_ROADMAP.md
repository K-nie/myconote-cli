# MycoNote competitiveness roadmap

Synthesis of a three-agent research effort (tool-engineering, ML/DL gene finders,
community/forum mining) + the 14-genome 4-tool benchmark + verification of what
each competitor arm actually ran. Internal planning doc: every external claim
carries its source; mechanism-justified gains with no isolated fungal point
estimate are marked "measure on our 14-genome set before any manuscript claim."

## 1. Diagnosis (verified, not inferred)

Benchmark headline (per-genome collapse, n=5): gene F1 loose — BRAKER3 0.930 >
funannotate 0.914 > MAKER 0.906 > **MycoNote 0.899**; strict — BRAKER 0.659 >
funannotate 0.636 > MAKER 0.625 > **MycoNote 0.609**.

Two facts reframe the whole problem:

- **The gap is small and it's a *strict*-F1 (CDS-structure) problem, not loose.**
  Loose: −1.5 pts vs funannotate, −3.1 vs BRAKER. Strict: −5.0 vs BRAKER. The
  priority is CDS-structure correctness / precision, not gene detection.
- **No tool used RNA-seq** (verified in the wrappers). BRAKER3 ran protein-only
  (OrthoDB Fungi, `--prot_seq`); funannotate ran BUSCO-seed per-genome Augustus
  training; MAKER ran pure ab-initio (`keep_preds=1`); MycoNote ran pure
  ab-initio consensus (Augustus 3.3.3 + SNAP + GeneMark-ES, agreement-weighted).
  So BRAKER3's famous ~+20 transcript-F1 RNA-seq lift is **not** in play here —
  the gap MycoNote must close is the 1.5–5 pt *structural* gap, not 20 pts.

Where each competitor's edge comes from, given what they ran:
- vs **BRAKER3** → protein evidence (OrthoDB) + TSEBRA transcript selection.
- vs **funannotate** → per-genome Augustus training (BUSCO-seeded) + **EVM
  weighted consensus** (not flat agreement).
- vs **MAKER** (also pure ab-initio, yet 0.906 > 0.899) → MycoNote's 3-way
  consensus isn't currently buying accuracy over one well-configured Augustus.

**The fungal ceiling caveat (most important for triage):** in compact,
intron-poor fungal genomes the ab-initio ceiling is already high, so
evidence-integration pays far less than in vertebrates. GeneMark-EP+ lifted
sensitivity on *N. crassa* (a fungus) by only ~2% even with the largest protein
set, "because accuracy of GeneMark-ES was high as typical for fungal genomes"
([GeneMark-EP+, NARGAB 2020](https://academic.oup.com/nargab/article/2/2/lqaa026/5836691)).
This is why MycoNote's `--protein-hints` A/B found no F1 lift (memory:
`accuracy_levers_closed`) — it is consistent with the literature, largely the
data, not just an integration bug.

## 2. Prioritized roadmap (ranked by accuracy-gain × feasibility × portability)

Portability-safe, Rust/config-only changes are at the top by design — they
preserve MycoNote's genuine edges (fastest prediction, best scaling, runs on any
CPU without AVX2).

### #1 — Post-prediction filtering → raise strict F1 (do first; Rust-native, zero portability cost)
Add a filter pass after the consensus: drop single-exon genes predicted purely
ab-initio with no cross-predictor agreement; require complete ORFs (start+stop,
no internal stop) unless flagged partial; drop single-predictor models below a
length/score floor. Expose as `--min-support` / `--keep-single-exon` (deterministic).
- Precedent: BRAKER3 explicitly removes single-exon purely-ab-initio genes to
  lift precision ([BRAKER3, Genome Res 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11216308/));
  MAKER's entire QC model is AED-based filtering (AED ≤0.5 "good", a good genome
  has >90% genes ≤0.5) ([MAKER2, BMC Bioinformatics 2011](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491)).
- Our AED analogue = the cross-predictor agreement we already compute.
- Risk: filtering trades sensitivity for precision — tune on the 14 genomes and
  report F1, not precision alone.

### #2 — Upgrade the combiner: weighted vote + Augustus-HiQ tier (the mechanism funannotate beats us with)
Generalize the combiner from "count agreement" to a weighted vote with a tunable
per-predictor weight and a bonus tier for Augustus calls that agree with another
predictor (our "high-quality Augustus" analogue).
- funannotate runs EVidenceModeler with default weights Augustus=1, **Augustus
  HiQ=2**, GeneMark=1, snap=1 — it trusts (high-quality) Augustus more than the
  rest, and beats us on the *same* ab-initio inputs ([funannotate predict docs](https://funannotate.readthedocs.io/en/latest/predict.html);
  [EVidenceModeler, Genome Biol 2008](https://link.springer.com/article/10.1186/gb-2008-9-1-r7)).
- **Do it Rust-native** (keeps the pure-Rust/run-anywhere edge). Wrapping EVM v2
  would match funannotate exactly but adds a Perl/bioconda dependency — only for
  a quick parity check.
- **Do NOT port TSEBRA for the ab-initio case** — it scores transcripts by
  RNA-seq/protein support (ε/weight thresholds over evidence); with no evidence
  it degenerates ([TSEBRA, BMC Bioinformatics 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8620231/)).
  TSEBRA is the right model for our *RNA-seq* mode, not this benchmark.
- Gain is mechanism-justified (no isolated fungal point estimate) → measure first.

### #3 — Pfam annotate speed: confirm threading/chunking (we already use hmmsearch)
MycoNote already prefers `hmmsearch` over `hmmscan` (src/annotate/pfam.rs), so the
big hmmscan→hmmsearch I/O win ([Eddy, cryptogenomicon](http://cryptogenomicon.org/hmmscan-vs-hmmsearch-speed-the-numerology.html))
is already banked. The remaining lever is re-parallelization: hmmsearch leaves
~65% of work serial, and wrappers that re-chunk it (pyhmmer) cut a 97 h run to
27 h on 6 cores ([PyHMMER, Bioinformatics 2023](https://academic.oup.com/bioinformatics/article/39/5/btad214/7131068)).
Shard proteins across cores + tune `--cpu`/`-Z`. Zero portability cost.
- Optional, flag-gated only: `nail` (Rust pHMM, ~26× HMMER, not SIMD-vectorized)
  ([nail, PMC10862755](https://pmc.ncbi.nlm.nih.gov/articles/PMC10862755/)) — but
  it **requires MMseqs2** (a C++ binary), which erodes the pure-Rust/any-CPU edge.
  Keep hmmsearch the default; expose `--pfam-backend {hmmsearch,nail}` only if #3's
  free parallelization isn't enough.

### #4 — RNA-seq benchmark arm + repositioning (where the real accuracy lives)
MycoNote already ships RNA-seq (v0.7.3); the benchmark just didn't use it.
Add a second benchmark arm **with** RNA-seq so the paper isn't comparing
MycoNote-ab-initio against tools in their evidence sweet spot, and reposition
docs to "if you have reads, use them." The ~20-pt lift is an RNA-seq effect
([BRAKER3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11216308/)). Harness work →
bench-harness-engineer. No portability change.

### #5 — Tiberius as an *optional* benchmark-gated predictor (the one DL bet)
For fungi, DL is not the obvious win it is for plants/vertebrates: Helixer only
*ties* Augustus/GeneMark-ES on fungi (benchmark-only). **Tiberius** is the one
worth testing — fungal ab-initio within ~3 pts of BRAKER3, dedicated fungal model
([Tiberius multi-clade, PMC13142403](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142403/)).
Integration shape has a precedent: funannotate2 runs Helixer standalone → takes
its GFF3 → feeds EVM as high-weight evidence; same adapter pattern fits our
combiner.
- **Portability cost:** pulls in a TensorFlow stack (CPU wheels need AVX); no hard
  GPU requirement at fungal sizes but it dents the single-binary/any-CPU story →
  ship as `--predictor tiberius`, never a default.
- **Data-leakage guard:** Tiberius/Helixer ship pretrained fungal models whose
  training species may overlap our panel — hold those genomes out when benchmarking.

### #6 — Document the known-nulls / hard limits (honesty, cheap)
- Protein hints give ~2% in fungi (ceiling above) → document as known-null; only
  revisit with evolutionarily-close proteins + confidence-tiered hints
  (miniprot + **miniprothint**, the GALBA recipe) and still expect low single
  digits ([GALBA 2023](https://www.biorxiv.org/content/10.1101/2023.04.10.536199v1.full.pdf);
  [miniprot](https://github.com/lh3/miniprot)).
- **tmel (125 Mb) is a hard limit for *all* tools** (~0.52–0.61 F1): >58% TEs +
  RIP (fungal C→T genome defense that destroys homology signal)
  ([T. melanosporum genome](https://www.nature.com/articles/nature08867);
  [RIP](https://pmc.ncbi.nlm.nih.gov/articles/PMC7882544/)). Frame it as a
  known-hard outlier — MycoNote collapses *like* the others, not worse.
- Repeat-rich genomes need a de-novo RepeatModeler library + **soft-mask**
  (never hard-mask — it buries exons) ([Galaxy masking](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/repeatmasker/tutorial.html)).
  An adaptive filter (lenient for compact genomes, aggressive for Tuber-like) is
  funannotate's documented design (`--repeat_filter` off by default, on for
  high-repeat genomes).

## 3. Portability guardrails (the edges we must NOT trade away)
- **Keep Augustus 3.3.3 the default** (runs on any CPU). Competitors need
  AVX2-compiled 3.5.0. At most offer 3.5.0 as an opt-in fast path for AVX2 hosts.
- Keep **hmmsearch** the portable Pfam default; `nail`/MMseqs2 behind a flag only.
- Keep the combiner **Rust-native** (avoid EVM's Perl / Tiberius's TF as hard deps).
- #1 (filtering), #2-route-1 (weighted combiner), and #3 (hmmsearch parallelism)
  are all portability-neutral — which is exactly why they lead.

## 4. The honest positioning (for the paper)
Not "faster" or "more accurate." The defensible story: **competitive accuracy
(within 1.5–5 F1 pts on ab-initio inputs), fastest prediction stage, best scaling
to large genomes, runs on any CPU, complete one-tool pipeline to NCBI submission.**
The roadmap closes the strict-F1 gap with portability-safe changes first.

## 5. Immediate next actions (suggested order)
1. #1 filtering + #3 hmmsearch parallelism — both free, portability-safe → implement, test on the 14 genomes.
2. #2 Rust-native weighted combiner (Augustus-HiQ tier) — the funannotate mechanism, kept pure-Rust.
3. #4 RNA-seq benchmark arm (bench-harness-engineer) — the honest accuracy story.
4. #5 Tiberius optional, benchmark-gated with leakage holdout.
5. #6 document the fungal protein-hint ceiling + the tmel hard limit.
