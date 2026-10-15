# Phase 4–5 execution status — production timing + accuracy matrix

Live reproducibility record for BENCHMARK_PLAN.md §7 Phases 4 (Stage-1 timing)
and 5 (production timing + accuracy). Written as the run proceeds; this is
manuscript Methods raw material. Last updated: 2026-10-15.

## What this run does and why
Four-tool fungal-annotation benchmark (MycoNote vs funannotate / MAKER /
BRAKER) on the full 14-genome panel. Two measurement axes, per §1–§3:
- **Timing** (real HPC-noise variance): Stage-1 n=5, production n=10; analysed
  on per-genome medians. Wall-time is the headline; CPU-seconds the robustness
  check.
- **Accuracy** (near-deterministic): 3-rep determinism audit, collapsed to the
  per-genome mean.

## Tool / data versions (verified 2026-10-15)
- MycoNote-CLI: `~/.local/bin/myconote-cli --version` → **v0.7.15**. Repo HEAD
  `d43caf7` (`git describe` = `v0.7.6-45-gd43caf7`).
  CAVEAT: `configs/tools.tsv` still records `0.7.14` / commit `046520d` — STALE.
  Needs an accurate update to the v0.7.15 build provenance before submission
  (§5 reproducibility item). Not changed here: the exact tag/SHA the installed
  binary was built from is not established on-cluster, and a guessed SHA in a
  reproducibility file is worse than a flagged gap.
- funannotate 1.8.17, MAKER 3.01.04, BRAKER 3.0.8 (per tools.tsv; conda envs
  `funannotate`, `maker`, `myconote_braker3`).
- Panel: 14 genomes in `configs/genomes.tsv`; all have `data/<g>/genome.fa` +
  `reference.gff3` (verified present + non-empty for all 14).

## Fairness: timing-node pinning decision (FLAGGED)
§3 calls for **same-node pinning** of the 4 tools per genome. On this shared
pool strict same-node pinning (`Machine == "<one host>"`) is **infeasible**:
- Only 15 AVX2 hosts exist, heterogeneous (6–102 cores); the pool carries
  500–700 idle jobs from other users.
- Pinning a genome's whole tool×rep block to one named host serialises ~20–40
  jobs onto a single contended machine and stalls BRAKER (16c/64 GB) indefinitely.

**Closest fair alternative adopted** (exactly the §-anticipated fallback, "pin to
one node class + report CPU-time"): pin **all four arms to one microarch
generation, x86_64-v3**, so wall-time is compared on identical-ISA hardware.
- The AVX2 pool splits into x86_64-v3 (9 hosts: scarcity-1,2,3,4,11,12,14,18,19)
  and x86_64-v4 (6 hosts). v3 chosen: more hosts (throughput, less contention)
  and more large-memory hosts able to seat BRAKER's 16c/64 GB slot
  (scarcity-11/12/18/19).
- Requirement on every timing sub:
  `requirements = (TARGET.has_avx2 =?= true) && (TARGET.Microarch =?= "x86_64-v3")`
  (matchmaking verified satisfiable: 27 slots match).
- **CPU-seconds** (`performance.json total_cpu_seconds`, captured by every
  wrapper via `/usr/bin/time -v`) is reported alongside wall-time as the
  co-tenant-contention-robust check. True single-tenant isolation is impossible
  on a partitionable-slot shared pool regardless of pinning, which is the
  independent reason CPU-seconds is reported.
- Accuracy runs are NOT pinned (deterministic, hardware-independent).

## Harness changes made on the cluster (MIRROR to canonical repo)
1. `jobs/myconote.sub`: queue block expanded 18→**42 rows** (14 genomes × 3
   reps); header comments updated. Backup: `jobs/myconote.sub.bak_pre42_<ts>`.
2. `scripts/gen_timing_subs.py`: NEW generator (expands myconote.sub; emits
   timing subs). Deterministic, idempotent.
3. NEW timing subs, x86_64-v3 pinned:
   - `jobs/{myconote,funannotate,maker,braker}_timing_s1.sub`   — 70 rows, n=5.
   - `jobs/{myconote,funannotate,maker,braker}_timing_prod.sub` — 140 rows, n=10.
   COSMETIC: the myconote timing subs inherited myconote.sub's "Submits 42 jobs"
   header comment; the queue block (70/140 rows) is authoritative.

## Pending edits BEFORE collect_results / make_figures (NOT yet made)
- `scripts/aggregate_metrics.py`: rep loop is hardcoded `for rep in [1,2,3]`;
  must glob `rep*` so n=5/n=10 timing dirs aggregate fully. (Accuracy 3-rep run
  is unaffected.)
- `scripts/make_figures.py`: `GENOMES` is the 6-genome list; expand to the 14
  panel. (It globs `rep*/metrics.json`, so it is already rep-count-flexible.)
- `collect_results.sh`: `REPS=(1 2 3)` + interactive `read -p` prompt — fine for
  the accuracy run; for timing dirs call `aggregate_metrics.py <dir>` directly
  (after the glob fix) to avoid the prompt.

## Submission log
### Phase 4 — Stage-1 timing (n=5, x86_64-v3 pinned) — SUBMITTED 2026-10-15 16:05
RESULTS_DIR = `results_timing_stage1/`
| tool | cluster | jobs |
|------|---------|------|
| myconote    | 146848 | 70 |
| funannotate | 146849 | 70 |
| maker       | 146850 | 70 |
| braker      | HELD — pending gate 146843 | (70 ready in braker_timing_s1.sub) |

Submit command (per arm):
`source configs/htcondor.conf && export RESULTS_DIR=$BENCHMARK_DIR/results_timing_stage1 && condor_submit -append "BENCHMARK_DIR=.." -append "DATA_DIR=.." -append "RESULTS_DIR=.." jobs/<tool>_timing_s1.sub`
Kill: `condor_rm 146848 146849 146850` (or `condor_rm narhmadey` for all).
Status at submit: 210 idle→running, 0 held; matchmaking satisfiable.

### BRAKER gate (146843, ccin header-fix re-run)
Still RUNNING on scarcity-18 at 16:10; no `braker.gtf`/`braker.gff3`/
`metrics.json`/`train.gb` produced yet. BRAKER (Stage-1 AND production) stays
HELD until 146843 exits 0 with braker.gff3/gtf + metrics.json + train.gb >0 genes.

### Phase 5 — production — NOT YET SUBMITTED
Gated on Stage-1 verification (per-tool×genome pass/fail). Will submit:
- n=10 timing (pinned): `*_timing_prod.sub` → `results_timing_prod/`.
- 3-rep accuracy: the 42-row `myconote.sub` + `funannotate.sub` / `maker.sub` /
  `braker.sub` → `results/` (or a dedicated accuracy dir).
