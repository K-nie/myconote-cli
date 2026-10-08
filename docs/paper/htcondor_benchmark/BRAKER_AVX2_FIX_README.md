# BRAKER arm — AVX2 scheduling fix (2026-10-07)

## What this documents and why
The BRAKER3 arm of the 4-tool gene-annotation benchmark (MycoNote vs MAKER vs
funannotate vs BRAKER) failed on the GLBRC HTCondor pool: a full 18-job run
(cluster 145888) produced `braker.gtf` for only the jobs that happened to land
on AVX2-capable execute nodes. This README records the root cause and the fix so
the result is reproducible and the Methods section is accurate.

## Root cause (diagnosed from cluster 145888 logs)
`augustus` 3.5.0 in the `myconote_braker3` conda env is compiled with AVX2
(1912 `ymm` instruction references; 0 `zmm`/AVX-512 — confirmed via
`objdump -d .../bin/augustus | grep -c '%ymm'`). On execute nodes WITHOUT AVX2
(microarch x86_64-v2, e.g. scarcity-10) augustus aborts with SIGILL —
"Illegal instruction (core dumped)" — at braker.pl line 2551 (the config-time
augustus self-test). braker.pl then reports "augustus not executable on this
machine" and exits 1 with no braker.gtf.

condor_history 145888 showed procs 3–17 (15 jobs) were all matched to
scarcity-10 and all exited 1; the jobs that landed on AVX2 nodes
(scarcity-20 = x86_64-v4, scarcity-3 = x86_64-v3) either scored (sce rep2) or ran
the full pipeline fine (GeneMark-EP completed, Augustus training proceeded).
So the failure is purely a CPU-instruction-set vs node-placement mismatch; it is
NOT a GeneMark-EP failure (GeneMark-EP succeeds on AVX2 nodes). The earlier
LD_LIBRARY_PATH fix addressed a different, missing-shared-lib case and does not
help a SIGILL.

Node CPU features (from condor_status -af Machine Microarch has_avx2):
  scarcity-10  x86_64-v2  has_avx2 = false  -> augustus SIGILL (all jobs failed)
  scarcity-3   x86_64-v3  has_avx2 = true   -> augustus OK, full pipeline runs
  scarcity-20  x86_64-v4  has_avx2 = true   -> augustus OK, sce rep2 scored

## The fix
1. jobs/braker.sub: add `requirements = (TARGET.has_avx2 =?= true)` so HTCondor
   only matches BRAKER jobs to AVX2-capable execute nodes (188 slots / 15
   machines available at fix time). The =?= operator treats an undefined ad as
   false, so nodes that do not publish the feature are safely excluded.
2. jobs/run_braker.sh: add a loud AVX2 preflight (grep -qw avx2 /proc/cpuinfo)
   that exits early with an actionable message if a job is ever mis-scheduled to
   a non-AVX2 node — so the failure is legible instead of a buried SIGILL.

Neither change alters the annotation algorithm or BRAKER parameters, so the tool
comparison stays fair. The only change is WHERE jobs are allowed to run.

## How to reproduce the run
    cd ~/myconote-cli/docs/paper/htcondor_benchmark
    source configs/htcondor.conf        # sets BENCHMARK_DIR/DATA_DIR/RESULTS_DIR, BRAKER_ENV=myconote_braker3
    BD=$BENCHMARK_DIR DD=$DATA_DIR RD=$RESULTS_DIR
    # smoke (2 previously-failed jobs):
    condor_submit -append "BENCHMARK_DIR = $BD" -append "DATA_DIR = $DD" -append "RESULTS_DIR = $RD" jobs/braker_smoke.sub
    # full 18 (clean uniform run), once smoke validates:
    condor_submit -append "BENCHMARK_DIR = $BD" -append "DATA_DIR = $DD" -append "RESULTS_DIR = $RD" jobs/braker.sub

BRAKER call (per job, unchanged):
    braker.pl --genome=genome.fa --prot_seq=OrthoDB/Fungi.fa --species=braker_<id>_rep<k> --workingdir=<out> --threads=16
GTF to GFF3 fallback: gtf2gff.pl; scoring: scripts/compare_annotations.py with a
seqid --rename-table built from the genome FASTA headers; "scored" = metrics.json.

## Inputs
- Genomes/refs: data/{sce,cal,ylp,ani,ncr,cne}/{genome.fa,reference.gff3}
- Protein evidence: ~/.myconote/dbs/orthodb/Fungi.fa (OrthoDB odb11 Fungi partition)
- 6 genomes x 3 replicates = 18 jobs

## Outputs
- results/braker/<id>/rep<k>/{braker.gtf,braker.gff3,metrics.json,performance.json}
- Target: 18/18 metrics.json

## Pinned versions
- braker.pl 3.0.8 (conda env myconote_braker3)
- augustus 3.5.0 (AVX2 build, dated Aug 20 2025)
- GeneMark-ES/EP suite: ~/.myconote/tools/gmes_linux_64_4, key ~/.gm_key
- HTCondor pool: scarcity (scarcity-ap-1 access point)

## Caveats / decisions
- AVX2 gate reduces the eligible node set but does not bias the comparison; it is
  a hard binary-compatibility requirement, not a performance tuning.
- Wall-time / RSS are measured with /usr/bin/time -v and are comparable across
  nodes only within the AVX2 class; note the node microarch mix in Methods if
  wall-time is reported.
