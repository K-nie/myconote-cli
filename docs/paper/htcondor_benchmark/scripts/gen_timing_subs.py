#!/usr/bin/env python3
"""Generate timing subs + expand myconote.sub for the 14-genome production matrix.

Run from the htcondor_benchmark/ directory. Idempotent: rewrites target files.

- Expands jobs/myconote.sub queue block to 14 genomes x 3 reps (42 rows).
- Creates jobs/{tool}_timing_s1.sub  (n=5 timing, x86_64-v3 pinned).
- Creates jobs/{tool}_timing_prod.sub (n=10 timing, x86_64-v3 pinned).

Timing fairness: strict same-node pinning (Machine == one named host) is
infeasible on this shared, heavily-loaded pool (one named host would serialise
the whole tool x rep block behind 500 competing idle jobs and stall BRAKER's
16c/64GB slot indefinitely). The closest fair alternative is to pin every arm
to ONE microarch generation so wall-time is compared on identical ISA hardware;
total_cpu_seconds (already captured by every wrapper via /usr/bin/time -v) is
the co-tenant-contention-robust check reported alongside wall-time.
"""
import re
import sys
from pathlib import Path

GENOMES = ["sce", "cal", "ylp", "ani", "ncr", "cne", "klac", "kpha",
           "fgra", "afum", "umay", "ccin", "spom", "tmel"]

V3_REQ = '(TARGET.has_avx2 =?= true) && (TARGET.Microarch =?= "x86_64-v3")'

JOBS = Path("jobs")
TOOLS = ["myconote", "funannotate", "maker", "braker"]


def queue_block(reps):
    lines = ["queue genome_id, rep from ("]
    for g in GENOMES:
        for r in range(1, reps + 1):
            lines.append(f"    {g} {r}")
    lines.append(")")
    return "\n".join(lines) + "\n"


def split_preamble(text):
    """Return everything before the 'queue genome_id' statement."""
    m = re.search(r"^queue\s+genome_id", text, re.MULTILINE)
    if not m:
        raise SystemExit("no 'queue genome_id' statement found")
    return text[:m.start()]


def swap_requirements(preamble, req):
    """Replace the active (non-comment) requirements line with `req`."""
    out, replaced = [], False
    for line in preamble.splitlines():
        if re.match(r"^\s*requirements\s*=", line) and not line.lstrip().startswith("#"):
            out.append(f"requirements = {req}")
            replaced = True
        else:
            out.append(line)
    if not replaced:
        # No active requirements line (shouldn't happen for these subs) -> add one.
        out.append(f"requirements = {req}")
    return "\n".join(out).rstrip() + "\n\n"


def main():
    # (A) Expand the accuracy/production myconote.sub to 14 x 3 = 42 rows.
    myco = JOBS / "myconote.sub"
    t = myco.read_text()
    pre = split_preamble(t)
    pre = pre.replace("Submits 18 jobs: 6 genomes x 3 replicates",
                      "Submits 42 jobs: 14 genomes x 3 replicates")
    pre = re.sub(r"# 6 genomes x 3 replicates = 18 jobs total\.",
                 "# 14 genomes x 3 replicates = 42 jobs total.", pre)
    myco.write_text(pre + queue_block(3))
    print(f"[A] myconote.sub expanded -> 42 rows "
          f"({(myco.read_text().count(chr(10)+'    ') )} queue entries approx)")

    # (B) Timing subs for all four tools.
    for tool in TOOLS:
        src = JOBS / f"{tool}.sub"
        pre = swap_requirements(split_preamble(src.read_text()), V3_REQ)
        hdr = (f"# TIMING sub (auto-generated) — {tool}, x86_64-v3 pinned.\n"
               f"# Fair wall-time: all arms pinned to one microarch generation\n"
               f"# (strict same-node pinning infeasible on this shared pool).\n"
               f"# CPU-seconds (performance.json total_cpu_seconds) is the\n"
               f"# contention-robust check. Submit with RESULTS_DIR pointed at a\n"
               f"# dedicated timing dir via condor_submit -append.\n\n")
        for tag, reps in (("s1", 5), ("prod", 10)):
            dst = JOBS / f"{tool}_timing_{tag}.sub"
            body = hdr + pre + queue_block(reps)
            dst.write_text(body)
            nrows = len(GENOMES) * reps
            print(f"[B] wrote {dst}  ({nrows} rows, n={reps})")


if __name__ == "__main__":
    if not JOBS.is_dir():
        sys.exit("run from htcondor_benchmark/ (no jobs/ dir here)")
    main()
