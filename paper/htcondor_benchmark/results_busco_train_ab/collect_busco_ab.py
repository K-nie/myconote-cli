#!/usr/bin/env python3
"""Aggregate the BUSCO-train A/B results into per-arm means + per-genome tables.

Reads <root>/<arm>/<genome>/{metrics.json, train_verdict.txt, time_predict.log}
for arms stock|self|busco over the 14-genome panel. Prints:
  (a) per-arm mean gene F1 (loose, strict) + mean predicted n_genes (n reps used),
  (b) per-genome trained-vs-fell-back table for the busco arm,
  (c) busco-vs-stock delta + predict wall-time (busco vs stock).
All means exclude genomes whose metrics.json is missing/unparseable (reported).
"""
import json, os, re, sys, statistics

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
GENOMES = "sce cal ylp ani ncr cne klac kpha fgra afum umay ccin spom tmel".split()
ARMS = ["stock", "self", "busco"]


def load_metrics(arm, g):
    p = os.path.join(ROOT, arm, g, "metrics.json")
    if not os.path.isfile(p):
        return None
    try:
        d = json.load(open(p))
        return {
            "loose": d["gene_level"]["loose"]["F1"],
            "strict": d["gene_level"]["strict"]["F1"],
            "n_genes": d["predicted_summary"]["n_genes"],
        }
    except Exception as e:
        print(f"  WARN parse {p}: {e}", file=sys.stderr)
        return None


def parse_time_log(arm, g):
    """Return (wall_s, max_rss_mb) from /usr/bin/time -v output."""
    p = os.path.join(ROOT, arm, g, "time_predict.log")
    wall = rss = None
    if not os.path.isfile(p):
        return (None, None)
    for line in open(p, errors="replace"):
        line = line.strip()
        if "Elapsed" in line and "wall clock" in line:
            tok = line.split()[-1]
            try:
                nums = [float(x) for x in tok.split(":")]
            except ValueError:
                continue
            if len(nums) == 3:
                wall = nums[0] * 3600 + nums[1] * 60 + nums[2]
            elif len(nums) == 2:
                wall = nums[0] * 60 + nums[1]
            elif len(nums) == 1:
                wall = nums[0]
        elif "Maximum resident set size" in line:
            try:
                rss = int(line.split(":")[1].strip()) / 1024.0  # KB -> MB (Linux)
            except Exception:
                pass
    return (wall, rss)


def verdict(arm, g):
    p = os.path.join(ROOT, arm, g, "train_verdict.txt")
    if not os.path.isfile(p):
        return "NO_VERDICT_FILE"
    return open(p).read().strip().replace("\n", " | ")


def mean(xs):
    xs = [x for x in xs if x is not None]
    return round(statistics.mean(xs), 4) if xs else None


data = {arm: {g: load_metrics(arm, g) for g in GENOMES} for arm in ARMS}
times = {arm: {g: parse_time_log(arm, g) for g in GENOMES} for arm in ARMS}

print("=" * 72)
print("(a) PER-ARM MEANS  (baseline to beat: loose 0.899 / strict 0.609)")
print("=" * 72)
print(f"{'arm':6s} {'n':>3s} {'meanF1_loose':>13s} {'meanF1_strict':>14s} {'mean_n_genes':>13s} {'mean_predict_wall_s':>20s}")
for arm in ARMS:
    present = [g for g in GENOMES if data[arm][g]]
    ml = mean([data[arm][g]["loose"] for g in present])
    ms = mean([data[arm][g]["strict"] for g in present])
    mn = mean([data[arm][g]["n_genes"] for g in present])
    mw = mean([times[arm][g][0] for g in GENOMES])
    print(f"{arm:6s} {len(present):>3d} {str(ml):>13s} {str(ms):>14s} {str(round(mn,1) if mn else mn):>13s} {str(round(mw,1) if mw else mw):>20s}")
    missing = [g for g in GENOMES if not data[arm][g]]
    if missing:
        print(f"       (missing/unscored: {','.join(missing)})")

print()
print("=" * 72)
print("(b) BUSCO ARM — per-genome trained vs fell-back")
print("=" * 72)
print(f"{'genome':7s} {'verdict':9s} {'detail'}")
n_trained = n_fellback = n_unknown = 0
for g in GENOMES:
    v = verdict("busco", g)
    tag = "TRAINED" if "verdict=TRAINED" in v else ("FELL_BACK" if "FELL_BACK" in v else "?")
    if tag == "TRAINED":
        n_trained += 1
    elif tag == "FELL_BACK":
        n_fellback += 1
    else:
        n_unknown += 1
    print(f"{g:7s} {tag:9s} {v}")
print(f"\nbusco summary: trained={n_trained} fell_back={n_fellback} unknown/missing={n_unknown} / {len(GENOMES)}")

print()
print("=" * 72)
print("(c) BUSCO vs STOCK — per-genome delta (F1) + predict wall-time")
print("=" * 72)
print(f"{'genome':7s} {'stock_loose':>11s} {'busco_loose':>11s} {'dL':>7s} {'stock_strict':>12s} {'busco_strict':>12s} {'dS':>7s} {'stock_wall':>10s} {'busco_wall':>10s}")
dL_all, dS_all = [], []
for g in GENOMES:
    s, b = data["stock"][g], data["busco"][g]
    sw, bw = times["stock"][g][0], times["busco"][g][0]
    if s and b:
        dL = round(b["loose"] - s["loose"], 4); dS = round(b["strict"] - s["strict"], 4)
        dL_all.append(dL); dS_all.append(dS)
    else:
        dL = dS = None
    def f(x): return f"{x:.4f}" if isinstance(x, float) else str(x)
    def w(x): return f"{x:.0f}" if isinstance(x, float) else str(x)
    print(f"{g:7s} {f(s['loose']) if s else 'NA':>11s} {f(b['loose']) if b else 'NA':>11s} {str(dL):>7s} "
          f"{f(s['strict']) if s else 'NA':>12s} {f(b['strict']) if b else 'NA':>12s} {str(dS):>7s} "
          f"{w(sw):>10s} {w(bw):>10s}")
if dL_all:
    print(f"\nmean delta busco-stock: loose {round(statistics.mean(dL_all),4)}  strict {round(statistics.mean(dS_all),4)}  (n={len(dL_all)} paired)")
