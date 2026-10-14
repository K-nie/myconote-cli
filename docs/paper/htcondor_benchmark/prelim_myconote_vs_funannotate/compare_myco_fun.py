#!/usr/bin/env python3
"""Preliminary MycoNote vs funannotate comparison from 1-rep smoke metrics.

MycoNote: results_smoke_expanded_v0715 (all 14, fresh, v0.7.15).
funannotate: 8 new genomes fresh (competitor smoke) + 6 originals rep1 from the
prior fair run (same v1.8.17, same references, same compare_annotations.py).
This is PRELIMINARY 1-rep data; the production matrix (n=10 timing, 3-rep
accuracy) supersedes it.
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = os.path.dirname(os.path.abspath(__file__))
CMP = os.path.join(S, "cmp")

# genome order: 6 originals then 8 new; label clades for readability
GENOMES = ["sce", "cal", "ylp", "ani", "ncr", "cne",
           "klac", "kpha", "fgra", "afum", "umay", "ccin", "spom", "tmel"]
NEW = {"klac", "kpha", "fgra", "afum", "umay", "ccin", "spom", "tmel"}


def load(path):
    with open(path) as fh:
        return json.load(fh)


def myco_paths(g):
    base = os.path.join(CMP, "myconote", g, "rep1")
    return os.path.join(base, "metrics.json"), os.path.join(base, "performance.json")


def fun_paths(g):
    if g in NEW:
        base = os.path.join(CMP, "funannotate_new", g, "rep1")
    else:
        base = os.path.join(CMP, "funannotate_orig", g, "rep1")
    return os.path.join(base, "metrics.json"), os.path.join(base, "performance.json")


def extract(mpath, ppath):
    m = load(mpath)
    p = load(ppath)
    st = p.get("stages", {})
    return {
        "n_genes": m["predicted_summary"]["n_genes"],
        "ref_genes": m["reference_summary"]["n_genes"],
        "gene_loose_F1": m["gene_level"]["loose"]["F1"],
        "gene_strict_F1": m["gene_level"]["strict"]["F1"],
        "exon_exact_F1": m["exon_level"]["exact_match"]["F1"],
        "nt_F1": m["nucleotide_level"]["F1"],
        "predict_min": st.get("predict", {}).get("wall_seconds", 0) / 60.0,
        "annotate_min": st.get("annotate", {}).get("wall_seconds", 0) / 60.0,
        "total_min": p.get("total_wall_seconds", 0) / 60.0,
        "peak_rss_gb": p.get("peak_rss_mb", 0) / 1024.0,
    }


rows = []
data = {}
for g in GENOMES:
    mm, mp = myco_paths(g)
    fm, fp = fun_paths(g)
    my = extract(mm, mp)
    fn = extract(fm, fp)
    data[g] = {"myconote": my, "funannotate": fn}
    rows.append((g, my, fn))

# ── Write comparison TSV ──────────────────────────────────────────────────
tsv = os.path.join(S, "comparison_myconote_vs_funannotate.tsv")
cols = ["n_genes", "gene_loose_F1", "gene_strict_F1", "exon_exact_F1", "nt_F1",
        "predict_min", "annotate_min", "total_min", "peak_rss_gb"]
with open(tsv, "w") as fh:
    hdr = ["genome", "ref_genes", "tool"] + cols
    fh.write("\t".join(hdr) + "\n")
    for g, my, fn in rows:
        for tool, d in (("myconote", my), ("funannotate", fn)):
            vals = [g, str(my["ref_genes"]), tool]
            for c in cols:
                v = d[c]
                vals.append(f"{v:.4f}" if isinstance(v, float) else str(v))
            fh.write("\t".join(vals) + "\n")

# ── Console table ──────────────────────────────────────────────────────────
print(f"{'genome':6} {'gLoose F1':>18} {'gStrict F1':>18} {'nt F1':>16} {'total min':>18}")
print(f"{'':6} {'myco / fun':>18} {'myco / fun':>18} {'myco / fun':>16} {'myco / fun':>18}")
for g, my, fn in rows:
    tag = "*" if g in NEW else " "
    print(f"{g+tag:6} "
          f"{my['gene_loose_F1']:.3f} / {fn['gene_loose_F1']:.3f}    "
          f"{my['gene_strict_F1']:.3f} / {fn['gene_strict_F1']:.3f}    "
          f"{my['nt_F1']:.3f} / {fn['nt_F1']:.3f}   "
          f"{my['total_min']:6.1f} / {fn['total_min']:6.1f}")

# win tallies
def tally(key, higher_better=True):
    mw = fw = tie = 0
    for g, my, fn in rows:
        a, b = my[key], fn[key]
        if abs(a - b) < 1e-9:
            tie += 1
        elif (a > b) == higher_better:
            mw += 1
        else:
            fw += 1
    return mw, fw, tie

print("\nWin tally (MycoNote / funannotate / tie) across 14 genomes:")
for k, hb in [("gene_loose_F1", True), ("gene_strict_F1", True),
              ("exon_exact_F1", True), ("nt_F1", True), ("total_min", False)]:
    mw, fw, tie = tally(k, hb)
    better = "higher" if hb else "lower"
    print(f"  {k:16} ({better} better): {mw} / {fw} / {tie}")

# ── Figure (Nature-style: sans-serif, lowercase bold panel letters, no titles) ─
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.6,
})
OK = {"myconote": "#0072B2", "funannotate": "#E69F00"}  # Okabe-Ito blue / orange
import numpy as np
x = np.arange(len(GENOMES))
w = 0.4
labels = [g + ("*" if g in NEW else "") for g in GENOMES]

panels = [
    ("gene_loose_F1", "Gene F1 (loose)", False),
    ("gene_strict_F1", "Gene F1 (strict)", False),
    ("nt_F1", "Nucleotide F1", False),
    ("total_min", "Total wall-clock (min)", True),
]
fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
for ax, (key, ylab, is_time) in zip(axes.flat, panels):
    my = [data[g]["myconote"][key] for g in GENOMES]
    fn = [data[g]["funannotate"][key] for g in GENOMES]
    ax.bar(x - w / 2, my, w, label="MycoNote", color=OK["myconote"])
    ax.bar(x + w / 2, fn, w, label="funannotate", color=OK["funannotate"])
    ax.set_ylabel(ylab)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=5.5)
    if not is_time:
        ax.set_ylim(0, 1)
    ax.tick_params(width=0.6, length=2.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
for ax, letter in zip(axes.flat, "abcd"):
    ax.text(-0.14, 1.04, letter, transform=ax.transAxes,
            fontweight="bold", fontsize=10, va="top", ha="left")
axes.flat[0].legend(frameon=False, fontsize=6, loc="lower right")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(S, f"figure_myconote_vs_funannotate.{ext}"),
                dpi=300, bbox_inches="tight")
print(f"\nWrote: {tsv}")
print("Wrote: figure_myconote_vs_funannotate.{pdf,png}")
print("(* = new expanded-panel genome; funannotate originals are prior-run rep1)")
