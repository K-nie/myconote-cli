#!/usr/bin/env python3
"""Per-test publication-grade benchmark figures + data files.

Produces ONE figure (PDF + 300 dpi PNG) and ONE data TSV per test, so each
comparison stands on its own and is work-matched (gene predictors and complete
annotation pipelines are never conflated on a speed axis):

  Test 1  accuracy             — loose/strict/nt gene F1, all 4 tools  (figure_accuracy)
  Test 2  prediction_speed     — gene-prediction wall-time, all 4 tools (figure_prediction_speed)
  Test 3  annotation_matched   — functional-annotation wall-time at MATCHED sources:
                                 MycoNote (no EggNog) vs funannotate, with MycoNote-full
                                 (with EggNog) shown for reference  (figure_annotation_matched)
  Test 4  completeness         — functional outputs each tool produces (table/TSV)
  Test 5  portability          — node-class availability (table/TSV)

Nature-style conventions: lowercase bold panel letters, sans-serif, Okabe-Ito
colourblind-safe palette, no titles on the figure, editable vector PDF + 300 dpi.

Inputs
------
--results-dir : benchmark results/ with <tool>/<genome>/rep<N>/metrics.json and
                timing logs (time_predict.log, time_annotate.log, time_maker.log),
                plus myconote_matched/<genome>/rep<N>/time_annotate.log for Test 3.
--out-dir     : where figures + *.tsv land (default: <results-dir>/../figures).
"""

import argparse
import glob
import json
import os
import re
import statistics as st
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TOOLS = ["myconote", "funannotate", "maker", "braker"]
TOOL_LABELS = {"myconote": "MycoNote", "funannotate": "funannotate", "maker": "MAKER", "braker": "BRAKER"}
TOOL_COLORS = {"myconote": "#0072B2", "funannotate": "#E69F00", "maker": "#009E73", "braker": "#CC79A7"}
GENOMES = ["sce", "cal", "ylp", "ani", "ncr", "cne", "klac", "kpha", "fgra", "afum", "umay", "ccin", "spom", "tmel"]
MM = 1.0 / 25.4
SINGLE_W, DOUBLE_W = 89 * MM, 183 * MM


def _rc():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "lines.linewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


def _save(fig, out_dir, stem):
    os.makedirs(out_dir, exist_ok=True)
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        p = os.path.join(out_dir, stem + ext)
        fig.savefig(p, bbox_inches="tight", **kw)
        print(f"  wrote {p}")
    plt.close(fig)


def _panel_letter(ax, letter):
    ax.text(-0.16, 1.04, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="right")


def _write_tsv(out_dir, name, header, rows):
    os.makedirs(out_dir, exist_ok=True)
    p = os.path.join(out_dir, name)
    with open(p, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    print(f"  wrote {p}")


# ── data loaders ─────────────────────────────────────────────────────────────
def _mins(path):
    try:
        t = open(path).read()
    except OSError:
        return None
    m = re.search(r"Elapsed.*?:\s*([0-9:.]+)", t)
    if not m:
        return None
    s = 0.0
    for x in m.group(1).split(":"):
        s = s * 60 + float(x)
    return s / 60.0


def load_accuracy(R):
    d = {m: {t: {g: [] for g in GENOMES} for t in TOOLS} for m in ("loose", "strict", "nt")}
    for t in TOOLS:
        for g in GENOMES:
            for mj in sorted(glob.glob(os.path.join(R, t, g, "rep*", "metrics.json"))):
                try:
                    j = json.load(open(mj)); gl = j["gene_level"]
                    d["loose"][t][g].append(gl["loose"]["F1"])
                    d["strict"][t][g].append(gl["strict"]["F1"])
                    d["nt"][t][g].append(j["nucleotide_level"]["F1"])
                except (KeyError, ValueError, OSError) as e:
                    print(f"  [skip] {mj}: {e}", file=sys.stderr)
    return d


def load_stage_minutes(R, subdir, fname):
    """All per-run minutes for one stage of one tool-dir (list)."""
    vals = []
    for log in glob.glob(os.path.join(R, subdir, "*", "rep*", fname)):
        v = _mins(log)
        if v is not None:
            vals.append(v)
    return vals


# ── Test 1: accuracy ─────────────────────────────────────────────────────────
def _grouped(ax, metric):
    present = [t for t in TOOLS if any(metric[t][g] for g in GENOMES)]
    n = len(present); w = 0.8 / max(n, 1); x = np.arange(len(GENOMES))
    for i, t in enumerate(present):
        means = [np.mean(metric[t][g]) if metric[t][g] else np.nan for g in GENOMES]
        sds = [np.std(metric[t][g]) if len(metric[t][g]) > 1 else 0 for g in GENOMES]
        ax.bar(x + (i - (n - 1) / 2) * w, means, w, yerr=sds, label=TOOL_LABELS[t],
               color=TOOL_COLORS[t], edgecolor="black", linewidth=0.4,
               error_kw={"elinewidth": 0.5, "capsize": 1.5, "capthick": 0.5})
    ax.set_xticks(x); ax.set_xticklabels(GENOMES); ax.set_ylim(0, 1.0)
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(length=2)


def test_accuracy(data, out_dir):
    _rc()
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(DOUBLE_W, DOUBLE_W * 0.34))
    _grouped(a, data["loose"]); a.set_ylabel("Gene-level loose F1"); _panel_letter(a, "a")
    _grouped(b, data["strict"]); b.set_ylabel("Gene-level strict F1"); _panel_letter(b, "b")
    _grouped(c, data["nt"]); c.set_ylabel("Nucleotide F1"); _panel_letter(c, "c")
    handles = [plt.Rectangle((0, 0), 1, 1, color=TOOL_COLORS[t], ec="black", lw=0.4) for t in TOOLS]
    fig.legend(handles, [TOOL_LABELS[t] for t in TOOLS], ncol=4, loc="lower center",
               frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    _save(fig, out_dir, "figure_accuracy")
    rows = []
    for t in TOOLS:
        for g in GENOMES:
            for m in ("loose", "strict", "nt"):
                vals = data[m][t][g]
                if vals:
                    rows.append([t, g, m, f"{np.mean(vals):.4f}", f"{(np.std(vals) if len(vals) > 1 else 0):.4f}", len(vals)])
    _write_tsv(out_dir, "accuracy.tsv", ["tool", "genome", "metric", "mean_F1", "sd", "n"], rows)


# ── Test 2: prediction speed (all 4; BRAKER timing unusable) ─────────────────
def test_prediction_speed(R, out_dir):
    _rc()
    # Each tool logs its prediction/whole-run wall-time to a tool-specific file
    # (BRAKER's /usr/bin/time wrapper writes time_braker.log; its internal
    # performance.json is all-zeros and is NOT used).
    pred = {
        "myconote": load_stage_minutes(R, "myconote", "time_predict.log"),
        "funannotate": load_stage_minutes(R, "funannotate", "time_predict.log"),
        "maker": load_stage_minutes(R, "maker", "time_maker.log"),
        "braker": load_stage_minutes(R, "braker", "time_braker.log"),
    }
    present = [t for t in TOOLS if pred[t]]
    fig, ax = plt.subplots(figsize=(SINGLE_W * 1.3, SINGLE_W * 1.0))
    means = [np.mean(pred[t]) for t in present]
    sds = [np.std(pred[t]) if len(pred[t]) > 1 else 0 for t in present]
    x = np.arange(len(present))
    ax.bar(x, means, 0.6, yerr=sds, color=[TOOL_COLORS[t] for t in present], edgecolor="black",
           linewidth=0.4, error_kw={"elinewidth": 0.5, "capsize": 1.5, "capthick": 0.5})
    ax.set_xticks(x); ax.set_xticklabels([TOOL_LABELS[t] for t in present])
    ax.set_ylabel("Gene-prediction wall-time (min)")
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(length=2)
    _save(fig, out_dir, "figure_prediction_speed")
    rows = [[t, f"{np.mean(pred[t]):.2f}", f"{(np.std(pred[t]) if len(pred[t]) > 1 else 0):.2f}", len(pred[t])] for t in present]
    for t in TOOLS:
        if not pred[t]:
            rows.append([t, "NA", "NA", "0  # timing log not found"])
    _write_tsv(out_dir, "prediction_speed.tsv", ["tool", "mean_min", "sd", "n"], rows)


# ── Test 3: functional-annotation speed at MATCHED sources ───────────────────
def test_annotation_matched(R, out_dir):
    _rc()
    matched = load_stage_minutes(R, "myconote_matched", "time_annotate.log")  # no EggNog
    fun = load_stage_minutes(R, "funannotate", "time_annotate.log")
    full = load_stage_minutes(R, "myconote", "time_annotate.log")  # with EggNog, for reference
    bars = []  # (label, color, vals, hatch)
    if matched:
        bars.append(("MycoNote\n(matched, no EggNog)", TOOL_COLORS["myconote"], matched, None))
    if fun:
        bars.append(("funannotate", TOOL_COLORS["funannotate"], fun, None))
    if full:
        bars.append(("MycoNote full\n(+ EggNog)", TOOL_COLORS["myconote"], full, "//"))
    fig, ax = plt.subplots(figsize=(SINGLE_W * 1.5, SINGLE_W * 1.05))
    x = np.arange(len(bars))
    for i, (lab, col, vals, hatch) in enumerate(bars):
        ax.bar(i, np.mean(vals), 0.6, yerr=(np.std(vals) if len(vals) > 1 else 0), color=col,
               edgecolor="black", linewidth=0.4, hatch=hatch,
               error_kw={"elinewidth": 0.5, "capsize": 1.5, "capthick": 0.5})
    ax.set_xticks(x); ax.set_xticklabels([b[0] for b in bars])
    ax.set_ylabel("Functional-annotation wall-time (min)")
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(length=2)
    _save(fig, out_dir, "figure_annotation_matched")
    rows = []
    for lab, _, vals, _ in bars:
        rows.append([lab.replace("\n", " "), f"{np.mean(vals):.2f}", f"{(np.std(vals) if len(vals) > 1 else 0):.2f}", len(vals)])
    _write_tsv(out_dir, "annotation_matched.tsv", ["arm", "mean_min", "sd", "n"], rows)


# ── Test 4 & 5: capability tables (static facts, written as TSV) ─────────────
def test_completeness(out_dir):
    # functional outputs each tool produces out of the box in this benchmark
    header = ["tool", "gene_structure", "swissprot", "pfam", "busco", "cazyme", "merops", "eggnog_orthology_GO_KEGG"]
    rows = [
        ["myconote", "yes", "yes", "yes", "yes", "yes", "yes", "yes"],
        ["funannotate", "yes", "yes", "yes", "yes", "yes", "yes", "no"],
        ["maker", "yes", "no", "no", "no", "no", "no", "no"],
        ["braker", "yes", "no", "no", "no", "no", "no", "no"],
    ]
    _write_tsv(out_dir, "completeness.tsv", header, rows)


def test_portability(out_dir):
    header = ["tool", "augustus_build", "requires_avx2", "runs_non_avx2_pool"]
    rows = [
        ["myconote", "3.3.3", "no", "yes"],
        ["maker", "n/a (SNAP/own)", "no", "yes"],
        ["funannotate", "3.5.0 (AVX2)", "yes", "no"],
        ["braker", "3.5.0 (AVX2)", "yes", "no"],
    ]
    _write_tsv(out_dir, "portability.tsv", header, rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args()
    out = a.out_dir or os.path.join(os.path.dirname(a.results_dir.rstrip("/")), "figures")
    print("Test 1 — accuracy"); test_accuracy(load_accuracy(a.results_dir), out)
    print("Test 2 — prediction speed"); test_prediction_speed(a.results_dir, out)
    print("Test 3 — annotation speed (matched sources)"); test_annotation_matched(a.results_dir, out)
    print("Test 4 — completeness"); test_completeness(out)
    print("Test 5 — portability"); test_portability(out)


if __name__ == "__main__":
    main()
