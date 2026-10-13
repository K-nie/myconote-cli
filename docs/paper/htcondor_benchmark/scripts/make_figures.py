#!/usr/bin/env python3
"""Publication-grade benchmark figures: MycoNote vs funannotate vs MAKER vs BRAKER.

Reads the per-tool/genome/rep ``metrics.json`` tree produced by the HTCondor
benchmark and renders a multi-panel accuracy + speed figure.

Nature-style conventions (fixed, not optional):
  * lowercase bold panel letters (a, b, c, d), top-left of each panel
  * sans-serif in-figure text; no titles rendered on the figure itself
  * colourblind-safe palette (Okabe-Ito); line weights >= 0.5 pt
  * composed to double-column width (183 mm)
  * output: editable vector PDF + 300 dpi PNG

Inputs
------
--results-dir : the benchmark `results/` dir containing <tool>/<genome>/rep<N>/metrics.json
                and (optional) time_predict.log / time_mask.log / time_sort.log per run.
--out-dir     : where to write figure_benchmark.pdf / .png (default: <results-dir>/../figures).

Usage
-----
    python3 make_figures.py --results-dir ../results --out-dir ../figures

Reps are averaged per (tool, genome); missing reps/tools are handled gracefully
(a tool with no data for a genome is simply absent from that group, not zero).
"""

import argparse
import glob
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Fixed conventions ────────────────────────────────────────────────────────
TOOLS = ["myconote", "funannotate", "maker", "braker"]
TOOL_LABELS = {
    "myconote": "MycoNote",
    "funannotate": "funannotate",
    "maker": "MAKER",
    "braker": "BRAKER",
}
# Okabe-Ito colourblind-safe palette (no red/green as the sole distinction).
TOOL_COLORS = {
    "myconote": "#0072B2",  # blue
    "funannotate": "#E69F00",  # orange
    "maker": "#009E73",  # bluish green
    "braker": "#CC79A7",  # reddish purple
}
# Genome display order (Saccharomycotina → Pezizomycotina → Basidiomycota).
GENOMES = ["sce", "cal", "ylp", "ani", "ncr", "cne"]

MM = 1.0 / 25.4  # mm → inches
DOUBLE_COL_W = 183 * MM


def _rc():
    """Nature-style rcParams: sans-serif, small legible text, thin clean lines."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "lines.linewidth": 0.8,
            "pdf.fonttype": 42,  # editable text in the PDF
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


# ── Data loading ─────────────────────────────────────────────────────────────
def load_metrics(results_dir):
    """Return nested dict: data[metric][tool][genome] = list of per-rep values.

    metric ∈ {loose, strict, nt}. Only genomes/tools with real metrics.json
    files contribute; everything else is simply absent.
    """
    data = {m: {t: {g: [] for g in GENOMES} for t in TOOLS} for m in ("loose", "strict", "nt")}
    for t in TOOLS:
        for g in GENOMES:
            for m in sorted(glob.glob(os.path.join(results_dir, t, g, "rep*", "metrics.json"))):
                try:
                    d = json.load(open(m))
                    gl = d["gene_level"]
                    data["loose"][t][g].append(gl["loose"]["F1"])
                    data["strict"][t][g].append(gl["strict"]["F1"])
                    data["nt"][t][g].append(d["nucleotide_level"]["F1"])
                except (KeyError, ValueError, OSError) as e:
                    print(f"  [skip] {m}: {e}", file=sys.stderr)
    return data


def _parse_time_log(path):
    """Minutes from a `/usr/bin/time -v` log (Elapsed wall-clock), else None."""
    try:
        txt = open(path).read()
    except OSError:
        return None
    # GNU time: "Elapsed (wall clock) time (h:mm:ss or m:ss): 1:23:45"
    m = re.search(r"Elapsed.*?:\s*([0-9:.]+)\s*$", txt, re.MULTILINE)
    if not m:
        return None
    parts = [float(x) for x in m.group(1).split(":")]
    secs = 0.0
    for p in parts:
        secs = secs * 60 + p
    return secs / 60.0


def load_timing(results_dir):
    """Return timing[tool][stage] = list of per-run minutes, for the predict stage.

    Speed is reported for the gene-prediction stage (time_predict.log), the one
    directly comparable across all four tools. Missing logs are skipped.
    """
    timing = {t: [] for t in TOOLS}
    for t in TOOLS:
        for log in glob.glob(os.path.join(results_dir, t, "*", "rep*", "time_predict.log")):
            mins = _parse_time_log(log)
            if mins is not None:
                timing[t].append(mins)
    return timing


# ── Plotting ─────────────────────────────────────────────────────────────────
def _grouped_bars(ax, data_metric, ylabel):
    """Grouped bars: x = genome, bars = tools (mean over reps, sd whiskers)."""
    present_tools = [t for t in TOOLS if any(data_metric[t][g] for g in GENOMES)]
    n = len(present_tools)
    width = 0.8 / max(n, 1)
    x = np.arange(len(GENOMES))
    for i, t in enumerate(present_tools):
        means, sds = [], []
        for g in GENOMES:
            vals = data_metric[t][g]
            means.append(np.mean(vals) if vals else np.nan)
            sds.append(np.std(vals) if len(vals) > 1 else 0.0)
        ax.bar(
            x + (i - (n - 1) / 2) * width,
            means,
            width,
            yerr=sds,
            label=TOOL_LABELS[t],
            color=TOOL_COLORS[t],
            edgecolor="black",
            linewidth=0.4,
            error_kw={"elinewidth": 0.5, "capsize": 1.5, "capthick": 0.5},
        )
    ax.set_xticks(x)
    ax.set_xticklabels([g for g in GENOMES])
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 1.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2)


def _speed_bars(ax, timing):
    """Per-tool gene-prediction wall-time (mean ± sd minutes). Lower is better."""
    present = [t for t in TOOLS if timing[t]]
    if not present:
        ax.text(0.5, 0.5, "prediction timing not available\n(run after jobs finish)",
                ha="center", va="center", transform=ax.transAxes, fontsize=6.5, color="#555")
        ax.set_axis_off()
        return
    means = [np.mean(timing[t]) for t in present]
    sds = [np.std(timing[t]) if len(timing[t]) > 1 else 0.0 for t in present]
    x = np.arange(len(present))
    ax.bar(x, means, 0.6, yerr=sds, color=[TOOL_COLORS[t] for t in present],
           edgecolor="black", linewidth=0.4,
           error_kw={"elinewidth": 0.5, "capsize": 1.5, "capthick": 0.5})
    ax.set_xticks(x)
    ax.set_xticklabels([TOOL_LABELS[t] for t in present], rotation=20, ha="right")
    ax.set_ylabel("Gene-prediction wall-time (min)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2)


def _panel_letter(ax, letter):
    ax.text(-0.14, 1.04, letter, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="right")


def make_figure(data, timing, out_dir):
    _rc()
    os.makedirs(out_dir, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL_W, DOUBLE_COL_W * 0.72))
    (a, b), (c, d) = axes

    _grouped_bars(a, data["loose"], "Gene-level loose F1")
    _panel_letter(a, "a")
    _grouped_bars(b, data["strict"], "Gene-level strict F1")
    _panel_letter(b, "b")
    _grouped_bars(c, data["nt"], "Nucleotide F1")
    _panel_letter(c, "c")
    _speed_bars(d, timing)
    _panel_letter(d, "d")

    # Single shared legend (tools), placed below the panels — no title on figure.
    handles = [plt.Rectangle((0, 0), 1, 1, color=TOOL_COLORS[t], ec="black", lw=0.4)
               for t in TOOLS]
    fig.legend(handles, [TOOL_LABELS[t] for t in TOOLS], ncol=4,
               loc="lower center", frameon=False, bbox_to_anchor=(0.5, 0.015))

    fig.tight_layout(rect=[0, 0.07, 1, 1])
    pdf = os.path.join(out_dir, "figure_benchmark.pdf")
    png = os.path.join(out_dir, "figure_benchmark.png")
    fig.savefig(pdf, bbox_inches="tight")  # vector
    fig.savefig(png, dpi=300, bbox_inches="tight")
    print(f"  wrote {pdf}")
    print(f"  wrote {png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()
    out_dir = args.out_dir or os.path.join(os.path.dirname(args.results_dir.rstrip("/")), "figures")
    data = load_metrics(args.results_dir)
    timing = load_timing(args.results_dir)
    make_figure(data, timing, out_dir)


if __name__ == "__main__":
    main()
