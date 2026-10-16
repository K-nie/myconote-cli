#!/usr/bin/env python3
"""Render the 12 non-redundant benchmark figures from the aggregated TSVs.
Run locally (needs matplotlib). Nature conventions: lowercase bold panel letters,
sans-serif, NO figure titles, Okabe-Ito colourblind-safe tool colours (consistent
across all figures), PDF + 300 dpi PNG. Each figure's descriptive title lives in
figures/README.md, never on the image.
"""
import csv
import os
import sys
import math
import statistics as st

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TSV = sys.argv[1]
OUT = sys.argv[2]
os.makedirs(OUT, exist_ok=True)

GENOMES = ["sce", "cal", "ylp", "ani", "ncr", "cne", "klac", "kpha",
           "fgra", "afum", "umay", "ccin", "spom", "tmel"]
TOOLS = ["myconote", "funannotate", "maker", "braker"]
TLAB = {"myconote": "MycoNote", "funannotate": "funannotate", "maker": "MAKER", "braker": "BRAKER"}
TCOL = {"myconote": "#0072B2", "funannotate": "#E69F00", "maker": "#009E73", "braker": "#CC79A7"}
CLADES = ["Saccharomycotina", "other-Ascomycota", "Basidiomycota"]
OKABE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442"]

MM = 1 / 25.4
SINGLE, DOUBLE = 89 * MM, 183 * MM

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "lines.linewidth": 0.9, "patch.linewidth": 0.4,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


def load(name):
    with open(os.path.join(TSV, name)) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fnum(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="right")


def save(fig, stem):
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        p = os.path.join(OUT, stem + ext)
        fig.savefig(p, bbox_inches="tight", **kw)
        print("wrote", p)
    plt.close(fig)


def tool_legend(fig, tools=TOOLS, y=-0.02):
    h = [plt.Rectangle((0, 0), 1, 1, color=TCOL[t], ec="black", lw=0.4) for t in tools]
    fig.legend(h, [TLAB[t] for t in tools], ncol=len(tools), loc="lower center",
               bbox_to_anchor=(0.5, y), frameon=False)


def grouped_bars(ax, rows, valkey, logy=False, tools=TOOLS):
    """x=genomes, grouped by tool. rows: list of dicts with genome/tool/<valkey>."""
    idx = {(r["genome"], r["tool"]): fnum(r.get(valkey)) for r in rows}
    n = len(tools)
    w = 0.8 / n
    x = np.arange(len(GENOMES))
    for i, t in enumerate(tools):
        vals = [idx.get((g, t)) for g in GENOMES]
        xs = x + (i - (n - 1) / 2) * w
        vals_plot = [v if (v is not None and (not logy or v > 0)) else np.nan for v in vals]
        ax.bar(xs, vals_plot, w, color=TCOL[t], edgecolor="black", linewidth=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels(GENOMES, rotation=90)
    if logy:
        ax.set_yscale("log")


# ---------- Fig 1: accuracy_by_genome (4 panels) ----------
def fig1():
    acc = load("accuracy.tsv")
    metrics = [("gene_f1_loose", "gene F1 (loose)"), ("gene_f1_strict", "gene F1 (strict)"),
               ("exon_f1_exact", "exon F1 (exact)"), ("nt_f1", "nucleotide F1")]
    fig, axs = plt.subplots(2, 2, figsize=(DOUBLE, DOUBLE * 0.62))
    for ax, (k, lab), L in zip(axs.flat, metrics, "abcd"):
        grouped_bars(ax, acc, k)
        ax.set_ylabel(lab)
        ax.set_ylim(0, 1.0)
        panel(ax, L)
    tool_legend(fig, y=-0.04)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, "fig1_accuracy_by_genome")


# ---------- Fig 2: accuracy_distribution (4 boxplot panels) ----------
def fig2():
    acc = load("accuracy.tsv")
    metrics = [("gene_f1_loose", "gene F1 (loose)"), ("gene_f1_strict", "gene F1 (strict)"),
               ("exon_f1_exact", "exon F1 (exact)"), ("nt_f1", "nucleotide F1")]
    fig, axs = plt.subplots(2, 2, figsize=(DOUBLE, DOUBLE * 0.62))
    for ax, (k, lab), L in zip(axs.flat, metrics, "abcd"):
        data = [[fnum(r[k]) for r in acc if r["tool"] == t and fnum(r[k]) is not None] for t in TOOLS]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showmeans=False,
                        medianprops=dict(color="black", linewidth=1.0),
                        flierprops=dict(marker="o", markersize=2, markerfacecolor="gray",
                                        markeredgecolor="none"))
        for patch, t in zip(bp["boxes"], TOOLS):
            patch.set_facecolor(TCOL[t])
            patch.set_edgecolor("black")
            patch.set_alpha(0.85)
        ax.set_xticks(range(1, len(TOOLS) + 1))
        ax.set_xticklabels([TLAB[t] for t in TOOLS], rotation=20, ha="right")
        ax.set_ylabel(lab)
        panel(ax, L)
    fig.tight_layout()
    save(fig, "fig2_accuracy_distribution")


# ---------- Fig 3: runtime_by_genome (wall, log y) ----------
def fig3():
    tw = load("timing_wall.tsv")
    fig, ax = plt.subplots(figsize=(DOUBLE, DOUBLE * 0.42))
    grouped_bars(ax, tw, "wall_s_median", logy=True)
    ax.set_ylabel("wall-clock time (s, log scale)")
    tool_legend(fig, y=-0.08)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, "fig3_runtime_by_genome")


# ---------- Fig 4: runtime_decomposition (predict vs annotate stacked) ----------
def fig4():
    tw = load("timing_wall.tsv")
    d = {(r["genome"], r["tool"]): r for r in tw}
    fig, ax = plt.subplots(figsize=(DOUBLE, DOUBLE * 0.45))
    n = len(TOOLS)
    w = 0.8 / n
    x = np.arange(len(GENOMES))
    for i, t in enumerate(TOOLS):
        xs = x + (i - (n - 1) / 2) * w
        pred, ann, tot, mono = [], [], [], []
        for g in GENOMES:
            r = d.get((g, t), {})
            p = fnum(r.get("predict_wall_s_median"))
            a = fnum(r.get("annotate_wall_s_median"))
            w_ = fnum(r.get("wall_s_median"))
            if p is not None and a is not None:
                pred.append(p); ann.append(a); tot.append(np.nan); mono.append(np.nan)
            else:
                pred.append(np.nan); ann.append(np.nan); tot.append(w_); mono.append(w_)
        c = TCOL[t]
        ax.bar(xs, pred, w, color=c, edgecolor="black", linewidth=0.3)
        ax.bar(xs, ann, w, bottom=np.nan_to_num(pred), color=c, alpha=0.5,
               edgecolor="black", linewidth=0.3)
        # monolithic tools (maker, braker): total bar, hatched, no stage split
        ax.bar(xs, mono, w, color=c, edgecolor="black", linewidth=0.3, hatch="///")
    ax.set_xticks(x); ax.set_xticklabels(GENOMES, rotation=90)
    ax.set_yscale("log")
    ax.set_ylabel("stage wall-clock (s, log scale)")
    # custom legend: tool colours + segment meaning
    h = [plt.Rectangle((0, 0), 1, 1, color=TCOL[t], ec="black", lw=0.4) for t in TOOLS]
    seg = [plt.Rectangle((0, 0), 1, 1, fc="gray", ec="black", lw=0.4),
           plt.Rectangle((0, 0), 1, 1, fc="gray", ec="black", lw=0.4, alpha=0.5),
           plt.Rectangle((0, 0), 1, 1, fc="gray", ec="black", lw=0.4, hatch="///")]
    leg1 = fig.legend(h, [TLAB[t] for t in TOOLS], ncol=4, loc="lower center",
                      bbox_to_anchor=(0.5, -0.02), frameon=False)
    fig.legend(seg, ["predict (solid)", "annotate (pale)", "total, no stage split"],
               ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.09), frameon=False)
    fig.add_artist(leg1)
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    save(fig, "fig4_runtime_decomposition")


# ---------- Fig 5: cpu_seconds (log y) ----------
def fig5():
    tc = load("timing_cpu.tsv")
    fig, ax = plt.subplots(figsize=(DOUBLE, DOUBLE * 0.42))
    grouped_bars(ax, tc, "cpu_s_median", logy=True)
    ax.set_ylabel("CPU time (s, log scale)")
    tool_legend(fig, y=-0.08)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, "fig5_cpu_seconds")


# ---------- Fig 6: peak_memory (a: peak log; b: predict vs annotate RSS) ----------
def fig6():
    rss = load("rss.tsv")
    d = {(r["genome"], r["tool"]): r for r in rss}
    fig, axs = plt.subplots(1, 2, figsize=(DOUBLE, DOUBLE * 0.42))
    grouped_bars(axs[0], rss, "peak_rss_mb_median", logy=True)
    axs[0].set_ylabel("peak RSS (MB, log scale)")
    panel(axs[0], "a")
    # panel b: predict vs annotate RSS, grouped (only tools that have the split)
    ax = axs[1]
    split_tools = ["myconote", "funannotate"]
    series = [(t, stage) for t in split_tools for stage in ("predict", "annotate")]
    n = len(series); w = 0.8 / n; x = np.arange(len(GENOMES))
    for i, (t, stage) in enumerate(series):
        xs = x + (i - (n - 1) / 2) * w
        key = f"{stage}_rss_mb_median"
        vals = [fnum(d.get((g, t), {}).get(key)) or np.nan for g in GENOMES]
        ax.bar(xs, vals, w, color=TCOL[t], alpha=1.0 if stage == "predict" else 0.5,
               edgecolor="black", linewidth=0.3,
               label=f"{TLAB[t]} {stage}")
    ax.set_xticks(x); ax.set_xticklabels(GENOMES, rotation=90)
    ax.set_ylabel("stage peak RSS (MB)")
    ax.legend(ncol=2, frameon=False, fontsize=5.5, loc="upper left")
    panel(ax, "b")
    tool_legend(fig, y=-0.06)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "fig6_peak_memory")


# ---------- Fig 7: myconote_annotate_sources (stacked, matched) ----------
def fig7():
    aps = [r for r in load("annotate_per_source.tsv") if r["variant"] == "matched"]
    srcs = ["mmseqs", "pfam", "busco", "cazyme", "merops"]
    scol = {s: OKABE[i] for i, s in enumerate(srcs)}
    d = {(r["genome"], r["source"]): fnum(r["seconds_median"]) for r in aps}
    fig, ax = plt.subplots(figsize=(DOUBLE, DOUBLE * 0.42))
    x = np.arange(len(GENOMES))
    bottom = np.zeros(len(GENOMES))
    for s in srcs:
        vals = np.array([d.get((g, s), 0.0) or 0.0 for g in GENOMES])
        ax.bar(x, vals, 0.65, bottom=bottom, color=scol[s], edgecolor="black", linewidth=0.3, label=s)
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels(GENOMES, rotation=90)
    ax.set_ylabel("annotate time by source (s)")
    ax.legend(ncol=5, loc="upper center", frameon=False, fontsize=6)
    fig.tight_layout()
    save(fig, "fig7_myconote_annotate_sources")


# ---------- Fig 8: scaling (log-log wall vs size, per-tool fit) ----------
def fig8():
    sc = load("scaling.tsv")
    fig, ax = plt.subplots(figsize=(SINGLE * 1.6, SINGLE * 1.15))
    for t in TOOLS:
        pts = [(fnum(r["size_mb"]), fnum(r["wall_min_median"])) for r in sc
               if r["tool"] == t and fnum(r["wall_min_median"]) and fnum(r["wall_min_median"]) > 0]
        if not pts:
            continue
        xs, ys = zip(*pts)
        ax.scatter(xs, ys, s=18, color=TCOL[t], edgecolor="black", linewidth=0.3, label=TLAB[t], zorder=3)
        lx, ly = np.log10(xs), np.log10(ys)
        b, a = np.polyfit(lx, ly, 1)
        xf = np.array([min(xs), max(xs)])
        ax.plot(xf, 10 ** (a + b * np.log10(xf)), color=TCOL[t], lw=1.0, alpha=0.8, zorder=2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("genome size (Mb, log scale)")
    ax.set_ylabel("wall-clock time (min, log scale)")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    save(fig, "fig8_scaling")


# ---------- Fig 9: accuracy_speed_tradeoff (one marker per tool) ----------
def fig9():
    acc = load("accuracy.tsv")
    tw = load("timing_wall.tsv")
    fig, ax = plt.subplots(figsize=(SINGLE * 1.6, SINGLE * 1.15))
    for t in TOOLS:
        f1 = [fnum(r["gene_f1_loose"]) for r in acc if r["tool"] == t and fnum(r["gene_f1_loose"]) is not None]
        wm = [fnum(r["wall_s_median"]) / 60 for r in tw if r["tool"] == t and fnum(r["wall_s_median"])]
        if not f1 or not wm:
            continue
        ax.scatter(st.median(wm), st.mean(f1), s=70, color=TCOL[t], edgecolor="black",
                   linewidth=0.5, zorder=3)
        ax.annotate(TLAB[t], (st.median(wm), st.mean(f1)), xytext=(4, 4),
                    textcoords="offset points", fontsize=7)
    ax.set_xlabel("median wall-clock time across genomes (min)")
    ax.set_ylabel("mean gene F1 (loose) across genomes")
    ax.set_ylim(0.88, 0.94)
    ax.grid(True, lw=0.3, alpha=0.4)
    fig.tight_layout()
    save(fig, "fig9_accuracy_speed_tradeoff")


# ---------- Fig 10: accuracy_by_clade (a loose, b strict) ----------
def fig10():
    acc = load("accuracy.tsv")
    fig, axs = plt.subplots(1, 2, figsize=(DOUBLE, DOUBLE * 0.42))
    for ax, (k, lab), L in zip(axs, [("gene_f1_loose", "gene F1 (loose)"),
                                      ("gene_f1_strict", "gene F1 (strict)")], "ab"):
        n = len(TOOLS); w = 0.8 / n; x = np.arange(len(CLADES))
        for i, t in enumerate(TOOLS):
            vals = []
            for cl in CLADES:
                vs = [fnum(r[k]) for r in acc if r["tool"] == t and r["clade"] == cl and fnum(r[k]) is not None]
                vals.append(st.mean(vs) if vs else np.nan)
            ax.bar(x + (i - (n - 1) / 2) * w, vals, w, color=TCOL[t], edgecolor="black", linewidth=0.3)
        ax.set_xticks(x)
        ax.set_xticklabels(["Saccharo-\nmycotina", "other\nAscomycota", "Basidio-\nmycota"], fontsize=6)
        ax.set_ylabel(lab); ax.set_ylim(0, 1.0)
        panel(ax, L)
    tool_legend(fig, y=-0.04)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "fig10_accuracy_by_clade")


# ---------- Fig 11: gene_count_calibration (pred vs ref, y=x) ----------
def fig11():
    gc = load("gene_counts.tsv")
    fig, ax = plt.subplots(figsize=(SINGLE * 1.6, SINGLE * 1.4))
    allv = []
    for t in TOOLS:
        xs = [fnum(r["n_ref_genes"]) for r in gc if r["tool"] == t and fnum(r["n_pred_genes_mean"]) is not None]
        ys = [fnum(r["n_pred_genes_mean"]) for r in gc if r["tool"] == t and fnum(r["n_pred_genes_mean"]) is not None]
        ax.scatter(xs, ys, s=20, color=TCOL[t], edgecolor="black", linewidth=0.3, label=TLAB[t], zorder=3)
        allv += xs + ys
    lo, hi = min(allv) * 0.9, max(allv) * 1.05
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.7, zorder=1, label="y = x")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("reference gene count")
    ax.set_ylabel("predicted gene count (mean)")
    ax.legend(frameon=False, loc="upper left", fontsize=6)
    fig.tight_layout()
    save(fig, "fig11_gene_count_calibration")


# ---------- Fig 12: eggnog_cost (per-source matched vs full, log y) ----------
def fig12():
    aps = load("annotate_per_source.tsv")
    srcs = ["mmseqs", "pfam", "busco", "cazyme", "merops", "eggnog"]
    agg = {}
    for variant in ("matched", "full"):
        for s in srcs:
            vals = [fnum(r["seconds_median"]) for r in aps
                    if r["variant"] == variant and r["source"] == s and fnum(r["seconds_median"])]
            agg[(variant, s)] = st.median(vals) if vals else 0.0
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.7, SINGLE * 1.2))
    x = np.arange(len(srcs)); w = 0.38
    ax.bar(x - w / 2, [agg[("matched", s)] or np.nan for s in srcs], w,
           color="#0072B2", edgecolor="black", linewidth=0.3, label="matched (no eggNOG, headline)")
    ax.bar(x + w / 2, [agg[("full", s)] or np.nan for s in srcs], w,
           color="#56B4E9", edgecolor="black", linewidth=0.3, label="full (with eggNOG, ref)")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(srcs, rotation=20, ha="right")
    ax.set_ylabel("annotate time by source (s, log scale)\nmedian across genomes")
    ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    save(fig, "fig12_eggnog_cost")


for fn in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12):
    fn()
print("all figures done.")
