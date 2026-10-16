#!/usr/bin/env python3
"""Aggregate the n=5 Stage-1 benchmark into figure-backing TSVs.

Authoritative aggregator for the 4-tool fungal-annotation benchmark. Handles the
data-layout facts the stock aggregate_metrics.py cannot:
  - MycoNote HEADLINE comes from results_timing_stage1_matched/ (no-eggNOG, v0.7.16);
    MycoNote FULL (with-eggNOG, v0.7.15) is a separate variant used ONLY for the
    eggNOG-cost figure and flagged contention-confounded.
  - BRAKER's performance.json is all-zero; its wall/CPU/RSS come from
    condor_history (dumped to figures/_braker_history_146852.tsv).
Collapse: accuracy -> per-genome MEAN (deterministic); timing -> per-genome MEDIAN.
Run on the cluster (light JSON parsing; submit-node safe).
"""
import json
import re
import statistics as st
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent  # htcondor_benchmark/
TIMING = BASE / "results_timing_stage1"
MATCHED = BASE / "results_timing_stage1_matched"
OUT = BASE / "figures"
OUT.mkdir(exist_ok=True)

GENOMES = ["sce", "cal", "ylp", "ani", "ncr", "cne", "klac", "kpha",
           "fgra", "afum", "umay", "ccin", "spom", "tmel"]
SIZE_MB = {"sce": 12, "cal": 14, "ylp": 21, "ani": 30, "ncr": 40, "cne": 19,
           "klac": 11, "kpha": 9, "fgra": 36, "afum": 29, "umay": 20,
           "ccin": 36, "spom": 13, "tmel": 125}
CLADE = {g: "Saccharomycotina" for g in ("sce", "cal", "ylp", "klac", "kpha")}
CLADE.update({g: "Basidiomycota" for g in ("cne", "umay", "ccin")})
CLADE.update({g: "other-Ascomycota" for g in ("ani", "ncr", "fgra", "afum", "spom", "tmel")})

TOOLS = ["myconote", "funannotate", "maker", "braker"]  # headline order
# variant -> result dir (braker handled via history)
DIRS = {"myconote": MATCHED / "myconote",
        "funannotate": TIMING / "funannotate",
        "maker": TIMING / "maker",
        "braker": TIMING / "braker",
        "myconote_full": TIMING / "myconote"}


def jload(p):
    try:
        return json.loads(Path(p).read_text())
    except Exception:
        return None


def reps_of(d):
    if not d.is_dir():
        return []
    return sorted(int(x.name[3:]) for x in d.glob("rep*") if x.name[3:].isdigit())


def braker_history():
    """(genome,rep) -> dict(wall,cpu,rss_mb) from the dumped condor_history."""
    out = {}
    hp = OUT / "_braker_history_146852.tsv"
    for line in hp.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        wall, ucpu, scpu, memMB = float(f[2]), float(f[3]), float(f[4]), float(f[5])
        args = f[8].split()
        g, rep = args[0], int(args[1])
        out[(g, rep)] = {"wall": wall, "cpu": ucpu + scpu, "rss_mb": memMB}
    return out


BRK = braker_history()
BUSCO_RE = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)")


def busco_summary(txt_path):
    try:
        t = Path(txt_path).read_text()
    except Exception:
        return None
    m = BUSCO_RE.search(t)
    if not m:
        return None
    lin = re.search(r"lineage dataset is:\s*(\S+)", t)
    return {"C": float(m.group(1)), "S": float(m.group(2)), "D": float(m.group(3)),
            "F": float(m.group(4)), "M": float(m.group(5)), "n": int(m.group(6)),
            "lineage": lin.group(1) if lin else "NA"}


def load_run(variant, g, rep):
    """One raw run -> flat dict (or None)."""
    d = DIRS[variant] / g / f"rep{rep}"
    mj = jload(d / "metrics.json")
    if mj is None:
        return None
    gl = mj.get("gene_level", {})
    ex = mj.get("exon_level", {})
    nt = mj.get("nucleotide_level", {})
    row = {
        "variant": variant, "tool": variant.replace("_full", ""),
        "genome": g, "rep": rep,
        "size_mb": SIZE_MB[g], "clade": CLADE[g],
        "gene_f1_loose": gl.get("loose", {}).get("F1"),
        "gene_f1_strict": gl.get("strict", {}).get("F1"),
        "exon_f1_exact": ex.get("exact_match", {}).get("F1"),
        "nt_f1": nt.get("F1"),
        "n_pred_genes": mj.get("predicted_summary", {}).get("n_genes"),
        "n_ref_genes": mj.get("reference_summary", {}).get("n_genes"),
        "wall_s": None, "cpu_s": None, "peak_rss_mb": None,
        "predict_wall_s": None, "annotate_wall_s": None,
        "predict_rss_mb": None, "annotate_rss_mb": None,
    }
    # timing
    if variant == "braker":
        h = BRK.get((g, rep))
        if h:
            row.update(wall_s=h["wall"], cpu_s=h["cpu"], peak_rss_mb=h["rss_mb"])
    else:
        pj = jload(d / "performance.json")
        if pj:
            row["wall_s"] = pj.get("total_wall_seconds")
            row["cpu_s"] = pj.get("total_cpu_seconds")
            row["peak_rss_mb"] = pj.get("peak_rss_mb")
            stg = pj.get("stages", {})
            if "predict" in stg:
                row["predict_wall_s"] = stg["predict"].get("wall_seconds")
                row["predict_rss_mb"] = round(stg["predict"].get("max_rss_kb", 0) / 1024, 1)
            if "annotate" in stg:
                row["annotate_wall_s"] = stg["annotate"].get("wall_seconds")
                row["annotate_rss_mb"] = round(stg["annotate"].get("max_rss_kb", 0) / 1024, 1)
    # per-source annotate timings (myconote variants only)
    row["_sources"] = {}
    tp = d / "annotate_out" / "annotate_timings.tsv"
    if tp.exists():
        for ln in tp.read_text().splitlines()[1:]:
            parts = ln.split("\t")
            if len(parts) == 2:
                row["_sources"][parts[0]] = float(parts[1])
    return row


# ---- load everything ----
raw = []
for variant in list(DIRS):
    for g in GENOMES:
        for rep in reps_of(DIRS[variant] / g):
            r = load_run(variant, g, rep)
            if r:
                raw.append(r)
print(f"loaded {len(raw)} raw runs")


def col(rows, key):
    return [r[key] for r in rows if r.get(key) is not None]


def by(variant_or_tool, g, use_tool=False):
    key = "tool" if use_tool else "variant"
    return [r for r in raw if r[key] == variant_or_tool and r["genome"] == g]


def mean(xs):
    xs = [x for x in xs if x is not None]
    return round(st.mean(xs), 4) if xs else ""


def med(xs):
    xs = [x for x in xs if x is not None]
    return round(st.median(xs), 2) if xs else ""


def write_tsv(name, header, rows):
    p = OUT / name
    with open(p, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if v is None else str(v) for v in r) + "\n")
    print(f"wrote {p}  ({len(rows)} rows)")


# ---- per_run_raw.tsv (all runs) ----
RAWCOLS = ["variant", "tool", "genome", "clade", "size_mb", "rep",
           "gene_f1_loose", "gene_f1_strict", "exon_f1_exact", "nt_f1",
           "n_pred_genes", "n_ref_genes", "wall_s", "cpu_s", "peak_rss_mb",
           "predict_wall_s", "annotate_wall_s", "predict_rss_mb", "annotate_rss_mb"]
write_tsv("per_run_raw.tsv", RAWCOLS, [[r[c] for c in RAWCOLS] for r in raw])

# ---- accuracy.tsv (per-genome MEAN, headline tools) ----
acc = []
for g in GENOMES:
    for t in TOOLS:
        rs = by(t, g)
        if not rs:
            continue
        acc.append([g, CLADE[g], t,
                    mean(col(rs, "gene_f1_loose")), mean(col(rs, "gene_f1_strict")),
                    mean(col(rs, "exon_f1_exact")), mean(col(rs, "nt_f1"))])
write_tsv("accuracy.tsv",
          ["genome", "clade", "tool", "gene_f1_loose", "gene_f1_strict",
           "exon_f1_exact", "nt_f1"], acc)

# clade_strata.tsv = accuracy.tsv (clade column already present) -- same backing, named per spec
write_tsv("clade_strata.tsv",
          ["genome", "clade", "tool", "gene_f1_loose", "gene_f1_strict",
           "exon_f1_exact", "nt_f1"], acc)

# ---- timing_wall.tsv (per-genome MEDIAN wall + stage medians) ----
tw = []
for g in GENOMES:
    for t in TOOLS:
        rs = by(t, g)
        if not rs:
            continue
        tw.append([g, CLADE[g], SIZE_MB[g], t,
                   med(col(rs, "wall_s")), med(col(rs, "predict_wall_s")),
                   med(col(rs, "annotate_wall_s"))])
write_tsv("timing_wall.tsv",
          ["genome", "clade", "size_mb", "tool", "wall_s_median",
           "predict_wall_s_median", "annotate_wall_s_median"], tw)

# ---- timing_cpu.tsv ----
tc = []
for g in GENOMES:
    for t in TOOLS:
        rs = by(t, g)
        if rs:
            tc.append([g, t, med(col(rs, "cpu_s"))])
write_tsv("timing_cpu.tsv", ["genome", "tool", "cpu_s_median"], tc)

# ---- rss.tsv ----
rr = []
for g in GENOMES:
    for t in TOOLS:
        rs = by(t, g)
        if rs:
            rr.append([g, t, med(col(rs, "peak_rss_mb")),
                       med(col(rs, "predict_rss_mb")), med(col(rs, "annotate_rss_mb"))])
write_tsv("rss.tsv",
          ["genome", "tool", "peak_rss_mb_median", "predict_rss_mb_median",
           "annotate_rss_mb_median"], rr)

# ---- scaling.tsv ----
sc = []
for g in GENOMES:
    for t in TOOLS:
        rs = by(t, g)
        if rs:
            w = med(col(rs, "wall_s"))
            sc.append([g, SIZE_MB[g], CLADE[g], t,
                       round(w / 60.0, 3) if w != "" else ""])
write_tsv("scaling.tsv",
          ["genome", "size_mb", "clade", "tool", "wall_min_median"], sc)

# ---- gene_counts.tsv ----
gc = []
for g in GENOMES:
    ref = ""
    for t in TOOLS:
        rs = by(t, g)
        if rs:
            if ref == "":
                ref = round(mean(col(rs, "n_ref_genes"))) if col(rs, "n_ref_genes") else ""
            gc.append([g, t, round(mean(col(rs, "n_pred_genes"))) if col(rs, "n_pred_genes") else "", ref])
write_tsv("gene_counts.tsv", ["genome", "tool", "n_pred_genes_mean", "n_ref_genes"], gc)

# ---- annotate_per_source.tsv (myconote matched + full), per-genome MEDIAN ----
aps = []
ALL_SRC = ["mmseqs", "pfam", "busco", "cazyme", "merops", "eggnog"]
for variant, label in (("myconote", "matched"), ("myconote_full", "full")):
    for g in GENOMES:
        rs = [r for r in raw if r["variant"] == variant and r["genome"] == g]
        if not rs:
            continue
        for s in ALL_SRC:
            vals = [r["_sources"].get(s) for r in rs if s in r["_sources"]]
            if vals:
                aps.append([g, label, s, round(st.median(vals), 3)])
write_tsv("annotate_per_source.tsv",
          ["genome", "variant", "source", "seconds_median"], aps)

# ---- completeness.tsv (BUSCO-of-proteins; ONLY myconote + funannotate) ----
comp = []
for g in GENOMES:
    mp = MATCHED / "myconote" / g / "rep1" / "annotate_out" / "busco"
    cand = list(mp.glob("*/short_summary.specific.*.txt")) if mp.is_dir() else []
    b = busco_summary(cand[0]) if cand else None
    if b:
        comp.append([g, "myconote", b["lineage"], b["C"], b["S"], b["D"], b["F"], b["M"], b["n"]])
    fp = TIMING / "funannotate" / g / "rep1" / "predict_out" / "annotate_misc" / "run_busco" / "short_summary_busco.txt"
    b = busco_summary(fp)
    if b:
        comp.append([g, "funannotate", b["lineage"], b["C"], b["S"], b["D"], b["F"], b["M"], b["n"]])
if comp:
    write_tsv("completeness.tsv",
              ["genome", "tool", "lineage", "C", "S", "D", "F", "M", "n"], comp)
    print("NOTE: completeness.tsv covers myconote+funannotate ONLY; different "
          "lineages/BUSCO versions -> NOT a fair cross-tool axis.")
else:
    print("completeness.tsv: no parseable BUSCO summaries -> skipped")

print("aggregation done.")
