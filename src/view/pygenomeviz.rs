//! pyGenomeViz backend for synteny rendering.
//!
//! Produces publication-quality static figures (PNG / PDF / SVG at 300 DPI)
//! by shelling out to the pyGenomeViz Python library. The Rust side collects
//! the per-genome inputs + plot parameters into a JSON spec and invokes a
//! generated Python driver that drives pyGenomeViz's built-in
//! MUMmer / BLAST / MMseqs aligners.
//!
//! This is checkpoint 1 of the pygenomeviz integration — focus is a clean
//! default output on 2-genome inputs. Gene-name cascades, KEGG coloring,
//! region zoom, and themes are wired through the JSON spec so later
//! checkpoints only touch the Python template.
//!
//! Python runtime: `myconote-pygv-python` wrapper at
//! `/Users/black_einstein/miniconda3/bin/myconote-pygv-python` points at a
//! dedicated conda env (`myconote-viz`) that ships pygenomeviz + MUMmer +
//! BLAST. Users install with:
//!     conda create -n myconote-viz -c bioconda -c conda-forge \
//!         python=3.12 pygenomeviz mummer4 blast matplotlib -y

use crate::utils::error::{MycoNoteError, Result};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use super::synteny::{PlotFormat, PlotTheme, PygvAligner, SyntenyConfig};

/// Render a synteny figure via pyGenomeViz. Assumes config.backend is
/// `PyGenomeViz`.
pub fn render_synteny(config: &SyntenyConfig) -> Result<()> {
    // ── Pre-flight checks ─────────────────────────────────────────────────
    if config.fasta1.is_none() || config.fasta2.is_none() {
        return Err(MycoNoteError::InvalidFormat(
            "pyGenomeViz backend requires --fasta1 and --fasta2 (the Python \
             library drives its own aligner on nucleotide FASTAs, not on pre-\
             computed PAF). Add --fasta1 <A.fa> --fasta2 <B.fa> and retry."
                .to_string(),
        ));
    }

    let python = locate_python()?;

    // Ensure output parent dir exists — users routinely pass nested paths.
    if let Some(parent) = config.output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(MycoNoteError::Io)?;
        }
    }

    // ── Build the JSON spec ──────────────────────────────────────────────
    let spec = build_spec_json(config);

    // Stage spec + driver script in a nanosecond-suffixed temp dir so
    // concurrent runs don't collide.
    let suffix = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp_dir = std::env::temp_dir().join(format!("myconote_pygv_{}", suffix));
    std::fs::create_dir_all(&tmp_dir).map_err(MycoNoteError::Io)?;

    let spec_path = tmp_dir.join("spec.json");
    std::fs::write(&spec_path, spec).map_err(MycoNoteError::Io)?;

    let driver_path = tmp_dir.join("driver.py");
    std::fs::write(&driver_path, PYTHON_DRIVER).map_err(MycoNoteError::Io)?;

    // ── Run Python driver ────────────────────────────────────────────────
    println!(
        "🎨 Rendering synteny figure via pyGenomeViz ({} aligner, {} theme)…",
        config.aligner.as_py_class(),
        config.theme.as_py_str()
    );
    let out = Command::new(&python)
        .arg(&driver_path)
        .arg(&spec_path)
        .output()
        .map_err(MycoNoteError::Io)?;

    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        let stdout = String::from_utf8_lossy(&out.stdout);
        let _ = std::fs::remove_dir_all(&tmp_dir);
        return Err(MycoNoteError::ExternalTool(format!(
            "pyGenomeViz driver failed (exit {}):\n--- stderr ---\n{}\n--- stdout ---\n{}",
            out.status.code().unwrap_or(-1),
            stderr.trim(),
            stdout.trim()
        )));
    }

    // Surface Python's stdout (useful diagnostics like "8 links plotted").
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        if !line.trim().is_empty() {
            println!("  {}", line);
        }
    }

    let _ = std::fs::remove_dir_all(&tmp_dir);

    println!("✓ Synteny figure: {}", config.output.display());
    Ok(())
}

/// Locate the pygenomeviz-capable Python interpreter. We install a wrapper
/// script in miniconda's base bin dir (see module docstring) to skirt the
/// shebang-resolution issue on fresh shells; if the wrapper is absent, fall
/// back to `python3` on PATH and hope pygenomeviz is importable.
fn locate_python() -> Result<PathBuf> {
    if let Ok(p) = which::which("myconote-pygv-python") {
        return Ok(p);
    }
    which::which("python3").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "myconote-pygv-python wrapper not found and `python3` not on PATH.\n\
             Install the viz env with:\n  \
             conda create -n myconote-viz -c bioconda -c conda-forge \\\n    \
                 python=3.12 pygenomeviz mummer4 blast matplotlib -y\n  \
             Then add a wrapper at ~/miniconda3/bin/myconote-pygv-python that\n  \
             exports the env's bin path and exec's its python."
                .to_string(),
        )
    })
}

/// Emit the JSON spec consumed by the Python driver. Hand-rolled to avoid
/// pulling in serde_json here — the structure is stable and small.
fn build_spec_json(config: &SyntenyConfig) -> String {
    let fa1 = config.fasta1.as_ref().unwrap().display().to_string();
    let fa2 = config.fasta2.as_ref().unwrap().display().to_string();
    let gff1 = config.gff1.display().to_string();
    let gff2 = config.gff2.display().to_string();
    let out = config.output.display().to_string();

    let region_json = match &config.region {
        Some(r) => escape_json(r),
        None => "null".to_string(),
    };
    // Implicit rule: --region turns on gene-feature arrows even when the
    // user didn't pass --show-features; the whole point of zooming is to
    // inspect features. show_features can still be explicitly set true for
    // whole-genome views if the user really wants it.
    let show_features = config.show_features || config.region.is_some();

    format!(
        r#"{{
  "genomes": [
    {{"name": {label1}, "fasta": {fa1}, "gff": {gff1}}},
    {{"name": {label2}, "fasta": {fa2}, "gff": {gff2}}}
  ],
  "output": {out},
  "format": {fmt},
  "aligner": {aligner},
  "theme": {theme},
  "min_identity": {min_id},
  "min_block_len": {min_block},
  "width_inches": {width},
  "track_height_inches": {track_h},
  "feature_track_ratio": 0.25,
  "dpi": 300,
  "show_features": {show_features},
  "show_gene_labels": {show_labels},
  "region": {region},
  "top_contigs": {top_contigs}
}}"#,
        label1 = escape_json(&config.label1),
        label2 = escape_json(&config.label2),
        fa1 = escape_json(&fa1),
        fa2 = escape_json(&fa2),
        gff1 = escape_json(&gff1),
        gff2 = escape_json(&gff2),
        out = escape_json(&out),
        fmt = escape_json(config.format.extension()),
        aligner = escape_json(config.aligner.as_py_class()),
        theme = escape_json(config.theme.as_py_str()),
        min_id = config.min_identity,
        min_block = config.min_block_len,
        width = config.width_inches,
        track_h = config.track_height_inches,
        show_features = show_features,
        show_labels = config.region.is_some(), // only label in region mode
        region = region_json,
        top_contigs = config.top_contigs,
    )
}

/// Escape a string so it can be dropped between `"…"` in JSON.
fn escape_json(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Python driver embedded as a static string. Reads the JSON spec from
/// argv[1], wires it to pyGenomeViz, saves the figure. Kept intentionally
/// small — feature-level customisation lands in later checkpoints.
const PYTHON_DRIVER: &str = r#"#!/usr/bin/env python3
"""
myconote pyGenomeViz driver — checkpoint 1.

Consumes a JSON spec from argv[1] produced by the Rust side and renders a
publication-quality synteny figure. Exits non-zero with a clear message on
missing deps or alignment failures; Rust surfaces both streams.
"""
import json
import os
import sys
from pathlib import Path


def fail(msg, code=1):
    print(msg, flush=True)
    sys.exit(code)


if len(sys.argv) < 2:
    fail("usage: driver.py <spec.json>", 2)

try:
    spec = json.loads(Path(sys.argv[1]).read_text())
except Exception as e:  # noqa: BLE001
    fail(f"failed to parse spec: {e}", 2)

try:
    from pygenomeviz import GenomeViz
    from pygenomeviz.align import Blast, MMseqs, MUMmer, AlignCoord
    from pygenomeviz.parser import Fasta, Gff
    from pygenomeviz.utils import ColorCycler
    import matplotlib

    matplotlib.use("Agg")  # headless — no Tk / Qt needed
    import matplotlib.pyplot as plt
except ImportError as e:  # noqa: BLE001
    fail(
        "pyGenomeViz import failed. Install with:\n  "
        "conda install -c bioconda pygenomeviz\n"
        f"  underlying error: {e}",
        3,
    )

aligner_name = spec["aligner"]
aligner_cls = {"MUMmer": MUMmer, "Blast": Blast, "MMseqs": MMseqs}.get(aligner_name)
if aligner_cls is None:
    fail(f"unknown aligner {aligner_name!r}", 2)

theme = spec.get("theme", "light")
width = float(spec.get("width_inches", 12.0))
track_h = float(spec.get("track_height_inches", 1.2))
feat_ratio = float(spec.get("feature_track_ratio", 0.25))
min_identity = float(spec.get("min_identity", 30))
min_block = int(spec.get("min_block_len", 1000))
show_features = bool(spec.get("show_features", False))
show_gene_labels = bool(spec.get("show_gene_labels", False))
region = spec.get("region")
top_contigs = int(spec.get("top_contigs", 0))
dpi = int(spec.get("dpi", 300))

ColorCycler.set_cmap("tab10")


def parse_region(s):
    """seqid:start-end → (seqid, start, end) or fail."""
    if ":" not in s:
        fail(f"invalid --region {s!r}; expected seqid:start-end", 2)
    seqid, coords = s.split(":", 1)
    if "-" not in coords:
        fail(f"invalid --region {s!r}; expected seqid:start-end", 2)
    lo, hi = coords.replace(",", "").split("-", 1)
    try:
        return seqid.strip(), int(lo), int(hi)
    except ValueError:
        fail(f"invalid --region {s!r}; coords must be integers", 2)


region_seqid, region_start, region_end = (None, None, None)
if region:
    region_seqid, region_start, region_end = parse_region(region)
    print(
        f"Region view: {region_seqid}:{region_start:,}-{region_end:,}",
        flush=True,
    )

# Parse the FASTA inputs so we can size chromosome tracks.
fasta_list = []
gff_list = []
for g in spec["genomes"]:
    fa_path = g["fasta"]
    gff_path = g.get("gff", "")
    if not Path(fa_path).exists():
        fail(f"FASTA not found: {fa_path}", 2)
    fasta_list.append(Fasta(fa_path))
    gff_list.append(Gff(gff_path) if gff_path and Path(gff_path).exists() else None)

# Optionally filter to top-N largest contigs per genome for clarity.
def filter_seqid2size(seqid2size, n):
    if n <= 0 or n >= len(seqid2size):
        return seqid2size
    items = sorted(seqid2size.items(), key=lambda kv: kv[1], reverse=True)[:n]
    return dict(items)


# ── Build GenomeViz figure ───────────────────────────────────────────────
gv = GenomeViz(
    fig_track_height=track_h,
    feature_track_ratio=feat_ratio,
    theme=theme,
    track_align_type="center",
)
# In region mode the scale is tight, use xticks; in whole-genome mode a
# single scale bar in the corner is less cluttered.
if region:
    gv.set_scale_xticks(ymargin=0.5)
else:
    gv.set_scale_bar(ymargin=0.5)

# Disambiguate parser names so AlignCoord refs don't collide on self-compare.
seen_stems = {}
for idx, fa in enumerate(fasta_list):
    base = fa.name
    seen_stems[base] = seen_stems.get(base, 0) + 1
    if seen_stems[base] > 1:
        fa.name = f"{base}#{idx + 1}"

track_colors = []
for g, fa, gff in zip(spec["genomes"], fasta_list, gff_list):
    color = ColorCycler()
    track_colors.append(color)

    # Decide the track's segment layout: either a single region range, or
    # (optionally top-N) contigs at full length.
    if region:
        segs = (region_start, region_end)
    else:
        segs = filter_seqid2size(fa.get_seqid2size(), top_contigs)

    track = gv.add_feature_track(
        fa.name,
        segs,
        label_kws={"color": color if theme == "dark" else "black"},
        align_label=False,
    )
    track.set_label(g["name"])

    for segment in track.segments:
        segment.add_feature(
            segment.start,
            segment.end,
            plotstyle="bigrbox",
            fc=color,
            ec="black",
            lw=0.4,
        )
        # Sublabels are useful in region mode and when few contigs; with
        # 20+ contigs the "0 - X bp" strings overlap into noise.
        seg_count = len(track.segments)
        if region or seg_count <= 6:
            segment.add_sublabel(ymargin=0.3)

    # Gene feature arrows only when explicitly requested — at
    # whole-chromosome × 6000-gene scale they compress to solid smears.
    if show_features and gff is not None:
        seqid_filter = {region_seqid} if region else None
        for seqid, features in gff.get_seqid2features("CDS").items():
            if seqid_filter is not None and seqid not in seqid_filter:
                continue
            try:
                segment = track.get_segment(seqid)
            except Exception:  # noqa: BLE001
                continue
            if region:
                features = [
                    f
                    for f in features
                    if not (
                        int(f.location.end) < region_start
                        or int(f.location.start) > region_end
                    )
                ]
            segment.add_features(
                features,
                plotstyle="bigarrow",
                fc="skyblue" if theme == "light" else "steelblue",
                lw=0.3,
                label_type="gene" if show_gene_labels else None,
            )

# ── Run aligner ──────────────────────────────────────────────────────────
print(f"Running {aligner_name} alignment on {len(fasta_list)} genomes…",
      flush=True)
aligner = aligner_cls(fasta_list)
align_coords = aligner.run()
align_coords = AlignCoord.filter(
    align_coords, length_thr=min_block, identity_thr=min_identity
)
print(f"{len(align_coords)} alignment link(s) passed filters "
      f"(min_len={min_block}, min_id={min_identity}).", flush=True)

# ── Plot alignment links ─────────────────────────────────────────────────
# --top-contigs and --region drop some contigs from the tracks; MUMmer
# may still return alignments touching those dropped contigs, and
# gv.add_link raises if it can't find the segment. Try each link and
# swallow the not-found exception rather than pre-filtering (which is
# hard to do 100%-correctly across pygenomeviz's internal name mangling).
from pygenomeviz.exception import (  # type: ignore
    FeatureTrackNotFoundError,
    SegmentNotFoundError,
)

if align_coords:
    idents = [ac.identity for ac in align_coords if ac.identity is not None]
    min_ident_seen = int(min(idents)) if idents else int(min_identity)
    # Degenerate case: all links share the same identity (e.g. self-
    # compare → 100 everywhere). With vmin=vmax=100 the colormap collapses
    # to the lightest shade and links become invisible. Drop the floor so
    # full-identity renders at the saturated end of the gradient.
    if min_ident_seen >= 99:
        min_ident_seen = max(int(min_identity) - 5, 50)
    # Saturated colors so links are visible against a light-grey link track.
    # Grey-on-grey was invisible; light-blue/red pops and still reads as
    # "forward = same orientation, red = inverted" to most readers.
    if theme == "dark":
        fwd_color, inv_color = "skyblue", "tomato"
    else:
        fwd_color, inv_color = "steelblue", "tomato"
    plotted = 0
    skipped = 0
    for ac in align_coords:
        # Region view: clip by coord range on the query side before trying
        # to plot (avoids draw-off-screen links + saves try/except cost).
        if region:
            q_seqid = ac.query_link[1]
            q_start = ac.query_link[2]
            q_end = ac.query_link[3]
            if q_seqid != region_seqid:
                skipped += 1
                continue
            if q_end < region_start or q_start > region_end:
                skipped += 1
                continue
        try:
            # Dropped `size=0.95` — that compresses link polygons vertically
            # to 95% of the track gap, which at small-track heights becomes
            # effectively invisible. Default (1.0, full track gap) is what
            # the pygenomeviz docs use.
            # Dropped `curve=False` — straight links are the default and the
            # flag wasn't behaving as expected on self-compare data (links
            # rendered but not visible against the track-gap background).
            gv.add_link(
                ac.query_link,
                ac.ref_link,
                color=fwd_color,
                inverted_color=inv_color,
                v=ac.identity,
                vmin=min_ident_seen,
            )
            plotted += 1
        except (FeatureTrackNotFoundError, SegmentNotFoundError):
            # Link references a contig that was filtered out
            # (--top-contigs). Silently skip.
            skipped += 1

    print(f"Plotted {plotted} link(s); {skipped} skipped (filtered contigs / region).",
          flush=True)
    if plotted > 0:
        gv.set_colorbar(
            [fwd_color, inv_color], vmin=min_ident_seen, bar_label="Identity (%)"
        )
else:
    print("⚠  No alignment links passed filters — figure shows genome tracks only.",
          flush=True)

# ── Size + save ──────────────────────────────────────────────────────────
n_taxa = len(spec["genomes"])
fig = gv.plotfig()
fig.set_size_inches(width, max(3.0, track_h * n_taxa + 1.5))

out_path = Path(spec["output"])
out_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out_path} ({dpi} DPI, {spec['format']}).")
"#;
