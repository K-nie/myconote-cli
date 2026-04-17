/// Synteny visualisation — 2-genome ribbon diagram
///
/// Workflow:
///   1. Run minimap2 in asm-to-asm mode to align two genomes
///   2. Parse the PAF output natively in Rust → syntenic blocks
///   3. Render a self-contained HTML file with a D3.js ribbon diagram
///      - Top panel: genome A chromosomes (to scale)
///      - Bottom panel: genome B chromosomes (to scale)
///      - Ribbons: one coloured band per syntenic block, shaded by identity
///      - Gene labels: shown on blocks when a names map is provided
///
/// The output is a single .html file with no external dependencies
/// (D3.js loaded from CDN, graceful offline error message).
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

pub struct SyntenyConfig {
    pub gff1: PathBuf,
    pub gff2: PathBuf,
    pub fasta1: Option<PathBuf>,
    pub fasta2: Option<PathBuf>,
    pub output: PathBuf,
    pub label1: String,
    pub label2: String,
    /// Minimum alignment block length to display (filters noise)
    pub min_block_len: u64,
    /// Gene ID → display name for labelling blocks
    pub names: HashMap<String, String>,
    /// NCBI taxon ID for online name resolution (optional)
    pub taxon_id: Option<u32>,
    /// Maximum gap (bp, on both query and target sides) between adjacent PAF
    /// hits to merge into a single chained synteny block. Set to 0 to disable
    /// chaining. Default 100_000 — tuned for fungal genomes where asm5
    /// typically fragments single colinear regions into 5–50 sub-hits.
    pub chain_gap: u64,
    /// Thread count passed to minimap2 (`-t`). Default 4.
    pub threads: u32,
    /// Keep the intermediate PAF instead of deleting the temp directory
    /// after rendering. Useful for debugging alignment issues or loading the
    /// raw hits into an external dotplot tool.
    pub keep_paf: bool,
}

impl Default for SyntenyConfig {
    fn default() -> Self {
        Self {
            gff1: PathBuf::new(),
            gff2: PathBuf::new(),
            fasta1: None,
            fasta2: None,
            output: PathBuf::from("synteny.html"),
            label1: "Genome A".to_string(),
            label2: "Genome B".to_string(),
            min_block_len: 1_000,
            names: HashMap::new(),
            taxon_id: None,
            chain_gap: 100_000,
            threads: 4,
            keep_paf: false,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PAF block representation
// ─────────────────────────────────────────────────────────────────────────────

/// One syntenic block parsed from a PAF line (possibly post-chaining).
#[derive(Debug, Clone)]
pub struct SyntenyBlock {
    pub query_name: String,
    pub query_start: u64,
    pub query_end: u64,
    pub query_len: u64,
    pub target_name: String,
    pub target_start: u64,
    pub target_end: u64,
    pub target_len: u64,
    pub strand: char,  // '+' or '-'
    pub identity: f64, // 0.0–1.0, weighted by block_len when chained
    pub block_len: u64,
    /// Column 10 of PAF — number of matching residues. Preserved so chaining
    /// can recompute identity as a length-weighted mean.
    pub residue_matches: u64,
}

/// Chromosome size map keyed by seqid.
type ChromSizes = HashMap<String, u64>;

/// A named gene interval used to label ribbons: (start, end, gene_id).
#[derive(Debug, Clone)]
struct GeneInterval {
    start: u64,
    end: u64,
    id: String,
}

/// Per-contig gene intervals, keyed by seqid.
type GenesByContig = HashMap<String, Vec<GeneInterval>>;

// ─────────────────────────────────────────────────────────────────────────────
// PAF parsing
// ─────────────────────────────────────────────────────────────────────────────

/// Parse a PAF file and return syntenic blocks above `min_len`.
fn parse_paf(paf_path: &Path, min_len: u64) -> Result<Vec<SyntenyBlock>> {
    use std::io::BufRead;
    let file = std::fs::File::open(paf_path).map_err(MycoNoteError::Io)?;
    let reader = std::io::BufReader::new(file);
    let mut blocks = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }

        let query_len: u64 = f[1].parse().unwrap_or(0);
        let query_start: u64 = f[2].parse().unwrap_or(0);
        let query_end: u64 = f[3].parse().unwrap_or(0);
        let strand_char = f[4].chars().next().unwrap_or('+');
        let target_len: u64 = f[6].parse().unwrap_or(0);
        let target_start: u64 = f[7].parse().unwrap_or(0);
        let target_end: u64 = f[8].parse().unwrap_or(0);
        let residue_matches: u64 = f[9].parse().unwrap_or(0);
        let block_len: u64 = f[10].parse().unwrap_or(0);

        if block_len < min_len {
            continue;
        }

        let identity = if block_len > 0 {
            residue_matches as f64 / block_len as f64
        } else {
            0.0
        };

        blocks.push(SyntenyBlock {
            query_name: f[0].to_string(),
            query_start,
            query_end,
            query_len,
            target_name: f[5].to_string(),
            target_start,
            target_end,
            target_len,
            strand: strand_char,
            identity,
            block_len,
            residue_matches,
        });
    }

    Ok(blocks)
}

// ─────────────────────────────────────────────────────────────────────────────
// Block chaining
// ─────────────────────────────────────────────────────────────────────────────

/// Merge adjacent colinear PAF hits into consolidated synteny blocks.
///
/// minimap2 `asm5` fragments a single colinear region into many short hits
/// when small indels or repeats break alignment. This pass groups blocks by
/// `(query_contig, target_contig, strand)`, sorts by query_start, then walks
/// the sorted list merging neighbours when both sides stay within `max_gap`
/// and coordinates remain monotonic (strictly increasing on query; on target
/// increasing for '+' strand, decreasing for '-' strand).
///
/// Merged identity is a length-weighted mean via summed `residue_matches /
/// block_len`. Setting `max_gap == 0` disables chaining and returns the input
/// unchanged.
fn chain_blocks(mut blocks: Vec<SyntenyBlock>, max_gap: u64) -> Vec<SyntenyBlock> {
    if max_gap == 0 || blocks.len() < 2 {
        return blocks;
    }

    // Group key: (query_contig, target_contig, strand).
    blocks.sort_by(|a, b| {
        a.query_name
            .cmp(&b.query_name)
            .then_with(|| a.target_name.cmp(&b.target_name))
            .then_with(|| a.strand.cmp(&b.strand))
            .then_with(|| a.query_start.cmp(&b.query_start))
    });

    let mut out: Vec<SyntenyBlock> = Vec::with_capacity(blocks.len());
    for b in blocks {
        let merged = out.last_mut().and_then(|cur| {
            if cur.query_name != b.query_name
                || cur.target_name != b.target_name
                || cur.strand != b.strand
            {
                return None;
            }

            // Query gap — blocks are sorted by query_start, so b.qs >= cur.qs.
            // Use saturating_sub so overlapping blocks (b.qs < cur.qe) map to
            // a zero gap rather than underflowing.
            let q_gap = b.query_start.saturating_sub(cur.query_end);
            if q_gap > max_gap {
                return None;
            }

            // Target-side monotonicity and gap depend on strand.
            let t_gap = match cur.strand {
                '+' => {
                    if b.target_start < cur.target_start {
                        return None; // order inversion — not colinear
                    }
                    b.target_start.saturating_sub(cur.target_end)
                }
                '-' => {
                    // Reverse-strand chains: target coords run backwards as
                    // query advances, so cur.target_start must sit AFTER
                    // b.target_end.
                    if b.target_end > cur.target_start {
                        return None;
                    }
                    cur.target_start.saturating_sub(b.target_end)
                }
                _ => return None,
            };
            if t_gap > max_gap {
                return None;
            }

            // Extend bounding box on both axes.
            cur.query_start = cur.query_start.min(b.query_start);
            cur.query_end = cur.query_end.max(b.query_end);
            cur.target_start = cur.target_start.min(b.target_start);
            cur.target_end = cur.target_end.max(b.target_end);
            cur.residue_matches += b.residue_matches;
            cur.block_len += b.block_len;
            cur.identity = if cur.block_len > 0 {
                cur.residue_matches as f64 / cur.block_len as f64
            } else {
                0.0
            };
            Some(())
        });

        if merged.is_none() {
            out.push(b);
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// minimap2 runner
// ─────────────────────────────────────────────────────────────────────────────

/// Run minimap2 in asm-to-asm mode and return path to the PAF output file.
fn run_minimap2(fasta1: &Path, fasta2: &Path, out_dir: &Path, threads: u32) -> Result<PathBuf> {
    // Check minimap2 is available
    let mm2 = which::which("minimap2").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "minimap2 not found in PATH. Install it with: conda install -c bioconda minimap2\n\
             or: brew install minimap2"
                .to_string(),
        )
    })?;

    let paf_path = out_dir.join("synteny_alignment.paf");
    let threads = threads.max(1);

    println!("  Running minimap2 asm-to-asm alignment (-t {})…", threads);
    let status = Command::new(&mm2)
        .args([
            "-cx",
            "asm5", // asm-to-asm preset (≥5% divergence)
            "--cs", // include cs tag for identity calculation
            "-t",
            &threads.to_string(),
            fasta1.to_str().unwrap_or(""),
            fasta2.to_str().unwrap_or(""),
        ])
        .stdout(std::fs::File::create(&paf_path).map_err(MycoNoteError::Io)?)
        .stderr(std::process::Stdio::inherit())
        .status()
        .map_err(|e| MycoNoteError::Io(e))?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "minimap2 alignment failed. Check the FASTA files are valid.".to_string(),
        ));
    }

    println!("  Alignment complete → {}", paf_path.display());
    Ok(paf_path)
}

// ─────────────────────────────────────────────────────────────────────────────
// Chromosome sizes from GFF3
// ─────────────────────────────────────────────────────────────────────────────

fn chrom_sizes_from_gff(gff_path: &Path) -> Result<ChromSizes> {
    use crate::parser::gff::GFFReader;
    let mut sizes: ChromSizes = HashMap::new();
    for result in GFFReader::from_path(gff_path)? {
        if let Ok(rec) = result {
            let max = sizes.entry(rec.seqid).or_insert(0);
            if rec.end > *max {
                *max = rec.end;
            }
        }
    }
    Ok(sizes)
}

/// Extract all gene intervals from a GFF3, sorted by start within each
/// contig. Used both for the on-track gene overlay (rendered as tick marks)
/// and for the PAF-block → overlapping-gene label lookup; label text is then
/// resolved through the `geneNames` map on the JS side, falling back to the
/// gene ID when no name is available.
fn genes_from_gff(gff_path: &Path) -> Result<GenesByContig> {
    use crate::parser::gff::GFFReader;
    let mut genes: GenesByContig = HashMap::new();
    for result in GFFReader::from_path(gff_path)? {
        let rec = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        if rec.feature_type != "gene" {
            continue;
        }
        let Some(id) = rec.id().cloned() else {
            continue;
        };
        genes.entry(rec.seqid).or_default().push(GeneInterval {
            start: rec.start,
            end: rec.end,
            id,
        });
    }
    for v in genes.values_mut() {
        v.sort_by_key(|g| g.start);
    }
    Ok(genes)
}

// ─────────────────────────────────────────────────────────────────────────────
// HTML renderer
// ─────────────────────────────────────────────────────────────────────────────

/// Serialise chromosome sizes to a JS object literal: `{chr: length, ...}`
fn chrom_sizes_to_js(sizes: &ChromSizes, label: &str) -> String {
    let mut sorted: Vec<(&String, &u64)> = sizes.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(b.0));
    let inner: String = sorted
        .iter()
        .map(|(k, v)| format!("\"{}\":{}", k, v))
        .collect::<Vec<_>>()
        .join(",");
    format!("const chrSizes_{} = {{{}}};", label, inner)
}

/// Serialise blocks to a JS array of objects.
fn blocks_to_js(blocks: &[SyntenyBlock]) -> String {
    let items: Vec<String> =
        blocks
            .iter()
            .map(|b| {
                format!(
            "{{qn:\"{}\",qs:{},qe:{},ql:{},tn:\"{}\",ts:{},te:{},tl:{},st:\"{}\",id:{:.4}}}",
            b.query_name, b.query_start, b.query_end, b.query_len,
            b.target_name, b.target_start, b.target_end, b.target_len,
            b.strand, b.identity
        )
            })
            .collect();
    format!("const synBlocks = [{}];", items.join(","))
}

/// Serialise the names map to JS: `{id: name, ...}`
fn names_to_js(names: &HashMap<String, String>) -> String {
    let inner: String = names
        .iter()
        .map(|(k, v)| {
            let v_esc = v.replace('"', "\\\"");
            format!("\"{}\":\"{}\"", k, v_esc)
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("const geneNames = {{{}}};", inner)
}

/// Serialise per-contig gene intervals to JS:
/// `{contig: [[start, end, "id"], ...], ...}` — sorted by start.
fn genes_to_js(genes: &GenesByContig, suffix: &str) -> String {
    let mut entries: Vec<(&String, &Vec<GeneInterval>)> = genes.iter().collect();
    entries.sort_by(|a, b| a.0.cmp(b.0));
    let body: String = entries
        .iter()
        .map(|(contig, ivs)| {
            let arr: String = ivs
                .iter()
                .map(|g| format!("[{},{},\"{}\"]", g.start, g.end, g.id))
                .collect::<Vec<_>>()
                .join(",");
            format!("\"{}\":[{}]", contig, arr)
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("const genes_{} = {{{}}};", suffix, body)
}

fn render_synteny_html(
    blocks: &[SyntenyBlock],
    sizes1: &ChromSizes,
    sizes2: &ChromSizes,
    genes1: &GenesByContig,
    genes2: &GenesByContig,
    config: &SyntenyConfig,
) -> String {
    let label1 = &config.label1;
    let label2 = &config.label2;

    let js_sizes1 = chrom_sizes_to_js(sizes1, "A");
    let js_sizes2 = chrom_sizes_to_js(sizes2, "B");
    let js_blocks = blocks_to_js(blocks);
    let js_names = names_to_js(&config.names);
    let js_genes1 = genes_to_js(genes1, "A");
    let js_genes2 = genes_to_js(genes2, "B");

    let block_count = blocks.len();
    let title = format!("MycoNote — Synteny: {} vs {}", label1, label2);

    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>{title}</title>
<style>
  :root{{--bg:#0f1117;--panel:#1a1d27;--accent:#00d4aa;--text:#e2e8f0;--muted:#718096;--ribbon-plus:#4299e1;--ribbon-minus:#fc8181;}}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{background:var(--bg);color:var(--text);font-family:'Inter',system-ui,sans-serif;height:100vh;display:flex;flex-direction:column}}
  header{{background:var(--panel);border-bottom:1px solid #2d3748;padding:12px 20px;display:flex;align-items:center;gap:12px}}
  header .logo{{font-size:1.3rem;font-weight:700;color:var(--accent);letter-spacing:-.5px}}
  header .subtitle{{color:var(--muted);font-size:.85rem}}
  #toolbar{{background:var(--panel);border-bottom:1px solid #2d3748;padding:8px 16px;display:flex;gap:16px;align-items:center;flex-wrap:wrap}}
  #toolbar label{{font-size:.8rem;color:var(--muted)}}
  #toolbar select,#toolbar input{{background:#2d3748;color:var(--text);border:1px solid #4a5568;border-radius:4px;padding:4px 8px;font-size:.8rem}}
  #toolbar button{{background:#2d3748;color:var(--accent);border:1px solid #4a5568;border-radius:4px;padding:4px 10px;font-size:.75rem;font-weight:600;cursor:pointer;font-family:inherit;transition:background .1s,color .1s}}
  #toolbar button:hover{{background:var(--accent);color:#0f1117}}
  #toolbar .sep{{width:1px;height:18px;background:#4a5568}}
  #stats{{font-size:.75rem;color:var(--muted);margin-left:auto}}
  #canvas-wrap{{flex:1;overflow:hidden;display:flex;flex-direction:column;align-items:center;justify-content:center;padding:20px}}
  svg{{width:100%;max-width:1400px;height:auto}}
  .tooltip{{position:fixed;background:#1a202c;border:1px solid #4a5568;border-radius:6px;padding:8px 12px;font-size:.8rem;pointer-events:none;opacity:0;transition:opacity .15s;max-width:300px;z-index:100}}
  #error-box{{display:none;background:#742a2a;border:1px solid #fc8181;border-radius:8px;padding:16px;margin:20px;text-align:center}}
</style>
</head>
<body>
<header>
  <span class="logo">🍄 MycoNote</span>
  <span class="subtitle">Synteny: <strong>{label1}</strong> vs <strong>{label2}</strong> — {block_count} blocks</span>
</header>
<div id="toolbar">
  <label>View: <select id="viewMode"><option value="ribbon">Ribbon</option><option value="dotplot">Dot plot</option></select></label>
  <label>Min identity: <input type="range" id="minId" min="0" max="100" value="70" step="1"/><span id="minIdVal">70%</span></label>
  <label>Colour by: <select id="colourMode"><option value="strand">Strand</option><option value="identity">Identity</option><option value="chrom">Chromosome</option></select></label>
  <label>Show gene names: <input type="checkbox" id="showNames" checked/></label>
  <label>Gene overlay: <input type="checkbox" id="showGenes"/></label>
  <span class="sep"></span>
  <button id="exportSvg" title="Download the current view as SVG">⬇ SVG</button>
  <button id="exportPng" title="Download the current view as PNG (2× resolution)">⬇ PNG</button>
  <div id="stats"></div>
</div>
<div id="canvas-wrap">
  <div id="error-box">⚠ Could not load D3.js from CDN. Please check your internet connection and reload the page.</div>
  <svg id="svg"></svg>
</div>
<div class="tooltip" id="tip"></div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js" crossorigin="anonymous"
  onerror="document.getElementById('error-box').style.display='block'"></script>
<script>
{js_sizes1}
{js_sizes2}
{js_blocks}
{js_names}
{js_genes1}
{js_genes2}

// ── Gene overlap lookup ─────────────────────────────────────────────────────
// For a block on `contig` spanning [qs, qe], return the display name of the
// gene with the largest overlap, or null if none overlap. genes_X contigs are
// pre-sorted by start, so we bail out as soon as a gene's start passes qe.
function labelForBlock(genesMap, contig, qs, qe) {{
  const ivs = genesMap[contig];
  if (!ivs) return null;
  let bestId = null;
  let bestOv = 0;
  for (const [gs, ge, gid] of ivs) {{
    if (gs > qe) break;
    if (ge < qs) continue;
    const ov = Math.min(ge, qe) - Math.max(gs, qs);
    if (ov > bestOv) {{ bestOv = ov; bestId = gid; }}
  }}
  return bestId ? (geneNames[bestId] || bestId) : null;
}}

// ── Colour helpers ──────────────────────────────────────────────────────────
const CHROM_COLOURS = d3.schemeTableau10.concat(d3.schemePastel1);
function chromColour(name, idx) {{ return CHROM_COLOURS[idx % CHROM_COLOURS.length]; }}
function identityColour(id) {{
  return d3.interpolateYlOrRd(id);  // 0=yellow, 1=red
}}
function strandColour(s) {{ return s === '+' ? 'var(--ribbon-plus)' : 'var(--ribbon-minus)'; }}

// ── Layout constants ────────────────────────────────────────────────────────
const MARGIN  = {{top:40, right:40, bottom:40, left:40}};
const TRACK_H = 22;
const RIBBON_AREA_H = 180;
const LABEL_PAD = 6;

function draw() {{
  const svg = d3.select('#svg');
  svg.selectAll('*').remove();

  const W  = parseInt(svg.style('width'))  || 1200;
  const innerW = W - MARGIN.left - MARGIN.right;

  // Filter blocks by min identity slider
  const minId = +document.getElementById('minId').value / 100;
  const colMode = document.getElementById('colourMode').value;
  const showNames = document.getElementById('showNames').checked;
  const showGenes = document.getElementById('showGenes').checked;
  const viewMode = document.getElementById('viewMode').value;
  const filtered = synBlocks.filter(b => b.id >= minId);
  document.getElementById('stats').textContent =
    `Showing ${{filtered.length}} / ${{synBlocks.length}} blocks (≥ ${{Math.round(minId*100)}}% identity)`;

  // ── Build chromosome order (sort by decreasing length) ──────────────────
  const chromsA = Object.entries(chrSizes_A).sort((a,b)=>b[1]-a[1]);
  const chromsB = Object.entries(chrSizes_B).sort((a,b)=>b[1]-a[1]);

  const totalA = chromsA.reduce((s,[,v])=>s+v, 0);
  const totalB = chromsB.reduce((s,[,v])=>s+v, 0);
  const GAP_PX = 4;

  const chromIdxA = Object.fromEntries(chromsA.map(([n],i) => [n, i]));
  const chromIdxB = Object.fromEntries(chromsB.map(([n],i) => [n, i]));
  const tip = document.getElementById('tip');

  // Ribbon and dot-plot paths share the chromosome ordering and colour
  // rules, but lay out on completely different coordinate systems, so each
  // mode builds its own offsets/scales below.

  if (viewMode === 'dotplot') {{
    drawDotPlot();
  }} else {{
    drawRibbon();
  }}

  // ══════════════════════════════════════════════════════════════════════════
  // Dot plot view
  // ══════════════════════════════════════════════════════════════════════════
  function drawDotPlot() {{
    // Left gutter for the Y axis (genome B contig labels), bottom gutter for
    // the X axis (genome A contig labels).
    const AXIS_LEFT = 110;
    const AXIS_BOTTOM = 60;
    const plotW = Math.max(100, innerW - AXIS_LEFT);
    const plotH = Math.max(100, Math.min(plotW, 640));  // cap for large screens

    const scaleX = (plotW - GAP_PX * (chromsA.length - 1)) / totalA;
    const scaleY = (plotH - GAP_PX * (chromsB.length - 1)) / totalB;

    const offX = {{}}; let xCursor = 0;
    for (const [name, len] of chromsA) {{ offX[name] = xCursor; xCursor += len * scaleX + GAP_PX; }}
    const offY = {{}}; let yCursor = 0;
    for (const [name, len] of chromsB) {{ offY[name] = yCursor; yCursor += len * scaleY + GAP_PX; }}

    const totalH = MARGIN.top + plotH + AXIS_BOTTOM + MARGIN.bottom;
    svg.attr('viewBox', `0 0 ${{W}} ${{totalH}}`);

    const root = svg.append('g')
      .attr('transform', `translate(${{MARGIN.left + AXIS_LEFT}},${{MARGIN.top}})`);

    // Plot frame
    root.append('rect')
      .attr('x', 0).attr('y', 0).attr('width', plotW).attr('height', plotH)
      .attr('fill', '#1a202c').attr('stroke', '#4a5568').attr('stroke-width', 0.5);

    // Contig gridlines + ticks — vertical (genome A) and horizontal (genome B)
    for (const [name, len] of chromsA) {{
      const x = offX[name] + len * scaleX;
      root.append('line')
        .attr('x1', x).attr('x2', x).attr('y1', 0).attr('y2', plotH)
        .attr('stroke', '#2d3748').attr('stroke-width', 0.5);
      // Diagonal label under the axis.
      root.append('text')
        .attr('x', offX[name] + len * scaleX / 2)
        .attr('y', plotH + 12)
        .attr('text-anchor', 'end').attr('font-size', 8).attr('fill', '#a0aec0')
        .attr('transform', `rotate(-45, ${{offX[name] + len * scaleX / 2}}, ${{plotH + 12}})`)
        .text(name);
    }}
    for (const [name, len] of chromsB) {{
      const y = offY[name] + len * scaleY;
      root.append('line')
        .attr('x1', 0).attr('x2', plotW).attr('y1', y).attr('y2', y)
        .attr('stroke', '#2d3748').attr('stroke-width', 0.5);
      root.append('text')
        .attr('x', -6).attr('y', offY[name] + len * scaleY / 2 + 3)
        .attr('text-anchor', 'end').attr('font-size', 8).attr('fill', '#a0aec0')
        .text(name);
    }}

    // Axis titles
    root.append('text')
      .attr('x', plotW / 2).attr('y', plotH + AXIS_BOTTOM - 10)
      .attr('text-anchor', 'middle').attr('font-size', 11).attr('font-weight', 'bold')
      .attr('fill', 'var(--accent)').text('{label1}');
    root.append('text')
      .attr('x', -AXIS_LEFT + 10).attr('y', plotH / 2)
      .attr('transform', `rotate(-90, ${{-AXIS_LEFT + 10}}, ${{plotH / 2}})`)
      .attr('text-anchor', 'middle').attr('font-size', 11).attr('font-weight', 'bold')
      .attr('fill', 'var(--accent)').text('{label2}');

    // One line segment per block. '+' strand → positive slope (query and
    // target both run forward); '-' strand → negative slope (target runs
    // backwards as query advances, the classic anti-diagonal).
    for (const b of filtered) {{
      if (!(b.qn in offX) || !(b.tn in offY)) continue;
      const x1 = offX[b.qn] + b.qs * scaleX;
      const x2 = offX[b.qn] + b.qe * scaleX;
      const y1 = offY[b.tn] + (b.st === '+' ? b.ts : b.te) * scaleY;
      const y2 = offY[b.tn] + (b.st === '+' ? b.te : b.ts) * scaleY;

      const colour = colMode === 'identity' ? identityColour(b.id)
                   : colMode === 'chrom'    ? chromColour(b.qn, chromIdxA[b.qn])
                   :                          strandColour(b.st);

      // Stroke width scales gently with identity so high-quality hits pop.
      const sw = 0.5 + 1.2 * b.id;

      root.append('line')
        .attr('x1', x1).attr('y1', y1).attr('x2', x2).attr('y2', y2)
        .attr('stroke', colour).attr('stroke-width', sw).attr('stroke-opacity', 0.9)
        .style('cursor', 'pointer')
        .on('mousemove', (event) => {{
          const nameA = labelForBlock(genes_A, b.qn, b.qs, b.qe);
          const nameB = labelForBlock(genes_B, b.tn, b.ts, b.te);
          const header = nameA || nameB || `${{b.qn}} ↔ ${{b.tn}}`;
          const geneLine = (nameA || nameB)
            ? `<span style="color:#a0aec0">Gene:</span> ${{nameA || '—'}} ↔ ${{nameB || '—'}}<br/>`
            : '';
          tip.style.opacity = 1;
          tip.style.left = (event.clientX + 12) + 'px';
          tip.style.top  = (event.clientY - 10) + 'px';
          tip.innerHTML  = `<b>${{header}}</b><br/>${{geneLine}}${{b.qn}}:${{b.qs.toLocaleString()}}–${{b.qe.toLocaleString()}}<br/>
            ↔ ${{b.tn}}:${{b.ts.toLocaleString()}}–${{b.te.toLocaleString()}}<br/>
            Strand: ${{b.st === '+' ? '➕ forward' : '➖ reverse'}}<br/>
            Identity: ${{(b.id*100).toFixed(1)}}%<br/>
            Length: ${{((b.qe-b.qs)/1000).toFixed(1)}} kb`;
        }})
        .on('mouseleave', () => {{ tip.style.opacity = 0; }});
    }}
  }}

  // ══════════════════════════════════════════════════════════════════════════
  // Ribbon view (original)
  // ══════════════════════════════════════════════════════════════════════════
  function drawRibbon() {{
  // Scale: pixel per base (use the larger genome to set scale)
  const scaleA = (innerW - GAP_PX * (chromsA.length - 1)) / totalA;
  const scaleB = (innerW - GAP_PX * (chromsB.length - 1)) / totalB;

  // x-offset for each chrom
  function buildOffsets(chroms, scale) {{
    const offsets = {{}};
    let x = 0;
    for (const [name, len] of chroms) {{
      offsets[name] = x;
      x += len * scale + GAP_PX;
    }}
    return offsets;
  }}
  const offA = buildOffsets(chromsA, scaleA);
  const offB = buildOffsets(chromsB, scaleB);

  // SVG total height
  const totalH = MARGIN.top + TRACK_H + RIBBON_AREA_H + TRACK_H + MARGIN.bottom + 30;
  svg.attr('viewBox', `0 0 ${{W}} ${{totalH}}`);

  const g = svg.append('g').attr('transform', `translate(${{MARGIN.left}},${{MARGIN.top}})`);

  const yA = 0;          // top of genome A track
  const yRibbon = TRACK_H + 8;
  const yB = TRACK_H + RIBBON_AREA_H + 16; // top of genome B track

  // ── Draw chromosome bars ────────────────────────────────────────────────
  function drawChroms(chroms, scale, offsets, yTop, label) {{
    for (const [name, len] of chroms) {{
      const x = offsets[name];
      const w = Math.max(1, len * scale);
      g.append('rect')
        .attr('x', x).attr('y', yTop)
        .attr('width', w).attr('height', TRACK_H)
        .attr('rx', 3)
        .attr('fill', '#2d3748').attr('stroke', '#4a5568').attr('stroke-width', 0.5);
      // Label (only if wide enough)
      if (w > 30) {{
        g.append('text')
          .attr('x', x + w/2).attr('y', yTop + TRACK_H/2 + 4)
          .attr('text-anchor','middle').attr('font-size', 9)
          .attr('fill', '#a0aec0').text(name);
      }}
    }}
    // Genome label on the left
    g.append('text')
      .attr('x', -LABEL_PAD).attr('y', yTop + TRACK_H/2 + 4)
      .attr('text-anchor','end').attr('font-size',11).attr('font-weight','bold')
      .attr('fill','var(--accent)').text(label);
  }}

  drawChroms(chromsA, scaleA, offA, yA, '{label1}');
  drawChroms(chromsB, scaleB, offB, yB, '{label2}');

  // ── Draw gene overlay ───────────────────────────────────────────────────
  // Tick marks on the chromosome bars for each annotated gene. Named genes
  // pop in the accent colour; unnamed ones get a dim fill. Ticks < 0.5 px
  // wide are skipped entirely (indistinguishable from a chromosome bar
  // edge at that zoom). A single SVG <g> wraps the overlay so it toggles
  // cleanly without re-walking the DOM.
  if (showGenes) {{
    function drawGenes(genesMap, offsets, scale, yTop, genomeLetter) {{
      const layer = g.append('g').attr('class', 'gene-overlay');
      for (const [contig, ivs] of Object.entries(genesMap)) {{
        const base = offsets[contig];
        if (base === undefined) continue;
        for (const [gs, ge, gid] of ivs) {{
          const w = Math.max(0.5, (ge - gs) * scale);
          if (w < 0.5) continue;
          const x = base + gs * scale;
          const named = geneNames[gid] !== undefined;
          layer.append('rect')
            .attr('x', x).attr('y', yTop + 2)
            .attr('width', w).attr('height', TRACK_H - 4)
            .attr('fill', named ? 'var(--accent)' : '#a0aec0')
            .attr('fill-opacity', named ? 0.85 : 0.35)
            .style('cursor', 'pointer')
            .on('mousemove', (event) => {{
              const name = geneNames[gid] || gid;
              tip.style.opacity = 1;
              tip.style.left = (event.clientX + 12) + 'px';
              tip.style.top  = (event.clientY - 10) + 'px';
              tip.innerHTML  = `<b>${{name}}</b><br/>${{genomeLetter}} · ${{contig}}:${{gs.toLocaleString()}}–${{ge.toLocaleString()}}<br/>
                Length: ${{(ge - gs).toLocaleString()}} bp`;
            }})
            .on('mouseleave', () => {{ tip.style.opacity = 0; }});
        }}
      }}
    }}
    drawGenes(genes_A, offA, scaleA, yA, '{label1}');
    drawGenes(genes_B, offB, scaleB, yB, '{label2}');
  }}

  // ── Draw ribbons ────────────────────────────────────────────────────────
  for (const b of filtered) {{
    if (!(b.qn in offA) || !(b.tn in offB)) continue;

    const x1 = offA[b.qn] + b.qs * scaleA;
    const x2 = offA[b.qn] + b.qe * scaleA;
    const x3 = offB[b.tn] + (b.st === '+' ? b.ts : b.te) * scaleB;
    const x4 = offB[b.tn] + (b.st === '+' ? b.te : b.ts) * scaleB;

    const y1 = yA + TRACK_H;
    const y2 = yB;
    const ymid = (y1 + y2) / 2;

    const colour = colMode === 'identity' ? identityColour(b.id)
                 : colMode === 'chrom'    ? chromColour(b.qn, chromIdxA[b.qn])
                 :                          strandColour(b.st);

    const path = `M${{x1}},${{y1}} C${{x1}},${{ymid}} ${{x3}},${{ymid}} ${{x3}},${{y2}}
                  L${{x4}},${{y2}} C${{x4}},${{ymid}} ${{x2}},${{ymid}} ${{x2}},${{y1}} Z`;

    g.append('path')
      .attr('d', path)
      .attr('fill', colour).attr('fill-opacity', 0.35)
      .attr('stroke', colour).attr('stroke-width', 0.4).attr('stroke-opacity', 0.7)
      .style('cursor','pointer')
      .on('mousemove', (event) => {{
        const nameA = labelForBlock(genes_A, b.qn, b.qs, b.qe);
        const nameB = labelForBlock(genes_B, b.tn, b.ts, b.te);
        const header = nameA || nameB || `${{b.qn}} ↔ ${{b.tn}}`;
        const geneLine = (nameA || nameB)
          ? `<span style="color:#a0aec0">Gene:</span> ${{nameA || '—'}} ↔ ${{nameB || '—'}}<br/>`
          : '';
        tip.style.opacity = 1;
        tip.style.left = (event.clientX + 12) + 'px';
        tip.style.top  = (event.clientY - 10) + 'px';
        tip.innerHTML  = `<b>${{header}}</b><br/>${{geneLine}}${{b.qn}}:${{b.qs.toLocaleString()}}–${{b.qe.toLocaleString()}}<br/>
          ↔ ${{b.tn}}:${{b.ts.toLocaleString()}}–${{b.te.toLocaleString()}}<br/>
          Strand: ${{b.st === '+' ? '➕ forward' : '➖ reverse'}}<br/>
          Identity: ${{(b.id*100).toFixed(1)}}%<br/>
          Length: ${{((b.qe-b.qs)/1000).toFixed(1)}} kb`;
      }})
      .on('mouseleave', () => {{ tip.style.opacity = 0; }});

    // Gene name label on large blocks — uses the top genome's overlapping gene.
    if (showNames && (x2 - x1) > 40) {{
      const mid = (x1 + x2) / 2;
      const labelText = labelForBlock(genes_A, b.qn, b.qs, b.qe);
      if (labelText) {{
        g.append('text')
          .attr('x', mid).attr('y', yA + TRACK_H - 3)
          .attr('text-anchor','middle').attr('font-size', 7).attr('fill','#e2e8f0')
          .attr('pointer-events','none')
          .text(labelText);
      }}
    }}
  }}
  }}  // end drawRibbon
}}  // end draw

// ── Controls ────────────────────────────────────────────────────────────────
document.getElementById('minId').addEventListener('input', function() {{
  document.getElementById('minIdVal').textContent = this.value + '%';
  draw();
}});
document.getElementById('viewMode').addEventListener('change', draw);
document.getElementById('colourMode').addEventListener('change', draw);
document.getElementById('showNames').addEventListener('change', draw);
document.getElementById('showGenes').addEventListener('change', draw);
window.addEventListener('resize', draw);

// ── Export helpers ──────────────────────────────────────────────────────────
// The in-page SVG uses CSS custom properties (var(--accent), var(--ribbon-plus),
// var(--ribbon-minus)) on fill/stroke presentation attributes. Those are
// resolved by the page's stylesheet at render time but do NOT survive a raw
// serialisation — a standalone .svg file has no :root and browsers won't
// resolve `var()` inside presentation attributes. So before export we rewrite
// every var() reference to the literal hex so the exported SVG/PNG is truly
// self-contained.
function serialiseSvg() {{
  const src = document.getElementById('svg');
  const clone = src.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  // Paint a background rect first so exported files aren't transparent.
  const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  bg.setAttribute('width', '100%');
  bg.setAttribute('height', '100%');
  bg.setAttribute('fill', '#0f1117');
  clone.insertBefore(bg, clone.firstChild);

  let xml = new XMLSerializer().serializeToString(clone);
  const palette = {{
    '--accent': '#00d4aa',
    '--ribbon-plus': '#4299e1',
    '--ribbon-minus': '#fc8181',
    '--text': '#e2e8f0',
    '--muted': '#718096',
  }};
  for (const [name, hex] of Object.entries(palette)) {{
    xml = xml.split('var(' + name + ')').join(hex);
  }}
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + xml;
}}

function triggerDownload(blob, filename) {{
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename; a.style.display = 'none';
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 500);
}}

document.getElementById('exportSvg').addEventListener('click', () => {{
  const xml = serialiseSvg();
  triggerDownload(new Blob([xml], {{type: 'image/svg+xml;charset=utf-8'}}), 'synteny.svg');
}});

document.getElementById('exportPng').addEventListener('click', () => {{
  const xml = serialiseSvg();
  const svgEl = document.getElementById('svg');
  const rect = svgEl.getBoundingClientRect();
  // 2× resolution for Retina-like sharpness; capped to avoid >16k canvas limit.
  const scale = Math.min(2, 16000 / Math.max(rect.width, rect.height, 1));
  const W = Math.round(rect.width  * scale);
  const H = Math.round(rect.height * scale);
  const img = new Image();
  const svgUrl = URL.createObjectURL(new Blob([xml], {{type: 'image/svg+xml;charset=utf-8'}}));
  img.onload = () => {{
    const canvas = document.createElement('canvas');
    canvas.width = W; canvas.height = H;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#0f1117';
    ctx.fillRect(0, 0, W, H);
    ctx.drawImage(img, 0, 0, W, H);
    URL.revokeObjectURL(svgUrl);
    canvas.toBlob(b => {{ if (b) triggerDownload(b, 'synteny.png'); }}, 'image/png');
  }};
  img.onerror = () => {{
    URL.revokeObjectURL(svgUrl);
    alert('PNG export failed — try SVG export instead.');
  }};
  img.src = svgUrl;
}});

// Initial draw
if (typeof d3 !== 'undefined') {{ draw(); }}
else {{ setTimeout(() => {{ if (typeof d3 !== 'undefined') draw(); }}, 500); }}
</script>
</body>
</html>
"#,
        title = title,
        label1 = label1,
        label2 = label2,
        block_count = block_count,
        js_sizes1 = js_sizes1,
        js_sizes2 = js_sizes2,
        js_blocks = js_blocks,
        js_names = js_names,
        js_genes1 = js_genes1,
        js_genes2 = js_genes2,
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn generate_synteny(config: &SyntenyConfig) -> Result<()> {
    // 1. Get chromosome sizes from GFF3 (works even without FASTA)
    println!("📐 Reading chromosome sizes from GFF3 files…");
    let sizes1 = chrom_sizes_from_gff(&config.gff1)?;
    let sizes2 = chrom_sizes_from_gff(&config.gff2)?;
    println!("   Genome A: {} sequences", sizes1.len());
    println!("   Genome B: {} sequences", sizes2.len());

    // 2. Run alignment if FASTA files are provided. Use a per-run temp
    //    directory (nanosecond-suffixed) so concurrent invocations don't
    //    clobber each other's PAF output, and clean it up unless the user
    //    asked to keep it for debugging.
    let tmp_dir = if config.fasta1.is_some() && config.fasta2.is_some() {
        let suffix = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        Some(std::env::temp_dir().join(format!("myconote_synteny_{}", suffix)))
    } else {
        None
    };

    let blocks = if let (Some(fa1), Some(fa2), Some(tmp)) =
        (&config.fasta1, &config.fasta2, &tmp_dir)
    {
        std::fs::create_dir_all(tmp).map_err(MycoNoteError::Io)?;
        let paf_path = run_minimap2(fa1, fa2, tmp, config.threads)?;
        let raw = parse_paf(&paf_path, config.min_block_len)?;
        println!(
            "   {} raw PAF hits ≥ {} bp",
            raw.len(),
            config.min_block_len
        );
        let chained = chain_blocks(raw, config.chain_gap);
        if config.chain_gap > 0 {
            println!(
                "   {} chained blocks (gap ≤ {} bp)",
                chained.len(),
                config.chain_gap
            );
        }
        chained
    } else {
        eprintln!("  ℹ  No FASTA files provided — rendering chromosome layout only (no ribbons).");
        eprintln!("     Add --fasta1 <genome_a.fa> --fasta2 <genome_b.fa> to enable alignment.");
        Vec::new()
    };

    // 3. Extract all gene intervals. Used both for the on-track gene overlay
    //    and to resolve PAF-block labels via genomic overlap (blocks key on
    //    contig, not gene ID — without this, `--names` does nothing).
    let genes1 = genes_from_gff(&config.gff1)?;
    let genes2 = genes_from_gff(&config.gff2)?;
    let total_genes: usize = genes1.values().map(|v| v.len()).sum::<usize>()
        + genes2.values().map(|v| v.len()).sum::<usize>();
    println!(
        "   {} genes extracted for overlay + block labels",
        total_genes
    );

    // 4. Render HTML
    println!("🎨 Rendering synteny diagram…");
    let html = render_synteny_html(&blocks, &sizes1, &sizes2, &genes1, &genes2, config);

    let mut out = std::fs::File::create(&config.output).map_err(MycoNoteError::Io)?;
    out.write_all(html.as_bytes()).map_err(MycoNoteError::Io)?;

    println!("✓ Synteny viewer: {}", config.output.display());
    println!("  Open in any modern browser — no server required.");

    // 5. Clean up the temp PAF unless the user asked to keep it.
    if let Some(tmp) = &tmp_dir {
        if config.keep_paf {
            println!("  PAF retained at: {}", tmp.display());
        } else if let Err(e) = std::fs::remove_dir_all(tmp) {
            eprintln!("  ⚠  Could not remove temp dir {}: {}", tmp.display(), e);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blk(
        qn: &str,
        qs: u64,
        qe: u64,
        tn: &str,
        ts: u64,
        te: u64,
        st: char,
        rm: u64,
    ) -> SyntenyBlock {
        let block_len = qe - qs;
        SyntenyBlock {
            query_name: qn.to_string(),
            query_start: qs,
            query_end: qe,
            query_len: 1_000_000,
            target_name: tn.to_string(),
            target_start: ts,
            target_end: te,
            target_len: 1_000_000,
            strand: st,
            identity: rm as f64 / block_len as f64,
            block_len,
            residue_matches: rm,
        }
    }

    #[test]
    fn chain_merges_two_close_plus_strand_hits() {
        let blocks = vec![
            blk("chr1", 100, 1_100, "ctgA", 500, 1_500, '+', 950),
            blk("chr1", 1_200, 2_200, "ctgA", 1_600, 2_600, '+', 980),
        ];
        let chained = chain_blocks(blocks, 10_000);
        assert_eq!(chained.len(), 1);
        let m = &chained[0];
        assert_eq!(m.query_start, 100);
        assert_eq!(m.query_end, 2_200);
        assert_eq!(m.target_start, 500);
        assert_eq!(m.target_end, 2_600);
        assert_eq!(m.residue_matches, 1_930);
        assert_eq!(m.block_len, 2_000);
        assert!((m.identity - 0.965).abs() < 1e-9);
    }

    #[test]
    fn chain_merges_two_close_minus_strand_hits() {
        // Reverse-strand chain: as query_start increases, target runs backwards.
        let blocks = vec![
            blk("chr1", 100, 1_100, "ctgA", 8_000, 9_000, '-', 900),
            blk("chr1", 1_200, 2_200, "ctgA", 6_500, 7_500, '-', 900),
        ];
        let chained = chain_blocks(blocks, 10_000);
        assert_eq!(chained.len(), 1);
        let m = &chained[0];
        assert_eq!(m.strand, '-');
        assert_eq!(m.query_start, 100);
        assert_eq!(m.query_end, 2_200);
        assert_eq!(m.target_start, 6_500);
        assert_eq!(m.target_end, 9_000);
    }

    #[test]
    fn chain_does_not_merge_across_large_gap() {
        let blocks = vec![
            blk("chr1", 100, 1_100, "ctgA", 500, 1_500, '+', 950),
            blk("chr1", 500_000, 501_000, "ctgA", 500_400, 501_400, '+', 950),
        ];
        let chained = chain_blocks(blocks, 10_000);
        assert_eq!(chained.len(), 2);
    }

    #[test]
    fn chain_does_not_merge_across_contigs_or_strands() {
        let blocks = vec![
            blk("chr1", 100, 1_100, "ctgA", 500, 1_500, '+', 950),
            blk("chr1", 1_200, 2_200, "ctgB", 1_600, 2_600, '+', 980),
            blk("chr1", 2_300, 3_300, "ctgA", 2_700, 3_700, '-', 970),
        ];
        let chained = chain_blocks(blocks, 10_000);
        assert_eq!(chained.len(), 3);
    }

    #[test]
    fn chain_rejects_non_monotonic_target_on_plus() {
        // Same query direction but target coords go backwards — not colinear.
        let blocks = vec![
            blk("chr1", 100, 1_100, "ctgA", 5_000, 6_000, '+', 950),
            blk("chr1", 1_200, 2_200, "ctgA", 1_000, 2_000, '+', 950),
        ];
        let chained = chain_blocks(blocks, 10_000);
        assert_eq!(chained.len(), 2);
    }

    #[test]
    fn chain_gap_zero_is_passthrough() {
        let blocks = vec![
            blk("chr1", 100, 1_100, "ctgA", 500, 1_500, '+', 950),
            blk("chr1", 1_200, 2_200, "ctgA", 1_600, 2_600, '+', 980),
        ];
        let chained = chain_blocks(blocks.clone(), 0);
        assert_eq!(chained.len(), 2);
    }

    #[test]
    fn parse_paf_reads_basic_record() {
        use std::io::Write;
        let mut tmp = std::env::temp_dir();
        tmp.push("myconote_synteny_parse_test.paf");
        let mut f = std::fs::File::create(&tmp).unwrap();
        // Two tab-separated PAF records + a short malformed line that must be skipped.
        writeln!(
            f,
            "qseq\t10000\t100\t1100\t+\ttseq\t20000\t500\t1500\t950\t1000\t60"
        )
        .unwrap();
        writeln!(
            f,
            "qseq\t10000\t2000\t2500\t-\ttseq\t20000\t3000\t3500\t480\t500\t60"
        )
        .unwrap();
        writeln!(f, "too\tshort").unwrap();
        drop(f);

        let blocks = parse_paf(&tmp, 100).unwrap();
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].query_name, "qseq");
        assert_eq!(blocks[0].strand, '+');
        assert_eq!(blocks[0].residue_matches, 950);
        assert!((blocks[0].identity - 0.95).abs() < 1e-9);
        assert_eq!(blocks[1].strand, '-');

        // min_len filter drops anything smaller than threshold.
        let filtered = parse_paf(&tmp, 750).unwrap();
        assert_eq!(filtered.len(), 1);

        let _ = std::fs::remove_file(&tmp);
    }
}
