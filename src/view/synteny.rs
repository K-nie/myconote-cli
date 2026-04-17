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
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PAF block representation
// ─────────────────────────────────────────────────────────────────────────────

/// One syntenic block parsed from a PAF line.
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
    pub identity: f64, // 0.0–1.0
    pub block_len: u64,
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
        });
    }

    Ok(blocks)
}

// ─────────────────────────────────────────────────────────────────────────────
// minimap2 runner
// ─────────────────────────────────────────────────────────────────────────────

/// Run minimap2 in asm-to-asm mode and return path to the PAF output file.
fn run_minimap2(fasta1: &Path, fasta2: &Path, out_dir: &Path) -> Result<PathBuf> {
    // Check minimap2 is available
    let mm2 = which::which("minimap2").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "minimap2 not found in PATH. Install it with: conda install -c bioconda minimap2\n\
             or: brew install minimap2"
                .to_string(),
        )
    })?;

    let paf_path = out_dir.join("synteny_alignment.paf");

    println!("  Running minimap2 asm-to-asm alignment…");
    let status = Command::new(&mm2)
        .args([
            "-cx",
            "asm5", // asm-to-asm preset (≥5% divergence)
            "--cs", // include cs tag for identity calculation
            "-t",
            "4", // 4 threads
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

/// Extract gene intervals from a GFF3, keeping only genes whose ID is present
/// in `name_keys` (so the emitted JSON stays small and every serialised gene
/// has a display name). Intervals are sorted by start per contig for the JS
/// binary-search-style overlap lookup.
fn named_genes_from_gff(gff_path: &Path, name_keys: &HashMap<String, String>) -> Result<GenesByContig> {
    use crate::parser::gff::GFFReader;
    let mut genes: GenesByContig = HashMap::new();
    if name_keys.is_empty() {
        return Ok(genes);
    }
    for result in GFFReader::from_path(gff_path)? {
        let rec = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        if rec.feature_type != "gene" {
            continue;
        }
        let Some(id) = rec.id().cloned() else { continue };
        if !name_keys.contains_key(&id) {
            continue;
        }
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
  <label>Min identity: <input type="range" id="minId" min="0" max="100" value="70" step="1"/><span id="minIdVal">70%</span></label>
  <label>Colour by: <select id="colourMode"><option value="strand">Strand</option><option value="identity">Identity</option><option value="chrom">Chromosome</option></select></label>
  <label>Show gene names: <input type="checkbox" id="showNames" checked/></label>
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
  const filtered = synBlocks.filter(b => b.id >= minId);
  document.getElementById('stats').textContent =
    `Showing ${{filtered.length}} / ${{synBlocks.length}} blocks (≥ ${{Math.round(minId*100)}}% identity)`;

  // ── Build chromosome order (sort by decreasing length) ──────────────────
  const chromsA = Object.entries(chrSizes_A).sort((a,b)=>b[1]-a[1]);
  const chromsB = Object.entries(chrSizes_B).sort((a,b)=>b[1]-a[1]);

  const totalA = chromsA.reduce((s,[,v])=>s+v, 0);
  const totalB = chromsB.reduce((s,[,v])=>s+v, 0);
  const GAP_PX = 4;

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

  const chromIdxA = Object.fromEntries(chromsA.map(([n],i) => [n, i]));
  const chromIdxB = Object.fromEntries(chromsB.map(([n],i) => [n, i]));

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

  // ── Draw ribbons ────────────────────────────────────────────────────────
  const tip = document.getElementById('tip');

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
}}

// ── Controls ────────────────────────────────────────────────────────────────
document.getElementById('minId').addEventListener('input', function() {{
  document.getElementById('minIdVal').textContent = this.value + '%';
  draw();
}});
document.getElementById('colourMode').addEventListener('change', draw);
document.getElementById('showNames').addEventListener('change', draw);
window.addEventListener('resize', draw);

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

    // 2. Run alignment if FASTA files are provided
    let blocks = if let (Some(fa1), Some(fa2)) = (&config.fasta1, &config.fasta2) {
        let tmp_dir = std::env::temp_dir().join("myconote_synteny");
        std::fs::create_dir_all(&tmp_dir).map_err(MycoNoteError::Io)?;
        let paf_path = run_minimap2(fa1, fa2, &tmp_dir)?;
        let blocks = parse_paf(&paf_path, config.min_block_len)?;
        println!(
            "   {} syntenic blocks ≥ {} bp",
            blocks.len(),
            config.min_block_len
        );
        blocks
    } else {
        eprintln!("  ℹ  No FASTA files provided — rendering chromosome layout only (no ribbons).");
        eprintln!("     Add --fasta1 <genome_a.fa> --fasta2 <genome_b.fa> to enable alignment.");
        Vec::new()
    };

    // 3. Extract named gene intervals so block labels can resolve gene IDs
    //    via genomic overlap (PAF blocks key on contig, not on gene ID).
    let genes1 = named_genes_from_gff(&config.gff1, &config.names)?;
    let genes2 = named_genes_from_gff(&config.gff2, &config.names)?;
    let labelable: usize = genes1.values().map(|v| v.len()).sum::<usize>()
        + genes2.values().map(|v| v.len()).sum::<usize>();
    if !config.names.is_empty() {
        println!(
            "   {} named genes available for ribbon labels",
            labelable
        );
    }

    // 4. Render HTML
    println!("🎨 Rendering synteny diagram…");
    let html = render_synteny_html(&blocks, &sizes1, &sizes2, &genes1, &genes2, config);

    let mut out = std::fs::File::create(&config.output).map_err(MycoNoteError::Io)?;
    out.write_all(html.as_bytes()).map_err(MycoNoteError::Io)?;

    println!("✓ Synteny viewer: {}", config.output.display());
    println!("  Open in any modern browser — no server required.");

    Ok(())
}
