/// Circular genome map — outputs a PNG file using the `plotters` crate.
///
/// Draws chromosomes/scaffolds as arcs around a circle, with gene density
/// shading and forward/reverse-strand gene tracks as coloured radial bars.
///
/// Layout (outermost → innermost):
///   1. Chromosome arcs — labelled, sized proportionally to sequence length
///   2. Gene density heat-map ring (window-based, blue gradient)
///   3. Forward-strand gene bars (blue)
///   4. Reverse-strand gene bars (orange)
///   5. CDS bars (green, slightly inner)

use crate::parser::gff::GFFReader;
use crate::utils::error::Result;
use crate::plot::PlotConfig;
use plotters::prelude::*;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn polar_to_cart(cx: f64, cy: f64, r: f64, angle_rad: f64) -> (f64, f64) {
    (cx + r * angle_rad.cos(), cy + r * angle_rad.sin())
}

struct ChromLayout {
    pub map:      HashMap<String, (f64, f64, u64)>, // id → (start_angle, end_angle, len)
    pub order:    Vec<String>,
    pub total_bp: u64,
}

impl ChromLayout {
    fn build(sizes: &HashMap<String, u64>) -> Self {
        let total_bp: u64 = sizes.values().sum();
        let n  = sizes.len() as f64;
        let gap = 0.002_f64 * n;
        let usable = 1.0 - gap;

        let mut order: Vec<String> = sizes.keys().cloned().collect();
        order.sort_by(|a, b| sizes[b].cmp(&sizes[a]));

        let mut map   = HashMap::new();
        let mut angle = -PI / 2.0;

        for chrom in &order {
            let len   = sizes[chrom];
            let frac  = len as f64 / total_bp as f64;
            let arc   = frac * usable * 2.0 * PI;
            map.insert(chrom.clone(), (angle, angle + arc, len));
            angle += arc + 0.002 * 2.0 * PI;
        }

        Self { map, order, total_bp }
    }

    fn pos_to_angle(&self, chrom: &str, pos: u64) -> Option<f64> {
        let (sa, ea, len) = *self.map.get(chrom)?;
        Some(sa + (pos as f64 / len as f64) * (ea - sa))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn plot<P: AsRef<Path>>(gff_path: P, config: &PlotConfig) -> Result<()> {
    println!("Generating circular genome map: {}", config.output_path);

    // ── Parse GFF ─────────────────────────────────────────────────────────
    let mut chrom_sizes: HashMap<String, u64> = HashMap::new();
    let mut genes:  Vec<(String, u64, u64, char)> = Vec::new(); // (seqid,start,end,strand)
    let mut cdss:   Vec<(String, u64, u64)>       = Vec::new();

    let reader = GFFReader::from_path(gff_path)?;
    for result in reader {
        let rec = result?;
        let e = chrom_sizes.entry(rec.seqid.clone()).or_insert(0);
        if rec.end > *e { *e = rec.end; }

        match rec.feature_type.as_str() {
            "gene" if config.show_genes => {
                genes.push((rec.seqid, rec.start, rec.end, rec.strand));
            }
            "CDS" if config.show_cds => {
                cdss.push((rec.seqid, rec.start, rec.end));
            }
            _ => {}
        }
    }

    if chrom_sizes.is_empty() {
        println!("No features found — cannot generate circular plot.");
        return Ok(());
    }

    let layout = ChromLayout::build(&chrom_sizes);
    let n_genes = genes.len();
    let arc_steps = 300usize;

    // ── Canvas ────────────────────────────────────────────────────────────
    let side = config.width.min(config.height);
    let root = BitMapBackend::new(&config.output_path, (config.width, config.height))
        .into_drawing_area();
    root.fill(&RGBColor(18, 18, 28))?;

    let cx = config.width  as f64 / 2.0;
    let cy = config.height as f64 / 2.0;
    let r_out     = side as f64 * 0.42;
    let r_chr_in  = side as f64 * 0.395;
    let r_dens_o  = side as f64 * 0.375;
    let r_dens_i  = side as f64 * 0.35;
    let r_fw_o    = side as f64 * 0.335;
    let r_fw_i    = side as f64 * 0.305;
    let r_rv_o    = side as f64 * 0.29;
    let r_rv_i    = side as f64 * 0.26;
    let r_cds_o   = side as f64 * 0.245;
    let r_cds_i   = side as f64 * 0.22;

    // ── Title ─────────────────────────────────────────────────────────────
    let title = config.title.as_deref().unwrap_or("Circular Genome Map");
    root.draw(&Text::new(
        title.to_string(),
        ((cx - title.len() as f64 * 5.5) as i32, 22),
        ("sans-serif", 22).into_font().color(&WHITE),
    ))?;
    let sub = format!("{} contigs  ·  {} genes  ·  {:.1} Mbp",
        chrom_sizes.len(), n_genes, layout.total_bp as f64 / 1e6);
    root.draw(&Text::new(
        sub,
        ((cx - 100.0) as i32, 50),
        ("sans-serif", 12).into_font().color(&RGBColor(160, 160, 175)),
    ))?;

    // ── Chromosome arcs ───────────────────────────────────────────────────
    for chrom in &layout.order {
        let (sa, ea, _) = layout.map[chrom];
        let col = chrom_color(chrom, &layout.order);

        let mut pts: Vec<(i32, i32)> = (0..=arc_steps)
            .map(|i| {
                let a = sa + (ea - sa) * i as f64 / arc_steps as f64;
                let (x, y) = polar_to_cart(cx, cy, r_out, a);
                (x as i32, y as i32)
            })
            .collect();
        pts.extend((0..=arc_steps).rev().map(|i| {
            let a = sa + (ea - sa) * i as f64 / arc_steps as f64;
            let (x, y) = polar_to_cart(cx, cy, r_chr_in, a);
            (x as i32, y as i32)
        }));
        root.draw(&Polygon::new(pts, col.filled()))?;

        // Label
        let mid = (sa + ea) / 2.0;
        if (ea - sa).to_degrees() > 2.5 {
            let (lx, ly) = polar_to_cart(cx, cy, r_out + 16.0, mid);
            let label = short_label(chrom);
            root.draw(&Text::new(
                label,
                (lx as i32 - 8, ly as i32 - 5),
                ("sans-serif", 9).into_font().color(&WHITE),
            ))?;
        }
    }

    // ── Gene density ring ─────────────────────────────────────────────────
    let window_bp: u64 = (layout.total_bp / 2000).max(5000);
    let mut density: HashMap<String, Vec<u32>> = HashMap::new();
    for chrom in &layout.order {
        let n = (chrom_sizes[chrom] / window_bp + 1) as usize;
        density.insert(chrom.clone(), vec![0u32; n]);
    }
    for (seqid, start, _, _) in &genes {
        if let Some(v) = density.get_mut(seqid) {
            let idx = (*start / window_bp) as usize;
            if idx < v.len() { v[idx] += 1; }
        }
    }
    let max_d = density.values().flat_map(|v| v.iter()).copied().max().unwrap_or(1).max(1) as f64;

    for (chrom, buckets) in &density {
        let (sa, ea, _) = layout.map[chrom];
        let n = buckets.len();
        for (i, &cnt) in buckets.iter().enumerate() {
            if cnt == 0 { continue; }
            let a0 = sa + (i as f64 / n as f64) * (ea - sa);
            let a1 = sa + ((i + 1) as f64 / n as f64) * (ea - sa);
            let steps = (((a1 - a0).abs() * arc_steps as f64 / (2.0 * PI)) as usize).max(2);
            let heat  = (cnt as f64 / max_d * 200.0) as u8 + 55;
            let col   = RGBColor(heat / 3, heat / 3, heat);
            let mut pts: Vec<(i32, i32)> = (0..=steps)
                .map(|s| { let a = a0 + (a1-a0)*s as f64/steps as f64; let (x,y)=polar_to_cart(cx,cy,r_dens_o,a); (x as i32,y as i32) })
                .collect();
            pts.extend((0..=steps).rev().map(|s| { let a = a0 + (a1-a0)*s as f64/steps as f64; let (x,y)=polar_to_cart(cx,cy,r_dens_i,a); (x as i32,y as i32) }));
            root.draw(&Polygon::new(pts, col.filled()))?;
        }
    }

    // ── Gene tracks ───────────────────────────────────────────────────────
    let draw_track = |root: &DrawingArea<BitMapBackend, plotters::coord::Shift>,
                      seqid: &str, start: u64, end: u64,
                      ro: f64, ri: f64, col: RGBColor| -> core::result::Result<(), Box<dyn std::error::Error + Send + Sync>> {
        let Some(a0) = layout.pos_to_angle(seqid, start) else { return Ok(()); };
        let Some(a1) = layout.pos_to_angle(seqid, end)   else { return Ok(()); };
        if (a1 - a0).abs() < 0.001 { return Ok(()); }
        let steps = (((a1-a0).abs() * arc_steps as f64 / (2.0*PI)) as usize).max(2);
        let mut pts: Vec<(i32,i32)> = (0..=steps).map(|s| { let a=a0+(a1-a0)*s as f64/steps as f64; let (x,y)=polar_to_cart(cx,cy,ro,a); (x as i32,y as i32) }).collect();
        pts.extend((0..=steps).rev().map(|s| { let a=a0+(a1-a0)*s as f64/steps as f64; let (x,y)=polar_to_cart(cx,cy,ri,a); (x as i32,y as i32) }));
        root.draw(&Polygon::new(pts, col.mix(0.75).filled()))?;
        Ok(())
    };

    for (seqid, start, end, strand) in &genes {
        let (ro, ri, col) = if *strand == '+' {
            (r_fw_o, r_fw_i, RGBColor(70, 130, 220))
        } else {
            (r_rv_o, r_rv_i, RGBColor(220, 110, 50))
        };
        let _ = draw_track(&root, seqid, *start, *end, ro, ri, col);
    }

    if config.show_cds {
        for (seqid, start, end) in &cdss {
            let _ = draw_track(&root, seqid, *start, *end, r_cds_o, r_cds_i, RGBColor(80, 200, 100));
        }
    }

    // ── Legend ────────────────────────────────────────────────────────────
    let items: &[(&str, RGBColor)] = &[
        ("+ strand genes", RGBColor(70, 130, 220)),
        ("− strand genes", RGBColor(220, 110, 50)),
        ("CDS",            RGBColor(80, 200, 100)),
        ("Gene density",   RGBColor(80, 80, 200)),
    ];
    let mut ly = config.height as i32 - 90;
    for (label, col) in items {
        root.draw(&Rectangle::new([(14, ly), (30, ly + 10)], col.filled()))?;
        root.draw(&Text::new(label.to_string(), (34, ly), ("sans-serif", 11).into_font().color(&WHITE)))?;
        ly += 17;
    }

    root.present()?;
    println!("✓ Circular plot saved: {} ({}×{} PNG)", config.output_path, config.width, config.height);
    Ok(())
}

fn chrom_color(chrom: &str, order: &[String]) -> RGBColor {
    const P: &[(u8,u8,u8)] = &[
        (86,180,233),(230,159,0),(0,158,115),(240,228,66),(0,114,178),
        (213,94,0),(204,121,167),(148,103,189),(23,190,207),(188,189,34),
        (31,119,180),(255,127,14),(44,160,44),(214,39,40),(148,103,189),
    ];
    let i = order.iter().position(|c| c == chrom).unwrap_or(0) % P.len();
    RGBColor(P[i].0, P[i].1, P[i].2)
}

fn short_label(s: &str) -> String {
    // "scaffold_001" → "001", "chr1" → "chr1", "NODE_1_length_..." → "1"
    if let Some(n) = s.split('_').last() { n.to_string() } else { s.to_string() }
}
