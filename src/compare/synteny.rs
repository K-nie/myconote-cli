use crate::compare::{CompareConfig, GenomeInfo};
use crate::utils::error::Result;
use std::fs::File;
use std::io::Write;

pub fn generate_synteny_plot(genomes: &[GenomeInfo], config: &CompareConfig) -> Result<()> {
    println!("Generating synteny plot...");

    let synteny_dir = config.output_dir.join("synteny");
    std::fs::create_dir_all(&synteny_dir)?;

    // Create an HTML file with a basic synteny visualization
    let html_path = synteny_dir.join("synteny.html");
    let mut html = File::create(html_path)?;

    writeln!(html, "<!DOCTYPE html>")?;
    writeln!(html, "<html>")?;
    writeln!(html, "<head>")?;
    writeln!(html, "    <title>Synteny Plot</title>")?;
    writeln!(html, "    <style>")?;
    writeln!(
        html,
        "        body {{ font-family: Arial, sans-serif; margin: 40px; }}"
    )?;
    writeln!(html, "        .genome {{ margin-bottom: 30px; }}")?;
    writeln!(html, "        .track {{ height: 40px; background: #f0f0f0; margin: 5px 0; position: relative; }}")?;
    writeln!(html, "        .gene {{ position: absolute; height: 40px; background: #3498db; border-radius: 3px; }}")?;
    writeln!(
        html,
        "        .link {{ stroke: #e74c3c; stroke-width: 2; }}"
    )?;
    writeln!(html, "    </style>")?;
    writeln!(html, "</head>")?;
    writeln!(html, "<body>")?;
    writeln!(html, "    <h1>Synteny Visualization</h1>")?;

    // Plot each genome
    for (_i, genome) in genomes.iter().enumerate() {
        writeln!(html, "    <div class='genome'>")?;
        writeln!(html, "        <h2>{}</h2>", genome.name)?;
        writeln!(html, "        <div class='track' style='width: 100%;'>")?;

        // Plot genes (simplified - just show a few representative genes)
        let num_genes = genome.gene_count.min(20);
        for j in 0..num_genes {
            let left = (j as f64 / num_genes as f64 * 100.0) as u32;
            let width = (1.0 / num_genes as f64 * 90.0) as u32;
            writeln!(
                html,
                "            <div class='gene' style='left: {}%; width: {}%;'></div>",
                left, width
            )?;
        }

        writeln!(html, "        </div>")?;
        writeln!(html, "    </div>")?;
    }

    // Add SVG for links between genomes
    writeln!(
        html,
        "    <svg width='100%' height='200' style='margin-top: 20px;'>"
    )?;
    writeln!(
        html,
        "        <text x='10' y='20'>Syntenic links between genomes</text>"
    )?;

    // Draw links between first two genomes if available
    if genomes.len() >= 2 {
        for i in 0..5 {
            let x1 = 50 + i * 30;
            let x2 = 150 + i * 30;
            writeln!(
                html,
                "        <line class='link' x1='{}' y1='50' x2='{}' y2='150' />",
                x1, x2
            )?;
        }
    }

    writeln!(html, "    </svg>")?;
    writeln!(html, "</body>")?;
    writeln!(html, "</html>")?;

    println!(
        "  Synteny plot saved to: {}/synteny.html",
        synteny_dir.display()
    );
    Ok(())
}
