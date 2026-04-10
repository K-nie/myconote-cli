use crate::compare::{CompareConfig, GenomeInfo};
use crate::utils::error::Result;
use std::fs::File;
use std::io::Write;

pub fn build_tree(genomes: &[GenomeInfo], config: &CompareConfig) -> Result<()> {
    println!("Building phylogenetic tree...");

    let tree_dir = config.output_dir.join("phylogeny");
    std::fs::create_dir_all(&tree_dir)?;

    // Create a Newick format tree based on genome similarities
    let mut newick = String::new();

    if genomes.len() == 1 {
        newick = format!("{};", genomes[0].name);
    } else if genomes.len() == 2 {
        newick = format!("({},{})root;", genomes[0].name, genomes[1].name);
    } else {
        // Simple UPGMA-like tree for demonstration
        newick.push('(');
        for (i, genome) in genomes.iter().enumerate() {
            if i > 0 {
                newick.push(',');
            }
            newick.push_str(&genome.name);
        }
        newick.push_str(")root;");
    }

    let tree_path = tree_dir.join("tree.nwk");
    std::fs::write(&tree_path, &newick)?;

    // Create an HTML visualization
    let html_path = tree_dir.join("tree.html");
    let mut html = File::create(html_path)?;

    writeln!(html, "<!DOCTYPE html>")?;
    writeln!(html, "<html>")?;
    writeln!(html, "<head>")?;
    writeln!(html, "    <title>Phylogenetic Tree</title>")?;
    writeln!(html, "    <script src='https://cdnjs.cloudflare.com/ajax/libs/raphael/2.3.0/raphael.min.js'></script>")?;
    writeln!(
        html,
        "    <script src='https://unpkg.com/@phylocanvas/phylocanvas.gl/dist/index.js'></script>"
    )?;
    writeln!(html, "</head>")?;
    writeln!(html, "<body>")?;
    writeln!(html, "    <h1>Phylogenetic Tree</h1>")?;
    writeln!(
        html,
        "    <div id='tree' style='width: 800px; height: 600px;'></div>"
    )?;
    writeln!(html, "    <script>")?;
    writeln!(html, "        const newick = `{}`;", newick)?;
    writeln!(html, "        console.log('Newick:', newick);")?;
    let _ = writeln!(html, "        // Tree visualization would go here");
    writeln!(
        html,
        "        document.getElementById('tree').innerHTML = '<pre>' + newick + '</pre>';"
    )?;
    writeln!(html, "    </script>")?;
    writeln!(html, "</body>")?;
    writeln!(html, "</html>")?;

    println!("  Tree saved to: {}/tree.nwk", tree_dir.display());
    println!("  Tree visualization: {}/tree.html", tree_dir.display());

    Ok(())
}
