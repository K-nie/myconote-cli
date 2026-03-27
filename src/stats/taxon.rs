use serde::{Serialize, Deserialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonReference {
    pub name: String,
    pub expected_gene_count: (usize, usize), // (min, max)
    pub expected_gene_length: (f64, f64),    // (min, max) in bp
    pub expected_exons_per_gene: (f64, f64), // (min, max)
    pub description: String,
}

impl TaxonReference {
    pub fn new(name: &str, 
               gene_min: usize, gene_max: usize,
               length_min: f64, length_max: f64,
               exons_min: f64, exons_max: f64,
               description: &str) -> Self {
        Self {
            name: name.to_string(),
            expected_gene_count: (gene_min, gene_max),
            expected_gene_length: (length_min, length_max),
            expected_exons_per_gene: (exons_min, exons_max),
            description: description.to_string(),
        }
    }
}

lazy_static::lazy_static! {
    pub static ref TAXON_DATABASE: HashMap<String, TaxonReference> = {
        let mut m = HashMap::new();
        
        // Fungi
        m.insert("fungi".to_string(), TaxonReference::new(
            "Fungi",
            5000, 15000,      // gene count range
            100.0, 2000.0,    // gene length range (bp)
            2.0, 5.0,         // exons per gene range
            "Fungal genomes typically have compact gene structures with few introns"
        ));
        
        // Ascomycota (a fungal subdivision)
        m.insert("ascomycota".to_string(), TaxonReference::new(
            "Ascomycota",
            6000, 12000,
            300.0, 1500.0,
            2.0, 4.0,
            "Ascomycete fungi like yeast and molds"
        ));
        
        // Basidiomycota
        m.insert("basidiomycota".to_string(), TaxonReference::new(
            "Basidiomycota",
            8000, 20000,
            400.0, 2000.0,
            3.0, 6.0,
            "Basidiomycete fungi like mushrooms and rusts"
        ));
        
        // Plants
        m.insert("plants".to_string(), TaxonReference::new(
            "Plants",
            20000, 50000,
            1000.0, 5000.0,
            4.0, 8.0,
            "Plant genomes often have larger genes with more exons"
        ));
        
        // Animals
        m.insert("animals".to_string(), TaxonReference::new(
            "Animals",
            15000, 25000,
            1000.0, 10000.0,
            5.0, 12.0,
            "Animal genomes have variable gene structures"
        ));
        
        // Mammals
        m.insert("mammals".to_string(), TaxonReference::new(
            "Mammals",
            19000, 23000,
            1500.0, 15000.0,
            6.0, 15.0,
            "Mammalian genes often have many exons due to alternative splicing"
        ));
        
        m
    };
}

pub fn get_taxon_warning(stats: &super::GenomeStatistics, taxon: &str) -> Option<String> {
    let db = TAXON_DATABASE.get(taxon)?;
    
    let mut warnings = Vec::new();
    
    // Check gene count
    if stats.total_genes < db.expected_gene_count.0 {
        warnings.push(format!(
            "Gene count ({}) is LOWER than expected for {} ({} - {})",
            stats.total_genes, db.name, db.expected_gene_count.0, db.expected_gene_count.1
        ));
    } else if stats.total_genes > db.expected_gene_count.1 {
        warnings.push(format!(
            "Gene count ({}) is HIGHER than expected for {} ({} - {})",
            stats.total_genes, db.name, db.expected_gene_count.0, db.expected_gene_count.1
        ));
    }
    
    // Check gene length
    let mean_len = stats.mean_gene_length();
    if mean_len < db.expected_gene_length.0 {
        warnings.push(format!(
            "Mean gene length ({:.1} bp) is SHORTER than expected for {} ({:.0} - {:.0} bp)",
            mean_len, db.name, db.expected_gene_length.0, db.expected_gene_length.1
        ));
    } else if mean_len > db.expected_gene_length.1 {
        warnings.push(format!(
            "Mean gene length ({:.1} bp) is LONGER than expected for {} ({:.0} - {:.0} bp)",
            mean_len, db.name, db.expected_gene_length.0, db.expected_gene_length.1
        ));
    }
    
    // Check exons per gene (if we have exon data)
    if stats.total_genes > 0 && stats.total_exons > 0 {
        let mean_exons = stats.total_exons as f64 / stats.total_genes as f64;
        if mean_exons < db.expected_exons_per_gene.0 {
            warnings.push(format!(
                "Mean exons per gene ({:.2}) is LOWER than expected for {} ({:.1} - {:.1})",
                mean_exons, db.name, db.expected_exons_per_gene.0, db.expected_exons_per_gene.1
            ));
        } else if mean_exons > db.expected_exons_per_gene.1 {
            warnings.push(format!(
                "Mean exons per gene ({:.2}) is HIGHER than expected for {} ({:.1} - {:.1})",
                mean_exons, db.name, db.expected_exons_per_gene.0, db.expected_exons_per_gene.1
            ));
        }
    }
    
    if warnings.is_empty() {
        Some(format!("✓ Genome statistics are within expected range for {}", db.name))
    } else {
        Some(format!(
            "⚠️  Taxonomic Benchmark Warning for {}:\n  {}\n  {}",
            db.name,
            warnings.join("\n  "),
            db.description
        ))
    }
}
