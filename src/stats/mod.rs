use crate::parser::gff::GFFReader;
use crate::parser::region::RegionSelector;
use crate::utils::error::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

pub mod taxon;

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct GenomeStatistics {
    pub total_features: usize,
    pub total_genes: usize,
    pub total_transcripts: usize,
    pub total_cds: usize,
    pub total_exons: usize,
    pub gene_lengths: Vec<u64>,
    pub chromosome_stats: HashMap<String, ChromosomeStats>,

    // Isoform tracking
    pub isoforms_per_gene: HashMap<String, usize>,
    pub primary_only_genes: Option<PrimaryGeneStats>,
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct ChromosomeStats {
    pub gene_count: usize,
    pub transcript_count: usize,
    pub cds_count: usize,
    pub exon_count: usize,
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct PrimaryGeneStats {
    pub total_genes: usize,
    pub total_transcripts: usize,
    pub total_cds: usize,
    pub total_exons: usize,
    pub gene_lengths: Vec<u64>,
    pub chromosome_stats: HashMap<String, ChromosomeStats>,
}

impl PrimaryGeneStats {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn mean_gene_length(&self) -> f64 {
        if self.gene_lengths.is_empty() {
            0.0
        } else {
            let sum: u64 = self.gene_lengths.iter().sum();
            sum as f64 / self.gene_lengths.len() as f64
        }
    }
}

impl GenomeStatistics {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_gff<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_gff_with_selector(path, &RegionSelector::new(), false)
    }

    pub fn from_gff_primary<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_gff_with_selector(path, &RegionSelector::new(), true)
    }

    pub fn from_gff_with_selector<P: AsRef<Path>>(
        path: P,
        selector: &RegionSelector,
        primary_only: bool,
    ) -> Result<Self> {
        let mut stats = GenomeStatistics::new();
        let reader = GFFReader::from_path(path)?;

        // Track relationships for isoform analysis
        let mut gene_to_transcripts: HashMap<String, Vec<String>> = HashMap::new();
        let mut transcript_to_cds: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut transcript_to_exons: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut transcript_lengths: HashMap<String, u64> = HashMap::new();
        // Per-gene chromosome (seqid). Needed by the primary_only path so it
        // can credit each collapsed isoform to its actual chromosome instead
        // of the arbitrary first chromosome in the HashMap iteration order.
        let mut gene_to_seqid: HashMap<String, String> = HashMap::new();

        let mut record_count = 0;

        // First pass: collect all features
        for result in reader {
            match result {
                Ok(record) => {
                    record_count += 1;

                    if !selector.should_include_feature(&record.seqid, record.start, record.end) {
                        continue;
                    }

                    stats.total_features += 1;

                    match record.feature_type.as_str() {
                        "gene" => {
                            stats.total_genes += 1;
                            stats.gene_lengths.push(record.length());

                            let chr_stats = stats
                                .chromosome_stats
                                .entry(record.seqid.clone())
                                .or_insert(ChromosomeStats::default());
                            chr_stats.gene_count += 1;

                            if let Some(id) = record.id() {
                                gene_to_transcripts.insert(id.clone(), Vec::new());
                                gene_to_seqid.insert(id.clone(), record.seqid.clone());
                            }
                        }
                        "mRNA" | "transcript" => {
                            stats.total_transcripts += 1;

                            if let (Some(parent), Some(id)) = (record.parent(), record.id()) {
                                if let Some(transcripts) = gene_to_transcripts.get_mut(parent) {
                                    transcripts.push(id.clone());
                                }
                                transcript_lengths.insert(id.clone(), record.length());
                            }

                            let chr_stats = stats
                                .chromosome_stats
                                .entry(record.seqid.clone())
                                .or_insert(ChromosomeStats::default());
                            chr_stats.transcript_count += 1;
                        }
                        "CDS" => {
                            stats.total_cds += 1;

                            if let Some(parent) = record.parent() {
                                let cds_coords = (record.start, record.end);
                                transcript_to_cds
                                    .entry(parent.clone())
                                    .or_insert(Vec::new())
                                    .push(cds_coords);
                            }

                            let chr_stats = stats
                                .chromosome_stats
                                .entry(record.seqid.clone())
                                .or_insert(ChromosomeStats::default());
                            chr_stats.cds_count += 1;
                        }
                        "exon" => {
                            stats.total_exons += 1;

                            if let Some(parent) = record.parent() {
                                let exon_coords = (record.start, record.end);
                                transcript_to_exons
                                    .entry(parent.clone())
                                    .or_insert(Vec::new())
                                    .push(exon_coords);
                            }

                            let chr_stats = stats
                                .chromosome_stats
                                .entry(record.seqid.clone())
                                .or_insert(ChromosomeStats::default());
                            chr_stats.exon_count += 1;
                        }
                        _ => {}
                    }
                }
                Err(e) => {
                    return Err(e);
                }
            }
        }

        let _ = record_count; // used only for counting; suppress unused warning

        // Calculate isoforms per gene
        for (gene_id, transcripts) in &gene_to_transcripts {
            stats
                .isoforms_per_gene
                .insert(gene_id.clone(), transcripts.len());
        }

        // If primary_only is true, build a separate PrimaryGeneStats object
        if primary_only {
            let mut primary_stats = PrimaryGeneStats::new();
            let mut primary_chromosome_stats: HashMap<String, ChromosomeStats> = HashMap::new();

            // For each gene, find the longest transcript
            for (gene_id, transcripts) in &gene_to_transcripts {
                if transcripts.is_empty() {
                    continue;
                }

                // Find the longest transcript for this gene
                let primary = transcripts
                    .iter()
                    .max_by_key(|t| transcript_lengths.get(*t).unwrap_or(&0))
                    .unwrap();

                primary_stats.total_genes += 1;
                primary_stats.total_transcripts += 1;

                if let Some(len) = transcript_lengths.get(primary) {
                    primary_stats.gene_lengths.push(*len);
                }

                // Add CDS counts for this transcript
                if let Some(cds_list) = transcript_to_cds.get(primary) {
                    primary_stats.total_cds += cds_list.len();
                }

                // Add exon counts for this transcript
                if let Some(exon_list) = transcript_to_exons.get(primary) {
                    primary_stats.total_exons += exon_list.len();
                }

                // Credit this gene to its actual chromosome (seqid captured
                // during the first pass) — the prior implementation did
                // `for (chr, _) in &stats.chromosome_stats { ...; break; }`
                // which silently credited every gene to the first chromosome
                // in HashMap iteration order, producing nonsense per-chrom
                // stats on any multi-chromosome input.
                if let Some(seqid) = gene_to_seqid.get(gene_id) {
                    let chr_stat = primary_chromosome_stats
                        .entry(seqid.clone())
                        .or_insert(ChromosomeStats::default());
                    chr_stat.gene_count += 1;
                    chr_stat.transcript_count += 1;
                    if let Some(cds_list) = transcript_to_cds.get(primary) {
                        chr_stat.cds_count += cds_list.len();
                    }
                    if let Some(exon_list) = transcript_to_exons.get(primary) {
                        chr_stat.exon_count += exon_list.len();
                    }
                }
            }

            primary_stats.chromosome_stats = primary_chromosome_stats;
            stats.primary_only_genes = Some(primary_stats);
        }

        Ok(stats)
    }

    pub fn mean_gene_length(&self) -> f64 {
        if self.gene_lengths.is_empty() {
            0.0
        } else {
            let sum: u64 = self.gene_lengths.iter().sum();
            sum as f64 / self.gene_lengths.len() as f64
        }
    }

    pub fn median_gene_length(&self) -> f64 {
        if self.gene_lengths.is_empty() {
            return 0.0;
        }
        let mut sorted = self.gene_lengths.clone();
        sorted.sort();
        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
        } else {
            sorted[mid] as f64
        }
    }

    pub fn min_gene_length(&self) -> u64 {
        if self.gene_lengths.is_empty() {
            0
        } else {
            *self.gene_lengths.iter().min().unwrap()
        }
    }

    pub fn max_gene_length(&self) -> u64 {
        if self.gene_lengths.is_empty() {
            0
        } else {
            *self.gene_lengths.iter().max().unwrap()
        }
    }

    pub fn n50(&self) -> u64 {
        if self.gene_lengths.is_empty() {
            return 0;
        }
        let mut lengths = self.gene_lengths.clone();
        lengths.sort_by(|a, b| b.cmp(a));
        let total: u64 = lengths.iter().sum();
        let half = total / 2;
        let mut cumulative = 0;
        for &len in &lengths {
            cumulative += len;
            if cumulative >= half {
                return len;
            }
        }
        0
    }

    pub fn get_taxon_warning(&self, taxon: &str) -> Option<String> {
        taxon::get_taxon_warning(self, taxon)
    }

    pub fn print_summary(&self) {
        println!("\n============================================================");
        println!("MYCONOTE GENOME STATISTICS SUMMARY");
        println!("============================================================");
        println!("\nFeature Counts:");
        println!("  Total features: {}", self.total_features);
        println!("  Genes: {}", self.total_genes);
        println!("  Transcripts: {}", self.total_transcripts);
        println!("  CDS: {}", self.total_cds);
        println!("  Exons: {}", self.total_exons);

        println!("\nIsoform Statistics:");
        let total_isoforms: usize = self.isoforms_per_gene.values().sum();
        if self.total_genes > 0 {
            let avg_isoforms = total_isoforms as f64 / self.total_genes as f64;
            println!("  Average isoforms per gene: {:.2}", avg_isoforms);
            println!("  Total isoforms: {}", total_isoforms);

            // Show distribution
            let mut single_copy = 0;
            let mut multi_copy = 0;
            for &count in self.isoforms_per_gene.values() {
                if count == 1 {
                    single_copy += 1;
                } else {
                    multi_copy += 1;
                }
            }
            println!(
                "  Single-isoform genes: {} ({:.1}%)",
                single_copy,
                (single_copy as f64 / self.total_genes as f64) * 100.0
            );
            println!(
                "  Multi-isoform genes: {} ({:.1}%)",
                multi_copy,
                (multi_copy as f64 / self.total_genes as f64) * 100.0
            );
        }

        // Show primary-only stats if available
        if let Some(primary) = &self.primary_only_genes {
            println!("\n[PRIMARY TRANSCRIPT ONLY]");
            println!("  Genes: {}", primary.total_genes);
            println!("  Transcripts: {}", primary.total_transcripts);
            println!("  CDS: {}", primary.total_cds);
            println!("  Exons: {}", primary.total_exons);
            println!("  Mean gene length: {:.2} bp", primary.mean_gene_length());

            println!("\n  Primary Chromosome Statistics:");
            for (chr, stats) in &primary.chromosome_stats {
                println!(
                    "    {}: {} genes, {} transcripts, {} CDS, {} exons",
                    chr,
                    stats.gene_count,
                    stats.transcript_count,
                    stats.cds_count,
                    stats.exon_count
                );
            }
        }

        println!("\nGene Length Statistics (bp):");
        println!("  Mean:   {:.2}", self.mean_gene_length());
        println!("  Median: {:.2}", self.median_gene_length());
        println!("  Min:    {}", self.min_gene_length());
        println!("  Max:    {}", self.max_gene_length());
        println!("  N50:    {}", self.n50());

        println!("\nChromosome Statistics:");
        for (chr, stats) in &self.chromosome_stats {
            println!(
                "  {}: {} genes, {} transcripts, {} CDS, {} exons",
                chr, stats.gene_count, stats.transcript_count, stats.cds_count, stats.exon_count
            );
        }
    }
}
