//! Static registry of every Y1000+ subset myconote can install.
//!
//! Source: Opulente et al. 2024 "Genomic factors shape carbon and nitrogen
//! metabolic niche breadth across Saccharomycotina yeasts." Science 384(6694).
//! Collection: <https://plus.figshare.com/collections/_/6714042>.
//!
//! Each entry points to its canonical figshare download and declares how to
//! normalise it into the on-disk layout that downstream subcommands expect.

use std::fmt;

/// Every downloadable Y1000+ subset. Keep the variant names in sync with
/// the user-facing `--include` strings (lowercase, no underscores).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Subset {
    Kegg,
    Busco,
    Codontable,
    Trna,
    Metabolism,
    Environment,
    Phenotypes,
    PhylogenyPlace,
    Annotations,
    Orthogroups,
    Proteomes,
    Cds,
    Genomes,
    Repeats,
    Domains,
}

impl Subset {
    /// All subsets in a stable order — drives `--list` output.
    pub const ALL: &'static [Subset] = &[
        Subset::Kegg,
        Subset::Busco,
        Subset::Codontable,
        Subset::Trna,
        Subset::Metabolism,
        Subset::Environment,
        Subset::Phenotypes,
        Subset::PhylogenyPlace,
        Subset::Annotations,
        Subset::Orthogroups,
        Subset::Proteomes,
        Subset::Cds,
        Subset::Genomes,
        Subset::Repeats,
        Subset::Domains,
    ];

    pub fn from_str(s: &str) -> Option<Subset> {
        match s.to_lowercase().as_str() {
            "kegg" => Some(Subset::Kegg),
            "busco" => Some(Subset::Busco),
            "codontable" => Some(Subset::Codontable),
            "trna" => Some(Subset::Trna),
            "metabolism" => Some(Subset::Metabolism),
            "environment" => Some(Subset::Environment),
            "phenotypes" => Some(Subset::Phenotypes),
            "phylogeny-place" | "phylogeny_place" | "phyloplace" => Some(Subset::PhylogenyPlace),
            "annotations" => Some(Subset::Annotations),
            "orthogroups" => Some(Subset::Orthogroups),
            "proteomes" => Some(Subset::Proteomes),
            "cds" => Some(Subset::Cds),
            "genomes" => Some(Subset::Genomes),
            "repeats" => Some(Subset::Repeats),
            "domains" => Some(Subset::Domains),
            _ => None,
        }
    }

    pub fn key(self) -> &'static str {
        match self {
            Subset::Kegg => "kegg",
            Subset::Busco => "busco",
            Subset::Codontable => "codontable",
            Subset::Trna => "trna",
            Subset::Metabolism => "metabolism",
            Subset::Environment => "environment",
            Subset::Phenotypes => "phenotypes",
            Subset::PhylogenyPlace => "phylogeny-place",
            Subset::Annotations => "annotations",
            Subset::Orthogroups => "orthogroups",
            Subset::Proteomes => "proteomes",
            Subset::Cds => "cds",
            Subset::Genomes => "genomes",
            Subset::Repeats => "repeats",
            Subset::Domains => "domains",
        }
    }

    /// One-line description shown by `setup --y1000plus --list`.
    pub fn summary(self) -> &'static str {
        match self {
            Subset::Kegg => "KEGG pathway transfer for annotate",
            Subset::Busco => "BUSCO percentile benchmarking for stats",
            Subset::Codontable => "Auto-detect CTG-Ser clade codon tables (predict)",
            Subset::Trna => "tRNA-count benchmarking vs 1,154 yeasts",
            Subset::Metabolism => "Carbon/nitrogen specialism prediction",
            Subset::Environment => "Isolation-source / ecological niche prediction",
            Subset::Phenotypes => "Growth-rate & phenotype prediction",
            Subset::PhylogenyPlace => "Place a new genome in the 1,154-yeast tree (EPA-ng)",
            Subset::Annotations => "Reference GFF3/GTF for comparative analyses",
            Subset::Orthogroups => "Pre-computed OrthoFinder orthogroups (skip OrthoFinder)",
            Subset::Proteomes => "Yeast-specific gene naming via DIAMOND",
            Subset::Cds => "CDS FASTAs for codon-usage / alignment work",
            Subset::Genomes => "Raw reference genome FASTAs",
            Subset::Repeats => "RepeatMasker-annotated repeats for mask lift-over",
            Subset::Domains => "InterProScan domain annotations (full bundle)",
        }
    }

    /// Which figshare file(s) back this subset. Most subsets map to a single
    /// archive; a few (phylogeny-place) bundle several small supporting files.
    pub fn files(self) -> &'static [FileSpec] {
        match self {
            Subset::Kegg => &[FileSpec {
                name: "y1000plus_annotations_pep_kegg.tar.gz",
                url: "https://ndownloader.figshare.com/files/40479305",
                size_bytes: 20_500_000, // ~20.5 MB
                extract: ExtractKind::TarGz,
                dest_subdir: "kegg",
            }],
            Subset::Busco => &[FileSpec {
                name: "Y1000p_BUCO_fulltable.tar.gz",
                url: "https://ndownloader.figshare.com/files/40535813",
                size_bytes: 78_900_000,
                extract: ExtractKind::TarGz,
                dest_subdir: "busco",
            }],
            Subset::Codontable => &[FileSpec {
                name: "y1000p_codetta_output.tar.gz",
                url: "https://ndownloader.figshare.com/files/40479242",
                size_bytes: 8_200_000,
                extract: ExtractKind::TarGz,
                dest_subdir: "codetta",
            }],
            Subset::Trna => &[FileSpec {
                name: "y1000p_tRNA_scan.tar.gz",
                url: "https://ndownloader.figshare.com/files/40479287",
                size_bytes: 45_500_000,
                extract: ExtractKind::TarGz,
                dest_subdir: "trna",
            }],
            Subset::Metabolism => &[
                // Species × KEGG-KO presence/absence matrix — the table that
                // powers metabolic niche-breadth reasoning in Opulente 2024.
                FileSpec {
                    name: "Y1000_KEGG_Annotations.xlsx",
                    url: "https://ndownloader.figshare.com/files/40809164",
                    size_bytes: 16_144_339,
                    extract: ExtractKind::Xlsx,
                    dest_subdir: "metabolism",
                },
                // Per-species generalist/specialist carbon & nitrogen labels.
                FileSpec {
                    name: "CarbonNitrogen_Classifications.zip",
                    url: "https://ndownloader.figshare.com/files/41164340",
                    size_bytes: 188_308,
                    extract: ExtractKind::Zip,
                    dest_subdir: "metabolism",
                },
            ],
            Subset::Environment => &[
                // OWL ontology of isolation environments (RDF/XML).
                FileSpec {
                    name: "budding_yeasts_isolations.owl",
                    url: "https://ndownloader.figshare.com/files/40637525",
                    size_bytes: 1_377_557,
                    extract: ExtractKind::Passthrough,
                    dest_subdir: "environment",
                },
                // Species-to-species distances from phylogenetic PCA of growth.
                FileSpec {
                    name: "PCA_distance_matrix.csv",
                    url: "https://ndownloader.figshare.com/files/40638095",
                    size_bytes: 515_597,
                    extract: ExtractKind::Passthrough,
                    dest_subdir: "environment",
                },
            ],
            Subset::Phenotypes => &[
                // Normalised per-strain growth rates across 24 conditions —
                // small (180 KB) but analytically richest for phenotype prediction.
                FileSpec {
                    name: "GrowthRates_Yeasts.xlsx",
                    url: "https://ndownloader.figshare.com/files/40705427",
                    size_bytes: 179_594,
                    extract: ExtractKind::Xlsx,
                    dest_subdir: "phenotypes",
                },
                // Growth at 37 °C (Y/N/W/V/S) from Figure 6 — tiny, drives
                // thermotolerance prediction.
                FileSpec {
                    name: "y1000p_growth_at_37.xlsx",
                    url: "https://ndownloader.figshare.com/files/40638881",
                    size_bytes: 51_409,
                    extract: ExtractKind::Xlsx,
                    dest_subdir: "phenotypes",
                },
            ],
            Subset::PhylogenyPlace => &[
                FileSpec {
                    name: "y1000p_1403_OG_files.tar.gz",
                    url: "https://ndownloader.figshare.com/files/40540178",
                    size_bytes: 724_686_853,
                    extract: ExtractKind::TarGz,
                    dest_subdir: "phylogeny/marker_ogs",
                },
                FileSpec {
                    name: "1175taxa_1403OGs_iqtree_ML_rooted.tre",
                    url: "https://ndownloader.figshare.com/files/43128517",
                    size_bytes: 63_791,
                    extract: ExtractKind::Passthrough,
                    dest_subdir: "phylogeny",
                },
                FileSpec {
                    name: "1175taxa_1403OGs_astral_length_rooted_lpp_added.tre",
                    url: "https://ndownloader.figshare.com/files/43128514",
                    size_bytes: 61_589,
                    extract: ExtractKind::Passthrough,
                    dest_subdir: "phylogeny",
                },
                FileSpec {
                    name: "1154yeasts_1403OGs_ml_timetree.tree",
                    url: "https://ndownloader.figshare.com/files/43128520",
                    size_bytes: 53_647,
                    extract: ExtractKind::Passthrough,
                    dest_subdir: "phylogeny",
                },
            ],
            Subset::Annotations => &[
                FileSpec {
                    name: "y1000p_gff3_files.tar.gz",
                    url: "https://ndownloader.figshare.com/files/40534991",
                    size_bytes: 376_532_012,
                    extract: ExtractKind::TarGz,
                    dest_subdir: "gff3",
                },
                FileSpec {
                    name: "y1000p_gtf_files.tar.gz",
                    url: "https://ndownloader.figshare.com/files/40534994",
                    size_bytes: 178_575_096,
                    extract: ExtractKind::TarGz,
                    dest_subdir: "gtf",
                },
            ],
            Subset::Orthogroups => &[FileSpec {
                name: "y1000p_orthofinder.tar.gz",
                url: "https://ndownloader.figshare.com/files/40479275",
                size_bytes: 1_300_000_000,
                extract: ExtractKind::TarGz,
                dest_subdir: "orthogroups",
            }],
            Subset::Proteomes => &[FileSpec {
                name: "y1000p_pep_files.tar.gz",
                url: "https://ndownloader.figshare.com/files/40534997",
                size_bytes: 1_953_997_436,
                extract: ExtractKind::TarGzThenDiamondIndex,
                dest_subdir: "pep",
            }],
            Subset::Cds => &[FileSpec {
                name: "y1000p_cds_files.tar.gz",
                url: "https://ndownloader.figshare.com/files/45439234",
                size_bytes: 2_938_548_138,
                extract: ExtractKind::TarGz,
                dest_subdir: "cds",
            }],
            Subset::Genomes => &[FileSpec {
                name: "y1000p_genome_files.zip",
                url: "https://ndownloader.figshare.com/files/40535165",
                size_bytes: 4_731_016_231,
                extract: ExtractKind::Zip,
                dest_subdir: "genomes",
            }],
            Subset::Repeats => &[FileSpec {
                name: "y1000p_repeatmasker_analysis.tar.gz",
                url: "https://ndownloader.figshare.com/files/40511831",
                size_bytes: 5_700_000_000,
                extract: ExtractKind::TarGz,
                dest_subdir: "repeatmasker",
            }],
            Subset::Domains => &[FileSpec {
                name: "annotations_pep_interproscan.zip",
                url: "https://ndownloader.figshare.com/files/43229655",
                size_bytes: 34_300_000_000,
                extract: ExtractKind::Zip,
                dest_subdir: "interproscan",
            }],
        }
    }

    /// Total bytes for this subset (sum across all its files).
    pub fn total_bytes(self) -> u64 {
        self.files().iter().map(|f| f.size_bytes).sum()
    }
}

impl fmt::Display for Subset {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.key())
    }
}

/// One figshare file that needs to be fetched + extracted.
#[derive(Debug, Clone)]
pub struct FileSpec {
    pub name: &'static str,
    pub url: &'static str,
    pub size_bytes: u64,
    pub extract: ExtractKind,
    /// Subpath inside `~/.myconote/y1000plus/` where extracted contents land.
    pub dest_subdir: &'static str,
}

/// How to unpack a downloaded file into the cache.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtractKind {
    /// Just copy the file into dest_subdir (for Newick, CSV, single FASTA).
    Passthrough,
    /// Standard `.tar.gz` → extract in place.
    TarGz,
    /// `.tar.gz` + post-process: concat all per-species FASTAs and build a
    /// DIAMOND index at `<root>/diamond_index.dmnd`.
    TarGzThenDiamondIndex,
    /// `.zip` → extract in place.
    Zip,
    /// `.xlsx` → parse via calamine into a `.tsv` inside dest_subdir.
    Xlsx,
}

/// Format a byte count for human display: 1.95 GB, 78.9 MB, 20.5 kB.
pub fn format_bytes(n: u64) -> String {
    const GB: u64 = 1_000_000_000;
    const MB: u64 = 1_000_000;
    const KB: u64 = 1_000;
    if n >= GB {
        format!("{:.2} GB", n as f64 / GB as f64)
    } else if n >= MB {
        format!("{:.1} MB", n as f64 / MB as f64)
    } else if n >= KB {
        format!("{:.1} kB", n as f64 / KB as f64)
    } else {
        format!("{} B", n)
    }
}
