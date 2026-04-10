/// Kingdom-aware parameter selection
///
/// Maps `--kingdom fungi|plant|animal|insect|protist` to:
///   - Recommended Augustus training species
///   - Expected gene density (genes / Mbp)
///   - Expected intron size range (min, max bp)
///   - BUSCO lineage dataset name
///   - Typical genome size range (Mbp) for sanity checks

#[derive(Debug, Clone, PartialEq)]
pub enum Kingdom {
    Fungi,
    Plant,
    Animal,
    Insect,
    Protist,
}

impl Kingdom {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().trim() {
            "fungi" | "fungal" | "fungus" => Kingdom::Fungi,
            "plant" | "plants" | "viridiplantae" => Kingdom::Plant,
            "animal" | "animals" | "metazoa" => Kingdom::Animal,
            "insect" | "insecta" => Kingdom::Insect,
            "protist" | "protista" => Kingdom::Protist,
            _ => Kingdom::Fungi, // safe default
        }
    }

    pub fn display_name(&self) -> &str {
        match self {
            Kingdom::Fungi => "Fungi",
            Kingdom::Plant => "Viridiplantae",
            Kingdom::Animal => "Metazoa",
            Kingdom::Insect => "Insecta",
            Kingdom::Protist => "Protista",
        }
    }

    /// Recommended Augustus species name, optionally overridden by `--species`
    pub fn default_augustus_species(&self) -> &str {
        match self {
            Kingdom::Fungi => "saccharomyces_cerevisiae_S288C",
            Kingdom::Plant => "arabidopsis",
            Kingdom::Animal => "human",
            Kingdom::Insect => "fly",
            Kingdom::Protist => "toxoplasma",
        }
    }

    /// Broader list of well-trained Augustus species for this kingdom.
    /// Users can pick one with `--species`.
    pub fn augustus_species_list(&self) -> &[&str] {
        match self {
            Kingdom::Fungi => &[
                "aspergillus_fumigatus",
                "aspergillus_nidulans",
                "aspergillus_oryzae",
                "botrytis_cinerea",
                "candida_albicans",
                "candida_guilliermondii",
                "candida_tropicalis",
                "chaetomium_globosum",
                "cryptococcus_neoformans_neoformans_JEC21",
                "fusarium_graminearum",
                "histoplasma_capsulatum",
                "magnaporthe_grisea",
                "neurospora_crassa",
                "phanerochaete_chrysosporium",
                "pichia_stipitis",
                "rhizopus_oryzae",
                "saccharomyces",
                "ustilago_maydis",
                "yarrowia_lipolytica",
            ],
            Kingdom::Plant => &[
                "arabidopsis",
                "rice",
                "maize",
                "tomato",
                "wheat",
                "sorghum",
                "populus",
                "vitis",
                "medicago",
            ],
            Kingdom::Animal => &[
                "human",
                "mouse",
                "chicken",
                "zebrafish",
                "c_elegans",
                "xenopus",
                "sea_anemon",
            ],
            Kingdom::Insect => &[
                "fly",
                "honeybee1",
                "tribolium2012",
                "nasonia",
                "pea_aphid",
                "aedes",
                "anopheles",
            ],
            Kingdom::Protist => &["toxoplasma", "leishmania", "tetrahymena"],
        }
    }

    /// Expected intron size range (min, max bp)
    pub fn intron_size_range(&self) -> (u64, u64) {
        match self {
            Kingdom::Fungi => (40, 2_000),
            Kingdom::Plant => (40, 50_000),
            Kingdom::Animal => (40, 500_000),
            Kingdom::Insect => (40, 50_000),
            Kingdom::Protist => (20, 1_000),
        }
    }

    /// BUSCO OrthoDB lineage dataset name
    pub fn busco_lineage(&self) -> &str {
        match self {
            Kingdom::Fungi => "fungi_odb10",
            Kingdom::Plant => "viridiplantae_odb10",
            Kingdom::Animal => "metazoa_odb10",
            Kingdom::Insect => "insecta_odb10",
            Kingdom::Protist => "eukaryota_odb10",
        }
    }

    /// Augustus UTR model flag.
    /// Disabled by default — not all species models ship with trained UTR
    /// parameters (e.g. saccharomyces_cerevisiae_S288C lacks _utr_probs.pbl).
    /// UTR prediction is enabled automatically only when a user-trained model
    /// that includes UTR parameters is detected, or via an explicit CLI flag.
    pub fn augustus_utr(&self) -> bool {
        false
    }
}
