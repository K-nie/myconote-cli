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

// ─────────────────────────────────────────────────────────────────────────────
// Clade-aware Augustus species selection (fungi)
// ─────────────────────────────────────────────────────────────────────────────
//
// `Kingdom::Fungi` previously defaulted *every* fungus to
// `saccharomyces_cerevisiae_S288C`.  S288C is a Saccharomycotina
// (budding-yeast) model with compact, intron-poor genes; on divergent fungi
// (Basidiomycota, Pezizomycotina) that gene model is badly wrong and
// prediction accuracy collapses silently — the benchmark harness recorded
// BUSCO completeness of 14.9 % and a near-zero strict gene-F1 on a
// Cryptococcus (Basidiomycota) genome run under S288C.
//
// The fix is a small curated registry mapping each pre-trained Augustus
// fungal model (the `Kingdom::Fungi` species list above) to its major clade,
// a per-clade default model, and a loud guard when a naive run would fall
// back to S288C off-Saccharomycotina.

/// Major fungal clades MycoNote recognises for Augustus model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FungalClade {
    /// Budding yeasts (e.g. *Saccharomyces*, *Candida*, *Pichia*, *Yarrowia*).
    Saccharomycotina,
    /// Filamentous ascomycetes (e.g. *Aspergillus*, *Neurospora*, *Fusarium*).
    Pezizomycotina,
    /// Basidiomycetes (e.g. *Cryptococcus*, *Ustilago*, *Phanerochaete*).
    Basidiomycota,
    /// Everything else (early-diverging lineages such as Mucoromycota). No
    /// strong pre-trained default — triggers the self-train recommendation.
    Other,
}

impl FungalClade {
    /// Parse a `--clade` value. Accepts the clade name, a short alias, and a
    /// couple of plain-language synonyms. Returns `None` for unknown input so
    /// the CLI can report a clear error listing the accepted values.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().trim() {
            "saccharomycotina" | "saccharomycetes" | "sacch" | "yeast" | "budding_yeast"
            | "budding-yeast" => Some(FungalClade::Saccharomycotina),
            "pezizomycotina"
            | "pezizo"
            | "filamentous"
            | "filamentous_ascomycete"
            | "euascomycete" => Some(FungalClade::Pezizomycotina),
            "basidiomycota" | "basidiomycete" | "basidio" => Some(FungalClade::Basidiomycota),
            "other" | "mucoromycota" | "zygomycete" | "early_diverging" => Some(FungalClade::Other),
            _ => None,
        }
    }

    /// Canonical name for display / logging.
    pub fn display_name(&self) -> &'static str {
        match self {
            FungalClade::Saccharomycotina => "Saccharomycotina",
            FungalClade::Pezizomycotina => "Pezizomycotina",
            FungalClade::Basidiomycota => "Basidiomycota",
            FungalClade::Other => "other fungi",
        }
    }

    /// The pre-trained Augustus model to use by default for this clade.
    /// `None` for `Other` — there is no single good stock model for
    /// early-diverging fungi, so a self-trained model is recommended instead.
    pub fn default_augustus_species(&self) -> Option<&'static str> {
        match self {
            // S288C is the right call *inside its own clade*.
            FungalClade::Saccharomycotina => Some("saccharomyces_cerevisiae_S288C"),
            // A. nidulans is the canonical well-trained filamentous ascomycete.
            FungalClade::Pezizomycotina => Some("aspergillus_nidulans"),
            // C. neoformans is the stock basidiomycete model.
            FungalClade::Basidiomycota => Some("cryptococcus_neoformans_neoformans_JEC21"),
            FungalClade::Other => None,
        }
    }
}

/// Curated registry: each pre-trained Augustus fungal model → its clade.
///
/// The species list mirrors `Kingdom::Fungi.augustus_species_list()` plus the
/// `saccharomyces_cerevisiae_S288C` default. It is the set of models that ship
/// with a standard Augustus `config/species/` install; MycoNote does not
/// enumerate the on-disk directory (it may be absent on this machine), so the
/// registry is defined from the standard Augustus fungal species and is the
/// documented source of truth for `--clade` defaults and species→clade lookup.
pub const FUNGAL_SPECIES_CLADES: &[(&str, FungalClade)] = &[
    // ── Saccharomycotina (budding yeasts) ────────────────────────────────────
    (
        "saccharomyces_cerevisiae_S288C",
        FungalClade::Saccharomycotina,
    ),
    ("saccharomyces", FungalClade::Saccharomycotina),
    ("candida_albicans", FungalClade::Saccharomycotina),
    ("candida_guilliermondii", FungalClade::Saccharomycotina),
    ("candida_tropicalis", FungalClade::Saccharomycotina),
    ("pichia_stipitis", FungalClade::Saccharomycotina),
    ("yarrowia_lipolytica", FungalClade::Saccharomycotina),
    // ── Pezizomycotina (filamentous ascomycetes) ─────────────────────────────
    ("aspergillus_fumigatus", FungalClade::Pezizomycotina),
    ("aspergillus_nidulans", FungalClade::Pezizomycotina),
    ("aspergillus_oryzae", FungalClade::Pezizomycotina),
    ("botrytis_cinerea", FungalClade::Pezizomycotina),
    ("chaetomium_globosum", FungalClade::Pezizomycotina),
    ("fusarium_graminearum", FungalClade::Pezizomycotina),
    ("histoplasma_capsulatum", FungalClade::Pezizomycotina),
    ("magnaporthe_grisea", FungalClade::Pezizomycotina),
    ("neurospora_crassa", FungalClade::Pezizomycotina),
    // ── Basidiomycota ────────────────────────────────────────────────────────
    (
        "cryptococcus_neoformans_neoformans_JEC21",
        FungalClade::Basidiomycota,
    ),
    ("phanerochaete_chrysosporium", FungalClade::Basidiomycota),
    ("ustilago_maydis", FungalClade::Basidiomycota),
    // ── Other (early-diverging) ──────────────────────────────────────────────
    ("rhizopus_oryzae", FungalClade::Other),
];

/// Look up the clade of a known Augustus fungal model by exact name.
pub fn clade_for_species(species: &str) -> Option<FungalClade> {
    FUNGAL_SPECIES_CLADES
        .iter()
        .find(|(name, _)| *name == species)
        .map(|(_, clade)| *clade)
}

/// Resolve the Augustus species name for a prediction run, applying
/// clade-aware defaults and returning an optional guard warning.
///
/// Precedence (highest first):
///   1. `explicit_species` — an explicit `--species <name>` (unchanged, no
///      warning). Callers resolve a self-trained model *before* this.
///   2. `clade` default — when `--clade` is given and the clade has a stock
///      model (Saccharomycotina / Pezizomycotina / Basidiomycota).
///   3. Fallback — for Fungi with no usable hint, S288C **plus** a prominent
///      warning recommending `--clade` / `--species` / self-training. For
///      non-fungal kingdoms, the kingdom default (unchanged, no warning).
///
/// Returns `(species, Some(warning))` when the fungal fallback is taken,
/// otherwise `(species, None)`. The function never hard-fails: a naive
/// `predict genome.fa` still runs, it is just no longer silent.
pub fn resolve_augustus_species(
    kingdom: &Kingdom,
    explicit_species: Option<&str>,
    clade: Option<FungalClade>,
) -> (String, Option<String>) {
    // 1. Explicit --species always wins and is left exactly as-is.
    if let Some(s) = explicit_species {
        return (s.to_string(), None);
    }

    // Clade-aware logic only applies to Fungi; other kingdoms are unchanged.
    if matches!(kingdom, Kingdom::Fungi) {
        // 2. Clade default (when the clade has a stock model).
        if let Some(c) = clade {
            if let Some(sp) = c.default_augustus_species() {
                return (sp.to_string(), None);
            }
            // `--clade other`: no stock model, fall through to the warning.
        }

        // 3. Fallback: S288C + loud guard.
        let hint = match clade {
            Some(FungalClade::Other) => "No stock Augustus model exists for early-diverging fungi",
            _ => "No --species or --clade given",
        };
        let warning = format!(
            "  ⚠  {hint}. Defaulting Augustus to 'saccharomyces_cerevisiae_S288C',\n\
             \x20     a Saccharomycotina (budding-yeast) model that is a POOR default for\n\
             \x20     other fungi — on divergent fungi (Basidiomycota, Pezizomycotina)\n\
             \x20     it silently collapses gene-prediction accuracy (BUSCO completeness\n\
             \x20     can fall below 15%). If this genome is not a budding yeast, pick one:\n\
             \x20       --clade <saccharomycotina|pezizomycotina|basidiomycota>\n\
             \x20       --species <augustus_model>      (see 'myconote-cli species')\n\
             \x20       --augustus-training self         (train a genome-specific model)"
        );
        return ("saccharomyces_cerevisiae_S288C".to_string(), Some(warning));
    }

    // Non-fungal kingdoms: unchanged behaviour, no warning.
    (kingdom.default_augustus_species().to_string(), None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_clade_has_expected_default() {
        assert_eq!(
            FungalClade::Saccharomycotina.default_augustus_species(),
            Some("saccharomyces_cerevisiae_S288C")
        );
        assert_eq!(
            FungalClade::Pezizomycotina.default_augustus_species(),
            Some("aspergillus_nidulans")
        );
        assert_eq!(
            FungalClade::Basidiomycota.default_augustus_species(),
            Some("cryptococcus_neoformans_neoformans_JEC21")
        );
        assert_eq!(FungalClade::Other.default_augustus_species(), None);
    }

    #[test]
    fn clade_from_str_accepts_names_and_aliases() {
        assert_eq!(
            FungalClade::from_str("saccharomycotina"),
            Some(FungalClade::Saccharomycotina)
        );
        assert_eq!(
            FungalClade::from_str("YEAST"),
            Some(FungalClade::Saccharomycotina)
        );
        assert_eq!(
            FungalClade::from_str("basidio"),
            Some(FungalClade::Basidiomycota)
        );
        assert_eq!(
            FungalClade::from_str(" Pezizomycotina "),
            Some(FungalClade::Pezizomycotina)
        );
        assert_eq!(FungalClade::from_str("other"), Some(FungalClade::Other));
        assert_eq!(FungalClade::from_str("bogus"), None);
    }

    #[test]
    fn registry_resolves_representative_species_to_right_clade() {
        assert_eq!(
            clade_for_species("saccharomyces_cerevisiae_S288C"),
            Some(FungalClade::Saccharomycotina)
        );
        assert_eq!(
            clade_for_species("aspergillus_nidulans"),
            Some(FungalClade::Pezizomycotina)
        );
        assert_eq!(
            clade_for_species("cryptococcus_neoformans_neoformans_JEC21"),
            Some(FungalClade::Basidiomycota)
        );
        assert_eq!(
            clade_for_species("ustilago_maydis"),
            Some(FungalClade::Basidiomycota)
        );
        assert_eq!(clade_for_species("not_a_real_model"), None);
    }

    #[test]
    fn registry_default_is_in_registry_and_maps_back_to_its_clade() {
        for clade in [
            FungalClade::Saccharomycotina,
            FungalClade::Pezizomycotina,
            FungalClade::Basidiomycota,
        ] {
            let def = clade
                .default_augustus_species()
                .expect("clade has a default");
            assert_eq!(
                clade_for_species(def),
                Some(clade),
                "default model {def} should be registered under its own clade"
            );
        }
    }

    #[test]
    fn explicit_species_overrides_clade_and_warns_never() {
        // --species wins even when a (different) --clade is also given, and no
        // warning is emitted — this is the benchmark's path, must be unchanged.
        let (sp, warn) = resolve_augustus_species(
            &Kingdom::Fungi,
            Some("fusarium_graminearum"),
            Some(FungalClade::Basidiomycota),
        );
        assert_eq!(sp, "fusarium_graminearum");
        assert!(warn.is_none());
    }

    #[test]
    fn clade_basidiomycota_selects_basidiomycete_not_s288c() {
        let (sp, warn) =
            resolve_augustus_species(&Kingdom::Fungi, None, Some(FungalClade::Basidiomycota));
        assert_eq!(sp, "cryptococcus_neoformans_neoformans_JEC21");
        assert_ne!(sp, "saccharomyces_cerevisiae_S288C");
        assert!(warn.is_none());
    }

    #[test]
    fn no_hint_fungus_still_returns_s288c_but_warns() {
        // Back-compat: a naive `predict genome.fa` still gets S288C (so results
        // are unchanged for anyone relying on it), but it is no longer silent.
        let (sp, warn) = resolve_augustus_species(&Kingdom::Fungi, None, None);
        assert_eq!(sp, "saccharomyces_cerevisiae_S288C");
        let w = warn.expect("naive fungal run must emit a guard warning");
        assert!(w.contains("saccharomyces_cerevisiae_S288C"));
        assert!(w.contains("--clade"));
    }

    #[test]
    fn clade_other_falls_back_to_s288c_with_warning() {
        let (sp, warn) = resolve_augustus_species(&Kingdom::Fungi, None, Some(FungalClade::Other));
        assert_eq!(sp, "saccharomyces_cerevisiae_S288C");
        assert!(warn.is_some());
    }

    #[test]
    fn non_fungal_kingdoms_are_unchanged_and_silent() {
        // Clade is fungal-only; a non-fungal kingdom ignores it and keeps its
        // own default with no warning.
        let (sp, warn) =
            resolve_augustus_species(&Kingdom::Plant, None, Some(FungalClade::Basidiomycota));
        assert_eq!(sp, "arabidopsis");
        assert!(warn.is_none());

        let (sp, warn) = resolve_augustus_species(&Kingdom::Animal, None, None);
        assert_eq!(sp, "human");
        assert!(warn.is_none());
    }
}
