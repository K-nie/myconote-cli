//! Named presets that compose common bundles of Y1000+ subsets.
//!
//! Presets layer additively: `phylogeny` is `starter` plus phylogenetic
//! placement; `compare` is `phylogeny` plus the proteomes/orthogroups needed
//! for comparative annotation. A user can always specify `--include a,b,c`
//! directly — presets are just convenience labels.
//!
//! Citation: Opulente et al. 2024, Science 384(6694): eadj4503.

use super::subsets::Subset;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    Starter,
    Phylogeny,
    Compare,
    Reference,
    Full,
}

impl Preset {
    pub fn from_str(s: &str) -> Option<Preset> {
        match s.to_lowercase().as_str() {
            "starter"   => Some(Preset::Starter),
            "phylogeny" => Some(Preset::Phylogeny),
            "compare"   => Some(Preset::Compare),
            "reference" => Some(Preset::Reference),
            "full"      => Some(Preset::Full),
            _ => None,
        }
    }

    pub fn key(self) -> &'static str {
        match self {
            Preset::Starter   => "starter",
            Preset::Phylogeny => "phylogeny",
            Preset::Compare   => "compare",
            Preset::Reference => "reference",
            Preset::Full      => "full",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            Preset::Starter   => "KEGG + BUSCO + Codetta + tRNA + metabolism + environment (~175 MB)",
            Preset::Phylogeny => "starter + phylogeny-place (~900 MB) — enables species placement",
            Preset::Compare   => "phylogeny + proteomes + orthogroups + annotations (~4 GB)",
            Preset::Reference => "compare + genomes + cds (~12 GB) — full reference genome bundle",
            Preset::Full      => "reference + repeats + domains + phenotypes (~52 GB) — everything",
        }
    }

    /// Which subsets this preset resolves to. Order matters for --list output.
    pub fn subsets(self) -> Vec<Subset> {
        match self {
            Preset::Starter => vec![
                Subset::Kegg, Subset::Busco, Subset::Codontable,
                Subset::Trna, Subset::Metabolism, Subset::Environment,
            ],
            Preset::Phylogeny => {
                let mut s = Preset::Starter.subsets();
                s.push(Subset::PhylogenyPlace);
                s
            }
            Preset::Compare => {
                let mut s = Preset::Phylogeny.subsets();
                s.extend([Subset::Proteomes, Subset::Orthogroups, Subset::Annotations]);
                s
            }
            Preset::Reference => {
                let mut s = Preset::Compare.subsets();
                s.extend([Subset::Genomes, Subset::Cds]);
                s
            }
            Preset::Full => {
                let mut s = Preset::Reference.subsets();
                s.extend([Subset::Repeats, Subset::Domains, Subset::Phenotypes]);
                s
            }
        }
    }

    pub const ALL: &'static [Preset] = &[
        Preset::Starter, Preset::Phylogeny, Preset::Compare,
        Preset::Reference, Preset::Full,
    ];
}
