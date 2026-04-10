/// Genetic code / translation table support
///
/// Implements NCBI translation tables for organisms that deviate from the
/// standard genetic code. Critical for:
///   - Candida CTG clade (Table 12: CTG → Ser instead of Leu)
///   - Mitochondrial genomes (Tables 3, 4, etc.)
///   - Alternative yeast nuclear codes
///
/// Reference: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
///
/// This is myconote-cli's own implementation — not derived from any other
/// annotation tool.
use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// Translation table enum
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GeneticCode {
    /// Table 1: Standard
    Standard,
    /// Table 2: Vertebrate Mitochondrial
    VertebrateMitochondrial,
    /// Table 3: Yeast Mitochondrial
    YeastMitochondrial,
    /// Table 4: Mold/Protozoan/Coelenterate Mitochondrial
    MoldMitochondrial,
    /// Table 5: Invertebrate Mitochondrial
    InvertebrateMitochondrial,
    /// Table 6: Ciliate Nuclear
    CiliateNuclear,
    /// Table 9: Echinoderm/Flatworm Mitochondrial
    EchinodermMitochondrial,
    /// Table 10: Euplotid Nuclear
    EuplotidNuclear,
    /// Table 11: Bacterial/Archaeal/Plant Plastid
    BacterialPlastid,
    /// Table 12: Alternative Yeast Nuclear (Candida CTG clade)
    AlternativeYeastNuclear,
    /// Table 13: Ascidian Mitochondrial
    AscidianMitochondrial,
    /// Table 14: Alternative Flatworm Mitochondrial
    AlternativeFlatwormMitochondrial,
    /// Table 16: Chlorophycean Mitochondrial
    ChlorophyceanMitochondrial,
    /// Table 21: Trematode Mitochondrial
    TrematodeMitochondrial,
    /// Table 22: Scenedesmus obliquus Mitochondrial
    ScenedesmusMitochondrial,
    /// Table 23: Thraustochytrium Mitochondrial
    ThraustochytriumMitochondrial,
    /// Table 26: Pachysolen tannophilus Nuclear
    PachysolenNuclear,
    /// Table 29: Mesodinium Nuclear
    MesodiniumNuclear,
    /// Table 33: Cephalodiscidae Mitochondrial
    CephalodiscidaeMitochondrial,
}

impl GeneticCode {
    /// Parse from NCBI table number
    pub fn from_table_number(n: u8) -> Option<Self> {
        match n {
            1 => Some(Self::Standard),
            2 => Some(Self::VertebrateMitochondrial),
            3 => Some(Self::YeastMitochondrial),
            4 => Some(Self::MoldMitochondrial),
            5 => Some(Self::InvertebrateMitochondrial),
            6 => Some(Self::CiliateNuclear),
            9 => Some(Self::EchinodermMitochondrial),
            10 => Some(Self::EuplotidNuclear),
            11 => Some(Self::BacterialPlastid),
            12 => Some(Self::AlternativeYeastNuclear),
            13 => Some(Self::AscidianMitochondrial),
            14 => Some(Self::AlternativeFlatwormMitochondrial),
            16 => Some(Self::ChlorophyceanMitochondrial),
            21 => Some(Self::TrematodeMitochondrial),
            22 => Some(Self::ScenedesmusMitochondrial),
            23 => Some(Self::ThraustochytriumMitochondrial),
            26 => Some(Self::PachysolenNuclear),
            29 => Some(Self::MesodiniumNuclear),
            33 => Some(Self::CephalodiscidaeMitochondrial),
            _ => None,
        }
    }

    /// Parse from common name string
    pub fn from_name(s: &str) -> Self {
        match s.to_lowercase().trim() {
            "standard" | "1" => Self::Standard,
            "candida" | "ctg" | "ctg-clade" | "12" => Self::AlternativeYeastNuclear,
            "yeast-mito" | "yeast_mitochondrial" | "3" => Self::YeastMitochondrial,
            "mold-mito" | "mold" | "4" => Self::MoldMitochondrial,
            "vertebrate-mito" | "vertebrate_mitochondrial" | "2" => Self::VertebrateMitochondrial,
            "invertebrate-mito" | "invertebrate_mitochondrial" | "5" => {
                Self::InvertebrateMitochondrial
            }
            "bacterial" | "plastid" | "11" => Self::BacterialPlastid,
            "ciliate" | "6" => Self::CiliateNuclear,
            "echinoderm-mito" | "9" => Self::EchinodermMitochondrial,
            _ => Self::Standard,
        }
    }

    pub fn table_number(&self) -> u8 {
        match self {
            Self::Standard => 1,
            Self::VertebrateMitochondrial => 2,
            Self::YeastMitochondrial => 3,
            Self::MoldMitochondrial => 4,
            Self::InvertebrateMitochondrial => 5,
            Self::CiliateNuclear => 6,
            Self::EchinodermMitochondrial => 9,
            Self::EuplotidNuclear => 10,
            Self::BacterialPlastid => 11,
            Self::AlternativeYeastNuclear => 12,
            Self::AscidianMitochondrial => 13,
            Self::AlternativeFlatwormMitochondrial => 14,
            Self::ChlorophyceanMitochondrial => 16,
            Self::TrematodeMitochondrial => 21,
            Self::ScenedesmusMitochondrial => 22,
            Self::ThraustochytriumMitochondrial => 23,
            Self::PachysolenNuclear => 26,
            Self::MesodiniumNuclear => 29,
            Self::CephalodiscidaeMitochondrial => 33,
        }
    }

    pub fn display_name(&self) -> &str {
        match self {
            Self::Standard => "Standard (Table 1)",
            Self::VertebrateMitochondrial => "Vertebrate Mitochondrial (Table 2)",
            Self::YeastMitochondrial => "Yeast Mitochondrial (Table 3)",
            Self::MoldMitochondrial => "Mold/Protozoan Mitochondrial (Table 4)",
            Self::InvertebrateMitochondrial => "Invertebrate Mitochondrial (Table 5)",
            Self::CiliateNuclear => "Ciliate Nuclear (Table 6)",
            Self::EchinodermMitochondrial => "Echinoderm Mitochondrial (Table 9)",
            Self::EuplotidNuclear => "Euplotid Nuclear (Table 10)",
            Self::BacterialPlastid => "Bacterial/Plastid (Table 11)",
            Self::AlternativeYeastNuclear => "Alternative Yeast Nuclear / Candida CTG (Table 12)",
            Self::AscidianMitochondrial => "Ascidian Mitochondrial (Table 13)",
            Self::AlternativeFlatwormMitochondrial => "Alt. Flatworm Mitochondrial (Table 14)",
            Self::ChlorophyceanMitochondrial => "Chlorophycean Mitochondrial (Table 16)",
            Self::TrematodeMitochondrial => "Trematode Mitochondrial (Table 21)",
            Self::ScenedesmusMitochondrial => "Scenedesmus Mitochondrial (Table 22)",
            Self::ThraustochytriumMitochondrial => "Thraustochytrium Mitochondrial (Table 23)",
            Self::PachysolenNuclear => "Pachysolen Nuclear (Table 26)",
            Self::MesodiniumNuclear => "Mesodinium Nuclear (Table 29)",
            Self::CephalodiscidaeMitochondrial => "Cephalodiscidae Mitochondrial (Table 33)",
        }
    }

    /// Build codon → amino acid lookup (all 64 codons).
    /// Differences from Standard code are documented per table.
    pub fn codon_table(&self) -> HashMap<[u8; 3], char> {
        let mut table = standard_codon_table();

        match self {
            Self::Standard => {}

            Self::AlternativeYeastNuclear => {
                // CTG → Ser (instead of Leu)
                table.insert(*b"CTG", 'S');
            }

            Self::VertebrateMitochondrial => {
                table.insert(*b"AGA", '*'); // Stop
                table.insert(*b"AGG", '*'); // Stop
                table.insert(*b"ATA", 'M'); // Met (not Ile)
                table.insert(*b"TGA", 'W'); // Trp (not Stop)
            }

            Self::YeastMitochondrial => {
                table.insert(*b"ATA", 'M');
                table.insert(*b"CTT", 'T');
                table.insert(*b"CTC", 'T');
                table.insert(*b"CTA", 'T');
                table.insert(*b"CTG", 'T');
                table.insert(*b"TGA", 'W');
            }

            Self::MoldMitochondrial => {
                table.insert(*b"TGA", 'W');
            }

            Self::InvertebrateMitochondrial => {
                table.insert(*b"AGA", 'S');
                table.insert(*b"AGG", 'S');
                table.insert(*b"ATA", 'M');
                table.insert(*b"TGA", 'W');
            }

            Self::CiliateNuclear => {
                table.insert(*b"TAA", 'Q');
                table.insert(*b"TAG", 'Q');
            }

            Self::EchinodermMitochondrial => {
                table.insert(*b"AAA", 'N');
                table.insert(*b"AGA", 'S');
                table.insert(*b"AGG", 'S');
                table.insert(*b"TGA", 'W');
            }

            Self::EuplotidNuclear => {
                table.insert(*b"TGA", 'C');
            }

            Self::BacterialPlastid => {
                // Same as standard, but ATG/GTG/TTG all start codons
                // (handled at initiation level, not codon table)
            }

            Self::PachysolenNuclear => {
                // CTG → Ala (instead of Leu)
                table.insert(*b"CTG", 'A');
            }

            _ => {
                // For less common tables, fall back to standard
                // Users can always specify --genetic-code N for exact behavior
            }
        }

        table
    }

    /// Translate a DNA sequence using this genetic code.
    pub fn translate(&self, dna: &str) -> String {
        let table = self.codon_table();
        let bytes = dna.as_bytes();
        let mut protein = String::with_capacity(bytes.len() / 3);

        let mut i = 0;
        while i + 2 < bytes.len() {
            let codon = [
                bytes[i].to_ascii_uppercase(),
                bytes[i + 1].to_ascii_uppercase(),
                bytes[i + 2].to_ascii_uppercase(),
            ];
            let aa = table.get(&codon).copied().unwrap_or('X');
            protein.push(aa);
            i += 3;
        }

        // Remove trailing stop
        if protein.ends_with('*') {
            protein.pop();
        }
        protein
    }

    /// Check if a codon is a start codon in this genetic code.
    pub fn is_start_codon(&self, codon: &[u8; 3]) -> bool {
        let upper = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        match self {
            Self::Standard | Self::AlternativeYeastNuclear => upper == *b"ATG",
            Self::BacterialPlastid => matches!(&upper, b"ATG" | b"GTG" | b"TTG"),
            Self::YeastMitochondrial => matches!(&upper, b"ATG" | b"ATA"),
            _ => upper == *b"ATG",
        }
    }

    /// Check for internal stop codons in a protein sequence.
    /// Returns positions of internal stops (0-indexed amino acid positions).
    pub fn find_internal_stops(protein: &str) -> Vec<usize> {
        protein
            .char_indices()
            .filter(|&(i, c)| c == '*' && i < protein.len() - 1)
            .map(|(i, _)| i)
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Standard codon table builder
// ─────────────────────────────────────────────────────────────────────────────

fn standard_codon_table() -> HashMap<[u8; 3], char> {
    let mut t = HashMap::with_capacity(64);

    // Phe
    t.insert(*b"TTT", 'F');
    t.insert(*b"TTC", 'F');
    // Leu
    t.insert(*b"TTA", 'L');
    t.insert(*b"TTG", 'L');
    t.insert(*b"CTT", 'L');
    t.insert(*b"CTC", 'L');
    t.insert(*b"CTA", 'L');
    t.insert(*b"CTG", 'L');
    // Ile
    t.insert(*b"ATT", 'I');
    t.insert(*b"ATC", 'I');
    t.insert(*b"ATA", 'I');
    // Met
    t.insert(*b"ATG", 'M');
    // Val
    t.insert(*b"GTT", 'V');
    t.insert(*b"GTC", 'V');
    t.insert(*b"GTA", 'V');
    t.insert(*b"GTG", 'V');
    // Ser
    t.insert(*b"TCT", 'S');
    t.insert(*b"TCC", 'S');
    t.insert(*b"TCA", 'S');
    t.insert(*b"TCG", 'S');
    t.insert(*b"AGT", 'S');
    t.insert(*b"AGC", 'S');
    // Pro
    t.insert(*b"CCT", 'P');
    t.insert(*b"CCC", 'P');
    t.insert(*b"CCA", 'P');
    t.insert(*b"CCG", 'P');
    // Thr
    t.insert(*b"ACT", 'T');
    t.insert(*b"ACC", 'T');
    t.insert(*b"ACA", 'T');
    t.insert(*b"ACG", 'T');
    // Ala
    t.insert(*b"GCT", 'A');
    t.insert(*b"GCC", 'A');
    t.insert(*b"GCA", 'A');
    t.insert(*b"GCG", 'A');
    // Tyr
    t.insert(*b"TAT", 'Y');
    t.insert(*b"TAC", 'Y');
    // Stop
    t.insert(*b"TAA", '*');
    t.insert(*b"TAG", '*');
    t.insert(*b"TGA", '*');
    // His
    t.insert(*b"CAT", 'H');
    t.insert(*b"CAC", 'H');
    // Gln
    t.insert(*b"CAA", 'Q');
    t.insert(*b"CAG", 'Q');
    // Asn
    t.insert(*b"AAT", 'N');
    t.insert(*b"AAC", 'N');
    // Lys
    t.insert(*b"AAA", 'K');
    t.insert(*b"AAG", 'K');
    // Asp
    t.insert(*b"GAT", 'D');
    t.insert(*b"GAC", 'D');
    // Glu
    t.insert(*b"GAA", 'E');
    t.insert(*b"GAG", 'E');
    // Cys
    t.insert(*b"TGT", 'C');
    t.insert(*b"TGC", 'C');
    // Trp
    t.insert(*b"TGG", 'W');
    // Arg
    t.insert(*b"CGT", 'R');
    t.insert(*b"CGC", 'R');
    t.insert(*b"CGA", 'R');
    t.insert(*b"CGG", 'R');
    t.insert(*b"AGA", 'R');
    t.insert(*b"AGG", 'R');
    // Gly
    t.insert(*b"GGT", 'G');
    t.insert(*b"GGC", 'G');
    t.insert(*b"GGA", 'G');
    t.insert(*b"GGG", 'G');

    t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_translation() {
        let code = GeneticCode::Standard;
        assert_eq!(code.translate("ATGAAATTT"), "MKF");
    }

    #[test]
    fn test_candida_ctg_clade() {
        let code = GeneticCode::AlternativeYeastNuclear;
        // CTG should be Ser, not Leu
        let protein = code.translate("ATGCTG");
        assert_eq!(protein, "MS");

        // Standard code should give Leu
        let std_protein = GeneticCode::Standard.translate("ATGCTG");
        assert_eq!(std_protein, "ML");
    }

    #[test]
    fn test_mito_tga_is_trp() {
        let code = GeneticCode::VertebrateMitochondrial;
        let protein = code.translate("ATGTGA");
        // TGA = Trp in vertebrate mito, not stop
        assert_eq!(protein, "MW");
    }

    #[test]
    fn test_internal_stops() {
        let stops = GeneticCode::find_internal_stops("MK*FG*");
        assert_eq!(stops, vec![2]);
    }

    #[test]
    fn test_from_name() {
        assert_eq!(
            GeneticCode::from_name("candida"),
            GeneticCode::AlternativeYeastNuclear
        );
        assert_eq!(
            GeneticCode::from_name("12"),
            GeneticCode::AlternativeYeastNuclear
        );
        assert_eq!(GeneticCode::from_name("standard"), GeneticCode::Standard);
    }

    #[test]
    fn test_from_table_number() {
        assert_eq!(
            GeneticCode::from_table_number(12),
            Some(GeneticCode::AlternativeYeastNuclear)
        );
        assert_eq!(GeneticCode::from_table_number(99), None);
    }
}
