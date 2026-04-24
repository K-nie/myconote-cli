//! Genetic code / translation table support.
//!
//! Implements the 25 NCBI translation tables currently defined at
//! <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>:
//! tables 1–6, 9–14, 16, 21–31, 33.
//!
//! The source of truth for every table is the `REGISTRY` static below.
//! Adding a newly-published NCBI table is a single registry entry:
//! table number, display name, codon deviations from Standard, and
//! (optional) alternative start codons. No match-arm maintenance across
//! the file.
//!
//! A `GeneticCode` value is a validated NCBI table number. Code outside
//! this module should construct one via `GeneticCode::from_table_number`
//! (returns `Option`) or use `GeneticCode::STANDARD` when defaulting.

use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// Public type: a validated NCBI table number
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GeneticCode {
    table: u8,
}

impl GeneticCode {
    /// NCBI table 1 (Standard). The default for any organism without a
    /// declared alternative code.
    pub const STANDARD: Self = Self { table: 1 };

    /// Validate `n` against the registry; return `Some(code)` if the
    /// number names a supported table, `None` otherwise.
    pub fn from_table_number(n: u8) -> Option<Self> {
        if REGISTRY.iter().any(|d| d.table == n) {
            Some(Self { table: n })
        } else {
            None
        }
    }

    /// Case-insensitive alias lookup. Accepts common spellings of the
    /// well-known tables (`candida`, `ctg`, `yeast-mito`, …) plus the
    /// literal table number as a string. Falls back to the Standard
    /// code on unknown input; callers that want strict validation
    /// should use `from_table_number`.
    pub fn from_name(s: &str) -> Self {
        let s = s.trim().to_ascii_lowercase();
        // Numeric alias path — try the number first, covers "1", "12", etc.
        if let Ok(n) = s.parse::<u8>() {
            if let Some(code) = Self::from_table_number(n) {
                return code;
            }
        }
        let table = match s.as_str() {
            "standard" => 1,
            "vertebrate-mito" | "vertebrate_mitochondrial" => 2,
            "yeast-mito" | "yeast_mitochondrial" => 3,
            "mold-mito" | "mold_mitochondrial" | "mold" | "protozoan" | "mycoplasma"
            | "spiroplasma" => 4,
            "invertebrate-mito" | "invertebrate_mitochondrial" => 5,
            "ciliate" | "ciliate_nuclear" | "dasycladacean" | "hexamita" => 6,
            "echinoderm-mito" | "flatworm-mito" => 9,
            "euplotid" | "euplotid_nuclear" => 10,
            "bacterial" | "archaeal" | "plastid" => 11,
            "candida" | "ctg" | "ctg-clade" | "alt-yeast" | "alternative_yeast_nuclear" => 12,
            "ascidian-mito" => 13,
            "alt-flatworm-mito" => 14,
            "chlorophycean-mito" => 16,
            "trematode-mito" => 21,
            "scenedesmus-mito" => 22,
            "thraustochytrium-mito" => 23,
            "rhabdopleuridae-mito" => 24,
            "sr1" | "gracilibacteria" => 25,
            "pachysolen" | "pachysolen_nuclear" => 26,
            "karyorelict" | "karyorelict_nuclear" => 27,
            "condylostoma" | "condylostoma_nuclear" => 28,
            "mesodinium" | "mesodinium_nuclear" => 29,
            "peritrich" | "peritrich_nuclear" => 30,
            "blastocrithidia" | "blastocrithidia_nuclear" => 31,
            "cephalodiscidae-mito" => 33,
            _ => 1,
        };
        Self { table }
    }

    /// The raw NCBI table number (for emitting `transl_table=N` in
    /// GFF3 / GenBank).
    pub fn table_number(&self) -> u8 {
        self.table
    }

    /// Human-readable name including the table number for log lines
    /// and reports.
    pub fn display_name(&self) -> String {
        let def = self.def();
        format!("{} (Table {})", def.name, def.table)
    }

    /// Non-`None` when this table has context-dependent codon behavior
    /// that we approximate rather than implement precisely (e.g.
    /// tables 27, 28, 31 where TAR / TGA can be either an amino acid
    /// or a stop depending on local context).
    pub fn caveat(&self) -> Option<&'static str> {
        self.def().caveat
    }

    /// Build the full 64-codon lookup for this table. Uses the
    /// Standard table as the baseline and applies this table's
    /// registered deviations.
    pub fn codon_table(&self) -> HashMap<[u8; 3], char> {
        let mut table = standard_codon_table();
        for (codon, aa) in self.def().deviations {
            table.insert(*codon, *aa);
        }
        table
    }

    /// Translate a DNA sequence using this genetic code. Trailing
    /// stops are trimmed; internal stops are preserved so callers can
    /// detect them via `find_internal_stops`.
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
        if protein.ends_with('*') {
            protein.pop();
        }
        protein
    }

    /// Is `codon` a valid initiation codon for this table?
    /// Every table permits ATG; some also permit GTG, TTG, ATA, etc.
    /// The extended list comes from `alt_start_codons` in the registry.
    pub fn is_start_codon(&self, codon: &[u8; 3]) -> bool {
        let upper = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        if upper == *b"ATG" {
            return true;
        }
        self.def().alt_start_codons.iter().any(|&alt| alt == upper)
    }

    /// Positions (0-indexed amino acid) of internal stops in a
    /// translated protein. Handy for flagging mistranslated CDS that
    /// slipped past the gene predictor.
    pub fn find_internal_stops(protein: &str) -> Vec<usize> {
        protein
            .char_indices()
            .filter(|&(i, c)| c == '*' && i < protein.len() - 1)
            .map(|(i, _)| i)
            .collect()
    }

    /// All supported tables in registry order. Used by CLI listing
    /// commands and by tests that want to sweep every table.
    pub fn all() -> impl Iterator<Item = GeneticCode> {
        REGISTRY.iter().map(|d| GeneticCode { table: d.table })
    }

    fn def(&self) -> &'static TableDef {
        REGISTRY
            .iter()
            .find(|d| d.table == self.table)
            .unwrap_or(&REGISTRY[0])
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Registry (the single source of truth)
// ─────────────────────────────────────────────────────────────────────────────

/// One NCBI translation table, as a pure-data record.
///
/// `deviations` lists every codon that differs from the Standard table
/// (Table 1). Each entry overrides the Standard assignment.
///
/// `alt_start_codons` lists additional permissible start codons
/// beyond ATG. ATG itself is always a start codon and is omitted here
/// to avoid duplication.
///
/// `caveat`, when present, documents a known limitation — typically
/// context-dependent codons (tables 27, 28, 31) that we resolve to a
/// single amino-acid interpretation rather than modelling the actual
/// context switch.
struct TableDef {
    table: u8,
    name: &'static str,
    deviations: &'static [([u8; 3], char)],
    alt_start_codons: &'static [[u8; 3]],
    caveat: Option<&'static str>,
}

/// The 25 NCBI translation tables currently defined (missing numbers
/// 7, 8, 15, 17–20, 32 are retired / withdrawn by NCBI).
///
/// Each entry is taken directly from
/// <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
/// (accessed 2026-04-24).
static REGISTRY: &[TableDef] = &[
    TableDef {
        table: 1,
        name: "Standard",
        deviations: &[],
        alt_start_codons: &[],
        caveat: None,
    },
    TableDef {
        table: 2,
        name: "Vertebrate Mitochondrial",
        deviations: &[
            (*b"AGA", '*'), // Arg → Stop
            (*b"AGG", '*'), // Arg → Stop
            (*b"ATA", 'M'), // Ile → Met
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[*b"ATT", *b"ATC", *b"ATA", *b"GTG"],
        caveat: None,
    },
    TableDef {
        table: 3,
        name: "Yeast Mitochondrial",
        deviations: &[
            (*b"ATA", 'M'), // Ile → Met
            (*b"CTT", 'T'), // Leu → Thr
            (*b"CTC", 'T'), // Leu → Thr
            (*b"CTA", 'T'), // Leu → Thr
            (*b"CTG", 'T'), // Leu → Thr
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[*b"ATA"],
        caveat: None,
    },
    TableDef {
        table: 4,
        name: "Mold / Protozoan / Coelenterate Mitochondrial; Mycoplasma / Spiroplasma",
        deviations: &[
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[
            *b"TTA", *b"TTG", *b"CTG", *b"ATT", *b"ATC", *b"ATA", *b"GTG",
        ],
        caveat: None,
    },
    TableDef {
        table: 5,
        name: "Invertebrate Mitochondrial",
        deviations: &[
            (*b"AGA", 'S'), // Arg → Ser
            (*b"AGG", 'S'), // Arg → Ser
            (*b"ATA", 'M'), // Ile → Met
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[*b"ATT", *b"ATC", *b"ATA", *b"GTG"],
        caveat: None,
    },
    TableDef {
        table: 6,
        name: "Ciliate / Dasycladacean / Hexamita Nuclear",
        deviations: &[
            (*b"TAA", 'Q'), // Stop → Gln
            (*b"TAG", 'Q'), // Stop → Gln
        ],
        alt_start_codons: &[],
        caveat: None,
    },
    TableDef {
        table: 9,
        name: "Echinoderm / Flatworm Mitochondrial",
        deviations: &[
            (*b"AAA", 'N'), // Lys → Asn
            (*b"AGA", 'S'), // Arg → Ser
            (*b"AGG", 'S'), // Arg → Ser
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[*b"GTG"],
        caveat: None,
    },
    TableDef {
        table: 10,
        name: "Euplotid Nuclear",
        deviations: &[
            (*b"TGA", 'C'), // Stop → Cys
        ],
        alt_start_codons: &[],
        caveat: None,
    },
    TableDef {
        table: 11,
        name: "Bacterial / Archaeal / Plant Plastid",
        deviations: &[],
        alt_start_codons: &[*b"GTG", *b"TTG", *b"ATT", *b"CTG"],
        caveat: None,
    },
    TableDef {
        table: 12,
        name: "Alternative Yeast Nuclear (Candida CTG clade)",
        deviations: &[
            (*b"CTG", 'S'), // Leu → Ser
        ],
        alt_start_codons: &[*b"CTG"],
        caveat: None,
    },
    TableDef {
        table: 13,
        name: "Ascidian Mitochondrial",
        deviations: &[
            (*b"AGA", 'G'), // Arg → Gly
            (*b"AGG", 'G'), // Arg → Gly
            (*b"ATA", 'M'), // Ile → Met
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[*b"ATA", *b"GTG", *b"TTG"],
        caveat: None,
    },
    TableDef {
        table: 14,
        name: "Alternative Flatworm Mitochondrial",
        deviations: &[
            (*b"AAA", 'N'), // Lys → Asn
            (*b"AGA", 'S'), // Arg → Ser
            (*b"AGG", 'S'), // Arg → Ser
            (*b"TAA", 'Y'), // Stop → Tyr
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[],
        caveat: None,
    },
    TableDef {
        table: 16,
        name: "Chlorophycean Mitochondrial",
        deviations: &[
            (*b"TAG", 'L'), // Stop → Leu
        ],
        alt_start_codons: &[],
        caveat: None,
    },
    TableDef {
        table: 21,
        name: "Trematode Mitochondrial",
        deviations: &[
            (*b"TGA", 'W'), // Stop → Trp
            (*b"ATA", 'M'), // Ile → Met
            (*b"AGA", 'S'), // Arg → Ser
            (*b"AGG", 'S'), // Arg → Ser
            (*b"AAA", 'N'), // Lys → Asn
        ],
        alt_start_codons: &[*b"GTG"],
        caveat: None,
    },
    TableDef {
        table: 22,
        name: "Scenedesmus obliquus Mitochondrial",
        deviations: &[
            (*b"TCA", '*'), // Ser → Stop
            (*b"TAG", 'L'), // Stop → Leu
        ],
        alt_start_codons: &[],
        caveat: None,
    },
    TableDef {
        table: 23,
        name: "Thraustochytrium Mitochondrial",
        deviations: &[
            (*b"TTA", '*'), // Leu → Stop
        ],
        alt_start_codons: &[*b"ATT", *b"GTG"],
        caveat: None,
    },
    TableDef {
        table: 24,
        name: "Rhabdopleuridae Mitochondrial",
        deviations: &[
            (*b"AGA", 'S'), // Arg → Ser
            (*b"AGG", 'K'), // Arg → Lys
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[*b"GTG", *b"TTG", *b"CTG"],
        caveat: None,
    },
    TableDef {
        table: 25,
        name: "Candidate Division SR1 / Gracilibacteria",
        deviations: &[
            (*b"TGA", 'G'), // Stop → Gly
        ],
        alt_start_codons: &[*b"TTG", *b"GTG"],
        caveat: None,
    },
    TableDef {
        table: 26,
        name: "Pachysolen tannophilus Nuclear",
        deviations: &[
            (*b"CTG", 'A'), // Leu → Ala
        ],
        alt_start_codons: &[*b"CTG"],
        caveat: None,
    },
    TableDef {
        table: 27,
        name: "Karyorelict Nuclear",
        deviations: &[
            (*b"TAA", 'Q'), // Stop → Gln (context-dependent)
            (*b"TAG", 'Q'), // Stop → Gln (context-dependent)
            (*b"TGA", 'W'), // Stop → Trp (context-dependent)
        ],
        alt_start_codons: &[],
        caveat: Some(
            "Table 27 has context-dependent stop codons: TAA/TAG/TGA may be \
             read-through amino acids OR terminators depending on local context. \
             MycoNote-CLI resolves them to their amino-acid reading (Q/Q/W); \
             proper context disambiguation is beyond this translator.",
        ),
    },
    TableDef {
        table: 28,
        name: "Condylostoma Nuclear",
        deviations: &[
            (*b"TAA", 'Q'), // Stop → Gln (context-dependent)
            (*b"TAG", 'Q'), // Stop → Gln (context-dependent)
            (*b"TGA", 'W'), // Stop → Trp (context-dependent)
        ],
        alt_start_codons: &[],
        caveat: Some("Table 28 has context-dependent stop codons — same caveat as Table 27."),
    },
    TableDef {
        table: 29,
        name: "Mesodinium Nuclear",
        deviations: &[
            (*b"TAA", 'Y'), // Stop → Tyr
            (*b"TAG", 'Y'), // Stop → Tyr
        ],
        alt_start_codons: &[],
        caveat: None,
    },
    TableDef {
        table: 30,
        name: "Peritrich Nuclear",
        deviations: &[
            (*b"TAA", 'E'), // Stop → Glu
            (*b"TAG", 'E'), // Stop → Glu
        ],
        alt_start_codons: &[],
        caveat: None,
    },
    TableDef {
        table: 31,
        name: "Blastocrithidia Nuclear",
        deviations: &[
            (*b"TAA", 'E'), // Stop → Glu (context-dependent)
            (*b"TAG", 'E'), // Stop → Glu (context-dependent)
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[],
        caveat: Some(
            "Table 31 has context-dependent TAA/TAG: may be read-through Glu or \
             terminators depending on local context. Resolved here as Glu.",
        ),
    },
    TableDef {
        table: 33,
        name: "Cephalodiscidae Mitochondrial",
        deviations: &[
            (*b"AGA", 'S'), // Arg → Ser
            (*b"AGG", 'K'), // Arg → Lys
            (*b"TAA", 'Y'), // Stop → Tyr
            (*b"TGA", 'W'), // Stop → Trp
        ],
        alt_start_codons: &[*b"GTG", *b"TTG", *b"CTG"],
        caveat: None,
    },
];

// ─────────────────────────────────────────────────────────────────────────────
// Standard codon table (Table 1) — the baseline everyone else overrides
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

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: translate a single codon and return the amino acid.
    fn tx1(code: GeneticCode, codon: &[u8; 3]) -> char {
        code.codon_table().get(codon).copied().unwrap_or('X')
    }

    // ── Registry integrity ───────────────────────────────────────────────────

    #[test]
    fn registry_has_25_tables() {
        assert_eq!(REGISTRY.len(), 25);
    }

    #[test]
    fn registry_covers_expected_ncbi_numbers() {
        let expected: [u8; 25] = [
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
            31, 33,
        ];
        let actual: Vec<u8> = REGISTRY.iter().map(|d| d.table).collect();
        assert_eq!(actual, expected.to_vec());
    }

    #[test]
    fn registry_numbers_are_unique_and_sorted() {
        let mut seen = std::collections::HashSet::new();
        let mut prev: i32 = -1;
        for d in REGISTRY {
            assert!(seen.insert(d.table), "duplicate table number {}", d.table);
            assert!(
                (d.table as i32) > prev,
                "registry not sorted at {}",
                d.table
            );
            prev = d.table as i32;
        }
    }

    // ── Lookup helpers ───────────────────────────────────────────────────────

    #[test]
    fn from_table_number_validates() {
        for n in [1u8, 2, 3, 12, 25, 33] {
            assert!(GeneticCode::from_table_number(n).is_some());
        }
        // 7, 8, 15, 17–20, 32 are retired by NCBI
        for n in [0u8, 7, 8, 15, 17, 18, 19, 20, 32, 99, 255] {
            assert!(
                GeneticCode::from_table_number(n).is_none(),
                "table {n} should not be supported"
            );
        }
    }

    #[test]
    fn from_name_covers_common_aliases() {
        assert_eq!(
            GeneticCode::from_name("candida").table_number(),
            12,
            "candida alias should map to table 12"
        );
        assert_eq!(GeneticCode::from_name("ctg").table_number(), 12);
        assert_eq!(GeneticCode::from_name("12").table_number(), 12);
        assert_eq!(GeneticCode::from_name("standard").table_number(), 1);
        assert_eq!(GeneticCode::from_name("yeast-mito").table_number(), 3);
        assert_eq!(GeneticCode::from_name("gracilibacteria").table_number(), 25);
        assert_eq!(
            GeneticCode::from_name("bogus-label").table_number(),
            1,
            "unknown aliases fall back to Standard (by design)"
        );
    }

    // ── Per-table codon deviations (sampled against NCBI reference) ──────────

    #[test]
    fn table_1_standard_has_no_deviations() {
        let code = GeneticCode::STANDARD;
        assert_eq!(tx1(code, b"CTG"), 'L');
        assert_eq!(tx1(code, b"ATA"), 'I');
        assert_eq!(tx1(code, b"TGA"), '*');
        assert_eq!(tx1(code, b"TAA"), '*');
    }

    #[test]
    fn table_2_vertebrate_mito() {
        let code = GeneticCode::from_table_number(2).unwrap();
        assert_eq!(tx1(code, b"AGA"), '*');
        assert_eq!(tx1(code, b"AGG"), '*');
        assert_eq!(tx1(code, b"ATA"), 'M');
        assert_eq!(tx1(code, b"TGA"), 'W');
    }

    #[test]
    fn table_3_yeast_mito() {
        let code = GeneticCode::from_table_number(3).unwrap();
        for ctn in [b"CTT", b"CTC", b"CTA", b"CTG"] {
            assert_eq!(tx1(code, ctn), 'T', "table 3: {:?} should be Thr", ctn);
        }
        assert_eq!(tx1(code, b"ATA"), 'M');
        assert_eq!(tx1(code, b"TGA"), 'W');
    }

    #[test]
    fn table_4_mold_mito() {
        let code = GeneticCode::from_table_number(4).unwrap();
        assert_eq!(tx1(code, b"TGA"), 'W');
        // No other deviations
        assert_eq!(tx1(code, b"CTG"), 'L');
    }

    #[test]
    fn table_5_invertebrate_mito() {
        let code = GeneticCode::from_table_number(5).unwrap();
        assert_eq!(tx1(code, b"AGA"), 'S');
        assert_eq!(tx1(code, b"AGG"), 'S');
        assert_eq!(tx1(code, b"ATA"), 'M');
        assert_eq!(tx1(code, b"TGA"), 'W');
    }

    #[test]
    fn table_6_ciliate() {
        let code = GeneticCode::from_table_number(6).unwrap();
        assert_eq!(tx1(code, b"TAA"), 'Q');
        assert_eq!(tx1(code, b"TAG"), 'Q');
        // TGA still stop
        assert_eq!(tx1(code, b"TGA"), '*');
    }

    #[test]
    fn table_9_echinoderm_mito() {
        let code = GeneticCode::from_table_number(9).unwrap();
        assert_eq!(tx1(code, b"AAA"), 'N');
        assert_eq!(tx1(code, b"AGA"), 'S');
        assert_eq!(tx1(code, b"AGG"), 'S');
        assert_eq!(tx1(code, b"TGA"), 'W');
    }

    #[test]
    fn table_10_euplotid() {
        let code = GeneticCode::from_table_number(10).unwrap();
        assert_eq!(tx1(code, b"TGA"), 'C');
        assert_eq!(tx1(code, b"TAA"), '*');
        assert_eq!(tx1(code, b"TAG"), '*');
    }

    #[test]
    fn table_11_bacterial_plastid_same_assignments_as_standard() {
        let code = GeneticCode::from_table_number(11).unwrap();
        // Bacterial table is IDENTICAL to Standard in codon assignments;
        // it only differs in start-codon tolerance (ATG/GTG/TTG/ATT/CTG).
        assert_eq!(tx1(code, b"CTG"), 'L');
        assert_eq!(tx1(code, b"TGA"), '*');
        assert!(code.is_start_codon(b"GTG"));
        assert!(code.is_start_codon(b"TTG"));
        assert!(code.is_start_codon(b"ATT"));
        assert!(code.is_start_codon(b"CTG"));
    }

    #[test]
    fn table_12_candida_ctg_clade() {
        let code = GeneticCode::from_table_number(12).unwrap();
        assert_eq!(tx1(code, b"CTG"), 'S');
        // Standard code should still give Leu — sanity of baseline.
        assert_eq!(tx1(GeneticCode::STANDARD, b"CTG"), 'L');
    }

    #[test]
    fn table_13_ascidian_mito() {
        let code = GeneticCode::from_table_number(13).unwrap();
        // WAS SILENTLY STANDARD in the pre-registry implementation.
        assert_eq!(tx1(code, b"AGA"), 'G');
        assert_eq!(tx1(code, b"AGG"), 'G');
        assert_eq!(tx1(code, b"ATA"), 'M');
        assert_eq!(tx1(code, b"TGA"), 'W');
    }

    #[test]
    fn table_14_alt_flatworm_mito() {
        let code = GeneticCode::from_table_number(14).unwrap();
        assert_eq!(tx1(code, b"AAA"), 'N');
        assert_eq!(tx1(code, b"AGA"), 'S');
        assert_eq!(tx1(code, b"AGG"), 'S');
        assert_eq!(tx1(code, b"TAA"), 'Y');
        assert_eq!(tx1(code, b"TGA"), 'W');
        // TAG still stop
        assert_eq!(tx1(code, b"TAG"), '*');
    }

    #[test]
    fn table_16_chlorophycean_mito() {
        let code = GeneticCode::from_table_number(16).unwrap();
        assert_eq!(tx1(code, b"TAG"), 'L');
        // TAA, TGA still stop
        assert_eq!(tx1(code, b"TAA"), '*');
        assert_eq!(tx1(code, b"TGA"), '*');
    }

    #[test]
    fn table_21_trematode_mito() {
        let code = GeneticCode::from_table_number(21).unwrap();
        assert_eq!(tx1(code, b"TGA"), 'W');
        assert_eq!(tx1(code, b"ATA"), 'M');
        assert_eq!(tx1(code, b"AGA"), 'S');
        assert_eq!(tx1(code, b"AGG"), 'S');
        assert_eq!(tx1(code, b"AAA"), 'N');
    }

    #[test]
    fn table_22_scenedesmus_mito() {
        let code = GeneticCode::from_table_number(22).unwrap();
        assert_eq!(tx1(code, b"TCA"), '*'); // Ser → Stop
        assert_eq!(tx1(code, b"TAG"), 'L'); // Stop → Leu
                                            // Other Ser codons unchanged
        assert_eq!(tx1(code, b"TCT"), 'S');
    }

    #[test]
    fn table_23_thraustochytrium_mito() {
        let code = GeneticCode::from_table_number(23).unwrap();
        assert_eq!(tx1(code, b"TTA"), '*'); // Leu → Stop
        assert_eq!(tx1(code, b"TTG"), 'L'); // other Leu unchanged
    }

    #[test]
    fn table_24_rhabdopleuridae_mito() {
        let code = GeneticCode::from_table_number(24).unwrap();
        assert_eq!(tx1(code, b"AGA"), 'S');
        assert_eq!(tx1(code, b"AGG"), 'K');
        assert_eq!(tx1(code, b"TGA"), 'W');
    }

    #[test]
    fn table_25_gracilibacteria() {
        let code = GeneticCode::from_table_number(25).unwrap();
        assert_eq!(tx1(code, b"TGA"), 'G');
        assert_eq!(tx1(code, b"TAA"), '*');
        assert_eq!(tx1(code, b"TAG"), '*');
    }

    #[test]
    fn table_26_pachysolen_nuclear() {
        let code = GeneticCode::from_table_number(26).unwrap();
        assert_eq!(tx1(code, b"CTG"), 'A'); // Leu → Ala
    }

    #[test]
    fn table_27_karyorelict_has_caveat() {
        let code = GeneticCode::from_table_number(27).unwrap();
        // Resolved to amino-acid reading, not stop.
        assert_eq!(tx1(code, b"TAA"), 'Q');
        assert_eq!(tx1(code, b"TAG"), 'Q');
        assert_eq!(tx1(code, b"TGA"), 'W');
        assert!(
            code.caveat().is_some(),
            "table 27 must flag its context-dependent limitation"
        );
    }

    #[test]
    fn table_28_condylostoma_has_caveat() {
        let code = GeneticCode::from_table_number(28).unwrap();
        assert_eq!(tx1(code, b"TAA"), 'Q');
        assert_eq!(tx1(code, b"TAG"), 'Q');
        assert_eq!(tx1(code, b"TGA"), 'W');
        assert!(code.caveat().is_some());
    }

    #[test]
    fn table_29_mesodinium() {
        let code = GeneticCode::from_table_number(29).unwrap();
        assert_eq!(tx1(code, b"TAA"), 'Y');
        assert_eq!(tx1(code, b"TAG"), 'Y');
        assert_eq!(tx1(code, b"TGA"), '*');
    }

    #[test]
    fn table_30_peritrich() {
        let code = GeneticCode::from_table_number(30).unwrap();
        assert_eq!(tx1(code, b"TAA"), 'E');
        assert_eq!(tx1(code, b"TAG"), 'E');
        assert_eq!(tx1(code, b"TGA"), '*');
    }

    #[test]
    fn table_31_blastocrithidia_has_caveat() {
        let code = GeneticCode::from_table_number(31).unwrap();
        assert_eq!(tx1(code, b"TAA"), 'E');
        assert_eq!(tx1(code, b"TAG"), 'E');
        assert_eq!(tx1(code, b"TGA"), 'W');
        assert!(code.caveat().is_some());
    }

    #[test]
    fn table_33_cephalodiscidae_mito() {
        let code = GeneticCode::from_table_number(33).unwrap();
        assert_eq!(tx1(code, b"AGA"), 'S');
        assert_eq!(tx1(code, b"AGG"), 'K');
        assert_eq!(tx1(code, b"TAA"), 'Y');
        assert_eq!(tx1(code, b"TGA"), 'W');
    }

    // ── Translation end-to-end ──────────────────────────────────────────────

    #[test]
    fn standard_translation() {
        assert_eq!(GeneticCode::STANDARD.translate("ATGAAATTT"), "MKF");
    }

    #[test]
    fn candida_ctg_vs_standard() {
        let ctg = GeneticCode::from_table_number(12).unwrap();
        assert_eq!(ctg.translate("ATGCTG"), "MS");
        assert_eq!(GeneticCode::STANDARD.translate("ATGCTG"), "ML");
    }

    #[test]
    fn vert_mito_tga_is_trp() {
        let code = GeneticCode::from_table_number(2).unwrap();
        assert_eq!(code.translate("ATGTGA"), "MW");
    }

    #[test]
    fn internal_stops_detected() {
        assert_eq!(
            GeneticCode::find_internal_stops("MK*FG*"),
            vec![2],
            "internal stop at AA position 2; trailing stop (position 5) ignored"
        );
    }

    // ── Start codons ────────────────────────────────────────────────────────

    #[test]
    fn atg_always_a_start_codon() {
        for code in GeneticCode::all() {
            assert!(
                code.is_start_codon(b"ATG"),
                "ATG must always start table {}",
                code.table_number()
            );
        }
    }

    #[test]
    fn bacterial_permits_extended_starts() {
        let code = GeneticCode::from_table_number(11).unwrap();
        assert!(code.is_start_codon(b"ATG"));
        assert!(code.is_start_codon(b"GTG"));
        assert!(code.is_start_codon(b"TTG"));
    }

    // ── Sweep-all correctness guard ─────────────────────────────────────────

    #[test]
    fn every_table_translates_atg_to_met() {
        for code in GeneticCode::all() {
            assert_eq!(
                code.translate("ATG"),
                "M",
                "table {} broke ATG=Met",
                code.table_number()
            );
        }
    }

    #[test]
    fn no_table_has_empty_display_name() {
        for code in GeneticCode::all() {
            let name = code.display_name();
            assert!(
                !name.trim().is_empty(),
                "table {} has empty display_name",
                code.table_number()
            );
        }
    }
}
