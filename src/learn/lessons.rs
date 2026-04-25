/// Lesson content for the interactive tutorial
///
/// Each lesson is a sequence of items: text blocks, questions,
/// code examples, live demos, and checkpoints. Lessons are ordered
/// from basic to advanced, following the pipeline order.
///
/// Design principles:
///   - Teach concepts BEFORE asking questions
///   - Explain WHY each tool was chosen (scientific rationale)
///   - Include references to key papers
///   - Show HOW to use each command with real examples
///   - Build understanding progressively

// ─────────────────────────────────────────────────────────────────────────────
// Data structures
// ─────────────────────────────────────────────────────────────────────────────

pub struct Lesson {
    pub title: &'static str,
    pub description: &'static str,
    pub est_minutes: u8,
    pub items: Vec<LessonItem>,
}

pub enum LessonItem {
    /// Explanatory text block (markdown-ish)
    Text(&'static str),

    /// Info panel with title and body (displayed in a box)
    Info(&'static str, &'static str),

    /// Code example: description + the command(s)
    CodeExample(&'static str, &'static str),

    /// Live demo — shows a command the user would run
    Demo(&'static str, &'static str),

    /// Interactive question
    Question {
        prompt: &'static str,
        kind: QuestionKind,
        hint: Option<&'static str>,
        explanation: Option<&'static str>,
    },

    /// Hands-on exercise — user types a command
    TryIt {
        instruction: &'static str,
        command_template: &'static str,
        hint: Option<&'static str>,
    },

    /// Progress checkpoint
    Checkpoint(&'static str),
}

pub enum QuestionKind {
    FreeText {
        answer: &'static str,
        accept_regex: Option<&'static str>,
    },
    MultipleChoice {
        choices: Vec<String>,
        correct_index: usize,
    },
    TrueFalse {
        answer: bool,
    },
    FillBlank {
        template: &'static str,
        answer: &'static str,
    },
    OrderSteps {
        steps: Vec<String>,
        correct_order: Vec<usize>,
    },
}

// ─────────────────────────────────────────────────────────────────────────────
// Lesson lookup
// ─────────────────────────────────────────────────────────────────────────────

pub fn find_lesson_by_name(name: &str) -> Option<usize> {
    let lower = name.to_lowercase();
    all_lessons()
        .iter()
        .enumerate()
        .find(|(_, l)| l.title.to_lowercase().contains(&lower))
        .map(|(i, _)| i + 1)
}

// ─────────────────────────────────────────────────────────────────────────────
// All lessons
// ─────────────────────────────────────────────────────────────────────────────

pub fn all_lessons() -> Vec<Lesson> {
    vec![
        // ── Foundation track (run in order if you're new) ────────────────────
        lesson_1_welcome(),
        lesson_2_setup(),
        lesson_3_sort_mask(),
        lesson_4_predict(),
        lesson_5_annotate(),
        lesson_6_advanced(),
        lesson_7_submit(),
        lesson_8_analysis(),
        // ── RNA-seq track (lessons 9-14, completable in any order) ──────────
        lesson_9_quant(),
        lesson_10_fetch_rna(),
        lesson_11_de_template(),
        lesson_12_ase(),
        lesson_13_ase_template(),
        lesson_14_go_template(),
        // ── Modern predictors track (lessons 15-16) ─────────────────────────
        lesson_15_braker(),
        lesson_16_genemark(),
        // ── Quality and reproducibility (lessons 17-20) ─────────────────────
        lesson_17_clean_contigs(),
        lesson_18_augustus_fungi(),
        lesson_19_reproducibility(),
        lesson_20_genetic_codes(),
    ]
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 1: Welcome & Orientation
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_1_welcome() -> Lesson {
    Lesson {
        title: "Welcome to MycoNote",
        description: "The annotation problem, how myconote-cli solves it, and core genomics concepts.",
        est_minutes: 12,
        items: vec![
            LessonItem::Text(
"THE ANNOTATION BOTTLENECK

As genome sequencing costs have plummeted (a fungal genome now costs
~$500 to sequence), the bottleneck in genomics has shifted from
generating sequence data to making sense of it. A raw genome assembly
is just a string of A, T, C, G — it tells you nothing about where the
genes are or what they do.

Genome annotation is the process of:
  1. Finding genes in the raw sequence (structural annotation)
  2. Determining what those genes do (functional annotation)

This is computationally intensive, requires multiple specialized tools,
and traditionally takes days of manual bioinformatics work. myconote-cli
automates this entire process into a single pipeline."
            ),

            LessonItem::Info("Why was myconote-cli created?",
"Existing annotation pipelines like funannotate (Palmer & Stajich, 2020)
are widely used but have limitations:
  - Written in Python, which is slow for I/O-heavy genomics tasks
  - Limited to fungi, with poor support for other kingdoms
  - No built-in output validation before NCBI submission
  - No genetic code support for organisms like the Candida CTG clade
  - No reproducibility tracking

myconote-cli addresses all of these gaps with a Rust implementation
that is faster, broader, and produces submission-ready outputs."
            ),

            LessonItem::Text(
"WHAT IS A GENOME ASSEMBLY?

When you sequence a genome, you don't get one long continuous sequence.
Instead, you get millions of short reads that are assembled into longer
contiguous sequences called 'contigs' or 'scaffolds'. These are stored
in FASTA format — the universal sequence format in bioinformatics.

A FASTA file looks like this:
  >scaffold_1
  ATCGATCGATCG...
  >scaffold_2
  GCTAGCTAGCTA...

Each '>' line is a header, followed by the DNA sequence."
            ),

            LessonItem::Question {
                prompt: "What file format does a genome assembly typically come in?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "GFF3 (Gene Feature Format)".into(),
                        "FASTA (sequence format with > headers)".into(),
                        "GenBank (NCBI annotated format)".into(),
                        "VCF (Variant Call Format)".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It's a simple format with > headers followed by sequences."),
                explanation: Some("FASTA (.fa, .fasta, .fna) stores raw sequences. GFF3 stores annotations ABOUT sequences. GenBank combines both. VCF stores variants."),
            },

            LessonItem::Text(
"THE MYCONOTE PIPELINE

myconote-cli organizes annotation into 7 sequential stages. Each stage
builds on the output of the previous one:

  1. sort     — Clean up contig names and order them by length
  2. mask     — Mark repetitive DNA so gene finders can ignore it
  3. train    — (Optional) Train gene predictors on your organism's RNA
  4. predict  — Find genes using multiple ab initio predictors
  5. update   — Refine gene boundaries using transcript evidence
  6. annotate — Assign function to each gene (what does it do?)
  7. submit   — Validate and prepare files for NCBI GenBank

This is the same conceptual workflow used by every genome project,
but myconote-cli handles it in a single tool instead of requiring
you to manually chain 20+ separate programs."
            ),

            LessonItem::Question {
                prompt: "Which step comes FIRST in the pipeline?",
                kind: QuestionKind::FreeText {
                    answer: "sort",
                    accept_regex: Some(r"(?i)^sort"),
                },
                hint: Some("We need clean, sorted contig names before anything else."),
                explanation: Some("'sort' standardizes contig names (scaffold_1, scaffold_2...) so all downstream tools see consistent identifiers."),
            },

            LessonItem::Text(
"WHAT IS SOFT-MASKING?

Eukaryotic genomes contain large amounts of repetitive DNA —
transposable elements, tandem repeats, and satellite sequences.
In fungi, 3-20% of the genome is repetitive; in plants, it can
be over 80%.

These repeats confuse gene prediction algorithms because they
look like real genes (they often contain open reading frames).

Soft-masking converts repetitive bases to lowercase letters:
  Before: ATCG ATCG ATCG (repeats look like real sequence)
  After:  atcg atcg atcg (gene finders know to skip these)

This is different from hard-masking (replacing with NNNN), which
permanently removes the sequence information."
            ),

            LessonItem::Question {
                prompt: "What does 'soft-masking' a genome mean?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Deleting repetitive sequences from the genome".into(),
                        "Replacing repeats with N characters (hard-masking)".into(),
                        "Converting repetitive bases to lowercase letters".into(),
                        "Compressing the FASTA file to save disk space".into(),
                    ],
                    correct_index: 2,
                },
                hint: Some("Think lowercase vs uppercase in a FASTA file."),
                explanation: Some("Soft-masking uses lowercase (atcg) for repeats. Gene predictors see them but know to be cautious. The sequence data is preserved, unlike hard-masking."),
            },

            LessonItem::Checkpoint("Core concepts understood"),

            LessonItem::Text(
"KINGDOM-AWARE DEFAULTS

Different eukaryotic lineages have dramatically different genome
architectures. myconote-cli adjusts its parameters automatically:

  Kingdom   Typical size   Intron range    Gene density     Example
  -------   ------------   ------------    ------------     -------
  fungi     12-40 Mb       40-2,000 bp     300-500/Mb       S. cerevisiae
  plant     100-3,000 Mb   40-50,000 bp    30-100/Mb        A. thaliana
  animal    500-3,000 Mb   40-500,000 bp   5-15/Mb          H. sapiens
  insect    100-2,000 Mb   40-50,000 bp    50-100/Mb        D. melanogaster
  protist   10-200 Mb      20-1,000 bp     200-600/Mb       T. gondii

Fungal genomes are compact with short introns — that's why myconote-cli
defaults to 'fungi'. Using the wrong kingdom causes gene finders to use
inappropriate intron size models, leading to missed or merged genes.

Ref: Stajich et al. (2009) 'The Fungi' Curr Biol 19(18):R840-5"
            ),

            LessonItem::Question {
                prompt: "For a Saccharomyces cerevisiae genome, which kingdom would you use?",
                kind: QuestionKind::FreeText {
                    answer: "fungi",
                    accept_regex: Some(r"(?i)^fung"),
                },
                hint: Some("Yeast is a fungus!"),
                explanation: Some("S. cerevisiae is a budding yeast — kingdom fungi. This is the default, so you don't even need to specify it."),
            },

            LessonItem::CodeExample(
                "See all available commands and their descriptions:",
                "myconote-cli --help"
            ),

            LessonItem::Question {
                prompt: "True or False: myconote-cli is focused primarily on fungal genomes.",
                kind: QuestionKind::TrueFalse { answer: true },
                hint: Some("Defaults, benchmarks, and validation panels target fungi."),
                explanation: Some("True. As of v0.5+, fungi is the primary supported use case — defaults, benchmarks, and the test fixtures are tuned for fungal annotation. The `--kingdom` flag also accepts plant, animal, insect, and protist values, but those paths are experimental and not validated on large genomes."),
            },
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 2: Setup & Installation
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_2_setup() -> Lesson {
    Lesson {
        title: "Setup & Installation",
        description: "Why each external tool is needed, how to install them, and database setup.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"THE TOOL ECOSYSTEM

Genome annotation requires many specialized algorithms, each designed
by different research groups over decades. myconote-cli orchestrates
30+ external tools. Here's WHY each major tool was chosen:

GENE PREDICTORS (finding genes in DNA):
  Augustus  — Generalized Hidden Markov Model with species-specific
              training. The gold standard for eukaryotic gene finding.
              Ref: Stanke et al. (2006) BMC Bioinformatics 7:62

  SNAP      — Semi-HMM-based predictor. Faster than Augustus, good
              as a secondary opinion for consensus building.
              Ref: Korf (2004) BMC Bioinformatics 5:59

  GeneMark  — Self-training predictor that doesn't need a pre-built
              model. Useful for novel organisms with no close relatives.
              Ref: Lomsadze et al. (2005) Nucleic Acids Res 33:6494

HOMOLOGY SEARCH (what are the genes?):
  MMseqs2   — 100x faster than BLAST with similar sensitivity.
              Used for protein homology against Swiss-Prot.
              Ref: Steinegger & Soding (2017) Nature Biotech 35:1026

  HMMER     — Profile hidden Markov models for domain detection.
              Used for Pfam domain annotation.
              Ref: Eddy (2011) PLoS Comput Biol 7:e1002195"
            ),

            LessonItem::Text(
"COMPLETENESS ASSESSMENT:
  BUSCO     — Benchmarking Universal Single-Copy Orthologs.
              Checks if expected conserved genes are present.
              A BUSCO score of 95%+ means your annotation is good.
              Ref: Manni et al. (2021) Mol Biol Evol 38:4647

REPEAT MASKING:
  RepeatModeler — Builds a de novo repeat library from your genome.
                  No reference database needed — discovers novel repeats.
  RepeatMasker  — Masks repeats using known repeat libraries.
                  Ref: Smit et al. (2013-2015) RepeatMasker Open-4.0

RNA-seq TRAINING:
  Trinity   — De novo transcript assembly from RNA-seq reads.
              Produces the transcript evidence needed to train predictors.
              Ref: Grabherr et al. (2011) Nature Biotech 29:644

  PASA      — Program to Assemble Spliced Alignments. Builds a
              transcript database and refines gene models with UTRs.
              Ref: Haas et al. (2003) Nucleic Acids Res 31:5654"
            ),

            LessonItem::CodeExample(
                "Check which tools are installed on your system:",
                "myconote-cli check"
            ),

            LessonItem::Question {
                prompt: "What command installs all missing external tools automatically?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli ___",
                    answer: "install",
                },
                hint: Some("It's a single word that means 'put these tools on my system'."),
                explanation: Some("'myconote-cli install' scans PATH for missing tools and installs them via conda/mamba. Tools requiring licences (GeneMark, SignalP) are flagged for manual download."),
            },

            LessonItem::Info("Package managers: conda vs mamba",
"Conda is the standard package manager for bioinformatics tools.
Mamba is a faster drop-in replacement that uses libsolv instead
of conda's Python-based dependency solver.

myconote-cli prefers mamba when available (10-50x faster installs)
but falls back to conda automatically. Both install from the
Bioconda channel — the community repository for bioinformatics
software with 8,000+ packages.

Ref: Gruning et al. (2018) Nature Methods 15:475 (Bioconda)"
            ),

            LessonItem::Text(
"REFERENCE DATABASES

Annotation doesn't happen in a vacuum — you're comparing your genes
against curated reference databases to infer function:

  Swiss-Prot    ~570,000 manually curated protein sequences with
                verified function. The highest-quality protein DB.
                Downloaded and indexed with MMseqs2 for fast search.

  Pfam-A        ~20,000 protein domain families represented as
                Hidden Markov Models. Tells you the domain architecture
                of each protein (e.g., 'has a kinase domain + SH2 domain').
                Ref: Mistry et al. (2021) Nucleic Acids Res 49:D412

  BUSCO lineages  Sets of single-copy orthologous genes expected
                  in a given lineage (fungi_odb10 has 758 genes).
                  Used to assess annotation completeness.

  dbCAN         CAZyme (Carbohydrate-Active Enzyme) HMM profiles.
                Critical for fungi that degrade plant biomass.

  MEROPS        Protease/peptidase classification database.
                Important for pathogenicity studies."
            ),

            LessonItem::Question {
                prompt: "What command downloads and indexes all reference databases?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli ___",
                    answer: "setup",
                },
                hint: Some("You're setting up the databases your annotations will search against."),
                explanation: Some("'myconote-cli setup' downloads Swiss-Prot, Pfam, dbCAN, MEROPS, and BUSCO lineages (~2.5 GB total). Each database is version-tracked for reproducibility."),
            },

            LessonItem::CodeExample(
                "Download only specific databases (saves time/space):",
                "myconote-cli setup --db swiss-prot pfam"
            ),

            LessonItem::Question {
                prompt: "Where are databases stored by default?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "/usr/local/share/myconote/".into(),
                        "~/.myconote/dbs/".into(),
                        "./databases/".into(),
                        "/tmp/myconote/".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It's in your home directory, hidden (starts with a dot)."),
                explanation: Some("~/.myconote/dbs/ is the default. Override with --db-dir. Version metadata is tracked so you know when each database was downloaded."),
            },

            LessonItem::Checkpoint("Environment setup understood"),

            LessonItem::Question {
                prompt: "Why does myconote-cli use MMseqs2 instead of BLAST for homology search?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "BLAST is no longer maintained".into(),
                        "MMseqs2 is ~100x faster with similar sensitivity".into(),
                        "BLAST cannot search protein databases".into(),
                        "MMseqs2 produces GFF3 output natively".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("Think about speed. Searching 5,000 proteins against Swiss-Prot takes ~90 seconds with MMseqs2."),
                explanation: Some("MMseqs2 uses a prefilter + k-mer matching strategy that achieves ~100x speedup over BLAST with comparable sensitivity (Steinegger & Soding, 2017)."),
            },
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 3: Sort & Mask
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_3_sort_mask() -> Lesson {
    Lesson {
        title: "Pre-processing: Sort & Mask",
        description: "Why pre-processing matters, how repeat masking works, and which engine to choose.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"WHY SORT?

Genome assemblers produce contigs with messy names like:
  >NODE_1_length_1856900_cov_45.2
  >k141_38927|size=500231

These cause problems downstream because:
  - Different tools truncate long headers differently
  - Special characters (|, =, spaces) break GFF3 parsers
  - NCBI rejects submissions with non-standard sequence IDs

The 'sort' command standardizes everything:
  >scaffold_001    (longest contig first)
  >scaffold_002
  ...

It also filters short contigs that are too small to contain genes
(typically < 500 bp) and writes a rename table so you can trace
the original names."
            ),

            LessonItem::CodeExample(
                "Sort contigs, remove anything under 500 bp, save rename mapping:",
                "myconote-cli sort assembly.fa --min-length 500 --rename-table rename.tsv"
            ),

            LessonItem::Question {
                prompt: "Why do we sort contigs by length (longest first)?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "It makes the FASTA file smaller".into(),
                        "Longest contigs usually contain the most genes".into(),
                        "It ensures consistent naming and helps downstream tools".into(),
                        "Gene predictors refuse unsorted input".into(),
                    ],
                    correct_index: 2,
                },
                hint: Some("Think about what happens when different tools see different contig names."),
                explanation: Some("Clean, consistent IDs prevent naming conflicts in GFF3 output, BAM alignments, and NCBI submission. The length-sorting is a convention that puts the most informative contigs first."),
            },

            LessonItem::Text(
"WHY MASK REPEATS?

Repetitive DNA is the enemy of gene prediction. Here's why:

  1. Transposable elements (TEs) contain their own genes (reverse
     transcriptase, integrases, etc.). Without masking, Augustus
     will predict thousands of spurious TE-derived 'genes'.

  2. Tandem repeats create frameshifts that confuse codon-based
     predictors, producing truncated or chimeric gene models.

  3. Segmental duplications cause the same gene to be predicted
     multiple times with conflicting coordinates.

In a typical 13 Mb fungal genome, 5-15% of the sequence is
repetitive. In plants, it can exceed 80%.

Ref: Wicker et al. (2007) Nature Reviews Genetics 8:973"
            ),

            LessonItem::Text(
"MASKING ENGINES — WHICH TO USE?

myconote-cli offers 5 masking strategies:

  self           Uses minimap2 to align the genome against itself.
                 Regions with multiple self-hits are masked.
                 PROS: No database needed, very fast.
                 CONS: Misses low-copy repeats. Best for quick checks.

  repeatmasker   Uses RepeatMasker with species-specific repeat
                 libraries from Dfam/RepBase.
                 PROS: Well-validated, fast.
                 CONS: Requires known repeat library for your species.

  repeatmodeler  Builds a DE NOVO repeat library using RepeatModeler,
                 then feeds it to RepeatMasker.
                 PROS: Best for novel genomes — discovers new repeats.
                 CONS: Slowest option (1-2 hours for a fungal genome).
                 THIS IS THE RECOMMENDED DEFAULT.

  both           RepeatMasker + self-alignment merged.
  full           RepeatModeler + RepeatMasker + self-alignment (most thorough).

Ref: Flynn et al. (2020) PNAS 117:9451 (RepeatModeler2)"
            ),

            LessonItem::TryIt {
                instruction: "Write the command to mask a genome using the repeatmodeler engine with 8 threads:",
                command_template: "myconote-cli mask genome.fa --engine ___ --threads ___",
                hint: Some("The engine is 'repeatmodeler' and threads is 8."),
            },

            LessonItem::Question {
                prompt: "What is the RECOMMENDED masking engine for a new, uncharacterized fungal genome?",
                kind: QuestionKind::FreeText {
                    answer: "repeatmodeler",
                    accept_regex: Some(r"(?i)repeatmodel"),
                },
                hint: Some("It builds a de novo repeat library specific to YOUR genome."),
                explanation: Some("repeatmodeler discovers novel repeat families without requiring a pre-existing library. Essential for non-model organisms where repeat databases are incomplete."),
            },

            LessonItem::Checkpoint("Pre-processing mastered"),

            LessonItem::Question {
                prompt: "True or False: Hard-masking (replacing repeats with N) is preferred for gene prediction.",
                kind: QuestionKind::TrueFalse { answer: false },
                hint: Some("If you replace sequence with Ns, you lose information permanently."),
                explanation: Some("False. Soft-masking (lowercase) is preferred because: (1) gene finders can still read the sequence if needed, (2) no information is lost, and (3) some real genes overlap with or are embedded in repeat regions."),
            },
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 4: Gene Prediction
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_4_predict() -> Lesson {
    Lesson {
        title: "Gene Prediction",
        description: "How ab initio prediction works, why we use multiple predictors, and evidence-based consensus.",
        est_minutes: 14,
        items: vec![
            LessonItem::Text(
"HOW GENE PREDICTION WORKS

Gene prediction algorithms scan the DNA sequence looking for signals
that indicate where genes begin and end:

  - Start codons (ATG) and stop codons (TAA, TAG, TGA)
  - Splice sites (GT...AG at intron boundaries)
  - Codon usage bias (real genes use codons differently than random DNA)
  - Promoter signals (TATA box, etc.)

These 'ab initio' (from the beginning) predictors use Hidden Markov
Models (HMMs) trained on known genes from your organism or a relative.

The key insight: NO SINGLE PREDICTOR IS PERFECT. Each has biases:
  - Augustus tends to predict longer genes (overshoots UTRs)
  - SNAP is fast but predicts more false positives
  - GeneMark works without training but misses short genes

That's why myconote-cli uses multiple predictors and merges them."
            ),

            LessonItem::Info("Why multiple predictors?",
"Studies consistently show that combining predictions from 2-3 tools
produces better gene sets than any single tool alone. This is called
'consensus gene calling' or 'evidence-based prediction'.

Hoff et al. (2016) showed that combining Augustus + GeneMark with
protein evidence found 15-25% more correct genes than either alone.

Ref: Hoff et al. (2016) Bioinformatics 32:767 (BRAKER)
Ref: Haas et al. (2008) Genome Biol 9:R7 (EvidenceModeler)"
            ),

            LessonItem::Text(
"THE PREDICTORS AND THEIR ROLES

  Augustus (weight: 10)
    The most accurate ab initio eukaryotic gene finder. Uses a
    generalized HMM with 'hints' from external evidence (protein
    alignments, RNA-seq). Has pre-trained models for ~100 species.
    Ref: Stanke et al. (2006) BMC Bioinformatics 7:62

  SNAP (weight: 3)
    Semi-HMM predictor. Much faster than Augustus but less accurate.
    Provides a 'second opinion' — genes found by both Augustus AND
    SNAP are very likely real.
    Ref: Korf (2004) BMC Bioinformatics 5:59

  GlimmerHMM (weight: 2)
    Interpolated Markov model predictor. Good at finding exons but
    can struggle with multi-exon genes. Optional third predictor.
    Ref: Majoros et al. (2004) Bioinformatics 20:2878

  GeneMark-ES (weight: 5)
    Self-training predictor that iteratively learns gene structure
    from the genome itself — no pre-trained model needed. Excellent
    for truly novel organisms with no close relatives.
    Ref: Lomsadze et al. (2005) Nucleic Acids Res 33:6494

  Protein evidence (weight: 20)
    Aligning known proteins to your genome using miniprot. This is
    the STRONGEST evidence — if a known protein aligns with introns
    to your genome, there's almost certainly a gene there.
    Ref: Li (2023) Bioinformatics 39:btad014 (miniprot)"
            ),

            LessonItem::Question {
                prompt: "Which evidence type gets the HIGHEST weight by default?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Augustus ab initio predictions".into(),
                        "SNAP predictions".into(),
                        "Protein alignment evidence".into(),
                        "GlimmerHMM predictions".into(),
                    ],
                    correct_index: 2,
                },
                hint: Some("Homology to known proteins is the strongest evidence that a gene exists."),
                explanation: Some("Protein evidence (weight 20) is strongest because if a known, curated protein aligns to your genome, it's near-certain a real gene is there. Ab initio predictions are probabilistic guesses."),
            },

            LessonItem::Text(
"THE EVIDENCE MODELER (EVM) CONSENSUS

After running all predictors, myconote-cli must decide: when two
predictors disagree about a gene's boundaries, which one is right?

The Evidence Modeler algorithm (Haas et al. 2008):
  1. Groups all predictions that overlap on the same strand
  2. Scores each prediction by its source weight
  3. For each overlap group, keeps the highest-scoring model
  4. Assigns clean sequential locus tags (e.g., MYORG_000001)

You can customize weights in a TOML file if you trust some predictors
more than others for your specific organism."
            ),

            LessonItem::CodeExample(
                "Basic gene prediction for a fungal genome:",
                "myconote-cli predict genome_masked.fa --kingdom fungi --locus-prefix MYORG"
            ),

            LessonItem::Question {
                prompt: "What does --locus-prefix control?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "The output directory name".into(),
                        "The gene ID naming scheme (e.g. MYORG_000001)".into(),
                        "The Augustus species model to use".into(),
                        "The BUSCO lineage for assessment".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It's a prefix that goes before each gene number."),
                explanation: Some("--locus-prefix MYORG gives genes IDs like MYORG_000001, MYORG_000002. NCBI requires a registered locus_tag prefix for genome submissions."),
            },

            LessonItem::Info("Custom evidence weights (TOML)",
"Create a weights.toml file to customize predictor trust:

  augustus = 10.0    # most trusted ab initio
  snap = 3.0        # secondary predictor
  protein = 25.0    # increase if you have good protein evidence
  est = 8.0         # transcript evidence
  genemark = 5.0    # self-training predictor
  glimmerhmm = 2.0  # third-tier predictor

Run with: myconote-cli predict genome.fa --weights weights.toml

Weights are saved to evidence_weights.toml in the output directory
for reproducibility."
            ),

            LessonItem::TryIt {
                instruction: "Write the command to predict genes with protein evidence from Swiss-Prot:",
                command_template: "myconote-cli predict masked.fa --kingdom fungi --protein-fasta ___",
                hint: Some("Point --protein-fasta to your Swiss-Prot FASTA file."),
            },

            LessonItem::Text(
"PLOIDY AND ALLELIC DUPLICATES

If your organism is diploid (2n) or polyploid, the genome assembly
may contain both copies of each chromosome. This means every gene
appears TWICE — inflating your gene count.

myconote-cli's --ploidy flag handles this:
  --ploidy 2   Expects allelic pairs. Adjusts overlap tolerance and
               reports allelic duplicates so you can collapse or keep them.

This is especially important for Candida albicans (diploid), many
plants (polyploid), and industrial yeast strains.

Ref: Todd et al. (2017) Nat Rev Microbiol 15:96 (Candida genomics)"
            ),

            LessonItem::Question {
                prompt: "What does the --ploidy flag do?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Sets the number of CPU threads".into(),
                        "Controls allelic duplicate handling in polyploid genomes".into(),
                        "Specifies the sequencing depth".into(),
                        "Selects the repeat masking engine".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("Think about what happens when you assemble a diploid genome — you get two copies of each gene."),
                explanation: Some("--ploidy 2 tells myconote-cli to expect allelic duplicates and adjusts the overlap tolerance in the evidence merger accordingly."),
            },

            LessonItem::Checkpoint("Gene prediction understood"),

            LessonItem::Question {
                prompt: "Put these pipeline steps in the correct order:",
                kind: QuestionKind::OrderSteps {
                    steps: vec![
                        "annotate".into(),
                        "predict".into(),
                        "mask".into(),
                        "sort".into(),
                    ],
                    correct_order: vec![3, 2, 1, 0],
                },
                hint: Some("Start with sort, end with annotate."),
                explanation: Some("sort -> mask -> predict -> annotate. Each step requires the output of the previous one."),
            },
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 5: Functional Annotation
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_5_annotate() -> Lesson {
    Lesson {
        title: "Functional Annotation",
        description: "How each annotation source works, why it matters, and genetic code biology.",
        est_minutes: 15,
        items: vec![
            LessonItem::Text(
"FROM GENES TO FUNCTION

At this point you have a GFF3 file with gene coordinates — but no idea
what any of them DO. Functional annotation assigns biological meaning
by comparing your proteins against curated databases.

myconote-cli integrates 15 annotation sources. Let's understand each:"
            ),

            LessonItem::Text(
"CORE ANNOTATIONS (always run):

  MMseqs2 vs Swiss-Prot
    What: Searches each predicted protein against UniProt/Swiss-Prot,
          a database of ~570,000 manually curated proteins.
    Why:  Gives you a product name ('Alcohol dehydrogenase 1') and
          a UniProt accession for each gene.
    How:  MMseqs2 uses precomputed k-mer indexes for ~100x speedup
          over BLAST, with comparable sensitivity.
    Ref:  Steinegger & Soding (2017) Nature Biotech 35:1026

  Pfam domains (hmmsearch)
    What: Searches proteins against ~20,000 Pfam-A domain models.
    Why:  Tells you the domain architecture — a protein with a
          'Kinase' domain + 'SH2' domain is likely a signaling kinase.
    How:  We use hmmsearch (not hmmscan) — it searches each HMM
          profile against all your proteins at once, giving a
          ~6x speedup for large proteomes.
    Ref:  Mistry et al. (2021) Nucleic Acids Res 49:D412

  BUSCO completeness
    What: Checks if expected single-copy orthologs are present.
    Why:  A BUSCO score tells you how complete your annotation is.
          95%+ is good; <80% suggests missing genes or assembly gaps.
    How:  Uses lineage-specific gene sets (e.g., fungi_odb10 = 758 genes).
    Ref:  Manni et al. (2021) Mol Biol Evol 38:4647

  GO terms
    What: Gene Ontology terms (biological process, molecular function,
          cellular component) assigned from UniProt matches.
    Why:  Enables enrichment analysis — 'are DNA repair genes
          overrepresented in my differentially expressed set?'
    Ref:  Gene Ontology Consortium (2021) Nucleic Acids Res 49:D325"
            ),

            LessonItem::Question {
                prompt: "Which tool does myconote-cli use for protein homology search against Swiss-Prot?",
                kind: QuestionKind::FreeText {
                    answer: "MMseqs2",
                    accept_regex: Some(r"(?i)mmseqs"),
                },
                hint: Some("It's ~100x faster than BLAST and starts with 'MM'."),
                explanation: Some("MMseqs2 achieves BLAST-like sensitivity at ~100x the speed using precomputed k-mer indexes and a two-stage search strategy."),
            },

            LessonItem::Text(
"OPTIONAL ANNOTATIONS (enable with flags):

  --trnascan   tRNA gene prediction via tRNAscan-SE
    What: Finds transfer RNA genes using covariance models.
    Why:  tRNAs are essential non-coding genes. Most fungi have
          150-500 tRNA genes. They're missed by protein-coding
          gene finders because they don't produce proteins.
    Ref:  Chan & Lowe (2019) Methods Mol Biol 1962:1

  --eggnog     EggNOG-mapper (COG/NOG functional categories)
    What: Assigns Clusters of Orthologous Groups and KEGG pathways.
    Why:  Provides broad functional categories (metabolism, signaling,
          etc.) and maps genes to metabolic pathways.
    Ref:  Cantalapiedra et al. (2021) Mol Biol Evol 38:5825

  --cazyme     CAZyme annotation via dbCAN
    What: Classifies carbohydrate-active enzymes (glycoside hydrolases,
          glycosyl transferases, etc.).
    Why:  Critical for fungi — CAZymes determine what substrates an
          organism can degrade. Key for biofuel, food, and pathogenicity.
    Ref:  Zheng et al. (2023) Nucleic Acids Res 51:D557

  --secretome  Signal peptide + transmembrane prediction
    What: Identifies secreted proteins (SignalP/DeepSig + DeepTMHMM).
    Why:  Secreted proteins include enzymes, toxins, and effectors.
          Key for understanding host-pathogen interactions.

  --antismash  Biosynthetic gene cluster detection
    What: Finds secondary metabolite clusters (PKS, NRPS, terpene).
    Why:  Fungi produce antibiotics, toxins, and pigments from BGCs.
    Ref:  Blin et al. (2023) Nucleic Acids Res 51:W46

  --merops     Protease classification
    What: Classifies peptidases into MEROPS families and clans.
    Why:  Proteases are drug targets and virulence factors.

  --interproscan  Comprehensive domain search (EBI REST API)
    What: Searches InterPro, which integrates 14 member databases
          (Pfam, TIGRFAM, SMART, Gene3D, Superfamily, etc.).
    Why:  The most thorough domain annotation available. Achieves
          98%+ annotation rate but requires internet access.
    Ref:  Blum et al. (2021) Nucleic Acids Res 49:D344"
            ),

            LessonItem::Question {
                prompt: "What does BUSCO measure?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Gene expression levels across tissues".into(),
                        "Genome completeness using conserved single-copy orthologs".into(),
                        "The number of tRNA genes in the genome".into(),
                        "Protein folding accuracy predictions".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It tells you how 'complete' your annotation is."),
                explanation: Some("BUSCO checks for 758 genes (in fungi_odb10) that should be present in single copy in every fungal genome. Missing or duplicated BUSCOs indicate assembly or annotation problems."),
            },

            LessonItem::CodeExample(
                "Run annotation with all sources enabled:",
                "myconote-cli annotate genes.gff3 --fasta genome.fa \\\n  --eggnog --cazyme --secretome --antismash --merops --trnascan \\\n  --interproscan --email you@email.edu"
            ),

            LessonItem::Text(
"GENETIC CODE: WHY IT MATTERS

The standard genetic code translates 64 codons into 20 amino acids.
But some organisms use DIFFERENT codes — and getting this wrong
silently corrupts every protein translation in your annotation.

The most important example for mycologists:

  CANDIDA CTG CLADE (Table 12)
    In ~400 yeast species, the codon CTG encodes SERINE instead of
    LEUCINE. If you annotate C. albicans with the standard code,
    every CTG in every gene gets the wrong amino acid.

    Species affected: C. albicans, C. tropicalis, C. parapsilosis,
    Debaryomyces hansenii, Meyerozyma guilliermondii, and ~400 others.

    Ref: Muhlhausen et al. (2016) Curr Opin Microbiol 32:16

  Other non-standard codes in fungi:
    Table 3  — Yeast mitochondrial (TGA = Trp, CTN = Thr)
    Table 26 — Pachysolen tannophilus nuclear (CTG = Ala)

myconote-cli supports 25 NCBI translation tables — tables 1–6, 9–14,
16, 21–31, and 33 — more than any other annotation pipeline."
            ),

            LessonItem::Question {
                prompt: "If you're annotating a Candida albicans genome, which genetic code should you use?",
                kind: QuestionKind::FreeText {
                    answer: "12",
                    accept_regex: Some(r"(?i)(12|candida|ctg)"),
                },
                hint: Some("Candida uses the CTG clade code where CTG encodes Serine, not Leucine."),
                explanation: Some("Table 12 (Alternative Yeast Nuclear Code). Without this flag, every CTG codon is mistranslated as Leu instead of Ser — silently corrupting your protein sequences."),
            },

            LessonItem::Checkpoint("Functional annotation mastered"),

            LessonItem::Question {
                prompt: "Which flag enables tRNA gene prediction?",
                kind: QuestionKind::FreeText {
                    answer: "--trnascan",
                    accept_regex: Some(r"(?i)--(trna|trnascan)"),
                },
                hint: Some("It wraps the tRNAscan-SE tool."),
                explanation: Some("--trnascan runs tRNAscan-SE, which uses covariance models to find tRNA genes. These are non-coding genes that standard protein-coding gene finders miss entirely."),
            },

            LessonItem::Question {
                prompt: "True or False: If a tool is missing, myconote-cli will crash.",
                kind: QuestionKind::TrueFalse { answer: false },
                hint: Some("Think about 'graceful degradation'."),
                explanation: Some("False. myconote-cli degrades gracefully — if hmmsearch isn't installed, Pfam is skipped; if tRNAscan-SE is missing, tRNA prediction is skipped. You get results from whatever IS available."),
            },
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 6: Advanced Features
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_6_advanced() -> Lesson {
    Lesson {
        title: "Advanced: Training & Evidence",
        description: "RNA-seq training theory, PASA updates, and protein-to-genome alignment.",
        est_minutes: 12,
        items: vec![
            LessonItem::Text(
"WHY TRAIN YOUR OWN MODELS?

Pre-trained Augustus models (e.g., 'saccharomyces_cerevisiae_S288C')
work well for organisms closely related to the training species. But
gene structure varies significantly even within fungi:

  - S. cerevisiae has very few introns (~5% of genes have introns)
  - Aspergillus has 2-3 introns per gene on average
  - Cryptococcus has 5-6 introns per gene

Using the wrong species model causes Augustus to either miss introns
(if your organism has more introns than the model expects) or predict
spurious introns (if your organism has fewer).

The solution: TRAIN a custom model on your own organism's gene
structures, derived from RNA-seq data."
            ),

            LessonItem::Text(
"THE TRAINING PIPELINE

myconote-cli's training workflow:

  1. Trinity assembly
     RNA-seq reads -> de novo transcript assembly
     Trinity uses a de Bruijn graph approach to reconstruct
     full-length transcripts from short reads.
     Ref: Grabherr et al. (2011) Nature Biotech 29:644

  2. minimap2 alignment
     Trinity transcripts -> aligned to genome
     minimap2 is used in splice-aware mode (-ax splice) to find
     where each transcript maps, including intron locations.
     Ref: Li (2018) Bioinformatics 34:3094

  3. PASA transcript database
     Aligned transcripts -> gene models with UTRs
     PASA (Program to Assemble Spliced Alignments) resolves
     overlapping transcript alignments and identifies complete
     gene models with 5' and 3' UTRs.
     Ref: Haas et al. (2003) Nucleic Acids Res 31:5654

  4. Augustus + SNAP training
     Complete gene models -> trained HMM parameters
     Uses the high-confidence PASA models to learn your organism's
     specific splice site signals, codon usage, and intron lengths."
            ),

            LessonItem::CodeExample(
                "Train predictors from paired-end RNA-seq:",
                "myconote-cli train genome_masked.fa \\\n  --left R1.fastq.gz --right R2.fastq.gz \\\n  --species my_organism --threads 8"
            ),

            LessonItem::Question {
                prompt: "What tool does myconote-cli use to assemble RNA-seq reads into transcripts?",
                kind: QuestionKind::FreeText {
                    answer: "Trinity",
                    accept_regex: Some(r"(?i)trinity"),
                },
                hint: Some("It's the most widely-used de novo RNA-seq assembler, named after a concept of three-in-one."),
                explanation: Some("Trinity (Grabherr et al. 2011) uses Inchworm, Chrysalis, and Butterfly modules to reconstruct transcripts from RNA-seq reads without a reference genome."),
            },

            LessonItem::Text(
"THE UPDATE STEP: ADDING UTRs

After gene prediction, the 'update' step refines gene models using
transcript evidence:

  - Extends 5' UTRs (untranslated regions before the start codon)
  - Extends 3' UTRs (after the stop codon)
  - Corrects intron/exon boundaries where transcript evidence disagrees
  - Discovers alternative isoforms (different mRNA splicing patterns)

This is done with PASA if available, or a lightweight minimap2-based
fallback that uses coverage data to extend UTR boundaries.

UTRs are important because they contain regulatory elements —
understanding UTR structure helps predict gene regulation."
            ),

            LessonItem::Question {
                prompt: "What is the 'update' step for?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Installing software updates for myconote-cli".into(),
                        "Refining gene models with UTRs and transcript evidence".into(),
                        "Downloading new versions of reference databases".into(),
                        "Updating the organism's taxonomy in NCBI".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It uses PASA to add 5'/3' UTRs and fix exon boundaries using RNA-seq evidence."),
                explanation: Some("'update' refines predicted gene models by adding UTRs, correcting splice sites, and discovering alternative isoforms using transcript alignments."),
            },

            LessonItem::Checkpoint("Advanced training understood"),

            LessonItem::Question {
                prompt: "What file format is used for custom evidence weights?",
                kind: QuestionKind::FreeText {
                    answer: "TOML",
                    accept_regex: Some(r"(?i)toml"),
                },
                hint: Some("It's a simple configuration format used in Rust projects (Cargo.toml uses it)."),
                explanation: Some("TOML (Tom's Obvious Minimal Language) is a human-readable config format. myconote-cli reads weights from .toml files and saves the weights used for each run for reproducibility."),
            },
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 7: NCBI Submission
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_7_submit() -> Lesson {
    Lesson {
        title: "NCBI Submission",
        description: "Why validation matters, common submission errors, and the table2asn workflow.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"WHY SUBMIT TO NCBI?

NCBI GenBank is the world's primary repository for annotated genome
sequences. Submitting your genome makes it:

  1. Publicly available for other researchers
  2. Searchable via BLAST (others can find homologs in your genome)
  3. Citable with an accession number in your publications
  4. Integrated into NCBI's taxonomy and gene databases

Most journals REQUIRE a GenBank accession before publication.

But NCBI has strict formatting requirements — and most annotation
pipelines produce output that needs manual cleanup before submission.
myconote-cli handles this automatically."
            ),

            LessonItem::Text(
"COMMON SUBMISSION ERRORS (AND HOW MYCONOTE CATCHES THEM)

myconote-cli validates your GFF3 BEFORE generating submission files:

  1. Duplicate feature IDs
     Problem: Two CDS features with the same ID='gene1-CDS'
     Cause:   Evidence merger bug (fixed in myconote-cli)
     Result:  NCBI rejects the submission silently

  2. Orphan Parent references
     Problem: A CDS references Parent='mRNA_42' but mRNA_42 doesn't exist
     Cause:   ID renaming during consensus merging
     Result:  table2asn crashes with a cryptic error

  3. Coordinate issues
     Problem: start > end for a feature
     Cause:   Strand confusion in GFF3 conversion
     Result:  NCBI validation error

  4. Missing qualifiers
     Problem: CDS without /product or gene without /locus_tag
     Cause:   Incomplete functional annotation
     Result:  NCBI returns warnings that delay acceptance

Other pipelines pass these errors to NCBI's tools and hope for the best.
myconote-cli catches them upfront with clear error messages."
            ),

            LessonItem::CodeExample(
                "Validate before committing to a full submission:",
                "myconote-cli submit annotated.gff3 --fasta genome.fa \\\n  --organism 'Aspergillus niger' --validate-only"
            ),

            LessonItem::Question {
                prompt: "What does myconote-cli check BEFORE generating submission files?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Whether the genome has been published before".into(),
                        "GFF3 compliance: duplicate IDs, orphan features, coordinate issues".into(),
                        "Whether you have a BioProject ID".into(),
                        "The organism's taxonomy in NCBI".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It validates the structural integrity of your GFF3 annotation."),
                explanation: Some("The validator checks for duplicate IDs, missing Parent references, coordinate issues, locus_tag presence, and cross-validates GFF3 seqids against the FASTA."),
            },

            LessonItem::Text(
"THE SUBMISSION WORKFLOW

  1. Register at NCBI (one-time)
     - Create a BioProject (groups all data for your study)
     - Create a BioSample (describes your organism/strain)
     - Request a locus_tag prefix (e.g., 'ASPNI')

  2. Validate
     myconote-cli submit genes.gff3 --fasta genome.fa \\
       --organism 'Aspergillus niger' --validate-only

  3. Generate submission files
     myconote-cli submit genes.gff3 --fasta genome.fa \\
       --organism 'Aspergillus niger' \\
       --strain 'CBS 513.88' \\
       --bioproject PRJNA123456 \\
       --biosample SAMN12345678 \\
       --locus-prefix ASPNI \\
       --genetic-code 1

  4. Upload to NCBI
     Submit the .sqn file at https://submit.ncbi.nlm.nih.gov/"
            ),

            LessonItem::Question {
                prompt: "What flag would you add to only validate without generating files?",
                kind: QuestionKind::FreeText {
                    answer: "--validate-only",
                    accept_regex: Some(r"(?i)--validate"),
                },
                hint: Some("You just want to check for errors, not produce output files."),
                explanation: Some("--validate-only runs all NCBI compliance checks and reports errors without generating .tbl or .sqn files. Always run this first."),
            },

            LessonItem::Info("BioProject & BioSample",
"Before submitting to NCBI, you must register:

  BioProject — Umbrella for your study. Gets an accession like PRJNA123456.
               Visit: https://submit.ncbi.nlm.nih.gov/subs/bioproject/

  BioSample  — Describes the biological material. Gets SAMN12345678.
               Visit: https://submit.ncbi.nlm.nih.gov/subs/biosample/

  locus_tag  — A unique prefix for your gene IDs. Register at:
               https://www.ncbi.nlm.nih.gov/genbank/genome_locus_tag/

These are free and take 1-2 business days to process."
            ),

            LessonItem::Checkpoint("Submission workflow complete"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 8: Stats, Comparative Genomics, and Format Conversion
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_8_analysis() -> Lesson {
    Lesson {
        title: "Stats, Comparative Genomics, and Conversion",
        description: "QC your annotation, run N-genome ortholog inference, convert between formats, and hand off to external visualization tools.",
        est_minutes: 12,
        items: vec![
            LessonItem::Text(
"BEYOND ANNOTATION: QC AND HANDOFF

After your pipeline finishes you have a GFF3, a FASTA, and a stack
of TSVs. Before you write the methods section you need to answer
three questions:

  1. Does this annotation look biologically reasonable? → `stats`
  2. How does it compare to related genomes? → `compare`
  3. How do I get my output into the tool my collaborator uses? → `convert`

This lesson covers all three. MycoNote-CLI deliberately does NOT
ship its own genome browser or circos plotter — that's a job for
best-in-class external tools (IGV, Proksee, JBrowse2, clinker)."
            ),

            LessonItem::Text(
"GENOME STATISTICS (myconote-cli stats)

The first question after annotation is: 'Does this look right?'
The `stats` command calculates:

  - Gene count, transcript count, CDS/exon counts
  - Gene length statistics (mean, median, N50)
  - GC content per contig
  - Exons per gene distribution
  - Per-chromosome breakdown
  - Optional benchmarking against taxon-specific expected ranges

With `--taxon fungi`, it warns when your numbers fall outside the
range typical for fungi:
  - 5,000-12,000 genes
  - Mean gene length of 1,200-1,800 bp
  - 1-3 exons per gene
  - 45-55% GC content

If your numbers are wildly different, something is likely wrong
with your assembly or prediction."
            ),

            LessonItem::CodeExample(
                "Generate statistics with taxonomic benchmarking:",
                "myconote-cli stats annotated.gff3 --taxon fungi"
            ),

            LessonItem::Question {
                prompt: "Which flag makes `stats` warn when your gene count looks unusual?",
                kind: QuestionKind::FreeText {
                    answer: "--taxon",
                    accept_regex: Some(r"(?i)--taxon"),
                },
                hint: Some("It compares your numbers against an expected range for a taxonomic group."),
                explanation: Some("`--taxon fungi` enables expected-range benchmarking. Other accepted values include ascomycota, basidiomycota."),
            },

            LessonItem::Text(
"N-GENOME ORTHOLOG INFERENCE (myconote-cli compare)

Once you have predictions for two or more genomes, you usually want
to ask: which genes are shared? which are lineage-specific? what's
the pan-genome shape?

`compare` wraps OrthoFinder (Emms & Kelly 2019, Genome Biol 20:238)
and emits:
  - orthogroups.tsv      Per-orthogroup gene IDs across genomes
  - per_genome_stats.tsv Total genes, orthogroup membership, singletons
  - pan_genome.tsv       Core / accessory / private gene counts
  - species_tree.nwk     Newick rooted species tree

With `--html`, all of the above plus an inline-SVG species tree are
bundled into a single self-contained `report.html` that opens
correctly from `file://` with no network — useful for sharing
with collaborators who don't run command-line tools."
            ),

            LessonItem::CodeExample(
                "Compare four fungal genomes and emit an HTML report:",
                "myconote-cli compare \\\n  genome_A.gff3 genome_B.gff3 genome_C.gff3 genome_D.gff3 \\\n  --fastas A.fa,B.fa,C.fa,D.fa \\\n  --output compare_out --html"
            ),

            LessonItem::Question {
                prompt: "What does `compare --html` produce that the plain TSV outputs do not?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "A high-resolution PDF for publication".into(),
                        "A self-contained interactive HTML report (sortable tables + inline-SVG tree)".into(),
                        "A circular Circos plot".into(),
                        "A JBrowse2 instance".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("Self-contained means no CDN, no network — opens directly from a file path."),
                explanation: Some("`--html` emits `compare_out/report.html` with CSS and JS inlined via `include_str!`, sortable per-genome stats and ortholog tables, and an inline-SVG species tree with branch lengths to scale."),
            },

            LessonItem::Text(
"PHYLOGENY: USE EXTERNAL TOOLS

MycoNote-CLI does NOT ship its own phylogenetic tree builder. The
upstream `phylogeny` subcommand was retired in v0.2.0. If you need
a maximum-likelihood species tree:

  1. Run `compare` to produce single-copy orthogroups.
  2. Find the directory `compare_out/Single_Copy_Orthologue_Sequences/`.
  3. Align each orthogroup with MAFFT (`mafft --auto`).
  4. Concatenate alignments into a supermatrix (e.g., FASconCAT-G).
  5. Run IQ-TREE 2 or RAxML-NG with appropriate partition models.

This is intentional separation of concerns: phylogenetic tree
inference is a deep field with its own tooling, and we'd rather
hand off to specialists than maintain a half-baked wrapper."
            ),

            LessonItem::Text(
"FORMAT CONVERSION (myconote-cli convert)

Bioinformatics uses dozens of file formats. `convert` handles:

  GFF3 conversions:
    --to gtf       Gene Transfer Format (Ensembl/GENCODE tools)
    --to bed       Browser Extensible Data (UCSC, bedtools)
    --to genbank   GenBank flat file (Geneious, SnapGene, Proksee)
    --to protein   Translated protein FASTA (BLAST, OrthoFinder)
    --to cds       Spliced nucleotide CDS (transcriptome for `quant`)

  FASTA conversions:
    --to fastq     FASTA -> FASTQ with dummy quality scores
    --to phylip    PHYLIP format (RAxML, PhyML)
    --to nexus     NEXUS format (MrBayes, PAUP*)

  VCF conversions:
    --to table     Tabular variant summary
    --to consensus Consensus FASTA with variants applied
    --to maf       Mutation Annotation Format

`convert --to cds` is the bridge into the RNA-seq stack: it
extracts the spliced nucleotide CDS that `quant` and `ase` use as
the salmon target transcriptome."
            ),

            LessonItem::Question {
                prompt: "How would you convert a GFF3 to GenBank format?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli convert genes.gff3 --to ___ --fasta genome.fa",
                    answer: "genbank",
                },
                hint: Some("The target format is GenBank (.gbk)."),
                explanation: Some("`--to genbank` produces a GenBank flat file compatible with Geneious, SnapGene, and Proksee. Requires `--fasta` because GenBank embeds the sequence."),
            },

            LessonItem::Text(
"VISUALIZATION: HAND OFF TO BEST-IN-CLASS TOOLS

MycoNote-CLI's strategy is convert-then-handoff:

  • IGV (desktop)
    Drop your GFF3 + genome FASTA in. Free, fast, works offline.
    Best for interactive browsing during analysis.

  • Proksee (web, https://proksee.ca/)
    Upload the `.gbk` from `convert --to genbank`. Excellent for
    publication-quality circular genome maps. No install needed.

  • JBrowse2
    Self-contained HTML browser. Use the published JSON converters
    to ingest your GFF3.

  • clinker (pip install clinker)
    Cross-species gene-cluster synteny diagrams from GenBank files.
    Great for comparing biosynthetic gene clusters across strains.

This separation lets you pick the best tool for the question
instead of being locked into a built-in viewer that does each
job 80%."
            ),

            LessonItem::Question {
                prompt: "Which external tool would you use to draw publication-quality circular genome maps from a `.gbk` produced by `convert --to genbank`?",
                kind: QuestionKind::FreeText {
                    answer: "Proksee",
                    accept_regex: Some(r"(?i)proksee"),
                },
                hint: Some("It's a free web tool; the name starts with 'P'."),
                explanation: Some("Proksee (https://proksee.ca/) is web-based, free, and produces publication-grade circular genome maps from GenBank files. No install required."),
            },

            LessonItem::Checkpoint("Analysis and handoff mastered"),

            LessonItem::Text(
"WHERE TO GO NEXT

Lessons 9-20 cover the modern features added since v0.5.0:

  9.  RNA-seq quantification (`quant`)
  10. SRA/ENA download (`fetch-rna`)
  11. Differential expression (`de-template`)
  12. Allele-specific expression (`ase`)
  13. ASE binomial test (`ase-template`)
  14. GO enrichment (`go-template`)
  15. BRAKER for maximum accuracy (`predict --use-braker`)
  16. GeneMark variants (`predict --genemark-mode`)
  17. Contig dedup (`clean --mode contigs`)
  18. Augustus fungal species bundle (`setup --db augustus-fungi`)
  19. Reproducibility manifests (`quant_bundle.json`, `ase_bundle.json`)
  20. NCBI codon tables (Candida CTG focus)

Each lesson is self-contained and assumes only the prerequisites
introduced in lessons 1-7. Skip ahead freely — the catalogue is
designed for non-linear use.

KEY REFERENCES:
  Stanke et al. (2006) BMC Bioinformatics 7:62 (Augustus)
  Steinegger & Soding (2017) Nature Biotech 35:1026 (MMseqs2)
  Manni et al. (2021) Mol Biol Evol 38:4647 (BUSCO)
  Haas et al. (2008) Genome Biol 9:R7 (EvidenceModeler)
  Emms & Kelly (2019) Genome Biol 20:238 (OrthoFinder)"
            ),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 9: RNA-seq Quantification (`quant`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_9_quant() -> Lesson {
    Lesson {
        title: "RNA-seq Quantification (quant)",
        description: "Salmon-based transcript quantification with fastp QC, decoy-aware indexing, tximport-ready outputs, and SHA256-keyed cache.",
        est_minutes: 12,
        items: vec![
            LessonItem::Text(
"FROM ANNOTATION TO EXPRESSION

Once you have a GFF3 annotation, the next question is usually:
'Which of these genes are actually expressed, and at what level?'

The `quant` subcommand runs the standard modern fungal RNA-seq
quantification recipe:

  raw FASTQs → fastp QC/trim → salmon (decoy-aware index) → counts + TPM

Outputs are tximport-ready, so you can drop them straight into
DESeq2 / edgeR / limma in R.

This lesson assumes you already have:
  - An annotated genome (GFF3 + FASTA from `annotate`)
  - Some RNA-seq FASTQs (paired or single-end)
  - salmon ≥ 1.10 and fastp ≥ 0.23 on PATH"
            ),

            LessonItem::Info("Check your tools first",
"Before starting, verify your environment:

  myconote-cli check
    Reports whether salmon and fastp are installed.

  myconote-cli install salmon fastp
    Installs them via conda/mamba if missing.

If you can't install external tools right now, you can still
complete this lesson — type 'skip' on the hands-on exercises and
read the explanations."
            ),

            LessonItem::Text(
"PREREQUISITE: A SPLICED-CDS TRANSCRIPTOME

salmon doesn't quantify against the genome FASTA — it quantifies
against a TRANSCRIPTOME (one sequence per transcript). For
eukaryotic genomes with introns, this means you need the spliced
nucleotide CDS, not the raw genome.

`convert --to cds` does exactly this: it walks each mRNA in the
GFF3, concatenates its child CDS features in genomic order,
reverse-complements on the `-` strand, and respects the phase
offset on the first CDS.

  myconote-cli convert annotated.gff3 --to cds \\
    --fasta genome.fa --output transcripts.fa

The resulting `transcripts.fa` is what `quant --transcripts` wants."
            ),

            LessonItem::Question {
                prompt: "What does `convert --to cds` produce that `convert --to protein` does not?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "A protein FASTA with stop codons stripped".into(),
                        "A spliced nucleotide FASTA (one entry per mRNA, suitable for salmon)".into(),
                        "A GenBank flat file with sequence + features".into(),
                        "A BED file of CDS coordinates".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("salmon quantifies nucleotides, not amino acids."),
                explanation: Some("`convert --to cds` emits spliced nucleotide CDS — concatenated child CDS features per mRNA, reverse-complemented on `-` strand, phase-respecting. `convert --to protein` translates to amino acids, which salmon can't use."),
            },

            LessonItem::Text(
"SAMPLE SHEET FORMAT

`quant` is sample-sheet-driven. The sheet is a TSV with at minimum:

  sample_id    fastq_r1                 fastq_r2          condition
  control_rep1 reads/control_R1.fq.gz   reads/control_R2.fq.gz   control
  control_rep2 reads/control_R2_R1.fq.gz reads/control_R2_R2.fq.gz control
  treated_rep1 reads/treated_R1.fq.gz   reads/treated_R2.fq.gz   treated
  treated_rep2 reads/treated_R2_R1.fq.gz reads/treated_R2_R2.fq.gz treated

Optional columns: `strandedness` (per-sample, default unstranded),
`batch`, and any others you want — unknown columns are preserved
verbatim and propagate to the output bundle.

Single-end? Leave `fastq_r2` empty for that row."
            ),

            LessonItem::CodeExample(
                "Quantify against your annotated genome:",
                "myconote-cli quant samples.tsv \\\n  --transcripts transcripts.fa \\\n  --decoys genome.fa \\\n  --output quant_out --threads 8"
            ),

            LessonItem::Text(
"DECOY-AWARE INDEXING (WHY YOU SET --decoys)

salmon's selective alignment can mistakenly map reads to
transcripts when the read's true source is intergenic / intronic
genomic DNA. The fix is to include the entire genome FASTA as
'decoy' sequence during index construction; salmon then assigns
reads to the genome instead of forcing them onto a transcript.

In MycoNote-CLI, `--decoys genome.fa` is the recommended default.
Skipping it is a footgun.

Ref: Srivastava et al. (2020) Genome Biol 21:239"
            ),

            LessonItem::Text(
"THE SHA256 INDEX CACHE

Building a salmon index takes minutes; it is wasteful to rebuild
when nothing has changed. `quant` keys its salmon index by
SHA256(transcripts.fa + decoys.fa + k-mer length). Identical
inputs → cache hit, no rebuild.

Cache directory precedence:
  --index-cache <dir> > $MYCONOTE_INDEX_CACHE > $XDG_CACHE_HOME
    > ~/.cache/myconote/ > ~/.myconote/

There is no silent fallback to /tmp — that would break bundle
provenance."
            ),

            LessonItem::Question {
                prompt: "What is the salmon index keyed on for the cache?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Filename + modification time".into(),
                        "SHA256 of transcripts.fa + decoys.fa + k-mer length".into(),
                        "The user's home directory".into(),
                        "The current date".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("Hash-based caching means content-addressed."),
                explanation: Some("SHA256 hashing means identical inputs are detected even if the file path or timestamp changes. This is what makes reruns deterministic and reproducible across machines."),
            },

            LessonItem::Text(
"OUTPUTS

After `quant` finishes you'll have:

  quant_out/
    counts.tsv             Wide matrix: transcript × sample → raw count
    tpm.tsv                Wide matrix: transcript × sample → TPM
    salmon/<sample>/quant.sf  tximport-native per-sample output
    fastp/<sample>.json    Per-sample QC report (Q30, dup rate, …)
    quant_bundle.json      Reproducibility manifest (covered in lesson 19)

For DE analysis, drop the `salmon/` directory into R via tximport
and pass to DESeq2 — see lesson 11 (`de-template`)."
            ),

            LessonItem::TryIt {
                instruction: "Write the command to extract a CDS transcriptome from your annotation, then quantify a sample sheet against it (4 threads):",
                command_template: "myconote-cli convert annotated.gff3 --to ___ --fasta genome.fa --output transcripts.fa\nmyconote-cli quant samples.tsv --transcripts transcripts.fa --decoys ___ --output quant_out --threads ___",
                hint: Some("The convert format is `cds`, the decoy is `genome.fa`, threads is `4`."),
            },

            LessonItem::Checkpoint("Quantification pipeline understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 10: SRA/ENA Download (`fetch-rna`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_10_fetch_rna() -> Lesson {
    Lesson {
        title: "SRA/ENA Download (fetch-rna)",
        description: "Download public RNA-seq FASTQs by accession with MD5 verification, ENA-default-with-sra-fallback, and quant-ready sample sheet emission.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"WHERE DOES PUBLIC RNA-SEQ LIVE?

If you're not generating your own RNA-seq, you're downloading it
from public archives. The two main mirrors:

  SRA   (NCBI)  — requires sra-toolkit (`prefetch` + `fasterq-dump`)
                  vdb-config + cloud-credential setup, slow, fiddly.

  ENA   (EMBL)  — REST API, no credentials, MD5-verifiable URLs.
                  This is MycoNote-CLI's default backend.

`fetch-rna` tries ENA first; if a run isn't mirrored there, falls
back to sra-toolkit (only if `--allow-sra` is set). It accepts:
  - Run IDs       (SRR1234567, ERR1234567)
  - Study IDs     (PRJNA123456, PRJEB123456)
  - Sample IDs    (SAMN12345678)
  - A text file with one accession per line"
            ),

            LessonItem::Info("Why ENA by default?",
"sra-toolkit is the official NCBI route, but in practice it has
real friction: vdb-config dialogs, cloud credentials, .ncbi/ dirs,
and slow `prefetch` against intermittent CDN endpoints.

ENA mirrors the same data with a stable REST API (`fastq_ftp` URL
field), per-file MD5s for verification, and no credentials. For
the 90% case where the run is mirrored to ENA, this is faster
and more reliable. The sra-toolkit fallback exists for the 10%
of runs ENA hasn't ingested yet."
            ),

            LessonItem::CodeExample(
                "Download an entire BioProject's RNA-seq, emit a quant-ready sample sheet:",
                "myconote-cli fetch-rna PRJNA123456 \\\n  --output rnaseq/ \\\n  --emit-sample-sheet samples.tsv \\\n  --threads 4"
            ),

            LessonItem::Text(
"MD5 VERIFICATION

ENA exposes per-file MD5 checksums via its REST API. `fetch-rna`
streams each FASTQ while computing MD5 on the fly and verifies
against the expected hash before declaring success. Mismatches
trigger one retry (transient network corruption is more common
than you'd think); a second mismatch is a hard error.

This is non-negotiable: a corrupted FASTQ can pass through your
entire pipeline and produce silently wrong DE results."
            ),

            LessonItem::Question {
                prompt: "What does `fetch-rna` do when the MD5 of a downloaded FASTQ doesn't match ENA's expected hash?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Silently keeps the file and continues".into(),
                        "Retries once, then fails loudly with the file path and expected/actual MD5s".into(),
                        "Switches to sra-toolkit automatically".into(),
                        "Skips that sample and moves on".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("Silent corruption is the worst possible outcome — fail loudly."),
                explanation: Some("One retry covers transient network corruption; a second failure is a hard error. The pipeline must not continue with corrupted data."),
            },

            LessonItem::Text(
"SAMPLE SHEET EMISSION

`--emit-sample-sheet samples.tsv` writes a TSV pre-populated for
`quant`:

  sample_id  fastq_r1                fastq_r2                condition
  SRR12345   rnaseq/SRR12345_1.fq.gz rnaseq/SRR12345_2.fq.gz unknown
  SRR12346   rnaseq/SRR12346_1.fq.gz rnaseq/SRR12346_2.fq.gz unknown

The `condition` column is filled with `unknown` — you edit the
file by hand to assign biological condition labels before passing
it to `quant`. This deliberate two-step process forces you to
read the SRA metadata and assign conditions correctly, instead of
silently treating all samples as one group."
            ),

            LessonItem::TryIt {
                instruction: "Write the command to download a single SRA run (SRR12345678) into ./rnaseq, emitting a sample sheet:",
                command_template: "myconote-cli fetch-rna ___ --output rnaseq/ --emit-sample-sheet samples.tsv",
                hint: Some("The accession goes immediately after the subcommand name."),
            },

            LessonItem::Checkpoint("Public-data ingestion mastered"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 11: Differential Expression (`de-template`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_11_de_template() -> Lesson {
    Lesson {
        title: "Differential Expression (de-template)",
        description: "Generate a self-contained DESeq2 R script from quant output: tximport, design formula, contrasts, apeglm shrinkage, MA + volcano plots.",
        est_minutes: 12,
        items: vec![
            LessonItem::Text(
"WHY A TEMPLATE INSTEAD OF AN R DEPENDENCY?

DESeq2 (Love et al. 2014) is the gold-standard fungal DE tool, but
running it requires R + Bioconductor. We deliberately do NOT take
R as a runtime dependency on MycoNote-CLI — that would explode the
install footprint and leak version conflicts into the Rust build.

Instead, `de-template` is Option 1D from the RNA-seq spec: it
emits a self-contained R script that:

  1. Loads tximport + DESeq2 + apeglm
  2. Reads `quant_out/salmon/<sample>/quant.sf`
  3. Builds the design formula and runs DESeq2
  4. Applies apeglm LFC shrinkage (with graceful fallback if missing)
  5. Writes per-contrast TSVs sorted by padj
  6. Renders MA + volcano plots with FDR + |LFC| threshold lines

You run the script with `Rscript de_analysis.R`. We write the
boilerplate; you provide R."
            ),

            LessonItem::Info("Prerequisite alerting in three layers",
"To avoid the 'tool wrote a script that won't run' failure mode:

  1. CLI `--help` lists R + Bioconductor packages required.
  2. Runtime stderr warning when `Rscript` isn't on PATH at all.
  3. The emitted R script itself starts with a `requireNamespace()`
     guard that prints a clear `BiocManager::install(...)` hint
     if any package is missing.

You'll always know what to install, regardless of which step
catches the gap first."
            ),

            LessonItem::Text(
"DESIGN FORMULAS AND CONTRASTS

`de-template` validates two things upfront before writing anything:

  1. Factor column existence
     `--factor condition` → must be a column in samples.tsv.

  2. Design-formula sanitization
     The string passed to `--design` is checked against `;`,
     backticks, newlines, `system(...)`, `eval(...)`. We reject
     anything that smells like an R injection — the user runs
     this script with their own privileges.

A simple two-condition contrast:
  --design '~ condition' --contrast condition,treated,control

A batch-corrected design:
  --design '~ batch + condition' --contrast condition,treated,control

Multiple contrasts in one invocation: repeat `--contrast`."
            ),

            LessonItem::CodeExample(
                "Emit a DE script for treated-vs-control with batch correction:",
                "myconote-cli de-template quant_out/ \\\n  --samples quant_out/samples.tsv \\\n  --design '~ batch + condition' \\\n  --contrast condition,treated,control \\\n  --output de_out/"
            ),

            LessonItem::Question {
                prompt: "Why does `de-template` reject a `--design` value containing a semicolon?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "DESeq2 doesn't accept semicolons".into(),
                        "Sanitizing prevents R injection in the emitted script".into(),
                        "It would crash R on Windows".into(),
                        "It's a syntax error in tximport".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("The user runs the emitted script with their own privileges."),
                explanation: Some("The string is interpolated into the R script we write. Tokens like `;`, backticks, `system(...)`, and `eval(...)` could let a malformed sample sheet inject arbitrary R code that runs as the user."),
            },

            LessonItem::Text(
"OUTPUT FILES (after Rscript)

  de_out/
    de_analysis.R                                 The script you Rscript
    de_<factor>_<level1>_vs_<level2>.tsv          DE results, sorted by padj
    de_<factor>_<level1>_vs_<level2>_MA.png       MA plot (mean expr vs LFC)
    de_<factor>_<level1>_vs_<level2>_volcano.png  Volcano (LFC vs -log10 padj)

Each plot includes horizontal/vertical threshold lines at the
default FDR (0.05) and |LFC| ≥ 1, configurable via flags. The
generated R script is fully readable — you can edit it after
emission to tweak plot aesthetics or add a custom contrast."
            ),

            LessonItem::Question {
                prompt: "After `de-template` writes the script, how do you actually run the analysis?",
                kind: QuestionKind::FreeText {
                    answer: "Rscript",
                    accept_regex: Some(r"(?i)rscript"),
                },
                hint: Some("It's a standard R command-line invocation."),
                explanation: Some("`Rscript de_out/de_analysis.R` runs the analysis. The tool stops at writing the script — execution is delegated to R so you control which R/Bioconductor versions are used."),
            },

            LessonItem::Checkpoint("DE templating understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 12: Allele-Specific Expression (`ase`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_12_ase() -> Lesson {
    Lesson {
        title: "Allele-Specific Expression (ase)",
        description: "Personalized per-haplotype salmon quantification on phased VCFs. Heterozygous diploid Candida, Saccharomyces F1 hybrids, lager allopolyploids.",
        est_minutes: 14,
        items: vec![
            LessonItem::Text(
"WHEN ALLELES MATTER

In a heterozygous diploid (Candida albicans), an F1 hybrid (yeast
crosses), or an allopolyploid (lager yeast), the two homologous
chromosomes carry different alleles. Standard RNA-seq
quantification merges reads from both alleles and reports the
SUM — losing the per-allele signal entirely.

ASE asks the per-allele question: 'For this transcript, are reads
distributed 70/30 across the two alleles, or 50/50?' Imbalance can
indicate cis-regulatory variation, imprinting, or selection.

`ase` is MycoNote-CLI's transcript-level ASE engine — not a
read-level WASP/STAR retrofit, but a personalized-transcriptome
salmon approach designed for the fungal use case."
            ),

            LessonItem::Text(
"THE PERSONALIZED-TRANSCRIPTOME APPROACH

Build a separate transcriptome per haplotype, then quantify each
sample against each haplotype with salmon:

  1. Read a phased VCF
  2. For each variant, apply ALT allele to the haplotype it's on.
  3. Build cds_hap0.fa and cds_hap1.fa per-haplotype CDS FASTAs.
  4. For each sample, run salmon × hap0, then salmon × hap1.
  5. Merge per-sample-per-haplotype counts into a transcript ×
     (sample.haplotype) wide matrix.
  6. Compute per-transcript informativeness (does it carry any
     phased differences between hap0 and hap1?).
  7. Detect mapping-rate asymmetry between haplotypes."
            ),

            LessonItem::Info("THE UNPHASED-HET HARD ERROR",
"`ase` REFUSES to operate on unphased heterozygous sites.

If you give it a VCF with a `0/1` (unphased) genotype, it errors
with the line number and a pointer at WhatsHap or HapCUT2 for
phasing. This is by design.

Why? An unphased het has unknown allele assignment — we don't
know which homologous chromosome carries the ALT. Picking
arbitrarily, or 'splitting evenly', would silently corrupt every
downstream count. Better to fail loudly than produce wrong numbers.

If you have only unphased calls, run a phasing tool first.
WhatsHap (Patterson et al. 2015) is the standard fungal choice."
            ),

            LessonItem::Question {
                prompt: "What does `ase` do when it encounters a `0/1` (unphased het) genotype in your VCF?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Silently picks one haplotype at random".into(),
                        "Splits the variant evenly between haplotypes".into(),
                        "Hard-errors with the line number and a pointer at WhatsHap".into(),
                        "Skips the variant and continues".into(),
                    ],
                    correct_index: 2,
                },
                hint: Some("Silent corruption is the worst possible outcome."),
                explanation: Some("Hard error. Unphased hets have unknown allele assignment, and picking arbitrarily would silently corrupt every downstream count. WhatsHap or HapCUT2 are the standard phasing tools."),
            },

            LessonItem::Text(
"INDEL SUPPORT AND SKIP CATEGORIES

`ase` handles SNVs, MNPs, and indels up to `--max-indel-size`
(default 50 bp). The personalization layer documents nine skip
categories so every variant lands on one side of the
applied/skipped divide with a reason:

  - exon_boundary_spanning   indel crosses an exon/intron junction
  - in_cis_overlap           two phased variants overlap on same hap
  - too_large                indel size > --max-indel-size
  - off_cds                  variant outside any CDS
  - non_acgt                 ALT contains ambiguity codes
  - … etc.

Both audit trails are emitted as `variants_applied.tsv` and
`variants_skipped.tsv` so you can review every decision."
            ),

            LessonItem::CodeExample(
                "Run ASE on a heterozygous Candida diploid:",
                "myconote-cli ase samples.tsv \\\n  --vcf phased.vcf.gz \\\n  --gff3 annotated.gff3 \\\n  --fasta genome.fa \\\n  --output ase_out --threads 8"
            ),

            LessonItem::Text(
"MAPPING-RATE ASYMMETRY WARNINGS

A subtle failure mode: if hap0's salmon mapping rate is 92% but
hap1's is 78% for the same sample, every count comparison is
biased — it looks like hap0 is over-expressed because more reads
mapped overall.

`ase` detects this with a configurable threshold (default 5
percentage points) and writes a warning to stderr plus a flag in
`ase_bundle.json`. It does NOT silently correct — that would mask
real biology. Investigate: it's usually a phasing error, a
contig-naming mismatch, or a contaminating second strain."
            ),

            LessonItem::Question {
                prompt: "What's the default mapping-rate-asymmetry threshold above which `ase` warns?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "1 percentage point".into(),
                        "5 percentage points".into(),
                        "10 percentage points".into(),
                        "There is no default — you must set it".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It's a small but meaningful number."),
                explanation: Some("5 percentage points is the default. A 5pp gap between hap0 and hap1 is small enough to occur from real biological imbalance but large enough to flag as worth investigating."),
            },

            LessonItem::Checkpoint("ASE quantification understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 13: ASE Binomial Test (`ase-template`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_13_ase_template() -> Lesson {
    Lesson {
        title: "ASE Binomial Test (ase-template)",
        description: "Per-sample-corrected binomial.test on ase output: sample-specific null, informativeness filter, base-R only (no Bioconductor).",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"FROM PER-HAPLOTYPE COUNTS TO SIGNIFICANCE

`ase` gives you per-transcript-per-haplotype counts. The next
question is: which transcripts are significantly imbalanced
between hap0 and hap1?

The standard fungal answer is the binomial exact test. For a
transcript with n0 reads on hap0 and n1 reads on hap1, ask: is
that ratio significantly different from the expected null?

`ase-template` emits a base-R (no Bioconductor) script that does
this with the right correction."
            ),

            LessonItem::Text(
"THE NULL-CORRECTION INSIGHT

Naive ASE testing uses `binom.test(c(n0, n1), p=0.5)` — testing
each transcript against a 50/50 null. But this is wrong: if the
sample's overall mapping rate is 92% on hap0 and 78% on hap1,
the null is NOT 50/50 — it's the per-sample library ratio.

`ase-template` computes a sample-specific null:

  null_p = total_hap0_counts / (total_hap0_counts + total_hap1_counts)

then tests each transcript against that null with binom.test:

  binom.test(c(n_hap0, n_hap1), p = null_p)

This is the most pedagogically valuable concept in the ASE stack:
testing against a 50/50 null when the sample-wide ratio is 60/40
will inflate your false-positive rate by ~10-20x."
            ),

            LessonItem::Info("Why is this concept the most important?",
"Because every undergraduate ASE tutorial uses 50/50 as the null.
Every. Single. One.

If you're working in a real fungal system — a heterozygous yeast
isolate, a hybrid, an allopolyploid — the 50/50 null is wrong.
Mapping bias, copy-number differences, and library prep biases
all push the actual null away from 50/50.

The fix is simple: use the sample's own total hap0:hap1 ratio
as the null. `ase-template` does this automatically. Without
this correction, you will have false positives that don't
replicate."
            ),

            LessonItem::Question {
                prompt: "What null proportion does `ase-template` use for the binom.test of each transcript?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Always 0.5 (50/50)".into(),
                        "The sample-specific total hap0:hap1 library ratio".into(),
                        "The genome-wide hap0:hap1 ratio".into(),
                        "A user-supplied constant".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("The null should reflect the actual library composition for that sample."),
                explanation: Some("Sample-specific null. If sample A maps 92% on hap0 and 78% on hap1, the per-sample null is 92/(92+78). Using 50/50 inflates false positives by 10-20x."),
            },

            LessonItem::Text(
"INFORMATIVENESS FILTERING

A transcript with no phased differences between hap0 and hap1
is non-informative for ASE — there's nothing to distinguish hap0
reads from hap1 reads. By default, `ase-template` filters to
informative transcripts only (set during `ase`). Override with
`--include-uninformative` if you want the full table.

Low-count rows (default `--min-reads 20` total) are also skipped
with a logged reason — binom.test on n=4 is noise."
            ),

            LessonItem::CodeExample(
                "Emit and run a per-sample-corrected ASE binomial test script:",
                "myconote-cli ase-template ase_out/ --output ase_de/\nRscript ase_de/ase_binomial.R"
            ),

            LessonItem::Text(
"OUTPUT FILES

  ase_de/
    ase_binomial.R                 The script you Rscript
    ase_results.tsv                Long-format: sample, transcript, n0, n1, p, padj
    <sample>_imbalance_hist.pdf    Per-sample log2(hap0/hap1) histogram

Each sample gets its own histogram so you can eyeball whether
ASE imbalance is a long-tail (a few real ASE genes) or
broadly-distributed (likely a phasing or mapping artefact)."
            ),

            LessonItem::Checkpoint("ASE binomial testing understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 14: GO Enrichment (`go-template`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_14_go_template() -> Lesson {
    Lesson {
        title: "GO Enrichment (go-template)",
        description: "topGO Fisher's exact test on DE results joined to annotate GO terms. BP/MF/CC ontologies, BH-adjusted within ontology, dot-plot output.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"FROM DE GENES TO BIOLOGICAL THEMES

After DE analysis you have a list of genes with adjusted p-values.
The biological question is: are those significant genes enriched
for any particular function — DNA repair, oxidative stress
response, cell-wall remodeling?

GO (Gene Ontology) enrichment answers this by testing whether your
DE-significant set carries more genes annotated with a given GO
term than expected by chance.

`go-template` joins three things:

  1. DE results (from `de-template`) — which genes are significant?
  2. Annotation (from `annotate`) — which GO terms does each gene have?
  3. The GO graph — for each term, count its descendants too.

Then runs topGO's classic Fisher's exact test."
            ),

            LessonItem::Text(
"THREE ONTOLOGIES, INDEPENDENTLY

GO has three orthogonal ontologies:

  BP   Biological Process    'response to oxidative stress'
  MF   Molecular Function    'ATP binding'
  CC   Cellular Component    'mitochondrial inner membrane'

Each is tested independently and BH-adjusted within ontology
(not across) — that's the convention because the three are
biologically distinct hierarchies. By default `go-template` runs
all three; pass `--ontology BP` to run just one.

Ref: Alexa et al. (2006) Bioinformatics 22:1600 (topGO)"
            ),

            LessonItem::Text(
"FLEXIBLE JOINING

DE results are usually keyed by transcript ID; `annotate` output
is usually keyed by gene/locus ID. `go-template` exposes both
join columns as flags:

  --de-id-col transcript_id
  --ann-id-col locus_tag

The GO column itself defaults to `go_terms` (matches `annotate`'s
output) but can be overridden. The separator within the GO cell
defaults to `,` but can be `;`, `|`, etc. depending on what your
annotation pipeline emits.

This means the same template handles transcript-level DE joined
to locus-tag annotations, gene-level DE joined to gene-id
annotations, or any other shared identifier."
            ),

            LessonItem::CodeExample(
                "Emit a GO-enrichment script for treated-vs-control DE:",
                "myconote-cli go-template \\\n  --de-results de_out/de_condition_treated_vs_control.tsv \\\n  --annotations annotate_out/annotations.tsv \\\n  --output go_out/\nRscript go_out/go_enrichment.R"
            ),

            LessonItem::Question {
                prompt: "How are p-values multiple-testing-corrected across the three GO ontologies?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "BH-adjusted across all three pooled together".into(),
                        "BH-adjusted within each ontology independently".into(),
                        "Bonferroni across all three pooled together".into(),
                        "No correction is applied".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("BP, MF, and CC are biologically independent hierarchies."),
                explanation: Some("Within each ontology independently. The three ontologies are biologically distinct; pooling for BH would mask BP signal under MF noise (or vice versa)."),
            },

            LessonItem::Text(
"OUTPUTS

  go_out/
    go_enrichment.R              The script you Rscript
    go_BP.tsv                    BP results, sorted by classic Fisher p
    go_MF.tsv
    go_CC.tsv
    go_combined_dotplot.pdf      Top-N terms across all 3 ontologies

The dot plot shows enrichment fold-change on x-axis, term on
y-axis, dot size = number of significant genes in the term, and
dot color = -log10(adjusted p). Standard publication-grade
output."
            ),

            LessonItem::Checkpoint("GO enrichment understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 15: BRAKER for Maximum Accuracy (`predict --use-braker`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_15_braker() -> Lesson {
    Lesson {
        title: "BRAKER for Maximum Accuracy",
        description: "When and how to bypass MycoNote's multi-tool consensus and run BRAKER 1/2/3 end-to-end via predict --use-braker. Mode auto-detection.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"BRAKER VS MYCONOTE'S DEFAULT CONSENSUS

The default `predict` runs Augustus + SNAP + GlimmerHMM + GeneMark
+ miniprot, then merges them through Evidence Modeler. That's a
robust consensus for organisms without extensive evidence.

But for organisms WITH extensive RNA-seq + protein databases, the
state of the art is BRAKER (Hoff et al. 2019, NAR Genom Bioinform
1:lqaa108): a tightly-integrated Augustus + GeneMark training
pipeline that produces higher-accuracy annotations than running
either alone.

`predict --use-braker` invokes braker.pl end-to-end. MycoNote does
not run Augustus/SNAP/GlimmerHMM/GeneMark/miniprot itself or the
EVM consensus stage when this flag is set; `braker.gff3` becomes
`consensus.gff3` directly. Downstream `update`, `annotate`,
`submit` all consume it unchanged."
            ),

            LessonItem::Info("When to choose BRAKER over the default",
"Choose BRAKER when:
  - You have RNA-seq aligned to the genome (BAM file).
  - You have a protein database from a related organism (FASTA).
  - Your organism is well-represented in OrthoDB or you have
    your own curated protein set.
  - You're doing a high-stakes, publication-quality annotation
    where extra accuracy matters.

Stick with the default consensus when:
  - You don't have substantial RNA-seq or protein evidence.
  - Your organism is novel and has no close OrthoDB representation.
  - You want a fast, reasonable first pass.
  - You need GeneMark-EP+ or ETP+ specifically (lesson 16)."
            ),

            LessonItem::Text(
"MODE AUTO-DETECTION

BRAKER has three modes that map to MycoNote's BrakerMode enum:

  BRAKER1  RNA-seq only             braker.pl --esmode
  BRAKER2  Protein only             braker.pl --epmode
  BRAKER3  RNA-seq + Protein        braker.pl --etpmode

`predict --use-braker` auto-detects from inputs:
  --braker-rna-bam set, --braker-proteins unset  →  BRAKER1
  --braker-rna-bam unset, --braker-proteins set  →  BRAKER2
  Both set                                       →  BRAKER3

Override with `--braker-mode <1|2|3>` if you want explicit
control."
            ),

            LessonItem::Question {
                prompt: "If you pass `--braker-rna-bam reads.bam` AND `--braker-proteins proteins.fa`, which BRAKER mode does MycoNote-CLI auto-select?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "BRAKER1".into(),
                        "BRAKER2".into(),
                        "BRAKER3".into(),
                        "It errors because the modes conflict".into(),
                    ],
                    correct_index: 2,
                },
                hint: Some("Both kinds of evidence → mode 3."),
                explanation: Some("BRAKER3 (--etpmode) is the joint RNA-seq + protein mode. Override with `--braker-mode 1` or `--braker-mode 2` if you want only one source despite supplying both."),
            },

            LessonItem::Text(
"FLAG-CONFLICT GUARD

`predict --use-braker` is incompatible with the standard predictor
flags:

  --genemark-mode      → conflict
  --genemark           → conflict (legacy alias)
  --genemark-hints     → conflict (legacy alias)
  --protein-evidence   → conflict
  --protein-fasta      → conflict
  --glimmerhmm         → conflict

If you set `--use-braker` plus any of the above, `check_braker_conflicts`
rejects the call with a SINGLE error listing every offending flag,
so the fix is one round-trip not whack-a-mole.

Note: `--no-snap` is intentionally NOT a conflict (disabling a
predictor that wouldn't run anyway is a no-op)."
            ),

            LessonItem::Text(
"GENETIC-CODE FORWARDING

`--genetic-code <n>` forwards through to BRAKER's `--translation_table`,
wiring the *Candida* CTG code (12) and ciliate-style codes (6)
end-to-end through Augustus + GeneMark inside BRAKER. This is the
fungal-genomics-relevant difference from running BRAKER directly:
you set the code once on the MycoNote command line and it
propagates correctly into BRAKER's internal calls."
            ),

            LessonItem::CodeExample(
                "Run BRAKER3 for a Candida CTG-clade species:",
                "myconote-cli predict masked.fa \\\n  --use-braker \\\n  --braker-rna-bam reads.bam \\\n  --braker-proteins proteins.fa \\\n  --genetic-code 12 \\\n  --output predict_out --threads 16"
            ),

            LessonItem::TryIt {
                instruction: "Write the command to run BRAKER1 (RNA-seq only) on a masked genome with 8 threads:",
                command_template: "myconote-cli predict masked.fa --use-braker --braker-rna-bam reads.bam --threads ___",
                hint: Some("Threads is 8."),
            },

            LessonItem::Checkpoint("BRAKER usage understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 16: GeneMark Variants (`predict --genemark-mode`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_16_genemark() -> Lesson {
    Lesson {
        title: "GeneMark Variants (--genemark-mode)",
        description: "ES vs ET vs EP+ vs ETP+: when each is appropriate, what evidence each consumes, the ProtHint dependency.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"FOUR GENEMARK MODES, FOUR USE CASES

GeneMark (Lomsadze et al., Nucleic Acids Res 33:6494) is a
self-training gene predictor that doesn't require a pre-built
species model — it iteratively learns gene structure from the
genome itself. Four operational modes ship today:

  ES     Self-training only. Pure ab initio. The fallback mode
         when no evidence is available.

  ET     Self-training + RNA-seq intron hints (`--ET`). Use when
         you have RNA-seq aligned to the genome.

  EP+    Self-training + protein evidence via ProtHint
         (`gmes_petap.pl --EP`). Use when you have a protein
         database but no RNA-seq.

  ETP+   Self-training + RNA-seq intron hints + protein evidence
         (`--ETP`). The strongest GeneMark mode, when you have
         both kinds of evidence."
            ),

            LessonItem::Info("Why EP+ matters for novel CTG-clade fungi",
"For a novel Candida CTG-clade isolate without RNA-seq:
  - No close trained Augustus species exists.
  - The alternative-yeast-nuclear code (Table 12) makes most
    published predictors error-prone if you forget to set it.
  - Self-training (ES mode) without anchoring will overfit weird
    CTG-encoded patterns.

EP+ solves this: protein evidence anchors the training, and the
genetic-code flag propagates through ProtHint correctly. This is
the highest-impact use case for the EP+ mode."
            ),

            LessonItem::Text(
"PROTHINT — THE EP+/ETP+ DEPENDENCY

EP+ and ETP+ each call ProtHint internally to convert the
genome + protein FASTA into the GFF hint files that GeneMark
consumes.

ProtHint installation reality (verified 2026-04-24 against
bioconda osx-64 and noarch):
  - There is NO standalone bioconda recipe for `prothint`.
    `conda install prothint` fails with PackagesNotFoundError.
  - ProtHint ships bundled with the GeneMark-ES installer
    tarball under `<install>/ProtHint/bin/`.
  - It also ships with the bioconda `braker3` package.

`myconote-cli check` flags ProtHint with a `manual_note` field
that explains this; `myconote-cli install` prints the same
guidance instead of trying a bioconda install that will fail."
            ),

            LessonItem::Question {
                prompt: "Which `--genemark-mode` value would you choose for a novel Candida species with NO RNA-seq but a curated protein FASTA from a related yeast?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "es".into(),
                        "et".into(),
                        "ep".into(),
                        "etp".into(),
                    ],
                    correct_index: 2,
                },
                hint: Some("No RNA-seq rules out modes that need RNA-seq hints."),
                explanation: Some("`ep` (EP+). Protein-only evidence. ETP requires both RNA-seq AND protein. ET requires RNA-seq. ES would self-train without evidence which is the weakest of the four."),
            },

            LessonItem::Text(
"LEGACY-FLAG TRANSLATION

The legacy v0.5.x flags `--genemark` and `--genemark-hints` are
auto-translated for backward compatibility:

  --genemark (alone)           → mode ES
  --genemark + --genemark-hints → mode ET
  --genemark-hints alone       → mode ET

EP/ETP have NO legacy alias because they did not exist in earlier
releases. New scripts should use `--genemark-mode` explicitly.

The mode-string parser is case-insensitive and accepts the
publication-style aliases too: `EP+` and `ETP+` both work."
            ),

            LessonItem::CodeExample(
                "Run ETP+ with both RNA-seq hints and a protein database:",
                "myconote-cli predict masked.fa \\\n  --genemark-mode etp \\\n  --genemark-hints rnaseq_hints.gff \\\n  --protein-fasta orthodb_fungi.fa \\\n  --output predict_out --threads 16"
            ),

            LessonItem::Question {
                prompt: "What does `--genemark-mode etp` require that the other modes do not?",
                kind: QuestionKind::FreeText {
                    answer: "both",
                    accept_regex: Some(r"(?i)(both|requires.*hints.*protein|protein.*hints|genemark-hints.*protein-fasta)"),
                },
                hint: Some("ETP needs two kinds of evidence."),
                explanation: Some("ETP requires both `--genemark-hints` and `--protein-fasta`. If either is missing, the run aborts with a named error — no silent fallback to a weaker mode."),
            },

            LessonItem::Checkpoint("GeneMark modes understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 17: Contig Dedup (`clean --mode contigs`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_17_clean_contigs() -> Lesson {
    Lesson {
        title: "Contig Dedup (clean --mode contigs)",
        description: "purge_dups-style contig deduplication via minimap2 self-alignment. Coverage/identity thresholds, deterministic tie-breaking, drop-report TSV.",
        est_minutes: 8,
        items: vec![
            LessonItem::Text(
"WHY DRAFT ASSEMBLIES HAVE DUPLICATE CONTIGS

Long-read fungal assemblies are great for contiguity but often
emit haplotigs — duplicate contigs representing the alternative
allele of a heterozygous region. If you don't purge them, you'll
double-count genes downstream.

The classic answer is `purge_dups` (Guan et al. 2020,
Bioinformatics 36:2896): align the assembly against itself and
drop contigs whose entire span is covered by a longer contig at
high identity.

`clean --mode contigs` is MycoNote-CLI's port of this idea: a
minimap2-based self-alignment, configurable thresholds, and
deterministic tie-breaking."
            ),

            LessonItem::Text(
"HOW IT WORKS

  1. minimap2 -X self-alignment of the genome FASTA.
  2. Parse the PAF output.
  3. For each contig, compute the longest-contig coverage at
     high-identity hits.
  4. Drop the contig if its span is ≥ --coverage covered by a
     longer contig at ≥ --identity (defaults: 0.95 / 0.95).
  5. Emit a cleaned FASTA + a TSV report listing every dropped
     contig and the contig that subsumed it.

Mutual subsumption is broken deterministically: longer contig
wins, lex-smaller name wins on ties. No randomness. Same input
→ same output bytes."
            ),

            LessonItem::CodeExample(
                "Purge haplotigs from a draft yeast assembly:",
                "myconote-cli clean --mode contigs assembly.fa \\\n  --coverage 0.95 --identity 0.95 \\\n  --output cleaned.fa --report drops.tsv"
            ),

            LessonItem::Question {
                prompt: "How does `clean --mode contigs` break mutual-subsumption ties (contig A covers B AND B covers A)?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Both are dropped".into(),
                        "Longer contig wins; lex-smaller name wins on size ties".into(),
                        "Random choice with a logged seed".into(),
                        "User must intervene manually".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("Reproducibility means deterministic tie-breaking."),
                explanation: Some("Longer wins, lex-smaller on size ties. Deterministic — same input, same output, every time."),
            },

            LessonItem::Text(
"WHEN TO RUN IT

`clean --mode contigs` is a PRE-annotation step for diploid /
hybrid / heterozygous draft assemblies. Run it after `sort` and
before `mask`:

  sort  →  clean --mode contigs  →  mask  →  predict  →  …

For high-quality reference assemblies (already polished and
purged upstream), it's a no-op and safe to skip.

For ASE work (lesson 12), do NOT purge haplotigs — the personalized-
transcriptome approach DEPENDS on having both haplotypes
represented in the genome FASTA + phased VCF."
            ),

            LessonItem::Question {
                prompt: "True or False: For ASE analysis (lesson 12), you should run `clean --mode contigs` before `predict` to purge haplotigs.",
                kind: QuestionKind::TrueFalse { answer: false },
                hint: Some("ASE depends on both haplotypes being present."),
                explanation: Some("False. ASE uses a personalized-transcriptome approach that needs both haplotypes preserved — purging would erase the very signal you're trying to measure. Use `clean --mode contigs` only when you intend to collapse to a single haploid representation."),
            },

            LessonItem::Checkpoint("Contig dedup understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 18: Augustus Fungal Species Bundle (`setup --db augustus-fungi`)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_18_augustus_fungi() -> Lesson {
    Lesson {
        title: "Augustus Fungal Species Bundle",
        description: "What `setup --db augustus-fungi` ships, AUGUSTUS_CONFIG_PATH, when to retrain vs. use a bundled species.",
        est_minutes: 8,
        items: vec![
            LessonItem::Text(
"49 CURATED FUNGAL AUGUSTUS SPECIES

A fresh Augustus install ships only ~10 fungal species models
(the classic ones — saccharomyces_cerevisiae_S288C, aspergillus_fumigatus,
neurospora_crassa). Funannotate downloads more from elsewhere;
those URLs aren't always stable.

`setup --db augustus-fungi` ships a curated manifest of 49 fungal
Augustus species, fetched directly from the upstream
`Gaius-Augustus/Augustus/config/species/` directory and verified
against the GitHub API. Coverage spans:

  - Saccharomycotina (budding yeasts: S. cerevisiae, Candida,
    Kluyveromyces, Yarrowia, …)
  - Taphrinomycotina (S. pombe, Pneumocystis)
  - Pezizomycotina (Aspergillus, Fusarium, Neurospora,
    Magnaporthe, Botrytis, Sclerotinia, …)
  - Basidiomycota (Cryptococcus, Coprinus, Laccaria,
    Phanerochaete, Ustilago)
  - Mucoromycota / Microsporidia (Rhizopus, Encephalitozoon)

Plus the bare-genus aliases that Augustus historically ships
(`saccharomyces`, `fusarium`, `cryptococcus`, etc.)."
            ),

            LessonItem::Info("Why 49 and not the spec's hoped-for ~100?",
"The original RNA-seq spec hoped for ~100 entries to match
funannotate's bundled tree. But upstream
`Augustus/config/species/` only ships ~50 fungal entries —
funannotate downloads beyond what's in the upstream repo.

We don't fabricate names that aren't upstream. Every name in
the curated 49 is verified against the GitHub API. If you need
a species that isn't bundled, train your own with
`myconote-cli train` (lesson 6)."
            ),

            LessonItem::Text(
"INSTALLING THE BUNDLE

  myconote-cli setup --db augustus-fungi --dry-run  # audit list first
  myconote-cli setup --db augustus-fungi            # actually fetch

After install, set the AUGUSTUS_CONFIG_PATH so Augustus picks the
configs up:

  export AUGUSTUS_CONFIG_PATH=$HOME/.myconote/augustus_config

Add this to your `~/.bashrc` or `~/.zshrc` for persistence.
`myconote-cli predict` reads this env var when invoking Augustus."
            ),

            LessonItem::Question {
                prompt: "After `setup --db augustus-fungi`, which environment variable must be set so Augustus finds the new species configs?",
                kind: QuestionKind::FreeText {
                    answer: "AUGUSTUS_CONFIG_PATH",
                    accept_regex: Some(r"(?i)augustus_config_path"),
                },
                hint: Some("It's an Augustus convention; the variable name spells out 'config path' for Augustus."),
                explanation: Some("AUGUSTUS_CONFIG_PATH=$HOME/.myconote/augustus_config. Set it in your shell rc file for persistence across terminal sessions."),
            },

            LessonItem::Text(
"WHEN TO RETRAIN VS USE A BUNDLED SPECIES

Use a bundled species when:
  - Your organism is closely related to one of the 49.
  - You don't have RNA-seq for proper training.
  - You're doing a quick first-pass annotation.

Retrain with `myconote-cli train` when:
  - Your organism is divergent (>20% intron-length difference) from
    every bundled species.
  - You have substantial RNA-seq from your strain.
  - You're producing a publication-quality reference annotation.

A common pattern: use `saccharomyces_cerevisiae_S288C` as a quick
first pass to get an initial gene set, then retrain on your own
organism for the final reference annotation."
            ),

            LessonItem::TryIt {
                instruction: "Write the dry-run command to audit which Augustus species would be downloaded:",
                command_template: "myconote-cli setup --db augustus-fungi --___",
                hint: Some("It's a common --flag for 'show me what you'd do without doing it'."),
            },

            LessonItem::Checkpoint("Augustus species bundle understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 19: Reproducibility Manifests
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_19_reproducibility() -> Lesson {
    Lesson {
        title: "Reproducibility Manifests",
        description: "What's in quant_bundle.json and ase_bundle.json. SHA256 input hashes. Re-running identically from a published bundle.",
        est_minutes: 8,
        items: vec![
            LessonItem::Text(
"REPRODUCIBILITY ISN'T OPTIONAL

If a paper reports 'salmon TPM 3.4 for transcript X in sample Y'
and a reviewer can't reproduce that number from the published
inputs, the paper has a problem.

MycoNote-CLI's RNA-seq + ASE stages each emit a JSON
reproducibility manifest at the end of the run:

  quant_out/quant_bundle.json    after `quant`
  ase_out/ase_bundle.json        after `ase`

Each captures everything a reviewer needs to re-run the analysis
identically: tool versions, parameters, inputs, and SHA256 hashes
of every input file."
            ),

            LessonItem::Text(
"WHAT'S IN quant_bundle.json

  {
    \"tool\": \"myconote-cli\",
    \"version\": \"0.7.1\",
    \"git_sha\": \"abc1234\",
    \"command\": \"myconote-cli quant samples.tsv --transcripts ...\",
    \"timestamp\": \"2026-04-25T19:42:13Z\",
    \"inputs\": {
      \"transcripts.fa\": \"sha256:e3a1b2…\",
      \"decoys.fa\":      \"sha256:7b4c5d…\",
      \"samples.tsv\":    \"sha256:1f2a3b…\"
    },
    \"per_sample\": {
      \"control_rep1\": {
        \"fastq_r1_sha256\": \"…\",
        \"fastq_r2_sha256\": \"…\",
        \"fastp_q30\": 0.96,
        \"salmon_mapping_rate\": 0.92,
        \"strandedness_inferred\": \"unstranded\"
      },
      …
    },
    \"versions\": {
      \"salmon\": \"1.10.2\",
      \"fastp\":  \"0.23.4\"
    }
  }

Every fastq is SHA256-hashed during streaming. So is every input
FASTA / TSV. If anyone re-runs this command with the same inputs,
they MUST get byte-identical bundle JSON (modulo timestamp)."
            ),

            LessonItem::Text(
"ase_bundle.json EXTENDS THE SCHEMA

`ase_bundle.json` carries everything in `quant_bundle.json` plus:

  per_haplotype:
    cds_hap0_sha256:  per-haplotype CDS FASTA hash
    cds_hap1_sha256:
  variants_applied: { snv: 4831, mnp: 73, ins: 412, del: 388 }
  variants_skipped: { exon_boundary_spanning: 12, in_cis_overlap: 4, … }
  per_sample_per_haplotype_mapping_rates:
    control_rep1.hap0: 0.92
    control_rep1.hap1: 0.88
    asymmetry_flag: false   # difference < 5pp
"
            ),

            LessonItem::Question {
                prompt: "What does the SHA256 hash of every input FASTQ in quant_bundle.json prevent?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Slow downloads".into(),
                        "Silent re-running on a different version of the same file (e.g. re-downloaded from SRA)".into(),
                        "Disk corruption during writes".into(),
                        "salmon's index cache misses".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("Filenames don't change but content can."),
                explanation: Some("Filename + size are not enough — a file can be re-downloaded from SRA with a different version, or a sample can be re-prepared. SHA256 catches content drift."),
            },

            LessonItem::Text(
"RE-RUNNING FROM A PUBLISHED BUNDLE

Suppose Paper X publishes its `quant_bundle.json`. As a reviewer
you can:

  1. Read the bundle's `command` field.
  2. Re-download the inputs (the bundle gives you SRA accessions
     in the per_sample block, plus expected SHA256s).
  3. Verify each downloaded SHA256 matches the bundle.
  4. Run the exact command. Compare your bundle to Paper X's.

If the SHA256s match and the command matches, your TPM matrix
should be byte-identical to Paper X's. If it isn't, the
divergence is a real bug — in your salmon version, your fastp
version, or somewhere in MycoNote's orchestration. The bundle
makes the divergence point findable."
            ),

            LessonItem::Question {
                prompt: "True or False: The bundle's `git_sha` field lets a reviewer reproduce the exact MycoNote-CLI build that ran.",
                kind: QuestionKind::TrueFalse { answer: true },
                hint: Some("Pinning the source commit pins the orchestration code."),
                explanation: Some("True. Combined with the version field, the git SHA pins the exact source tree. `cargo install --git https://github.com/K-nie/myconote-cli --rev <sha>` reproduces the build."),
            },

            LessonItem::Checkpoint("Reproducibility understood"),
        ],
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LESSON 20: NCBI Codon Tables (Candida CTG focus)
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_20_genetic_codes() -> Lesson {
    Lesson {
        title: "NCBI Codon Tables (Candida CTG focus)",
        description: "25 NCBI tables, Table 12 for Candida CTG-clade species, the silent-mistranslation failure mode under Standard, propagation through to BRAKER.",
        est_minutes: 12,
        items: vec![
            LessonItem::Text(
"WHEN THE STANDARD CODE IS SILENTLY WRONG

The 'standard' genetic code (NCBI Table 1) is what every
introductory textbook teaches: 64 codons → 20 amino acids, with
TAA/TAG/TGA as stop. It's correct for most eukaryotes.

It is WRONG for several lineages of fungal interest. The most
important for fungal genomics:

  Table 12 — Alternative Yeast Nuclear Code
    Used by ~400 species in the Candida CTG clade.
    The codon CTG, which Standard reads as Leucine, is read as
    SERINE in this lineage.
    Affects: C. albicans, C. tropicalis, C. parapsilosis,
    Debaryomyces hansenii, Meyerozyma guilliermondii, Lodderomyces,
    Spathaspora, …

  Table 26 — Pachysolen Tannophilus Nuclear
    CTG is read as Alanine in P. tannophilus.

  Table 3  — Yeast Mitochondrial
    TGA → Trp; CTN → Thr (different from nuclear).

  Tables 4, 5, 9, 10, 14, 16  — various mitochondrial / ciliate /
    flatworm reassignments.

Refs:
  Santos et al. (2011) Mol Biol Evol 28:2185 (CTG-Ser evolution)
  Mühlhausen et al. (2016) Curr Opin Microbiol 32:16 (CTG-clade overview)"
            ),

            LessonItem::Info("THE SILENT-MISTRANSLATION FAILURE MODE",
"Annotate Candida albicans with the Standard code and your
pipeline will:
  - Run without errors or warnings.
  - Produce protein FASTAs that look fine syntactically.
  - Pass NCBI submission validation.
  - Mistranslate every CTG codon (and there are thousands in a
    typical Candida genome) as Leucine instead of Serine.

Every downstream BLAST hit will be against subtly-wrong
sequences. Every Pfam domain hit will use slightly-wrong
amino acid composition. Every functional prediction will
inherit the error. And nothing will fail loudly until a
reviewer notices the protein doesn't match published Candida
sequences.

This is what 'silent mistranslation' means, and it is the
single most common way to wreck a Candida annotation."
            ),

            LessonItem::Question {
                prompt: "What happens if you annotate a Candida albicans genome with the Standard genetic code (Table 1) instead of Table 12?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Annotation fails with a clear error".into(),
                        "All proteins are silently mistranslated at every CTG codon (CTG read as Leu instead of Ser)".into(),
                        "Augustus refuses to run".into(),
                        "Only mitochondrial genes are affected".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("'Silent' is a key word here."),
                explanation: Some("Silent mistranslation. Every CTG codon (typically 1,000s in a Candida genome) is read as Leu instead of Ser. The pipeline runs to completion with no error, producing subtly-wrong proteins that propagate through every downstream analysis."),
            },

            LessonItem::Text(
"USING TABLE 12 IN MYCONOTE-CLI

The `--genetic-code <n>` flag is consistent across subcommands:

  myconote-cli annotate genes.gff3 --fasta genome.fa --genetic-code 12
  myconote-cli predict masked.fa --use-braker --genetic-code 12
  myconote-cli submit annotated.gff3 --fasta genome.fa --genetic-code 12

The flag forwards through to:
  - The internal CDS translator (annotate / submit)
  - Augustus (via species model selection)
  - GeneMark (via `--gcode 12` on the gmes_petap.pl call)
  - BRAKER (via `--translation_table 12` when --use-braker is set)

This is the 'first-class flag' contribution: you set it ONCE on
the MycoNote command line, and it propagates correctly into every
predictor and every translator. No per-component config files,
no easy-to-forget knobs."
            ),

            LessonItem::Text(
"THE 25-TABLE REGISTRY

MycoNote-CLI's genetic-code registry is data-driven. Tables 1-6,
9-14, 16, 21-31, and 33 are all implemented and tested
per-table against NCBI's reference codon usage. Three tables (27,
28, 31) use context-dependent stop-codon reassignments that
cannot be resolved from codon identity alone — for these, the
tool resolves to the amino-acid reading and flags the limitation
in the annotation report rather than silently falling back to
the Standard code.

Adding a future NCBI-defined table is a single registry entry
(table number, display name, codon deviations from Standard,
optional alternative starts). No match-arm maintenance across
the file."
            ),

            LessonItem::Question {
                prompt: "Which `--genetic-code <n>` value would you pass for an annotation of Pachysolen tannophilus?",
                kind: QuestionKind::FreeText {
                    answer: "26",
                    accept_regex: Some(r"^\s*26\s*$"),
                },
                hint: Some("Pachysolen has its own table — not 12."),
                explanation: Some("Table 26. Pachysolen reads CTG as Alanine, distinct from the Candida CTG clade's Serine reassignment (Table 12). Both deviate from Standard but in different directions."),
            },

            LessonItem::CodeExample(
                "End-to-end Candida CTG-clade annotation, codon-correct:",
                "myconote-cli predict masked.fa --use-braker \\\n  --braker-rna-bam reads.bam --braker-proteins proteins.fa \\\n  --genetic-code 12 --output predict_out\nmyconote-cli annotate predict_out/consensus.gff3 \\\n  --fasta genome.fa --genetic-code 12 \\\n  --interproscan --email you@email.edu\nmyconote-cli submit annotate_out/annotated.gff3 \\\n  --fasta genome.fa --organism 'Candida albicans' \\\n  --genetic-code 12 --locus-prefix CALB"
            ),

            LessonItem::Checkpoint("Genetic-code biology mastered"),

            LessonItem::Text(
"YOU'VE COMPLETED THE FULL CATALOGUE.

20 lessons covering:
  Foundation track  (1-8): pipeline basics through NCBI submission
  RNA-seq track     (9-14): quant, fetch-rna, de-template, ase
  Predictors track  (15-16): BRAKER, GeneMark variants
  Quality track     (17-20): clean, augustus-fungi, reproducibility, codon tables

Where to go from here:
  - Read the manuscript: docs/paper/myconote_manuscript.md
  - Browse the doc site: https://k-nie.github.io/myconote-cli/
  - Try `myconote-cli explain <stage>` to interpret your own outputs

Good luck with your annotations."
            ),
        ],
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Catalogue tests — run with `cargo test --lib lessons::tests`
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalogue_has_twenty_lessons() {
        let all = all_lessons();
        assert_eq!(
            all.len(),
            20,
            "expected 20 lessons in the catalogue (8 foundation + 12 v0.7.1 expansion)"
        );
    }

    #[test]
    fn every_lesson_has_a_title_and_at_least_one_item() {
        for (i, lesson) in all_lessons().iter().enumerate() {
            assert!(
                !lesson.title.is_empty(),
                "lesson {} has an empty title",
                i + 1
            );
            assert!(
                !lesson.description.is_empty(),
                "lesson {} ('{}') has an empty description",
                i + 1,
                lesson.title
            );
            assert!(
                !lesson.items.is_empty(),
                "lesson {} ('{}') has zero items",
                i + 1,
                lesson.title
            );
        }
    }

    #[test]
    fn every_lesson_ends_with_a_verification_step() {
        // Pedagogically: every lesson should end with either a Checkpoint
        // (silent verification) or a Text block summarising what was learned
        // (terminal narrative). A lesson that ends mid-question is malformed.
        for (i, lesson) in all_lessons().iter().enumerate() {
            let last = lesson.items.last().expect("non-empty by previous test");
            let ok = matches!(
                last,
                LessonItem::Checkpoint(_)
                    | LessonItem::Text(_)
                    | LessonItem::Question { .. }
                    | LessonItem::Info(_, _)
            );
            assert!(
                ok,
                "lesson {} ('{}') ends with an unsupported terminal item",
                i + 1,
                lesson.title
            );
        }
    }

    #[test]
    fn lesson_titles_are_unique() {
        let all = all_lessons();
        let mut titles: Vec<&'static str> = all.iter().map(|l| l.title).collect();
        titles.sort();
        let dedup_len = titles
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        assert_eq!(titles.len(), dedup_len, "duplicate lesson titles detected");
    }

    #[test]
    fn rna_seq_track_lessons_present_by_keyword() {
        // Lessons 9-14 are the RNA-seq track. Verify each is wired into
        // the catalogue by searching titles, so a future refactor that
        // accidentally drops one fails this test instead of silently
        // shrinking the tutorial.
        let titles: Vec<&'static str> = all_lessons().iter().map(|l| l.title).collect();
        let must_contain = [
            "quant",
            "fetch-rna",
            "de-template",
            "ase",
            "Binomial",
            "go-template",
        ];
        for needle in must_contain {
            let found = titles.iter().any(|t| t.contains(needle));
            assert!(
                found,
                "expected at least one lesson title to contain '{}', got titles: {:?}",
                needle, titles
            );
        }
    }

    #[test]
    fn modern_predictor_lessons_present() {
        let titles: Vec<&'static str> = all_lessons().iter().map(|l| l.title).collect();
        for needle in ["BRAKER", "GeneMark"] {
            let found = titles.iter().any(|t| t.contains(needle));
            assert!(
                found,
                "expected lesson title containing '{}', got titles: {:?}",
                needle, titles
            );
        }
    }

    #[test]
    fn quality_track_lessons_present() {
        let titles: Vec<&'static str> = all_lessons().iter().map(|l| l.title).collect();
        for needle in [
            "Contig Dedup",
            "Augustus Fungal",
            "Reproducibility",
            "Codon Tables",
        ] {
            let found = titles.iter().any(|t| t.contains(needle));
            assert!(
                found,
                "expected lesson title containing '{}', got titles: {:?}",
                needle, titles
            );
        }
    }

    #[test]
    fn find_lesson_by_name_works_for_new_lessons() {
        // Sanity: name-based lookup should resolve each new lesson.
        for needle in [
            "quant",
            "fetch-rna",
            "de-template",
            "ase",
            "BRAKER",
            "GeneMark",
            "Contig Dedup",
            "Reproducibility",
            "Codon",
        ] {
            assert!(
                find_lesson_by_name(needle).is_some(),
                "find_lesson_by_name('{}') returned None",
                needle
            );
        }
    }

    #[test]
    fn lesson_durations_are_reasonable() {
        // Each lesson should be 5..=20 minutes — anything outside that
        // suggests a stale or runaway lesson.
        for (i, lesson) in all_lessons().iter().enumerate() {
            assert!(
                lesson.est_minutes >= 5 && lesson.est_minutes <= 20,
                "lesson {} ('{}') has est_minutes {} outside [5,20]",
                i + 1,
                lesson.title,
                lesson.est_minutes
            );
        }
    }

    #[test]
    fn no_lesson_references_removed_subcommands() {
        // Phylogeny / place / Y1000+ were removed in v0.5.0. No lesson
        // body should mention them as live subcommands.
        let dead = ["myconote-cli phylogeny", "myconote-cli place", "y1000plus"];
        for lesson in all_lessons() {
            for item in &lesson.items {
                let body: String = match item {
                    LessonItem::Text(t) => (*t).to_string(),
                    LessonItem::Info(title, body) => format!("{}\n{}", title, body),
                    LessonItem::CodeExample(d, c) => format!("{}\n{}", d, c),
                    LessonItem::Demo(d, c) => format!("{}\n{}", d, c),
                    LessonItem::TryIt {
                        instruction,
                        command_template,
                        ..
                    } => format!("{}\n{}", instruction, command_template),
                    LessonItem::Checkpoint(m) => (*m).to_string(),
                    LessonItem::Question {
                        prompt,
                        explanation,
                        ..
                    } => format!("{}\n{}", prompt, explanation.unwrap_or("")),
                };
                let lower = body.to_lowercase();
                for needle in dead {
                    assert!(
                        !lower.contains(&needle.to_lowercase()),
                        "lesson '{}' still references removed subcommand/term '{}'",
                        lesson.title,
                        needle
                    );
                }
            }
        }
    }
}
