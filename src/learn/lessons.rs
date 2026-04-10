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
        lesson_1_welcome(),
        lesson_2_setup(),
        lesson_3_sort_mask(),
        lesson_4_predict(),
        lesson_5_annotate(),
        lesson_6_advanced(),
        lesson_7_submit(),
        lesson_8_analysis(),
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
                prompt: "True or False: myconote-cli can only annotate fungal genomes.",
                kind: QuestionKind::TrueFalse { answer: false },
                hint: Some("Remember the 5 kingdoms we just discussed."),
                explanation: Some("False. It supports fungi, plants, animals, insects, and protists — each with kingdom-specific parameter tuning."),
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

myconote-cli supports 18 NCBI translation tables — more than any
other annotation pipeline."
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
// LESSON 8: Analysis & Visualization
// ═══════════════════════════════════════════════════════════════════════════════

fn lesson_8_analysis() -> Lesson {
    Lesson {
        title: "Analysis & Visualization",
        description: "Statistics, genome maps, phylogenetics, comparative genomics, and format conversion.",
        est_minutes: 12,
        items: vec![
            LessonItem::Text(
"BEYOND ANNOTATION: ANALYSIS TOOLS

myconote-cli includes analysis capabilities that most annotation
pipelines lack. These tools let you explore, compare, and publish
your results without switching to separate software."
            ),

            LessonItem::Text(
"GENOME STATISTICS (myconote-cli stats)

After annotation, the first question is: 'Does this look right?'
The stats command calculates:

  - Gene count, transcript count, CDS/exon counts
  - Gene length statistics (mean, median, N50)
  - GC content per contig
  - Exons per gene distribution
  - Per-chromosome breakdown

With --taxon fungi, it benchmarks your numbers against expected
ranges. For example, a typical fungal genome has:
  - 5,000-12,000 genes
  - Mean gene length of 1,200-1,800 bp
  - 1-3 exons per gene
  - 45-55% GC content

If your numbers are wildly different, something may be wrong with
your assembly or prediction."
            ),

            LessonItem::CodeExample(
                "Generate statistics with taxonomic benchmarking:",
                "myconote-cli stats annotated.gff3 --taxon fungi"
            ),

            LessonItem::Text(
"PHYLOGENETICS (myconote-cli phylogeny)

Build maximum-likelihood phylogenetic trees using IQ-TREE 2, the
state-of-the-art ML tree inference program.

Key features:
  - ModelFinder Plus (MFP): automatically selects the best-fit
    substitution model from 286 candidates
  - Ultrafast Bootstrap (UFBoot2): 1000 replicates in minutes
    instead of hours
  - Partitioned analysis for multi-gene datasets

Ref: Minh et al. (2020) Mol Biol Evol 37:1530 (IQ-TREE 2)
Ref: Kalyaanamoorthy et al. (2017) Nature Methods 14:587 (ModelFinder)"
            ),

            LessonItem::Question {
                prompt: "What tool does myconote-cli use for phylogenetic inference?",
                kind: QuestionKind::FreeText {
                    answer: "IQ-TREE",
                    accept_regex: Some(r"(?i)iq.?tree"),
                },
                hint: Some("It's the most-cited ML tree builder. Starts with 'IQ'."),
                explanation: Some("IQ-TREE 2 (Minh et al. 2020) performs maximum-likelihood phylogenetic inference with automatic model selection and ultrafast bootstrapping."),
            },

            LessonItem::Text(
"VISUALIZATION

  plot — Genome maps
    Linear: horizontal tracks showing gene density per contig
    Circular: ideogram-style whole-genome view (like Circos but built-in)
    Output: PNG or SVG for publication

  view — Interactive genome browser
    JBrowse2: generates a self-contained HTML file with an interactive
              browser — no server needed, works offline
    UCSC: generates custom track URLs for the UCSC Genome Browser

  synteny — Comparative ribbon diagrams
    Compares two genomes using minimap2 whole-genome alignment and
    draws ribbons connecting syntenic (conserved order) regions"
            ),

            LessonItem::Question {
                prompt: "What command creates a circular genome map?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli ___ genes.gff3 --type circular --output map.png",
                    answer: "plot",
                },
                hint: Some("You want to plot the genome."),
                explanation: Some("'myconote-cli plot' generates publication-quality genome maps. Use --type linear for horizontal tracks or --type circular for Circos-style rings."),
            },

            LessonItem::Text(
"FORMAT CONVERSION (myconote-cli convert)

Bioinformatics uses dozens of file formats. myconote-cli converts
between 15+ formats:

  GFF3 conversions:
    --to gtf       Gene Transfer Format (Ensembl/GENCODE tools)
    --to bed       Browser Extensible Data (UCSC, bedtools)
    --to genbank   GenBank flat file (Geneious, SnapGene)
    --to protein   Translated protein FASTA (BLAST, OrthoFinder)

  Sequence conversions:
    --to phylip    PHYLIP format (RAxML, PhyML)
    --to nexus     NEXUS format (MrBayes, PAUP*)
    --to fasta     FASTQ -> FASTA (strip quality scores)

  VCF conversions:
    --to table     Tabular variant summary
    --to consensus Consensus FASTA with variants applied"
            ),

            LessonItem::Question {
                prompt: "How would you convert a GFF3 to GenBank format?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli convert genes.gff3 --to ___ --fasta genome.fa",
                    answer: "genbank",
                },
                hint: Some("The target format is GenBank (.gbk)."),
                explanation: Some("--to genbank produces a GenBank flat file compatible with Geneious, SnapGene, and NCBI tools. Requires --fasta because GenBank includes the sequence."),
            },

            LessonItem::Checkpoint("Analysis tools mastered"),

            LessonItem::Text(
"CONGRATULATIONS! You've completed all 8 lessons.

Here's the full workflow for annotating a new genome:

  # 1. Install tools and databases (one-time)
  myconote-cli install --yes
  myconote-cli setup

  # 2. Pre-process
  myconote-cli sort assembly.fa --min-length 500
  myconote-cli mask assembly_sorted.fa --engine repeatmodeler

  # 3. Train (optional — if you have RNA-seq)
  myconote-cli train assembly_masked.fa --left R1.fq --right R2.fq

  # 4. Predict genes
  myconote-cli predict assembly_masked.fa --kingdom fungi

  # 5. Annotate
  myconote-cli annotate predict_out/consensus.gff3 --fasta assembly.fa \\
    --trnascan --interproscan --email you@email.edu

  # 6. Submit to NCBI
  myconote-cli submit annotate_out/annotated.gff3 --fasta assembly.fa \\
    --organism 'My organism' --validate-only

KEY REFERENCES:
  Stanke et al. (2006) BMC Bioinformatics 7:62 (Augustus)
  Steinegger & Soding (2017) Nature Biotech 35:1026 (MMseqs2)
  Manni et al. (2021) Mol Biol Evol 38:4647 (BUSCO)
  Minh et al. (2020) Mol Biol Evol 37:1530 (IQ-TREE 2)
  Haas et al. (2008) Genome Biol 9:R7 (EvidenceModeler)

You're ready to annotate real genomes. Good luck!"
            ),
        ],
    }
}
