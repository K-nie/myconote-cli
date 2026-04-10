/// Lesson content for the interactive tutorial
///
/// Each lesson is a sequence of items: text blocks, questions,
/// code examples, live demos, and checkpoints. Lessons are ordered
/// from basic to advanced, following the pipeline order.

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
        description: "Get oriented: what myconote-cli does, how the pipeline works, and key concepts.",
        est_minutes: 8,
        items: vec![
            LessonItem::Text(
"Welcome! myconote-cli is a genome annotation pipeline written in Rust.
It takes a genome assembly (FASTA) and produces a fully annotated gene set
(GFF3) with functional descriptions, domain annotations, GO terms, and more.

Think of it as an assembly line for genome annotation:
  assembly.fa → sort → mask → train → predict → update → annotate → submit

Each step builds on the previous one. Let's learn what each does."
            ),

            LessonItem::Info("Why Rust?",
"myconote-cli is written in Rust for speed and safety.
It uses parallel processing (rayon), memory-mapped I/O (memmap2),
and zero-copy parsing for blazing performance on large genomes.
External tools (Augustus, SNAP, etc.) are called as subprocesses."
            ),

            LessonItem::Question {
                prompt: "What file format does a genome assembly typically come in?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "GFF3".into(),
                        "FASTA".into(),
                        "GenBank".into(),
                        "VCF".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It's a simple format with > headers followed by sequences."),
                explanation: Some("FASTA (.fa, .fasta, .fna) is the standard format for genome assemblies."),
            },

            LessonItem::Text(
"The myconote-cli pipeline has 7 main stages:

  1. sort     — Sort contigs by length, clean up headers
  2. mask     — Identify and soft-mask repetitive regions
  3. train    — Train gene predictors using RNA-seq data
  4. predict  — Predict gene models (Augustus, SNAP, EVM consensus)
  5. update   — Refine predictions with transcript evidence (UTRs)
  6. annotate — Add functional annotations (Swiss-Prot, Pfam, GO, etc.)
  7. submit   — Prepare files for NCBI GenBank submission"
            ),

            LessonItem::Question {
                prompt: "Which step comes FIRST in the pipeline?",
                kind: QuestionKind::FreeText {
                    answer: "sort",
                    accept_regex: Some(r"(?i)^sort"),
                },
                hint: Some("We need clean, sorted contig names before anything else."),
                explanation: Some("'sort' renames and orders contigs, giving you a clean starting point."),
            },

            LessonItem::Question {
                prompt: "What does 'soft-masking' a genome mean?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Deleting repetitive sequences from the genome".into(),
                        "Replacing repeats with N characters".into(),
                        "Converting repetitive bases to lowercase letters".into(),
                        "Compressing the FASTA file".into(),
                    ],
                    correct_index: 2,
                },
                hint: Some("Think lowercase vs uppercase in a FASTA file."),
                explanation: Some("Soft-masking uses lowercase letters (atcg) for repeats while keeping the sequence intact. Gene predictors can then ignore these regions."),
            },

            LessonItem::Checkpoint("Basic concepts understood"),

            LessonItem::Text(
"myconote-cli supports 5 kingdoms, each with tailored defaults:

  fungi    — Default. Compact genomes, short introns (~40-2000 bp).
  plant    — Large genomes, many introns (up to 50 kb).
  animal   — Very large genomes, huge introns (up to 500 kb).
  insect   — Medium genomes with moderate introns.
  protist  — Diverse group with variable genome sizes."
            ),

            LessonItem::Question {
                prompt: "For a Saccharomyces cerevisiae genome, which kingdom would you use?",
                kind: QuestionKind::FreeText {
                    answer: "fungi",
                    accept_regex: Some(r"(?i)^fung"),
                },
                hint: Some("Yeast is a fungus!"),
                explanation: Some("S. cerevisiae is a fungus. Use --kingdom fungi (the default)."),
            },

            LessonItem::CodeExample(
                "To see all available commands:",
                "myconote-cli --help"
            ),

            LessonItem::Question {
                prompt: "True or False: myconote-cli can only annotate fungal genomes.",
                kind: QuestionKind::TrueFalse { answer: false },
                hint: Some("Remember the 5 kingdoms we just discussed."),
                explanation: Some("False — it supports fungi, plants, animals, insects, and protists."),
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
        description: "Install dependencies, download databases, and verify your environment.",
        est_minutes: 6,
        items: vec![
            LessonItem::Text(
"Before annotating a genome, you need:
  1. External tools (Augustus, SNAP, HMMER, MMseqs2, etc.)
  2. Reference databases (Swiss-Prot, Pfam, BUSCO lineages)

myconote-cli can install everything automatically!"
            ),

            LessonItem::CodeExample(
                "Check which tools are already installed:",
                "myconote-cli check"
            ),

            LessonItem::Question {
                prompt: "What command installs all missing external tools?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli ___",
                    answer: "install",
                },
                hint: Some("It's a single word that means 'put these tools on my system'."),
                explanation: Some("'myconote-cli install' scans for missing tools and installs them via conda/mamba."),
            },

            LessonItem::Info("Package managers",
"myconote-cli prefers mamba (faster) but falls back to conda.
If neither is available, it tells you exactly what to install manually.
Some tools (GeneMark) require a licence and can't be auto-installed."
            ),

            LessonItem::Question {
                prompt: "What command downloads and indexes the reference databases?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli ___",
                    answer: "setup",
                },
                hint: Some("You're setting up the databases your annotations will search against."),
                explanation: Some("'myconote-cli setup' downloads Swiss-Prot, Pfam, dbCAN, MEROPS, and BUSCO lineages."),
            },

            LessonItem::CodeExample(
                "Download only specific databases:",
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
                explanation: Some("Databases go to ~/.myconote/dbs/ by default. Override with --db-dir."),
            },

            LessonItem::Checkpoint("Environment ready"),

            LessonItem::Question {
                prompt: "What command shows available pre-trained Augustus species models?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli ___",
                    answer: "species",
                },
                hint: Some("You want to see what species are available for gene prediction."),
                explanation: Some("'myconote-cli species' lists all Augustus species on your system, grouped by kingdom."),
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
        description: "Clean up your assembly and mask repetitive regions before gene prediction.",
        est_minutes: 7,
        items: vec![
            LessonItem::Text(
"Before predicting genes, we prepare the genome in two steps:

  1. SORT — Rename contigs to clean IDs (scaffold_1, scaffold_2, ...),
             sort by length, and optionally filter short scaffolds.

  2. MASK — Identify and soft-mask repetitive elements so gene
             predictors don't get confused by transposons and satellites."
            ),

            LessonItem::CodeExample(
                "Sort contigs and remove anything under 500 bp:",
                "myconote-cli sort assembly.fa --min-length 500"
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
                explanation: Some("Consistent, clean IDs prevent naming conflicts throughout the pipeline."),
            },

            LessonItem::Text(
"Masking has several engines — choose based on your needs:

  self          Fast minimap2 self-alignment (no database needed)
  repeatmasker  RepeatMasker with species repeat library
  repeatmodeler De novo repeat library then RepeatMasker [RECOMMENDED]
  full          All methods combined (most thorough)"
            ),

            LessonItem::TryIt {
                instruction: "Write the command to mask a genome using the repeatmodeler engine with 8 threads:",
                command_template: "myconote-cli mask genome.fa --engine ___ --threads ___",
                hint: Some("The engine is 'repeatmodeler' and threads is 8."),
            },

            LessonItem::Question {
                prompt: "What is the RECOMMENDED masking engine for a new fungal genome?",
                kind: QuestionKind::FreeText {
                    answer: "repeatmodeler",
                    accept_regex: Some(r"(?i)repeatmodel"),
                },
                hint: Some("It builds a de novo repeat library specific to your genome."),
                explanation: Some("repeatmodeler creates a custom repeat library, then feeds it to RepeatMasker. Best for novel genomes."),
            },

            LessonItem::Checkpoint("Pre-processing mastered"),

            LessonItem::Question {
                prompt: "True or False: Hard-masking (replacing repeats with N) is preferred for gene prediction.",
                kind: QuestionKind::TrueFalse { answer: false },
                hint: Some("If you replace sequence with Ns, you lose information."),
                explanation: Some("False — soft-masking (lowercase) is preferred. Predictors can still read the sequence but know to be cautious."),
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
        description: "Predict genes using Augustus, SNAP, and the Evidence Modeler consensus.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"Gene prediction is the heart of the pipeline. myconote-cli uses
multiple ab initio predictors and merges them with evidence:

  Augustus    — Most accurate ab initio predictor (weight: 10)
  SNAP        — Fast secondary predictor (weight: 3)
  GlimmerHMM  — Optional third predictor (weight: 2)
  GeneMark-ES — Self-training predictor (weight: 5)
  Protein     — miniprot/exonerate protein alignment (weight: 20)

The Evidence Modeler merges overlapping predictions, keeping the
highest-scoring model at each locus."
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
                hint: Some("Homology to known proteins is the strongest evidence."),
                explanation: Some("Protein evidence has weight 20 — homology to real proteins is stronger than any single ab initio prediction."),
            },

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
                explanation: Some("--locus-prefix MYORG gives genes IDs like MYORG_000001, MYORG_000002, etc."),
            },

            LessonItem::Info("Evidence Weights",
"You can customize weights with a TOML file:

  # weights.toml
  augustus = 10.0
  snap = 3.0
  protein = 25.0      # increase protein weight
  genemark = 5.0

  myconote-cli predict genome.fa --weights weights.toml"
            ),

            LessonItem::TryIt {
                instruction: "Write the command to predict genes with protein evidence from Swiss-Prot:",
                command_template: "myconote-cli predict masked.fa --kingdom fungi --protein-fasta ___",
                hint: Some("Point --protein-fasta to your Swiss-Prot FASTA file."),
            },

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
                hint: Some("Think about what happens when you assemble a diploid genome."),
                explanation: Some("--ploidy 2 tells myconote-cli to expect allelic duplicates in a diploid assembly and adjust accordingly."),
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
                    correct_order: vec![3, 2, 1, 0],  // sort, mask, predict, annotate
                },
                hint: Some("Start with sort, end with annotate."),
                explanation: Some("The order is: sort → mask → predict → annotate."),
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
        description: "Add biological meaning to predicted genes: Swiss-Prot, Pfam, GO terms, BUSCO, tRNA, and more.",
        est_minutes: 10,
        items: vec![
            LessonItem::Text(
"After predicting genes, we need to figure out what they DO.
myconote-cli runs multiple annotation sources in parallel:

  Core (always on):
    MMseqs2     → Swiss-Prot homology (product names)
    Pfam        → Domain architecture
    BUSCO       → Completeness assessment
    GO terms    → Gene Ontology from UniProt hits

  Optional (enable with flags):
    --trnascan     → tRNA gene prediction
    --eggnog       → COG/NOG functional categories + KEGG
    --cazyme       → Carbohydrate-active enzymes (dbCAN)
    --secretome    → Signal peptides + transmembrane domains
    --antismash    → Biosynthetic gene clusters
    --merops       → Protease families
    --interproscan → Comprehensive domain search (EBI API)"
            ),

            LessonItem::Question {
                prompt: "Which tool does myconote-cli use for protein homology search against Swiss-Prot?",
                kind: QuestionKind::FreeText {
                    answer: "MMseqs2",
                    accept_regex: Some(r"(?i)mmseqs"),
                },
                hint: Some("It's much faster than BLAST and starts with 'MM'."),
                explanation: Some("MMseqs2 is used for fast, sensitive protein homology search."),
            },

            LessonItem::CodeExample(
                "Run annotation with all the bells and whistles:",
                "myconote-cli annotate genes.gff3 --fasta genome.fa \\\n  --eggnog --cazyme --secretome --antismash --merops --trnascan"
            ),

            LessonItem::Question {
                prompt: "What does BUSCO measure?",
                kind: QuestionKind::MultipleChoice {
                    choices: vec![
                        "Gene expression levels".into(),
                        "Genome completeness using conserved single-copy orthologs".into(),
                        "The number of tRNA genes".into(),
                        "Protein folding accuracy".into(),
                    ],
                    correct_index: 1,
                },
                hint: Some("It tells you how 'complete' your annotation is."),
                explanation: Some("BUSCO checks whether expected universal single-copy genes are present, giving a completeness score."),
            },

            LessonItem::Info("Genetic Code",
"For organisms with non-standard genetic codes, use --genetic-code:

  1  = Standard (default)
  12 = Candida CTG clade (CTG = Ser instead of Leu)
  3  = Yeast mitochondrial
  4  = Mold mitochondrial

  myconote-cli annotate genes.gff3 --fasta genome.fa --genetic-code 12"
            ),

            LessonItem::Question {
                prompt: "If you're annotating a Candida albicans genome, which genetic code should you use?",
                kind: QuestionKind::FreeText {
                    answer: "12",
                    accept_regex: Some(r"(?i)(12|candida|ctg)"),
                },
                hint: Some("Candida uses the CTG clade code where CTG encodes Serine."),
                explanation: Some("Table 12 (Alternative Yeast Nuclear Code) — CTG = Ser instead of Leu."),
            },

            LessonItem::Checkpoint("Annotation pipeline mastered"),

            LessonItem::Question {
                prompt: "Which flag enables tRNA gene prediction?",
                kind: QuestionKind::FreeText {
                    answer: "--trnascan",
                    accept_regex: Some(r"(?i)--(trna|trnascan)"),
                },
                hint: Some("It wraps the tRNAscan-SE tool."),
                explanation: Some("--trnascan runs tRNAscan-SE to find tRNA genes in your genome."),
            },

            LessonItem::Question {
                prompt: "True or False: If a tool is missing, myconote-cli will crash.",
                kind: QuestionKind::TrueFalse { answer: false },
                hint: Some("Think about 'graceful degradation'."),
                explanation: Some("False — myconote-cli degrades gracefully. If a tool is missing, it skips that step and continues."),
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
        description: "RNA-seq training, custom evidence weights, and protein-to-genome alignment.",
        est_minutes: 8,
        items: vec![
            LessonItem::Text(
"For the best gene models, TRAIN your own Augustus/SNAP models
using RNA-seq data from your organism. This is optional but strongly
recommended if you have RNA-seq.

The training pipeline:
  RNA-seq reads → Trinity assembly → minimap2 alignment →
  PASA transcript database → Extract complete gene models →
  Train Augustus + SNAP"
            ),

            LessonItem::CodeExample(
                "Train predictors from paired-end RNA-seq:",
                "myconote-cli train genome_masked.fa \\\n  --left R1.fastq.gz --right R2.fastq.gz \\\n  --species my_organism"
            ),

            LessonItem::Question {
                prompt: "What tool does myconote-cli use to assemble RNA-seq reads into transcripts?",
                kind: QuestionKind::FreeText {
                    answer: "Trinity",
                    accept_regex: Some(r"(?i)trinity"),
                },
                hint: Some("It's named after a concept of three-in-one."),
                explanation: Some("Trinity performs de novo RNA-seq assembly, producing transcript sequences for training."),
            },

            LessonItem::Text(
"After training, use your custom species model in prediction:

  myconote-cli predict genome_masked.fa \\
    --species my_organism \\
    --snap-hmm train_out/snap_training.hmm

This gives much better results than using a generic pre-trained model!"
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
                hint: Some("It uses PASA to add 5'/3' UTRs and fix exon boundaries."),
                explanation: Some("'update' refines gene models — adding UTRs, correcting intron/exon boundaries, and finding alternative isoforms via PASA."),
            },

            LessonItem::Checkpoint("Advanced training understood"),

            LessonItem::Question {
                prompt: "What file format is used for custom evidence weights?",
                kind: QuestionKind::FreeText {
                    answer: "TOML",
                    accept_regex: Some(r"(?i)toml"),
                },
                hint: Some("It's a simple configuration format used in Rust projects (Cargo uses it)."),
                explanation: Some("TOML files (e.g. weights.toml) let you customize the weight given to each predictor."),
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
        description: "Prepare your annotation for GenBank submission using the submit command.",
        est_minutes: 6,
        items: vec![
            LessonItem::Text(
"Once your genome is annotated, you'll want to submit it to NCBI GenBank.
myconote-cli automates the entire submission prep:

  1. Validates your GFF3 for NCBI compliance
  2. Generates a feature table (.tbl format)
  3. Runs table2asn to create the submission file (.sqn)
  4. Produces all required metadata files"
            ),

            LessonItem::CodeExample(
                "Prepare a submission:",
                "myconote-cli submit annotated.gff3 --fasta genome.fa \\\n  --organism 'Aspergillus niger' \\\n  --strain CBS 513.88 \\\n  --locus-prefix ASPNI \\\n  --email your@email.edu"
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
                hint: Some("It validates the GFF3 for common structural problems."),
                explanation: Some("The validator checks for duplicate IDs, missing Parent references, coordinate issues, and required attributes."),
            },

            LessonItem::Question {
                prompt: "What flag would you add to only validate without generating files?",
                kind: QuestionKind::FreeText {
                    answer: "--validate-only",
                    accept_regex: Some(r"(?i)--validate"),
                },
                hint: Some("You just want to validate, nothing else."),
                explanation: Some("--validate-only runs NCBI compliance checks without generating submission files."),
            },

            LessonItem::Info("BioProject & BioSample",
"Before submitting to NCBI, register your project:
  1. Create a BioProject at https://submit.ncbi.nlm.nih.gov/
  2. Create a BioSample for your organism/strain
  3. Pass the accessions to myconote-cli:
     --bioproject PRJNA123456 --biosample SAMN12345678"
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
        description: "Statistics, plots, phylogenetics, synteny, and format conversion.",
        est_minutes: 8,
        items: vec![
            LessonItem::Text(
"myconote-cli includes powerful analysis tools beyond annotation:

  stats     — Gene counts, lengths, GC%, exon stats, N50
  plot      — Linear and circular genome maps
  phylogeny — ML trees with IQ-TREE (ModelFinder + UFBoot)
  compare   — Multi-genome comparison
  synteny   — Ribbon diagrams between two genomes
  view      — Interactive browser (JBrowse2 / UCSC)
  convert   — Format conversion (GFF3 ↔ GTF/BED/GenBank/FASTA)"
            ),

            LessonItem::CodeExample(
                "Generate statistics with taxonomic benchmarking:",
                "myconote-cli stats annotated.gff3 --taxon fungi"
            ),

            LessonItem::Question {
                prompt: "What command creates a circular genome map?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli ___ genes.gff3 --type circular --output map.png",
                    answer: "plot",
                },
                hint: Some("You want to plot the genome."),
                explanation: Some("'myconote-cli plot' generates publication-quality genome maps."),
            },

            LessonItem::Question {
                prompt: "How would you convert a GFF3 to GenBank format?",
                kind: QuestionKind::FillBlank {
                    template: "myconote-cli convert genes.gff3 --to ___ --fasta genome.fa",
                    answer: "genbank",
                },
                hint: Some("The target format is GenBank (.gbk)."),
                explanation: Some("--to genbank converts GFF3 to GenBank format (requires --fasta for the sequences)."),
            },

            LessonItem::TryIt {
                instruction: "Write a command to build a phylogenetic tree from a multiple sequence alignment:",
                command_template: "myconote-cli phylogeny alignment.fa --model MFP --bootstrap 1000",
                hint: Some("MFP = ModelFinder Plus (auto-selects the best model)."),
            },

            LessonItem::Question {
                prompt: "What tool does myconote-cli use for phylogenetic inference?",
                kind: QuestionKind::FreeText {
                    answer: "IQ-TREE",
                    accept_regex: Some(r"(?i)iq.?tree"),
                },
                hint: Some("It's a modern ML tree builder. Starts with IQ."),
                explanation: Some("IQ-TREE 2 is used for maximum-likelihood phylogenetic inference with ultrafast bootstrapping."),
            },

            LessonItem::Checkpoint("Analysis tools mastered"),

            LessonItem::Text(
"Congratulations! You've completed all lessons.

Here's the full workflow for annotating a new genome:

  # 1. Install tools and databases
  myconote-cli install
  myconote-cli setup

  # 2. Pre-process
  myconote-cli sort assembly.fa --min-length 500
  myconote-cli mask assembly_sorted.fa --engine repeatmodeler

  # 3. Train (if you have RNA-seq)
  myconote-cli train assembly_masked.fa --left R1.fq --right R2.fq --species myorg

  # 4. Predict genes
  myconote-cli predict assembly_masked.fa --kingdom fungi --species myorg

  # 5. Annotate
  myconote-cli annotate predict_out/consensus.gff3 --fasta assembly.fa \\
    --eggnog --cazyme --trnascan --secretome

  # 6. Submit to NCBI
  myconote-cli submit annotate_out/annotated.gff3 --fasta assembly.fa \\
    --organism 'My organism'

You're ready to annotate real genomes. Good luck!"
            ),
        ],
    }
}
