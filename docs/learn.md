# Interactive Tutorial (`myconote-cli learn`)

Learn myconote-cli step by step, right in your terminal. Inspired by R's **swirl** package — no browser needed, no external tools required.

---

## Getting Started

```bash
myconote-cli learn          # list all lessons and your progress
myconote-cli learn 1        # start lesson 1
myconote-cli learn basics   # start a lesson by name
myconote-cli learn --resume # pick up where you left off
myconote-cli learn --reset  # clear all progress
```

Aliases also work: `myconote-cli tutorial` or `myconote-cli swirl`

---

## Available Lessons

The catalogue is split into four tracks. Lessons within each track can be completed in any order, but new users should walk the foundation track first.

### Foundation track (run in order if you're new)

| # | Lesson | Duration | Topics |
|---|--------|----------|--------|
| 1 | **Welcome to MycoNote** | ~12 min | Pipeline overview, kingdoms, file formats, key concepts |
| 2 | **Setup & Installation** | ~10 min | `install`, `setup`, `check`, `species` commands |
| 3 | **Pre-processing: Sort & Mask** | ~10 min | `sort`, masking engines, soft vs hard masking |
| 4 | **Gene Prediction** | ~14 min | Augustus, SNAP, EVM consensus, protein evidence, ploidy, weights |
| 5 | **Functional Annotation** | ~15 min | Swiss-Prot, Pfam, BUSCO, tRNA, genetic codes, InterProScan |
| 6 | **Advanced: Training & Evidence** | ~12 min | RNA-seq training, Trinity, PASA, custom evidence weights |
| 7 | **NCBI Submission** | ~10 min | GFF3 validation, `.tbl` generation, table2asn, BioProject |
| 8 | **Stats, Comparative Genomics, and Conversion** | ~12 min | `stats`, `compare`, `convert`, handoff to Proksee/IGV/clinker |

### RNA-seq track (lessons 9-14, completable in any order)

| # | Lesson | Duration | Topics |
|---|--------|----------|--------|
| 9 | **RNA-seq Quantification (`quant`)** | ~12 min | salmon decoy-aware index, fastp QC, tximport-ready outputs, SHA256 cache |
| 10 | **SRA/ENA Download (`fetch-rna`)** | ~10 min | ENA REST default, sra-toolkit fallback, MD5 verification, sample-sheet emission |
| 11 | **Differential Expression (`de-template`)** | ~12 min | tximport + DESeq2 + apeglm template, design formulas, contrasts |
| 12 | **Allele-Specific Expression (`ase`)** | ~14 min | Phased VCF + personalized transcriptome, unphased-het hard error, asymmetry detection |
| 13 | **ASE Binomial Test (`ase-template`)** | ~10 min | Per-sample-corrected null, informativeness filter, base-R only |
| 14 | **GO Enrichment (`go-template`)** | ~10 min | topGO Fisher's exact test, BP/MF/CC, BH within ontology, dot plot |

### Modern predictors track (lessons 15-16)

| # | Lesson | Duration | Topics |
|---|--------|----------|--------|
| 15 | **BRAKER for Maximum Accuracy** | ~10 min | `predict --use-braker`, BRAKER1/2/3 mode auto-detection, conflict guard |
| 16 | **GeneMark Variants (`--genemark-mode`)** | ~10 min | ES/ET/EP+/ETP+, ProtHint, when each is appropriate |

### Quality and reproducibility track (lessons 17-20)

| # | Lesson | Duration | Topics |
|---|--------|----------|--------|
| 17 | **Contig Dedup (`clean --mode contigs`)** | ~8 min | minimap2 self-alignment, coverage/identity thresholds, deterministic ties |
| 18 | **Augustus Fungal Species Bundle** | ~8 min | The 49-species manifest, `AUGUSTUS_CONFIG_PATH`, retrain vs bundled |
| 19 | **Reproducibility Manifests** | ~8 min | `quant_bundle.json` / `ase_bundle.json`, SHA256 hashing, byte-identical re-runs |
| 20 | **NCBI Codon Tables (Candida CTG focus)** | ~12 min | 25 NCBI tables, Table 12, silent-mistranslation failure mode |

**Total: ~217 minutes (~3.6 hours)** covering the full pipeline plus modern features added in v0.5+.

---

## During a Lesson

Type your answer and press Enter. Special commands:

| Command | What it does |
|---------|-------------|
| `hint` | Get a contextual hint for the current question |
| `info` | Peek at the expected answer |
| `skip` | Move to the next question without answering |
| `quit` / `bye` | Save progress and exit |

After 3 wrong attempts, the correct answer is revealed automatically — the goal is learning, not testing.

---

## Question Types

The tutorial uses 5 question types:

**Multiple choice** — pick a, b, c, or d
```
  Which tool does myconote-cli use for protein homology search?
    a) BLAST
    b) MMseqs2
    c) DIAMOND
    d) HMMER

  > b
  Correct!
```

**Free text** — type the answer (fuzzy-matched, typos forgiven)
```
  What command installs all missing external tools?

  > install
  Correct!
```

**True/False**
```
  True or False: myconote-cli is focused primarily on fungal genomes.

  > true
  Correct! Defaults, benchmarks, and the test fixtures are all tuned
  for fungal annotation. The `--kingdom` flag accepts plant/animal/
  insect/protist too, but those paths are experimental and not
  validated on large genomes.
```

**Fill-in-the-blank** — complete a command template
```
  Complete the command:
  myconote-cli ___ genes.gff3 --to gtf

  > convert
  Correct!
```

**Order steps** — arrange pipeline stages correctly
```
  Put these steps in the correct pipeline order:
    1) annotate
    2) predict
    3) mask
    4) sort

  > 4,3,2,1
  Correct! sort -> mask -> predict -> annotate
```

---

## Progress Tracking

Your progress is saved to `~/.myconote/learn_progress.json` and persists between sessions.

- Completed lessons show a green checkmark in the lesson list
- `--resume` picks up exactly where you left off
- `--reset` clears all progress to start fresh

---

## For Workshop Instructors

The `learn` system can be used in a classroom setting:

1. Students install myconote-cli on their laptops
2. Each student works through lessons at their own pace
3. No internet required (all lessons are built into the binary)
4. Progress is per-user, so students can resume between sessions

For a longer workshop format, see the [Workshop Lesson](lesson.md) page which provides a 3-hour instructor-guided curriculum.
