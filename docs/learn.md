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

| # | Lesson | Duration | Topics |
|---|--------|----------|--------|
| 1 | **Welcome to MycoNote** | ~8 min | Pipeline overview, kingdoms, file formats, key concepts |
| 2 | **Setup & Installation** | ~6 min | `install`, `setup`, `check`, `species` commands |
| 3 | **Pre-processing: Sort & Mask** | ~7 min | `sort`, masking engines, soft vs hard masking |
| 4 | **Gene Prediction** | ~10 min | Augustus, SNAP, EVM consensus, protein evidence, ploidy, weights |
| 5 | **Functional Annotation** | ~10 min | Swiss-Prot, Pfam, BUSCO, tRNA, genetic codes, InterProScan |
| 6 | **Advanced: Training & Evidence** | ~8 min | RNA-seq training, Trinity, PASA, custom evidence weights |
| 7 | **NCBI Submission** | ~6 min | GFF3 validation, `.tbl` generation, table2asn, BioProject |
| 8 | **Analysis & Visualization** | ~8 min | Stats, plots, phylogenetics, synteny, format conversion |

**Total: ~63 minutes** covering the full pipeline from basics to NCBI submission.

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
  True or False: myconote-cli can only annotate fungal genomes.

  > false
  Correct! It supports fungi, plants, animals, insects, and protists.
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
