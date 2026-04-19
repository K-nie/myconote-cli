---
marp: true
theme: default
paginate: true
size: 16:9
header: "MycoNote-CLI · GSA 2026 · Narh-Madey"
footer: "Hittinger Lab · Laboratory of Genetics · UW–Madison"
style: |
  section { font-family: 'Inter', sans-serif; }
  h1, h2 { color: #00573F; }
  code { background: #F3F4F6; padding: 2px 4px; border-radius: 3px; }
  .small { font-size: 70%; color: #555; }
  .cite { font-size: 55%; color: #666; font-style: italic; }
---

<!-- _paginate: false -->
<!-- _header: "" -->
<!-- _footer: "" -->

# MycoNote-CLI

## From raw contigs to NCBI-ready submission, in one binary

### With built-in functional placement and AI-assisted interpretation

---

**Benjamin Narh-Madey**
Hittinger Lab · Laboratory of Genetics
University of Wisconsin–Madison
`narhmadey@wisc.edu`

<span class="small">github.com/K-nie/myconote-cli · MIT · v0.2.0</span>

---

## The bottleneck has moved

Long-read sequencing is **cheap**. Annotation is not.

- **funannotate** — 3 fungi / day on a workstation; manual cleanup required
- **MAKER** — flexible but configuration-heavy
- **BRAKER** — RNA-seq gold standard, but a single pipeline step

**Result:** annotation, not sequencing, dictates paper timelines.

MycoNote-CLI is a single Rust binary that takes a eukaryotic genome assembly from raw contigs to NCBI-ready submission — with reproducibility baked in.

---

## One binary, seven commands

```
sort → mask → train → predict → update → annotate → submit
```

Each step is a separate subcommand. Stop anywhere. Resume anywhere. Standard formats (GFF3, FASTA, GenBank) throughout.

- **Rust** — blazing-fast parallelism, no Python dependency hell
- **Per-contig parallelism** — Augustus, repeat-masking, EVM all rayon-backed
- **JSON reproducibility reports** — tool versions, parameters, checksums
- **Docker + Singularity** — pinned dependency graphs for HPC and cloud

---

## Multi-predictor consensus, not single-predictor hope

```
         ┌─ Augustus   ─┐
         ├─ SNAP       ─┤
FASTA →  ├─ GlimmerHMM ─┼→ EVM consensus → GFF3 → annotate
         ├─ GeneMark   ─┤
         └─ miniprot   ─┘
              protein evidence
```

- Weighted evidence integration — one bad predictor does not sink the call
- Ploidy-aware allelic duplicate detection
- **18 NCBI genetic code tables** supported

---

## Why genetic-code correctness matters

| Table | Name | When it applies |
|------:|------|-----------------|
| 1  | Standard | Most eukaryotes (default) |
| 12 | Alternative Yeast Nuclear | **Candida CTG clade (CTG = Ser)** |
| 3  | Yeast Mitochondrial | Yeast mt genomes |

The CTG clade reassigned **CTG → Serine** ~170–270 Mya.

Annotate *C. albicans* with table 1 and every Ser encoded by CTG reads as Leu in silico → wrong Pfam hits, wrong domain architecture, NCBI validation errors.

<span class="cite">Santos & Tuite (1995) NAR · Krassowski et al. (2018) Nat Commun · Opulente et al. (2024) Science</span>

---

## `place` — Y1000+ functional placement

**Input:** your `annotate_out/annotations.tsv` (KEGG-KO per gene)
**Reference:** Y1000+ — 1,154 sequenced yeasts (Opulente et al. 2024)
**Method:** KEGG-KO Jaccard similarity → ranked neighbours

Auto-predictions when reference subsets are installed:

- **Codon table** — from top hit → prevents the CTG trap
- **C/N metabolic lifestyle** — weighted vote across neighbours
- **Growth at 37 °C** — thermotolerance / pathogenic potential
- **Ecological niche** — insect gut, flower nectar, rotting wood, …

<span class="cite">Opulente DA et al. (2024). Science 384(6694): eadj4503</span>

---

## `explain` — AI interpretation with scientific guardrails

```bash
myconote-cli explain predict
```

- **Local Ollama backend** — your data never leaves the machine
- **RAM-aware model selection** — 8 GB baseline → 48 GB top tier
- **Deterministic fallback** — `--no-llm` mode needs no model

**Guardrails** — where most LLM-bio tools fail:

- Every claim must cite `[knowledge:…]`, `[paper:<doi>]`, `[data:<file>:<line>]`, or `[rule:…]`
- **Post-hoc citation validator** strips tags the model invented
- **Ethics gate** — pattern-based refusal for out-of-scope queries
- `--trace` shows prompt, retrieval scores, validator decisions
- `--dry-run` assembles the prompt without calling the LLM

---

## `batch` — many genomes, one command

```bash
# From a sample sheet
myconote-cli batch samples.tsv --parallel 4

# HTCondor submit-file generation for HPC
myconote-cli batch genomes/ --condor --condor-mem 64G
```

- **Directory mode** or **sample sheet** (TSV with per-genome settings)
- **Resume** after partial failure — no recomputation
- **HTCondor first-class support** — generates `condor.sub` + per-genome wrappers
- Live dashboard for local parallel runs

---

## `compare` — real N-genome orthology

OrthoFinder wrapping — not a toy.

- **DendroBLAST** default (fast) or **MSA** mode (MAFFT + trimAl)
- **Supermatrix bridge into `phylogeny`** — concatenated single-copy orthologs → IQ-TREE species tree
- **Tiered genome-count caps** — 5 fungi / 3 Arabidopsis-class dicots / 2 maize-class monocots (because OrthoFinder scales as *n*²)
- Disambiguates duplicate genome names before feeding to OrthoFinder

```bash
myconote-cli compare --kingdom fungi proteomes/
```

---

## Validation — two real genomes

**Brettanomyces bruxellensis** (novel assembly, 30 contigs, 12.9 Mb)

| Step | Time | Result |
|------|-----:|--------|
| predict (Augustus) | ~3 min | 5,218 genes |
| annotate (MMseqs2) | ~1.5 min | 75% with product names |
| annotate (Pfam, hmmsearch) | ~7 min | 85.6% with domains |
| NCBI validation | <1 s | **PASSED (0 errors)** |

**Candida tropicalis** — CTG clade; runs cleanly with `--genetic-code 12` auto-suggested by `place`.

---

## Reproducibility, built in

- **Docker + Singularity** — pinned tool versions, single-command deployment
- **CITATION.cff** — GitHub renders a one-click citation
- **JSON run reports** — every parameter, tool version, checksum
- **Semantic versioning** — v0.2.0 changelog explicit about prototypes that did not ship
- **MIT licence** — reuse freely with attribution
- **CI** — GitHub Actions (fmt + clippy + test + cargo-audit)

<span class="small">github.com/K-nie/myconote-cli/releases/tag/v0.2.0</span>

---

## What MycoNote-CLI does NOT do

Honest limitations matter for a methods paper.

- Functional assignments are **hypotheses**, not facts — a 30 % identity MMseqs2 hit is a best guess
- Gene prediction still makes mistakes, especially at UTRs and in repeat-rich regions
- LLM outputs are **suggestive, not definitive** — verify against your organism and aims
- v0.2.0 has **weeks of community testing**, not years — validate critical calls against a second tool
- Comprehensive 8-genome benchmark vs funannotate / MAKER / BRAKER is *planned*, not yet run

---

## Roadmap

- **Comprehensive benchmark** — 8 reference genomes × 4 tools × 3 replicates on the GLBRC HTCondor pool
- **Interactive web UI** — upload FASTA, watch pipeline progress, download GenBank
- **More Augustus training species** — expanding beyond the default fungal set
- **Integration with existing community pipelines** where it adds value, not replaces

---

## Take-home

1. **One binary**, seven commands — annotation without the dependency labyrinth
2. **Y1000+ placement** catches silent errors (codon table, lifestyle mismatch) before they ruin a paper
3. **AI interpretation with real guardrails** — citation validation, ethics gate, deterministic rules fallback — not another hallucination engine
4. **MIT, reproducible, open** — clone, cite, contribute

---

## Acknowledgements

- **Chris Todd Hittinger** and the **Hittinger Lab**
- **Laboratory of Genetics**, University of Wisconsin–Madison
- **GLBRC** — HTCondor pool access
- **Dana Opulente** and the **Y1000+** consortium for the 1,154-yeast reference
- The OrthoFinder, IQ-TREE, Augustus, SNAP, EvidenceModeler, and PASA communities

---

<!-- _paginate: false -->

## Thank you

**github.com/K-nie/myconote-cli**

```bash
curl -fsSL https://raw.githubusercontent.com/K-nie/myconote-cli/main/quick-install.sh | bash
```

Benjamin Narh-Madey — `narhmadey@wisc.edu`

<span class="cite">Narh-Madey B. (2026). MycoNote-CLI v0.2.0. https://github.com/K-nie/myconote-cli</span>
