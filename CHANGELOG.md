# Changelog

All notable changes to MycoNote-CLI are recorded here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **`de-template` subcommand** — generate a self-contained R script
  that runs DESeq2 differential-expression analysis on `quant`
  output. The tool writes the script; the user runs it with
  `Rscript`. Implementation is Option 1D from
  `scratch/rnaseq_spec_decisions.md` — no R runtime dependency on
  our side, we just emit tximport + DESeq2 + apeglm boilerplate
  with the user's design and contrasts substituted in.
  - Prerequisite alerting in three layers: CLI `--help` block,
    stderr warning when `Rscript` isn't on PATH, and a
    `requireNamespace()` guard at the top of the emitted R script
    with a clear `BiocManager::install(...)` hint on failure.
  - Input validation before writing the script: factor column
    existence, level values present, design-formula sanitization
    against `;`, backticks, newlines, `system(...)`, `eval(...)`.
  - Multiple contrasts per invocation via repeated `--contrast`.
  - Emits per-contrast TSV (sorted by padj) + MA plot + volcano
    plot with FDR and |LFC| threshold lines.
  - `apeglm` LFC shrinkage on when the package is installed, with
    a graceful fallback to raw LFCs when it isn't.
  - 23 unit tests + 7 CI-safe integration tests + 1 ignored
    live `Rscript parse()` syntactic validation (verified against
    R 4.4.2).
- **`fetch-rna` subcommand** — download public RNA-seq FASTQs by
  SRA/ENA accession. ENA REST is the default backend (no credentials,
  no `vdb-config`); sra-toolkit (`prefetch` + `fasterq-dump`) is an
  opt-in fallback for runs ENA hasn't mirrored. Accepts run IDs,
  study IDs, project IDs, sample IDs, experiment IDs, or a text
  file with one accession per line. Streams files with on-the-fly
  MD5 verification and retry on transient failures. Emits a
  `samples.tsv` pre-populated for `quant`.
  - Dep added: `md-5 = "0.10"` for ENA MD5 checksum verification.
  - `sra-tools` registered in `install` / `check` with
    `used_by="fetch-rna"` (optional).
- **`quant` subcommand** — RNA-seq expression quantification against
  an annotated fungal genome using salmon with a decoy-aware index.
  Reads → fastp (QC/trim) → salmon quant → wide `counts.tsv` +
  `tpm.tsv` matrices + per-sample `quant.sf` (tximport-native) + a
  fixed-schema `quant_bundle.json` reproducibility manifest. DE
  analysis (DESeq2 / edgeR / limma) stays in R.
  - Sample-sheet driven: TSV with sample_id, fastq_r1 (paired/SE),
    optional fastq_r2, condition, strandedness, batch. Unknown
    columns are preserved verbatim.
  - SHA256-keyed salmon-index cache — reruns on the same inputs
    skip the expensive rebuild. Cache dir follows precedence
    `--index-cache` > `MYCONOTE_INDEX_CACHE` > `$XDG_CACHE_HOME` >
    `~/.cache/myconote/` > `~/.myconote/`; no silent fallback to
    `/tmp` (that would break bundle provenance).
  - `salmon` 1.10+ and `fastp` 0.23+ registered in `install` / `check`
    with `used_by = "quant"`.
  - Parallelism: ships with serial `--jobs 1` in 0.3.0 pending an
    A. niger 6-sample benchmark on two machine shapes (see
    `scratch/rnaseq_spec_decisions.md` §5).
  - `explain quant` integrates the new stage with thresholds for
    mapping rate, fastp Q30, duplication rate, insert-size peak.
- **`convert --to cds`** — extract spliced nucleotide CDS per mRNA
  (concatenated child CDS, reverse-complemented on `-` strand, phase
  offset honored). Feeds directly into `quant` as the salmon target
  transcriptome.

### Changed
- **Scope narrowed to fungi.** Tool description, README, mkdocs, and
  `learn` tutorial now lead with fungal genome annotation. Non-fungal
  eukaryotes still work but are no longer a first-class target.

### Removed
- **`phylogeny` subcommand** — IQ-TREE wrapper has been dropped. Users
  who need an ML tree should align compare's
  `Single_Copy_Orthologue_Sequences/` with MAFFT and run IQ-TREE or
  RAxML-NG externally. This also drops the internal `align` module
  (MAFFT driver + supermatrix builder) and `compare --species-tree`.
- **`place` subcommand** and the entire **Y1000+ reference bundle**
  (`setup --y1000plus`, `stats --benchmark y1000plus`,
  `--annotated` KEGG-KO percentile rank). The 1,154-yeast Y1000+ work
  is reference data; hosting the downloader, Excel/TSV parsers, KO
  Jaccard, codon-table/phenotype/niche predictors, and tarball
  extractors inside an annotation tool was scope creep. The Opulente
  et al. (2024) dataset remains publicly accessible for external use.
- Dependencies dropped now that the above modules are gone: `ndarray`,
  `ndarray-stats`, `tar`, `zip`, `calamine`. (`sha2` was dropped here
  then re-added for the `quant` index cache — see Added above.)
- External-tool registrations removed from `install` / `check`:
  `iqtree`, `mafft`, `muscle`.

### Chat / `explain` validator
- Replaced the citation-tag validator (rarely fired) with a
  subcommand validator: lines that invoke a non-existent
  `myconote-cli <subcommand>` (e.g. `myconote-cli repeatmasker`) are
  stripped from LLM output before the user sees them. The
  authoritative subcommand list lives in
  `src/chat/validator.rs::VALID_SUBCOMMANDS` and is mirrored into
  the system prompt in `src/chat/prompts.rs`.

## [0.2.0] — 2026-04-18

Major release: AI-assisted interpretation, multi-genome batch orchestration,
Y1000+ 1,154-yeast placement with functional predictions, and real
N-genome ortholog inference via OrthoFinder.

### Added
- **`explain` subcommand** — LLM-powered interpreter for every pipeline stage.
  - Local Ollama backend (default `llama3.1`, `--endpoint` configurable).
  - Rules-only `--no-llm` fallback; auto-detects available RAM to pick a model.
  - RAG over four source types: `knowledge:`, `paper:<doi>`, `data:<file>:<line>`, `rule:`.
  - Post-hoc citation validator strips any tag the model invents.
  - Pattern-based ethics gate on incoming queries.
  - `--trace` / `--dry-run` for transparency; `--paste` reads stdin (e.g. NCBI errors).
  - Reproducibility bundles written to `explain_*/` (gitignored).
  - Integration tests + eval harness with fixtures and scoreboard.
  - Paper corpus setup with manifest validation and BM25 retrieval.
  - Auto-suggests relevant `learn` lesson on first run.
- **`batch` subcommand** — annotate many genomes from a directory or sample
  sheet, with first-class HTCondor support.
- **`place` subcommand** — place a user genome into the Y1000+ 1,154-yeast
  reference bundle.
  - Codon-table advice for the top hit when the `codontable` subset is installed.
  - Carbon/nitrogen lifestyle prediction from the `metabolism` subset.
  - Growth-at-37 °C prediction from the `phenotypes` subset.
  - Ecological-niche prediction from the `environment` OWL subset.
- **`stats --benchmark`** — percentile-of-user gene-count and tRNA-count
  against the Y1000+ reference distribution.
- **`compare` (real)** — N-genome ortholog inference wrapping OrthoFinder
  (default `dendroblast`; MAFFT species-tree bridge available).
- **Ollama auto-setup** — install check, version probe, model pull.
- **Scientific disclaimer** — every `explain` output carries a
  "suggestive, not definitive" banner.

### Changed
- `compare` now wraps OrthoFinder properly instead of returning the stale
  synteny placeholder.
- `predict`: real per-contig Augustus parallelism.
- `predict`: SNAP ZOE lookup fixed.
- `predict`: `--glimmerhmm` / `--genemark` dispatch wired into
  `run_prediction`.
- `annotate`: `--eggnog` / `--cazyme` / `--secretome` / `--antismash`
  dispatch wired.
- `submit`: multi-exon CDS emitted as a single joined feature block (not N
  separate ones); `product` and `Dbxref` propagation fixed.
- `phylogeny`: path-like `--prefix` handled correctly; output parents auto-
  created.
- All subcommands: output-parent directories auto-created.
- Release profile hardening: LTO, single codegen unit, stripped, opt-level 3
  (already in v0.1.0, re-confirmed).

### Removed
- `plot` (circular/linear) and `view` (JBrowse2/UCSC/synteny) subcommands
  — their outputs fell short of publication-grade; users are now directed
  to external tools (Proksee, IGV, clinker) that already do this well.
  A pyGenomeViz synteny backend was prototyped during the cycle but not
  retained.
- Fake `compare` subcommand (redirects users to the real implementation).

### Fixed
- `compare`: `.fas` extension handling.
- `train`: model-completeness check.
- `stats`: per-chromosome accounting bug.
- Several dead modules and stale files removed.
- CI: formatting + new RUSTSEC advisory ignores; pre-push hook guards the
  CI budget.

### Repository hygiene
- Python bytecode caches gitignored.
- Temporary branch-cleanup files removed.

---

## [0.1.0] — 2026-04-09

Initial public release. Pipeline: `sort → mask → train → predict → update
→ annotate → submit`. Eight analysis commands (`stats`, `phylogeny`,
`compare`, `convert`, `clean`, `fix`) and five utilities (`install`,
`check`, `setup`, `species`, `learn`).

[0.2.0]: https://github.com/K-nie/myconote-cli/releases/tag/v0.2.0
[0.1.0]: https://github.com/K-nie/myconote-cli/releases/tag/v0.1.0
