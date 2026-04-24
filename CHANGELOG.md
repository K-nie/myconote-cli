# Changelog

All notable changes to MycoNote-CLI are recorded here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
  `ndarray-stats`, `tar`, `zip`, `sha2`, `calamine`.
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
