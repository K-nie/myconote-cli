# Changelog

All notable changes to MycoNote-CLI are recorded here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.1] — 2026-04-25

### Fixed
- **Audit pass over user-facing claims.** Comprehensive sweep across
  README, CHANGELOG, help text, learn lessons, and source after three
  rapid feature releases (v0.5.0 → v0.7.0). Full report in
  `docs/audit/v0.7.1_audit.md`. Specific corrections:
  - README "18 NCBI translation tables supported in total" → "25 NCBI
    translation tables supported in total (tables 1–6, 9–14, 16,
    21–31, 33)" — matches the data-driven registry in
    `src/annotate/genetic_code.rs::REGISTRY`.
  - README "8 lessons" → "20 lessons" with the new four-track
    catalogue description (foundation 1-8, RNA-seq 9-14, predictors
    15-16, quality 17-20).
  - `learn` lesson 5 narrative "myconote-cli supports 18 NCBI
    translation tables" → "25 NCBI translation tables — tables 1–6,
    9–14, 16, 21–31, and 33". Same drift, same fix.
  - `learn` lesson 1 final True/False reframed: "myconote-cli can only
    annotate fungal genomes" → "myconote-cli is focused primarily on
    fungal genomes" (matches the v0.5.0 fungi-narrow scope: defaults
    target fungi; non-fungal kingdom flags remain as experimental
    scaffold).
  - `learn` lesson 8 ("Analysis & Visualization") rewritten end-to-end
    from `plot` / `view` / `synteny` (none of which exist as
    subcommands; they were aspirational in v0.1.0 and removed before
    v0.2.0) to `stats` + `compare` + `convert` with explicit handoff
    guidance to IGV / Proksee / JBrowse2 / clinker as external
    visualization tools. Now titled "Stats, Comparative Genomics, and
    Conversion".
  - `docs/learn.md` rewritten to reflect the 20-lesson catalogue with
    four tracks; the previous version still listed "8 \| Analysis &
    Handoff \| ~8 min \| Stats, phylogenetics, Y1000+ placement, …"
    which referenced removed features.
  - `docs/paper/myconote_manuscript.md` annotated with a v0.7.1 audit
    TODO block at the top covering two count drifts ("fifteen
    functional annotation sources" — code has 11; "ships ~35 curated
    fungal Augustus species" — manifest is 49 since v0.7.0). Code is
    the source of truth; user revises framing.
- **Removed dead `src/cli/` clap scaffold.** `src/cli/commands.rs` and
  `src/cli/mod.rs` were unused pre-v0.2.0 clap-derived enum +
  argument structs (`Phylogeny(PhylogenyArgs)`, `Blast(BlastArgs)`,
  `Align(AlignArgs)`, `CompareArgs.build_tree`, `CompareArgs.plot_type`)
  for subcommands that were never wired up — main.rs uses manual arg
  parsing throughout. Deleted entire directory; removed `pub mod cli;`
  from `src/main.rs` and `src/lib.rs`. Build remains warning-clean.
  Without this removal, repo-wide greps for `phylogeny` / `Blast` /
  `Align` would keep returning false positives that obscure real
  drift in future audits.

### Changed
- **`learn` tutorial expanded from 8 to 20 lessons** (~127 additional
  minutes of content) covering the modern features added in v0.5.0
  through v0.7.0. Lessons 9-20 are completable in any order — the
  catalogue is split into four tracks:
  - Foundation track (1-8, sequential): pipeline basics through NCBI
    submission. Lessons 1, 5, and 8 cleaned up as part of the audit
    (see Fixed above).
  - RNA-seq track (9-14): `quant` (salmon decoy-aware index, fastp QC,
    SHA256 cache), `fetch-rna` (ENA REST default, MD5 verification,
    sample-sheet emission), `de-template` (tximport + DESeq2 +
    apeglm), `ase` (phased VCF, personalized transcriptome,
    unphased-het hard error), `ase-template` (per-sample-corrected
    binomial test — the most pedagogically valuable concept in the
    ASE stack), `go-template` (topGO Fisher's exact, BP/MF/CC).
  - Modern predictors track (15-16): `predict --use-braker` (BRAKER
    1/2/3 mode auto-detection, conflict guard, genetic-code
    forwarding), `predict --genemark-mode` (ES/ET/EP+/ETP+, ProtHint
    dependency).
  - Quality and reproducibility track (17-20): `clean --mode contigs`
    (purge_dups-style haplotig dedup), `setup --db augustus-fungi`
    (the 49-species manifest, AUGUSTUS_CONFIG_PATH), reproducibility
    manifests (`quant_bundle.json` / `ase_bundle.json` schema, SHA256
    input hashing, byte-identical re-runs from a published bundle),
    NCBI codon tables (Candida CTG-clade focus, the
    silent-mistranslation failure mode under Standard, Table 12
    propagation through to BRAKER's `--translation_table`).
- Each new lesson follows the existing pattern (objective →
  narrative → command → verify → next-step nudge), routes to
  `myconote-cli check` / `install` first when external tools are
  needed, and offers `skip` for users without the dependencies.
- 10 new lib tests in `src/learn/lessons.rs::tests` guard the
  catalogue: count (must be 20), title uniqueness, every lesson has
  a verification terminal step, lesson durations within
  `[5, 20]` minutes, RNA-seq / predictors / quality tracks present
  by keyword, name-based lookup integrity for new lessons, and a
  regression test that no lesson body re-introduces removed
  subcommands (`myconote-cli phylogeny`, `myconote-cli place`,
  `Y1000+`).

## [0.7.0] — 2026-04-24

### Added
- **`update --kallisto`** — funannotate-style abundance pre-filter
  for the `update` step. When set, MycoNote-CLI builds a kallisto
  index over the transcript set, runs `kallisto quant` per
  RNA-seq sample, drops transcripts whose TPM falls below
  `--kallisto-min-tpm` (default **1.0**), and passes only the
  surviving transcript IDs through to PASA's update step. The
  filtered ID list is written to
  `update_out/pasa_kallisto_filter.txt` for auditability. Pass
  semantics are "in any sample" (max TPM ≥ threshold) so a
  transcript that's expressed in one condition isn't silently
  dropped because it's quiet in another. New `update::kallisto`
  module owns the driver: header-row-aware abundance.tsv parser
  that tolerates extra columns and rejects missing required
  columns with named errors, deterministic per-transcript
  aggregation via BTreeMap, and a `KallistoFilter` config struct.
  Kallisto CLI verified against pachterlab.github.io/kallisto/manual
  on 2026-04-24 (kallisto 0.50.0):
  `kallisto index -i out.idx <fa>`,
  `kallisto quant -i idx -o out r1 r2`. When `--kallisto` is set
  but `kallisto` is missing from PATH, the run aborts with a
  pointer at `myconote-cli install update`; when set without an
  RNA-seq input or with only `--rna-bam` (which kallisto cannot
  re-quantify), the run aborts with an actionable message. No
  silent fallback. `kallisto` registered in install/check
  catalogues with `conda_pkg=kallisto`, `conda_chan=bioconda`,
  `used_by=update`. 9 unit tests in src/update/kallisto.rs cover
  the parser and aggregation; 4 in tests/unit_tests.rs drive the
  CLI parser path (valid/garbage/negative `--kallisto-min-tpm`,
  `--kallisto` without RNA-seq).
- **`compare --html`** — self-contained interactive HTML report
  alongside the existing TSV outputs. Single file at
  `compare_out/report.html`, opens correctly from `file://` with
  no network access — CSS and JS are inlined via `include_str!`,
  no CDN, no Google Fonts, no JS frameworks. Carries five panels:
  pan-genome shape (table + stacked-bar SVG), per-genome stats
  (sortable table), ortholog table (sortable + text-filterable,
  capped at 5,000 rows with a pointer at the full TSV), species
  tree (inline SVG with branch lengths to scale), and a
  reproducibility footer with tool version + git SHA + run
  command + SHA256 of every input GFF/FASTA. Templating via
  `tinytemplate` 1.2 (runtime, no proc-macros — picked over
  askama for footprint and over handlebars for dep-tree size).
  Newick parser is 56 lines hand-rolled, rejects bare strings
  without `(` or `;` so truncated files don't silently parse as
  single-leaf trees. Every user-supplied gene ID, orthogroup ID,
  and genome name (including SVG leaf labels) runs through an
  HTML-escape pass that catches all five XSS-relevant characters,
  defending against malicious input from upstream pipelines.
  Client-side JS is 70 lines of vanilla — no React, no jQuery. 12
  unit tests cover the escape helper, full report rendering,
  Newick parsing, SVG emission, and the per-genome stats
  arithmetic; the XSS guard is tested with a literal
  `<script>alert('pwned')</script>` payload as a gene ID.
- **`setup --db augustus-fungi`** expansion — curated fungal
  Augustus species manifest grows from 35 to 49 entries by
  filtering the upstream `Gaius-Augustus/Augustus/config/species/`
  directory (167 dirs, queried 2026-04-24) against a curated
  fungal-genus catalogue plus the bare-genus aliases Augustus
  historically ships (`saccharomyces`, `fusarium`, `cryptococcus`,
  `coprinus`, `histoplasma`, `neurospora`, `pneumocystis`,
  `ustilago`, `pchrysosporium`, `anidulans`). Every name in the
  expanded manifest verifies against the upstream GitHub API. The
  spec hoped for ~100 entries to match funannotate's bundled tree,
  but `Augustus/config/species/` only ships ~50 fungal entries —
  funannotate downloads beyond what's in the upstream repo. We
  don't fabricate names that aren't upstream. The species-list
  test bound widens to 45..=60 to accommodate the bare-genus
  aliases plus `verticillium_longisporum1` and
  `cryptococcus_neoformans_neoformans_JEC21`. Catalogue
  description and `setup --db augustus-fungi` help text bump from
  "~30" to "~50".

### Fixed
- **`myconote-cli --help`** now lists the `explain` subcommand
  under a new "AI-powered interpretation" section between Analysis
  and Utility, mirroring the README structure. The arm itself
  worked since 0.6.0 but users discovering MycoNote-CLI through
  `--help` never saw the LLM interpreter exists.

## [0.6.0] — 2026-04-24

### Added
- **`predict --genemark-mode <es|et|ep|etp>`** — adds the two
  protein-evidence-guided GeneMark variants that were previously
  unreachable (EP+ and ETP+) and replaces the dual `--genemark` /
  `--genemark-hints` flags with a single mode selector. EP+
  (`gmes_petap.pl --EP`) is the highest-impact addition for novel
  CTG-clade fungi without RNA-seq because the alternative-yeast
  nuclear code makes published predictors error-prone, and protein
  evidence anchors the training. EP/ETP+ each call ProtHint
  internally to convert the genome + protein FASTA into the GFF
  hint files that GeneMark consumes; ETP+ additionally takes the
  RNA-seq intron hints. ETP requires both `--genemark-hints` and
  `--protein-fasta` and fails loudly if either is missing — no
  silent fallback. The legacy flags still work and are translated:
  `--genemark` (alone) → ES; `--genemark` + hints → ET; hints alone
  → ET. EP/ETP have no legacy alias because they did not exist in
  earlier releases. `prothint.py` registered in `install` / `check`
  with `manual_note` flagging that ProtHint has no standalone
  bioconda recipe (verified 2026-04-24 against bioconda osx-64 and
  noarch — `PackagesNotFoundError` for both `prothint` and the
  `prothint*` glob); ProtHint ships bundled with the GeneMark-ES
  installer tarball under `<install>/ProtHint/bin/` and with the
  bioconda `braker3` package. 11 unit tests cover mode-string
  parsing (`es|et|ep|etp` plus the publication-style `EP+`/`ETP+`
  aliases, case-insensitive), the legacy-alias resolution table,
  and the canonical long-name strings used in user-facing logs.
  The end-to-end live ETP test is gated under `#[ignore]` like
  `de-template`'s `Rscript parse()` check, because actually
  running `gmes_petap.pl --ETP` needs a Georgia-Tech-licensed
  GeneMark binary plus ProtHint together.
- **`predict --use-braker`** — BRAKER 1/2/3 as a first-class
  predictor. When set, MycoNote-CLI does not run Augustus / SNAP /
  GlimmerHMM / GeneMark / miniprot itself or the EVM consensus
  stage; `braker.pl` is invoked end-to-end and `braker.gff3` is
  published as `consensus.gff3` so the downstream `update`,
  `annotate`, `submit` stages consume it unchanged. New
  `src/predict/braker.rs` exposes a `BrakerMode` enum
  (`Braker1` / `Braker2` / `Braker3`) that maps to BRAKER's
  `--esmode` / `--epmode` / `--etpmode` switches; flag names
  verified against `Gaius-Augustus/BRAKER` `scripts/braker.pl`
  GetOptions on 2026-04-24. Mode auto-detects from inputs:
  RNA-only → BRAKER1, protein-only → BRAKER2, both → BRAKER3.
  `--braker-mode <1|2|3>` overrides auto-detection;
  `--braker-rna-bam <file>` is repeatable (and accepts a
  comma-separated list); `--braker-proteins <fa>` carries the
  protein database. New `--genetic-code <n>` flag forwards through
  to BRAKER as `--translation_table`, wiring the *Candida* CTG
  code (12) and ciliate-style codes (6) end-to-end through
  Augustus + GeneMark inside BRAKER. New
  `check_braker_conflicts` rejects `--use-braker` combined with
  any of the standard predictor flags (`--genemark-mode`,
  `--genemark`, `--genemark-hints`, `--protein-evidence`,
  `--protein-fasta`, `--glimmerhmm`) and lists every offending
  flag in a single error so the fix is one round-trip, not
  whack-a-mole. `--no-snap` is intentionally not a conflict
  (disabling a predictor that wouldn't run anyway is a no-op).
  `braker.pl` registered in `install` / `check` with
  `conda_pkg=braker3` (verified 2026-04-24: bioconda osx-64 lists
  `braker3` 3.0.3 → 3.0.8; the binary it provides is `braker.pl`,
  not `braker3.pl`). 25 tests added (10 in
  `src/predict/braker.rs`, 15 in `tests/unit_tests.rs`): mode-string
  parsing (numeric + named), the mode-flag string contract against
  upstream BRAKER CLI, the four-cell auto-detection truth table,
  default-config invariants, and the eight conflict cases. Live
  BRAKER smoke test gated under `#[ignore]`.

### Documentation
- **`docs/pipeline/predict.md`** rewritten from a stub. New
  "Choosing a GeneMark mode" subsection lays out the four GeneMark
  modes side by side with their inputs, quality tier, and when to
  reach for each, with an OrthoDB-fungi example for EP+. New
  "Using BRAKER for maximum accuracy" subsection covers the
  auto-detection table, the genetic-code forwarding, and the
  full dependency chain (braker3 + Augustus +
  `AUGUSTUS_CONFIG_PATH` + licensed GeneMark + ProtHint).

## [0.5.1] — 2026-04-25

### Added
- **`clean --mode contigs`** — extends the existing GFF3-cleanup
  subcommand with a FASTA mode that runs `minimap2 -X` self-alignment
  on a genome assembly and drops contigs whose entire span is covered
  by a longer contig at high identity. Useful for purging
  near-duplicate haplotigs from draft assemblies (purge_dups /
  purge_haplotigs territory). Defaults: `--coverage 0.95
  --identity 0.95`. Emits a cleaned FASTA plus a TSV report listing
  each dropped contig and the contig that subsumed it. The PAF
  intermediate is removed on success and kept on failure for
  debugging. Mutual subsumption is broken deterministically (longer
  wins, lex-smaller on ties). `clean` is now registered as a `used_by`
  consumer for `minimap2` in both `check` and `install` catalogues.
  11 unit tests cover the FASTA reader, PAF parser, drop decisions,
  and FASTA + report writers.
- **`setup --db augustus-fungi`** — curated bundle of ~35 fungal
  Augustus species fetched from the upstream GitHub mirror and
  installed into `~/.myconote/augustus_config/species/<name>/` so a
  fresh MycoNote-CLI install can predict on most fungi out of the
  box. List spans Saccharomycotina, Taphrinomycotina,
  Pezizomycotina, Basidiomycota, Mucoromycota, and Microsporidia;
  every name matches the upstream `config/species/<name>/` directory
  verbatim. Six config files per species
  (parameters/metapars/{exon,intron,igenic}_probs/weightmatrix). A
  new `--dry-run` flag prints what would be downloaded without
  fetching, so the curated list can be audited; a new singular `--db
  <name>` flag complements the existing `--dbs` plural form. After
  install, the user is told to set
  `AUGUSTUS_CONFIG_PATH=$HOME/.myconote/augustus_config` so Augustus
  picks the configs up. 9 unit tests cover catalogue presence, list
  invariants (size, uniqueness, breadth across major clades), URL
  construction, and the `species_present` completeness check.
- **`go-template` subcommand** — same Option 1D pattern as
  `de-template` / `ase-template`. Takes a `de-template` results TSV
  plus an `annotate` `annotations.tsv` (carrying a `go_terms` column)
  and emits a self-contained R script that builds the gene → GO
  mapping, defines the universe + foreground from the DE results, and
  runs topGO's classic Fisher's exact test for each ontology in
  BP / MF / CC. Outputs one TSV per ontology (sorted by classic
  Fisher *p*, BH-adjusted within ontology) plus a combined dot plot
  PDF. Three-layer prereq alerting (CLI help / runtime stderr /
  in-script `requireNamespace("topGO")` guard) mirrors the other two
  template subcommands. Configurable join columns
  (`--de-id-col` / `--ann-id-col`) and GO-column / separator so the
  same template handles transcript-level DE joined to locus-tag
  annotations or any other shared identifier. 18 unit tests cover
  argv parsing, ontology choices, render embedding, R-string
  escaping, and the rfc3339 helper. New
  `docs/analysis/go-template.md` at parity with the other two
  template docs; mkdocs nav and README analysis-commands table
  updated.

## [0.5.0] — 2026-04-24

### Added
- **`ase` subcommand** — allele-specific expression for heterozygous,
  hybrid, or polyploid fungal genomes. Builds personalized
  transcriptomes per haplotype from a phased VCF, then quantifies
  each sample against each haplotype with salmon. Emits
  transcript × (sample.haplotype) counts + TPM matrices, a
  per-transcript informativeness summary, variants-applied and
  variants-skipped audit TSVs, and an `ase_bundle.json`
  reproducibility manifest (extends the `quant_bundle.json` schema
  with per-haplotype transcriptome SHA256s, variant-application
  counts by category, and per-sample per-haplotype mapping rates +
  asymmetry flag). Full design in `scratch/ase_spec.md`; all 10
  open decisions locked on 2026-04-24 before implementation.
  - Custom phased-VCF parser (`src/ase/vcf.rs`) — zero new crates.
    Unphased heterozygous sites are a hard error with line number,
    not a silent skip. Supports `.vcf` and `.vcf.gz`. Handles SNVs,
    MNPs, and indels up to `--max-indel-size` (default 50 bp).
  - Personalization (`src/ase/personalize.rs`) — strand-aware
    variant → CDS mapping, reverse-complement on `-` strand,
    in-cis variant-overlap rejection, exon-boundary-spanning
    indels skipped. Nine documented skip categories; every variant
    lands on one side of the applied/skipped divide with a reason.
  - Per-sample-per-haplotype salmon driver (`src/ase/quant.rs`) —
    reuses `quant`'s SHA256-keyed index cache (different CDS FASTA
    → different cache slot automatically), one fastp per sample
    shared across both salmon runs, mapping-rate asymmetry
    detector (default threshold 5 %).
  - Merge + summary (`src/ase/merge.rs`) — long-to-wide join with
    deterministic column ordering (`<sample>.<hap>`), per-transcript
    informativeness flag computed from variant counts, max
    read-count asymmetry across samples.
  - Requires `salmon` (≥1.10) and `fastp` (≥0.23) on PATH; both
    registered in `install` / `check` already via `quant`.
  - 54 unit tests covering VCF parsing, variant application,
    asymmetry detection, merge logic, bundle round-tripping, and
    argument parsing.
- **`ase-template` subcommand** — companion R-script generator for
  `ase` output, same Option 1D pattern as `de-template`. Emits a
  base-R script (no Bioconductor) that runs a per-transcript
  binomial exact test on `(count_hap0, count_hap1)` with a
  sample-specific null proportion derived from total hap0:hap1
  library-size ratio (corrects for global mapping-rate asymmetry).
  BH-adjusts per sample. Outputs a long-format results TSV and a
  per-sample imbalance-histogram PDF. Handles column parsing that
  allows dots in sample IDs, filters to informative transcripts by
  default (overridable via `--include-uninformative`), and skips
  low-count rows below `--min-reads` (default 20) with a logged
  reason. 21 unit tests.
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

[0.7.1]: https://github.com/K-nie/myconote-cli/releases/tag/v0.7.1
[0.7.0]: https://github.com/K-nie/myconote-cli/releases/tag/v0.7.0
[0.6.0]: https://github.com/K-nie/myconote-cli/releases/tag/v0.6.0
[0.5.1]: https://github.com/K-nie/myconote-cli/releases/tag/v0.5.1
[0.5.0]: https://github.com/K-nie/myconote-cli/releases/tag/v0.5.0
[0.2.0]: https://github.com/K-nie/myconote-cli/releases/tag/v0.2.0
[0.1.0]: https://github.com/K-nie/myconote-cli/releases/tag/v0.1.0
