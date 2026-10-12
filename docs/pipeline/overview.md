# Pipeline Overview

MycoNote-CLI organises genome annotation into seven sequential commands. Each command produces standard-format output files that feed directly into the next step.

```
┌──────┐  ┌──────┐  ┌───────┐  ┌─────────┐  ┌────────┐  ┌──────────┐  ┌────────┐
│ sort │→ │ mask │→ │ train │→ │ predict │→ │ update │→ │ annotate │→ │ submit │
└──────┘  └──────┘  └───────┘  └─────────┘  └────────┘  └──────────┘  └────────┘
```

| Step | Command | What it does | Key tools |
|------|---------|--------------|-----------|
| 1 | `sort` | Reorder + length-filter contigs (lossless by default; `--ncbi-clean` to rename) | Built-in |
| 2 | `mask` | Identify and soft-mask repetitive elements | RepeatModeler2, RepeatMasker |
| 3 | `train` | Train gene predictors on RNA evidence | PASA, Trinity, Augustus training |
| 4 | `predict` | Predict gene models ab initio + consensus | Augustus, GlimmerHMM, SNAP, GeneMark, miniprot |
| 5 | `update` | Refine models with RNA-seq | PASA update |
| 6 | `annotate` | Assign functional annotations | MMseqs2, eggNOG, dbCAN, InterProScan |
| 7 | `submit` | Prepare NCBI GenBank submission | table2asn, built-in validation |

You can start the pipeline at any step if you already have intermediate files. For example, if you have a masked genome and trained Augustus parameters, you can jump straight to `predict`.

---

## Multi-genome mode

For projects with multiple genomes, use `batch` to run the full pipeline across all of them:

```bash
myconote-cli batch genomes/ --kingdom fungi --threads 8
myconote-cli batch samples.tsv --condor   # HTCondor HPC submission
```

See [batch documentation](../analysis/batch.md) for details on sample sheets, HTCondor, and resume.

---

## AI-powered interpretation

After any stage, use `explain` to interpret the results:

```bash
myconote-cli explain predict              # full: rules + local LLM
myconote-cli explain predict --no-llm     # rules only
```

See [explain documentation](../analysis/explain.md) for details.

---

## Common flags

Pipeline commands take the input as a **positional argument** (not a `--genome`
flag) and write to a path given by `--output`/`-o`. Common flags:

| Flag | Description |
|------|-------------|
| `<input>` (positional) | Input genome FASTA / GFF3 for that stage |
| `--output`, `-o` | Output path (file or directory, per command) |
| `--threads` | CPU threads where the stage uses them (e.g. `predict`, `annotate`; `sort` has none) |
| `--kingdom` | `fungi` (default, primary use case). `plant` / `animal` / `insect` / `protist` accepted but experimental — see the multi-kingdom note on the home page. |
| `--species` / `--clade` | Override the Augustus species model (`predict`); `--clade` picks a clade-appropriate fungal default |

---

## Environment variables

| Variable | Description |
|----------|-------------|
| `NCBI_EMAIL` | Used by BLAST remote queries |
| `AUGUSTUS_CONFIG_PATH` | Override Augustus config directory |
| `EVM_HOME` | Path to EvidenceModeler installation |
| `MYCONOTE_CHAT_ENDPOINT` | Override Ollama endpoint for `explain` |
| `MYCONOTE_CHAT_MODEL` | Override LLM model for `explain` |
| `MYCONOTE_CHAT_API_KEY` | Optional API key for LLM endpoint |
