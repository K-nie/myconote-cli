# Pipeline Overview

Myconote_CLI organises genome annotation into six sequential commands. Each command produces standard-format output files that feed directly into the next step.

```
┌─────────┐   ┌─────────┐   ┌─────────┐   ┌─────────────┐   ┌─────────┐   ┌──────────┐
│  sort   │ → │  mask   │ → │  train  │ → │   predict   │ → │ update  │ → │ annotate │
└─────────┘   └─────────┘   └─────────┘   └─────────────┘   └─────────┘   └──────────┘
```

| Step | Command | What it does | Key tools |
|------|---------|--------------|-----------|
| 1 | `sort` | Standardise scaffold names, filter short sequences | Built-in |
| 2 | `mask` | Identify and soft-mask repetitive elements | RepeatModeler2, RepeatMasker |
| 3 | `train` | Train gene predictors on RNA evidence | PASA, Trinity, Augustus training |
| 4 | `predict` | Predict gene models ab initio | Augustus, GlimmerHMM, SNAP, EVM |
| 5 | `update` | Refine models with RNA-seq | PASA update |
| 6 | `annotate` | Assign functional annotations | BLAST, eggNOG, dbCAN, InterProScan |

You can start the pipeline at any step if you already have intermediate files. For example, if you have a masked genome and trained Augustus parameters, you can jump straight to `predict`.

---

## Common flags

Most pipeline commands accept:

| Flag | Description |
|------|-------------|
| `--genome` | Input genome FASTA |
| `--out` | Output directory (created if absent) |
| `--threads` | CPU threads (default: 4) |
| `--kingdom` | `fungi` / `plant` / `animal` / `protist` (sets species defaults) |
| `--species` | Override Augustus species model |
| `--config` | Path to a TOML config file |

---

## Environment variables

| Variable | Description |
|----------|-------------|
| `NCBI_EMAIL` | Used by BLAST remote queries |
| `AUGUSTUS_CONFIG_PATH` | Override Augustus config directory |
| `EVM_HOME` | Path to EvidenceModeler installation |
