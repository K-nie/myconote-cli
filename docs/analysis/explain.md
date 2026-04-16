# myconote explain

AI-powered interpretation of genome annotation pipeline results. Uses a local LLM via Ollama -- no data leaves your machine.

```bash
myconote-cli explain <stage> [options]
myconote-cli explain --paste [options]
```

## How it works

```
User input
    |
[Ethics classifier] -- refuse biosecurity/fabrication/embargo requests
    |
[Context builder]   -- scan stage output directory, build StageContext
    |
[Rule engine]       -- deterministic findings from TOML-defined rules
    |
[BM25 retrieval]    -- search knowledge base + paper corpus
    |
[Command recommender] -- match findings to tool catalog
    |
  --no-llm? --> print findings + commands, exit
    |
[Prompt assembler]  -- system prompt + findings + retrieved context
    |
[Local LLM (Ollama)] -- generate interpretation
    |
[Citation validator] -- strip ungrounded claims
    |
[Renderer]          -- short / verbose / trace output
    |
[Bundle writer]     -- save reproducibility bundle
```

## Stages

| Stage | What `explain` interprets |
|-------|--------------------------|
| `sort` | Contig count, name mapping, filtered scaffolds |
| `mask` | Repeat content percentage, repeat family counts |
| `train` | Training report quality, model convergence |
| `predict` | Gene count, predictor agreement, expected ranges for kingdom |
| `update` | UTR extension stats, gene model changes |
| `annotate` | Functional coverage per source, missing annotations |
| `submit` | NCBI validation errors and warnings |

## Options

| Flag | Description |
|------|-------------|
| `<stage>` | Pipeline stage to interpret: sort, mask, train, predict, update, annotate, submit |
| `--dir <path>` | Stage output directory (default: auto-discovered) |
| `--paste` | Read from stdin instead of scanning a directory |
| `--stdin-format <fmt>` | Override format detection: auto, gff, vcf, log, error |
| `--model <name>` | Override LLM model |
| `--endpoint <url>` | Override Ollama endpoint |
| `--no-llm` | Skip LLM call, print rule-based findings only |
| `--verbose` | Show all findings with evidence |
| `--trace` | Show full prompt + retrieval scores |
| `--dry-run` | Print assembled prompt without calling LLM |

## Modes

### Stage mode (default)

Scans the output directory for a pipeline stage and interprets the results:

```bash
myconote-cli explain predict
myconote-cli explain predict --dir path/to/predict_out/
myconote-cli explain annotate --verbose
```

### Paste mode

Reads from stdin and auto-detects the format:

```bash
echo "ERROR: SEQ_FEAT.NoStop" | myconote-cli explain --paste
cat predict_summary.txt | myconote-cli explain --paste
head -5 genes.gff3 | myconote-cli explain --paste
```

Detected formats: GFF3, FASTA headers, NCBI validation errors, tool logs, TSV/CSV.

### Rules-only mode

Works without Ollama. Prints deterministic findings and command recommendations:

```bash
myconote-cli explain predict --no-llm
```

If Ollama is not running and `--no-llm` is not set, the tool prints rule-based findings with a note: "LLM unavailable -- showing rule-based findings only."

## Smart model selection

The tool auto-detects system RAM and selects the most capable model:

| Available RAM | Model | Quality |
|---------------|-------|---------|
| 48+ GB | llama3.3:70b-instruct-q4_K_M | Best quality |
| 24+ GB | qwen2.5:32b-instruct-q4_K_M | Excellent quality |
| 16+ GB | mistral-small:22b | Strong quality |
| 12+ GB | qwen2.5:14b | Good quality |
| 8+ GB | llama3.1:8b | Baseline |

Override with `--model <name>` or `MYCONOTE_CHAT_MODEL` environment variable.

## Configuration

### Config file (`~/.myconote/config.toml`)

```toml
[chat]
enabled = true                       # false disables explain entirely
endpoint = "http://localhost:11434"   # Ollama endpoint
model = "mistral-small:22b"          # override auto-selected model
timeout_s = 120                      # LLM request timeout
```

### Environment variables

| Variable | Description |
|----------|-------------|
| `MYCONOTE_CHAT_ENDPOINT` | Override Ollama endpoint |
| `MYCONOTE_CHAT_MODEL` | Override LLM model |
| `MYCONOTE_CHAT_API_KEY` | Optional API key |

Precedence: CLI flags > environment variables > config file > defaults.

## Setup

```bash
myconote-cli setup ollama        # install Ollama + pull best model
myconote-cli setup chat-corpus   # download Q1 papers for grounded citations (optional)
```

## Citation policy

- Every factual claim in the LLM output must cite a source: knowledge base, paper, user data, or rule
- Citations that cannot be resolved to a retrieved source are stripped automatically
- Paper citations come only from Q1 open-access journals (CC-BY licensed)
- The paper corpus is optional -- without it, explain works fully but prints no paper citations

## Ethics

The explain command refuses requests involving:

- **Biosecurity**: pathogen enhancement, select-agent engineering
- **Fabrication**: generating fake results, inventing statistics
- **Embargo violation**: circumventing NCBI embargo/hold rules
- **Academic misconduct**: ghostwriting methods sections, plagiarism

Legitimate biology questions (AMR genes, virulence factors, toxins, genome editing) are always answered.

## Reproducibility bundles

Every `explain` call writes a bundle to `./explain_<stage>_<timestamp>/`:

| File | Contents |
|------|----------|
| `prompt.txt` | The full prompt sent to the LLM |
| `context.json` | StageContext with artifacts and stats |
| `findings.json` | Deterministic rule engine findings |
| `response.txt` | LLM response (if applicable) |
| `env.json` | Model, endpoint, timestamp, myconote version |

## Scientific disclaimer

All interpretations are suggestive, not definitive. Findings must be independently verified in the context of your specific organism, assembly, and research goals. Do not rely solely on this tool for scientific conclusions or publication-ready claims.

## Examples

```bash
# Quick interpretation after prediction
myconote-cli explain predict

# Detailed view with all evidence
myconote-cli explain predict --verbose

# Debug the prompt (privacy audit)
myconote-cli explain predict --dry-run

# Analyze an NCBI validation error
echo "ERROR: SEQ_FEAT.NoStop at feature lcl|scaffold_1:gene-MYORG_001234" \
  | myconote-cli explain --paste

# Rules only on a specific directory
myconote-cli explain annotate --dir /data/project/annotate_out --no-llm
```
