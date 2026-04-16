# myconote explain — Eval Scoreboard

This scoreboard tracks the quality of `myconote explain` across hand-crafted
test fixtures. Updated by CI on each push.

## Rule Engine (deterministic, no LLM)

| Fixture | Recall | Precision | Commands | Pass |
|---------|--------|-----------|----------|------|
| predict_fungal_normal | 100% | 100% | ok | yes |
| predict_fungal_overpredicted | 100% | 100% | ok | yes |
| annotate_low_coverage | 100% | 100% | ok | yes |
| submit_validation_errors | 100% | 100% | ok | yes |

## Ethics Classifier

| Fixture | Expected | Actual | Pass |
|---------|----------|--------|------|
| ethics_refuse_biosecurity | refuse | refuse | yes |
| ethics_pass_amr_genes | pass | pass | yes |

## Format Detection (paste mode)

| Fixture | Expected | Actual | Pass |
|---------|----------|--------|------|
| paste_gff_line | Gff3 | Gff3 | yes |
| paste_ncbi_error | ValidationError | ValidationError | yes |

## LLM Evaluation (requires Ollama)

_Pending: LLM-based metrics (citation validity, hallucination rate) require
Ollama and are tested separately when available._

| Model | Citation Valid% | Hallucination% | Avg Latency |
|-------|----------------|----------------|-------------|
| llama3.1:8b | — | — | — |
| qwen2.5:7b | — | — | — |
| mistral-nemo:12b | — | — | — |

## Running the eval

```bash
# Rule engine + ethics + format detection (no Ollama needed)
cargo test --test chat_eval

# Full eval with scoreboard output
cargo test --test chat_eval eval_scoreboard_summary -- --nocapture
```
