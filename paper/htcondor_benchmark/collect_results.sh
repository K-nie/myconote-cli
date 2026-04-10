#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# collect_results.sh
# Verify all HTCondor jobs completed, then aggregate into manuscript tables.
#
# Usage:
#   bash collect_results.sh
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BENCHMARK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${DATA_DIR:-$BENCHMARK_DIR/data}"
RESULTS_DIR="${RESULTS_DIR:-$BENCHMARK_DIR/results}"
CONFIG="$BENCHMARK_DIR/configs/genomes.tsv"

TOOLS=(myconote funannotate maker braker)
REPS=(1 2 3)

echo "═══════════════════════════════════════════════════════════"
echo "MycoNote-CLI Benchmark: Result Collection"
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "Results dir : $RESULTS_DIR"
echo ""

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "ERROR: Results directory not found: $RESULTS_DIR"
    exit 1
fi

# ── Step 1: Enumerate expected runs ──────────────────────────────────────
GENOMES=()
while IFS=$'\t' read -r id rest; do
    [[ "$id" == "id" ]] && continue
    [[ -z "$id" ]] && continue
    GENOMES+=("$id")
done < "$CONFIG"

EXPECTED=$(( ${#TOOLS[@]} * ${#GENOMES[@]} * ${#REPS[@]} ))
echo "Expected runs: $EXPECTED (${#TOOLS[@]} tools x ${#GENOMES[@]} genomes x ${#REPS[@]} replicates)"
echo ""

# ── Step 2: Check completion status ──────────────────────────────────────
echo "── Step 2: Checking completion status ──"
declare -A tool_complete
declare -A tool_missing

total_complete=0
total_missing=0
missing_list=()

for tool in "${TOOLS[@]}"; do
    tool_complete[$tool]=0
    tool_missing[$tool]=0
    for genome in "${GENOMES[@]}"; do
        for rep in "${REPS[@]}"; do
            metrics="$RESULTS_DIR/$tool/$genome/rep$rep/metrics.json"
            if [[ -f "$metrics" ]]; then
                tool_complete[$tool]=$(( ${tool_complete[$tool]} + 1 ))
                total_complete=$((total_complete + 1))
            else
                tool_missing[$tool]=$(( ${tool_missing[$tool]} + 1 ))
                total_missing=$((total_missing + 1))
                missing_list+=("$tool/$genome/rep$rep")
            fi
        done
    done
done

for tool in "${TOOLS[@]}"; do
    printf "  %-12s : %2d complete, %2d missing\n" \
        "$tool" "${tool_complete[$tool]}" "${tool_missing[$tool]}"
done
echo ""
printf "  TOTAL        : %2d / %d complete\n" "$total_complete" "$EXPECTED"
echo ""

# ── Step 3: Report missing runs ──────────────────────────────────────────
if [[ "$total_missing" -gt 0 ]]; then
    echo "── Missing runs ──"
    for m in "${missing_list[@]}"; do
        echo "  $m"
    done
    echo ""
    echo "Missing runs can be re-submitted manually via condor_submit, or"
    echo "you can proceed with aggregation on the partial set (not recommended"
    echo "for the final manuscript)."
    echo ""
    read -r -p "Proceed with partial aggregation? [y/N] " reply
    if [[ ! "$reply" =~ ^[Yy]$ ]]; then
        echo "Aborting. Re-run after all jobs complete."
        exit 1
    fi
fi

# ── Step 4: Run aggregation ──────────────────────────────────────────────
echo "── Step 4: Aggregating metrics ──"
python3 "$BENCHMARK_DIR/scripts/aggregate_metrics.py" "$RESULTS_DIR"
echo ""

# ── Step 5: Summary of outputs ───────────────────────────────────────────
AGG_DIR="$RESULTS_DIR/aggregated"
echo "── Aggregated outputs ──"
for f in per_run_metrics.tsv per_tool_summary.tsv comparison_table.tsv statistical_tests.tsv; do
    if [[ -f "$AGG_DIR/$f" ]]; then
        lines=$(wc -l < "$AGG_DIR/$f")
        printf "  %-28s (%d lines)\n" "$f" "$lines"
    else
        printf "  %-28s MISSING\n" "$f"
    fi
done
echo ""

echo "═══════════════════════════════════════════════════════════"
echo "Collection complete."
echo ""
echo "Manuscript-ready tables: $AGG_DIR"
echo ""
echo "Recommended next steps:"
echo "  1. Inspect comparison_table.tsv for the Results section"
echo "  2. Inspect statistical_tests.tsv for p-values and effect sizes"
echo "  3. Feed per_run_metrics.tsv into R/Python for figure generation"
echo "═══════════════════════════════════════════════════════════"
