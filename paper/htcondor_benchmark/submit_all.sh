#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# submit_all.sh
# Submit benchmark jobs to HTCondor.
#
# Usage:
#   bash submit_all.sh                  # all 4 tools (72 jobs)
#   bash submit_all.sh --myconote-only  # MycoNote-CLI only (18 jobs)
#   bash submit_all.sh --tools m,f      # comma list: m=myconote f=funannotate
#                                         #            b=braker    k=maker
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BENCHMARK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${DATA_DIR:-$BENCHMARK_DIR/data}"
RESULTS_DIR="${RESULTS_DIR:-$BENCHMARK_DIR/results}"

# ── Parse tool selection ──────────────────────────────────────────────────
TOOLS_TO_SUBMIT="myconote,funannotate,braker,maker"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --myconote-only)
            TOOLS_TO_SUBMIT="myconote"
            shift
            ;;
        --tools)
            # Expand short codes: m=myconote, f=funannotate, b=braker, k=maker
            raw="$2"
            expanded=""
            IFS=',' read -ra codes <<< "$raw"
            for c in "${codes[@]}"; do
                case "$c" in
                    m|myconote)     expanded+=",myconote" ;;
                    f|funannotate)  expanded+=",funannotate" ;;
                    b|braker)       expanded+=",braker" ;;
                    k|maker)        expanded+=",maker" ;;
                    *) echo "ERROR: unknown tool '$c'" >&2; exit 1 ;;
                esac
            done
            TOOLS_TO_SUBMIT="${expanded#,}"
            shift 2
            ;;
        -h|--help)
            sed -n '2,14p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument '$1' (try --help)" >&2
            exit 1
            ;;
    esac
done

# Make these available to HTCondor submit files
export BENCHMARK_DIR
export DATA_DIR
export RESULTS_DIR

echo "═══════════════════════════════════════════════════════════"
echo "Submitting MycoNote-CLI Benchmark to HTCondor"
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "Benchmark dir : $BENCHMARK_DIR"
echo "Data dir      : $DATA_DIR"
echo "Results dir   : $RESULTS_DIR"
echo "Tools         : $TOOLS_TO_SUBMIT"
echo ""

# Verify HTCondor is available
if ! command -v condor_submit &>/dev/null; then
    echo "ERROR: condor_submit not found. Are you on a submit node?"
    exit 1
fi

submit_tool() {
    local tool="$1"
    local job_count="$2"
    local label="$3"
    echo "── Submitting $label jobs ($job_count) ──"
    condor_submit \
        -append "BENCHMARK_DIR = $BENCHMARK_DIR" \
        -append "DATA_DIR = $DATA_DIR" \
        -append "RESULTS_DIR = $RESULTS_DIR" \
        "$BENCHMARK_DIR/jobs/${tool}.sub"
    echo ""
}

# Submit in priority order (fastest / most informative first)
total_jobs=0
for tool in myconote funannotate braker maker; do
    if [[ ",$TOOLS_TO_SUBMIT," == *",$tool,"* ]]; then
        case "$tool" in
            myconote)    submit_tool myconote 18 "MycoNote-CLI" ;;
            funannotate) submit_tool funannotate 18 "funannotate" ;;
            braker)      submit_tool braker 18 "BRAKER" ;;
            maker)       submit_tool maker 18 "MAKER" ;;
        esac
        total_jobs=$((total_jobs + 18))
    fi
done

echo "═══════════════════════════════════════════════════════════"
echo "$total_jobs jobs submitted."
echo ""
echo "Monitor with:"
echo "  condor_q \$(whoami)"
echo "  condor_q -nobatch \$(whoami)"
echo ""
echo "When all jobs complete, run:"
echo "  bash collect_results.sh"
echo "═══════════════════════════════════════════════════════════"
