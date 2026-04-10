#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# submit_all.sh
# Submit all 96 benchmark jobs to HTCondor.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BENCHMARK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${DATA_DIR:-$BENCHMARK_DIR/data}"
RESULTS_DIR="${RESULTS_DIR:-$BENCHMARK_DIR/results}"

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
echo ""

# Verify HTCondor is available
if ! command -v condor_submit &>/dev/null; then
    echo "ERROR: condor_submit not found. Are you on a submit node?"
    exit 1
fi

# ── Submit MycoNote-CLI jobs first (smallest) ─────────────────────────────
echo "── Submitting MycoNote-CLI jobs (24) ──"
condor_submit \
    -append "BENCHMARK_DIR = $BENCHMARK_DIR" \
    -append "DATA_DIR = $DATA_DIR" \
    -append "RESULTS_DIR = $RESULTS_DIR" \
    "$BENCHMARK_DIR/jobs/myconote.sub"
echo ""

# ── Submit funannotate jobs ───────────────────────────────────────────────
echo "── Submitting funannotate jobs (24) ──"
condor_submit \
    -append "BENCHMARK_DIR = $BENCHMARK_DIR" \
    -append "DATA_DIR = $DATA_DIR" \
    -append "RESULTS_DIR = $RESULTS_DIR" \
    "$BENCHMARK_DIR/jobs/funannotate.sub"
echo ""

# ── Submit BRAKER jobs ────────────────────────────────────────────────────
echo "── Submitting BRAKER jobs (24) ──"
condor_submit \
    -append "BENCHMARK_DIR = $BENCHMARK_DIR" \
    -append "DATA_DIR = $DATA_DIR" \
    -append "RESULTS_DIR = $RESULTS_DIR" \
    "$BENCHMARK_DIR/jobs/braker.sub"
echo ""

# ── Submit MAKER jobs (longest, last) ─────────────────────────────────────
echo "── Submitting MAKER jobs (24) ──"
condor_submit \
    -append "BENCHMARK_DIR = $BENCHMARK_DIR" \
    -append "DATA_DIR = $DATA_DIR" \
    -append "RESULTS_DIR = $RESULTS_DIR" \
    "$BENCHMARK_DIR/jobs/maker.sub"
echo ""

echo "═══════════════════════════════════════════════════════════"
echo "All 96 jobs submitted."
echo ""
echo "Monitor with:"
echo "  condor_q \$(whoami)"
echo "  condor_q -nobatch \$(whoami)"
echo ""
echo "When all jobs complete, run:"
echo "  bash collect_results.sh"
echo "═══════════════════════════════════════════════════════════"
