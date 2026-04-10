#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# setup.sh
# One-time setup for the MycoNote-CLI benchmark on HTCondor.
#
# - Downloads all 8 reference genomes and curated annotations
# - Verifies that conda environments for funannotate, maker, braker exist
# - Creates the results directory structure
# - Reports any missing dependencies
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BENCHMARK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${DATA_DIR:-$BENCHMARK_DIR/data}"
RESULTS_DIR="${RESULTS_DIR:-$BENCHMARK_DIR/results}"

echo "═══════════════════════════════════════════════════════════"
echo "MycoNote-CLI Benchmark Setup"
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "Benchmark dir : $BENCHMARK_DIR"
echo "Data dir      : $DATA_DIR"
echo "Results dir   : $RESULTS_DIR"
echo ""

# ── Step 1: Download reference genomes ────────────────────────────────────
echo "── Step 1: Downloading reference genomes ──"
bash "$BENCHMARK_DIR/scripts/download_references.sh" "$DATA_DIR"
echo ""

# ── Step 2: Verify conda environments ─────────────────────────────────────
echo "── Step 2: Verifying conda environments ──"

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || \
    source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || \
    { echo "ERROR: Could not find conda installation"; exit 1; }

check_env() {
    local env_name="$1"
    if conda env list | grep -q "^$env_name "; then
        echo "  ✓ Environment '$env_name' exists"
    else
        echo "  ✗ Environment '$env_name' MISSING - install with:"
        case "$env_name" in
            funannotate) echo "      conda create -n funannotate -c bioconda funannotate" ;;
            maker)       echo "      conda create -n maker -c bioconda maker" ;;
            braker)      echo "      conda create -n braker -c bioconda braker3" ;;
        esac
    fi
}

check_env funannotate
check_env maker
check_env braker

# Check myconote-cli is in PATH
if command -v myconote-cli &>/dev/null; then
    echo "  ✓ myconote-cli found: $(myconote-cli --version 2>&1 | head -1 || echo 'unknown')"
else
    echo "  ✗ myconote-cli NOT in PATH - install with:"
    echo "      curl -fsSL https://raw.githubusercontent.com/K-nie/myconote-cli/main/quick-install.sh | bash"
fi
echo ""

# ── Step 3: Create results directory structure ────────────────────────────
echo "── Step 3: Creating results directory structure ──"
mkdir -p "$RESULTS_DIR"/{myconote,funannotate,maker,braker}
echo "  ✓ Created $RESULTS_DIR"
echo ""

# ── Step 4: Make all scripts executable ───────────────────────────────────
echo "── Step 4: Making scripts executable ──"
chmod +x "$BENCHMARK_DIR/jobs"/*.sh
chmod +x "$BENCHMARK_DIR/scripts"/*.sh
chmod +x "$BENCHMARK_DIR/scripts"/*.py
echo "  ✓ All scripts executable"
echo ""

# ── Step 5: Verify download integrity ─────────────────────────────────────
echo "── Step 5: Verifying downloaded data ──"
total_genomes=0
ok_genomes=0
while IFS=$'\t' read -r id rest; do
    if [[ "$id" == "id" ]]; then continue; fi
    total_genomes=$((total_genomes + 1))
    if [[ -f "$DATA_DIR/$id/genome.fa" ]] && [[ -f "$DATA_DIR/$id/reference.gff3" ]]; then
        ok_genomes=$((ok_genomes + 1))
    fi
done < "$BENCHMARK_DIR/configs/genomes.tsv"
echo "  $ok_genomes / $total_genomes genomes ready"
echo ""

# ── Final check ──────────────────────────────────────────────────────────
echo "═══════════════════════════════════════════════════════════"
if [[ "$ok_genomes" -eq "$total_genomes" ]]; then
    echo "✓ Setup complete. Ready to submit benchmark jobs."
    echo ""
    echo "Next steps:"
    echo "  1. Verify environment variables in submit files match your setup"
    echo "  2. Submit jobs with: bash submit_all.sh"
else
    echo "⚠ Setup incomplete. Please address missing items above."
fi
echo "═══════════════════════════════════════════════════════════"
