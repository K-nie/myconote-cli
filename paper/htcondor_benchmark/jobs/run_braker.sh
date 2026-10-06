#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_braker.sh
# Wrapper for BRAKER3 benchmark runs.
# Note: BRAKER works best with RNA-seq evidence; for this benchmark we use
# BRAKER3 in protein-only mode (BRAKER2-style) using OrthoDB.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GENOME_ID="$1"
REP="$2"
BENCHMARK_DIR="$3"
DATA_DIR="$4"
RESULTS_DIR="$5"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate braker

GENOME_DIR="$DATA_DIR/$GENOME_ID"
GENOME_FA="$GENOME_DIR/genome.fa"
REFERENCE_GFF="$GENOME_DIR/reference.gff3"

OUT_DIR="$RESULTS_DIR/braker/$GENOME_ID/rep$REP"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

LOG="$OUT_DIR/run.log"

# Use OrthoDB protein database appropriate for the kingdom
CONFIG="$BENCHMARK_DIR/configs/genomes.tsv"
KINGDOM=$(awk -F'\t' -v id="$GENOME_ID" '$1==id {print $3}' "$CONFIG")

case "$KINGDOM" in
    fungi)   PROTEIN_DB="/staging/braker/odb11/Fungi.fa" ;;
    # Non-fungal kingdoms intentionally omitted — the benchmark panel
    # (configs/genomes.tsv) is fungi-only. If a future release validates
    # the --kingdom plant|animal|insect|protist paths, restore these
    # cases: Viridiplantae.fa, Metazoa.fa, Arthropoda.fa, Eukaryota.fa.
    *)       echo "ERROR: unsupported kingdom '$KINGDOM' for this benchmark panel" >&2
             exit 1 ;;
esac

/usr/bin/time -v -o "$OUT_DIR/time_braker.log" \
    braker.pl \
    --genome="$GENOME_FA" \
    --prot_seq="$PROTEIN_DB" \
    --species=braker_${GENOME_ID}_rep${REP} \
    --workingdir="$OUT_DIR" \
    --threads=16 \
    --useexisting \
    >> "$LOG" 2>&1

PREDICTED_GFF="$OUT_DIR/braker.gff3"
[[ -f "$PREDICTED_GFF" ]] || PREDICTED_GFF=$(find "$OUT_DIR" -name "*.gff3" -print -quit 2>/dev/null || true)

python3 "$BENCHMARK_DIR/scripts/compare_annotations.py" \
    "$PREDICTED_GFF" \
    "$REFERENCE_GFF" \
    --label "braker_${GENOME_ID}_rep${REP}" \
    --output "$OUT_DIR/metrics.json"

python3 - << PYEOF > "$OUT_DIR/performance.json"
import json

def parse_time_log(path):
    """Parse /usr/bin/time -v output. Byte-identical Elapsed/CPU parsing to
    run_myconote.sh so every tool's wall and CPU time are measured the same
    way. The old regex required fractional seconds and so dropped any stage
    >=1h (h:mm:ss, integer seconds), under-reporting wall. Split the trailing
    token on ':' and sum 1/2/3 fields instead.
    """
    data = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if 'Elapsed' in line and 'wall clock' in line:
                    tok = line.split()[-1]
                    try:
                        nums = [float(p) for p in tok.split(':')]
                    except ValueError:
                        nums = []
                    if len(nums) == 3:
                        data['wall_seconds'] = round(nums[0]*3600 + nums[1]*60 + nums[2], 2)
                    elif len(nums) == 2:
                        data['wall_seconds'] = round(nums[0]*60 + nums[1], 2)
                    elif len(nums) == 1:
                        data['wall_seconds'] = round(nums[0], 2)
                elif 'Maximum resident set size' in line:
                    data['max_rss_kb'] = int(line.split(':')[1].strip())
                elif 'User time' in line:
                    data['user_seconds'] = float(line.split(':')[1].strip())
                elif 'System time' in line:
                    data['system_seconds'] = float(line.split(':')[1].strip())
    except FileNotFoundError:
        pass
    return data

s = parse_time_log('$OUT_DIR/time_braker.log')
print(json.dumps({
    'tool': 'braker3',
    'genome': '$GENOME_ID',
    'replicate': $REP,
    'total_wall_seconds': round(s.get('wall_seconds', 0), 1),
    'total_cpu_seconds': round(s.get('user_seconds', 0) + s.get('system_seconds', 0), 1),
    'total_user_seconds': round(s.get('user_seconds', 0), 1),
    'total_system_seconds': round(s.get('system_seconds', 0), 1),
    'peak_rss_mb': round(s.get('max_rss_kb', 0) / 1024, 1),
}, indent=2))
PYEOF

echo "Done. Metrics: $OUT_DIR/metrics.json"
