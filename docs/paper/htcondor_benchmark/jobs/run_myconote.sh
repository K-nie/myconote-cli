#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_myconote.sh
# Wrapper to run a single MycoNote-CLI benchmark and compute metrics.
#
# Arguments: $1=genome_id $2=replicate $3=benchmark_dir $4=data_dir $5=results_dir
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GENOME_ID="$1"
REP="$2"
BENCHMARK_DIR="$3"
DATA_DIR="$4"
RESULTS_DIR="$5"

GENOME_DIR="$DATA_DIR/$GENOME_ID"
GENOME_FA="$GENOME_DIR/genome.fa"
REFERENCE_GFF="$GENOME_DIR/reference.gff3"

OUT_DIR="$RESULTS_DIR/myconote/$GENOME_ID/rep$REP"
mkdir -p "$OUT_DIR"

LOG="$OUT_DIR/run.log"
TIME_LOG="$OUT_DIR/time.log"
METRICS_JSON="$OUT_DIR/metrics.json"

# ── Read genome metadata from configs/genomes.tsv ──────────────────────────
CONFIG="$BENCHMARK_DIR/configs/genomes.tsv"
GENOME_LINE=$(awk -F'\t' -v id="$GENOME_ID" '$1==id {print; exit}' "$CONFIG")

if [[ -z "$GENOME_LINE" ]]; then
    echo "ERROR: Genome ID not found: $GENOME_ID"
    exit 1
fi

ORGANISM=$(echo "$GENOME_LINE" | cut -f2)
KINGDOM=$(echo "$GENOME_LINE" | cut -f3)
GENETIC_CODE=$(echo "$GENOME_LINE" | cut -f5)
BUSCO_LINEAGE=$(echo "$GENOME_LINE" | cut -f6)

LOCUS_PREFIX=$(echo "$GENOME_ID" | tr '[:lower:]' '[:upper:]')

echo "════════════════════════════════════════════════════════"
echo "MycoNote-CLI Benchmark"
echo "Genome:     $ORGANISM ($GENOME_ID)"
echo "Replicate:  $REP"
echo "Kingdom:    $KINGDOM"
echo "Code:       $GENETIC_CODE"
echo "Threads:    8"
echo "Output:     $OUT_DIR"
echo "════════════════════════════════════════════════════════"

# ── Step 1: Sort ───────────────────────────────────────────────────────────
echo "[$(date +%T)] Step 1: sort"
/usr/bin/time -v -o "$OUT_DIR/time_sort.log" \
    myconote-cli sort "$GENOME_FA" \
    --output "$OUT_DIR/sorted.fa" \
    --min-length 500 \
    >> "$LOG" 2>&1

# ── Step 2: Mask ───────────────────────────────────────────────────────────
# Use the default SelfAlign engine: pure-Rust tandem-repeat finder plus
# minimap2 self-alignment (gracefully skips the minimap2 step if the
# binary isn't in PATH). RepeatModeler2 and RepeatMasker are NOT required,
# which matters on HPC clusters without those tools pre-installed.
echo "[$(date +%T)] Step 2: mask"
/usr/bin/time -v -o "$OUT_DIR/time_mask.log" \
    myconote-cli mask "$OUT_DIR/sorted.fa" \
    --output "$OUT_DIR/masked.fa" \
    --threads 8 \
    >> "$LOG" 2>&1

# ── Step 3: Predict ────────────────────────────────────────────────────────
echo "[$(date +%T)] Step 3: predict"
/usr/bin/time -v -o "$OUT_DIR/time_predict.log" \
    myconote-cli predict "$OUT_DIR/masked.fa" \
    --kingdom "$KINGDOM" \
    --locus-prefix "$LOCUS_PREFIX" \
    --output "$OUT_DIR/predict_out" \
    --threads 8 \
    >> "$LOG" 2>&1

# ── Step 4: Annotate ───────────────────────────────────────────────────────
echo "[$(date +%T)] Step 4: annotate"
/usr/bin/time -v -o "$OUT_DIR/time_annotate.log" \
    myconote-cli annotate "$OUT_DIR/predict_out/consensus.gff3" \
    --fasta "$OUT_DIR/masked.fa" \
    --output "$OUT_DIR/annotate_out" \
    --kingdom "$KINGDOM" \
    --locus-prefix "$LOCUS_PREFIX" \
    --genetic-code "$GENETIC_CODE" \
    --threads 8 \
    >> "$LOG" 2>&1

# ── Step 5: Validate for NCBI ──────────────────────────────────────────────
echo "[$(date +%T)] Step 5: submit validation"
myconote-cli submit "$OUT_DIR/annotate_out/annotated.gff3" \
    --fasta "$OUT_DIR/masked.fa" \
    --organism "$ORGANISM" \
    --validate-only \
    > "$OUT_DIR/submission_validation.txt" 2>&1 || true

# ── Step 6: Compute accuracy metrics vs reference ──────────────────────────
echo "[$(date +%T)] Step 6: comparing to reference annotation"
python3 "$BENCHMARK_DIR/scripts/compare_annotations.py" \
    "$OUT_DIR/annotate_out/annotated.gff3" \
    "$REFERENCE_GFF" \
    --label "myconote_${GENOME_ID}_rep${REP}" \
    --output "$METRICS_JSON" \
    2>> "$LOG"

# ── Step 7: Aggregate timing data ──────────────────────────────────────────
echo "[$(date +%T)] Step 7: aggregating performance data"
python3 - << PYEOF > "$OUT_DIR/performance.json"
import json
import re

def parse_time_log(path):
    """Parse /usr/bin/time -v output."""
    data = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if 'Elapsed' in line and 'wall clock' in line:
                    m = re.search(r'(\d+):?(\d+):(\d+\.\d+)', line)
                    if m:
                        h, mm, ss = m.groups()
                        if h:
                            data['wall_seconds'] = int(h)*3600 + int(mm)*60 + float(ss)
                        else:
                            data['wall_seconds'] = int(mm)*60 + float(ss)
                elif 'Maximum resident set size' in line:
                    data['max_rss_kb'] = int(line.split(':')[1].strip())
                elif 'User time' in line:
                    data['user_seconds'] = float(line.split(':')[1].strip())
                elif 'System time' in line:
                    data['system_seconds'] = float(line.split(':')[1].strip())
    except FileNotFoundError:
        pass
    return data

result = {
    'tool': 'myconote-cli',
    'genome': '$GENOME_ID',
    'replicate': $REP,
    'stages': {
        'sort':     parse_time_log('$OUT_DIR/time_sort.log'),
        'mask':     parse_time_log('$OUT_DIR/time_mask.log'),
        'predict':  parse_time_log('$OUT_DIR/time_predict.log'),
        'annotate': parse_time_log('$OUT_DIR/time_annotate.log'),
    }
}

# Compute totals
total_wall = sum(s.get('wall_seconds', 0) for s in result['stages'].values())
max_rss = max((s.get('max_rss_kb', 0) for s in result['stages'].values()), default=0)
result['total_wall_seconds'] = round(total_wall, 1)
result['peak_rss_mb'] = round(max_rss / 1024, 1)

print(json.dumps(result, indent=2))
PYEOF

echo "[$(date +%T)] Done. Output: $OUT_DIR"
echo "Metrics: $METRICS_JSON"
