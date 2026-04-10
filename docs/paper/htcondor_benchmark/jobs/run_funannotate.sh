#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_funannotate.sh
# Wrapper for funannotate benchmark runs.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GENOME_ID="$1"
REP="$2"
BENCHMARK_DIR="$3"
DATA_DIR="$4"
RESULTS_DIR="$5"

# Activate funannotate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate funannotate

GENOME_DIR="$DATA_DIR/$GENOME_ID"
GENOME_FA="$GENOME_DIR/genome.fa"
REFERENCE_GFF="$GENOME_DIR/reference.gff3"

OUT_DIR="$RESULTS_DIR/funannotate/$GENOME_ID/rep$REP"
mkdir -p "$OUT_DIR"

LOG="$OUT_DIR/run.log"
METRICS_JSON="$OUT_DIR/metrics.json"

CONFIG="$BENCHMARK_DIR/configs/genomes.tsv"
GENOME_LINE=$(awk -F'\t' -v id="$GENOME_ID" '$1==id {print; exit}' "$CONFIG")
ORGANISM=$(echo "$GENOME_LINE" | cut -f2 | tr '_' ' ')
KINGDOM=$(echo "$GENOME_LINE" | cut -f3)

LOCUS_PREFIX=$(echo "$GENOME_ID" | tr '[:lower:]' '[:upper:]')

cd "$OUT_DIR"

# ── Step 1: clean ──────────────────────────────────────────────────────────
/usr/bin/time -v -o "$OUT_DIR/time_clean.log" \
    funannotate clean -i "$GENOME_FA" -o cleaned.fa --minlen 500 \
    >> "$LOG" 2>&1

# ── Step 2: sort ───────────────────────────────────────────────────────────
/usr/bin/time -v -o "$OUT_DIR/time_sort.log" \
    funannotate sort -i cleaned.fa -o sorted.fa -b scaffold \
    >> "$LOG" 2>&1

# ── Step 3: mask ───────────────────────────────────────────────────────────
/usr/bin/time -v -o "$OUT_DIR/time_mask.log" \
    funannotate mask -i sorted.fa -o masked.fa --cpus 8 \
    >> "$LOG" 2>&1

# ── Step 4: predict ────────────────────────────────────────────────────────
/usr/bin/time -v -o "$OUT_DIR/time_predict.log" \
    funannotate predict -i masked.fa -o predict_out \
    --species "$ORGANISM" --strain rep$REP \
    --name "$LOCUS_PREFIX" --cpus 8 --augustus_species saccharomyces_cerevisiae_S288C \
    >> "$LOG" 2>&1

# ── Step 5: annotate ───────────────────────────────────────────────────────
/usr/bin/time -v -o "$OUT_DIR/time_annotate.log" \
    funannotate annotate -i predict_out --cpus 8 \
    >> "$LOG" 2>&1

# Find the final annotated GFF3. Use `find ... -print -quit` rather than
# `find | head -1` so the pipeline can't trip set -o pipefail via SIGPIPE.
PREDICTED_GFF=$(find predict_out -path "*annotate*" -name "*.gff3" -print -quit 2>/dev/null || true)
if [[ -z "$PREDICTED_GFF" ]]; then
    PREDICTED_GFF=$(find predict_out -name "*.gff3" -print -quit 2>/dev/null || true)
fi

# ── Step 6: Compute metrics ────────────────────────────────────────────────
python3 "$BENCHMARK_DIR/scripts/compare_annotations.py" \
    "$PREDICTED_GFF" \
    "$REFERENCE_GFF" \
    --label "funannotate_${GENOME_ID}_rep${REP}" \
    --output "$METRICS_JSON"

# ── Step 7: Performance ────────────────────────────────────────────────────
python3 - << PYEOF > "$OUT_DIR/performance.json"
import json, re

def parse_time(path):
    data = {}
    try:
        with open(path) as f:
            for line in f:
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
    except FileNotFoundError:
        pass
    return data

result = {
    'tool': 'funannotate',
    'genome': '$GENOME_ID',
    'replicate': $REP,
    'stages': {
        'clean':    parse_time('$OUT_DIR/time_clean.log'),
        'sort':     parse_time('$OUT_DIR/time_sort.log'),
        'mask':     parse_time('$OUT_DIR/time_mask.log'),
        'predict':  parse_time('$OUT_DIR/time_predict.log'),
        'annotate': parse_time('$OUT_DIR/time_annotate.log'),
    }
}

total_wall = sum(s.get('wall_seconds', 0) for s in result['stages'].values())
max_rss = max((s.get('max_rss_kb', 0) for s in result['stages'].values()), default=0)
result['total_wall_seconds'] = round(total_wall, 1)
result['peak_rss_mb'] = round(max_rss / 1024, 1)
print(json.dumps(result, indent=2))
PYEOF

echo "Done. Metrics: $METRICS_JSON"
