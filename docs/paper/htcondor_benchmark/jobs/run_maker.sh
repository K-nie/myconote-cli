#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_maker.sh
# Wrapper for MAKER benchmark runs.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GENOME_ID="$1"
REP="$2"
BENCHMARK_DIR="$3"
DATA_DIR="$4"
RESULTS_DIR="$5"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate maker

GENOME_DIR="$DATA_DIR/$GENOME_ID"
GENOME_FA="$GENOME_DIR/genome.fa"
REFERENCE_GFF="$GENOME_DIR/reference.gff3"

OUT_DIR="$RESULTS_DIR/maker/$GENOME_ID/rep$REP"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

LOG="$OUT_DIR/run.log"

# Generate maker control files
maker -CTL >> "$LOG" 2>&1

# Edit maker_opts.ctl to point at our genome
sed -i "s|^genome=.*|genome=$GENOME_FA|" maker_opts.ctl
sed -i "s|^model_org=.*|model_org=fungi|" maker_opts.ctl
sed -i "s|^augustus_species=.*|augustus_species=saccharomyces_cerevisiae_S288C|" maker_opts.ctl
sed -i "s|^cpus=.*|cpus=8|" maker_opts.ctl

/usr/bin/time -v -o "$OUT_DIR/time_maker.log" \
    maker maker_opts.ctl maker_bopts.ctl maker_exe.ctl \
    >> "$LOG" 2>&1

# Find MAKER's output GFF3. Use `find -print -quit` so the pipeline
# can't trip set -o pipefail via SIGPIPE on early head exit.
PREDICTED_GFF=$(find . -name "*.all.gff" -print -quit 2>/dev/null || true)
if [[ -z "$PREDICTED_GFF" ]]; then
    PREDICTED_GFF=$(find . -name "*.gff" -print -quit 2>/dev/null || true)
fi

# Compute metrics
python3 "$BENCHMARK_DIR/scripts/compare_annotations.py" \
    "$PREDICTED_GFF" \
    "$REFERENCE_GFF" \
    --label "maker_${GENOME_ID}_rep${REP}" \
    --output "$OUT_DIR/metrics.json"

# Performance summary
python3 - << PYEOF > "$OUT_DIR/performance.json"
import json, re
data = {}
try:
    with open('$OUT_DIR/time_maker.log') as f:
        for line in f:
            if 'Elapsed' in line and 'wall clock' in line:
                m = re.search(r'(\d+):?(\d+):(\d+\.\d+)', line)
                if m:
                    h, mm, ss = m.groups()
                    data['wall_seconds'] = (int(h)*3600 if h else 0) + int(mm)*60 + float(ss)
            elif 'Maximum resident set size' in line:
                data['max_rss_kb'] = int(line.split(':')[1].strip())
except FileNotFoundError:
    pass

result = {
    'tool': 'maker',
    'genome': '$GENOME_ID',
    'replicate': $REP,
    'total_wall_seconds': round(data.get('wall_seconds', 0), 1),
    'peak_rss_mb': round(data.get('max_rss_kb', 0) / 1024, 1),
}
print(json.dumps(result, indent=2))
PYEOF

echo "Done. Metrics: $OUT_DIR/metrics.json"
