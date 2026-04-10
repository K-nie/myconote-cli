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
    plant)   PROTEIN_DB="/staging/braker/odb11/Viridiplantae.fa" ;;
    animal)  PROTEIN_DB="/staging/braker/odb11/Metazoa.fa" ;;
    insect)  PROTEIN_DB="/staging/braker/odb11/Arthropoda.fa" ;;
    protist) PROTEIN_DB="/staging/braker/odb11/Eukaryota.fa" ;;
    *)       PROTEIN_DB="/staging/braker/odb11/Eukaryota.fa" ;;
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
[[ -f "$PREDICTED_GFF" ]] || PREDICTED_GFF=$(find "$OUT_DIR" -name "*.gff3" | head -1)

python3 "$BENCHMARK_DIR/scripts/compare_annotations.py" \
    "$PREDICTED_GFF" \
    "$REFERENCE_GFF" \
    --label "braker_${GENOME_ID}_rep${REP}" \
    --output "$OUT_DIR/metrics.json"

python3 - << PYEOF > "$OUT_DIR/performance.json"
import json, re
data = {}
try:
    with open('$OUT_DIR/time_braker.log') as f:
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

print(json.dumps({
    'tool': 'braker3',
    'genome': '$GENOME_ID',
    'replicate': $REP,
    'total_wall_seconds': round(data.get('wall_seconds', 0), 1),
    'peak_rss_mb': round(data.get('max_rss_kb', 0) / 1024, 1),
}, indent=2))
PYEOF

echo "Done. Metrics: $OUT_DIR/metrics.json"
