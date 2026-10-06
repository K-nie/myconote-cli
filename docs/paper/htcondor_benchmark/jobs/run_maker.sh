#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_maker.sh
# Wrapper for MAKER benchmark runs.
# Strategy (native + clade Augustus): MAKER is given the clade-appropriate
# pre-trained Augustus species per genome via augustus_species in maker_opts.ctl.
# MAKER runs on the raw genome, so predicted seqids keep the original NCBI
# accessions and match reference.gff3 directly — no rename lift-over.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GENOME_ID="$1"
REP="$2"
BENCHMARK_DIR="$3"
DATA_DIR="$4"
RESULTS_DIR="$5"

# ── Conda env ────────────────────────────────────────────────────────────
MAKER_ENV="${MAKER_ENV:-maker}"
for _p in "/opt/bifxapps/miniconda3/etc/profile.d/conda.sh" \
          "$HOME/miniconda3/etc/profile.d/conda.sh" \
          "$HOME/anaconda3/etc/profile.d/conda.sh"; do
    [[ -f "$_p" ]] && source "$_p" && break
done
conda activate "$MAKER_ENV" 2>/dev/null \
    || conda activate "$HOME/.conda/envs/$MAKER_ENV" 2>/dev/null \
    || echo "WARN: 'conda activate $MAKER_ENV' failed; relying on bin PATH"
if [[ -d "$HOME/.conda/envs/$MAKER_ENV/bin" ]]; then
    export PATH="$HOME/.conda/envs/$MAKER_ENV/bin:$PATH"
fi

GENOME_DIR="$DATA_DIR/$GENOME_ID"
GENOME_FA="$GENOME_DIR/genome.fa"
REFERENCE_GFF="$GENOME_DIR/reference.gff3"

OUT_DIR="$RESULTS_DIR/maker/$GENOME_ID/rep$REP"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

LOG="$OUT_DIR/run.log"

# Per-genome clade Augustus species (pre-trained models shipped with Augustus).
case "$GENOME_ID" in
    sce) AUGUSTUS_SPECIES="saccharomyces_cerevisiae_S288C" ;;
    cal) AUGUSTUS_SPECIES="candida_albicans" ;;
    ylp) AUGUSTUS_SPECIES="yarrowia_lipolytica" ;;
    ani) AUGUSTUS_SPECIES="aspergillus_nidulans" ;;
    ncr) AUGUSTUS_SPECIES="neurospora_crassa" ;;
    cne) AUGUSTUS_SPECIES="cryptococcus_neoformans_neoformans_JEC21" ;;
    *)   echo "ERROR: no Augustus species mapped for '$GENOME_ID'" >&2; exit 1 ;;
esac

# Generate maker control files
maker -CTL >> "$LOG" 2>&1

# Edit maker_opts.ctl to point at our genome and the clade model.
# Repeat masking is disabled (model_org blank + RepeatMasker removed from the
# exe list): MAKER 3.01.04 cannot parse RepeatMasker 4.2.3 / FamDB-3.9 and
# aborts at startup with "Could not determine if DFam is installed" (the env
# ships only the Dfam root partition; the fungi partition is absent). The
# benchmark panel is repeat-poor ascomycete/basidiomycete genomes, so skipping
# repeat masking has minimal impact on de-novo gene prediction, and it keeps the
# MAKER arm reproducible without a curated Dfam/RepBase partition.
sed -i "s|^genome=.*|genome=$GENOME_FA|" maker_opts.ctl
sed -i "s|^model_org=.*|model_org=|" maker_opts.ctl
sed -i "s|^augustus_species=.*|augustus_species=$AUGUSTUS_SPECIES|" maker_opts.ctl
sed -i "s|^cpus=.*|cpus=8|" maker_opts.ctl
# keep_preds=1: no EST/protein evidence is supplied, so MAKER would otherwise
# drop all unsupported ab-initio models and emit an empty annotation. Keeping
# predictions makes MAKER a clade-Augustus ab-initio predictor — the direct
# analogue of the native+clade Augustus strategy used for the other arms.
sed -i "s|^keep_preds=.*|keep_preds=1|" maker_opts.ctl
# Remove RepeatMasker from the exe list so MAKER skips the broken FamDB probe.
sed -i "s|^RepeatMasker=.*|RepeatMasker=|" maker_exe.ctl

/usr/bin/time -v -o "$OUT_DIR/time_maker.log" \
    maker maker_opts.ctl maker_bopts.ctl maker_exe.ctl \
    >> "$LOG" 2>&1

# MAKER writes per-contig GFFs into the output datastore; the final genome-wide
# annotation only exists after gff3_merge walks the datastore index. Skipping
# this leaves no *.all.gff, so a bare `find *.gff` would match a single-contig
# evidence fragment (0 genes → every metric 0). -n drops the embedded FASTA so
# the merged GFF is pure annotation for the comparator.
DS_INDEX=$(find . -name "*_master_datastore_index.log" -print -quit 2>/dev/null || true)
PREDICTED_GFF="$OUT_DIR/genome.all.gff"
gff3_merge -n -d "$DS_INDEX" -o "$PREDICTED_GFF" >> "$LOG" 2>&1

# Fall back to any merged *.all.gff if the explicit output path is missing.
if [[ ! -s "$PREDICTED_GFF" ]]; then
    PREDICTED_GFF=$(find . -name "*.all.gff" -print -quit 2>/dev/null || true)
fi

# Compute metrics
python3 "$BENCHMARK_DIR/scripts/compare_annotations.py" \
    "$PREDICTED_GFF" \
    "$REFERENCE_GFF" \
    --label "maker_${GENOME_ID}_rep${REP}" \
    --output "$OUT_DIR/metrics.json"

# Performance summary
python3 - << PYEOF > "$OUT_DIR/performance.json"
import json

data = {}
try:
    with open('$OUT_DIR/time_maker.log') as f:
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

result = {
    'tool': 'maker',
    'genome': '$GENOME_ID',
    'replicate': $REP,
    'total_wall_seconds': round(data.get('wall_seconds', 0), 1),
    'total_cpu_seconds': round(data.get('user_seconds', 0) + data.get('system_seconds', 0), 1),
    'total_user_seconds': round(data.get('user_seconds', 0), 1),
    'total_system_seconds': round(data.get('system_seconds', 0), 1),
    'peak_rss_mb': round(data.get('max_rss_kb', 0) / 1024, 1),
}
print(json.dumps(result, indent=2))
PYEOF

echo "Done. Metrics: $OUT_DIR/metrics.json"
