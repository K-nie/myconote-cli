#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_braker.sh
# Wrapper for BRAKER3 benchmark runs.
# Note: BRAKER works best with RNA-seq evidence; for this benchmark we use
# BRAKER3 in protein-only mode (BRAKER2-style) using OrthoDB.
# BRAKER trains its own Augustus/GeneMark models de novo per genome, so there
# is no --augustus_species to set — this is the fair, as-designed workflow.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GENOME_ID="$1"
REP="$2"
BENCHMARK_DIR="$3"
DATA_DIR="$4"
RESULTS_DIR="$5"

# ── Conda env ────────────────────────────────────────────────────────────
# BRAKER3 is installed as the `myconote_braker3` env on this pool (not `braker`).
# Source the base conda profile from wherever it lives, activate, and also
# prepend the env bin so braker.pl resolves even if activation no-ops.
BRAKER_ENV="${BRAKER_ENV:-myconote_braker3}"
for _p in "/opt/bifxapps/miniconda3/etc/profile.d/conda.sh" \
          "$HOME/miniconda3/etc/profile.d/conda.sh" \
          "$HOME/anaconda3/etc/profile.d/conda.sh"; do
    [[ -f "$_p" ]] && source "$_p" && break
done
conda activate "$BRAKER_ENV" 2>/dev/null \
    || conda activate "$HOME/.conda/envs/$BRAKER_ENV" 2>/dev/null \
    || echo "WARN: 'conda activate $BRAKER_ENV' failed; relying on bin PATH"
if [[ -d "$HOME/.conda/envs/$BRAKER_ENV/bin" ]]; then
    export PATH="$HOME/.conda/envs/$BRAKER_ENV/bin:$PATH"
fi
# The env's augustus (3.5.0) is reported "not executable on this machine" by
# braker.pl on some execute nodes — a missing shared lib (it runs fine on nodes
# that have it), not a bad binary. Prepend the env's lib so augustus resolves
# its libs everywhere. Without this, all 18 braker jobs failed at braker.pl:2551.
if [[ -d "$HOME/.conda/envs/$BRAKER_ENV/lib" ]]; then
    export LD_LIBRARY_PATH="$HOME/.conda/envs/$BRAKER_ENV/lib:${LD_LIBRARY_PATH:-}"
fi

# ── GeneMark ─────────────────────────────────────────────────────────────
# BRAKER's protein pipeline (GeneMark-EP+) requires the GeneMark-ES suite,
# which is license-gated and not shipped via conda. It is provisioned under
# ~/.myconote/tools with the academic key at ~/.gm_key. Export GENEMARK_PATH
# and prepend its bin so braker.pl resolves gmes_petap.pl.
export GENEMARK_PATH="${GENEMARK_PATH:-$HOME/.myconote/tools/gmes_linux_64_4}"
if [[ -d "$GENEMARK_PATH" ]]; then
    export PATH="$GENEMARK_PATH:$PATH"
fi
if [[ ! -s "$HOME/.gm_key" ]]; then
    echo "ERROR: GeneMark key not found at ~/.gm_key" >&2
    exit 1
fi
# ProtHint feeds protein splice hints into GeneMark-EP+; it ships inside the
# GeneMark-ES suite. braker.pl needs PROTHINT_PATH and prothint.py on PATH.
export PROTHINT_PATH="${PROTHINT_PATH:-$GENEMARK_PATH/ProtHint/bin}"
if [[ -d "$PROTHINT_PATH" ]]; then
    export PATH="$PROTHINT_PATH:$PATH"
fi

GENOME_DIR="$DATA_DIR/$GENOME_ID"
GENOME_FA="$GENOME_DIR/genome.fa"
REFERENCE_GFF="$GENOME_DIR/reference.gff3"

OUT_DIR="$RESULTS_DIR/braker/$GENOME_ID/rep$REP"
# braker.pl aborts if its --workingdir already holds its subdirs (e.g. a stale
# GeneMark-ES/ from a prior run): "Failed to create directory .../GeneMark-ES!".
# Wipe the dir first so re-runs are idempotent; first-time runs are unaffected.
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

LOG="$OUT_DIR/run.log"

# Per-job writable Augustus config. BRAKER trains a new Augustus species into
# AUGUSTUS_CONFIG_PATH/species and aborts if that species dir already exists
# (e.g. a prior run of this same genome×rep). Pointing every job at the shared
# env config also races when the 18-job panel runs in parallel. Give each job
# its own fresh copy under OUT_DIR so training writes are isolated and reruns
# are clean. --fresh ensures a stale copy from an aborted run is replaced.
export AUGUSTUS_CONFIG_PATH="$OUT_DIR/augustus_config"
rm -rf "$AUGUSTUS_CONFIG_PATH"
cp -r "$HOME/.conda/envs/$BRAKER_ENV/config" "$AUGUSTUS_CONFIG_PATH"
rm -rf "$AUGUSTUS_CONFIG_PATH/species/braker_${GENOME_ID}_rep${REP}"

# OrthoDB protein partition (odb11). /staging is execute-node scratch and is not
# present on every node, so the partition is provisioned on cephfs home instead.
PROTEIN_DB="${BRAKER_PROTEIN_DB:-$HOME/.myconote/dbs/orthodb/Fungi.fa}"
CONFIG="$BENCHMARK_DIR/configs/genomes.tsv"
KINGDOM=$(awk -F'\t' -v id="$GENOME_ID" '$1==id {print $3}' "$CONFIG")
if [[ "$KINGDOM" != "fungi" ]]; then
    echo "ERROR: unsupported kingdom '$KINGDOM' for this benchmark panel" >&2
    exit 1
fi
if [[ ! -s "$PROTEIN_DB" ]]; then
    echo "ERROR: OrthoDB protein partition not found: $PROTEIN_DB" >&2
    exit 1
fi

/usr/bin/time -v -o "$OUT_DIR/time_braker.log" \
    braker.pl \
    --genome="$GENOME_FA" \
    --prot_seq="$PROTEIN_DB" \
    --species=braker_${GENOME_ID}_rep${REP} \
    --workingdir="$OUT_DIR" \
    --threads=16 \
    >> "$LOG" 2>&1

PREDICTED_GFF="$OUT_DIR/braker.gff3"
[[ -f "$PREDICTED_GFF" ]] || PREDICTED_GFF=$(find "$OUT_DIR" -name "*.gff3" -print -quit 2>/dev/null || true)

# BRAKER runs on the raw genome, so predicted seqids keep the original
# accessions (NC_...) and match reference.gff3 directly — no rename lift-over.
python3 "$BENCHMARK_DIR/scripts/compare_annotations.py" \
    "$PREDICTED_GFF" \
    "$REFERENCE_GFF" \
    --label "braker_${GENOME_ID}_rep${REP}" \
    --output "$OUT_DIR/metrics.json"

python3 - << PYEOF > "$OUT_DIR/performance.json"
import json

def parse_time(path):
    data = {}
    try:
        for line in open(path):
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

t = parse_time('$OUT_DIR/time_braker.log')
print(json.dumps({
    'tool': 'braker3',
    'genome': '$GENOME_ID',
    'replicate': $REP,
    'total_wall_seconds': round(t.get('wall_seconds', 0), 1),
    'total_cpu_seconds': round(t.get('user_seconds', 0) + t.get('system_seconds', 0), 1),
    'total_user_seconds': round(t.get('user_seconds', 0), 1),
    'total_system_seconds': round(t.get('system_seconds', 0), 1),
    'peak_rss_mb': round(t.get('max_rss_kb', 0) / 1024, 1),
}, indent=2))
PYEOF

echo "Done. Metrics: $OUT_DIR/metrics.json"
