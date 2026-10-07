#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_funannotate.sh
# Wrapper for funannotate benchmark runs.
# Strategy (native + clade Augustus): funannotate trains its own Augustus model
# de novo per genome, seeded with the clade-appropriate BUSCO species via
# --busco_seed_species. This is funannotate's as-designed workflow — we do NOT
# hand it a pre-trained S288C model, which would be unfair across clades.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GENOME_ID="$1"
REP="$2"
BENCHMARK_DIR="$3"
DATA_DIR="$4"
RESULTS_DIR="$5"

# ── Conda env ────────────────────────────────────────────────────────────
# Source the base conda profile from wherever it lives, activate, and also
# prepend the env bin so funannotate resolves even if activation no-ops on an
# execute node where the base install isn't mounted.
FUN_ENV="${FUNANNOTATE_ENV:-funannotate}"
for _p in "/opt/bifxapps/miniconda3/etc/profile.d/conda.sh" \
          "$HOME/miniconda3/etc/profile.d/conda.sh" \
          "$HOME/anaconda3/etc/profile.d/conda.sh"; do
    [[ -f "$_p" ]] && source "$_p" && break
done
conda activate "$FUN_ENV" 2>/dev/null \
    || conda activate "$HOME/.conda/envs/$FUN_ENV" 2>/dev/null \
    || echo "WARN: 'conda activate $FUN_ENV' failed; relying on bin PATH"
if [[ -d "$HOME/.conda/envs/$FUN_ENV/bin" ]]; then
    export PATH="$HOME/.conda/envs/$FUN_ENV/bin:$PATH"
fi

# Augustus writes trained species into its config tree, so AUGUSTUS_CONFIG_PATH
# must point at a writable copy; the funannotate env ships one under config/.
export AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH:-$HOME/.conda/envs/$FUN_ENV/config}"

# funannotate annotate needs its reference DBs (Pfam, dbCAN, MEROPS, UniProt,
# InterPro, BUSCO), provisioned once via `funannotate setup -d`.
export FUNANNOTATE_DB="${FUNANNOTATE_DB:-$HOME/.myconote/dbs/funannotate}"
if [[ ! -d "$FUNANNOTATE_DB" ]]; then
    echo "ERROR: FUNANNOTATE_DB not found: $FUNANNOTATE_DB (run funannotate setup -d)" >&2
    exit 1
fi

GENOME_DIR="$DATA_DIR/$GENOME_ID"
GENOME_FA="$GENOME_DIR/genome.fa"
REFERENCE_GFF="$GENOME_DIR/reference.gff3"

OUT_DIR="$RESULTS_DIR/funannotate/$GENOME_ID/rep$REP"
# funannotate predict aborts if predict_out/ already exists; wipe OUT_DIR so
# re-runs are idempotent and don't collide with a prior run's output.
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

LOG="$OUT_DIR/run.log"
METRICS_JSON="$OUT_DIR/metrics.json"

CONFIG="$BENCHMARK_DIR/configs/genomes.tsv"
GENOME_LINE=$(awk -F'\t' -v id="$GENOME_ID" '$1==id {print; exit}' "$CONFIG")
ORGANISM=$(echo "$GENOME_LINE" | cut -f2 | tr '_' ' ')
KINGDOM=$(echo "$GENOME_LINE" | cut -f3)

LOCUS_PREFIX=$(echo "$GENOME_ID" | tr '[:lower:]' '[:upper:]')

# Per-genome clade Augustus species used to seed de-novo BUSCO training.
case "$GENOME_ID" in
    sce) AUGUSTUS_SPECIES="saccharomyces_cerevisiae_S288C" ;;
    cal) AUGUSTUS_SPECIES="candida_albicans" ;;
    ylp) AUGUSTUS_SPECIES="yarrowia_lipolytica" ;;
    ani) AUGUSTUS_SPECIES="aspergillus_nidulans" ;;
    ncr) AUGUSTUS_SPECIES="neurospora_crassa" ;;
    cne) AUGUSTUS_SPECIES="cryptococcus_neoformans_neoformans_JEC21" ;;
    *)   echo "ERROR: no Augustus seed species mapped for '$GENOME_ID'" >&2; exit 1 ;;
esac

cd "$OUT_DIR"

# ── Step 1: clean ──────────────────────────────────────────────────────────
/usr/bin/time -v -o "$OUT_DIR/time_clean.log" \
    funannotate clean -i "$GENOME_FA" -o cleaned.fa --minlen 500 \
    >> "$LOG" 2>&1

# ── Step 2: sort ───────────────────────────────────────────────────────────
# --minlen 0 is passed explicitly: funannotate 1.8.17's sort.py does not apply
# its own argparse default, so minlen stays None and `if minlen > 0` throws
# TypeError. clean (Step 1) already dropped <500 bp contigs, so 0 keeps the rest.
/usr/bin/time -v -o "$OUT_DIR/time_sort.log" \
    funannotate sort -i cleaned.fa -o sorted.fa -b scaffold --minlen 0 \
    >> "$LOG" 2>&1

# ── Rename lift-over ─────────────────────────────────────────────────────
# funannotate sort renames contigs to scaffold_N (longest first), so the
# predicted GFF3 uses scaffold_N while reference.gff3 keeps NCBI accessions.
# Sort does not alter sequence content, so md5-match each sorted contig back
# to its pre-sort id to build the original_id/new_id/length table that
# compare_annotations.py --rename-table expects (its first line is a header).
python3 - "$OUT_DIR/cleaned.fa" "$OUT_DIR/sorted.fa" "$OUT_DIR/rename_table.tsv" << 'PYEOF'
import sys, hashlib

def read_fa(path):
    out, name, seq = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        out[name] = ''.join(seq)
    return out

orig = read_fa(sys.argv[1])
new = read_fa(sys.argv[2])
md5_orig = {hashlib.md5(s.upper().encode()).hexdigest(): (n, len(s)) for n, s in orig.items()}
with open(sys.argv[3], 'w') as fh:
    fh.write("original_id\tnew_id\tlength\n")
    for new_id, s in new.items():
        h = hashlib.md5(s.upper().encode()).hexdigest()
        if h in md5_orig:
            orig_id, ln = md5_orig[h]
            fh.write(f"{orig_id}\t{new_id}\t{ln}\n")
PYEOF

# ── Step 3: mask ───────────────────────────────────────────────────────────
/usr/bin/time -v -o "$OUT_DIR/time_mask.log" \
    funannotate mask -i sorted.fa -o masked.fa --cpus 8 \
    >> "$LOG" 2>&1

# ── Step 4: predict (de-novo Augustus, clade-seeded) ─────────────────────────
/usr/bin/time -v -o "$OUT_DIR/time_predict.log" \
    funannotate predict -i masked.fa -o predict_out \
    --species "$ORGANISM" --strain rep$REP \
    --name "$LOCUS_PREFIX" --cpus 8 --busco_seed_species "$AUGUSTUS_SPECIES" \
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
    --rename-table "$OUT_DIR/rename_table.tsv" \
    --output "$METRICS_JSON"

# ── Step 7: Performance ────────────────────────────────────────────────────
python3 - << PYEOF > "$OUT_DIR/performance.json"
import json

def parse_time(path):
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
total_user = sum(s.get('user_seconds', 0) for s in result['stages'].values())
total_system = sum(s.get('system_seconds', 0) for s in result['stages'].values())
max_rss = max((s.get('max_rss_kb', 0) for s in result['stages'].values()), default=0)
result['total_wall_seconds'] = round(total_wall, 1)
result['total_cpu_seconds'] = round(total_user + total_system, 1)
result['total_user_seconds'] = round(total_user, 1)
result['total_system_seconds'] = round(total_system, 1)
result['peak_rss_mb'] = round(max_rss / 1024, 1)
print(json.dumps(result, indent=2))
PYEOF

echo "Done. Metrics: $METRICS_JSON"
