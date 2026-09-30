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

# ── Database locations ─────────────────────────────────────────────────────
# annotate DB-backed sources (--eggnog/--cazyme/--merops) plus Swiss-Prot and
# Pfam auto-resolve from --db-dir. Default matches the tool's own default and
# htcondor.conf; override MYCONOTE_DB_DIR (e.g. to shared storage) upstream.
MYCONOTE_DB_DIR="${MYCONOTE_DB_DIR:-$HOME/.myconote/dbs}"
# predict/train read Augustus species from this FIXED path (ignores --db-dir).
export AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH:-$HOME/.myconote/augustus_config}"

# annotate --eggnog resolves the emapper.py DB from --eggnog-db, NOT --db-dir
# (the tool never derives it from --db-dir). Point it at the copy provisioned
# by setup.sh under $MYCONOTE_DB_DIR/eggnog.
EGGNOG_DB_DIR="${EGGNOG_DB_DIR:-$MYCONOTE_DB_DIR/eggnog}"

# ── Runtime tool environment ───────────────────────────────────────────────
# predict/annotate shell out to external binaries that MycoNote locates by a
# bare PATH lookup (augustus, snap, diamond, hmmsearch, minimap2, run_dbcan,
# emapper.py). These live in conda envs on this cluster, so put them on PATH:
#   - base "myconote" env: snap, diamond, hmmsearch, minimap2, run_dbcan and
#     emapper.py (all working). Its bundled augustus links the wrong boost, so
#     we do NOT use it for augustus.
#   - "myconote_augustus" env: a working augustus 3.3.3, prepended so it wins.
# The metrics scripts (Steps 6-7) are stdlib-only, so the env's python3 is fine.
MYCONOTE_ENV="${MYCONOTE_ENV:-myconote}"
AUGUSTUS_ENV="${AUGUSTUS_ENV:-myconote_augustus}"
if [[ -z "${CONDA_PROFILE:-}" ]]; then
    for _p in "/opt/bifxapps/miniconda3/etc/profile.d/conda.sh" \
              "$HOME/miniconda3/etc/profile.d/conda.sh" \
              "$HOME/anaconda3/etc/profile.d/conda.sh"; do
        [[ -f "$_p" ]] && CONDA_PROFILE="$_p" && break
    done
fi
if [[ -n "${CONDA_PROFILE:-}" && -f "$CONDA_PROFILE" ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_PROFILE"
    conda activate "$MYCONOTE_ENV" 2>/dev/null \
        || conda activate "$HOME/.conda/envs/$MYCONOTE_ENV" 2>/dev/null \
        || echo "WARN: 'conda activate $MYCONOTE_ENV' failed; using bin PATH fallback"
fi
# Robust fallback: `conda activate` needs the base conda's profile.d/conda.sh,
# which lives under /opt/bifxapps and is NOT mounted on every execute node — so
# activation silently no-ops there, and then snap/mmseqs/hmmscan/run_dbcan/
# emapper.py vanish from PATH (only augustus+diamond survive via the env prepend
# below), which makes annotate skip Swiss-Prot/Pfam/eggnog/CAZyme and predict
# fall back to Augustus-only. The env's bin dir is on shared cephfs and its
# binaries are RPATH-linked to ../lib, so prepending it directly gives every
# tool regardless of whether activation worked. Verified: mmseqs/hmmscan/snap/
# run_dbcan/emapper.py all run correctly from bin alone.
if [[ -d "$HOME/.conda/envs/$MYCONOTE_ENV/bin" ]]; then
    export PATH="$HOME/.conda/envs/$MYCONOTE_ENV/bin:$PATH"
fi
# Prepend the working augustus (self-contained via its own RPATH-linked boost).
for _aug in "$HOME/.conda/envs/$AUGUSTUS_ENV/bin" \
            "/opt/bifxapps/miniconda3/envs/$AUGUSTUS_ENV/bin"; do
    if [[ -x "$_aug/augustus" ]]; then export PATH="$_aug:$PATH"; break; fi
done
# Ensure the freshly built myconote-cli wins over anything the env put on PATH.
if [[ -x "$HOME/.local/bin/myconote-cli" ]]; then
    export PATH="$HOME/.local/bin:$PATH"
fi

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

# ── Clade-appropriate Augustus species model ────────────────────────────────
# predict defaults --kingdom fungi to saccharomyces_cerevisiae_S288C for EVERY
# genome (src/predict/kingdom.rs). That budding-yeast model is near-intron-less,
# so it wrecks gene structure on intron-rich fungi (e.g. Cryptococcus scored
# BUSCO 14.9% / strict gene F1 ~0 under S288C). funannotate/BRAKER train
# species-specific models, so a fair comparison must give MycoNote the matching
# model too. Each name below is an installed Augustus species dir on this pool;
# strain-exact where available (sce=S288C, cne=JEC21).
case "$GENOME_ID" in
    sce) AUGUSTUS_SPECIES="saccharomyces_cerevisiae_S288C" ;;
    cal) AUGUSTUS_SPECIES="candida_albicans" ;;
    ylp) AUGUSTUS_SPECIES="yarrowia_lipolytica" ;;
    ani) AUGUSTUS_SPECIES="aspergillus_nidulans" ;;
    ncr) AUGUSTUS_SPECIES="neurospora_crassa" ;;
    cne) AUGUSTUS_SPECIES="cryptococcus_neoformans_neoformans_JEC21" ;;
    *)   echo "ERROR: no Augustus species mapping for genome '$GENOME_ID'"; exit 1 ;;
esac

echo "════════════════════════════════════════════════════════"
echo "MycoNote-CLI Benchmark"
echo "Genome:     $ORGANISM ($GENOME_ID)"
echo "Replicate:  $REP"
echo "Kingdom:    $KINGDOM"
echo "Species:    $AUGUSTUS_SPECIES"
echo "Code:       $GENETIC_CODE"
echo "Threads:    8"
echo "Output:     $OUT_DIR"
echo "════════════════════════════════════════════════════════"

# ── Step 1: Sort ───────────────────────────────────────────────────────────
echo "[$(date +%T)] Step 1: sort"
# --rename-table records the original_id→scaffold_N contig mapping. sort renames
# contigs to scaffold_N (longest first) for NCBI-clean output, so the predicted
# GFF3 ends up in the scaffold_N namespace while reference.gff3 keeps the original
# accessions (e.g. NC_001133.9). Step 6 uses this table to lift predicted seqids
# back before comparing — without it every metric is 0 (no seqid overlaps).
/usr/bin/time -v -o "$OUT_DIR/time_sort.log" \
    myconote-cli sort "$GENOME_FA" \
    --output "$OUT_DIR/sorted.fa" \
    --min-length 500 \
    --rename-table "$OUT_DIR/rename_table.tsv" \
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
    --species "$AUGUSTUS_SPECIES" \
    --locus-prefix "$LOCUS_PREFIX" \
    --output "$OUT_DIR/predict_out" \
    --threads 8 \
    >> "$LOG" 2>&1

# ── Step 4: Annotate ───────────────────────────────────────────────────────
# Full functional-annotation suite, matched to funannotate's default workflow:
# Swiss-Prot + Pfam (always on) plus EggNog COG/NOG, CAZyme (dbCAN), and MEROPS
# protease families. All resolve from --db-dir. antiSMASH and InterProScan are
# deliberately excluded — they require network access from execute nodes.
echo "[$(date +%T)] Step 4: annotate"
/usr/bin/time -v -o "$OUT_DIR/time_annotate.log" \
    myconote-cli annotate "$OUT_DIR/predict_out/consensus.gff3" \
    --fasta "$OUT_DIR/masked.fa" \
    --output "$OUT_DIR/annotate_out" \
    --kingdom "$KINGDOM" \
    --locus-prefix "$LOCUS_PREFIX" \
    --genetic-code "$GENETIC_CODE" \
    --db-dir "$MYCONOTE_DB_DIR" \
    --eggnog \
    --eggnog-db "$EGGNOG_DB_DIR" \
    --cazyme \
    --merops \
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
    --rename-table "$OUT_DIR/rename_table.tsv" \
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
