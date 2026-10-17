#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_busco_ab.sh -- A/B harness for Augustus training modes (v0.9.0).
#
# Arms: stock | self | busco. Identical inputs; the ONLY variable is
# --augustus-training. Reuses the Stage-1 matched rep1 masked FASTA + rename
# table (so sort/mask are not re-run and not re-measured), runs ONLY `predict`,
# and scores predict_out/consensus.gff3 against the reference (coding-only).
#
# Args: $1=genome_id  $2=arm(stock|self|busco)  $3=benchmark_dir  $4=results_dir
#
# Race-safety: each job copies AUGUSTUS_CONFIG_PATH to its own Condor scratch
# (48M) so concurrent self/busco jobs never collide on the shared user config,
# and the shared config is left pristine for other benchmarks.
# -----------------------------------------------------------------------------
set -euo pipefail

GENOME_ID="$1"
ARM="$2"
BENCHMARK_DIR="$3"
RESULTS_DIR="$4"

# -- Locate reused Stage-1 inputs (matched first, plain fallback) ------------
S1M="$BENCHMARK_DIR/results_timing_stage1_matched/myconote/$GENOME_ID/rep1"
S1P="$BENCHMARK_DIR/results_timing_stage1/myconote/$GENOME_ID/rep1"
if [[ -f "$S1M/masked.fa" ]]; then
    SRC="$S1M"
elif [[ -f "$S1P/masked.fa" ]]; then
    SRC="$S1P"
else
    echo "ERROR: no masked.fa for $GENOME_ID under matched or plain Stage-1" >&2
    exit 1
fi
MASKED="$SRC/masked.fa"
RENAME="$SRC/rename_table.tsv"
REFERENCE="$BENCHMARK_DIR/data/$GENOME_ID/reference.gff3"
for f in "$MASKED" "$RENAME" "$REFERENCE"; do
    [[ -f "$f" ]] || { echo "ERROR: missing input $f" >&2; exit 1; }
done

# -- Runtime tool environment (mirrors jobs/run_myconote.sh) -----------------
export AUGUSTUS_CONFIG_PATH_SRC="${AUGUSTUS_CONFIG_PATH:-$HOME/.myconote/augustus_config}"
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
    set +u                      # conda activation hooks read unset vars
    # shellcheck source=/dev/null
    source "$CONDA_PROFILE"
    conda activate "$MYCONOTE_ENV" 2>/dev/null \
        || conda activate "$HOME/.conda/envs/$MYCONOTE_ENV" 2>/dev/null \
        || echo "WARN: conda activate $MYCONOTE_ENV failed; using bin PATH fallback" >&2
    set -u
fi
# Robust PATH fallback (activation silently no-ops on nodes without the base
# profile): prepend the env bins directly. myconote env first (busco, blast,
# hmmsearch, snap, diamond, minimap2), then the working augustus env, then the
# freshly built v0.9.0 binary.
if [[ -d "$HOME/.conda/envs/$MYCONOTE_ENV/bin" ]]; then
    export PATH="$HOME/.conda/envs/$MYCONOTE_ENV/bin:$PATH"
fi
for _aug in "$HOME/.conda/envs/$AUGUSTUS_ENV/bin" \
            "/opt/bifxapps/miniconda3/envs/$AUGUSTUS_ENV/bin"; do
    if [[ -x "$_aug/augustus" ]]; then
        export PATH="$_aug:$PATH"
        export AUGUSTUS_SCRIPTS_PATH="$_aug"   # new_species.pl / etraining live here
        break
    fi
done
if [[ -x "$HOME/.local/bin/myconote-cli" ]]; then
    export PATH="$HOME/.local/bin:$PATH"
fi
# GeneMark-ES -- same provisioning as the baseline run_myconote.sh, so all three
# arms get GeneMark identically and differ ONLY in the Augustus model.
export GENEMARK_PATH="${GENEMARK_PATH:-$HOME/.myconote/tools/gmes_linux_64_4}"
if [[ -d "$GENEMARK_PATH" ]]; then
    export PATH="$GENEMARK_PATH:$PATH"
fi
# BBMap for BUSCO genome mode. The `busco` wrapper in the myconote env execs the
# myconote_busco env's python+busco; BUSCO's BBTools runner shells out to
# `stats.sh` (name="bbtools", cmd="stats.sh") by bare PATH lookup. stats.sh lives
# ONLY in the myconote_busco env, which was NOT on PATH, so busco genome mode
# aborted with "bbtools tool cannot be found" and `--augustus-training busco`
# FELL BACK to stock on every genome. Append the env at the TAIL so stats.sh and
# its java resolve, WITHOUT shadowing augustus (that env ships 3.5.0; stock/self
# use myconote_augustus 3.3.3 -- the trained and stock arms must share one augustus
# for a fair comparison) or tblastn/busco/snap, which stay in the higher-priority
# envs prepended above. No java exists in any higher-priority env, so stats.sh
# uses myconote_busco's java (BBTools 39.81, verified).
if [[ -d "$HOME/.conda/envs/myconote_busco/bin" ]]; then
    export PATH="$PATH:$HOME/.conda/envs/myconote_busco/bin"
fi

# -- Per-job writable Augustus config copy (race-free) -----------------------
SCRATCH="${_CONDOR_SCRATCH_DIR:-/tmp}/busco_ab_${ARM}_${GENOME_ID}_$$"
mkdir -p "$SCRATCH"
cp -a "$AUGUSTUS_CONFIG_PATH_SRC" "$SCRATCH/augustus_config"
export AUGUSTUS_CONFIG_PATH="$SCRATCH/augustus_config"
cleanup() { rm -rf "$SCRATCH"; }
trap cleanup EXIT

# -- Genome metadata from configs/genomes.tsv --------------------------------
CONFIG="$BENCHMARK_DIR/configs/genomes.tsv"
GENOME_LINE=$(awk -F'\t' -v id="$GENOME_ID" '$1==id {print; exit}' "$CONFIG")
[[ -n "$GENOME_LINE" ]] || { echo "ERROR: genome '$GENOME_ID' not in $CONFIG" >&2; exit 1; }
ORGANISM=$(echo "$GENOME_LINE" | cut -f2)
KINGDOM=$(echo "$GENOME_LINE" | cut -f3)
GENETIC_CODE=$(echo "$GENOME_LINE" | cut -f5)
BUSCO_LINEAGE_TSV=$(echo "$GENOME_LINE" | cut -f6)
LOCUS_PREFIX=$(echo "$GENOME_ID" | tr '[:lower:]' '[:upper:]')

# Clade-appropriate Augustus species -- verbatim from jobs/run_myconote.sh.
case "$GENOME_ID" in
    sce)  AUGUSTUS_SPECIES="saccharomyces_cerevisiae_S288C" ;;
    cal)  AUGUSTUS_SPECIES="candida_albicans" ;;
    ylp)  AUGUSTUS_SPECIES="yarrowia_lipolytica" ;;
    ani)  AUGUSTUS_SPECIES="aspergillus_nidulans" ;;
    ncr)  AUGUSTUS_SPECIES="neurospora_crassa" ;;
    cne)  AUGUSTUS_SPECIES="cryptococcus_neoformans_neoformans_JEC21" ;;
    klac) AUGUSTUS_SPECIES="kluyveromyces_lactis" ;;
    kpha) AUGUSTUS_SPECIES="pichia_stipitis" ;;
    fgra) AUGUSTUS_SPECIES="fusarium_graminearum" ;;
    afum) AUGUSTUS_SPECIES="aspergillus_fumigatus" ;;
    umay) AUGUSTUS_SPECIES="ustilago_maydis" ;;
    ccin) AUGUSTUS_SPECIES="coprinus_cinereus" ;;
    spom) AUGUSTUS_SPECIES="schizosaccharomyces_pombe" ;;
    tmel) AUGUSTUS_SPECIES="aspergillus_nidulans" ;;
    *)    echo "ERROR: no Augustus species mapping for '$GENOME_ID'" >&2; exit 1 ;;
esac

# -- Arm -> training flag -----------------------------------------------------
TRAIN_ARGS=()
case "$ARM" in
    stock) : ;;                                   # default stock pretrained model
    self)  TRAIN_ARGS=(--augustus-training self) ;;
    busco) TRAIN_ARGS=(--augustus-training busco) ;;
    *)     echo "ERROR: unknown arm '$ARM' (stock|self|busco)" >&2; exit 1 ;;
esac

OUT="$RESULTS_DIR/$ARM/$GENOME_ID"
rm -rf "$OUT"; mkdir -p "$OUT"          # condor .log/.out/.err live in RESULTS_DIR/logs, untouched
POUT="$SCRATCH/predict_out"             # heavy BUSCO/Augustus I/O on local scratch
LOG="$OUT/run.log"
METRICS="$OUT/metrics.json"

echo "========================================================"
echo "BUSCO-train A/B -- arm=$ARM  genome=$ORGANISM ($GENOME_ID)"
echo "Kingdom=$KINGDOM  Species=$AUGUSTUS_SPECIES  Code=$GENETIC_CODE"
echo "busco_lineage(tsv, informational only)=$BUSCO_LINEAGE_TSV"
echo "masked=$MASKED"
echo "AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH (per-job copy)"
echo "myconote-cli: $(myconote-cli --version 2>/dev/null | grep -oE 'v[0-9.]+' | head -1)"
echo "busco: $(busco --version 2>/dev/null)"
echo "========================================================"

# -- Predict (the only measured stage) ---------------------------------------
echo "[$(date +%T)] predict (arm=$ARM)"
set +e
/usr/bin/time -v -o "$OUT/time_predict.log" \
    myconote-cli predict "$MASKED" \
        --kingdom "$KINGDOM" \
        --species "$AUGUSTUS_SPECIES" \
        --locus-prefix "$LOCUS_PREFIX" \
        --genetic-code "$GENETIC_CODE" \
        --output "$POUT" \
        --threads 8 \
        ${TRAIN_ARGS[@]+"${TRAIN_ARGS[@]}"} \
        > "$LOG" 2>&1
PREDICT_RC=$?
set -e
echo "predict rc=$PREDICT_RC" | tee -a "$LOG"

# -- Copy artifacts + BUSCO-training evidence back to shared results ---------
if [[ -f "$POUT/consensus.gff3" ]]; then
    cp -f "$POUT/consensus.gff3" "$OUT/consensus.gff3"
fi
if [[ -d "$POUT/busco_train" ]]; then
    find "$POUT/busco_train" -name "short_summary*.txt" -exec cp -f {} "$OUT/" \; 2>/dev/null || true
    find "$POUT/busco_train" -name "full_table.tsv" -exec cp -f {} "$OUT/busco_full_table.tsv" \; 2>/dev/null || true
    if [[ -f "$POUT/busco_train/busco_training.gff3" ]]; then
        cp -f "$POUT/busco_train/busco_training.gff3" "$OUT/busco_training.gff3"
    fi
fi

# -- Score predict output vs reference (coding-only; lift seqids via rename) --
echo "[$(date +%T)] score consensus.gff3 vs reference"
if [[ -f "$OUT/consensus.gff3" ]]; then
    python3 "$BENCHMARK_DIR/scripts/compare_annotations.py" \
        "$OUT/consensus.gff3" \
        "$REFERENCE" \
        --rename-table "$RENAME" \
        --label "${ARM}_${GENOME_ID}" \
        --output "$METRICS" \
        2>> "$LOG" || echo "WARN: scorer failed" | tee -a "$LOG"
else
    echo "ERROR: no consensus.gff3 produced -- cannot score" | tee -a "$LOG"
fi

# -- Training-engagement verdict (the half-point of this A/B) ----------------
VERDICT="$OUT/train_verdict.txt"
{
  echo "genome=$GENOME_ID arm=$ARM"
  if [[ "$ARM" == "busco" ]]; then
      if grep -q "Trained species (BUSCO-seeded):" "$LOG"; then
          nseed=$(grep -oE "Built BUSCO training set: [0-9]+" "$LOG" | grep -oE "[0-9]+" | head -1)
          sp=$(grep "Trained species (BUSCO-seeded):" "$LOG" | tail -1 | sed 's/.*: //')
          echo "verdict=TRAINED busco_training_genes=${nseed:-NA} trained_species=${sp:-NA}"
      elif grep -q "BUSCO-SEEDED TRAINING DID NOT RUN" "$LOG"; then
          reason=$(grep -E "^\s*reason\s*:" "$LOG" | tail -1 | sed 's/.*reason[[:space:]]*:[[:space:]]*//')
          echo "verdict=FELL_BACK reason=${reason:-unknown}"
      else
          echo "verdict=UNKNOWN (no train/fallback marker in run.log)"
      fi
  elif [[ "$ARM" == "self" ]]; then
      if grep -q "Trained species:" "$LOG"; then
          echo "verdict=TRAINED (self)"
      elif grep -q "SELF-TRAINING DID NOT RUN" "$LOG"; then
          reason=$(grep -E "^\s*reason\s*:" "$LOG" | tail -1 | sed 's/.*reason[[:space:]]*:[[:space:]]*//')
          echo "verdict=FELL_BACK reason=${reason:-unknown}"
      else
          echo "verdict=UNKNOWN"
      fi
  else
      echo "verdict=STOCK (no training requested)"
  fi
  echo "predict_rc=$PREDICT_RC"
} > "$VERDICT"
cat "$VERDICT"

echo "[$(date +%T)] done -- $OUT"
exit 0
