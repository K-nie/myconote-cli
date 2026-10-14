#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# build_cazyme_db.sh
# One-time provisioning of a REAL dbCAN CAZyme DIAMOND database for the
# MycoNote-CLI arm of the HTCondor benchmark.
#
# WHY: the CAZyme DB previously shipped under ~/.myconote/dbs/dbcan/ was corrupt
# — CAZyDB.fa and dbCAN-HMMdb-V12.txt were both 18,618-byte HTML error pages
# (a 404/redirect captured by the July-9 2025 download), 0 sequences / 0 HMMs.
# With no valid DB and no --cazyme-db flag, `annotate --cazyme` silently no-oped.
#
# WHAT: downloads the authoritative dbCAN v5-2 CAZyDB protein FASTA from the
# public dbCAN AWS S3 bucket (CAZyDB.07242025.fa, CAZyDB release 2025-07-24),
# reformats the FASTA headers so the CAZy family is the leading '|'-delimited
# token (required by MycoNote's parse_diamond_cazyme), and builds a DIAMOND
# database CAZyDB.dmnd that the `annotate --cazyme --cazyme-db` fallback uses.
#
# COMPARABILITY: funannotate's CAZyme step runs hmmsearch against dbCAN HMMs
# (HMMdb release 14.0, 2025-08-26, derived from CAZyDB 2025-07-10). This builds
# the SAME-era CAZyDB release (2025-07-24) but drives it with DIAMOND-blastp,
# a different (sequence-similarity, not profile-HMM) method. See the benchmark
# README CAZyme caveat.
#
# MUST run on an HTCondor execute node (never the submit node).
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────────────
DBCAN_DIR="${DBCAN_DIR:-$HOME/.myconote/dbs/dbcan}"
S3_BASE="https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2_9-13-2025"
CAZYDB_FASTA_NAME="CAZyDB.07242025.fa"
CAZYDB_URL="$S3_BASE/$CAZYDB_FASTA_NAME"

RAW_FASTA="$DBCAN_DIR/$CAZYDB_FASTA_NAME"
DMND="$DBCAN_DIR/CAZyDB.dmnd"
BACKUP_DIR="$DBCAN_DIR/corrupt_html_backup"

# ── Tool environment (match run_myconote.sh: prepend myconote env bin) ─────────
MYCONOTE_ENV="${MYCONOTE_ENV:-myconote}"
if [[ -d "$HOME/.conda/envs/$MYCONOTE_ENV/bin" ]]; then
    export PATH="$HOME/.conda/envs/$MYCONOTE_ENV/bin:$PATH"
fi

echo "════════════════════════════════════════════════════════"
echo "Build dbCAN CAZyme DIAMOND DB for MycoNote benchmark"
echo "Host:     $(hostname)"
echo "Date:     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "DBCAN_DIR: $DBCAN_DIR"
echo "Source:   $CAZYDB_URL"
echo "diamond:  $(command -v diamond) ($(diamond --version 2>&1 | head -1))"
echo "════════════════════════════════════════════════════════"

mkdir -p "$DBCAN_DIR"

# ── Step 0: probe S3 reachability from this execute node ───────────────────────
echo "[$(date +%T)] Step 0: probe S3 reachability"
if ! curl -sI --max-time 30 "$CAZYDB_URL" | grep -qiE "HTTP/.* 200"; then
    echo "ERROR: cannot reach $CAZYDB_URL from this execute node (no HTTP 200)."
    echo "       Execute-node outbound internet may be blocked. STOP."
    curl -sI --max-time 30 "$CAZYDB_URL" 2>&1 | head -10 || true
    exit 3
fi
echo "  S3 reachable (HTTP 200)."

# ── Step 1: back up the corrupt HTML DB files (keep as evidence) ──────────────
echo "[$(date +%T)] Step 1: back up corrupt HTML files"
mkdir -p "$BACKUP_DIR"
for f in CAZyDB.fa dbCAN-HMMdb-V12.txt; do
    if [[ -f "$DBCAN_DIR/$f" ]]; then
        if head -c 20 "$DBCAN_DIR/$f" | grep -qi "<!DOCTYPE\|<html"; then
            mv -f "$DBCAN_DIR/$f" "$BACKUP_DIR/$f"
            echo "  moved corrupt $f -> $BACKUP_DIR/"
        else
            echo "  WARNING: $f does not look like HTML; leaving in place"
        fi
    fi
done

# ── Step 2: download the real CAZyDB protein FASTA (idempotent) ───────────────
echo "[$(date +%T)] Step 2: download $CAZYDB_FASTA_NAME"
if [[ -s "$RAW_FASTA" ]] && ! head -c 20 "$RAW_FASTA" | grep -qi "<!DOCTYPE\|<html" \
   && [[ "$(head -1 "$RAW_FASTA" | cut -c1)" == ">" ]]; then
    echo "  valid FASTA already present; skipping download"
else
    curl -L --fail --max-time 3600 -o "$RAW_FASTA" "$CAZYDB_URL"
fi
ls -l "$RAW_FASTA"
echo "  md5: $(md5sum "$RAW_FASTA" | awk '{print $1}')"

# ── Step 3: verify it is a real FASTA (not HTML) ──────────────────────────────
echo "[$(date +%T)] Step 3: verify real FASTA"
if head -c 20 "$RAW_FASTA" | grep -qi "<!DOCTYPE\|<html"; then
    echo "ERROR: downloaded file is HTML, not FASTA. STOP."
    head -c 300 "$RAW_FASTA"; echo
    exit 2
fi
echo "  head:"; head -3 "$RAW_FASTA"
NSEQ=$(grep -c "^>" "$RAW_FASTA")
echo "  sequence count (grep -c '^>'): $NSEQ"
if [[ "$NSEQ" -le 0 ]]; then
    echo "ERROR: 0 sequences in FASTA. STOP."; exit 2
fi

# ── Step 4: reformat headers CAZy-family-first for parse_diamond_cazyme ────────
# dbCAN CAZyDB mixes header formats: most are '>ACCESSION|FAMILY' (e.g.
# >QDZ53590.1|GH31_13), but the bundle also includes JGI/other entries whose
# CAZy family is NOT in field 2 (naively moving field 1 to the end then yields a
# bogus numeric "family"). MycoNote's parse_diamond_cazyme reads the family as
# the FIRST '|'-token (stripping any '_subfamily'), so we SCAN all '|'-fields for
# the first token matching a CAZy-family pattern (GH/GT/PL/CE/AA/CBM + digits)
# and make THAT the leading token:
#   >QDZ53590.1|GH31_13     -> >GH31_13|QDZ53590.1   (parser -> family "GH31")
#   >ACC|123456|GH5         -> >GH5|ACC              (parser -> family "GH5")
# Sequences with no recognizable CAZy-family token in any field are DROPPED
# (they cannot yield a valid family call anyway); counts are logged.
echo "[$(date +%T)] Step 4: reformat headers (CAZy-family-first) + diamond makedb"
REFMT="$DBCAN_DIR/CAZyDB.reformatted.fa"
awk '
    /^>/ { sub(/^>/,""); n=split($0,a,"|"); fam="";
           for(i=1;i<=n;i++){ if(a[i] ~ /^(GH|GT|PL|CE|AA|CBM)[0-9]/){ fam=a[i]; break } }
           if(fam==""){ skip=1; dropped++; next }
           skip=0; kept++; print ">"fam"|"a[1]; next }
    skip!=1 { print }
    END { printf("  reformat: kept=%d dropped=%d (no CAZy-family token)\n", kept, dropped) > "/dev/stderr" }
' "$RAW_FASTA" > "$REFMT"
echo "  sample reformatted headers:"; { grep "^>" "$REFMT" | head -5 ; } || true
echo "  leading-token family-class tally (reformatted):"
{ grep "^>" "$REFMT" | sed -E 's/^>//; s/[0-9_|].*//' | sort | uniq -c | sort -rn ; } || true
RKEPT=$(grep -c "^>" "$REFMT")
echo "  reformatted sequence count (grep -c '^>'): $RKEPT"
diamond makedb --in "$REFMT" --db "$DMND" --threads 8
rm -f "$REFMT"   # keep only the raw FASTA + the .dmnd on shared storage

# ── Step 5: verify the DIAMOND DB ─────────────────────────────────────────────
echo "[$(date +%T)] Step 5: verify DIAMOND DB"
ls -l "$DMND"
{ diamond dbinfo --db "$DMND" 2>&1 | head -20 ; } || true

# ── Step 6: write provenance stamp ────────────────────────────────────────────
echo "[$(date +%T)] Step 6: provenance stamp"
cat > "$DBCAN_DIR/dbcan.version" <<EOF
# dbCAN CAZyme DIAMOND database for MycoNote-CLI benchmark
# Provisioned: $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname)
source_url        $CAZYDB_URL
cazydb_release    2025-07-24
dbcan_db_bundle   db_v5-2_9-13-2025 (dbCAN v5-2, 2025-09-13)
raw_fasta         $RAW_FASTA
raw_fasta_md5     $(md5sum "$RAW_FASTA" | awk '{print $1}')
raw_fasta_nseq    $NSEQ
retained_nseq     ${RKEPT:-NA} (sequences with a recognizable CAZy-family token)
diamond_db        $DMND
diamond_db_md5    $(md5sum "$DMND" | awk '{print $1}')
diamond_version   $(diamond --version 2>&1 | head -1)
method            DIAMOND blastp (family from best hit); header reformatted family-first
comparability     funannotate uses dbCAN HMMdb 14.0 (CAZyDB 2025-07-10) via HMMER;
                  this DB is the same-era CAZyDB release driven by DIAMOND (different method)
EOF
cat "$DBCAN_DIR/dbcan.version"

echo "[$(date +%T)] DONE. CAZyDB.dmnd ready at $DMND"
