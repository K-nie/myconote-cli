#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# download_references.sh
# Downloads all reference genomes + curated annotations from NCBI RefSeq.
#
# Uses NCBI accessions (GCF_*) from configs/genomes.tsv column 7, auto-
# discovers the assembly folder from the NCBI FTP directory index (folder
# names include an assembly label that varies per assembly), and validates
# each downloaded file for plausible content and size before accepting it.
#
# Usage:
#   bash download_references.sh [data_dir]
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

DATA_DIR="${1:-./benchmark_data}"
CONFIG="$(cd "$(dirname "$0")/.." && pwd)/configs/genomes.tsv"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG" >&2
    exit 1
fi

mkdir -p "$DATA_DIR"

NCBI_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all"

# ── Helpers ──────────────────────────────────────────────────────────────

file_size() {
    # Portable file size (works on Linux + BSD/macOS)
    if stat -c%s "$1" &>/dev/null; then
        stat -c%s "$1"
    else
        stat -f%z "$1"
    fi
}

# Parse GCF_000146045.2 -> GCF/000/146/045
accession_path() {
    local acc="$1"
    local prefix="${acc:0:3}"
    local digits="${acc:4}"
    digits="${digits%%.*}"         # strip .version
    # Left-pad to 9 digits just in case
    while [[ "${#digits}" -lt 9 ]]; do digits="0${digits}"; done
    echo "${prefix}/${digits:0:3}/${digits:3:3}/${digits:6:3}"
}

# Scrape NCBI parent directory index to find the full assembly folder name
find_ncbi_subdir() {
    local acc="$1"
    local parent_url="${NCBI_BASE}/$(accession_path "$acc")/"
    local html
    html=$(curl -fsSL --retry 3 --retry-delay 5 "$parent_url" 2>/dev/null) || return 1
    # Folders look like: <a href="GCF_000146045.2_R64/">
    echo "$html" \
        | grep -oE "${acc}_[A-Za-z0-9._-]+/" \
        | head -1 \
        | sed 's:/$::'
}

validate_fasta() {
    local f="$1"
    local min_bytes="$2"
    # Use wc -c for size — portable, follows symlinks, immune to NFS
    # attribute-cache lag that can confuse `stat` and `du` on GLBRC.
    local sz
    sz=$(wc -c < "$f" 2>/dev/null || echo 0)
    if [[ "$sz" -lt "$min_bytes" ]]; then
        echo "    ✗ FASTA too small: $sz bytes (expected >= $min_bytes)"
        return 1
    fi
    if ! head -1 "$f" | grep -q '^>'; then
        echo "    ✗ Not a valid FASTA (first line does not start with '>')"
        echo "    ↳ first line: $(head -1 "$f" | cut -c1-80)"
        return 1
    fi
    return 0
}

validate_gff3() {
    local f="$1"
    local min_bytes="$2"
    local sz
    sz=$(wc -c < "$f" 2>/dev/null || echo 0)
    if [[ "$sz" -lt "$min_bytes" ]]; then
        echo "    ✗ GFF3 too small: $sz bytes (expected >= $min_bytes)"
        return 1
    fi
    # Note: `head -N | grep` is a trap with `set -o pipefail`. When grep
    # finds a match early and exits, head gets SIGPIPE and the pipeline is
    # reported as failed — even though the file is perfectly valid. Instead,
    # run grep directly on the file with -m1 so it stops after the first
    # match with exit 0 and no pipe involved.
    if ! grep -qm1 '^##gff-version' "$f"; then
        echo "    ✗ Not a valid GFF3 (no ##gff-version header anywhere in file)"
        echo "    ↳ first line: $(head -1 "$f" | cut -c1-80)"
        return 1
    fi
    if ! grep -qm1 '^[^#]' "$f"; then
        echo "    ✗ GFF3 contains only comment lines (no feature records)"
        return 1
    fi
    return 0
}

# Download a URL with retries, resume support, and cleanup on failure.
download_with_retry() {
    local url="$1"
    local dest="$2"
    rm -f "$dest"
    if curl -fSL --retry 5 --retry-delay 10 --connect-timeout 30 \
            -o "$dest" "$url"; then
        return 0
    else
        rm -f "$dest"
        return 1
    fi
}

# ── Main loop ────────────────────────────────────────────────────────────

total=0
ok_count=0
failed_ids=()

# Read config into arrays so we can iterate without a subshell
ids=()
organisms=()
size_mbs=()
accs=()
while IFS=$'\t' read -r id organism kingdom size_mb code busco acc; do
    [[ "$id" == "id" ]] && continue
    [[ "$id" == \#* ]] && continue   # skip comment lines
    [[ -z "$id" ]] && continue
    ids+=("$id")
    organisms+=("$organism")
    size_mbs+=("$size_mb")
    accs+=("$acc")
done < "$CONFIG"

for i in "${!ids[@]}"; do
    id="${ids[$i]}"
    organism="${organisms[$i]}"
    size_mb="${size_mbs[$i]}"
    acc="${accs[$i]}"
    total=$((total + 1))

    GENOME_DIR="$DATA_DIR/$id"
    mkdir -p "$GENOME_DIR"
    GENOME_FILE="$GENOME_DIR/genome.fa"
    ANNO_FILE="$GENOME_DIR/reference.gff3"

    # Minimum expected sizes (bytes): genome ~25% of expected uncompressed size,
    # GFF3 >= 100 kB (NCBI RefSeq annotations with thousands of features easily
    # exceed this; anything smaller is probably an error page).
    min_genome_bytes=$(( size_mb * 256 * 1024 ))   # 0.25 * size_mb * 1MB
    min_gff_bytes=$(( 100 * 1024 ))

    # Skip if already present AND valid
    if [[ -f "$GENOME_FILE" && -f "$ANNO_FILE" ]]; then
        if validate_fasta "$GENOME_FILE" "$min_genome_bytes" &>/dev/null \
           && validate_gff3 "$ANNO_FILE" "$min_gff_bytes" &>/dev/null; then
            echo "[skip] $id ($organism) — already present and valid"
            ok_count=$((ok_count + 1))
            continue
        else
            echo "[revalidate] $id — existing files invalid, re-downloading"
            rm -f "$GENOME_FILE" "$ANNO_FILE"
        fi
    fi

    echo ""
    echo "── $id ($organism) — $acc ──"

    # Resolve the full NCBI assembly folder
    echo "  Resolving NCBI folder..."
    subdir=$(find_ncbi_subdir "$acc" || true)
    if [[ -z "$subdir" ]]; then
        echo "  ✗ Could not locate $acc on NCBI"
        failed_ids+=("$id")
        continue
    fi
    echo "  → $subdir"

    parent="${NCBI_BASE}/$(accession_path "$acc")/${subdir}"
    genome_url="${parent}/${subdir}_genomic.fna.gz"
    anno_url="${parent}/${subdir}_genomic.gff.gz"

    # ── Genome ──
    echo "  Downloading genome from $genome_url"
    TMP_GENOME="$GENOME_DIR/genome.fna.gz"
    if ! download_with_retry "$genome_url" "$TMP_GENOME"; then
        echo "  ✗ curl failed for genome"
        failed_ids+=("$id")
        continue
    fi
    if ! gunzip -t "$TMP_GENOME" &>/dev/null; then
        echo "  ✗ genome .gz file is corrupt"
        rm -f "$TMP_GENOME"
        failed_ids+=("$id")
        continue
    fi
    gunzip -c "$TMP_GENOME" > "$GENOME_FILE"
    rm -f "$TMP_GENOME"
    if ! validate_fasta "$GENOME_FILE" "$min_genome_bytes"; then
        rm -f "$GENOME_FILE"
        failed_ids+=("$id")
        continue
    fi
    echo "    ✓ genome.fa ($(wc -c < "$GENOME_FILE") bytes)"

    # ── Annotation ──
    echo "  Downloading annotation from $anno_url"
    TMP_ANNO="$GENOME_DIR/reference.gff.gz"
    if ! download_with_retry "$anno_url" "$TMP_ANNO"; then
        echo "  ✗ curl failed for annotation"
        rm -f "$GENOME_FILE"
        failed_ids+=("$id")
        continue
    fi
    if ! gunzip -t "$TMP_ANNO" &>/dev/null; then
        echo "  ✗ annotation .gz file is corrupt"
        rm -f "$TMP_ANNO" "$GENOME_FILE"
        failed_ids+=("$id")
        continue
    fi
    gunzip -c "$TMP_ANNO" > "$ANNO_FILE"
    rm -f "$TMP_ANNO"
    if ! validate_gff3 "$ANNO_FILE" "$min_gff_bytes"; then
        rm -f "$GENOME_FILE" "$ANNO_FILE"
        failed_ids+=("$id")
        continue
    fi
    echo "    ✓ reference.gff3 ($(wc -c < "$ANNO_FILE") bytes)"

    ok_count=$((ok_count + 1))
done

echo ""
echo "──────────────────────────────────────────"
echo "Downloaded and validated: $ok_count / $total"
if [[ "${#failed_ids[@]}" -gt 0 ]]; then
    echo "Failed: ${failed_ids[*]}"
fi
echo "Data dir: $DATA_DIR"
if command -v du &>/dev/null; then
    echo "Total size: $(du -sh "$DATA_DIR" 2>/dev/null | cut -f1)"
fi

if [[ "${#failed_ids[@]}" -gt 0 ]]; then
    exit 2
fi
