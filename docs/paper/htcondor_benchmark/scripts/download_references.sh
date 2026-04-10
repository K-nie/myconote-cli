#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# download_references.sh
# Downloads all 8 reference genomes and curated annotations for the benchmark.
#
# Usage:
#   bash download_references.sh [data_dir]
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

DATA_DIR="${1:-./benchmark_data}"
CONFIG="$(dirname "$0")/../configs/genomes.tsv"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

mkdir -p "$DATA_DIR"

# Skip header row
tail -n +2 "$CONFIG" | while IFS=$'\t' read -r id organism kingdom size code busco genome_url anno_url; do
    GENOME_DIR="$DATA_DIR/$id"
    mkdir -p "$GENOME_DIR"

    GENOME_FILE="$GENOME_DIR/genome.fa"
    ANNO_FILE="$GENOME_DIR/reference.gff3"

    if [[ -f "$GENOME_FILE" ]] && [[ -f "$ANNO_FILE" ]]; then
        echo "[skip] $id ($organism) — already present"
        continue
    fi

    echo ""
    echo "── $id ($organism) ──"

    # Download genome
    if [[ ! -f "$GENOME_FILE" ]]; then
        echo "  Downloading genome from $genome_url"
        TMP_GENOME="$GENOME_DIR/genome.tmp"
        if curl -fsSL "$genome_url" -o "$TMP_GENOME" 2>&1; then
            # Decompress if gzipped
            if file "$TMP_GENOME" | grep -q gzip; then
                gunzip -c "$TMP_GENOME" > "$GENOME_FILE"
                rm "$TMP_GENOME"
            else
                mv "$TMP_GENOME" "$GENOME_FILE"
            fi
            echo "  ✓ Genome saved to $GENOME_FILE"
        else
            echo "  ✗ Failed to download genome for $id"
            continue
        fi
    fi

    # Download annotation
    if [[ ! -f "$ANNO_FILE" ]]; then
        echo "  Downloading annotation from $anno_url"
        TMP_ANNO="$GENOME_DIR/anno.tmp"
        if curl -fsSL "$anno_url" -o "$TMP_ANNO" 2>&1; then
            if file "$TMP_ANNO" | grep -q gzip; then
                gunzip -c "$TMP_ANNO" > "$ANNO_FILE"
                rm "$TMP_ANNO"
            else
                mv "$TMP_ANNO" "$ANNO_FILE"
            fi
            echo "  ✓ Annotation saved to $ANNO_FILE"
        else
            echo "  ✗ Failed to download annotation for $id"
            continue
        fi
    fi

    # Report sizes
    echo "  Genome size: $(du -h "$GENOME_FILE" | cut -f1)"
    echo "  Annotation size: $(du -h "$ANNO_FILE" | cut -f1)"
done

echo ""
echo "──────────────────────────────────────────"
echo "Download complete. Data dir: $DATA_DIR"
echo "Total size: $(du -sh "$DATA_DIR" | cut -f1)"
