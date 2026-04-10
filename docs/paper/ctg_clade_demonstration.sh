#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# CTG Clade Genetic Code Demonstration
#
# Demonstrates that annotating Candida tropicalis with the standard genetic
# code (Table 1) versus the alternative yeast nuclear code (Table 12)
# produces measurably different protein sequences.
#
# This addresses Reviewer #2 Concern #6: "Show concretely that annotating
# Candida albicans with --genetic-code 12 produces different (better)
# results than with the standard code."
#
# Usage:
#   bash ctg_clade_demonstration.sh
#
# Output:
#   ctg_demo/standard_code/  - Annotation with Table 1 (incorrect)
#   ctg_demo/ctg_code/       - Annotation with Table 12 (correct)
#   ctg_demo/comparison.txt  - Quantitative comparison
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

OUTDIR="ctg_demo"
GFF="tests/data/candida_tropicalis.final.gff3"
FASTA="tests/data/candida_tropicalis.fas"

# Verify inputs
if [[ ! -f "$GFF" ]] || [[ ! -f "$FASTA" ]]; then
    echo "ERROR: Test data not found. Run from the myconote-cli repository root."
    exit 1
fi

mkdir -p "$OUTDIR"

# ── Run 1: Standard code (Table 1) - INCORRECT for CTG clade ────────────
echo "=========================================="
echo "Run 1: Standard genetic code (Table 1)"
echo "Expected: Each CTG translated as Leucine (L)"
echo "=========================================="

myconote-cli annotate "$GFF" \
    --fasta "$FASTA" \
    --output "$OUTDIR/standard_code" \
    --kingdom fungi \
    --locus-prefix CTRO_STD \
    --genetic-code 1 \
    --no-busco --no-pfam \
    --threads 4 2>&1 | tail -20

echo ""

# ── Run 2: Alternative yeast nuclear code (Table 12) - CORRECT ─────────
echo "=========================================="
echo "Run 2: Alternative yeast nuclear (Table 12)"
echo "Expected: Each CTG translated as Serine (S)"
echo "=========================================="

myconote-cli annotate "$GFF" \
    --fasta "$FASTA" \
    --output "$OUTDIR/ctg_code" \
    --kingdom fungi \
    --locus-prefix CTRO_CTG \
    --genetic-code 12 \
    --no-busco --no-pfam \
    --threads 4 2>&1 | tail -20

echo ""
echo "=========================================="
echo "Quantitative comparison"
echo "=========================================="

PROTEINS_STD="$OUTDIR/standard_code/proteins.fa"
PROTEINS_CTG="$OUTDIR/ctg_code/proteins.fa"

# Verify both protein files exist
if [[ ! -f "$PROTEINS_STD" ]] || [[ ! -f "$PROTEINS_CTG" ]]; then
    echo "ERROR: Protein output files not generated"
    exit 1
fi

REPORT="$OUTDIR/comparison.txt"

{
    echo "CTG CLADE GENETIC CODE DEMONSTRATION"
    echo "======================================"
    echo "Generated: $(date)"
    echo ""
    echo "Test organism: Candida tropicalis (CTG clade)"
    echo "Input GFF3:    $GFF"
    echo "Input FASTA:   $FASTA"
    echo ""

    # Total proteins
    n_std=$(grep -c "^>" "$PROTEINS_STD")
    n_ctg=$(grep -c "^>" "$PROTEINS_CTG")
    echo "Proteins (Table 1):  $n_std"
    echo "Proteins (Table 12): $n_ctg"
    echo ""

    # Check that protein counts match
    if [[ "$n_std" != "$n_ctg" ]]; then
        echo "WARNING: Protein counts differ between runs"
    fi

    # Total residues
    res_std=$(grep -v "^>" "$PROTEINS_STD" | tr -d '\n' | wc -c | tr -d ' ')
    res_ctg=$(grep -v "^>" "$PROTEINS_CTG" | tr -d '\n' | wc -c | tr -d ' ')
    echo "Total residues (Table 1):  $res_std"
    echo "Total residues (Table 12): $res_ctg"
    echo ""

    # Count L vs S in each
    l_std=$(grep -v "^>" "$PROTEINS_STD" | tr -d '\n' | tr -dc 'L' | wc -c | tr -d ' ')
    s_std=$(grep -v "^>" "$PROTEINS_STD" | tr -d '\n' | tr -dc 'S' | wc -c | tr -d ' ')
    l_ctg=$(grep -v "^>" "$PROTEINS_CTG" | tr -d '\n' | tr -dc 'L' | wc -c | tr -d ' ')
    s_ctg=$(grep -v "^>" "$PROTEINS_CTG" | tr -d '\n' | tr -dc 'S' | wc -c | tr -d ' ')

    echo "Amino acid counts:"
    echo "                Table 1   Table 12   Diff"
    echo "  Leucine  (L): $l_std    $l_ctg    $((l_std - l_ctg))"
    echo "  Serine   (S): $s_std    $s_ctg    $((s_ctg - s_std))"
    echo ""
    echo "Expected: Table 12 should have FEWER L and MORE S compared to Table 1."
    echo "          The difference equals the number of CTG codons in CDSs."
    echo ""

    # Per-protein comparison
    if command -v python3 &>/dev/null; then
        echo "Per-protein analysis:"
        python3 << PYEOF
def read_fasta(path):
    seqs = {}
    cur = None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line
    return seqs

std = read_fasta("$PROTEINS_STD")
ctg = read_fasta("$PROTEINS_CTG")

common = set(std) & set(ctg)
diffs = []
for k in common:
    if len(std[k]) != len(ctg[k]):
        continue
    n_diff = sum(1 for a, b in zip(std[k], ctg[k]) if a != b)
    diffs.append(n_diff)

if diffs:
    diffs.sort()
    n = len(diffs)
    mean_diff = sum(diffs) / n
    median_diff = diffs[n // 2]
    max_diff = max(diffs)
    n_changed = sum(1 for d in diffs if d > 0)
    print(f"  Proteins compared:        {n}")
    print(f"  Proteins with changes:    {n_changed} ({n_changed/n*100:.1f}%)")
    print(f"  Mean changes per protein: {mean_diff:.2f}")
    print(f"  Median changes:           {median_diff}")
    print(f"  Maximum changes:          {max_diff}")
    print()
    print("  Distribution of changes per protein:")
    bins = [0, 1, 5, 10, 25, 50, 1000]
    for lo, hi in zip(bins, bins[1:]):
        count = sum(1 for d in diffs if lo <= d < hi)
        if count > 0:
            print(f"    {lo:>4}-{hi:<4}: {count:>5} proteins")
PYEOF
    fi

    echo ""
    echo "INTERPRETATION:"
    echo "  - Each difference is a CTG codon translated as L (Table 1) vs S (Table 12)"
    echo "  - L is hydrophobic (Kyte-Doolittle: 3.8)"
    echo "  - S is polar, hydrogen-bond donor (Kyte-Doolittle: -0.8)"
    echo "  - The substitution is biologically significant: protein folding,"
    echo "    domain architecture, and functional predictions all change"
    echo ""
    echo "CONCLUSION:"
    echo "  Annotating Candida species with the standard code (Table 1) silently"
    echo "  produces incorrect protein sequences. MycoNote-CLI's --genetic-code 12"
    echo "  flag prevents this error. Other tools support non-standard codes through"
    echo "  configuration files, but MycoNote-CLI exposes this through a single CLI flag."

} | tee "$REPORT"

echo ""
echo "Full report written to: $REPORT"
