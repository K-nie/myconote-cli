#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# myconote-cli test runner
# Run from the project root:  bash test.sh
# ─────────────────────────────────────────────────────────────────────────────

set -e
BIN="./target/release/myconote-cli"
GFF="tests/data/candida_tropicalis.final.gff3"
FA="tests/data/candida_tropicalis.fas"
TMP=$(mktemp -d)
PASS=0; FAIL=0

GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; NC='\033[0m'

ok()   { echo -e "${GREEN}✓${NC} $1"; PASS=$((PASS+1)); }
fail() { echo -e "${RED}✗${NC} $1"; FAIL=$((FAIL+1)); }
hdr()  { echo -e "\n${YELLOW}── $1 ──────────────────────────────────${NC}"; }

# ── 0. Build ──────────────────────────────────────────────────────────────────
hdr "Build"
if cargo build --release 2>&1; then
    ok "cargo build --release"
else
    echo -e "${RED}Build failed — aborting tests.${NC}"
    exit 1
fi

# ── 1. Help / top-level ───────────────────────────────────────────────────────
hdr "Help & routing"

$BIN 2>&1 | grep -q "Genome Annotation Pipeline\|Usage:"  && ok "help: no args shows usage" || fail "help: no args"
$BIN check 2>&1 | grep -qi "tool\|OK\|MISSING"  && ok "check: runs"       || fail "check"
$BIN setup --list 2>&1 | grep -q "swiss-prot"   && ok "setup --list"      || fail "setup --list"
$BIN species --grouped 2>&1 | grep -qi "fungi"  && ok "species --grouped" || fail "species --grouped"

# ── 2. stats ──────────────────────────────────────────────────────────────────
hdr "stats"

$BIN stats "$GFF" 2>&1 | grep -qE "[0-9]+"                && ok "stats: human output"    || fail "stats: human output"
$BIN stats "$GFF" --format json 2>&1 | grep -q "gene_count" && ok "stats: json"           || fail "stats: json"
$BIN stats "$GFF" --format csv  2>&1 | grep -q ","          && ok "stats: csv"            || fail "stats: csv"
$BIN stats "$GFF" --taxon fungi 2>&1 | grep -qE "[0-9]+"   && ok "stats: --taxon fungi"  || fail "stats: --taxon fungi"

# ── 3. sort ───────────────────────────────────────────────────────────────────
hdr "sort"

$BIN sort "$FA" -o "$TMP/sorted.fa" && \
    grep -c "^>" "$TMP/sorted.fa" | grep -q "^24$" && \
    ok "sort: all 24 scaffolds preserved" || fail "sort: sequence count"

$BIN sort "$FA" --prefix chr -o "$TMP/sorted_chr.fa" && \
    grep -q "^>chr_" "$TMP/sorted_chr.fa" && \
    ok "sort: --prefix chr" || fail "sort: --prefix"

$BIN sort "$FA" --min-length 500000 -o "$TMP/sorted_long.fa" && \
    N=$(grep -c "^>" "$TMP/sorted_long.fa") && \
    [ "$N" -lt 24 ] && \
    ok "sort: --min-length filters short scaffolds ($N of 24 kept)" || fail "sort: --min-length"

$BIN sort "$FA" --rename-table "$TMP/rename.tsv" -o "$TMP/sorted_tbl.fa" && \
    [ -f "$TMP/rename.tsv" ] && \
    ok "sort: --rename-table created" || fail "sort: --rename-table"

# ── 4. clean ──────────────────────────────────────────────────────────────────
hdr "clean"

$BIN clean "$GFF" -o "$TMP/clean.gff3" && \
    grep -c $'\tgene\t' "$TMP/clean.gff3" | awk '{exit ($1 < 6000)}' && \
    ok "clean: gene count preserved" || fail "clean: gene count"

$BIN clean "$GFF" --remove-orphans -o "$TMP/clean_orphans.gff3" && \
    ok "clean: --remove-orphans" || fail "clean: --remove-orphans"

# ── 5. convert ────────────────────────────────────────────────────────────────
hdr "convert"

$BIN convert "$GFF" --to gtf -o "$TMP/out.gtf" && \
    grep -q "transcript_id" "$TMP/out.gtf" && \
    ok "convert: GFF3 → GTF" || fail "convert: GFF3 → GTF"

$BIN convert "$GFF" --to bed -o "$TMP/out.bed" && \
    awk '{exit (NF != 6)}' "$TMP/out.bed" && \
    ok "convert: GFF3 → BED6 (6 fields)" || fail "convert: GFF3 → BED6"

$BIN convert "$GFF" --to bed12 -o "$TMP/out.bed12" && \
    ok "convert: GFF3 → BED12" || fail "convert: GFF3 → BED12"

$BIN convert "$GFF" --to table -o "$TMP/out.tsv" && \
    grep -q $'\t' "$TMP/out.tsv" && \
    ok "convert: GFF3 → TSV table" || fail "convert: GFF3 → TSV table"

$BIN convert "$GFF" --to protein --fasta "$FA" -o "$TMP/proteins.faa" && \
    N=$(grep -c "^>" "$TMP/proteins.faa") && \
    [ "$N" -gt 100 ] && \
    ok "convert: GFF3 → protein ($N sequences)" || fail "convert: GFF3 → protein"

$BIN convert "$GFF" --to genbank --fasta "$FA" -o "$TMP/out.gbk" && \
    grep -q "^LOCUS" "$TMP/out.gbk" && \
    ok "convert: GFF3 → GenBank (LOCUS line present)" || fail "convert: GFF3 → GenBank"

# ── 6. fix ────────────────────────────────────────────────────────────────────
hdr "fix"

# Need a GenBank to fix; use the one we just created
if [ -f "$TMP/out.gbk" ]; then
    $BIN fix "$TMP/out.gbk" -o "$TMP/fixed.gbk" --report "$TMP/fix_report.txt" && \
        [ -f "$TMP/fixed.gbk" ] && grep -q "Records processed" "$TMP/fix_report.txt" && \
        ok "fix: GenBank repaired, report created" || fail "fix: GenBank repair"

    $BIN fix "$TMP/out.gbk" -o "$TMP/dry.gbk" --dry-run && \
        [ ! -f "$TMP/dry.gbk" ] && \
        ok "fix: --dry-run creates no output file" || fail "fix: --dry-run"
else
    echo "  (skipped fix tests — GenBank not created)"
fi

# ── 7. plot ───────────────────────────────────────────────────────────────────
hdr "plot"

$BIN plot "$GFF" --output "$TMP/linear.png" --type linear && \
    python3 -c "
f=open('$TMP/linear.png','rb')
h=f.read(4)
assert h==b'\x89PNG', f'Not a PNG: {h}'
print('valid PNG')
" && ok "plot: linear PNG (valid magic bytes)" || fail "plot: linear PNG"

$BIN plot "$GFF" --output "$TMP/circular.png" --type circular && \
    python3 -c "
f=open('$TMP/circular.png','rb')
h=f.read(4)
assert h==b'\x89PNG', f'Not a PNG: {h}'
print('valid PNG')
" && ok "plot: circular PNG (valid magic bytes)" || fail "plot: circular PNG"

# ── 8. view ───────────────────────────────────────────────────────────────────
hdr "view"

$BIN view "$GFF" -o "$TMP/view.html" && \
    grep -qi "jbrowse\|html" "$TMP/view.html" && \
    ok "view: JBrowse2 HTML created" || fail "view: JBrowse2 HTML"

$BIN view "$GFF" --browser ucsc -o "$TMP/ucsc.html" && \
    ok "view: UCSC HTML created" || fail "view: UCSC HTML"

# ── 9. Help text for all subcommands ─────────────────────────────────────────
hdr "Subcommand help"

for CMD in sort mask train predict update annotate remote fix convert clean synteny view; do
    $BIN $CMD 2>&1 | grep -qi "usage\|options\|error" && \
        ok "help: $CMD" || fail "help: $CMD"
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "─────────────────────────────────────────────────────────────────"
echo -e "  ${GREEN}Passed: $PASS${NC}   ${RED}Failed: $FAIL${NC}"
echo "─────────────────────────────────────────────────────────────────"
echo "  Temporary files: $TMP"
echo ""

[ "$FAIL" -eq 0 ] && exit 0 || exit 1
