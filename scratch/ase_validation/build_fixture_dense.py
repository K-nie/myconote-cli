#!/usr/bin/env python3
"""
Same as build_fixture.py but injects ~20 phased SNVs per transcript
(roughly one per 60–80 bp), matching the variant density of real F1
yeast hybrids (e.g. Scer × Spar at ~10–15% divergence in non-coding
regions, ~1–2% in coding). Single-variant simulations don't give
salmon enough seed-level discrimination — we need many variants per
read for the haplotype-specific quant signal to emerge.
"""

from pathlib import Path
import random

GENOME_SRC = Path("/Users/black_einstein/Desktop/PhD_Coursework/SPRING2026/PHYLOGENETICS/Bioinformatics_Tools_/MycoNote_CLI/myconote-cli/scratch/e2e/genome.fa")
GFF3_SRC   = Path("/Users/black_einstein/Desktop/PhD_Coursework/SPRING2026/PHYLOGENETICS/Bioinformatics_Tools_/MycoNote_CLI/myconote-cli/scratch/e2e/genes_clean.gff3")
OUT_DIR    = Path("/tmp/ase_validation/fixture")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CONTIG     = "NW_003020038.1"
TX_CHOICES = [
    ("g000005.m1",  70),
    ("g000010.m1",  30),
    ("g000012.m1",  50),
]
SNVS_PER_TX = 20   # ≈ 1 per 60-80 bp on a 1.0-1.5 kb CDS
GENE_IDS = {m.rsplit(".", 1)[0] for m, _ in TX_CHOICES}

# ── 1. read genome ───────────────────────────────────────────────────────────
def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

genome = read_fasta(GENOME_SRC)
contig_seq = genome[CONTIG]

# ── 2. subset GFF3 ──────────────────────────────────────────────────────────
target_mrnas = {m for m, _ in TX_CHOICES}
kept_lines = []
cds_by_mrna = {}
seen = set()
for raw in GFF3_SRC.read_text().splitlines():
    if not raw.strip() or raw.startswith("#"):
        if raw.startswith("##gff-version") and not kept_lines: kept_lines.append(raw)
        continue
    cols = raw.split("\t")
    if len(cols) < 9 or cols[0] != CONTIG: continue
    attrs = dict(a.split("=", 1) for a in cols[8].split(";") if "=" in a)
    t = cols[2]
    if t == "gene"  and attrs.get("ID") in GENE_IDS:           kept_lines.append(raw)
    elif t == "mRNA"  and attrs.get("ID") in target_mrnas:     kept_lines.append(raw); seen.add(attrs["ID"])
    elif t == "exon"  and attrs.get("Parent") in target_mrnas: kept_lines.append(raw)
    elif t == "CDS"   and attrs.get("Parent") in target_mrnas:
        kept_lines.append(raw)
        cds_by_mrna.setdefault(attrs["Parent"], []).append((int(cols[3]), int(cols[4]), cols[6]))

assert seen == target_mrnas
(OUT_DIR / "ase.gff3").write_text("\n".join(kept_lines) + "\n")

# ── 3. one-contig genome FASTA ──────────────────────────────────────────────
sub_fa = OUT_DIR / "ase_genome.fa"
with sub_fa.open("w") as fh:
    fh.write(f">{CONTIG}\n")
    for i in range(0, len(contig_seq), 60): fh.write(contig_seq[i:i+60] + "\n")

# ── 4. dense phased VCF, ~SNVS_PER_TX per CDS ───────────────────────────────
rng = random.Random(42)
nts = list("ACGT")
vcf_rows = []
design_rows = ["transcript_id\tcds_start\tcds_end\tn_snvs\thap0_pct"]
all_variants_meta = []  # for ground-truth read simulation
for mrna_id, pct in TX_CHOICES:
    blocks = cds_by_mrna[mrna_id]
    assert len(blocks) == 1, mrna_id
    s, e, strand = blocks[0]
    # Pick SNV positions evenly spaced through the CDS, jittered.
    cds_len = e - s + 1
    step = cds_len // (SNVS_PER_TX + 1)
    used_pos = set()
    snv_positions = []
    for k in range(1, SNVS_PER_TX + 1):
        p = s + k * step + rng.randint(-step // 4, step // 4)
        # Stay safely inside the CDS so reads can still span the variants.
        p = max(s + 30, min(e - 30, p))
        if p in used_pos: continue
        used_pos.add(p)
        snv_positions.append(p)
    snv_positions.sort()
    for p in snv_positions:
        ref = contig_seq[p - 1].upper()
        if ref not in nts: continue
        alt = rng.choice([n for n in nts if n != ref])
        vcf_rows.append(f"{CONTIG}\t{p}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t0|1")
        all_variants_meta.append((mrna_id, s, e, p, ref, alt))
    design_rows.append(f"{mrna_id}\t{s}\t{e}\t{len(snv_positions)}\t{pct}")
    print(f"  {mrna_id}: CDS {s}-{e} (+), {len(snv_positions)} SNVs, hap0_pct={pct}")

vcf_header = [
    "##fileformat=VCFv4.2",
    "##INFO=<ID=.,Number=.,Type=String,Description=\"dummy\">",
    "##FILTER=<ID=PASS,Description=\"passed\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    f"##contig=<ID={CONTIG},length={len(contig_seq)}>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1",
]
(OUT_DIR / "phased.vcf").write_text("\n".join(vcf_header + vcf_rows) + "\n")
(OUT_DIR / "design.tsv").write_text("\n".join(design_rows) + "\n")

# Detail file used by simulate_reads_dense.py to apply the *same* variants
# we wrote into the VCF, in the same order.
import json
(OUT_DIR / "variants_meta.json").write_text(json.dumps(all_variants_meta))
print(f"wrote VCF with {len(vcf_rows)} variants")
