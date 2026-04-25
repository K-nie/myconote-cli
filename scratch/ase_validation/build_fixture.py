#!/usr/bin/env python3
"""
Build a controlled ASE validation fixture.

Pipeline:
1. Subset the Candida tropicalis genome fixture to one contig.
2. Subset the GFF3 to 3 plus-strand single-block CDSs on that contig.
3. Inject 3 phased heterozygous SNVs, one inside each CDS, picking a
   reference base that can be changed without creating a stop codon
   in the CDS reading frame (arbitrary — we only care about counts).
4. Write the phased VCF.
5. Emit a design TSV recording the true hap0:hap1 read ratio we will
   simulate for each transcript.
"""

from pathlib import Path
import random

GENOME_SRC = Path("/Users/black_einstein/Desktop/PhD_Coursework/SPRING2026/PHYLOGENETICS/Bioinformatics_Tools_/MycoNote_CLI/myconote-cli/scratch/e2e/genome.fa")
GFF3_SRC   = Path("/Users/black_einstein/Desktop/PhD_Coursework/SPRING2026/PHYLOGENETICS/Bioinformatics_Tools_/MycoNote_CLI/myconote-cli/scratch/e2e/genes_clean.gff3")
OUT_DIR    = Path("/tmp/ase_validation/fixture")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CONTIG    = "NW_003020038.1"
TX_CHOICES = [
    # (mrna_id, ratio_hap0_percent)
    ("g000005.m1",  70),   # hap0-biased  — should read as strong ASE toward hap0
    ("g000010.m1",  30),   # hap1-biased  — strong ASE toward hap1
    ("g000012.m1",  50),   # balanced — no ASE
]
# Gene IDs are mRNA IDs without the trailing `.m1`.
GENE_IDS = {m.rsplit(".", 1)[0] for m, _ in TX_CHOICES}

# ── 1. read genome ───────────────────────────────────────────────────────────
def read_fasta(p):
    seqs = {}
    name = None
    buf = []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs

genome = read_fasta(GENOME_SRC)
contig_seq = genome[CONTIG]
print(f"contig {CONTIG}: {len(contig_seq):,} bp")

# ── 2. subset GFF3 to target gene set ───────────────────────────────────────
target_mrnas = {m for m, _ in TX_CHOICES}
kept_lines = []
cds_by_mrna = {}
seen_mrnas = set()
for raw in GFF3_SRC.read_text().splitlines():
    if raw.startswith("#") or not raw.strip():
        if raw.startswith("##gff-version") and not kept_lines:
            kept_lines.append(raw)
        continue
    cols = raw.split("\t")
    if len(cols) < 9 or cols[0] != CONTIG:
        continue
    attrs = dict(a.split("=", 1) for a in cols[8].split(";") if "=" in a)
    if cols[2] == "gene":
        gid = attrs.get("ID", "")
        if gid in GENE_IDS:
            kept_lines.append(raw)
    elif cols[2] == "mRNA":
        mid = attrs.get("ID", "")
        if mid in target_mrnas:
            seen_mrnas.add(mid)
            kept_lines.append(raw)
    elif cols[2] == "exon":
        pid = attrs.get("Parent", "")
        if pid in target_mrnas:
            kept_lines.append(raw)
    elif cols[2] == "CDS":
        pid = attrs.get("Parent", "")
        if pid in target_mrnas:
            kept_lines.append(raw)
            cds_by_mrna.setdefault(pid, []).append((int(cols[3]), int(cols[4]), cols[6]))

assert seen_mrnas == target_mrnas, f"expected {target_mrnas}, got {seen_mrnas}"

sub_gff = OUT_DIR / "ase.gff3"
sub_gff.write_text("\n".join(kept_lines) + "\n")
print(f"wrote {sub_gff}: {len(kept_lines)} lines")

# ── 3. subset genome to one contig so salmon index stays small ──────────────
sub_fa = OUT_DIR / "ase_genome.fa"
with sub_fa.open("w") as fh:
    fh.write(f">{CONTIG}\n")
    for i in range(0, len(contig_seq), 60):
        fh.write(contig_seq[i:i+60] + "\n")
print(f"wrote {sub_fa}")

# ── 4. build phased VCF with 1 SNV per CDS ──────────────────────────────────
# Rule: pick a position near the middle of each CDS on the + strand, swap
# REF base for one different nt. Emit GT=0|1 (het, phased).
random.seed(42)
nts = list("ACGT")
vcf_rows = []
design_rows = ["transcript_id\tmrna_id\tcds_start\tcds_end\tsnv_pos\tsnv_ref\tsnv_alt\thap0_pct"]
for mrna_id, pct in TX_CHOICES:
    cds_blocks = cds_by_mrna[mrna_id]
    assert len(cds_blocks) == 1, f"{mrna_id} has {len(cds_blocks)} CDS blocks; expected 1"
    s, e, strand = cds_blocks[0]
    snv_pos = (s + e) // 2  # 1-based
    ref = contig_seq[snv_pos - 1].upper()
    if ref not in nts:
        ref = "A"
    alt = next(n for n in nts if n != ref)
    vcf_rows.append(f"{CONTIG}\t{snv_pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t0|1")
    design_rows.append(f"{mrna_id}\t{mrna_id}\t{s}\t{e}\t{snv_pos}\t{ref}\t{alt}\t{pct}")
    print(f"  mRNA={mrna_id}: CDS {s}-{e} +strand, SNV {ref}>{alt} @ {snv_pos}, hap0_pct={pct}")

vcf_header = [
    "##fileformat=VCFv4.2",
    "##INFO=<ID=.,Number=.,Type=String,Description=\"dummy\">",
    "##FILTER=<ID=PASS,Description=\"passed\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    f"##contig=<ID={CONTIG},length={len(contig_seq)}>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1",
]
vcf_path = OUT_DIR / "phased.vcf"
vcf_path.write_text("\n".join(vcf_header + vcf_rows) + "\n")
print(f"wrote {vcf_path}")

design_path = OUT_DIR / "design.tsv"
design_path.write_text("\n".join(design_rows) + "\n")
print(f"wrote {design_path}")
