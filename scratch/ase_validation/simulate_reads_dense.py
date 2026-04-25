#!/usr/bin/env python3
"""Simulate reads for the dense-variant fixture."""

import gzip
import json
import random
import subprocess
from pathlib import Path

OUT_DIR = Path("/tmp/ase_validation/fixture")
DESIGN  = OUT_DIR / "design.tsv"
CDS_FA  = OUT_DIR / "ase_cds.fa"
META    = OUT_DIR / "variants_meta.json"

TOTAL_READS_PER_TX = 4000
READ_LEN = 100
INSERT   = 200

def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name is not None: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None: seqs[name] = "".join(buf)
    return seqs

cds = read_fasta(CDS_FA)
variants = json.loads(META.read_text())

# Apply REF→ALT to hap1 at each variant position (CDS-frame; all + strand
# single-block so cds_idx = genome_pos − cds_start).
hap0_seqs = {tx: cds[tx] for tx in cds}
hap1_seqs = {tx: list(cds[tx]) for tx in cds}
for (tx, cds_s, cds_e, gpos, ref, alt) in variants:
    cds_idx = gpos - cds_s
    cur = hap1_seqs[tx][cds_idx].upper()
    assert cur == ref, f"{tx}@{gpos}: expected {ref}, got {cur}"
    hap1_seqs[tx][cds_idx] = alt
hap1_seqs = {tx: "".join(s) for tx, s in hap1_seqs.items()}

# Verify hap0 / hap1 differ at the right number of sites per transcript.
counts = {}
for tx in cds:
    counts[tx] = sum(a != b for a, b in zip(hap0_seqs[tx], hap1_seqs[tx]))
print("hap0/hap1 diffs per tx:", counts)

# Per-tx target read counts.
design = []
for line in DESIGN.read_text().splitlines()[1:]:
    f = line.split("\t")
    design.append({"mrna": f[0], "hap0_pct": int(f[4])})

# Simulate via wgsim.
all_r1 = OUT_DIR / "sample1_R1.fastq"
all_r2 = OUT_DIR / "sample1_R2.fastq"
for f in (all_r1, all_r2, Path(str(all_r1) + ".gz"), Path(str(all_r2) + ".gz")):
    if f.exists(): f.unlink()

read_counter = 0
for d in design:
    tx = d["mrna"]
    n_total = TOTAL_READS_PER_TX
    n_hap0 = int(round(n_total * d["hap0_pct"] / 100))
    n_hap1 = n_total - n_hap0
    print(f"  {tx}: hap0_reads={n_hap0}, hap1_reads={n_hap1}")
    for hap_name, n_reads, src in [
        ("hap0", n_hap0, hap0_seqs[tx]),
        ("hap1", n_hap1, hap1_seqs[tx]),
    ]:
        if n_reads == 0: continue
        single = OUT_DIR / f"sim_{tx}_{hap_name}.fa"
        single.write_text(f">{tx}\n{src}\n")
        r1 = OUT_DIR / f"sim_{tx}_{hap_name}_R1.fq"
        r2 = OUT_DIR / f"sim_{tx}_{hap_name}_R2.fq"
        cmd = ["wgsim",
               "-1", str(READ_LEN), "-2", str(READ_LEN),
               "-d", str(INSERT), "-s", "20",
               "-e", "0.001", "-r", "0", "-R", "0", "-X", "0",
               "-N", str(n_reads),
               "-S", str(read_counter + 1),
               str(single), str(r1), str(r2)]
        read_counter += 1
        subprocess.run(cmd, check=True, capture_output=True)
        with all_r1.open("a") as o, r1.open() as i: o.write(i.read())
        with all_r2.open("a") as o, r2.open() as i: o.write(i.read())
        single.unlink(); r1.unlink(); r2.unlink()

for fq in (all_r1, all_r2):
    subprocess.run(["gzip", "-f", str(fq)], check=True)

(OUT_DIR / "samples.tsv").write_text(
    "sample_id\tfastq_r1\tfastq_r2\tcondition\n"
    f"sample1\t{all_r1}.gz\t{all_r2}.gz\tF1_hybrid\n"
)
print(f"wrote {OUT_DIR / 'samples.tsv'}")
