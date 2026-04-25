#!/usr/bin/env python3
"""
Build hap0 + hap1 personalized transcriptomes manually (ground truth),
then simulate paired-end reads from each at the configured ratio so we
have a sample whose true hap0:hap1 mix is known.

We do this independently from src/ase/personalize.rs — that's the unit
we're validating. The whole point of the fixture is that we know the
right answer because we constructed the reads ourselves.

Output: fixture/sample1_R1.fastq.gz + fixture/sample1_R2.fastq.gz
"""

import gzip
import random
import subprocess
from pathlib import Path

OUT_DIR = Path("/tmp/ase_validation/fixture")
DESIGN  = OUT_DIR / "design.tsv"
CDS_FA  = OUT_DIR / "ase_cds.fa"

# Simulation parameters chosen so the binomial test has enough power per
# transcript while staying small enough that wgsim finishes in seconds.
TOTAL_READS_PER_TX = 4000   # combined hap0+hap1 reads per transcript
READ_LEN = 100
INSERT   = 200

# ── 1. read the reference CDS FASTA ─────────────────────────────────────────
def read_fasta(p):
    seqs = {}
    name = None
    buf = []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs

cds = read_fasta(CDS_FA)
print(f"loaded {len(cds)} reference transcripts")

# ── 2. read design + apply each variant to hap1 CDS ─────────────────────────
design = []
for line in DESIGN.read_text().splitlines()[1:]:
    f = line.split("\t")
    design.append({
        "mrna_id":   f[0],
        "cds_start": int(f[2]),
        "cds_end":   int(f[3]),
        "snv_pos":   int(f[4]),
        "snv_ref":   f[5],
        "snv_alt":   f[6],
        "hap0_pct":  int(f[7]),
    })

# For each design row, swap REF→ALT in hap1's copy of the transcript.
# All chosen CDSs are plus-strand single-block, so CDS-frame position =
# snv_pos − cds_start (0-based).
hap0_seqs = {tx: cds[tx] for tx in cds}
hap1_seqs = {tx: list(cds[tx]) for tx in cds}
for d in design:
    tx = d["mrna_id"]
    cds_idx = d["snv_pos"] - d["cds_start"]   # 0-based offset into CDS
    cur = hap1_seqs[tx][cds_idx].upper()
    assert cur == d["snv_ref"], f"{tx}: expected {d['snv_ref']} at CDS pos {cds_idx}, got {cur}"
    hap1_seqs[tx][cds_idx] = d["snv_alt"]
hap1_seqs = {tx: "".join(s) for tx, s in hap1_seqs.items()}

# Sanity: hap0 != hap1 at exactly 1 position per transcript.
for d in design:
    tx = d["mrna_id"]
    diffs = sum(a != b for a, b in zip(hap0_seqs[tx], hap1_seqs[tx]))
    assert diffs == 1, f"{tx}: expected 1 diff, got {diffs}"
print("hap0 / hap1 differ at exactly 1 position per transcript ✓")

# Write each haplotype's transcripts as a separate FASTA so wgsim can
# simulate from one haplotype at a time.
def write_fa(path, seqs):
    with path.open("w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")

hap0_fa = OUT_DIR / "truth_hap0.fa"
hap1_fa = OUT_DIR / "truth_hap1.fa"
write_fa(hap0_fa, hap0_seqs)
write_fa(hap1_fa, hap1_seqs)
print(f"wrote {hap0_fa}, {hap1_fa}")

# ── 3. simulate reads via wgsim, separately per (transcript, haplotype) ─────
# wgsim simulates from a FASTA. To enforce per-transcript ratios, simulate
# per single-transcript FASTA. wgsim parameters:
#   -1 / -2 read length
#   -d / -s mean / sd insert size
#   -e error rate (we keep low so true SNV is recoverable)
#   -r mutation rate 0 (we want only the SNV we injected)
#   -R indel fraction 0
#   -N total read pair count
random.seed(7)
all_r1 = OUT_DIR / "sample1_R1.fastq"
all_r2 = OUT_DIR / "sample1_R2.fastq"
for f in (all_r1, all_r2):
    if f.exists():
        f.unlink()

read_counter = 0
for d in design:
    tx = d["mrna_id"]
    n_total = TOTAL_READS_PER_TX
    n_hap0  = int(round(n_total * d["hap0_pct"] / 100))
    n_hap1  = n_total - n_hap0
    print(f"  {tx}: hap0_reads={n_hap0}, hap1_reads={n_hap1}")
    for hap_name, n_reads, src_seq in [
        ("hap0", n_hap0, hap0_seqs[tx]),
        ("hap1", n_hap1, hap1_seqs[tx]),
    ]:
        if n_reads == 0:
            continue
        single_fa = OUT_DIR / f"sim_{tx}_{hap_name}.fa"
        single_fa.write_text(f">{tx}\n{src_seq}\n")
        r1 = OUT_DIR / f"sim_{tx}_{hap_name}_R1.fq"
        r2 = OUT_DIR / f"sim_{tx}_{hap_name}_R2.fq"
        cmd = [
            "wgsim",
            "-1", str(READ_LEN), "-2", str(READ_LEN),
            "-d", str(INSERT), "-s", "20",
            "-e", "0.001", "-r", "0", "-R", "0", "-X", "0",
            "-N", str(n_reads),
            "-S", str(read_counter + 1),
            str(single_fa), str(r1), str(r2),
        ]
        read_counter += 1
        subprocess.run(cmd, check=True, capture_output=True)
        # Append to combined FASTQ
        with all_r1.open("a") as out_fh, r1.open() as in_fh:
            out_fh.write(in_fh.read())
        with all_r2.open("a") as out_fh, r2.open() as in_fh:
            out_fh.write(in_fh.read())
        single_fa.unlink()
        r1.unlink()
        r2.unlink()

# gzip combined FASTQs (fastp/salmon are happy either way; gzipped saves
# disk and exercises the gz path).
for fq in (all_r1, all_r2):
    subprocess.run(["gzip", "-f", str(fq)], check=True)
print(f"wrote {all_r1}.gz, {all_r2}.gz")

# ── 4. write sample sheet ───────────────────────────────────────────────────
sheet = OUT_DIR / "samples.tsv"
sheet.write_text(
    "sample_id\tfastq_1\tfastq_2\tcondition\n"
    f"sample1\t{all_r1}.gz\t{all_r2}.gz\tF1_hybrid\n"
)
print(f"wrote {sheet}")
