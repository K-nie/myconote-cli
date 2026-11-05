#!/usr/bin/env python3
"""Q5: lost consensus votes.
On ref loci where funannotate strict-matches but MycoNote consensus does NOT,
determine whether any MycoNote raw predictor (augustus/snap/genemark) HAD a gene
model whose CDS structure strict-matches the reference. If so, the combiner had
the correct answer and discarded it.

All inputs lifted to reference (NC_) namespace. Strict = identical ordered CDS
interval set (same definition as compare_annotations.py / the scored benchmark).
"""
import sys, importlib.util
from collections import defaultdict
SCORER="/Users/black_einstein/Desktop/PhD_Coursework/SPRING2026/PHYLOGENETICS/Bioinformatics_Tools_/MycoNote_CLI/myconote-cli/docs/paper/htcondor_benchmark/scripts/compare_annotations.py"
spec=importlib.util.spec_from_file_location("ca",SCORER); ca=importlib.util.module_from_spec(spec); spec.loader.exec_module(ca)

g=sys.argv[1]
ref_path=sys.argv[2]
myc_rename=sys.argv[3]
fun_gff=sys.argv[4]
fun_rename=sys.argv[5]
cons=sys.argv[6]
aug=sys.argv[7]
snap=sys.argv[8] if sys.argv[8]!="none" else None
gm=sys.argv[9] if sys.argv[9]!="none" else None

mycmap=ca.load_rename_map(myc_rename)
funmap=ca.load_rename_map(fun_rename)

ref =ca.parse_gff3(ref_path, coding_only=True)
funp=ca.parse_gff3(fun_gff, seqid_map=funmap, coding_only=True)
consp=ca.parse_gff3(cons, seqid_map=mycmap, coding_only=True)
augp=ca.parse_gff3(aug, seqid_map=mycmap, coding_only=True)
snapp=ca.parse_gff3(snap, seqid_map=mycmap, coding_only=True) if snap else {}
gmp  =ca.parse_gff3(gm, seqid_map=mycmap, coding_only=True) if gm else {}

def cds_index(genes):
    """map (seqid,strand) -> set of CDS-interval-tuples -> representative gene"""
    idx=defaultdict(dict)
    for s,gl in genes.items():
        for gg in gl:
            cis=ca.cds_interval_set(gg)
            if cis: idx[(s,gg.strand)][cis]=gg
    return idx

ref_idx=cds_index(ref)
fun_idx=cds_index(funp)
cons_idx=cds_index(consp)
aug_idx=cds_index(augp)
snap_idx=cds_index(snapp)
gm_idx=cds_index(gmp)

def strict_has(idx, s, strand, cis): return cis in idx.get((s,strand),{})

# iterate every reference gene
total_ref=0
fun_ok=0                 # fun strict-matches ref
myc_miss_fun_ok=0        # fun ok but myconote consensus NOT strict
recoverable=0            # of those, a myc raw predictor HAD the exact structure
by_pred=defaultdict(int) # which predictor held the lost-but-correct model
examples=[]
for s,gl in ref.items():
    for rg in gl:
        cis=ca.cds_interval_set(rg)
        if not cis: continue
        total_ref+=1
        strand=rg.strand
        f=strict_has(fun_idx,s,strand,cis)
        c=strict_has(cons_idx,s,strand,cis)
        if f: fun_ok+=1
        if f and not c:
            myc_miss_fun_ok+=1
            held=[]
            if strict_has(aug_idx,s,strand,cis): held.append("augustus")
            if snap and strict_has(snap_idx,s,strand,cis): held.append("snap")
            if gm and strict_has(gm_idx,s,strand,cis): held.append("genemark")
            if held:
                recoverable+=1
                key="+".join(held)
                by_pred[key]+=1
                if len(examples)<6: examples.append((s,rg.start,rg.end,strand,key))

print(f"\n==== Q5 {g} ====")
print(f"reference coding genes                              : {total_ref}")
print(f"funannotate strict-correct                          : {fun_ok}")
print(f"  of which MycoNote consensus is NOT strict-correct : {myc_miss_fun_ok}")
print(f"  of those, a MycoNote RAW predictor HAD exact model: {recoverable}  ({100*recoverable/myc_miss_fun_ok:.1f}% of the gap-loci)" if myc_miss_fun_ok else "  n/a")
print(f"  --> combiner discarded a correct model this often : {recoverable} loci")
print(f"  predictor that held the discarded-correct model   :")
for k,v in sorted(by_pred.items(), key=lambda x:-x[1]):
    print(f"       {k:22s} {v}")
if examples:
    print("  examples (seqid,start,end,strand,held-by):")
    for e in examples: print(f"       {e}")
