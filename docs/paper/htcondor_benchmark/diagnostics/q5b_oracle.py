#!/usr/bin/env python3
"""Q5b: oracle-combiner headroom. For ALL reference coding genes, how many does
the shipped consensus get strict-correct, vs how many a *perfect* combiner could
get by always choosing the raw predictor (augustus/snap/genemark) that has the
exact-correct CDS structure when any one does. Quantifies strict-F1 left on the
table by the current weighted-vote combiner, independent of funannotate.
"""
import sys, importlib.util
from collections import defaultdict
SCORER="/Users/black_einstein/Desktop/PhD_Coursework/SPRING2026/PHYLOGENETICS/Bioinformatics_Tools_/MycoNote_CLI/myconote-cli/docs/paper/htcondor_benchmark/scripts/compare_annotations.py"
spec=importlib.util.spec_from_file_location("ca",SCORER); ca=importlib.util.module_from_spec(spec); spec.loader.exec_module(ca)
g,ref_path,myc_rename,cons,aug,snap,gm=sys.argv[1:8]
mycmap=ca.load_rename_map(myc_rename)
ref=ca.parse_gff3(ref_path,coding_only=True)
def idx(p):
    d=defaultdict(set)
    if p=="none": return d
    gg=ca.parse_gff3(p,seqid_map=mycmap,coding_only=True)
    for s,gl in gg.items():
        for x in gl:
            cis=ca.cds_interval_set(x)
            if cis: d[(s,x.strand)].add(cis)
    return d
ci=idx(cons); ai=idx(aug); si=idx(snap); gi=idx(gm)
tot=cons_ok=oracle_ok=0
only_raw=0  # raw has it, consensus doesn't
for s,gl in ref.items():
    for rg in gl:
        cis=ca.cds_interval_set(rg)
        if not cis: continue
        tot+=1
        key=(s,rg.strand)
        c = cis in ci[key]
        r = (cis in ai[key]) or (cis in si[key]) or (cis in gi[key])
        if c: cons_ok+=1
        if c or r: oracle_ok+=1
        if r and not c: only_raw+=1
print(f"{g:5s} ref={tot:6d}  consensus_strict={cons_ok:6d} ({100*cons_ok/tot:5.1f}%)  "
      f"oracle_strict={oracle_ok:6d} ({100*oracle_ok/tot:5.1f}%)  "
      f"recoverable_by_better_combiner={only_raw:5d} (+{100*only_raw/tot:4.1f} pts sens)")
