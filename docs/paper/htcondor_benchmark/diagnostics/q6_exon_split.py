#!/usr/bin/env python3
"""Q6: gene-level loose+strict F1 split by SINGLE-exon vs MULTI-exon genes.
Partition each side (pred, ref) by its own coding-exon (CDS) count:
single = 1 CDS, multi = >=2 CDS. Then compute sens/spec/F1 within each class.
Strict = identical ordered CDS structure (Eilbeck). Loose = same-strand span overlap.
"""
import sys, importlib.util
from collections import defaultdict
SCORER="/Users/black_einstein/Desktop/PhD_Coursework/SPRING2026/PHYLOGENETICS/Bioinformatics_Tools_/MycoNote_CLI/myconote-cli/docs/paper/htcondor_benchmark/scripts/compare_annotations.py"
spec=importlib.util.spec_from_file_location("ca",SCORER); ca=importlib.util.module_from_spec(spec); spec.loader.exec_module(ca)

pred_path,ref_path,rename,label=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
seqmap=ca.load_rename_map(rename) if rename not in("-","none") else None
pred=ca.parse_gff3(pred_path,seqid_map=seqmap,coding_only=True)
ref =ca.parse_gff3(ref_path,coding_only=True)

def ncds(g): return len(ca.get_cds_features(g))

def metrics(pred,ref,mode):
    # returns dict class-> (tp_sens, n_ref, tp_spec, n_pred)
    res={'single':[0,0,0,0],'multi':[0,0,0,0]}
    # specificity: iterate pred genes, find a ref match of ANY class for loose,
    # or exact-structure match for strict; classify by PRED's own exon count
    ref_by=defaultdict(list)
    for s,gl in ref.items():
        for g in gl: ref_by[s].append(g)
    for s in ref_by: ref_by[s].sort(key=lambda g:g.start)
    for s,gl in pred.items():
        for pg in gl:
            cls='single' if ncds(pg)==1 else 'multi'
            res[cls][3]+=1
            matched=False
            for rg in ref_by.get(s,[]):
                if pg.strand!=rg.strand: continue
                if rg.end<pg.start: continue
                if rg.start>pg.end: break
                if mode=='loose':
                    if ca.overlap_length(pg.start,pg.end,rg.start,rg.end)>0: matched=True;break
                else:
                    if ca.genes_match_strict(pg,rg): matched=True;break
            if matched: res[cls][2]+=1
    # sensitivity: iterate ref genes, classify by REF's own exon count
    pred_by=defaultdict(list)
    for s,gl in pred.items():
        for g in gl: pred_by[s].append(g)
    for s in pred_by: pred_by[s].sort(key=lambda g:g.start)
    for s,gl in ref.items():
        for rg in gl:
            cls='single' if ncds(rg)==1 else 'multi'
            res[cls][1]+=1
            matched=False
            for pg in pred_by.get(s,[]):
                if pg.strand!=rg.strand: continue
                if pg.end<rg.start: continue
                if pg.start>rg.end: break
                if mode=='loose':
                    if ca.overlap_length(pg.start,pg.end,rg.start,rg.end)>0: matched=True;break
                else:
                    if ca.genes_match_strict(pg,rg): matched=True;break
            if matched: res[cls][0]+=1
    return res

print(f"\n==== {label} ====")
for mode in ('loose','strict'):
    r=metrics(pred,ref,mode)
    print(f" [{mode}]")
    for cls in('single','multi'):
        tpS,nref,tpP,npred=r[cls]
        sens=tpS/nref if nref else 0; spec=tpP/npred if npred else 0
        f1=2*sens*spec/(sens+spec) if (sens+spec) else 0
        print(f"   {cls:7s}  F1={f1:.4f}  sens={sens:.4f}(ref {tpS}/{nref})  spec={spec:.4f}(pred {tpP}/{npred})")
