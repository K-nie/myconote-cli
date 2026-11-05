#!/usr/bin/env python3
"""Q3 error taxonomy: classify MycoNote (or any predictor) gene-level
disagreements vs a reference into missed / extra / correct(exact) /
correct-locus-wrong-boundaries / merged / split / tangled.

Reuses compare_annotations.py parsing so the gene set, coding-only filter and
seqid lifting are IDENTICAL to the scored benchmark.
"""
import sys, importlib.util
from collections import defaultdict

SCORER = "/Users/black_einstein/Desktop/PhD_Coursework/SPRING2026/PHYLOGENETICS/Bioinformatics_Tools_/MycoNote_CLI/myconote-cli/docs/paper/htcondor_benchmark/scripts/compare_annotations.py"
spec = importlib.util.spec_from_file_location("ca", SCORER)
ca = importlib.util.module_from_spec(spec); spec.loader.exec_module(ca)

pred_path, ref_path, rename = sys.argv[1], sys.argv[2], sys.argv[3]
seqmap = ca.load_rename_map(rename) if rename not in ("-", "none") else None
pred = ca.parse_gff3(pred_path, seqid_map=seqmap, coding_only=True)
ref  = ca.parse_gff3(ref_path, coding_only=True)

# assign stable indices
pred_nodes, ref_nodes = [], []
for s, gl in pred.items():
    for g in gl: pred_nodes.append((s, g))
for s, gl in ref.items():
    for g in gl: ref_nodes.append((s, g))

# build adjacency: pred i <-> ref j if same strand and gene-span overlap (loose)
adj_p = defaultdict(set); adj_r = defaultdict(set)
ref_by_seq = defaultdict(list)
for j,(s,g) in enumerate(ref_nodes): ref_by_seq[s].append(j)
for s in ref_by_seq: ref_by_seq[s].sort(key=lambda j: ref_nodes[j][1].start)

for i,(s,pg) in enumerate(pred_nodes):
    for j in ref_by_seq.get(s, []):
        rg = ref_nodes[j][1]
        if pg.strand != rg.strand: continue
        if rg.end < pg.start: continue
        if rg.start > pg.end: break
        if ca.overlap_length(pg.start,pg.end,rg.start,rg.end) > 0:
            adj_p[i].add(j); adj_r[j].add(i)

# connected components over bipartite graph
seen_p=set(); seen_r=set()
cats = defaultdict(int)          # category -> count of pred genes (or ref genes for missed)
locus_examples = defaultdict(list)
# per-category also track ref-gene counts for the "missed" side
def comp(i0, ispred):
    stack=[('p',i0)] if ispred else [('r',i0)]
    ps=set(); rs=set()
    while stack:
        t,x=stack.pop()
        if t=='p':
            if x in ps: continue
            ps.add(x)
            for j in adj_p[x]:
                if j not in rs: stack.append(('r',j))
        else:
            if x in rs: continue
            rs.add(x)
            for k in adj_r[x]:
                if k not in ps: stack.append(('p',k))
    return ps, rs

n_pred=len(pred_nodes); n_ref=len(ref_nodes)
missed=extra=exact=wrongb=merged=split=tangled=0
merged_ref=split_pred=tangled_pred=tangled_ref=0

for i in range(n_pred):
    if i in seen_p: continue
    if not adj_p[i]:
        extra+=1; seen_p.add(i); continue
    ps, rs = comp(i, True)
    seen_p|=ps; seen_r|=rs
    np_, nr_ = len(ps), len(rs)
    if np_==1 and nr_==1:
        p=next(iter(ps)); r=next(iter(rs))
        if ca.genes_match_strict(pred_nodes[p][1], ref_nodes[r][1]):
            exact+=1
        else:
            wrongb+=1
    elif np_==1 and nr_>=2:
        merged+=1; merged_ref+=nr_
    elif np_>=2 and nr_==1:
        split+=1; split_pred+=np_
    else:
        tangled+=1; tangled_pred+=np_; tangled_ref+=nr_

# missed = ref genes with no pred overlap
for j in range(n_ref):
    if not adj_r[j]:
        missed+=1

g=pred_path.split('/')
print(f"\n==== {sys.argv[4] if len(sys.argv)>4 else pred_path} ====")
print(f"pred genes={n_pred}  ref genes={n_ref}")
print(f"{'category':32s}{'count':>8s}  {'%pred':>7s} {'%ref':>7s}  note")
def row(name,c,note=""):
    print(f"{name:32s}{c:8d}  {100*c/n_pred:6.1f}% {100*c/n_ref:6.1f}%  {note}")
row("exact (correct CDS structure)", exact)
row("correct-locus-wrong-boundaries", wrongb)
row("merged (1 myc spans >=2 ref)", merged, f"covers {merged_ref} ref genes")
row("split (>=2 myc in 1 ref)", split, f"from {split_pred} myc genes")
row("tangled (>=2 x >=2)", tangled, f"{tangled_pred} myc / {tangled_ref} ref")
row("extra (myc, no ref overlap)", extra)
print(f"{'missed (ref, no myc overlap)':32s}{missed:8d}  {'':>7s} {100*missed/n_ref:6.1f}%")
# accounting check
acct_pred = exact+wrongb+merged+split_pred+tangled_pred+extra
print(f"pred accounted: {acct_pred} (vs {n_pred})")
# among 1:1 matched loci, fraction exact vs wrong-boundary
oneone = exact+wrongb
if oneone:
    print(f"\nAmong 1:1-matched loci (n={oneone}): exact={exact} ({100*exact/oneone:.1f}%)  wrong-boundaries={wrongb} ({100*wrongb/oneone:.1f}%)")
