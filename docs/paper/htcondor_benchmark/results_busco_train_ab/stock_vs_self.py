import json, sys, statistics
LOCAL = sys.argv[1]
G = "sce cal ylp ani ncr cne klac kpha fgra afum umay ccin spom tmel".split()
def m(a, g):
    try:
        d = json.load(open(f"{LOCAL}/{a}/{g}/metrics.json"))
        return d["gene_level"]["loose"]["F1"], d["gene_level"]["strict"]["F1"], d["predicted_summary"]["n_genes"]
    except Exception:
        return None
hdr = ["g", "stk_L", "slf_L", "dL", "stk_S", "slf_S", "dS", "stk_ng", "slf_ng"]
print("".join(f"{h:>8}" for h in hdr))
dls, dss = [], []
for g in G:
    s, f = m("stock", g), m("self", g)
    if s and f:
        dl, ds = f[0] - s[0], f[1] - s[1]
        dls.append(dl); dss.append(ds)
        print(f"{g:>8}{s[0]:8.4f}{f[0]:8.4f}{dl:+8.4f}{s[1]:8.4f}{f[1]:8.4f}{ds:+8.4f}{s[2]:8d}{f[2]:8d}")
print(f"\nmean self-stock delta: loose {statistics.mean(dls):+.4f}  strict {statistics.mean(dss):+.4f}  (n={len(dls)})")
print(f"self LOSES on loose {sum(1 for d in dls if d<0)}/14, on strict {sum(1 for d in dss if d<0)}/14")
