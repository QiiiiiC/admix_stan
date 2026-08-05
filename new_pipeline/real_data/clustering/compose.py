"""How close is each selected cluster to a reported 1000G subpopulation?

Diagnostic only -- the clustering never sees the panel, and nothing downstream
depends on this file.  It exists to answer one question: did "most connected
cluster >= min size" land on a real group, or on an arbitrary slice?

Three overlap numbers, because "% overlap" is ambiguous for two sets of
different size and each answers a different question:

  purity  |C n S| / |C|   what fraction of the CLUSTER is that subpopulation.
                          Low purity = the cluster pooled several groups, which
                          is the substructure we were trying to remove.
  recall  |C n S| / |S|   what fraction of the SUBPOPULATION was captured.
                          Low recall = the cluster is a fragment of a real group.
  Jaccard |C n S| / |C u S|   the symmetric one; 1.0 only if C == S exactly.

Purity alone is not enough: a cluster of 5 GWD is 100% pure and 4% of GWD.
"""
import os
import json
from collections import Counter, defaultdict

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROUTES = [("SNP kmeans", "highcov_snp.txt"),
          ("IBD louvain", "highcov_ibd_louvain.txt"),
          ("IBD spectral", "highcov_ibd_spectral.txt"),
          ("IBD dbscan", "highcov_ibd_dbscan.txt")]
POPS = ["AFR", "EAS", "EUR", "SAS"]

panel = {}
with open(os.path.join(B, "integrated_call_samples_v3.20130502.ALL.panel")) as f:
    next(f)
    for line in f:
        c = line.split()
        if len(c) >= 3:
            panel[c[0]] = c[1]

# subpopulation sizes WITHIN our unrelated 4-pop panel, not 1000G at large
full = defaultdict(Counter)
for line in open(os.path.join(B, "high_cov", "samples", "labels_4pop.txt")):
    s, p = line.split()[:2]
    if s in panel:
        full[p][panel[s]] += 1

lines = [__doc__.strip(), "",
         f"{'route':13s} {'pop':4s} {'n':>4s} {'subpop':>7s} {'purity':>7s} "
         f"{'recall':>7s} {'Jacc':>6s}  composition of the selected cluster",
         "-" * 118]
rows = []
for rname, fn in ROUTES:
    path = os.path.join(B, "clustering", fn)
    if not os.path.exists(path):
        continue
    sel = defaultdict(Counter)
    for line in open(path):
        s, p = line.split()[:2]
        if s in panel:
            sel[p][panel[s]] += 1
    for p in POPS:
        cnt = sel[p]
        n = sum(cnt.values())
        if not n:
            continue
        top, hit = cnt.most_common(1)[0]
        nS = full[p][top]
        purity = hit / n
        recall = hit / nS if nS else 0.0
        jacc = hit / (n + nS - hit) if (n + nS - hit) else 0.0
        comp = ", ".join(f"{k} {v}/{full[p][k]}" for k, v in cnt.most_common(5))
        lines.append(f"{rname:13s} {p:4s} {n:>4d} {top:>7s} {purity:6.0%} "
                     f"{recall:6.0%} {jacc:6.2f}  {comp}")
        rows.append(dict(route=rname, pop=p, n=n, subpop=top,
                         purity=round(purity, 4), recall=round(recall, 4),
                         jaccard=round(jacc, 4),
                         composition={k: [v, full[p][k]] for k, v in cnt.most_common()}))
    lines.append("")

text = "\n".join(lines)
print(text)
out_txt = os.path.join(B, "clustering", "cluster_overlap.txt")
out_json = os.path.join(B, "clustering", "cluster_overlap.json")
open(out_txt, "w").write(text + "\n")
json.dump(rows, open(out_json, "w"), indent=2)
print(f"[save] {out_txt}\n[save] {out_json}")
