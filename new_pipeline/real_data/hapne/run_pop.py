"""HapNe-IBD for one population from the high_cov pipeline.

Usage:  python run_pop.py FIN                     (1000G subpopulation)
        python run_pop.py EAS --super             (1000G superpopulation)
        python run_pop.py AFR --labels ../clustering/highcov_ibd.txt \
                              --name AFR_ibd --out-root ../clustering/hapne

Reads high_cov/ibd_all_masked.ibd.gz, keeps within-population pairs, splits into
the 39 GRCh38 chromosome arms, and fits HapNe-IBD.  Any label file of the form
"<sample> <pop>" works, so the same code path serves the reported 1000G groups
and the clusters found by cluster_populations.py.
"""
import os, sys, gzip, argparse
from collections import defaultdict

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HC = os.path.join(B, "high_cov")
REG = os.path.join(B, "tools", "HapNe", "src", "hapne", "files", "regions_grch38masked.txt")

ap = argparse.ArgumentParser()
ap.add_argument("pop")
ap.add_argument("--super", action="store_true", help="pop is a superpopulation")
ap.add_argument("--labels", default=None,
                help="'<sample> <pop>' file; overrides the built-in 1000G labels")
ap.add_argument("--name", default=None, help="output subfolder (default: pop)")
ap.add_argument("--out-root", default=os.path.join(B, "hapne"))
ap.add_argument("--mask", default=None, help="BED-like file of regions to drop")
a = ap.parse_args()

if a.labels:
    labfile = a.labels if os.path.isabs(a.labels) else os.path.abspath(a.labels)
else:
    labfile = os.path.join(HC, "samples",
                           "labels_4pop.txt" if a.super else "labels_subpop.txt")
members = {c[0] for c in (l.split() for l in open(labfile)) if len(c) > 1 and c[1] == a.pop}
if not members:
    sys.exit(f"no samples for {a.pop} in {labfile}")
name = a.name or a.pop
print(f"{name}: {len(members)} individuals  ({os.path.basename(labfile)})")

mask = defaultdict(list)
if a.mask:
    for line in open(a.mask):
        if line.strip() and not line.startswith("#"):
            c = line.split(); mask[c[0]].append((float(c[1]), float(c[2])))

regions = []
with open(REG) as f:
    next(f)
    for line in f:
        c = line.split("\t")
        regions.append(("chr" + c[0], int(c[1]), int(c[2]), c[3]))
by_chr = defaultdict(list)
for ch, lo, hi, nm in regions:
    by_chr[ch].append((lo, hi, nm))

out = os.path.join(a.out_root, name)
os.makedirs(f"{out}/IBD", exist_ok=True)
H = {n: gzip.open(f"{out}/IBD/{n}.ibd.gz", "wt") for _, _, _, n in regions}
kept = defaultdict(int); n_pair = n_mask = n_span = 0
with gzip.open(f"{HC}/ibd_all_masked.ibd.gz", "rt") as fh:
    for line in fh:
        line = line.rstrip("\n")
        c = line.split("\t")
        if c[0] not in members or c[2] not in members:
            continue
        n_pair += 1
        ch, s, e = c[4], int(c[5]), int(c[6])
        if any(s < hi and e > lo for lo, hi in mask.get(ch, [])):
            n_mask += 1; continue
        hit = False
        for lo, hi, nm in by_chr.get(ch, []):
            if s >= lo and e <= hi:
                H[nm].write(line + "\n"); kept[nm] += 1; hit = True; break
        if not hit:
            n_span += 1
for h in H.values():
    h.close()
empty = [n for _, _, _, n in regions if kept[n] == 0]
for n in empty:
    os.remove(f"{out}/IBD/{n}.ibd.gz")
print(f"  within-pop segs={n_pair} masked={n_mask} out-of-arm={n_span} "
      f"written={sum(kept.values())} regions={len(regions)-len(empty)}")

with open(f"{out}/config.ini", "w") as f:
    f.write(f"""[CONFIG]
output_folder={out}
nb_samples={len(members)}
population_name={name}
ibd_files={out}/IBD
column_cm_length=8
genome_build=grch38masked
""")

sys.path.insert(0, os.path.join(B, "tools", "HapNe", "src"))
from configparser import ConfigParser
from hapne.ibd import build_hist
from hapne import hapne_ibd
import logging
logging.basicConfig(level=logging.WARNING)
cfg = ConfigParser(); cfg.read(f"{out}/config.ini")
build_hist(cfg)
hapne_ibd(cfg)
import pandas as pd
r = pd.read_csv(f"{out}/HapNe/hapne.csv")
print(f"  Ne: " + "  ".join(f"gen{int(r.TIME[i])}={r['Q0.5'][i]:,.0f}"
                            for i in (0, 9, 24, 49, 99) if i < len(r)))
