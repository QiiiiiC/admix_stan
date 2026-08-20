"""Apply bad_regions_grch38.bed to ibd_all.ibd.gz -> ibd_all_masked.ibd.gz.

A segment is dropped if it OVERLAPS any masked window, not merely if its
midpoint falls inside one.  The artifact is a spurious haplotype match seeded in
the bad region; hap-ibd then extends it outward, so a segment that starts inside
and runs out is contaminated over its whole reported length, and its length --
which is the quantity every downstream model fits -- cannot be trusted.  This
also matches how generate_stan_data.py masks, so the two stay consistent.
"""
import gzip, os, sys
from collections import defaultdict

B = os.path.dirname(os.path.abspath(__file__))
bed = sys.argv[1] if len(sys.argv) > 1 else f"{B}/bad_regions_grch38.bed"
src = f"{B}/ibd_all.ibd.gz"
dst = f"{B}/ibd_all_masked.ibd.gz"

mask = defaultdict(list)
for line in open(bed):
    if line.strip() and not line.startswith("#"):
        c = line.split()
        mask[c[0]].append((int(c[1]), int(c[2])))
print(f"[mask] {sum(len(v) for v in mask.values())} regions from {os.path.basename(bed)}")

n_in = n_out = 0
per_chr_in = defaultdict(int); per_chr_out = defaultdict(int)
with gzip.open(src, "rt") as fh, gzip.open(dst, "wt") as out:
    for line in fh:
        line = line.rstrip("\n")
        c = line.split("\t")
        n_in += 1
        ch, s, e = c[4], int(c[5]), int(c[6])
        per_chr_in[ch] += 1
        if any(s < hi and e > lo for lo, hi in mask.get(ch, [])):
            continue
        out.write(line + "\n")
        n_out += 1
        per_chr_out[ch] += 1
print(f"[mask] {n_in:,} -> {n_out:,} segments ({100*(n_in-n_out)/n_in:.1f}% removed)")
print(f"[save] {dst}")
print(f"\n{'chr':>6} {'before':>10} {'after':>10} {'kept':>7}")
for ch in sorted(per_chr_in, key=lambda x: int(x.replace("chr", ""))):
    print(f"{ch:>6} {per_chr_in[ch]:>10,} {per_chr_out[ch]:>10,} "
          f"{100*per_chr_out[ch]/per_chr_in[ch]:>6.1f}%")
