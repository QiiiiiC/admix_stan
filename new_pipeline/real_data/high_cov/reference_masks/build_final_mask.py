"""Production IBD mask = three layers, in order of how well-referenced they are.

  1. outside-arm   HapNe's 39 GRCh38 chromosome arms (Fournier et al. 2023,
                   Nat Commun 14:7107, src/hapne/files/regions_grch38.txt).
                   A POSITIVE definition of analysable genome; the complement
                   (centromeres, acrocentric short arms, telomere caps) is
                   excluded by construction.  These regions carry 33 markers/Mb
                   against 257/Mb elsewhere -- 8x sparser -- which is the
                   mechanism: with no markers there is no evidence of
                   discordance, so hap-ibd extends straight through and emits
                   spurious long segments.

  2. extreme-ibd   IBDNe's rule (Browning & Browning 2015, AJHG 97:404-418):
                   loci carrying an extreme number of IBD segments.  Threshold
                   6x the median per-locus coverage, measured per cM on
                   within-subpopulation pairs.  Chosen because the flagged
                   REGION SET is identical from 5x to 10x -- only its edges
                   move -- so the answer is insensitive across a factor of two.
                   IBDNe's second clause (also drop surviving fragments < 50 cM
                   between anchors) is NOT applied: measured here it buys 1.7%
                   on the exposure ratio for 197 cM, and only 0.5% of the map
                   lies in fragments that short.  These loci have NORMAL marker
                   density (273/Mb), so they are a different artifact from (1) --
                   false haplotype matching from duplication or mismapping.

  3. crosspop      Cross-population enrichment >= 4x the genome-wide cross/within
                   ratio (high_cov/derive_mask.py criterion (c)).  Ours, not
                   published, but 9 of its 10 regions coincide with Li et al.
                   2014 (PLoS Genet 10:e1004144) Table 3 or fall outside the
                   arms, so it is reproducing known artifact regions by an
                   independent route.  Catches what (2) cannot: a locus can be
                   enriched in cross-continental pairs -- which share almost no
                   true IBD -- without its total coverage being extreme.

Layers (2) and (3) are complementary by construction: (2) is a within-population
coverage test, (3) is a between-population ratio test.  8p23 (chr8:12-13 Mb) and
chr1:133-135 Mb are found only by (3); chr2:2-6 Mb only by (2).
"""
import os, argparse
import numpy as np

B = os.path.dirname(os.path.abspath(__file__))
HC = os.path.dirname(B)

ap = argparse.ArgumentParser()
ap.add_argument("--arm-ibdne", default=f"{B}/bad_regions_arm_ibdne_nofill_grch38.bed")
ap.add_argument("--crosspop", default=f"{B}/crosspop_enrich_grch38.bed")
ap.add_argument("--genmap-dir", default=f"{HC}/genmap")
ap.add_argument("--out", default=f"{HC}/bad_regions_grch38.bed")
a = ap.parse_args()

GM = {}
def cm(ch, s, e):
    c = ch[3:] if ch.startswith("chr") else ch
    if c not in GM:
        d = np.loadtxt(f"{a.genmap_dir}/plink.chr{c}.GRCh38.map", usecols=(2, 3))
        GM[c] = d[np.argsort(d[:, 1])]
    d = GM[c]
    return float(np.interp(e, d[:, 1], d[:, 0]) - np.interp(s, d[:, 1], d[:, 0]))


def read(path, default_tag):
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        c = line.split()
        out.append((c[0], int(c[1]), int(c[2]),
                    c[3] if len(c) > 3 and not c[3][0].isdigit() else default_tag))
    return out


regs = read(a.arm_ibdne, "arm-ibdne") + read(a.crosspop, "crosspop")
by = {}
for ch, s, e, tag in regs:
    by.setdefault(tag, 0)
    by[tag] += cm(ch, s, e)
print("input layers (cM, before merging overlaps):")
for t, v in sorted(by.items()):
    print(f"   {t:<14} {v:7.2f} cM")

merged = []
for ch, s, e, tag in sorted(regs, key=lambda r: (int(r[0][3:]), r[1])):
    if merged and merged[-1][0] == ch and s <= merged[-1][2]:
        merged[-1][2] = max(merged[-1][2], e)
        if tag not in merged[-1][3].split("+"):
            merged[-1][3] += "+" + tag
    else:
        merged.append([ch, s, e, tag])

with open(a.out, "w") as f:
    f.write("# Production IBD mask, GRCh38.  Built by reference_masks/build_final_mask.py\n")
    f.write("# outside-arm : complement of HapNe's 39 chromosome arms "
            "(Fournier et al. 2023, Nat Commun 14:7107)\n")
    f.write("# extreme-ibd : IBDNe extreme-locus rule (Browning & Browning 2015, "
            "AJHG 97:404-418), 6x median per-cM coverage, no 50 cM fill\n")
    f.write("# crosspop    : cross-population enrichment >= 4x "
            "(high_cov/derive_mask.py criterion c)\n")
    for ch, s, e, tag in merged:
        f.write(f"{ch}\t{s}\t{e}\t{tag}\n")

tot = sum(cm(c, s, e) for c, s, e, _ in merged)
print(f"\n[save] {a.out}: {len(merged)} regions, {tot:.2f} cM "
      f"({100*tot/3472:.1f}% of the 3472 cM inside the arms)")
print("\nregions that are NOT simply outside-arm:")
for ch, s, e, tag in merged:
    if "outside-arm" not in tag:
        print(f"   {ch}:{s/1e6:>7.1f}-{e/1e6:<7.1f} Mb {cm(ch,s,e):6.2f} cM  {tag}")
