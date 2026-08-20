#!/usr/bin/env bash
# HapNe-IBD on every 1000G subpopulation and superpopulation, from the rebuilt
# high_cov/ibd_all_masked.ibd.gz.
#
# The mask is now bad_regions_grch38.bed built by
# high_cov/reference_masks/build_final_mask.py: HapNe's 39 chromosome arms,
# plus IBDNe's extreme-locus rule (6x median per-cM coverage, no 50 cM fill),
# plus our cross-population enrichment criterion.  96.59 cM.
#
# Failures are expected and tolerated: a population with too few segments in too
# few arms cannot be fit, and that is itself a result -- see the summary table.
set -u
cd "$(dirname "$0")"
PY=/opt/miniconda3/envs/genetics_env/bin/python
mkdir -p logs

SUB="GWD YRI IBS TSI CHS JPT CHB GIH STU ITU ESN LWK KHV FIN CEU ACB PJL CDX GBR BEB MSL ASW"
SUP="AFR EAS EUR SAS"

for p in $SUB; do
  echo "=== $p (subpopulation) ==="
  $PY run_pop.py "$p" > logs/$p.log 2>&1 \
    && tail -2 logs/$p.log \
    || echo "  FAILED (see logs/$p.log): $(grep -iEm1 'error|exception' logs/$p.log | cut -c1-90)"
done

for p in $SUP; do
  echo "=== $p (superpopulation) ==="
  $PY run_pop.py "$p" --super > logs/$p.log 2>&1 \
    && tail -2 logs/$p.log \
    || echo "  FAILED (see logs/$p.log): $(grep -iEm1 'error|exception' logs/$p.log | cut -c1-90)"
done
echo "ALL_DONE"
