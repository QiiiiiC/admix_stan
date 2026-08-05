#!/usr/bin/env bash
# Substructure clustering on the high_cov pipeline output.
#
# Why this exists: HapNe on the pooled 1000G superpopulations gives curves that
# no demography produces -- EAS falls only 6x between generation 4 and 100
# (essentially flat), AFR 4x -- while FIN alone falls 76x and reproduces the
# published curve.  A pooled superpopulation is not one panmictic population,
# and both HapNe and the coalescent leaf in the Stan model assume it is.  So:
# find the homogeneous core of each superpopulation and refit on that.
#
#   IBD  : high_cov/ibd_all_masked.ibd.gz  (merged, masked, min-seed 1.5)
#   SNP  : merged_pruned/vcf/*             (LD-pruned genotypes, phase irrelevant)
#
# Selection is --select cohesive: keep clusters of >= --min-cluster-size, then
# take the most internally cohesive one.  See cluster_populations.py's docstring.
set -euo pipefail
cd "$(dirname "$0")/.."

PY=/opt/miniconda3/envs/genetics_env/bin/python
OUT=clustering/highcov

$PY cluster_populations.py \
    --labels   high_cov/samples/labels_4pop.txt \
    --pop-order AFR EAS EUR SAS \
    --ibd-glob 'high_cov/ibd_all_masked.ibd.gz' \
    --vcf-glob 'merged_pruned/vcf/*_pruned.vcf.gz' \
    --panel    integrated_call_samples_v3.20130502.ALL.panel \
    --ibd-min-cm 2.0 \
    --ibd-cluster spectral \
    --select cohesive \
    --min-cluster-size 50 \
    --n-pc 4 --kmax 8 --knn 20 \
    --out "$OUT" \
    2>&1 | tee clustering/run.log
