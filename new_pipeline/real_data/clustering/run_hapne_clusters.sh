#!/usr/bin/env bash
# HapNe-IBD on the cluster-selected subpopulations, both routes.
#
# The test: if the pooled superpopulation curves are flat/non-monotonic because
# the pool is not panmictic, then refitting on a homogeneous cluster drawn from
# the same segments should restore a real trajectory -- without changing the IBD
# data, the mask, the arms, or any HapNe setting.  Only the sample set moves.
set -euo pipefail
cd "$(dirname "$0")/.."

PY=/opt/miniconda3/envs/genetics_env/bin/python
for m in ibd snp; do
  for p in AFR EAS EUR SAS; do
    echo "=== ${p}_${m} ==="
    $PY hapne/run_pop.py "$p" \
        --labels "clustering/highcov_${m}.txt" \
        --name "${p}_${m}" \
        --out-root clustering/hapne
  done
done
