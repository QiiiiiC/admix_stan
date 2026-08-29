#!/usr/bin/env bash
# End-to-end sweep for one leaf set: Stan inputs -> 21 topology fits -> ranking.
#
#   usage:  run_one.sh <tag> <POP1> <POP2> <POP3> [seeds...]
#   e.g.    run_one.sh gbr_ceu_ibs GBR CEU IBS
#           run_one.sh gbr_ibs_yri GBR IBS YRI 1 7 13
#
# Everything lands in 3pop/<tag>/.  Two tags are fully independent -- separate
# output dirs and separate CmdStan scratch -- so they can run concurrently.
#
# POP order must be identical in both steps: the saved matrix rows are indexed
# by pop_order and the topologies address those rows positionally.
set -euo pipefail
cd "$(dirname "$0")"

if [ "$#" -lt 4 ]; then
    echo "usage: $0 <tag> <POP1> <POP2> <POP3> [seeds...]" >&2
    exit 2
fi
TAG=$1; P1=$2; P2=$3; P3=$4; shift 4
SEEDS=${*:-"1 7 13 23 31 47"}
PY=/opt/miniconda3/envs/genetics_env/bin/python

if [ ! -f "$TAG/stan_data/stan_$TAG.npz" ]; then
    bash make_stan_data.sh "$TAG" "$P1" "$P2" "$P3"
else
    echo "[skip] $TAG/stan_data/stan_$TAG.npz already exists"
fi

$PY fit_topologies.py --tag "$TAG" --pops "$P1" "$P2" "$P3" --seeds $SEEDS \
    2>&1 | tee "$TAG/fit_all.log"
$PY compare_elbo.py --tag "$TAG" 2>&1 | tee "$TAG/comparison/ranking.txt"
