#!/usr/bin/env bash
# Stan inputs for one three-population run.
#
#   usage:  make_stan_data.sh <tag> <POP1> <POP2> <POP3>
#   e.g.    make_stan_data.sh gbr_ceu_ibs GBR CEU IBS
#
# Writes 3pop/<tag>/stan_data/stan_<tag>.{npz,json} plus the labels file it used.
# One npz serves all 21 topologies for that leaf set: the saved matrices are
# indexed by --pop-order, and a topology only decides which EVENTS connect those
# rows, never their order.  So build once per leaf set, fit 21 times.
#
#   IBD  : high_cov/ibd_all_masked.ibd.gz  (merged 0.6cM/1 discord, then masked
#          with the referenced 96.59 cM mask)
#   SNP  : merged_pruned/vcf/* -- chr9 excluded, its gzip is truncated (`gzip -t`
#          fails).  Only the SNP side loses it; the IBD side reads high_cov where
#          chr9 is fine.  Harmless: the two products are normalised independently.
#   bins : 37 uniform 0.5 cM bins, 2.0 -> 20.5 cM (HapNe-style).
#
# --mask-regions is passed even though ibd_all_masked has already had those
# segments removed: the second pass drops nothing, but it is what subtracts the
# masked genetic length from the IBD normaliser, which still has to happen.
set -euo pipefail
cd "$(dirname "$0")/.."

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <tag> <POP1> <POP2> <POP3>" >&2
    exit 2
fi
TAG=$1; P1=$2; P2=$3; P3=$4
PY=/opt/miniconda3/envs/genetics_env/bin/python
DIR=3pop/$TAG
mkdir -p "$DIR/stan_data"

LAB=$DIR/labels_$TAG.txt
awk -v a="$P1" -v b="$P2" -v c="$P3" '$2==a||$2==b||$2==c' \
    high_cov/samples/labels_subpop.txt > "$LAB"
echo "[labels] $LAB"
awk '{print $2}' "$LAB" | sort | uniq -c

$PY generate_stan_data.py \
    --labels    "$LAB" \
    --pop-order "$P1" "$P2" "$P3" \
    --ibd-method file \
    --ibd-glob  'high_cov/ibd_all_masked.ibd.gz' \
    --vcf-glob  'merged_pruned/vcf/merged_all.chr[!9]*_pruned.vcf.gz' \
    --genmap-dir high_cov/genmap \
    --mask-regions high_cov/bad_regions_grch38.bed \
    --bins-uniform 2.0 20.5 0.5 \
    --min-maf 0.05 \
    --out "$DIR/stan_data/stan_$TAG" \
    2>&1 | tee "$DIR/make_stan_data.log"
