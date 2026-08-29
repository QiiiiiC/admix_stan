#!/usr/bin/env bash
# Stan inputs for the cluster-selected leaves, from the high_cov pipeline.
#
#   leaves : clustering/highcov_ibd_spectral.txt   (AFR/EAS/EUR/SAS = 99/94/101/87,
#            i.e. LWK / CDX / FIN / GIH, found blind by connectivity clustering)
#   IBD    : high_cov/ibd_all_masked.ibd.gz  -- already merged (0.6cM/1 discord)
#            and masked, one file for all 22 autosomes.  The chromosome comes
#            from the segment's own column 5, so a single merged file is fine.
#   SNP    : merged_pruned/vcf/*  -- LD-pruned genotypes.  Phase is irrelevant to
#            the covariance, so the old Beagle phasing does not matter here.
#            chr9 is EXCLUDED: merged_all.chr9_pruned.vcf.gz is a truncated gzip
#            (`gzip -t` fails; every other chromosome passes).  The [!9] class
#            drops it and keeps chr19 (whose first character is 1).  Only the SNP
#            side loses chr9 -- the IBD side reads high_cov, where chr9 is fine --
#            which is harmless because the two products are normalised
#            independently: the covariance is a per-site average and the IBD
#            fraction is divided by the genetic length of the chromosomes it
#            actually saw.
#   bins   : 37 uniform 0.5 cM bins from 2.0 to 20.5 cM (HapNe-style).
#
# The mask is passed even though ibd_all_masked.ibd.gz has already had those
# segments removed: re-running it drops nothing, but it is what subtracts the
# masked genetic length from the IBD normalizer, which still has to happen.
set -euo pipefail
cd "$(dirname "$0")/.."

PY=/opt/miniconda3/envs/genetics_env/bin/python

$PY generate_stan_data.py \
    --labels    clustering/highcov_ibd_spectral.txt \
    --pop-order AFR EAS EUR SAS \
    --ibd-method file \
    --ibd-glob  'high_cov/ibd_all_masked.ibd.gz' \
    --vcf-glob  'merged_pruned/vcf/merged_all.chr[!9]*_pruned.vcf.gz' \
    --genmap-dir high_cov/genmap \
    --mask-regions high_cov/bad_regions_grch38.bed \
    --bins-uniform 2.0 20.5 0.5 \
    --min-maf 0.05 \
    --out clustering/results_noadmix/stan_data/stan_4pop_clustered \
    2>&1 | tee clustering/make_stan_data.log
