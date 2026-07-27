#!/bin/bash
# Reproduces the HapNe paper's IBD pipeline exactly:
#   per-population sample subset -> hap-ibd DEFAULT params -> merge-ibd-segments defaults
set -e
export PATH=/opt/miniconda3/envs/genetics_env/bin:$PATH
BASE=/Users/qichen/github/admix_stan/new_pipeline/real_data
D=$BASE/hapne_fin_exact
JAR=$BASE/tools/merge-ibd-segments.jar
for c in 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22; do
  V=$D/vcf/chr$c.vcf.gz
  M=$BASE/merged_all/genmap/plink.chr$c.GRCh37.map
  [ -f "$V" ] || bcftools view -S $D/fin.samples --force-samples -Oz -o $V \
      $BASE/merged_all/vcf/merged_all.phased.chr$c.vcf.gz
  # hap-ibd with DEFAULT parameters (min-seed=2.0, min-output=2.0, min-extend=1.0,
  # min-markers=100, min-mac=2, max-gap=1000) -- exactly as the paper states
  [ -f "$D/ibd/chr$c.ibd.gz" ] || hap-ibd gt=$V map=$M out=$D/ibd/chr$c nthreads=4 > /dev/null 2>&1
  # merge-ibd-segments with DEFAULT parameters (gap 0.6 cM, 1 discordant homozygote)
  if [ ! -f "$D/merged/chr$c.ibd.gz" ]; then
    gunzip -c $D/ibd/chr$c.ibd.gz \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,"1.0",$8}' \
      | java -jar $JAR $V $M 0.6 1 \
      | awk 'BEGIN{OFS="\t"}{h2=($2==0?1:$2); h4=($4==0?1:$4); print $1,h2,$3,h4,$5,$6,$7,$9}' \
      | gzip > $D/merged/chr$c.ibd.gz
  fi
  echo "chr$c raw=$(gunzip -c $D/ibd/chr$c.ibd.gz|wc -l) merged=$(gunzip -c $D/merged/chr$c.ibd.gz|wc -l)"
done
echo PIPELINE_DONE
