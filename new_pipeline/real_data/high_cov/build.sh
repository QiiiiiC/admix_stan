#!/bin/bash
# ---------------------------------------------------------------------------
# high_cov IBD pipeline (GRCh38, consortium SHAPEIT2-duoHMM phased)
#   1. drop relatives  -> 2504 unrelated panel, restricted to AFR/EAS/EUR/SAS (2157)
#   2. global MAF>=0.05 across ALL retained samples (one marker set for every pop)
#   3. uniform thinning to ~700k genome-wide (array density; NOT LD pruning)
#   4. hap-ibd  min-seed=1.5
#   5. merge-ibd-segments 0.6 / 1
#   6. concatenate to a single ibd_all.ibd.gz for downstream use
# ---------------------------------------------------------------------------
set -e
export PATH=/opt/miniconda3/envs/genetics_env/bin:$PATH
H=/Users/qichen/github/admix_stan/new_pipeline/real_data/high_cov
JAR=/Users/qichen/github/admix_stan/new_pipeline/real_data/tools/merge-ibd-segments.jar
P=CCDG_14151_B01_GRM_WGS_2020-08-05
S=$H/samples/unrelated_4pop.txt
mkdir -p $H/{array_vcf,ibd,merged,logs}

# ---- pass A: count common variants genome-wide to set the thinning stride ----
if [ ! -f $H/logs/maf_counts.txt ]; then
  : > $H/logs/maf_counts.txt
  for c in $(seq 1 22); do
    n=$(bcftools view -S $S --force-samples -Ou $H/vcf/${P}_chr${c}.filtered.shapeit2-duohmm-phased.vcf.gz \
        | bcftools view -m2 -M2 -v snps -q 0.05:minor -H | wc -l)
    echo "$c $n" >> $H/logs/maf_counts.txt
    echo "  counted chr$c: $n"
  done
fi
TOT=$(awk '{s+=$2} END{print s}' $H/logs/maf_counts.txt)
STRIDE=$(( (TOT + 350000) / 700000 )); [ $STRIDE -lt 1 ] && STRIDE=1
echo "MAF>=0.05 total=$TOT  stride=$STRIDE  -> ~$((TOT/STRIDE)) markers"
echo $STRIDE > $H/logs/stride.txt

# ---- pass B: subset + MAF + thin -> hap-ibd -> merge ----
for c in $(seq 1 22); do
  V=$H/array_vcf/chr${c}.vcf.gz
  M=$H/genmap/plink.chr${c}.GRCh38.map
  if [ ! -f "$V" ]; then
    bcftools view -S $S --force-samples -Ou $H/vcf/${P}_chr${c}.filtered.shapeit2-duohmm-phased.vcf.gz \
      | bcftools view -m2 -M2 -v snps -q 0.05:minor \
      | awk -v s=$STRIDE '/^#/{print;next} {n++; if(n%s==0) print}' | bgzip > "$V"
    bcftools index -t "$V"
  fi
  [ -f "$H/ibd/chr${c}.ibd.gz" ] || hap-ibd gt=$V map=$M out=$H/ibd/chr${c} \
      min-seed=1.5 nthreads=8 > $H/logs/hapibd_chr${c}.log 2>&1
  if [ ! -f "$H/merged/chr${c}.ibd.gz" ]; then
    gunzip -c $H/ibd/chr${c}.ibd.gz \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,"1.0",$8}' \
      | java -jar $JAR $V $M 0.6 1 \
      | awk 'BEGIN{OFS="\t"}{h2=($2==0?1:$2); h4=($4==0?1:$4); print $1,h2,$3,h4,$5,$6,$7,$9}' \
      | gzip > $H/merged/chr${c}.ibd.gz
  fi
  echo "chr$c markers=$(gunzip -c $V|grep -vc '^#') raw=$(gunzip -c $H/ibd/chr${c}.ibd.gz|wc -l) merged=$(gunzip -c $H/merged/chr${c}.ibd.gz|wc -l)"
done

# ---- single combined IBD file for downstream (our framework + hapne subsetting) ----
gunzip -c $H/merged/chr*.ibd.gz | gzip > $H/ibd_all.ibd.gz
echo "ibd_all.ibd.gz segments: $(gunzip -c $H/ibd_all.ibd.gz | wc -l)"
echo BUILD_DONE
