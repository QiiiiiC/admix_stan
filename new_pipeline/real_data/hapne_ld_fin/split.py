"""Replicates hapne.convert.tools.split_convert_vcf_in_parallel, but sources the
per-chromosome VCF (HapNe assumes one all-chromosome file)."""
import os, sys, pandas as pd
BASE="/Users/qichen/github/admix_stan/new_pipeline/real_data"
OUT=f"{BASE}/hapne_ld_fin/DATA/GENOTYPES"; os.makedirs(OUT, exist_ok=True)
reg=pd.read_csv(f"{BASE}/tools/HapNe/src/hapne/src_regions.txt" if False else
                f"{BASE}/tools/HapNe/src/hapne/files/regions_grch37noc9.txt", sep="\t")
for _,r in reg.iterrows():
    c=r["CHR"]; name=r["NAME"]
    if os.path.exists(f"{OUT}/{name}.bed"): continue
    vcf=f"{BASE}/merged_all/vcf/merged_all.phased.chr{c}.vcf.gz"
    cmd=(f"plink --vcf {vcf} --make-bed --keep {BASE}/hapne_ld_fin/fin.keep "
         f"--out {OUT}/{name} --cm-map {BASE}/hapne_ld_fin/map/genetic_map_chr@.txt "
         f"--chr {c} --from-bp {r['FROM_BP']} --to-bp {r['TO_BP']} --const-fid "
         f"--threads 2 --memory 4096 --maf 0.249 --snps-only --geno 0.8 > /dev/null 2>&1")
    os.system(cmd)
    ok=os.path.exists(f"{OUT}/{name}.bed")
    print(f"{'ok ' if ok else 'FAIL'} {name}", flush=True)
