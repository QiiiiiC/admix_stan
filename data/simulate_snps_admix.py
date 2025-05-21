import subprocess
from methods import sim_snp_admix, var_matrix
import msprime
import tskit
import numpy as np
import pandas as pd
import json

seed = 42
length = 1e7
n= 20
m = 2
N = [3000,3000,3000,3000,3000,3000]
T = [30,150,800]

ts = sim_snp_admix(N,T,m,length,n,seed)
mut_ts = msprime.sim_mutations(ts,discrete_genome=False, rate=1e-8, random_seed=1234)


anc_t = T[-1] 
ancestral_mutations = []

for mut in mut_ts.mutations():
    node = mut.node  
    time = mut_ts.node(node).time  
    
    if time > anc_t:  
        ancestral_mutations.append(mut.site)
filtered_ts = mut_ts.delete_sites([mut.site for mut in mut_ts.mutations() if mut.id not in ancestral_mutations])
with open("snp_ts.vcf", "w") as vcf_file:
    filtered_ts.write_vcf(vcf_file)

vcf_file = "snp_ts.vcf"
output_prefix = "plink_output"
pruned_prefix = "pruned"


def run_cmd(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

# Convert VCF to PLINK binary
run_cmd(f"plink --vcf {vcf_file} --make-bed --out {output_prefix} --double-id")

# LD pruning
run_cmd(f"plink --bfile {output_prefix} --indep-pairwise 100 5 0.2 --out {pruned_prefix}")

# Extract pruned SNPs
run_cmd(f"plink --bfile {output_prefix} --extract {pruned_prefix}.prune.in --make-bed --out {pruned_prefix}_final")

print("LD pruning complete.")

with open(f"{pruned_prefix}.prune.in", "r") as file:
    pruned_snps_list = [line.strip() for line in file]
pruned_snps_list = [int(i) for i in pruned_snps_list]


pop1_samples = np.arange(2*n)  
pop2_samples = np.arange(2*n, 4*n)  
pop3_samples = np.arange(4*n,6*n)
pop4_samples = np.arange(6*n,8*n)


freq_array = []
freq_ = []
final_pruned = []
last_position = 0
for variant in filtered_ts.variants():
    if variant.index in pruned_snps_list:
          
        genotypes = variant.genotypes  
        freq_pop1 = np.sum(genotypes[pop1_samples]) / len(pop1_samples)
        freq_pop2 = np.sum(genotypes[pop2_samples]) / len(pop2_samples)
        freq_pop3 = np.sum(genotypes[pop3_samples]) / len(pop3_samples)
        freq_pop4 = np.sum(genotypes[pop4_samples]) / len(pop4_samples)
        freq_pop = np.sum(genotypes)/(8*n)
        if freq_pop > 0.1:
            final_pruned.append(variant.index)
            freq_.append(freq_pop)
            freq_array.append([freq_pop1, freq_pop2,freq_pop3,freq_pop4])

freq_array = np.array(freq_array)

adjusted_ratio = sum([i*(1-i) for i in freq_])/len(freq_)
print(adjusted_ratio)
# Using the block method
# b is the number of SNPs per block
output_mean, output_var = var_matrix(freq_array, b = 20)
# m is the number of populations
m = len(output_mean)

cols = ['group','d_mean','d_var','adjusted_factor']
data = pd.DataFrame(columns = cols)
for i in range(m):
    for j in range(i,m):
        new_row = [[i,j],output_mean[i][j], output_var[i][j],adjusted_ratio]
        data.loc[len(data)] = new_row
json_snp = data.to_dict(orient='list')
with open('json_snp.json', 'w') as f:
    json.dump(json_snp, f)

