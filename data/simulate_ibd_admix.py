# This file contains the code that simulates SNPs data and IBD data for given population structure parameters.
# The population structure in this case contains three populations and one outer group. The most ancestral time and the outer group's population size will be known in inference.
# For the SNPs data, follow the procedure in 'Notes on SNPs simulation'.
from methods import sim_ibd_admix
import json 

seed = 42
length = 1e6
n= 20
m = 0.3
N = [5000,5000,5000,5000,5000,5000]
L = [0.5,0.6,0.7,0.8,0.9,1,1.5,2,10,300]
T = [50,500,1000]
alpha = 0.75

data = sim_ibd_admix(N,T,L,m,length,n,alpha,seed)
with open(f'ibd_admix_data.json','w') as json_file:
    json.dump(data,json_file)