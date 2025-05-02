# This file contains the code that simulates SNPs data and IBD data for given population structure parameters.
# The population structure in this case contains three populations and one outer group. The most ancestral time and the outer group's population size will be known in inference.
# For the SNPs data, follow the procedure in 'Notes on SNPs simulation'.
from methods import sim_ibd_admix
import json 

seed = 42
length = 1e6
n= 10
m = 2
N = [5000,5000,5000,5000,10000,10000]
L =  [0,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1,1.5,2.5,5,20,200] 
T = [30,100,1000]

data = sim_ibd_admix(N,T,L,m,length,n,seed)
with open(f'ibd_admix_data.json','w') as json_file:
    json.dump(data,json_file)