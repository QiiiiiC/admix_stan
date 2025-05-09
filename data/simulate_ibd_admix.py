# This file contains the code that simulates SNPs data and IBD data for given population structure parameters.
# The population structure in this case contains three populations and one outer group. The most ancestral time and the outer group's population size will be known in inference.
# For the SNPs data, follow the procedure in 'Notes on SNPs simulation'.
from methods import sim_ibd_admix
import json 

seed = 42
length = 1e7
n= 20
m = 2
N = [3000,3000,3000,3000,7000,7000]
L = [0.7,0.75,0.8,0.85,0.9,0.95,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,10]
T = [30,150,400]

data = sim_ibd_admix(N,T,L,m,length,n,seed)
with open(f'ibd_admix_data.json','w') as json_file:
    json.dump(data,json_file)