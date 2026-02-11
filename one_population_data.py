import msprime
import numpy as np
import pandas as pd
import math
from collections import defaultdict

# def all_ibd_segments(ts):
#     n = ts.num_samples
#     trees_iter = ts.trees()
#     tree = next(trees_iter)
#     last_mrca_m = np.zeros((n,n))
#     for i in range(n):
#         for j in range(i,n):
#             last_mrca_m[i][j] = tree.mrca(i,j)
#     last_left_m = np.zeros((n,n))
#     segment_lengths_m = [[[]for x in range(n)]for y in range(n)]
#     for tree in trees_iter:
#         for i in range(n):
#             for j in range(i,n):
#                 mrca = tree.mrca(i,j)
#                 last_mrca = last_mrca_m[i][j]
#                 if mrca!= last_mrca:
#                     left = tree.interval[0]
#                     last_left = last_left_m[i][j]
#                     segment_lengths_m[i][j].append((left-last_left)/ts.sequence_length)
#                     last_mrca_m[i][j] = mrca
#                     last_left_m[i][j] = left
#     for i in range(n):
#         for j in range(i,n):
#             segment_lengths_m[i][j].append((ts.sequence_length-last_left_m[i][j])/ts.sequence_length)
#     return segment_lengths_m

def all_ibd_segments(ts):
    """
    Manually extracts IBD segments by iterating trees, 
    but MERGES adjacent trees if the TMRCA for the pair is unchanged.
    """
    # Initialize a matrix to store current segment starts for every pair
    n = ts.num_samples
    current_mrca = np.zeros((n, n)) - 1
    current_start = np.zeros((n, n))
    
    # Store completed segments: segments[u][v] = [len1, len2, ...]
    segments = defaultdict(lambda: defaultdict(list))
    
    # Iterate continuously through the genome
    for tree in ts.trees():
        interval_end = tree.interval.right
        
        # This is efficient for small n, slow for large n
        # For every pair, check their MRCA time/node
        for i in range(n):
            for j in range(i + 1, n):
                mrca_node = tree.mrca(i, j)
                
                # If this is the first tree, initialize
                if current_mrca[i, j] == -1:
                    current_mrca[i, j] = mrca_node
                    current_start[i, j] = tree.interval.left
                
                # If MRCA changed, the segment ended. Record it.
                elif current_mrca[i, j] != mrca_node:
                    # Calculate length (fraction of genome or bp)
                    # Here we store fraction to match your 'l' definition
                    seg_len = (tree.interval.left - current_start[i, j]) / ts.sequence_length
                    segments[i][j].append(seg_len)
                    
                    # Start new segment
                    current_mrca[i, j] = mrca_node
                    current_start[i, j] = tree.interval.left
                    
    # Flush the final segments at the end of the chromosome
    for i in range(n):
        for j in range(i + 1, n):
            seg_len = (ts.sequence_length - current_start[i, j]) / ts.sequence_length
            segments[i][j].append(seg_len)
            
    return segments

def simulate(N,L,m,length,n,s):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
     # bb is the number of bits simulated.
    bb = m*length

    # kk is measured in centiMorgans.
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": n}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb,
        random_seed = s,
        model = 'smc'
    )
    all = all_ibd_segments(ts)
    out = {'index':[],'fraction':[],'u':[],'v':[],'group':[]}
    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        #for within data, have (n*2)*(n*2-1)/2 pairs and for between data we have (n*2)*(n*2) pairs.
        #within A
        for j in range(n*2):
            for k in range(j+1,n*2):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,0]]
    return out

def simulate2(N,T,L,m,length,n):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name = 'B', initial_size = N[0])
    demography.add_population(name = 'anc', initial_size = N[0])
    demography.add_population_split(time = T, derived = ['A','B'], ancestral = 'anc')

     # bb is the number of bits simulated.
    bb = m*length

    # kk is measured in centiMorgans.
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": n,'B':n}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb,
    )
    all = all_ibd_segments(ts)
    out = {'index':[],'fraction':[],'u':[],'v':[],'group':[]}
    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        #for within data, have (n*2)*(n*2-1)/2 pairs and for between data we have (n*2)*(n*2) pairs.
        #within A
        for j in range(n*2):
            for k in range(n*2,n*4):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,1]]
    return out