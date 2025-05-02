import msprime
import numpy as np
import pandas as pd
import math

# create a list of migration matrices, whenever there's a event(split/merge), create a migration matrix for that particular event.
# the starting element is always an identity matrix.
# each migration matrix is of dimension k*k where k is the total number of all populations(instead of starting populations).
def mig_matrix(nodes, events, nodes_map):
    n = len(nodes)
    events_num = int(len(events)-1-len([e for e in events if e[-1]==0])/2)
    out = [np.diag([1.0 for i in range(n)]).tolist() for i in range(events_num + 1)]
    
    # cnt is for loop indices.
    cnt = 1
    for i in range(1,len(events)):
        # when there's a split backward in time, two elements are created in 'events', so skip one of them.
        if events[i][1] != events[i-1][1]:
            # this is a merge backward in time.
            if events[i][-1] == 1:
                dest = events[i][0]
                change = nodes[dest]['children']
                for j in change:
                    out[cnt][nodes_map[j]][nodes_map[j]]=0
                    out[cnt][nodes_map[j]][nodes_map[dest]]=1
                cnt += 1
            # this is a split backward in time
            else:
                change = nodes[events[i][0]]['children'][0]
                dest = nodes[change]['parents']
                out[cnt][nodes_map[change]][nodes_map[change]]=0
                for j in dest:
                    out[cnt][nodes_map[change]][nodes_map[j]]=nodes[j]['frac']
                cnt += 1
    return out

def all_ibd_segments(ts):
    n = ts.num_samples
    trees_iter = ts.trees()
    tree = next(trees_iter)
    last_mrca_m = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            last_mrca_m[i][j] = tree.mrca(i,j)
    last_left_m = np.zeros((n,n))
    segment_lengths_m = [[[]for x in range(n)]for y in range(n)]
    for tree in trees_iter:
        for i in range(n):
            for j in range(i,n):
                mrca = tree.mrca(i,j)
                last_mrca = last_mrca_m[i][j]
                if mrca!= last_mrca:
                    left = tree.interval[0]
                    last_left = last_left_m[i][j]
                    segment_lengths_m[i][j].append((left-last_left)/ts.sequence_length)
                    last_mrca_m[i][j] = mrca
                    last_left_m[i][j] = left
    for i in range(n):
        for j in range(i,n):
            segment_lengths_m[i][j].append((ts.sequence_length-last_left_m[i][j])/ts.sequence_length)
    return segment_lengths_m

def sim_ibd_admix(N,T,L,m,length,n,seed):
    #Set the follwing model as the default model.
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name="ADMIX", initial_size=N[1])
    demography.add_population(name="B", initial_size=N[2])
    demography.add_population(name="OUTER", initial_size=N[3])
    demography.add_population(name ='SUB_ANC', initial_size = N[4])
    demography.add_population(name = 'ANC', initial_size = N[5])
    demography.add_admixture(time=T[0], derived='ADMIX', ancestral=["A","B"],proportions = [0.25,0.75])
    demography.add_population_split(time=T[1], derived=["A", "B"], ancestral="SUB_ANC")
    demography.add_population_split(time = T[2], derived = ['SUB_ANC','OUTER'], ancestral = 'ANC')
    
    # bb is the number of bits simulated.
    bb = m*length

    # kk is measured in centiMorgans.
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": n, "ADMIX": n, "B": n, 'OUTER':n}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb,
        random_seed = seed
    )
    all = all_ibd_segments(ts)
    out = {'index':[],'number':[],'fraction':[],'u':[],'v':[],'group':[],'true_N':[],'true_T1':[],'true_T2':[],'true_T3':[]}

    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        #for within data, have (n*2)*(n*2-1)/2 pairs and for between data we have (n*2)*(n*2) pairs.
        #within A
        for j in range(n*2):
            for k in range(j+1,n*2):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,0]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]

       #between A and ADMIX
        for j in range(n*2):
            for k in range(n*2,n*4):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,1]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]

        #between A and B
        for j in range(n*2):
            for k in range(n*4,n*6):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,2]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]
        
        #between A and OUTER
        for j in range(n*2):
            for k in range(n*6,n*8):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,3]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]
        
        #within ADMIX and ADMIX
        for j in range(n*2,n*4):
            for k in range(j+1,n*4):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[1,1]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]

        #between ADMIX and B
        for j in range(n*2,n*4):
            for k in range(n*4,n*6):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[1,2]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]
        
        #between ADMIX and OUTER
        for j in range(n*2,n*4):
            for k in range(n*6,n*8):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[1,3]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]

        #within B and B
        for j in range(n*4,n*6):
            for k in range(j+1,n*6):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[2,2]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]
        
        #between B and OUTER
        for j in range(n*4,n*6):
            for k in range(n*6,n*8):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[2,3]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]

        #within OUTER and OUTER
        for j in range(n*6,n*8):
            for k in range(j+1,n*8):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[3,3]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
                out['true_T3'] += [T[2]]
                
    return out

def sim_snp_admix(N,T,m,length,n,seed):
    #Set the follwing model as the default model.
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name="ADMIX", initial_size=N[1])
    demography.add_population(name="B", initial_size=N[2])
    demography.add_population(name="OUTER", initial_size=N[3])
    demography.add_population(name ='SUB_ANC', initial_size = N[4])
    demography.add_population(name = 'ANC', initial_size = N[5])
    demography.add_admixture(time=T[0], derived='ADMIX', ancestral=["A","B"],proportions = [0.25,0.75])
    demography.add_population_split(time=T[1], derived=["A", "B"], ancestral="SUB_ANC")
    demography.add_population_split(time = T[2], derived = ['SUB_ANC','OUTER'], ancestral = 'ANC')
    
    # bb is the number of bits simulated.
    bb = m*length

    ts = msprime.sim_ancestry(
        samples={"A": n, "ADMIX": n, "B": n, 'OUTER':n}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb,
        random_seed = seed
    )
    return ts

def var_matrix(freq_array,b):
    n = len(freq_array)
    m = len(freq_array[0])
    bb = n//b
    output_mean_total = [np.zeros((m,m)) for i in range(bb)]
    output_var = np.zeros((m,m))
    for l in range(bb):
        for i in range(m):
            for j in range(i,m):
                total_sum = 0
                for k in range(l*10,(l*10+10)):
                    mean = sum(freq_array[k])/m
                    total_sum += (freq_array[k][i]-mean)*(freq_array[k][j]-mean)
                    
                output_mean_total[l][i][j] = total_sum/b
    output_mean = sum(output_mean_total)/len(output_mean_total)

    for i in range(m):
        for j in range(i,m):
            inc = 0
            for l in range(bb):
                inc += (output_mean_total[l][i][j] - output_mean[i][j])**2
            output_var[i][j] = math.sqrt(inc/(bb*(bb-1)))

    return (output_mean,output_var)
