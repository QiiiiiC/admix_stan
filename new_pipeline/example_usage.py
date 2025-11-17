from methods import UnifiedDemography, simulate_msprime, simulate_snp_pruning, simulate_snp_data
from pathlib import Path


############################## Parameter block ##############################

ud = UnifiedDemography()
ud.add_population("R", 3000) # node times (gens before present)
ud.add_population("E", 1500)
ud.add_population("A", 0)
ud.add_population("Admix", 0)
ud.add_population("C", 0)
ud.add_population('D', 0)

"""
Population splits:
tree edges with edge-specific Ne, 
from child to parent with specified effective population size.
"""
ud.add_edge('E', 'R', Ne=1.0e4)
ud.add_edge('D', 'R', Ne=1.0e4)
ud.add_edge('A', 'E', Ne=1.0e4)
ud.add_edge('C', 'E', Ne=1.0e4)

"""
Admixture event: 
Admix splits goes to A and C at time 600 gens ago,
with proportion 0.25 to A and 0.75 to C.
"""
ud.add_admixture(time=600, child_present="Admix", p1="A", p2="C", alpha=0.25)

ud.set_leaves(['A', 'Admix', 'C', 'D'])

ud.set_ancestral_population_size(1.0e5)


############################## Simulation block ##############################

# n_per_leaf is a shortcut to specify equal number of samples per leaf population
# also support samples_per_pop dictionary to specify different sample sizes.
# diploid samples are created by default
mts, ud = simulate_msprime(
    ud,
    n_per_leaf=20,
    sequence_length=5e6,
    recombination_rate=1e-8,
    mutation_rate=1e-8,
    seed=42
)

# mts, dem = simulate_msprime(
#     ud,
#     samples_per_pop={"A":30, "Admix":20, "C":25, "D":15},
#     n_per_leaf=20,
#     sequence_length=5e6,
#     recombination_rate=1e-8,
#     mutation_rate=1e-8,
#     seed=42
# )

# LD prunning via PLINK
"""
Output the 


"""
pruned_sites = simulate_snp_pruning(mts, ud, temporal_folder='./plink_output')

# Get the observed sample covariance matrix of the pruned SNPs