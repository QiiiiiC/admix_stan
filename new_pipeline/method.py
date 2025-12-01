import msprime
import tskit
import numpy as np
from typing import Dict, List, Optional
from demography import DemographicTopology
def simulate_msprime_from_topology(
    topology: DemographicTopology,
    n_samples_per_leaf: int = 20,
    seq_length: float = 1e7,
    recombination_rate: float = 1e-8,
    mutation_rate: float = 1e-8,
    global_ne: float = 10000,
    seed: int = 42
):
    """
    Translates the custom DemographicTopology into an msprime simulation.
    
    Args:
        topology: The instantiated DemographicTopology object.
        n_samples_per_leaf: Number of diploid samples to take from each leaf node (time 0).
        global_ne: Fallback effective population size if a Node.ne is None.
    
    Returns:
        mts: The mutated TreeSequence.
        demography: The msprime.Demography object used (for debugging/plotting).
    """
    
    # 1. Initialize msprime Demography
    demography = msprime.Demography()

    # 2. Register Populations (Nodes)
    # In msprime, every Node in your graph acts as a distinct population 
    # that exists for a specific duration.
    for name, node in topology.nodes.items():
        # Use node.ne if set, otherwise fallback to global_ne
        pop_ne = node.ne if node.ne is not None else global_ne
        
        # We define the population. 
        # Note: We don't strictly need to set initially_active=False because
        # msprime infers activity from events, but defining them all is safe.
        demography.add_population(
            name=name, 
            initial_size=pop_ne
        )

    # 3. Translate Events
    # We iterate through your ordered_events list to apply transitions.
    for event in topology.ordered_events:
        e_type = event['type']
        
        if e_type == 'MERGE':
            # LOGIC: Two children merge backwards into one parent.
            # MSPrime Perspective: A Parent splits forwards into two children.
            child_1 = event['children'][0]
            child_2 = event['children'][1]
            parent = event['parent']
            
            # The event happens at the Parent's start time (which is Children's end time)
            time = topology.get_node(parent).time_start
            
            if time is None:
                raise ValueError(f"Time not set for MERGE event creating {parent}")

            demography.add_population_split(
                time=time,
                derived=[child_1, child_2],
                ancestral=parent
            )

        elif e_type == 'ADMIXTURE':
            # LOGIC: One child splits backwards into two parents.
            # MSPrime Perspective: Child 'derived' from mixture of Ancestors.
            child_name = event['child']
            parent_names = event['parents'] # List of [p1, p2]
            
            child_node = topology.get_node(child_name)
            time = child_node.time_end
            
            if time is None:
                raise ValueError(f"Time not set for ADMIXTURE event involving {child_name}")

            # Retrieve fractions stored in the child node
            # Ensure the order of proportions matches the order of parent_names
            proportions = []
            for p_name in parent_names:
                frac = child_node.admixture_fractions.get(p_name)
                if frac is None:
                    raise ValueError(f"Admixture fraction missing for {p_name} in {child_name}")
                proportions.append(frac)

            demography.add_admixture(
                time=time,
                derived=child_name,
                ancestral=parent_names,
                proportions=proportions
            )

    # 4. Sort Events
    # msprime requires events to be time-ordered. 
    demography.sort_events()

    # 5. Define Samples
    # We take samples from the 'initial_leaves' defined in your topology at time 0
    sample_list = []
    for leaf_name in topology.initial_leaves:
        sample_list.append(msprime.SampleSet(
            num_samples=n_samples_per_leaf, 
            population=leaf_name, 
            time=0
        ))

    # 6. Run Simulation
    print("Starting msprime simulation...")
    ts = msprime.sim_ancestry(
        samples=sample_list,
        demography=demography,
        sequence_length=seq_length,
        recombination_rate=recombination_rate,
        random_seed=seed
    )
    
    print("Adding mutations...")
    mts = msprime.sim_mutations(
        ts, 
        rate=mutation_rate, 
        random_seed=seed
    )
    
    return mts, demography