import msprime
import tskit
import numpy as np
import subprocess
from pathlib import Path
from typing import Dict, Optional, List, Tuple
from collections import defaultdict
import itertools

def to_msprime_demography(topology):
    """
    Converts the custom DemographicTopology into an msprime.Demography object.
    """
    # 1. Initialize msprime Demography
    demography = msprime.Demography()

    # 2. Define Populations (Nodes)
    # msprime requires us to register populations before using them in events.
    # We use the 'ne' attribute from your Node as 'initial_size'.
    for name, node in topology.nodes.items():
        if node.ne is None:
            # Fallback if Ne isn't set (though it should be for valid sims)
            print(f"Warning: Node {name} has no Ne set. Defaulting to 1.")
            ne_val = 1
        else:
            ne_val = node.ne
        
        demography.add_population(name=name, initial_size=ne_val)

    # 3. Process Events in Order
    # Your events are ordered. We iterate through them and map them to msprime calls.
    for event in topology.ordered_events:
        
        if event['type'] == 'MERGE':
            # Logic: Two populations (children) merge into one (parent) backwards in time.
            # In msprime: add_population_split(time, derived=[c1, c2], ancestral=parent)
            c1, c2 = event['children']
            parent = event['parent']
            
            # The event time is the START of the parent node
            time = topology.nodes[parent].time_start
            
            if time is None:
                raise ValueError(f"Time not set for merge event into {parent}")
                
            demography.add_population_split(
                time=time,
                derived=[c1, c2],
                ancestral=parent
            )

        elif event['type'] == 'ADMIXTURE':
            # Logic: One population (child) splits into two (parents) backwards in time.
            # In msprime: add_admixture(time, derived=child, ancestral=[p1, p2], proportions=...)
            child = event['child']
            p1, p2 = event['parents']
            
            # The event time is the END of the child node (when it looks back and splits)
            time = topology.nodes[child].time_end
            
            if time is None:
                raise ValueError(f"Time not set for admixture event of {child}")

            # Retrieve fractions
            # Note: msprime expects proportions in the order of 'ancestral' list [p1, p2]
            frac1 = topology.nodes[child].admixture_fractions.get(p1)
            frac2 = topology.nodes[child].admixture_fractions.get(p2)
            
            if frac1 is None or frac2 is None:
                raise ValueError(f"Admixture fractions not set for {child}")

            demography.add_admixture(
                time=time,
                derived=child,
                ancestral=[p1, p2],
                proportions=[frac1, frac2]
            )

    # 4. Sort events
    # msprime requires events to be added, but sometimes internal sorting is needed 
    # if the input order wasn't strictly chronological. 
    # However, Demography.sort_events() usually handles this automatically before simulation.
    demography.sort_events()
    
    return demography

def simulate_msprime(
    topology, # this is the class DemographicTopology
    samples_per_pop= None,
    n_per_leaf=20,
    sequence_length = 1e7,
    recombination_rate = 1e-8,
    mutation_rate = 1e-8,
    seed = 42
):
    """
    Simulates ancestry and mutations using msprime based on the topology.
    """
    # Default: sample from initial leaves
    leaves = topology.initial_leaves
    
    if samples_per_pop is None:
        samples_per_pop = {pop: n_per_leaf for pop in leaves}

    # Convert to msprime demography
    dem = to_msprime_demography(topology)
    
    # Create SampleSet
    # Note: We assume leaves are at time 0. If your topology has ancient leaves, 
    # you might need to check node.time_start.
    samples = []
    for p in leaves:
        if p in samples_per_pop:
            samples.append(msprime.SampleSet(samples_per_pop[p], population=p, time=0))

    # Simulate Ancestry
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=dem,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=seed,
        model = 'smc'
    )
    
    # Simulate Mutations
    mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seed, discrete_genome=True)
    
    return mts


def simulate_snp_pruning(
        mts,
        topology,
        temporal_folder,
        cutoff_time = None, 
        cutoff_freq = 0.1,
        pruning_blocksize = 100,
        pruning_blockstep = 5,
        pruning_r2 = 0.2
):
    """
    Filters for ancient mutations, exports to VCF, runs PLINK for LD pruning, 
    and returns the pruned SNP indices and their frequency array.
    """
    p = Path(temporal_folder)
    p.mkdir(parents=True, exist_ok=True)
    
    # 1. Filter for Ancestral Mutations
    if cutoff_time is None:
        # Default: Keep everything
        filtered_ts = mts
    else:
        # Keep only mutations OLDER than cutoff_time
        sites_to_remove = []
        for mut in mts.mutations():
            if mts.node(mut.node).time < cutoff_time:
                sites_to_remove.append(mut.site)
        
        if sites_to_remove:
            filtered_ts = mts.delete_sites(sites_to_remove)
        else:
            filtered_ts = mts

    # 2. Export VCF for PLINK
    vcf_path = p / "snp_ts.vcf"
    with open(vcf_path, "w") as vcf_file:
        filtered_ts.write_vcf(vcf_file)

    output_prefix = p / "snp_ts"
    pruned_prefix = p / "pruned"

    def run_cmd(cmd):
        print(f"Running: {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    # 3. Run PLINK Commands
    # Convert VCF to BED
    run_cmd(f"plink --vcf {vcf_path} --make-bed --out {output_prefix} --double-id --allow-extra-chr --silent")

    # LD pruning
    run_cmd(f"plink --bfile {output_prefix} --indep-pairwise {pruning_blocksize} {pruning_blockstep} {pruning_r2} --out {pruned_prefix} --allow-extra-chr --silent")

    # Extract pruned SNPs list
    prune_in_file = p / "pruned.prune.in"
    if not prune_in_file.exists():
        raise FileNotFoundError(f"PLINK did not generate {prune_in_file}. Check if plink is installed.")
        
    with open(prune_in_file, "r") as file:
        pruned_ids = set(line.strip() for line in file)

    print(f"LD pruning complete. Kept {len(pruned_ids)} SNPs.")

    # 4. Compute Frequencies for Pruned SNPs
    # Map population names to sample indices in the tree sequence
    leaves = topology.initial_leaves
    
    # Build dictionary of {pop_name: [sample_ids]}
    pop_sample_ids = {}
    for pop_name in leaves:
        # Get the ID of the population from the tree sequence
        # Note: msprime demography keeps the names
        pop_id = -1
        for p in filtered_ts.populations():
            if p.metadata.get('name') == pop_name:
                pop_id = p.id
                break
        
        if pop_id != -1:
            pop_sample_ids[pop_name] = filtered_ts.samples(population=pop_id)
        else:
            pop_sample_ids[pop_name] = []

    final_pruned_indices = []
    freq_matrix = []

    # Iterate variants and check if they are in the pruned set
    for variant in filtered_ts.variants():
        # The VCF ID is usually the variant.site.id or similar. 
        # tskit write_vcf uses site IDs by default.
        # We assume the ID in pruned_ids matches str(variant.site.id)
        # Note: variant.index is the index in the list, variant.site.id is the permanent ID.
        
        if str(variant.index) in pruned_ids:
            
            genotypes = variant.genotypes
            total_samples = len(genotypes)
            global_freq = np.sum(genotypes) / total_samples
            
            # Optional: Global MAF filter (e.g. > 0.1)
            if global_freq > cutoff_freq:
                # Calculate frequency per population
                row = []
                for pop in leaves:
                    s_indices = pop_sample_ids[pop]
                    if len(s_indices) > 0:
                        # Extract genotypes for this population
                        pop_genos = genotypes[s_indices]
                        freq = np.sum(pop_genos) / len(pop_genos)
                        row.append(freq)
                    else:
                        row.append(0.0)
                
                freq_matrix.append(row)
                final_pruned_indices.append(variant.site.id)

    freq_array = np.array(freq_matrix) 

    return final_pruned_indices, freq_array


def calculate_ibd_fractions_path(ts, bins, cm_per_unit=1e-6, num_bootstraps=1000):    

    
    sample_nodes = ts.samples()
    node_to_pop = ts.nodes_population[sample_nodes]
    pop_ids = np.unique(node_to_pop)
    num_pops = len(pop_ids)
    

    pop_samples = defaultdict(list)
    for u in sample_nodes:
        pop_samples[node_to_pop[u]].append(u)

    results = {b_i: defaultdict(lambda: defaultdict(float)) 
               for b_i in range(len(bins))}
    
    genome_length = ts.sequence_length * cm_per_unit
    
    # Filter tiny segments
    min_bin_val = min(b[0] for b in bins)
    min_span_ts_units = min_bin_val / cm_per_unit
    
    ibd_iter = ts.ibd_segments(
        store_pairs=True, 
        store_segments=True,
        min_span=min_span_ts_units  
    )

    print(f"Iterating IBD segments (min_span={min_span_ts_units:.2f})...")
    
    for (u, v), segments in ibd_iter.items(): 
        p_u = ts.nodes_population[u]
        p_v = ts.nodes_population[v]
        p_i, p_j = sorted((p_u, p_v))
        
        # Unique identifier for this specific pair of individuals
        pair_key = tuple(sorted((u, v)))
        
        for seg in segments:
            seg_len = (seg.right - seg.left) * cm_per_unit
            
            for b_i, (min_len, max_len) in enumerate(bins):
                if min_len < seg_len <= max_len:
                    # FIX 3: Store fraction for the pair using tuple key (p_i, p_j)
                    # This matches how the bootstrap loop tries to retrieve it later.
                    results[b_i][(p_i, p_j)][pair_key] += (seg_len / genome_length)
                    break 
    
    final_mean_matrix = {}
    final_var_matrix = {}

    for b_i in results:
        mean_matrix = np.zeros((num_pops, num_pops))
        var_matrix = np.zeros((num_pops, num_pops))

        for i in range(num_pops):
            for j in range(i, num_pops):
                # Calculate total theoretical pairs (N)
                # Now this works because pop_samples[i] is a list
                if i == j:
                    n = len(pop_samples[i]) #
                    num_pairs = n * (n - 1) // 2
                else:
                    num_pairs = len(pop_samples[i]) * len(pop_samples[j])
                
                if num_pairs == 0:
                    continue

                # Retrieve observed non-zero fractions
                # This works now because we stored data with key (i, j)
                observed_dict = results[b_i].get((i, j), {})
                observed_values = np.array(list(observed_dict.values()))
                
                # The rest are zeros
                count_zeros = num_pairs - len(observed_values)
                
                # Construct the full population of pairs
                full_population = np.concatenate([
                    observed_values, 
                    np.zeros(count_zeros)
                ])

                # A. Original Mean
                original_mean = np.mean(full_population)
                
                # B. Bootstrap Variance
                if num_bootstraps > 0:
                    boot_samples = np.random.choice(full_population, size=(num_bootstraps, num_pairs), replace=True)
                    boot_means = np.mean(boot_samples, axis=1)
                    boot_var = np.var(boot_means)
                else:
                    boot_var = 0.0

                # Fill Matrices
                mean_matrix[i, j] = mean_matrix[j, i] = original_mean
                var_matrix[i, j] = var_matrix[j, i] = boot_var

        final_mean_matrix[b_i] = mean_matrix
        final_var_matrix[b_i] = var_matrix

    return final_mean_matrix, final_var_matrix


def calculate_ibd_fractions_mrca(ts, bins, cm_per_unit=1e-6, num_bootstraps=1000):
    """
    Calculates IBD fractions using the strict 'MRCA-span' definition.
    Segments are defined by continuous TMRCA, ignoring changes in topology 
    that do not alter the most recent common ancestor.
    """
    
    # 1. Setup Population Mappings
    sample_nodes = ts.samples()
    num_samples = len(sample_nodes)
    node_to_pop = ts.nodes_population[sample_nodes]
    pop_ids = np.unique(node_to_pop)
    num_pops = len(pop_ids)

    pop_samples = defaultdict(list)
    for u in sample_nodes:
        pop_samples[node_to_pop[u]].append(u)

    # 2. Initialize Results Container
    # structure: results[bin_index][(pop_i, pop_j)][pair_key] = total_fraction
    results = {b_i: defaultdict(lambda: defaultdict(float)) 
               for b_i in range(len(bins))}
    
    genome_length_cm = ts.sequence_length * cm_per_unit

    # 3. Iterate Trees to Extract MRCA-span Segments
    # We maintain the state of the current segment for every pair of samples.
    # state[(u, v)] = (current_mrca_node, start_coordinate_bp)
    current_state = {}
    
    # Initialize with the first tree
    tree_iter = ts.trees()
    first_tree = next(tree_iter)
    
    # Pre-calculate pairs to iterate (u, v)
    # Using indices to access the sample_nodes array is faster
    pairs = []
    for i in range(num_samples):
        for j in range(i + 1, num_samples):
            u, v = sample_nodes[i], sample_nodes[j]
            pairs.append((u, v))
            
            # Initialize state
            mrca = first_tree.mrca(u, v)
            current_state[(u, v)] = (mrca, first_tree.interval.left)

    print(f"Scanning trees for MRCA segments ({len(pairs)} pairs)...")

    # Helper to process a finished segment
    def process_segment(u, v, length_bp):
        length_cm = length_bp * cm_per_unit
        
        # Binning logic
        for b_i, (min_len, max_len) in enumerate(bins):
            # Note: Strict inequality matching your request (min < len <= max)
            # Adjust if you need inclusive lower bound
            if min_len < length_cm < max_len:
                p_u = node_to_pop[u]
                p_v = node_to_pop[v]
                p_i, p_j = sorted((p_u, p_v))
                pair_key = tuple(sorted((u, v)))
                
                results[b_i][(p_i, p_j)][pair_key] += (length_cm / genome_length_cm)
                break

    # Iterate through the rest of the trees
    for tree in tree_iter:
        current_left = tree.interval.left
        
        for (u, v) in pairs:
            new_mrca = tree.mrca(u, v)
            old_mrca, start_pos = current_state[(u, v)]
            
            # MRCA changed? Segment ended.
            if new_mrca != old_mrca:
                seg_len = current_left - start_pos
                process_segment(u, v, seg_len)
                
                # Start new segment
                current_state[(u, v)] = (new_mrca, current_left)

    # 4. Flush the final segments (end of genome)
    final_pos = ts.sequence_length
    for (u, v) in pairs:
        old_mrca, start_pos = current_state[(u, v)]
        seg_len = final_pos - start_pos
        process_segment(u, v, seg_len)

    # 5. Matrix Construction & Bootstrapping (Same as your original code)
    final_mean_matrix = {}
    final_var_matrix = {}

    for b_i in results:
        mean_matrix = np.zeros((num_pops, num_pops))
        var_matrix = np.zeros((num_pops, num_pops))

        for i in range(num_pops):
            for j in range(i, num_pops):
                # Calculate total theoretical pairs
                n_i = len(pop_samples[i])
                n_j = len(pop_samples[j])
                
                if i == j:
                    num_pairs = n_i * (n_i - 1) // 2
                else:
                    num_pairs = n_i * n_j
                
                if num_pairs == 0:
                    continue

                # Retrieve observed non-zero fractions
                observed_dict = results[b_i].get((i, j), {})
                observed_values = np.array(list(observed_dict.values()))
                
                # The rest are zeros
                count_zeros = num_pairs - len(observed_values)
                
                # Construct the full population of pairs
                full_population = np.concatenate([
                    observed_values, 
                    np.zeros(count_zeros)
                ])

                # A. Original Mean
                original_mean = np.mean(full_population)
                
                # B. Bootstrap Variance
                if num_bootstraps > 0:
                    boot_samples = np.random.choice(full_population, size=(num_bootstraps, num_pairs), replace=True)
                    boot_means = np.mean(boot_samples, axis=1)
                    boot_var = np.var(boot_means)
                else:
                    boot_var = 0.0

                # Fill Matrices
                mean_matrix[i, j] = mean_matrix[j, i] = original_mean
                var_matrix[i, j] = var_matrix[j, i] = boot_var

        final_mean_matrix[b_i] = mean_matrix
        final_var_matrix[b_i] = var_matrix

    return final_mean_matrix, final_var_matrix


def _resolve_n_samples(dem, n_samples_per_pop):
    """Resolve per-population haploid sample counts into a list ordered to
    match the IBD/SNP matrix rows (i.e. dem.initial_leaves order).

    The IBD zero-segment likelihood normalises by the number of sample pairs,
    which depends on the number of haploid samples in each population:
    within a population (i == j) there are n*(n-1)/2 pairs, and between two
    populations (i != j) there are n_i * n_j pairs.

    Parameters
    ----------
    n_samples_per_pop : int | dict | sequence
        Number of *haploid* samples per population. Accepts
          - dict {leaf_name: count}  (looked up by dem.initial_leaves),
          - int                      (same count for every population), or
          - sequence of length n_leaves (already in leaf order).
    """
    leaves = dem.initial_leaves
    n_leaves = len(leaves)

    if n_samples_per_pop is None:
        raise ValueError(
            "n_samples_per_pop must be provided (haploid sample count per "
            "population) so the IBD pair counts reflect the actual sample sizes."
        )

    if isinstance(n_samples_per_pop, dict):
        return [int(n_samples_per_pop[p]) for p in leaves]
    if np.isscalar(n_samples_per_pop):
        return [int(n_samples_per_pop)] * n_leaves

    seq = [int(x) for x in n_samples_per_pop]
    if len(seq) != n_leaves:
        raise ValueError(
            f"n_samples_per_pop has length {len(seq)} but the topology has "
            f"{n_leaves} leaves."
        )
    return seq


def build_ibd_stan_data(
    dem,
    ibd_mean: dict,
    ibd_var: dict,
    bins: list,
    n_samples_per_pop=None,
    T_max: float = None,
    se_floor: float = 1e-8,
    cm: float = None
) -> dict:
    """
    Build a Stan data dictionary for the CORRECTED IBD model that tracks
    joint lineage paths instead of marginal locations.

    The key change from the old version: we now provide precomputed leaf-pair
    indices (pair_i, pair_j) so the Stan model can iterate over unique leaf
    pairs and track joint (lineage_1_location, lineage_2_location) states.

    Parameters
    ----------
    dem : DemographicTopology
        Topology of populations.

    ibd_mean : dict
        {bin_index: np.ndarray (n_leaves, n_leaves)} — mean IBD fractions
        per population pair per bin.

    ibd_var : dict
        {bin_index: np.ndarray (n_leaves, n_leaves)} — bootstrap *variances*.

    bins : list of [min_cM, max_cM]
        Bin boundaries in centiMorgans.

    T_max : float, optional
        Upper bound for the final open-ended epoch (generations).

    se_floor : float
        Minimum SE to avoid zero-variance likelihoods.

    Returns
    -------
    ibd_data : dict
        Ready to pass to cmdstanpy model.sample(data=ibd_data).
    """
    ibd_data = {}

    # ------------------------------------------------------------------
    # 1. Topology block (unchanged)
    # ------------------------------------------------------------------
    matrices, events, admixture_map, admixture_map_id = \
        dem.get_topology_matrix_representation()

    n_leaves = len(dem.initial_leaves)
    n_nodes = len(dem.nodes)
    n_events = len(dem.ordered_events)

    ibd_data['n_leaves']    = n_leaves
    ibd_data['n_nodes']     = n_nodes
    ibd_data['n_events']    = n_events
    ibd_data['n_admixture'] = dem.n_admix
    ibd_data['migration_matrices'] = matrices

    # admixture_map: [event_index, src, tgt1, tgt2]
    m = []
    for i in admixture_map_id.keys():
        m.append(
            [i, admixture_map_id[i]['child']] +
            admixture_map_id[i]['parents']
        )
    ibd_data['admixture_map'] = m

    # Index arrays: 1-based to match Stan
    ibd_data['admixture_indices']     = [i + 2 for i in admixture_map_id.keys()]
    ibd_data['fixed_indices_shifted'] = [
        i for i in range(2, n_events + 2)
        if i not in ibd_data['admixture_indices']
    ]
    ibd_data['fixed_indices'] = [i - 1 for i in ibd_data['fixed_indices_shifted']]

    # ------------------------------------------------------------------
    # 2. T_max
    # ------------------------------------------------------------------
    ibd_data['T_max'] = T_max

    # ------------------------------------------------------------------
    # 3. Bin boundaries
    # ------------------------------------------------------------------
    n_bins = len(bins)
    ibd_data['n_bins']      = n_bins
    ibd_data['bin_length']  = [list(b) for b in bins]

    # ------------------------------------------------------------------
    # 4.: Precompute leaf pair indices for the corrected model
    #    We enumerate all unique pairs (i, j) with i <= j (1-based).
    #    This includes self-pairs (i, i) for diagonal IBD entries.
    # ------------------------------------------------------------------
    pair_i_list = []
    pair_j_list = []
    for i in range(1, n_leaves + 1):
        for j in range(i, n_leaves + 1):
            pair_i_list.append(i)
            pair_j_list.append(j)

    n_leaf_pairs = len(pair_i_list)
    assert n_leaf_pairs == n_leaves * (n_leaves + 1) // 2

    ibd_data['n_leaf_pairs'] = n_leaf_pairs
    ibd_data['pair_i'] = pair_i_list
    ibd_data['pair_j'] = pair_j_list

    # Haploid sample count per population, for the IBD zero-segment pair counts
    ibd_data['n_samples'] = _resolve_n_samples(dem, n_samples_per_pop)

    # ------------------------------------------------------------------
    # 5. IBD observations: ibd_hat and ibd_se
    # ------------------------------------------------------------------
    assert len(ibd_mean) == n_bins, \
        f"ibd_mean has {len(ibd_mean)} entries but bins has {n_bins}"
    assert ibd_mean[0].shape == (n_leaves, n_leaves), \
        f"Expected ibd_mean matrices of shape ({n_leaves},{n_leaves}), " \
        f"got {ibd_mean[0].shape}."

    ibd_hat = np.zeros((n_bins, n_leaves, n_leaves))
    ibd_se  = np.zeros((n_bins, n_leaves, n_leaves))

    for b in range(n_bins):
        ibd_hat[b] = ibd_mean[b]
        se = np.sqrt(np.maximum(ibd_var[b], 0.0))
        ibd_se[b]  = np.maximum(se, se_floor)

    ibd_data['ibd_hat'] = ibd_hat.tolist()
    ibd_data['ibd_se']  = ibd_se.tolist()

    # ------------------------------------------------------------------
    # 6. cM simulated
    # ------------------------------------------------------------------    
    ibd_data['cm'] = cm

    # ------------------------------------------------------------------
    # 7. Summary
    # ------------------------------------------------------------------
    print(f"[build_ibd_stan_data] Stan data assembled (corrected joint-path model):")
    print(f"  n_leaves     : {n_leaves}")
    print(f"  n_nodes      : {n_nodes}")
    print(f"  n_events     : {n_events}")
    print(f"  n_admixture  : {ibd_data['n_admixture']}")
    print(f"  n_samples    : {ibd_data['n_samples']}")
    print(f"  n_leaf_pairs : {n_leaf_pairs}")
    print(f"  n_bins       : {n_bins}")
    print(f"  T_max        : {T_max:.0f} generations")
    print(f"  ibd_hat      : [{ibd_hat.min():.3e}, {ibd_hat.max():.3e}]")
    print(f"  ibd_se       : [{ibd_se.min():.3e},  {ibd_se.max():.3e}]")

    return ibd_data


def build_snp_stan_data(
    dem,
    w_hat: np.ndarray,
    w_se: np.ndarray,
) -> dict:
    """
    Build a Stan data dictionary for snp_model_Nfixed.stan.

    Parameters
    ----------
    dem : DemographicTopology
    w_hat : np.ndarray (n_leaves, n_leaves)
        Observed double-centered SNP covariance matrix.
    w_se : np.ndarray (n_leaves, n_leaves)
        Standard errors for w_hat entries.

    Returns
    -------
    dict : ready for cmdstanpy.
    """
    data = {}

    # -- Topology (same as IBD) --
    matrices, events, admixture_map, admixture_map_id = \
        dem.get_topology_matrix_representation()

    n_leaves = len(dem.initial_leaves)
    n_nodes = len(dem.nodes)
    n_events = len(dem.ordered_events)

    data['n_leaves']    = n_leaves
    data['n_nodes']     = n_nodes
    data['n_events']    = n_events
    data['n_admixture'] = dem.n_admix
    data['migration_matrices'] = matrices

    m = []
    for i in admixture_map_id.keys():
        m.append(
            [i, admixture_map_id[i]['child']] +
            admixture_map_id[i]['parents']
        )
    data['admixture_map'] = m

    data['admixture_indices']     = [i + 2 for i in admixture_map_id.keys()]
    data['fixed_indices_shifted'] = [
        i for i in range(2, n_events + 2)
        if i not in data['admixture_indices']
    ]
    data['fixed_indices'] = [i - 1 for i in data['fixed_indices_shifted']]

    # -- Leaf pairs --
    pair_i_list = []
    pair_j_list = []
    for i in range(1, n_leaves + 1):
        for j in range(i, n_leaves + 1):
            pair_i_list.append(i)
            pair_j_list.append(j)

    data['n_leaf_pairs'] = len(pair_i_list)
    data['pair_i'] = pair_i_list
    data['pair_j'] = pair_j_list

    # -- Observations --
    assert w_hat.shape == (n_leaves, n_leaves), \
        f"w_hat shape {w_hat.shape} != ({n_leaves},{n_leaves})"
    data['w_hat'] = w_hat.tolist()
    data['w_se']  = w_se.tolist()

    print(f"[build_snp_stan_data] Stan data assembled:")
    print(f"  n_leaves     : {n_leaves}")
    print(f"  n_nodes      : {n_nodes}")
    print(f"  n_events     : {n_events}")
    print(f"  n_leaf_pairs : {data['n_leaf_pairs']}")
    print(f"  w_hat range  : [{w_hat.min():.3e}, {w_hat.max():.3e}]")
    print(f"  w_se  range  : [{w_se.min():.3e}, {w_se.max():.3e}]")

    return data


def build_mixed_stan_data(
    dem,
    ibd_mean: dict,
    ibd_var: dict,
    bins: list,
    w_hat: np.ndarray,
    w_se: np.ndarray,
    n_samples_per_pop=None,
    T_max: float = None,
    se_floor: float = 1e-8,
    cm: float = None,
    ibd_count: dict = None,
) -> dict:
    """
    Build a Stan data dictionary for mixed_model.stan (IBD + SNP composite
    likelihood with effective_N as a parameter).

    Parameters
    ----------
    dem : DemographicTopology
    ibd_mean : dict
        {bin_index: np.ndarray (n_leaves, n_leaves)} — mean IBD fractions.
    ibd_var : dict
        {bin_index: np.ndarray (n_leaves, n_leaves)} — IBD variances.
    bins : list of [min_cM, max_cM]
    w_hat : np.ndarray (n_leaves, n_leaves)
        Observed double-centered SNP covariance matrix.
    w_se : np.ndarray (n_leaves, n_leaves)
        Standard errors for w_hat entries.
    T_max : float
        Upper bound for the final IBD epoch.
    se_floor : float
        Minimum SE for IBD observations.
    cm : float
        Genome length in cM for Poisson zero-IBD likelihood.

    Returns
    -------
    dict : ready for cmdstanpy.
    """
    data = {}

    # -- Topology (shared) --
    matrices, events, admixture_map, admixture_map_id = \
        dem.get_topology_matrix_representation()

    n_leaves = len(dem.initial_leaves)
    n_nodes = len(dem.nodes)
    n_events = len(dem.ordered_events)

    data['n_leaves']    = n_leaves
    data['n_nodes']     = n_nodes
    data['n_events']    = n_events
    data['n_admixture'] = dem.n_admix
    data['migration_matrices'] = matrices

    m = []
    for i in admixture_map_id.keys():
        m.append(
            [i, admixture_map_id[i]['child']] +
            admixture_map_id[i]['parents']
        )
    data['admixture_map'] = m

    data['admixture_indices']     = [i + 2 for i in admixture_map_id.keys()]
    data['fixed_indices_shifted'] = [
        i for i in range(2, n_events + 2)
        if i not in data['admixture_indices']
    ]
    data['fixed_indices'] = [i - 1 for i in data['fixed_indices_shifted']]

    # -- Leaf pairs (shared) --
    pair_i_list = []
    pair_j_list = []
    for i in range(1, n_leaves + 1):
        for j in range(i, n_leaves + 1):
            pair_i_list.append(i)
            pair_j_list.append(j)

    data['n_leaf_pairs'] = len(pair_i_list)
    data['pair_i'] = pair_i_list
    data['pair_j'] = pair_j_list

    # Haploid sample count per population, for the IBD zero-segment pair counts
    data['n_samples'] = _resolve_n_samples(dem, n_samples_per_pop)

    # -- IBD-specific data --
    data['T_max'] = T_max

    n_bins = len(bins)
    data['n_bins']     = n_bins
    data['bin_length'] = [list(b) for b in bins]

    assert len(ibd_mean) == n_bins, \
        f"ibd_mean has {len(ibd_mean)} entries but bins has {n_bins}"

    ibd_hat = np.zeros((n_bins, n_leaves, n_leaves))
    ibd_se  = np.zeros((n_bins, n_leaves, n_leaves))

    for b in range(n_bins):
        ibd_hat[b] = ibd_mean[b]
        se = np.sqrt(np.maximum(ibd_var[b], 0.0))
        ibd_se[b]  = np.maximum(se, se_floor)

    data['ibd_hat'] = ibd_hat.tolist()
    data['ibd_se']  = ibd_se.tolist()
    # Raw segment counts, for the Poisson likelihood.  ibd_hat/ibd_se are still
    # passed: the generated quantities report both parameterisations so the
    # Poisson fit can be read against the old normal one.
    if ibd_count is not None:
        data['ibd_count'] = np.stack([np.asarray(ibd_count[b], int)
                                      for b in range(n_bins)]).tolist()
    else:
        data['ibd_count'] = np.zeros((n_bins, n_leaves, n_leaves), int).tolist()
    data['cm'] = cm

    # -- SNP-specific data --
    assert w_hat.shape == (n_leaves, n_leaves), \
        f"w_hat shape {w_hat.shape} != ({n_leaves},{n_leaves})"
    data['w_hat'] = w_hat.tolist()
    data['w_se']  = w_se.tolist()

    # -- Summary --
    print(f"[build_mixed_stan_data] Stan data assembled (IBD + SNP composite):")
    print(f"  n_leaves     : {n_leaves}")
    print(f"  n_samples    : {data['n_samples']}")
    print(f"  n_nodes      : {n_nodes}")
    print(f"  n_events     : {n_events}")
    print(f"  n_leaf_pairs : {data['n_leaf_pairs']}")
    print(f"  n_bins       : {n_bins}")
    print(f"  T_max        : {T_max:.0f} generations")
    print(f"  cm           : {cm}")
    print(f"  ibd_hat      : [{ibd_hat.min():.3e}, {ibd_hat.max():.3e}]")
    print(f"  ibd_se       : [{ibd_se.min():.3e}, {ibd_se.max():.3e}]")
    print(f"  w_hat range  : [{w_hat.min():.3e}, {w_hat.max():.3e}]")
    print(f"  w_se  range  : [{w_se.min():.3e}, {w_se.max():.3e}]")

    return data