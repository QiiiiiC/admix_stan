# demography_edgewise.py
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import numpy as np
import msprime
import subprocess
import tskit
from pathlib import Path


@dataclass
class Node:
    name: str
    time: float  

@dataclass
class Edge:
    child: str
    parent: str
    Ne: float    

@dataclass
class Admixture:
    time: float                
    child_present: str         
    parents: Tuple[str, str]   
    alpha: float               
    Ne_child: float            

@dataclass
class UnifiedDemography:
    nodes: Dict[str, Node] = field(default_factory=dict)
    edges: List[Edge] = field(default_factory=list)
    admixtures: List[Admixture] = field(default_factory=list)
    leaves: List[str] = field(default_factory=list)
    ancestral_population_size: float = field(default_factory=lambda: 1e5)

    def add_population(self, name: str, time: float):
        self.nodes[name] = Node(name, float(time))

    def add_edge(self, child: str, parent: str, Ne: float):
        tp, tc = self.nodes[parent].time, self.nodes[child].time
        assert tp > tc, f"edge {parent}->{child}: parent({tp}) must be older than child({tc})"
        assert Ne > 0
        self.edges.append(Edge(child, parent, float(Ne)))

    def add_admixture(self, time: float, child_present: str, p1: str, p2: str,
                      alpha: float, Ne_child: float):
        assert 0.0 <= alpha <= 1.0
        assert Ne_child > 0
        self.admixtures.append(
            Admixture(float(time), child_present, (p1, p2), float(alpha), float(Ne_child))
        )

    def set_leaves(self, leaves: List[str]):
        self.leaves = list(leaves)

    def set_ancestral_population_size(self, Ne: float):
        self.ancestral_population_size = float(Ne)

    @property
    def ancestral_time(self) -> float:
        return max((n.time for n in self.nodes.values()), default=0.0)


def to_msprime(ud: UnifiedDemography) -> msprime.Demography:
    dem = msprime.Demography()

    # 1) add all populations with baseline size
    for name in ud.nodes:
        dem.add_population(name=name, initial_size=ud.ancestral_population_size)

    for e in ud.edges:
        t_child = ud.nodes[e.child].time
        dem.add_population_parameters_change(time=t_child, population=e.child, initial_size=e.Ne)

    for p in ud.admixtures:
        t_child = ud.nodes[p.child_present].time
        dem.add_population_parameters_change(time=t_child, population=p.child_present,
                                             initial_size=p.Ne_child)

    parent_children: Dict[Tuple[str, float], List[str]] = {}
    for e in ud.edges:
        parent_children.setdefault((e.parent, ud.nodes[e.parent].time), []).append(e.child)
    for (parent, t_parent), kids in parent_children.items():
        dem.add_population_split(time=t_parent, derived=kids, ancestral=parent)

    for p in ud.admixtures:
        dem.add_admixture(time=p.time, derived=p.child_present,
                          ancestral=[p.parents[0], p.parents[1]],
                          proportions=[p.alpha, 1.0 - p.alpha])

    dem.sort_events()
    return dem


def simulate_msprime(
    ud: UnifiedDemography,
    samples_per_pop: Optional[Dict[str, int]] = None,
    n_per_leaf: int = 20,
    sequence_length: float = 1e7,
    recombination_rate: float = 1e-8,
    mutation_rate: float = 1e-8,
    seed: Optional[int] = 1
):
    # default: sample from time-0 nodes if leaves not provided
    if not ud.leaves:
        ud.set_leaves([name for name, node in ud.nodes.items() if node.time == 0.0])

    if samples_per_pop is None:
        samples_per_pop = {pop: n_per_leaf for pop in ud.leaves}

    dem = to_msprime(ud)
    samples = [msprime.SampleSet(samples_per_pop[p], population=p, time=0) for p in ud.leaves]

    ts = msprime.sim_ancestry(
        samples=samples,
        demography=dem,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=seed,
    )
    mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seed, discrete_genome=False)
    return mts, ud


def simulate_snp_pruning(
        mts: tskit.TreeSequence,
        ud: UnifiedDemography,
        temporal_folder: str,
        pruning_blocksize:int = 100,
        pruning_blockstep:int = 5,
        pruning_r2:float = 0.2
):
    p = Path.cwd() / temporal_folder
    p.mkdir(parents=True, exist_ok=True)
    # In SNP only consider ancestral mutations before the ancient time.
    anc_t = ud.ancestral_time

    ancestral_mutations = []

    for mut in mts.mutations():
        node = mut.node  
        time = mts.node(node).time  
        
        if time > anc_t:  
            ancestral_mutations.append(mut.site)
    filtered_ts = mts.delete_sites([mut.site for mut in mts.mutations() if mut.id not in ancestral_mutations])
    with open(f"{temporal_folder}/snp_ts.vcf", "w") as vcf_file:
        filtered_ts.write_vcf(vcf_file)

    vcf_file = f"{temporal_folder}/snp_ts.vcf"
    output_prefix = f"{temporal_folder}/plink_output"
    pruned_prefix = f"{temporal_folder}/pruned"

    def run_cmd(cmd):
        print(f"Running: {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    # Convert VCF to PLINK binary
    run_cmd(f"plink --vcf {vcf_file} --make-bed --out {output_prefix} --double-id")

    # LD pruning
    run_cmd(f"plink --bfile {output_prefix} --indep-pairwise {pruning_blocksize} {pruning_blockstep} {pruning_r2} --out {pruned_prefix}")

    # Extract pruned SNPs
    run_cmd(f"plink --bfile {output_prefix} --extract {pruned_prefix}.prune.in --make-bed --out {pruned_prefix}_final")

    print("LD pruning complete.")

    with open(f"{pruned_prefix}.prune.in", "r") as file:
        pruned_snps_list = [line.strip() for line in file]
    pruned_snps_list = [int(i) for i in pruned_snps_list]



    for variant in filtered_ts.variants():
        if variant.index in pruned_snps_list:
            
            genotypes = variant.genotypes  
            freq_pop = np.sum(genotypes)/filtered_ts.num_individuals
            if freq_pop > 0.1:
                freq_array = [np.sum(genotypes[sample]) / len(sample) for sample in popn_samples]

            if freq_pop > 0.1:
                final_pruned.append(variant.index)
                freq_.append(freq_pop)
                freq_array.append([freq_pop1, freq_pop2,freq_pop3,freq_pop4])

    freq_array = np.array(freq_array)

    return pruned_snps_list




def snp_data(
        ud: UnifiedDemography,
        freq_array: np.ndarray
):
    n_popn, n_snp = len(freq_array[0]), len(freq_array)
    M = np.eye(n_popn) - np.ones((n_popn, n_popn)) / n_popn # centering matrix
    W_hat = (freq_array @ M).T @ (freq_array @ M) / n_snp
   
    
    return W_hat

def snp_data_blocked(
        ud: UnifiedDemography,
        freq_array: np.ndarray,
        block_size: int
):
    n_popn, n_snp = len(freq_array[0]), len(freq_array)
    M = np.eye(n_popn) - np.ones((n_popn, n_popn)) / n_popn # centering matrix
    W_blocks = []
    n_blocks = n_snp // block_size
    for i in range(n_blocks):
        start = i * block_size
        end = start + block_size
        block = freq_array[:, start:end]
        W_block = (block @ M).T @ (block @ M) / block_size
        W_blocks.append(W_block)
    W_blocks = np.array(W_blocks)
    W_hat = np.mean(W_blocks, axis=0)
    
    return W_hat, W_blocks
