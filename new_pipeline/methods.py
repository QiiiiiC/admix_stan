# demography_edgewise.py
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import numpy as np
import msprime
import subprocess
import tskit
from pathlib import Path
from itertools import combinations_with_replacement


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

