from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import numpy as np
import msprime
import subprocess
import tskit



class Node:
    """
    Represents a population.
    
    Terminology:
    - Ancestors (Parents): Nodes backwards in time.
    - Descendants (Children): Nodes closer to the present.
    """
    def __init__(self, name):
        self.name = name
        
        self.descendants = [] 
        
        self.ancestors = []  
        
        self.admixture_fractions = {} 

        # Parameters
        self.ne = None          # Effective population size
        self.time_start = None  # Time (generations ago) this node appears backwards
        self.time_end = None    # Time (generations ago) this node merges/splits

    def add_descendant(self, node):
        """Register a node that comes AFTER this one (closer to present)."""
        self.descendants.append(node)

    def add_ancestor(self, node, fraction = 1.0):
        """Register a node that comes BEFORE this one (further in past)."""
        self.ancestors.append(node)
        self.admixture_fractions[node.name] = fraction

    @property
    def is_root(self):
        return len(self.ancestors) == 0

    @property
    def is_leaf(self):
        return len(self.descendants) == 0
    
    @property
    def is_admixed(self):
        return len(self.ancestors) > 1
    

class DemographicTopology:
    def __init__(self, leaf_names):
        self.nodes = {}
        self.ordered_events = [] # Stores events in the order they are added (1, 2, 3...)
        self.initial_leaves = leaf_names # Keep track of starting nodes for validation
        
        # Step 1: Initialize present-day samples (Leaves)
        for name in leaf_names:
            node = Node(name)
            node.time_start = 0.0 
            self.nodes[name] = node

    def get_node(self, name):
        if name not in self.nodes:
            raise ValueError(f"Node '{name}' not found.")
        return self.nodes[name]

    # ==========================================
    # PHASE 2: TOPOLOGY (Define Order of Events)
    # ==========================================

    def add_merge_event(self, child_1_name, child_2_name, new_parent_name):
        """
        Event Type: MERGE (Backwards)
        Two lineages join into one ancestor.
        """
        if new_parent_name in self.nodes:
            raise ValueError(f"Node '{new_parent_name}' already exists.")

        child1 = self.get_node(child_1_name)
        child2 = self.get_node(child_2_name)
        parent = Node(new_parent_name)

        # Link Topology
        # 1. Child looks back at Parent (Ancestry)
        child1.add_ancestor(parent, fraction=1.0)
        child2.add_ancestor(parent, fraction=1.0)
        
        # 2. Parent looks forward at Children (Descendants)
        parent.add_descendant(child1)
        parent.add_descendant(child2)

        self.nodes[new_parent_name] = parent
        
        # Record Event
        event_id = len(self.ordered_events) + 1
        self.ordered_events.append({
            'id': event_id,
            'type': 'MERGE',
            'children': [child_1_name, child_2_name],
            'parent': new_parent_name
        })
        print(f"[Event {event_id}] Merge: {child_1_name} & {child_2_name} -> {new_parent_name}")

    def add_admixture_event(self, child_name, parent_1_name, parent_2_name):
        """
        Event Type: ADMIXTURE/SPLIT (Backwards)
        One lineage splits into two ancestral sources.
        Note: Fractions are not set here; they are parameters added in Step 3.
        """
        if parent_1_name in self.nodes or parent_2_name in self.nodes:
            raise ValueError("New parent nodes must have unique names.")

        child = self.get_node(child_name)
        p1 = Node(parent_1_name)
        p2 = Node(parent_2_name)

        # Link Topology
        # 1. Child looks back at Parents (Ancestry)
        # Initialize with None/Placeholder fraction, to be set in Step 3
        child.add_ancestor(p1, fraction=None)
        child.add_ancestor(p2, fraction=None)

        # 2. Parents look forward at Child (Descendants)
        p1.add_descendant(child)
        p2.add_descendant(child)

        self.nodes[parent_1_name] = p1
        self.nodes[parent_2_name] = p2

        # Record Event
        event_id = len(self.ordered_events) + 1
        self.ordered_events.append({
            'id': event_id,
            'type': 'ADMIXTURE',
            'child': child_name,
            'parents': [parent_1_name, parent_2_name]
        })
        print(f"[Event {event_id}] Admixture: {child_name} -> {parent_1_name} & {parent_2_name}")

    # ==========================================
    # PHASE 3: PARAMETERS
    # ==========================================

    def set_node_ne(self, node_name, ne):
        """Set Effective Population Size for a node."""
        self.get_node(node_name).ne = ne

    def set_merge_time(self, parent_name, time):
        """
        Set parameters for a Merge event.
        - parent_name: The node created by the merge.
        - time: When the merge happened (generations ago).
        """
        parent = self.get_node(parent_name)
        parent.time_start = time
        
        # The children end exactly when the parent starts
        for child in parent.descendants:
            child.time_end = time

    def set_admixture_parameters(self, child_name, time, fraction_parent_1, parent_1_name):
        """
        Set parameters for an Admixture event.
        - child_name: The node that splits backwards.
        - time: When the admixture happened.
        - fraction_parent_1: Fraction of ancestry from Parent 1.
        - parent_1_name: Name of Parent 1 (to distinguish which fraction is which).
        """
        child = self.get_node(child_name)
        child.time_end = time
        
        fraction_parent_2 = 1.0 - fraction_parent_1
        
        # Identify the parents from the topology
        p1_found = False
        for ancestor in child.ancestors:
            ancestor.time_start = time # Parents start when child splits
            
            if ancestor.name == parent_1_name:
                child.admixture_fractions[ancestor.name] = fraction_parent_1
                p1_found = True
            else:
                # Assume the other parent gets the remainder
                child.admixture_fractions[ancestor.name] = fraction_parent_2
        
        if not p1_found:
             raise ValueError(f"Parent {parent_1_name} is not an ancestor of {child_name}")

    def finalize_root(self):
        """
        Finds the unique root of the topology and sets its end time to infinity.
        """
        roots = [n for n in self.nodes.values() if n.is_root]
        
        # We don't raise error here if multiple roots exist, as is_valid handles that check.
        # We just set inf for any node that currently looks like a root.
        for r in roots:
             r.time_end = float('inf')
             print(f"[Finalize] Set {r.name} time_end to infinity.")

    # ==========================================
    # UTILITIES & VALIDATION
    # ==========================================

    def is_valid(self):
        """
        Validates the topology by simulating the events in order.
        Checks:
        1. Connectivity: Leaves must reduce to exactly one Root.
        2. Event Logic: Inputs to events must be active available lineages.
        3. Time Consistency: Parents must start >= Children end.
        """
        print("\n--- Validating Topology ---")
        
        # Simulation State: Set of active lineages (nodes currently existing at this step)
        active_lineages = set(self.initial_leaves)
        
        for event in self.ordered_events:
            e_id = event['id']
            e_type = event['type']
            
            if e_type == 'MERGE':
                c1, c2 = event['children']
                p = event['parent']
                
                if c1 not in active_lineages or c2 not in active_lineages:
                    raise ValueError(f"Event {e_id} Invalid: Children {c1} or {c2} are not active lineages at this step.")
                
                active_lineages.remove(c1)
                active_lineages.remove(c2)
                active_lineages.add(p)
                
            elif e_type == 'ADMIXTURE':
                c = event['child']
                p1, p2 = event['parents']
                
                if c not in active_lineages:
                     raise ValueError(f"Event {e_id} Invalid: Child {c} is not an active lineage at this step.")
                
                active_lineages.remove(c)
                active_lineages.add(p1)
                active_lineages.add(p2)

        # 1. Check Connectivity
        if len(active_lineages) != 1:
            raise ValueError(f"Invalid Topology: Ends with {len(active_lineages)} roots ({active_lineages}). Must be exactly 1.")
        
        root_name = list(active_lineages)[0]
        print(f"✓ Connectivity Check Passed. Final Root: {root_name}")

        # 2. Check Time Consistency (if times are set)
        for node in self.nodes.values():
            if node.time_start is not None and node.time_end is not None:
                if node.time_start > node.time_end:
                     # Exception: If end is inf, it's fine
                     if node.time_end != float('inf'):
                        raise ValueError(f"Time Error in {node.name}: Starts at {node.time_start} but ends earlier at {node.time_end}.")
        
        print("✓ Time Consistency Check Passed.")
        return True

    def get_topology_matrix_representation(self):
        """
        Returns:
        1. migration_matrices: List of NxN matrices (one per event).
        2. event_types: List of strings ('MERGE' or 'ADMIXTURE').
        3. admixture_map: Dict mapping event index to admixture details.
        
        Matrix Logic:
        - Rows: Source nodes (Children)
        - Cols: Dest nodes (Parents)
        - If Merge: M[child][parent] = 1, M[other][other] = 1 (Identity for others)
        - If Admixture: Identity Matrix (as requested)
        """
        # 1. Map all nodes to indices
        # We sort to ensure deterministic ordering
        all_node_names = sorted(self.nodes.keys())
        node_to_idx = {name: i for i, name in enumerate(all_node_names)}
        n_nodes = len(all_node_names)
        
        migration_matrices = []
        event_types = []
        admixture_map = {}

        for idx, event in enumerate(self.ordered_events):
            e_type = event['type']
            event_types.append(e_type)
            
            # Create base matrix (Identity)
            # Default behavior: Nodes map to themselves (persistence) unless changed
            mat = np.eye(n_nodes)
            
            if e_type == 'MERGE':
                c1 = event['children'][0]
                c2 = event['children'][1]
                p = event['parent']
                
                i_c1 = node_to_idx[c1]
                i_c2 = node_to_idx[c2]
                i_p = node_to_idx[p]
                
                # Zero out the children self-mapping (they are merging, not persisting as themselves)
                mat[i_c1, i_c1] = 0
                mat[i_c2, i_c2] = 0
                
                # Map children to parent
                mat[i_c1, i_p] = 1
                mat[i_c2, i_p] = 1
                
                # Note: The parent 'p' implicitly starts existing here. 
                # If 'p' was not active before, M[p][p] being 1 is fine (0->0 transition effectively).
                
            elif e_type == 'ADMIXTURE':
                # User Request: "If there's an admixture event, simply put an identity matrix there."
                # We store the details in the map instead.
                admixture_map[idx] = {
                    'child': event['child'],
                    'parents': event['parents']
                }
                # Matrix remains Identity
            
            migration_matrices.append(mat)

        return migration_matrices, event_types, admixture_map

    def print_summary(self):
        print(f"\n{'Node':<15} {'Time Range':<15} {'Ne':<8} {'Ancestors (Frac)':<30} {'Descendants'}")
        print("-" * 90)
        
        sorted_nodes = sorted(self.nodes.values(), key=lambda x: x.time_start if x.time_start is not None else -1)
        
        for n in sorted_nodes:
            t_start = str(n.time_start) if n.time_start is not None else "?"
            t_end = str(n.time_end) if n.time_end is not None else "Past"
            
            anc_str = ""
            if n.ancestors:
                parts = []
                for anc in n.ancestors:
                    frac = n.admixture_fractions.get(anc.name)
                    frac_str = f"{frac:.2f}" if frac is not None else "?"
                    parts.append(f"{anc.name}({frac_str})")
                anc_str = ", ".join(parts)
            else:
                anc_str = "ROOT"

            children = ", ".join([c.name for c in n.descendants]) or "LEAF"
            
            print(f"{n.name:<15} {t_start}-{t_end:<9} {str(n.ne):<8} {anc_str:<30} {children}")