from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import numpy as np
import msprime
import subprocess
import tskit
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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
        self.n_admix = 0
        
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

        # Record Event (1-indexed to match MERGE numbering)
        event_id = len(self.ordered_events) + 1
        self.ordered_events.append({
            'id': event_id,
            'type': 'ADMIXTURE',
            'child': child_name,
            'parents': [parent_1_name, parent_2_name]
        })

        self.n_admix += 1
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

    def plot_demography(self, scale=True, ax=None):
        """
        Plots the demographic topology.
        
        Args:
            scale (bool): If True, Y-axis is generations (requires parameters).
                          If False, Y-axis is event rank (topology only).
            ax: Optional matplotlib axis.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))

        # =========================================================
        # 1. Determine Y-coordinates (Vertical Placement)
        # =========================================================
        node_y_intervals = {} # Stores [y_start, y_end] for each node

        if scale:
            # --- REAL TIME MODE ---
            # Use the actual time parameters set by the user
            # Collect all time points to determine margins
            all_times = {0.0}
            for n in self.nodes.values():
                if n.time_start is not None: all_times.add(n.time_start)
                if n.time_end is not None and n.time_end != float('inf'): all_times.add(n.time_end)
            
            max_time = max(all_times) if all_times else 0
            top_margin = max_time * 1.1 if max_time > 0 else 1.0

            for name, node in self.nodes.items():
                start = node.time_start if node.time_start is not None else 0.0
                end = node.time_end if node.time_end is not None else 0.0
                if end == float('inf'): end = top_margin
                node_y_intervals[name] = [start, end]
                
            ax.set_ylabel("Time (Generations ago)")

        else:
            # --- TOPOLOGY ONLY MODE ---
            # Ignore time parameters. Use Event Rank (1, 2, 3...)
            # Leaves start at 0.
            
            # Initialize all nodes to start at 0 and end at "undefined"
            for name in self.nodes:
                node_y_intervals[name] = [0.0, None] # [start, end]

            # Iterate through events to build the vertical structure
            # Each event bumps the Y-axis up by 1 unit
            max_rank = len(self.ordered_events) + 1
            
            for i, event in enumerate(self.ordered_events):
                rank = float(i + 1)
                
                if event['type'] == 'MERGE':
                    c1, c2 = event['children']
                    parent = event['parent']
                    
                    # Children end at this rank
                    node_y_intervals[c1][1] = rank
                    node_y_intervals[c2][1] = rank
                    
                    # Parent starts at this rank
                    node_y_intervals[parent][0] = rank
                    
                elif event['type'] == 'ADMIXTURE':
                    child = event['child']
                    p1, p2 = event['parents']
                    
                    # Child ends at this rank
                    node_y_intervals[child][1] = rank
                    
                    # Parents start at this rank
                    node_y_intervals[p1][0] = rank
                    node_y_intervals[p2][0] = rank

            # Close off any nodes that didn't end (Roots)
            for name in node_y_intervals:
                if node_y_intervals[name][1] is None:
                    node_y_intervals[name][1] = max_rank
            
            ax.set_ylabel("Event Rank (Unscaled)")


        # =========================================================
        # 2. Determine X-coordinates (Horizontal Placement)
        # =========================================================
        node_x = {}
        
        # Initialize leaves evenly spaced
        for i, name in enumerate(self.initial_leaves):
            node_x[name] = float(i)

        def _find_first_merge_partner(pop, start_idx):
            """Walk the event list from start_idx forward; return the
            partner pop in pop's first MERGE, or None if pop is absorbed
            by another ADMIX first."""
            for ev in self.ordered_events[start_idx:]:
                if ev['type'] == 'MERGE' and pop in ev['children']:
                    others = [c for c in ev['children'] if c != pop]
                    return others[0] if others else None
                if ev['type'] == 'ADMIXTURE' and ev['child'] == pop:
                    return None
            return None

        # Process events to assign X to ancestors. Both rules are
        # MIDPOINT-OF-PARTNERS:
        #   MERGE   parent -> midpoint(child_1, child_2)
        #   ADMIX   parent -> midpoint(admix_child, eventual_merge_partner)
        # An ADMIX output that gets absorbed by another ADMIX (no
        # merge partner) uses a small fallback offset from its child.
        for event_idx, event in enumerate(self.ordered_events):
            if event['type'] == 'MERGE':
                c1, c2 = event['children']
                parent = event['parent']
                if c1 in node_x and c2 in node_x:
                    node_x[parent] = (node_x[c1] + node_x[c2]) / 2.0
                else:
                    node_x[parent] = 0.0

            elif event['type'] == 'ADMIXTURE':
                child = event['child']
                p1, p2 = event['parents']
                if child in node_x:
                    base_x = node_x[child]
                    for p, default_sign in [(p1, -1), (p2, +1)]:
                        partner = _find_first_merge_partner(p, event_idx + 1)
                        if partner is not None and partner in node_x:
                            node_x[p] = (base_x + node_x[partner]) / 2.0
                        else:
                            node_x[p] = base_x + default_sign * 0.3

        # =========================================================
        # 3. Drawing
        # =========================================================
        colors = plt.cm.get_cmap('tab20', len(self.nodes))
        
        # Draw Nodes (Vertical Pipes)
        for i, (name, interval) in enumerate(node_y_intervals.items()):
            if name not in node_x: continue
            
            y_bottom, y_top = interval
            x = node_x[name]
            
            ax.plot([x, x], [y_bottom, y_top], lw=6, color=colors(i), solid_capstyle='butt')
            
            # Label the node (midpoint)
            label_y = (y_bottom + y_top) / 2
            ax.text(x + 0.1, label_y, name, fontsize=9, color='black', ha='left')

        # Draw Connectors (Horizontal Lines)
        for i, event in enumerate(self.ordered_events):
            # Rank for unscaled is simply i+1, for scaled we need to lookup times
            if scale:
                # Need to lookup the specific time based on event type
                if event['type'] == 'MERGE':
                    time_val = self.nodes[event['parent']].time_start
                else:
                    time_val = self.nodes[event['child']].time_end
                if time_val is None: time_val = 0
                y_pos = time_val
            else:
                y_pos = float(i + 1)

            if event['type'] == 'MERGE':
                c1, c2 = event['children']
                x1, x2 = node_x[c1], node_x[c2]
                ax.plot([x1, x2], [y_pos, y_pos], 'k-', lw=2)

            elif event['type'] == 'ADMIXTURE':
                child = event['child']
                p1, p2 = event['parents']
                x1, x2 = node_x[p1], node_x[p2]
                ax.plot([x1, x2], [y_pos, y_pos], 'k--', lw=2)
                
                # Only print fractions if scaled (and thus presumably parameters exist)
                if scale:
                    frac1 = self.nodes[child].admixture_fractions.get(p1)
                    frac2 = self.nodes[child].admixture_fractions.get(p2)
                    
                    l1 = f"{frac1:.2f}" if isinstance(frac1, (int, float)) else "?"
                    l2 = f"{frac2:.2f}" if isinstance(frac2, (int, float)) else "?"

                    ax.text(node_x[p1], y_pos, l1, ha='center', va='bottom', fontsize=8)
                    ax.text(node_x[p2], y_pos, l2, ha='center', va='bottom', fontsize=8)

        ax.set_title("Demographic Topology")
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.set_xticks([])
        
        if not ax:
            plt.show()

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
        all_node_names = self.nodes.keys()
        node_to_idx = {name: i for i, name in enumerate(all_node_names)}
        n_nodes = len(all_node_names)
        
        migration_matrices = []
        event_types = []
        admixture_map = {}
        admixture_map_id = {}

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
                admixture_map_id[idx] = {
                    'child': node_to_idx[event['child']],
                    'parents': [node_to_idx[i] for i in event['parents']]
                }
                # Matrix remains Identity
            
            migration_matrices.append(mat)

        return migration_matrices, event_types, admixture_map, admixture_map_id

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