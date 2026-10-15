"""Checks of the graph space and its moves.

    python -B -m unittest test_graphs -v
"""
import random
import unittest

import graphs as gr
from methods.demography import DemographicTopology

SPACES = [  # (leaves, space); small enough to enumerate
    ("abc", gr.Space(max_admix=1, bubbles=True, ghost=True)),
    ("abc", gr.Space(max_admix=2)),
    ("abcd", gr.Space(max_admix=1)),
    ("abcd", gr.Space(max_admix=1, targets="leaves")),
]


def four_leaf_truth():
    d = DemographicTopology(list("abcd"))
    d.add_admixture_event("b", "b1", "b2")
    d.add_merge_event("a", "b1", "ab")
    d.add_merge_event("b2", "c", "cb")
    d.add_merge_event("cb", "d", "cbd")
    d.add_merge_event("ab", "cbd", "root")
    return d


class GraphTests(unittest.TestCase):
    def test_tree_counts(self):
        """Rooted binary trees on L leaves: (2L - 3)!! = 3, 15, 105."""
        for leaves, n in (("abc", 3), ("abcd", 15), ("abcde", 105)):
            seen, _ = gr.enumerate_space(gr.Space(max_admix=0), leaves)
            self.assertEqual(len(seen), n)

    def test_moves_are_reversible_and_stay_valid(self):
        """Every neighbour relation is mutual, every graph satisfies the constraints and
        has a valid listing, and the per-branch model accepts it."""
        import topology as tp
        for leaves, sp in SPACES:
            seen, nbrs = gr.enumerate_space(sp, leaves)
            for k, n in nbrs.items():
                for k2 in n:
                    self.assertIn(k, nbrs[k2], (leaves, sp, k, k2))
            for g in seen.values():
                self.assertTrue(gr.allowed(g, sp))
                tp.topology_data(gr.to_topology(g))

    def test_relocate_is_reversible(self):
        for leaves, sp in SPACES[1:3]:
            _, nbrs = gr.enumerate_space(sp, leaves, relocate=True)
            for k, n in nbrs.items():
                for k2 in n:
                    self.assertIn(k, nbrs[k2], (leaves, sp, k, k2))

    def test_key_ignores_names_listing_and_source_sides(self):
        d = four_leaf_truth()
        swapped = DemographicTopology(list("abcd"))       # sources swapped, other names and listing
        swapped.add_admixture_event("b", "x2", "x1")
        swapped.add_merge_event("x2", "c", "q")
        swapped.add_merge_event("q", "d", "r")
        swapped.add_merge_event("a", "x1", "p")
        swapped.add_merge_event("r", "p", "top")
        self.assertEqual(gr.from_topology(d).key(), gr.from_topology(swapped).key())
        self.assertEqual(gr.from_topology(d).key(), gr.from_topology(gr.to_topology(gr.from_topology(d))).key())

    def test_truth_is_in_the_four_leaf_space(self):
        seen, _ = gr.enumerate_space(gr.Space(max_admix=1), "abcd")
        self.assertIn(gr.from_topology(four_leaf_truth()).key(), seen)

    def test_moves_reach_every_graph(self):
        """From any graph, a random walk of moves eventually visits the whole space."""
        seen, nbrs = gr.enumerate_space(gr.Space(max_admix=1), "abc")
        rng = random.Random(0)
        for start in seen:
            k, visited = start, {start}
            for _ in range(2000):
                k = rng.choice(list(nbrs[k])); visited.add(k)
            self.assertEqual(visited, set(seen))


if __name__ == "__main__":
    unittest.main()
