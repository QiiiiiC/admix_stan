import json
from pathlib import Path
import unittest

import numpy as np

import study as s

HERE = s.HERE
VARYING = s.ROOT / "simulation_study" / "four_leaf_b_admixture_varying_ne"


class DesignTests(unittest.TestCase):
    def setUp(self):
        self.c = s.load_config()
        self.v = s.load_config(VARYING / "pipeline_config.json")

    def test_growth_rule(self):
        ne = s.branch_ne(self.c)
        self.assertEqual(ne["root"], 10000.0)
        self.assertAlmostEqual(ne["a"], 150000 ** (1 - 30 / 400) * 10000 ** (30 / 400), places=3)
        for leaf, parent in (("a", "ab"), ("ab", "root"), ("d", "cbd"), ("cbd", "root"), ("c", "cb"), ("cb", "cbd")):
            self.assertGreater(ne[leaf], ne[parent], f"{leaf} should be larger than {parent} (growth)")
        leaves = [ne[x] for x in s.POPS]
        self.assertEqual(len(set(round(x) for x in leaves)), 4)
        self.assertGreater(min(leaves), 20000)

    def test_msprime_sizes_equal_truth(self):
        for c in (self.c, self.v):
            d = s.msprime_demography(c)
            truth = {p["name"][3:]: p["truth"] for p in s.parameters(c) if p["name"].startswith("Ne_")}
            self.assertEqual({p.name: 2 * p.initial_size for p in d.populations}, truth)
            self.assertTrue(all(abs(truth[n] - s.branch_ne(c)[n]) <= 2 for n in s.NODES))
            self.assertEqual(sorted(e.time for e in d.events), [c["times"][k] for k in s.EVENTS])

    def test_varying_design_blocks_are_its_own_simulation(self):
        """The varying design's saved trees and caches carry the seeds this pipeline
        would use (block_seed(seed, 0, k)), so they are its own simulation."""
        from methods.simulate import block_seed
        pool = s.ROOT / self.v["trees_from"]
        caches = sorted((VARYING / "runs" / "pipeline" / "pool_000").glob("block_*.npz"))
        self.assertEqual(len(caches), self.v["blocks"])
        for k, cache in enumerate(caches):
            self.assertEqual(int(np.load(cache)["seed"]), block_seed(self.v["seed"], 0, k))
            self.assertTrue((pool / f"block_{k:03d}.trees").exists())

    def test_hapibd_matches_truth_above_its_detection_limit(self):
        """hap-IBD on one simulated block: the same segments as the genealogy from 2 cM,
        and not more than the truth in the lowest bin (it only misses short ones)."""
        import tempfile
        import tskit
        from methods import ibd
        c = self.c
        ts = tskit.load(HERE / "runs" / "default" / "pool_000" / "block_000.trees")
        with tempfile.TemporaryDirectory() as tmp:
            hap, _ = s.hapibd(ts, c, tmp, ibd.hapibd_command())
        true, _ = s.true_ibd(ts, c)
        above = s.edges(c)[:-1] >= 2.0
        np.testing.assert_array_equal(hap[above], true[above])
        self.assertLessEqual(hap[0].sum(), true[0].sum())
        self.assertGreater(hap[0].sum(), 0.7 * true[0].sum())

    def test_tasks_and_data(self):
        c = self.c
        for cfg in (c, self.v):
            self.assertEqual(cfg["sources"], ["true", "hapibd"])
            self.assertEqual(len(s.fit_tasks(cfg, cfg["sources"], list(s.MODEL_NAMES), list(s.GRAPHS))), 8 * 100 * 2 * 2 * 2)
            self.assertEqual({t[2] for t in s.nuts_tasks(cfg, list(s.GRAPHS))}, {"true"})
        self.assertEqual(len(s.nuts_tasks(c, list(s.GRAPHS))), 8 * 10 * 2)
        self.assertEqual(len(c["bin_edges_cm"]) - 1, 22)
        nb, L = 22, 4
        obs = dict(ibd_hat=np.zeros((nb, L, L)), ibd_se=np.ones((nb, L, L)), ibd_count=np.ones((nb, L, L), int),
                   w_hat=np.eye(L) * .01, w_se=np.ones((L, L)) * .01, cm=100.0)
        for model in s.MODEL_NAMES:
            for graph, build in s.GRAPHS.items():
                data = s.stan_data(graph, obs, c, model)
                self.assertEqual(data["n_events"], len(build().ordered_events))
                self.assertEqual(data["n_bins"], nb)
                self.assertEqual("w_hat" in data, model == "mixed")


if __name__ == "__main__":
    unittest.main()
