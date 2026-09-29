"""Fast invariants of the methods package (no Stan compilation).

    python -B -m unittest methods.test_methods -v      # from the repository root
"""
import re
import unittest

import numpy as np

from methods import blocks, diagnostics, ibd, models, snp, stan_data, topologies
from methods.demography import DemographicTopology, replay, valid_event_orders, clades
from methods.simulate import to_msprime, size_at


def loop_plus_ancient():
    """Four leaves: recent admix-through-b loop on a, ancient admixture on c."""
    d = DemographicTopology(["a", "b", "c", "d"])
    d.add_admixture_event("a", "aP1", "aP2")
    d.add_merge_event("aP2", "b", "bP")
    d.add_merge_event("bP", "aP1", "ab")
    d.add_admixture_event("c", "cP1", "cP2")
    d.add_merge_event("ab", "cP1", "left")
    d.add_merge_event("d", "cP2", "right")
    d.add_merge_event("left", "right", "root")
    d.set_uniform_ne(10000)
    d.set_admixture_parameters("a", 20, 0.5, "aP1")
    d.set_merge_time("bP", 60); d.set_merge_time("ab", 100)
    d.set_admixture_parameters("c", 500, 0.7, "cP1")
    d.set_merge_time("left", 700); d.set_merge_time("right", 900); d.set_merge_time("root", 1000)
    d.finalize_root()
    return d


def data_block_names(path):
    block = models._block(path.read_text(), "data")
    block = re.sub(r"//[^\n]*", "", block)
    return set(re.findall(r"\b(\w+)\s*;", block))


class MethodsTests(unittest.TestCase):
    def test_dsl_truth_and_validity(self):
        d = loop_plus_ancient()
        self.assertTrue(d.is_valid())
        self.assertEqual(d.event_times(), [20, 60, 100, 500, 700, 900, 1000])
        self.assertEqual(d.admixture_fractions(), [0.5, 0.7])
        names = [n for n, _ in diagnostics.true_parameters(d)]
        self.assertEqual(names[:3], ["t_a", "t_bP", "t_ab"])
        self.assertEqual(clades(d.initial_leaves, d.ordered_events)["ab"], {"aP1", "aP2", "b"})

    def test_event_orders(self):
        d = loop_plus_ancient()
        orders = valid_event_orders(d.initial_leaves, d.ordered_events)
        self.assertIn(d.ordered_events, orders)
        for order in orders:
            replay(d.initial_leaves, order).is_valid()
        # the c admixture may sit anywhere before `left`; `right` after it
        self.assertEqual(len(orders), len({tuple(e["id"] for e in o) for o in orders}))

    def test_three_leaf_search(self):
        self.assertEqual(len(topologies.enumerate_all(["a", "b", "c"])), 21)
        cands = topologies.ordered_candidates(["a", "b", "c"])
        self.assertEqual(len(cands), 27)
        self.assertEqual(len({x["graph"] for x in cands}), 21)

    def test_msprime_growth(self):
        d = loop_plus_ancient()
        demo = to_msprime(d, sizes={"a": (100000.0, 10000.0)})
        dbg = demo.debug()
        names = [p.name for p in demo.populations]
        hap = lambda name, t: 2 * dbg.population_size_trajectory(np.array([t]))[0][names.index(name)]
        self.assertAlmostEqual(hap("a", 0), 100000, delta=1)
        self.assertAlmostEqual(hap("a", 20 - 1e-6), 10000, delta=10)
        self.assertAlmostEqual(hap("b", 0), 10000, delta=1)
        self.assertAlmostEqual(float(size_at(100000, 10000, 0, 20, 10)), np.sqrt(1e9), places=3)

    def test_pairs_and_snp_covariance(self):
        np.testing.assert_array_equal(ibd.pair_counts([4, 6]), [[6, 24], [24, 15]])
        rng = np.random.default_rng(0)
        numer = rng.normal(size=(8, 4, 4)); numer = numer + numer.transpose(0, 2, 1)
        w, se = snp.snp_covariance(numer, rng.uniform(1, 2, 8), [40] * 4)
        self.assertEqual(w.shape, (4, 4)); self.assertTrue(np.all(se > 0))
        np.testing.assert_allclose(w, w.T)

    def test_block_selection_and_aggregate(self):
        s = blocks.selections(1, 0, 50, 10, 3)
        self.assertEqual(len(s), 3); self.assertTrue(all(len(x) == 10 and x == sorted(x) for x in s))
        self.assertNotEqual(s, blocks.selections(1, 0, 50, 10, 3, 500))   # a length key draws fresh sets
        L, nb = 4, 5
        fake = [dict(true_count=np.ones((nb, L, L), int), true_length=np.full((nb, L, L), 3.0),
                     snp_numer=np.eye(L)[None].repeat(5, 0) * (k + 1), snp_denom=np.ones(5))
                for k in range(3)]
        obs = blocks.aggregate(fake, [0, 2], "true", 50.0, [10] * L)
        self.assertEqual(obs["cm"], 100.0)
        np.testing.assert_allclose(obs["ibd_hat"][0], 6.0 / (100.0 * ibd.pair_counts([10] * L)))

    def test_normaliser_counts(self):
        expected = {"ibd_Nfixed": 3, "snp_Nfixed": 4, "mixed_Nfixed": 4, "ibd_Nfixed_poisson": 3,
                    "mixed_Nfixed_poisson": 4, "ibd_Nvarying_poisson": 5, "snp_Nvarying": 6,
                    "mixed_Nvarying_poisson": 6, "ibd_Nsmooth": 6, "ibd_Nsmooth_poisson": 6,
                    "snp_Nsmooth": 7, "mixed_Nsmooth_poisson": 7, "mixed_Nsmooth_poisson_separate": 11}
        for name, spec in models.MODELS.items():
            source, n = models.normalized_source(spec)
            self.assertEqual(n, expected[name], name)
            self.assertIn("n_events * 0.01", source)
            self.assertNotRegex(models._block(source, "model"), r"(?m)^\s*\w+(\[[^\]]*\])?\s*~")

    def test_stan_data_matches_every_model(self):
        d = loop_plus_ancient()
        L, e = 4, ibd.bin_edges(1.0, 5.0, 1.0)
        obs = dict(ibd_hat=np.zeros((4, L, L)), ibd_se=np.ones((4, L, L)), ibd_count=np.zeros((4, L, L), int),
                   w_hat=np.zeros((L, L)), w_se=np.ones((L, L)), cm=100.0)
        for name, spec in models.MODELS.items():
            data = stan_data.stan_data(d, obs, spec, edges=e, n_haploid=[40] * L)
            self.assertEqual(set(data), data_block_names(spec.path), name)
            self.assertEqual(np.asarray(data["admixture_map"]).shape, (2, 4))
            if spec.likelihood == "poisson":
                self.assertIn("ibd_count", data)
        self.assertNotIn("ibd_se", stan_data.stan_data(d, obs, models.MODELS["mixed_Nfixed_poisson"],
                                                       edges=e, n_haploid=[40] * L))

    def test_poisson_models_are_derived(self):
        """Every *_poisson IBD/mixed model is exactly what make_poisson.py derives from its source."""
        import importlib.util
        spec = importlib.util.spec_from_file_location("make_poisson", models.STAN / "make_poisson.py")
        mp = importlib.util.module_from_spec(spec); spec.loader.exec_module(mp)
        for name, mixed in mp.SOURCES:
            self.assertEqual((models.STAN / f"{name}_poisson.stan").read_text(),
                             mp.convert((models.STAN / f"{name}.stan").read_text(), mixed), name)

    def test_bridge_recovers_known_normaliser(self):
        """Unnormalised N(0, S) in 3-d: log Z = 0.5 log det(2 pi S) + offset, recovered to its error."""
        from methods.evidence import bridge
        rng = np.random.default_rng(0)
        S = np.array([[1.0, 0.6, 0.0], [0.6, 2.0, 0.3], [0.0, 0.3, 0.5]]); Si = np.linalg.inv(S)
        logp = lambda x: -0.5 * np.einsum("ij,jk,ik->i", x, Si, x) + 7.0
        truth = 0.5 * np.linalg.slogdet(2 * np.pi * S)[1] + 7.0
        post = rng.multivariate_normal(np.zeros(3), S, 4000)
        q_cov = 1.5 * S; q_ci = np.linalg.inv(q_cov)
        logq = lambda x: -0.5 * np.einsum("ij,jk,ik->i", x, q_ci, x) - 0.5 * np.linalg.slogdet(2 * np.pi * q_cov)[1]
        prop = rng.multivariate_normal(np.zeros(3), q_cov, 4000)
        est, err = bridge(logp(post), logq(post), logp(prop), logq(prop))
        self.assertLess(abs(est - truth), 4 * err + 1e-3)
        self.assertLess(err, 0.05)

    def test_trajectory_reference(self):
        G, n0 = 10.0, 15000.0
        ref = diagnostics.trajectory_reference(lambda t: size_at(n0 * G, n0, 0, 50, t), 0, 50)
        self.assertAlmostEqual(ref["harmonic"], n0 * np.log(G) / (1 - 1 / G), delta=n0 * 1e-3)
        self.assertIsNone(diagnostics.trajectory_reference(lambda t: n0 + 0 * t, 300, np.inf)["t1"])


if __name__ == "__main__":
    unittest.main()
