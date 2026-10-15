"""Correctness of the relaxed-graph models.

    python -B -m unittest test_relaxed -v
"""
import logging
import sys
import unittest

import numpy as np

import data as dt
from models import compile_model, at
from skeleton import Skeleton

sys.path.insert(0, str(dt.ROOT))
from methods import stan_data as sd
from methods.models import MODELS

logging.getLogger("cmdstanpy").setLevel(logging.WARNING)
C = dt.load_config()
SK = Skeleton(C["leaves"])
THETA = dict(times=[20.0, 35.0, 60.0, 200.0], f=0.7, N=15000.0)


def synthetic_obs(seed=0):
    """Any data will do: the tests compare two implementations of one likelihood."""
    r = np.random.default_rng(seed)
    L, B = SK.L, len(C["bin_edges_cm"]) - 1
    counts = r.poisson(30, size=(B, L, L)); counts = np.triu(counts) + np.triu(counts, 1).transpose(0, 2, 1)
    w = r.normal(0, 0.01, (L, L)); w = (w + w.T) / 2
    return dict(ibd_count=counts, w_hat=w, w_se=np.full((L, L), 0.004), cm=500.0,
                ibd_hat=np.zeros((B, L, L)), ibd_se=np.ones((B, L, L)))


def vertex_data(base, vertices):
    return dict(base, V=len(vertices), vertex_admix=[v[0] + 1 for v in vertices],
                vertex_merge=[[p + 1 for p in v[1:]] for v in vertices], log_prior=[0.0] * len(vertices))


class RelaxedGraphTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.obs = synthetic_obs()
        cls.base = dt.stan_data(C, SK, cls.obs)
        cls.discrete = compile_model("discrete")
        cls.estep = compile_model("relaxed_estep")
        cls.standard = __import__("cmdstanpy").CmdStanModel(stan_file=str(MODELS["mixed_Nfixed_poisson"].path))

    def test_skeleton_counts(self):
        from collections import Counter
        g = SK.graphs()
        self.assertEqual((len(SK.vertices), len(g)), (108, 45))
        # (binary?, vertices per graph): a binary graph has 1 or 2 vertices (its mirror
        # with the admixture sources swapped); a three-way root 3 or 6
        self.assertEqual(Counter((SK.is_binary(SK.vertices[ix[0]]), len(ix)) for ix in g.values()),
                         {(True, 1): 12, (True, 2): 21, (False, 3): 6, (False, 6): 6})

    def test_graph_keys_do_not_depend_on_the_process(self):
        """Graph identity must not depend on set iteration order (string hashing is
        randomised per process)."""
        import subprocess
        code = "from skeleton import Skeleton; print(sorted(Skeleton(['a','b','c']).graphs()))"
        runs = {subprocess.run([sys.executable, "-B", "-c", code], capture_output=True, text=True, check=True,
                               env={"PYTHONHASHSEED": str(k), "PATH": ""}).stdout for k in range(4)}
        self.assertEqual(len(runs), 1)

    def test_relaxed_matrices_are_row_stochastic(self):
        r = np.random.default_rng(1)
        for _ in range(20):
            wa = r.dirichlet(np.ones(SK.L)); wm = [r.dirichlet(np.ones(len(SK.pairs))) for _ in range(SK.n_merge)]
            for M in SK.matrices(wa, wm, r.uniform()):
                self.assertTrue(np.all(M >= -1e-15)); np.testing.assert_allclose(M.sum(1), 1)

    def test_vertices_equal_the_standard_model(self):
        """Deme representation = node representation of methods' mixed_Nfixed_poisson,
        for every binary vertex: same IBD and SNP log-likelihood at the same theta."""
        binary = [v for v in SK.vertices if SK.is_binary(v)]
        ours = at(self.discrete, vertex_data(self.base, binary), THETA)["loglik"]
        spec = MODELS["mixed_Nfixed_poisson"]
        e = np.asarray(C["bin_edges_cm"])
        for v, mine in zip(binary, ours):
            d = sd.stan_data(SK.topology(v), self.obs, spec, edges=e, n_haploid=dt.n_haploid(C), T_max=C["T_max"])
            ref = at(self.standard, d, dict(times=THETA["times"], admixture_fractions=[THETA["f"]],
                                            effective_N=THETA["N"]))
            self.assertAlmostEqual(mine, ref["lp_ibd"] + ref["lp_snp"], delta=1e-6 * abs(mine), msg=SK.label(v))

    def test_relaxed_at_a_vertex_equals_discrete(self):
        lam = 0.5
        for v in (SK.vertices[0], SK.vertices[40], SK.vertices[107]):
            y_a = np.zeros(SK.L); y_a[v[0]] = 40 * lam
            y_m = []
            for p in v[1:]:
                y = np.zeros(len(SK.pairs)); y[p] = 40 * lam; y_m.append(y.tolist())
            relaxed = at(self.estep, dict(self.base, **{"lambda": lam, "times": THETA["times"],
                                                                     "f": THETA["f"], "N": THETA["N"]}),
                         dict(y_admix=y_a.tolist(), y_merge=y_m))["loglik"]
            discrete = at(self.discrete, vertex_data(self.base, [v]), THETA)["loglik"][0]
            self.assertAlmostEqual(relaxed, discrete, delta=1e-8 * abs(discrete), msg=SK.label(v))

    def test_discrete_posterior_is_the_vertex_softmax(self):
        out = at(self.discrete, vertex_data(self.base, SK.vertices), THETA)
        np.testing.assert_allclose(out["posterior"].sum(), 1)
        ll = out["loglik"]
        np.testing.assert_allclose(out["posterior"], np.exp(ll - ll.max()) / np.exp(ll - ll.max()).sum(), rtol=1e-9, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
