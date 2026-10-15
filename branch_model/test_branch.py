"""The per-branch model against the current event-order model.

    python -B -m unittest test_branch -v
"""
import logging
import unittest

import numpy as np

import topology as tp
from methods import stan_data as sd
from methods.demography import DemographicTopology
from methods.models import MODELS

logging.getLogger("cmdstanpy").setLevel(logging.ERROR)
EDGES = np.array([1.0, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25])
T_MAX = 1e5
CURRENT = MODELS["mixed_Nsmooth_poisson"]

# Graphs as {event: (kind, inputs, outputs)}; any consistent heights define a history.
FOUR_LEAF = dict(A=("admix", ["b"], ["b1", "b2"]), E1=("merge", ["a", "b1"], ["ab"]), J=("merge", ["b2", "c"], ["cb"]),
                 E2=("merge", ["cb", "d"], ["cbd"]), R=("merge", ["ab", "cbd"], ["root"]))
TREE = dict(E1=("merge", ["a", "b"], ["ab"]), E2=("merge", ["c", "d"], ["cd"]), R=("merge", ["ab", "cd"], ["root"]))
LOOP = dict(A1=("admix", ["a"], ["aP1", "aP2"]), M1=("merge", ["aP2", "b"], ["bP"]), M2=("merge", ["bP", "aP1"], ["ab"]),
            A2=("admix", ["c"], ["cP1", "cP2"]), M3=("merge", ["ab", "cP1"], ["left"]),
            M4=("merge", ["d", "cP2"], ["right"]), R=("merge", ["left", "right"], ["root"]))
CASES = [  # (graph, heights): each row puts the unrelated events in a different time order
    (FOUR_LEAF, dict(A=20, E1=60, J=100, E2=200, R=400)),
    (FOUR_LEAF, dict(A=20, E1=150, J=100, E2=200, R=400)),
    (FOUR_LEAF, dict(A=20, E1=300, J=100, E2=200, R=400)),
    (TREE, dict(E1=50, E2=52, R=300)),
    (TREE, dict(E1=52, E2=50, R=300)),
    (LOOP, dict(A1=20, M1=60, M2=100, A2=500, M3=700, M4=900, R=1000)),
    (LOOP, dict(A2=10, A1=20, M4=40, M1=60, M2=100, M3=700, R=1000)),
    (LOOP, dict(A1=20, A2=30, M1=60, M4=80, M2=100, M3=120, R=150)),
]


def build(graph, listing):
    d = DemographicTopology(list("abcd"))
    for e in listing:
        kind, ins, outs = graph[e]
        if kind == "merge":
            d.add_merge_event(*ins, *outs)
        else:
            d.add_admixture_event(*ins, *outs)
    return d


def at(model, data, params):
    fit = model.sample(data=data, inits=params, fixed_param=True, chains=1, iter_sampling=1, seed=1,
                       show_progress=False, show_console=False, sig_figs=17)
    return {k: v[0] for k, v in fit.stan_variables().items()}


def synthetic_obs(seed=0, L=4):
    r = np.random.default_rng(seed); B = len(EDGES) - 1
    counts = r.poisson(30, (B, L, L)); counts = np.triu(counts) + np.triu(counts, 1).transpose(0, 2, 1)
    w = r.normal(0, 0.01, (L, L)); w = (w + w.T) / 2
    return dict(ibd_count=counts, ibd_hat=np.zeros((B, L, L)), ibd_se=np.ones((B, L, L)), w_hat=w,
                w_se=np.full((L, L), 0.004), cm=500.0)


class BranchModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from cmdstanpy import CmdStanModel
        cls.branch = CmdStanModel(stan_file=str(tp.STAN_FILE))
        cls.current = CmdStanModel(stan_file=str(CURRENT.path))
        cls.obs = synthetic_obs()

    def compare(self, graph, heights, listing, seed):
        """Same parameters in both models: the current one built in the time order of
        `heights`, the per-branch one from `listing` (any listing respecting the partial order)."""
        r = np.random.default_rng(seed)
        raw = {x: r.normal(0, 0.6) for x in "abcd"}
        for _, _, outs in graph.values():
            raw.update({x: r.normal(0, 0.6) for x in outs})
        n_admix = sum(kind == "admix" for kind, _, _ in graph.values())
        common = dict(admixture_fractions=r.uniform(0.15, 0.85, n_admix).tolist(), mu_log=float(np.log(12000)),
                      sigma_log=0.3, tau=0.25)
        order = sorted(heights, key=heights.get)
        old = build(graph, order)
        # admixture fractions are indexed in the order the admixtures are added
        fr = dict(zip([e for e in graph if graph[e][0] == "admix"], common["admixture_fractions"]))
        ref = at(self.current, sd.stan_data(old, self.obs, CURRENT, edges=EDGES, n_haploid=[40] * 4, T_max=T_MAX),
                 dict(common, admixture_fractions=[fr[e] for e in order if graph[e][0] == "admix"],
                      times=np.diff([0] + [heights[e] for e in order]).tolist(), Ne_raw=[raw[x] for x in old.nodes]))
        new = build(graph, listing)
        out = at(self.branch, tp.stan_data(new, self.obs, EDGES, [40] * 4, T_MAX),
                 dict(common, admixture_fractions=[fr[e] for e in listing if graph[e][0] == "admix"],
                      times=tp.heights_to_times(new, [heights[e] for e in listing]), Ne_raw=[raw[x] for x in new.nodes]))
        ne_old = dict(zip(old.nodes, ref["Ne"])); ne_new = dict(zip(new.nodes, out["Ne"]))
        np.testing.assert_allclose([ne_new[x] for x in ne_old], list(ne_old.values()), rtol=1e-12)
        np.testing.assert_allclose(out["ibd_number"], ref["ibd_number"], rtol=1e-10)
        np.testing.assert_allclose(out["W_centered"], ref["W_centered"], rtol=1e-9, atol=1e-15)
        self.assertAlmostEqual(out["lp_ibd"], ref["lp_ibd"], delta=1e-7 * abs(ref["lp_ibd"]))
        self.assertAlmostEqual(out["lp_snp"], ref["lp_snp"], delta=1e-7 * abs(ref["lp_snp"]))

    def test_equals_the_current_model_in_every_time_order(self):
        for k, (graph, heights) in enumerate(CASES):
            with self.subTest(case=k):
                self.compare(graph, heights, list(graph), seed=k)

    def test_any_listing_gives_the_same_answer(self):
        """The listing of events is bookkeeping: one that contradicts the time order
        (but respects the partial order) changes nothing."""
        self.compare(TREE, dict(E1=50, E2=52, R=300), ["E2", "E1", "R"], seed=13)
        self.compare(FOUR_LEAF, dict(A=20, E1=60, J=100, E2=200, R=400), ["A", "J", "E2", "E1", "R"], seed=11)
        self.compare(LOOP, dict(A1=20, M1=60, M2=100, A2=500, M3=700, M4=900, R=1000),
                     ["A2", "A1", "M4", "M1", "M2", "M3", "R"], seed=12)

    def test_listing_must_respect_the_partial_order(self):
        with self.assertRaises(Exception):
            tp.topology_data(build(FOUR_LEAF, ["E1", "A", "J", "E2", "R"]))


if __name__ == "__main__":
    unittest.main()
