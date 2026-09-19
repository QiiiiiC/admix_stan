"""Scientific invariants; no long simulations or Stan compilation in unit tests."""
import copy
import gzip
import itertools
from pathlib import Path
import sys
import tempfile
import unittest

import msprime
import numpy as np
import tskit

from study import (load_config, candidates, selections, simulation_demography, quiet,
                   to_msprime_demography, true_ibd, pair_counts, snp_covariance,
                   aggregate, stan_data, replay, smooth_data, hapibd, edges)
from run_study import truth_for, parameter_summary


class StudyTests(unittest.TestCase):
    def setUp(self):
        self.c = load_config()

    def test_design_and_all_event_orders(self):
        xs = candidates()
        self.assertEqual(len(xs), 27)
        self.assertEqual(len(set(x["graph"] for x in xs)), 21)
        correct = [x for x in xs if x["correct"]]
        self.assertEqual(len(correct), 2)
        self.assertEqual(sum(truth_for(x, self.c)["correct_order"] for x in correct), 1)
        self.assertTrue(any(x["admixed"] and x["events"][0]["type"] == "MERGE" for x in xs))
        self.assertEqual(len(edges(self.c))-1, 37)
        for scenario, times in (("hidden_loop", [12, 16, 40, 70, 100, 300]),
                                ("no_loop", [16, 40, 100, 300])):
            d = quiet(simulation_demography, self.c, scenario)
            dem = quiet(to_msprime_demography, d)
            dem.validate()
            self.assertEqual([ev.time for ev in dem.events], times)
            self.assertTrue(all(p.initial_size == 7500 for p in dem.populations))
        for x in xs:
            d = quiet(replay, x["events"])
            sd = smooth_data(d)
            self.assertTrue(all(p == 0 or p > i+1 for i, p in enumerate(sd["ne_parent"])))

    def test_shared_selections_and_exposure(self):
        first = selections(self.c, 0)
        self.assertEqual(first, selections(self.c, 0))
        self.assertNotEqual(first, selections(self.c, 1))
        self.assertTrue(all(len(x) == len(set(x)) == 30 and min(x) >= 0 and max(x) < 50 for x in first))
        np.testing.assert_array_equal(np.diag(pair_counts(self.c)), [435]*3)
        self.assertEqual(pair_counts(self.c)[0, 1], 900)

    def test_mrca_ignores_intermediate_path_changes(self):
        # The paths change at 5 cM but every pair has the same MRCA on all 10 cM.
        tab = tskit.TableCollection(sequence_length=10_000_000)
        tab.populations.metadata_schema = tskit.MetadataSchema.permissive_json()
        for name in ("a", "b", "c"):
            tab.populations.add_row(metadata={"name": name})
        for p in range(3):
            tab.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=p)
        for time in (1, 1, 2):
            tab.nodes.add_row(time=time, population=0)
        for left, right, mid in ((0, 5_000_000, 3), (5_000_000, 10_000_000, 4)):
            tab.edges.add_row(left, right, mid, 0)
            tab.edges.add_row(left, right, 5, mid)
            tab.edges.add_row(left, right, 5, 1)
            tab.edges.add_row(left, right, 5, 2)
        tab.sort(); ts = tab.tree_sequence()
        count, length = true_ibd(ts, self.c)
        b = np.searchsorted(edges(self.c), 10.0, side="left")-1
        self.assertEqual(count[b, 0, 1], 1)
        self.assertEqual(length[b, 0, 1], 10.0)
        self.assertEqual(count[:, 0, 1].sum(), 1)

    def test_fast_mrca_matches_reference_scan(self):
        c = copy.deepcopy(self.c)
        c["bin_edges_cm"] = [0.01, 1.0, 0.01]
        d = quiet(to_msprime_demography, quiet(simulation_demography, c, "hidden_loop"))
        ts = msprime.sim_ancestry(samples={p: 2 for p in ("a", "b", "c")}, demography=d,
                                 sequence_length=1e6, recombination_rate=1e-8, random_seed=17)
        count, length = true_ibd(ts, c)
        expected_c = np.zeros_like(count); expected_l = np.zeros_like(length)
        e = edges(c)
        for u, v in itertools.combinations(ts.samples(), 2):
            i, j = sorted((ts.node(u).population, ts.node(v).population))
            old, start = None, 0
            segments = []
            for tree in ts.trees():
                m = tree.mrca(u, v)
                if old is not None and m != old:
                    segments.append(tree.interval.left-start)
                    start = tree.interval.left
                old = m
            segments.append(ts.sequence_length-start)
            for span in segments:
                size = span*1e-6; b = np.searchsorted(e, size, side="left")-1
                if 0 <= b < len(count):
                    expected_c[b, i, j] += 1; expected_l[b, i, j] += size
        for i in range(3):
            for j in range(i, 3):
                np.testing.assert_array_equal(count[:, i, j], expected_c[:, i, j])
                np.testing.assert_allclose(length[:, i, j], expected_l[:, i, j], atol=1e-12)

    def test_aggregation_keeps_raw_counts_and_pair_normalization(self):
        c = self.c
        shape = (37, 3, 3)
        blocks = []
        for i in range(50):
            blocks.append(dict(true_count=np.ones(shape, np.int64)*(i+1),
                true_length=np.ones(shape)*(i+1)*3,
                snp_numer=np.repeat((np.eye(3)*(i+1))[None], 5, axis=0),
                snp_denom=np.ones(5)*(i+2)))
        chosen = selections(c, 0)[0]
        obs = aggregate(blocks, chosen, "true", c)
        total = sum(i+1 for i in chosen)
        self.assertEqual(obs["cm"], 1500)
        self.assertEqual(obs["ibd_count"][0, 0, 1], total)
        self.assertAlmostEqual(obs["ibd_hat"][0, 0, 1], total*3/(1500*900))
        self.assertAlmostEqual(obs["ibd_hat"][0, 0, 0], total*3/(1500*435))
        # Every candidate receives observations in a,b,c order, including trees.
        for candidate in candidates():
            data = stan_data(candidate, obs, c)
            np.testing.assert_array_equal(data["ibd_count"], obs["ibd_count"])
            self.assertEqual(data["n_samples"], [30]*3)
            self.assertEqual(data["admixture_map"].shape, (data["n_admixture"], 4))

    def test_snp_ratio_and_delete_block_se(self):
        numer = np.asarray([np.eye(3)*x for x in (1, 4, 9)])
        denom = np.asarray([2, 3, 4])
        w, se = snp_covariance(numer, denom, self.c)
        center = np.eye(3)-np.ones((3, 3))/3
        np.testing.assert_allclose(w, numer.sum(0)/denom.sum()-center/30)
        jack = np.asarray([(numer.sum(0)-numer[i])/(denom.sum()-denom[i]) for i in range(3)])
        expected = np.maximum(np.sqrt(2/3*np.sum((jack-jack.mean(0))**2, axis=0)), 1e-8)
        np.testing.assert_allclose(se, expected)

    def test_hapibd_includes_hbd(self):
        # Fake caller tests file parsing/exposure, including a within-individual segment.
        d = quiet(to_msprime_demography, quiet(simulation_demography, self.c, "no_loop"))
        ts = msprime.sim_ancestry(samples={p: 2 for p in ("a", "b", "c")}, demography=d,
                                 sequence_length=10e6, recombination_rate=0, random_seed=19)
        ts = msprime.sim_mutations(ts, rate=1e-9, random_seed=31)
        with tempfile.TemporaryDirectory() as td:
            script = Path(td)/"caller.py"
            script.write_text("import gzip,sys\n"
                "out=next(x[4:] for x in sys.argv if x.startswith('out='))\n"
                "with gzip.open(out+'.ibd.gz','wt') as f: f.write('sample_0 1 sample_2 1 1 1 3000001 3.0\\n')\n"
                "with gzip.open(out+'.hbd.gz','wt') as f: f.write('sample_0 1 sample_0 2 1 1 2500001 2.5\\n')\n")
            count, lengths = hapibd(ts, self.c, Path(td)/"output", [sys.executable, str(script)])
            self.assertEqual(count[:, 0, 0].sum(), 1)
            self.assertEqual(count[:, 0, 1].sum(), 1)
            self.assertEqual(lengths[:, 0, 0].sum(), 2.5)
            self.assertEqual(lengths[:, 0, 1].sum(), 3.)

    def test_weighted_intervals(self):
        result = parameter_summary(np.asarray([[1., 2.], [3., 4.]]), np.asarray([0.25, 0.75]))
        np.testing.assert_allclose(result["mean"], [2.5, 3.5])
        empty = parameter_summary(np.empty((2, 0)), np.asarray([0.25, 0.75]))
        self.assertEqual(empty["mean"], [])


if __name__ == "__main__":
    unittest.main()
