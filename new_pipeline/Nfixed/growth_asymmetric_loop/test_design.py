import unittest
import numpy as np
import study as s
import run_study as runner


class DesignTests(unittest.TestCase):
    def setUp(self):
        self.c = s.load_config()
        self.candidates = s.quiet(s.candidates)

    def test_generating_sizes_in_msprime(self):
        """Sizes at branch ends, from the DemographyDebugger, both scenarios."""
        t = self.c["times"]; n0 = self.c["haploid_ne"]; g = self.c["ne"]
        for scenario in s.SCENARIOS:
            d = s.msprime_demography(self.c, scenario)
            dbg = d.debug()
            names = [p.name for p in d.populations]
            def haploid(name, time):
                return 2*dbg.population_size_trajectory(np.array([time]))[0][names.index(name)]
            grow = scenario == "growth"
            self.assertAlmostEqual(haploid("a", 0), n0*(g["growth_fold_a"] if grow else 1), delta=1)
            self.assertAlmostEqual(haploid("a", t["left"]-1e-6), n0, delta=n0*0.01)
            self.assertAlmostEqual(haploid("b", 0), n0*(g["growth_fold_b"] if grow else 1), delta=1)
            self.assertAlmostEqual(haploid("b", t["loop_open"]-1e-6), n0, delta=n0*0.01)
            mid = (t["loop_open"]+t["loop_close"])/2
            self.assertAlmostEqual(haploid("loop1", mid), n0, delta=1)
            self.assertAlmostEqual(haploid("loop2", mid), n0/(g["ratio_loop"] if grow else 1), delta=1)
            self.assertAlmostEqual(haploid("b1", t["b_split"]+1), n0, delta=1)
            self.assertAlmostEqual(haploid("b2", t["b_split"]+1), n0/(g["ratio_split"] if grow else 1), delta=1)
            self.assertAlmostEqual(haploid("root", t["root"]+50), n0, delta=1)
            self.assertEqual(sorted(e.time for e in d.events if not hasattr(e, "growth_rate")),
                             [t[k] for k in ("loop_open", "loop_close", "b_split", "left", "right", "root")])

    def test_effective_ne_composite_and_references(self):
        c = self.c; t = c["times"]; f = c["fractions"]["loop"]; n0 = c["haploid_ne"]
        # loop interval: two lineages coalesce only if they fall in the same branch
        loop_t = (t["loop_open"]+t["loop_close"])/2
        for scenario, n2 in (("growth", n0/c["ne"]["ratio_loop"]), ("constant", n0)):
            expect = 1/(f**2/n0 + (1-f)**2/n2)
            self.assertAlmostEqual(float(s.effective_ne(c, scenario, "b", loop_t)), expect, places=6)
        self.assertAlmostEqual(float(s.effective_ne(c, "growth", "b", loop_t)), 20000.0, places=6)
        self.assertAlmostEqual(float(s.effective_ne(c, "constant", "b", loop_t)), 24000.0, places=6)
        # harmonic mean of an exponential branch: N0 ln G / (1 - 1/G)
        G = c["ne"]["growth_fold_a"]
        ref = s.reference(c, "growth", "a")
        self.assertAlmostEqual(ref["harmonic"], n0*np.log(G)/(1-1/G), delta=n0*1e-3)
        self.assertAlmostEqual(ref["young"], n0*G); self.assertAlmostEqual(ref["old"], n0)
        self.assertEqual((ref["t0"], ref["t1"]), (0, t["left"]))
        const = s.reference(c, "constant", "a")
        self.assertTrue(all(abs(const[k]-n0) < 1e-6 for k in ("young", "old", "harmonic", "arithmetic")))
        # explicit-loop naming: b is the leaf only
        leaf = s.reference(c, "growth", "b", collapsed=False)
        self.assertEqual((leaf["t0"], leaf["t1"]), (0, t["loop_open"]))
        self.assertEqual(s.reference(c, "growth", "b", collapsed=True)["t1"], t["b_split"])
        self.assertIsNone(s.reference(c, "growth", "root")["t1"])

    def test_search_and_branch_map(self):
        self.assertEqual(len(self.candidates), 29)
        self.assertEqual(len({x["graph"] for x in self.candidates}), 22)
        self.assertEqual(sum(x["explicit_loop"] for x in self.candidates), 2)
        self.assertEqual(sum(x["correct"] for x in self.candidates), 4)
        for x in self.candidates:
            d = s.quiet(s.replay, x["events"]); data = s.smooth_data(d)
            self.assertTrue(all(p == 0 or p > i+1 for i, p in enumerate(data["ne_parent"])))
            bm = s.branch_map(x)
            self.assertEqual(bool(bm), x["correct"])
        g15 = next(x for x in self.candidates if x["id"] == "g15_o1")
        self.assertEqual(sorted(s.branch_map(g15).values()), sorted(s.BACKBONE))
        g22 = next(x for x in self.candidates if x["id"] == "g22_o1")
        self.assertEqual(sorted(s.branch_map(g22).values()),
                         sorted(["a", "b", "c", "loop1", "loop2", "anc_b", "b.1", "b.2", "left", "right", "root"]))
        mapping = {m["name"]: m for m in s.event_parameters(g22, self.c)}
        self.assertEqual(mapping["loop_major_fraction"]["index"], 0)
        self.assertEqual(mapping["b_fraction"]["index"], 1)
        self.assertTrue(runner.truth_for(g22, self.c)["correct_order"])
        slices = [[x["id"] for i, x in enumerate(self.candidates) if i % 4 == k] for k in range(4)]
        self.assertEqual(sorted(sum(slices, [])), sorted(x["id"] for x in self.candidates))

    def test_residual_tilt(self):
        flat = np.zeros((37, 3, 3)); self.assertEqual(s.residual_tilt(flat, self.c)["0,0"]["slope"], 0.0)
        rising = np.zeros((37, 3, 3)); rising[:, 1, 1] = np.linspace(-2, 2, 37)
        tilt = s.residual_tilt(rising, self.c)
        self.assertAlmostEqual(tilt["1,1"]["slope"], 4/18.0, places=6)   # 4 units over 18 cM
        self.assertAlmostEqual(tilt["1,1"]["correlation"], 1.0, places=6)
        self.assertEqual(tilt["0,0"]["slope"], 0.0)

    def test_summaries_carry_references(self):
        x = next(x for x in self.candidates if x["id"] == "g15_o1")
        count = np.ones((37, 3, 3), int)*2
        obs = dict(ibd_hat=np.ones((37, 3, 3))*.001, ibd_se=np.ones((37, 3, 3))*.0001,
                   ibd_count=count, w_hat=np.eye(3)*.01, w_se=np.ones((3, 3))*.01, cm=1500)
        n = len(x["nodes"])
        sv = dict(cumulative_times=np.array([[30, 50, 100, 300]]*2, float),
                  admixture_fractions=np.array([[.7], [.8]]),
                  Ne_ibd=np.ones((2, n))*20000, Ne_snp=np.ones((2, n))*10000,
                  ibd_fraction=np.ones((2, 37, 3, 3))*.001,
                  ibd_number=np.repeat((count/(1500*s.pair_counts(self.c)))[None], 2, axis=0),
                  W_centered=np.repeat(obs["w_hat"][None], 2, axis=0))
        r = runner.summaries(sv, np.array([.5, .5]), obs, x, self.c, "growth")
        np.testing.assert_allclose(r["residuals"]["ibd_pearson"], 0)
        self.assertEqual([q["node"] for q in r["ne_reference"]], x["nodes"])
        by = {q["branch"]: q for q in r["ne_reference"]}
        self.assertAlmostEqual(by["a"]["young"], 150000.0)
        self.assertEqual(by["b"]["t1"], 30)
        self.assertAlmostEqual(r["log_ne_ibd_over_snp"]["mean"][0], np.log(2), places=6)
        self.assertEqual(r["residual_tilt"]["0,0"]["slope"], 0.0)


if __name__ == "__main__":
    unittest.main()
