import unittest

import numpy as np

import study as s
import run_study as runner


class DesignTests(unittest.TestCase):
    def setUp(self):
        self.c = s.load_config()
        self.cands = {x["id"]: x for x in s.candidates()}

    def test_generating_history_in_msprime(self):
        d = s.msprime_demography(self.c)
        t = self.c["times"]
        self.assertEqual(sorted(e.time for e in d.events),
                         [t[k] for k in ("loop_split", "loop_donor_merge", "loop_close", "ancient_admix",
                                         "left", "right", "root")])
        self.assertTrue(all(p.initial_size == self.c["haploid_ne"] / 2 and p.growth_rate == 0
                            for p in d.populations))
        admix = {e.derived: e for e in d.events if hasattr(e, "proportions")}
        np.testing.assert_allclose(admix["a"].proportions, [0.5, 0.5])
        np.testing.assert_allclose(admix["c"].proportions, [0.7, 0.3])
        self.assertEqual(list(admix["c"].ancestral), ["cP1", "cP2"])

    def test_candidates(self):
        true = s.true_topology(self.c)
        self.assertEqual(self.cands["T_true"]["events"], true.ordered_events)
        shape = {k: (len(x["events"]), x["n_admixture"]) for k, x in self.cands.items()}
        self.assertEqual(shape, {"T_true": (7, 2), "T_alt1": (5, 1), "T_alt2": (5, 1)})
        for x in self.cands.values():
            s.replay(x["events"]).is_valid()

    def test_semantic_parameters(self):
        c = self.c
        m = {k: {p["name"]: p for p in s.event_parameters(x, c)} for k, x in self.cands.items()}
        truth = s.true_topology(c)
        self.assertEqual([m["T_true"][n]["truth"] for n in ("loop_split", "loop_donor_merge", "loop_close",
                                                           "ancient_admix", "left", "right", "root")],
                         truth.event_times())
        self.assertEqual((m["T_true"]["loop_fraction"]["index"], m["T_true"]["ancient_fraction"]["index"]), (0, 1))
        self.assertEqual([m["T_true"][n]["truth"] for n in ("loop_fraction", "ancient_fraction")],
                         truth.admixture_fractions())
        self.assertEqual(set(m["T_alt1"]), {"loop_close", "ancient_admix", "ancient_fraction", "left", "right", "root"})
        self.assertEqual(m["T_alt1"]["ancient_fraction"]["index"], 0)
        self.assertEqual(m["T_alt1"]["ancient_admix"]["index"], 1)
        self.assertEqual(set(m["T_alt2"]), {"loop_split", "loop_fraction", "loop_donor_merge", "loop_close"})

    def test_tasks_and_selections(self):
        c = self.c
        work = s.tasks(c, s.SOURCES, list(s.MODEL_NAMES))
        per_rep = (2 * 2 + 1) * 3       # ibd and mixed per IBD source, SNP-only once; 3 candidates
        self.assertEqual(len(work), len(c["genome_cm"]) * c["replicates"] * per_rep)
        self.assertEqual(len(set(work)), len(work))
        self.assertTrue(all(src == "none" for _, _, src, model, _ in work if model == "snp"))
        slices = [[t for i, t in enumerate(work) if i % 8 == k] for k in range(8)]
        self.assertEqual(sorted(sum(slices, [])), sorted(work))
        pilot = s.tasks(c, s.SOURCES, list(s.MODEL_NAMES), replicate_limit=1)
        self.assertEqual(len(pilot), len(c["genome_cm"]) * per_rep)
        seen = set()
        for cm in c["genome_cm"]:
            sel = s.selections(c, 0, cm)
            self.assertEqual(len(sel), c["replicates"])
            self.assertTrue(all(len(x) == round(cm / c["block_cm"]) == len(set(x)) for x in sel))
            self.assertEqual(sel, s.selections(c, 0, cm))
            seen.add(tuple(sel[0]))
        self.assertEqual(len(seen), len(c["genome_cm"]))

    def test_stan_data_and_summaries(self):
        c = self.c; L = 4; nb = len(s.edges(c)) - 1
        obs = dict(ibd_hat=np.full((nb, L, L), 1e-4), ibd_se=np.full((nb, L, L), 1e-5),
                   ibd_count=np.ones((nb, L, L), int), w_hat=np.eye(L) * .01, w_se=np.ones((L, L)) * .01, cm=100.0)
        for x in self.cands.values():
            for model in s.MODEL_NAMES:
                data = s.stan_data(x, obs, c, model)
                self.assertEqual(data["n_events"], len(x["events"]))
                self.assertNotIn("ibd_se", data)   # Poisson IBD: no standard error enters
                self.assertEqual(("ibd_count" in data, "w_hat" in data),
                                 (model != "snp", model != "ibd"))
        x = self.cands["T_true"]
        sv = dict(cumulative_times=np.array([[20, 60, 100, 500, 700, 900, 1000]] * 2, float),
                  admixture_fractions=np.array([[.5, .7], [.5, .7]]), effective_N=np.array([1e4, 1e4]),
                  W_centered=np.repeat(obs["w_hat"][None], 2, axis=0))
        r = runner.summaries(sv, np.array([.5, .5]), obs, x, c)
        errors = {k: v["mean"] - v["truth"] for k, v in r["semantic_parameters"].items()}
        self.assertEqual(len(errors), 10)
        self.assertTrue(all(abs(e) < 1e-9 for e in errors.values()))


if __name__ == "__main__":
    unittest.main()
