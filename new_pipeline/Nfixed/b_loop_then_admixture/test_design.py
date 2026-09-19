import unittest
import numpy as np
import study as s
import run_study as runner


class DesignTests(unittest.TestCase):
    def setUp(self):
        self.c=s.load_config()
        self.candidates=s.quiet(s.candidates)

    def test_twenty_subsamples_not_twenty_individuals(self):
        self.assertEqual(self.c["replicates"],20)
        self.assertEqual(self.c["diploid_samples"],15)
        samples=s.selections(self.c,0)
        self.assertEqual(len(samples),20)
        self.assertTrue(all(len(set(x))==30 for x in samples))
        self.assertEqual(samples,s.selections(self.c,0))

    def test_generating_events_and_sizes(self):
        for scenario,times in (("b_loop",[2,15,20,40,100,300]),("no_loop",[20,40,100,300])):
            d=s.quiet(s.simulation_demography,self.c,scenario)
            ms=s.quiet(s.to_msprime_demography,d); ms.validate()
            self.assertEqual([e.time for e in ms.events],times)
            self.assertTrue(all(p.initial_size==7500 for p in ms.populations))
            self.assertFalse(any(e.get("child")=="c" for e in d.ordered_events))

    def test_whole_search_plus_true_loop(self):
        self.assertEqual(len(self.candidates),29)
        self.assertEqual(len({x["graph"] for x in self.candidates}),22)
        self.assertEqual(sum(x["explicit_loop"] for x in self.candidates),2)
        for x in self.candidates:
            d=s.quiet(s.replay,x["events"]); data=s.smooth_data(d)
            self.assertTrue(all(p==0 or p>i+1 for i,p in enumerate(data["ne_parent"])))
        self.assertEqual(list(s.VARIANTS),["poisson_shared","poisson_separate"])

    def test_loop_and_external_admixture_are_distinct(self):
        x=next(x for x in self.candidates if x["id"]=="g22_o1")
        mapping={x["name"]:x for x in s.event_parameters(x,self.c)}
        self.assertEqual(mapping["loop_major_fraction"]["index"],0)
        self.assertEqual(mapping["b_fraction"]["index"],1)
        self.assertEqual(mapping["b_split"]["index"],2)
        self.assertTrue(mapping["loop_major_fraction"]["fold"])
        self.assertFalse(mapping["b_fraction"]["fold"])
        self.assertTrue(runner.truth_for(x,self.c)["correct_order"])
        other=next(x for x in self.candidates if x["id"]=="g22_o2")
        self.assertFalse(runner.truth_for(other,self.c)["correct_order"])

    def test_semantic_summaries_and_residuals(self):
        x=next(x for x in self.candidates if x["id"]=="g22_o1")
        count=np.ones((37,3,3),int)*2
        obs=dict(ibd_hat=np.ones((37,3,3))*.001,ibd_se=np.ones((37,3,3))*.0001,
                 ibd_count=count,w_hat=np.eye(3)*.01,w_se=np.ones((3,3))*.01,cm=1500)
        sv=dict(cumulative_times=np.array([[2,15,20,40,100,300]]*2),
                admixture_fractions=np.array([[.25,.6],[.75,.8]]),
                Ne=np.ones((2,len(x["nodes"]))) *15000,
                ibd_fraction=np.ones((2,37,3,3))*.001,
                ibd_number=np.repeat((count/(1500*s.pair_counts(self.c)))[None],2,axis=0),
                W_centered=np.repeat(obs["w_hat"][None],2,axis=0))
        r=runner.summaries(sv,np.array([.5,.5]),obs,x,self.c)
        self.assertAlmostEqual(r["semantic_parameters"]["b_fraction"]["mean"],.7)
        self.assertAlmostEqual(r["semantic_parameters"]["loop_major_fraction"]["mean"],.75)
        np.testing.assert_allclose(r["residuals"]["ibd_pearson"],0)
        np.testing.assert_allclose(r["residuals"]["snp_raw"],0)


if __name__=="__main__":
    unittest.main()
