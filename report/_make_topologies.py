"""Generate topology figures for the APM report using the project's
DemographicTopology.plot_demography. Run with the genetics_env python.
Reconstructs (does not import) each experiment's true topology so importing
the experiment scripts does not trigger their pipelines."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# shim: plot_demography uses the removed plt.cm.get_cmap(name, lut)
if not hasattr(plt.cm, "get_cmap"):
    plt.cm.get_cmap = plt.get_cmap

_PARENT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "new_pipeline")
sys.path.insert(0, _PARENT)
from demography import DemographicTopology

OUT = os.path.dirname(os.path.abspath(__file__))


def save(dem, fname, title):
    fig, ax = plt.subplots(figsize=(10, 6.5))
    dem.plot_demography(scale=True, ax=ax)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", fname)


# ---- §1 / §4.1.1 worked example: single recent admixture (Ne uniform 6000 haploid) ----
d1 = DemographicTopology(['a', 'admix', 'b', 'c'])
d1.add_admixture_event('admix', 'admixPa1', 'admixPa2')
d1.add_merge_event('a', 'admixPa1', 'aAdmix')
d1.add_merge_event('b', 'admixPa2', 'bAdmix')
d1.add_merge_event('aAdmix', 'bAdmix', 'subAnc')
d1.add_merge_event('subAnc', 'c', 'anc')
for n in d1.nodes:
    d1.set_node_ne(n, 3000)            # diploid 3000 = haploid 6000
d1.set_admixture_parameters('admix', 20, 0.75, 'admixPa1')
d1.set_merge_time('aAdmix', 40)
d1.set_merge_time('bAdmix', 50)
d1.set_merge_time('subAnc', 200)
d1.set_merge_time('anc', 500)
d1.finalize_root()
save(d1, "topology_example.png", "One-admixture topology (uniform Ne)")


# ---- §4.1.2 two-admixture: recent admix-through-b on a + ancient admix on c ----
d2 = DemographicTopology(['a', 'b', 'c', 'd'])
d2.add_admixture_event('a', 'aP1', 'aP2')
d2.add_merge_event('aP2', 'b', 'bP')
d2.add_merge_event('bP', 'aP1', 'ab')
d2.add_admixture_event('c', 'cP1', 'cP2')
d2.add_merge_event('ab', 'cP1', 'left')
d2.add_merge_event('d', 'cP2', 'right')
d2.add_merge_event('left', 'right', 'root')
for n in d2.nodes:
    d2.set_node_ne(n, 5000)            # diploid 5000 = haploid 10000
d2.set_admixture_parameters('a', time=10, fraction_parent_1=0.5, parent_1_name='aP1')
d2.set_merge_time('bP', 55)
d2.set_merge_time('ab', 75)
d2.set_admixture_parameters('c', time=500, fraction_parent_1=0.70, parent_1_name='cP1')
d2.set_merge_time('left', 700)
d2.set_merge_time('right', 900)
d2.set_merge_time('root', 1000)
d2.finalize_root()
save(d2, "topology_sim2.png", "Recent admixture-through-b on a + ancient admixture on c")


# ---- §6 v2 (varying Ne, success): weak (0.1) recent on b + old on d ----
v2 = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
v2.add_admixture_event('b', 'bP1', 'bP2')
v2.add_merge_event('a', 'bP1', 'ab')
v2.add_merge_event('c', 'bP2', 'cb')
v2.add_admixture_event('d', 'dP1', 'dP2')
v2.add_merge_event('cb', 'dP1', 'cbd')
v2.add_merge_event('ab', 'cbd', 'left')
v2.add_merge_event('e', 'dP2', 'right')
v2.add_merge_event('left', 'right', 'root')
for n in v2.nodes:
    v2.set_node_ne(n, 15000 // 2)
v2.set_node_ne('b', 5000 // 2)
v2.set_node_ne('d', 5000 // 2)
v2.set_admixture_parameters('b', time=20, fraction_parent_1=0.9, parent_1_name='bP1')
v2.set_admixture_parameters('d', time=500, fraction_parent_1=0.7, parent_1_name='dP1')
v2.set_merge_time('ab', 50)
v2.set_merge_time('cb', 100)
v2.set_merge_time('cbd', 700)
v2.set_merge_time('left', 1000)
v2.set_merge_time('right', 1200)
v2.set_merge_time('root', 1500)
v2.finalize_root()
save(v2, "topology_mixed_v2.png", "Varying Ne (success): recent on b + old on d")


# ---- §6 v4 (varying Ne, failure): two recent pulses on b + old on d ----
GAP = 5
v4 = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
v4.add_admixture_event('b', 'bP1_r', 'bM')
v4.add_merge_event('a', 'bP1_r', 'ab')
v4.add_admixture_event('bM', 'bP2_r', 'bMc')
v4.add_merge_event('ab', 'bP2_r', 'abP')
v4.add_merge_event('c', 'bMc', 'bC')
v4.add_admixture_event('d', 'dP1', 'dP2')
v4.add_merge_event('bC', 'dP1', 'cbd')
v4.add_merge_event('abP', 'cbd', 'left')
v4.add_merge_event('e', 'dP2', 'right')
v4.add_merge_event('left', 'right', 'root')
for n in v4.nodes:
    v4.set_node_ne(n, 10000 // 2)
v4.set_node_ne('d', 5000 // 2)
v4.set_admixture_parameters('b', time=10, fraction_parent_1=0.25, parent_1_name='bP1_r')
v4.set_merge_time('ab', 10 + GAP)
v4.set_admixture_parameters('bM', time=80, fraction_parent_1=(0.5 - 0.25) / (1 - 0.25), parent_1_name='bP2_r')
v4.set_merge_time('abP', 80 + GAP)
v4.set_merge_time('bC', 80 + 2 * GAP)
v4.set_admixture_parameters('d', time=700, fraction_parent_1=0.30, parent_1_name='dP1')
v4.set_merge_time('cbd', 900)
v4.set_merge_time('left', 1100)
v4.set_merge_time('right', 1200)
v4.set_merge_time('root', 1500)
v4.finalize_root()
save(v4, "topology_mixed_v4.png", "Varying Ne (failure): two recent pulses on b + old on d")

print("done")
