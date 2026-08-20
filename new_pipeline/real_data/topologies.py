"""
Candidate admixture-graph topologies over the five 1000-Genomes superpopulations.

A topology here is purely structural: it lists the MERGE / ADMIXTURE events.  The
event *parameters* (times, admixture fractions, effective sizes) are NOT set --
they are exactly what the Stan model infers in ``infer_topology.py``.

IMPORTANT
---------
A topology's leaves MUST be in the same order as the ``pop_order`` saved by
``generate_stan_data.py`` (default ``AFR, AMR, EAS, EUR, SAS``), because
``assemble_stan_data`` matches the IBD/SNP matrix rows to ``dem.initial_leaves``
by position.  All builders below take ``pop_order`` and use it verbatim as the
leaf order, so just keep the two consistent.

Add your own topology by writing a ``build_*(pop_order)`` function that returns a
``DemographicTopology`` and registering it in ``TOPOLOGIES`` at the bottom.  You
can also point ``infer_topology.py --topology /path/to/file.py`` at any module
that defines ``build_topology(pop_order)``.
"""

import os
import sys
import importlib.util

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_THIS_DIR)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from demography import DemographicTopology  # noqa: E402

# Must match generate_stan_data.DEFAULT_POP_ORDER.
DEFAULT_POP_ORDER = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def _check_leaves(pop_order):
    """Leaves are taken verbatim from pop_order; just sanity-check membership."""
    need = set(DEFAULT_POP_ORDER)
    have = set(pop_order)
    if have != need:
        raise ValueError(
            f"pop_order {list(pop_order)} must be a permutation of "
            f"{DEFAULT_POP_ORDER} for the built-in topologies."
        )
    return list(pop_order)


def build_amr_admixture(pop_order=DEFAULT_POP_ORDER):
    """Out-of-Africa tree with AMR as a 2-way admixed population.

    Structure (backwards in time):
      * AMR splits into a Eurasian source and an African source (1 ADMIXTURE).
      * The Eurasian source joins EUR; (EUR+src) joins SAS, then EAS  -> OOA.
      * The African source joins AFR.
      * OOA and (AFR + African source) coalesce at the root.

    This is the natural minimal model of New-World admixture in the 1000G
    superpopulations (a 2-way proxy for the real EUR/AFR/Native-American mix,
    since there is no Native-American leaf here).
    """
    leaves = _check_leaves(pop_order)
    d = DemographicTopology(leaves)

    d.add_admixture_event("AMR", "AMR_eur", "AMR_afr")   # AMR = Eurasian + African
    d.add_merge_event("EUR", "AMR_eur", "eur_src")        # Eurasian source -> EUR
    d.add_merge_event("eur_src", "SAS", "weur")           # West-Eurasian clade
    d.add_merge_event("weur", "EAS", "ooa")               # out-of-Africa clade
    d.add_merge_event("AFR", "AMR_afr", "afr_src")        # African source -> AFR
    d.add_merge_event("ooa", "afr_src", "root")           # root
    return d


def build_tree_no_admix(pop_order=DEFAULT_POP_ORDER):
    """Plain out-of-Africa tree, no admixture (a null comparison topology).

      ((( (EUR, SAS), EAS), AMR), AFR)
    AMR is treated as a normal Eurasian sister here (no New-World admixture).
    """
    leaves = _check_leaves(pop_order)
    d = DemographicTopology(leaves)

    d.add_merge_event("EUR", "SAS", "eur_sas")
    d.add_merge_event("eur_sas", "EAS", "weur_eas")
    d.add_merge_event("weur_eas", "AMR", "eurasia")
    d.add_merge_event("eurasia", "AFR", "root")
    return d


def build_tree_2pop(pop_order=("AFR", "EUR")):
    """Two leaves, one merge.  The minimal graph: a single split.

      MERGE  AFR + EUR -> root          1 event, 3 nodes, 0 admixtures

    Diagnostic topology.  There is essentially nothing here to misspecify -- no
    admixture, no branching order, no internal branches -- so if the SNP and IBD
    likelihoods still disagree on this graph, the disagreement is mechanical
    (estimator scale, normalisation, priors) rather than a wrong demography.

    Worth knowing before reading any SNP result on two leaves: the centered
    covariance is rank 1.  Both lineages sit in their own leaf until the root,
    so the cross term snp_mean[AFR,EUR] is exactly 0 and double-centering gives

        W_centered = g * [[1, -1], [-1, 1]],   g = (t/4)(1/Ne_AFR + 1/Ne_EUR)

    i.e. the SNP data is ONE number and cannot separate t from Ne, nor Ne_AFR
    from Ne_EUR.  IBD identifies all three from the length spectrum, so the
    clean test is whether the IBD-implied g matches the observed one.
    """
    leaves = list(pop_order)
    if len(leaves) != 2:
        raise ValueError(f"pop_order {leaves} must have exactly 2 populations.")
    d = DemographicTopology(leaves)
    d.add_merge_event(leaves[0], leaves[1], "root")
    return d


def build_tree_3pop(pop_order):
    """((A, B), C) -- two sister leaves plus an outgroup.  2 merges, 5 nodes.

    Three leaves is the SMALLEST graph on which SNP drift is identifiable per
    branch, and that is the whole reason to prefer it over build_tree_2pop.  A
    double-centered symmetric L x L covariance carries L(L-1)/2 free numbers; an
    unrooted tree on L leaves has 2L-3 branches:

        L = 2 :  1 statistic,  1 branch  -> only x_A + x_B, never the split
        L = 3 :  3 statistics, 3 branches -> x_A, x_B, x_C each identified
        L = 4 :  6 statistics, 5 branches -> over-determined, topology testable

    With A and B sister leaves their branch DURATIONS are equal, so the SNP-side
    ratio x_A / x_B estimates Ne_B / Ne_A with no reference to t at all -- a
    t-free quantity that can be compared directly against the IBD-side Ne ratio.
    On two leaves that comparison is impossible in principle, not just noisy.

    Admixture is tested by the model itself: `3pop/enumerate_3pop.py` enumerates
    every graph on these leaves with at most one admixture event and
    `3pop/fit_topologies.py` ranks them, so a leaf that is really a mixture shows
    up as an admixture graph winning rather than as a pre-fit exclusion.
    """
    leaves = list(pop_order)
    if len(leaves) != 3:
        raise ValueError(f"pop_order {leaves} must have exactly 3 populations.")
    d = DemographicTopology(leaves)
    d.add_merge_event(leaves[0], leaves[1], "sister")   # the two close leaves
    d.add_merge_event("sister", leaves[2], "root")      # outgroup joins
    return d


def build_tree_no_admix_4pop(pop_order=("AFR", "EAS", "EUR", "SAS")):
    """Plain out-of-Africa tree on 4 populations, NO admixture.

      ((EUR, SAS), EAS), AFR)

    Backwards in time:
      MERGE  EUR + SAS      -> eur_sas      # West Eurasian ancestor
      MERGE  eur_sas + EAS  -> eurasia      # out-of-Africa ancestor
      MERGE  eurasia + AFR  -> root

    3 merges, 0 admixtures, 7 nodes.  This is the null topology: it is what the
    two-admixture graph reduces to if neither admixture edge is needed.  Worth
    fitting first now that the leaves are single subpopulations rather than
    pooled superpopulations -- the two admixture edges in build_two_admix_4pop
    were added to absorb ASW/ACB European ancestry inside the pooled AFR
    superpopulation and the SAS West/East mixture, and the first of those
    problems is removed by construction once AFR is a single subpopulation.
    """
    leaves = list(pop_order)
    need = {"AFR", "EAS", "EUR", "SAS"}
    if set(leaves) != need:
        raise ValueError(f"pop_order {leaves} must be a permutation of {sorted(need)}.")
    d = DemographicTopology(leaves)
    d.add_merge_event("EUR", "SAS", "eur_sas")      # 1
    d.add_merge_event("eur_sas", "EAS", "eurasia")  # 2
    d.add_merge_event("eurasia", "AFR", "root")     # 3
    return d


def build_amr_2way_recent(pop_order=DEFAULT_POP_ORDER):
    """Corrected AMR 2-way admixture: recent EUR+AFR sources, deep AFR split.

    Fixes the event-ordering bug in ``build_amr_admixture``: there the African
    source rejoined AFR at an event placed AFTER out-of-Africa, so the model's
    monotonic cumulative_times forced AMR's African ancestry to be older than
    OOA.  Here the recent admixture-source rejoins (EUR+AMR_eur, AFR+AMR_afr)
    are the FIRST events, and the deep continental splits come later, so the
    chronology is respected:

      recent:  AMR = EUR-source + AFR-source ; each source rejoins its continent
      deep:    (EUR2, SAS) -> West Eurasia ; +EAS -> OOA ; +AFR2 -> root

    Still a 2-way model (no Native-American leaf among the 5 superpops), so the
    AMR admixture fraction is only an EUR/AFR proxy of the true 3-way mix.
    """
    leaves = _check_leaves(pop_order)
    d = DemographicTopology(leaves)

    d.add_admixture_event("AMR", "AMR_eur", "AMR_afr")   # recent New-World mix
    d.add_merge_event("EUR", "AMR_eur", "EUR2")           # EUR source -> EUR (recent)
    d.add_merge_event("AFR", "AMR_afr", "AFR2")           # AFR source -> AFR (recent)
    d.add_merge_event("EUR2", "SAS", "weur")              # West-Eurasian clade
    d.add_merge_event("weur", "EAS", "ooa")               # out-of-Africa clade
    d.add_merge_event("AFR2", "ooa", "root")              # deep AFR-vs-rest split
    return d


def build_amr_3way(pop_order=DEFAULT_POP_ORDER):
    """AMR as a 3-way mix (EUR + AFR + Native-American proxied by EAS).

    The Native-American source has no leaf in the 5 superpopulations, but it is
    East-Asian-related, so we model it as a lineage that diverges from EAS.  AMR
    is a nested admixture: first EUR-source vs (everything else), then within the
    rest African-source vs Native-American-source.  Each source rejoins its
    continent at a biologically appropriate depth (EUR/AFR recent, NAT~EAS at the
    ~20-26 kya NAT-EAS divergence), which also lets the EUR/AFR admixture
    fractions be identified instead of collapsing.
    """
    leaves = _check_leaves(pop_order)
    d = DemographicTopology(leaves)

    d.add_admixture_event("AMR", "AMR_eur", "AMR_ne")     # European vs non-European
    d.add_admixture_event("AMR_ne", "AMR_afr", "AMR_na")  # African vs Native-American
    d.add_merge_event("EUR", "AMR_eur", "EUR2")           # EUR source -> EUR (recent)
    d.add_merge_event("AFR", "AMR_afr", "AFR2")           # AFR source -> AFR (recent)
    d.add_merge_event("EAS", "AMR_na", "EAS2")            # NAT source -> EAS (~20-26 kya)
    d.add_merge_event("EUR2", "SAS", "weur")              # West-Eurasian clade
    d.add_merge_event("weur", "EAS2", "ooa")              # out-of-Africa clade
    d.add_merge_event("AFR2", "ooa", "root")              # deep AFR-vs-rest split
    return d


def build_sas_admix_4pop(pop_order=("AFR", "EAS", "EUR", "SAS")):
    """4-population topology (AMR dropped) with SAS as a 2-way admixture.

    Leaves: AFR, EAS, EUR, SAS  (must match the data pop_order).

    Backwards in time:
      * SAS splits into two sources, SAS_a and SAS_b (1 ADMIXTURE).
      * SAS_a joins EUR        -> eur_sas    (West-Eurasian / ANI-like source)
      * eur_sas joins EAS      -> eurasia
      * eurasia joins SAS_b    -> ooa        (deep AASI-like source)
      * ooa joins AFR          -> root       (out-of-Africa split)

    Motivation: AMR is the most internally structured superpop (it pools several
    New-World groups), which inflates apparent recent Ne; dropping it and treating
    SAS as genuinely admixed isolates that substructure question.
    """
    leaves = list(pop_order)
    need = {"AFR", "EAS", "EUR", "SAS"}
    if set(leaves) != need:
        raise ValueError(f"pop_order {leaves} must be a permutation of {sorted(need)}.")
    d = DemographicTopology(leaves)
    d.add_admixture_event("SAS", "SAS_a", "SAS_b")   # SAS = source_a + source_b
    d.add_merge_event("EUR", "SAS_a", "eur_sas")      # ANI-like source -> EUR
    d.add_merge_event("eur_sas", "EAS", "eurasia")    # + EAS
    d.add_merge_event("eurasia", "SAS_b", "ooa")      # + deep SAS source
    d.add_merge_event("ooa", "AFR", "root")           # out-of-Africa split
    return d


def build_sas_bridge_4pop(pop_order=("AFR", "EAS", "EUR", "SAS")):
    """4-pop topology with SAS as a EUR x EAS bridge (literature-grounded).

    Modern South Asians = ANI (West-Eurasian: Iranian-farmer + Steppe) + AASI
    (a deep East-Eurasian lineage, sister to East Asians, diverged ~45 kya;
    Narasimhan 2019).  So SAS is a recent 2-way admixture whose two sources sit
    on OPPOSITE sides of the West/East-Eurasian split:

      ADMIX  SAS -> SAS_ani, SAS_aasi          # recent (Bronze-Age) mixture
      MERGE  EUR + SAS_ani -> west_eurasian      # ANI joins West Eurasia
      MERGE  EAS + SAS_aasi -> east_eurasian     # AASI joins East Eurasia (~45 kya)
      MERGE  west_eurasian + east_eurasian -> ooa  # West-East split (~40 kya)
      MERGE  ooa + AFR -> root                    # out-of-Africa (~60-120 kya)

    Refines build_sas_admix_4pop, where the second SAS source was placed basal
    to OOA; here it correctly sits on the East-Eurasian (EAS) side.
    """
    leaves = list(pop_order)
    need = {"AFR", "EAS", "EUR", "SAS"}
    if set(leaves) != need:
        raise ValueError(f"pop_order {leaves} must be a permutation of {sorted(need)}.")
    d = DemographicTopology(leaves)
    d.add_admixture_event("SAS", "SAS_ani", "SAS_aasi")   # ANI (West) + AASI (East)
    d.add_merge_event("EUR", "SAS_ani", "west_eurasian")   # ANI -> West Eurasia
    d.add_merge_event("EAS", "SAS_aasi", "east_eurasian")  # AASI -> East Eurasia
    d.add_merge_event("west_eurasian", "east_eurasian", "ooa")  # West-East split
    d.add_merge_event("ooa", "AFR", "root")                # out-of-Africa
    return d


def build_two_admix_4pop(pop_order=("AFR", "EAS", "EUR", "SAS")):
    """4-pop graph with TWO admixture events (AFR and the EUR+SAS ancestor).

    Leaves: AFR, EAS, EUR, SAS  (must match the data pop_order).

    Backwards in time:
      MERGE  EUR + SAS   -> a
      ADMIX  AFR         -> b, c         # AFR is 2-way admixed (b + c)
      ADMIX  a           -> d, e         # the EUR+SAS ancestor is 2-way admixed (d + e)
      MERGE  c + d        -> cd
      MERGE  e + EAS      -> e_eas
      MERGE  b + cd       -> b_cd
      MERGE  b_cd + e_eas -> root

    2 admixture events + 5 merges = 7 events, 13 nodes.  admixture_fractions[0]
    is AFR's fraction from source b (1-frac from c); admixture_fractions[1] is
    the EUR+SAS ancestor's fraction from source d (1-frac from e).
    """
    leaves = list(pop_order)
    need = {"AFR", "EAS", "EUR", "SAS"}
    if set(leaves) != need:
        raise ValueError(f"pop_order {leaves} must be a permutation of {sorted(need)}.")
    d = DemographicTopology(leaves)
    d.add_merge_event("EUR", "SAS", "a")          # 1
    d.add_admixture_event("AFR", "b", "c")         # 2: AFR = b + c
    d.add_admixture_event("a", "d", "e")           # 3: a = d + e
    d.add_merge_event("c", "d", "cd")              # 4
    d.add_merge_event("e", "EAS", "e_eas")         # 5
    d.add_merge_event("b", "cd", "b_cd")           # 6
    d.add_merge_event("b_cd", "e_eas", "root")     # 7
    return d


# Registry of built-in topologies.  Keys are the --topology names.
TOPOLOGIES = {
    "two_admix_4pop": build_two_admix_4pop,
    "amr_admixture": build_amr_admixture,
    "amr_2way_recent": build_amr_2way_recent,
    "amr_3way": build_amr_3way,
    "tree_no_admix": build_tree_no_admix,
    "sas_admix_4pop": build_sas_admix_4pop,
    "sas_bridge_4pop": build_sas_bridge_4pop,
}


def load_topology(name_or_path, pop_order=DEFAULT_POP_ORDER):
    """Resolve --topology to a DemographicTopology.

    ``name_or_path`` is either a key in ``TOPOLOGIES`` or a path to a .py file
    that defines ``build_topology(pop_order)``.
    """
    if name_or_path in TOPOLOGIES:
        return TOPOLOGIES[name_or_path](pop_order)

    if name_or_path.endswith(".py") and os.path.exists(name_or_path):
        spec = importlib.util.spec_from_file_location("user_topology", name_or_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        if not hasattr(mod, "build_topology"):
            raise AttributeError(
                f"{name_or_path} must define build_topology(pop_order)."
            )
        return mod.build_topology(pop_order)

    raise ValueError(
        f"Unknown topology {name_or_path!r}. Choose one of "
        f"{sorted(TOPOLOGIES)} or pass a .py file defining build_topology()."
    )


if __name__ == "__main__":
    # Quick structural sanity check of every built-in topology.
    for nm, builder in TOPOLOGIES.items():
        d = builder()
        print(f"\n=== {nm} ===")
        d.is_valid()
        d.print_summary()
