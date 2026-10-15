"""Topology data for the per-branch model, from a `DemographicTopology`.

The graph language records events in the order they were added; here that order is
used only as ONE listing that respects "must come after".  What the model receives
is the partial order itself: for every event the events that create its inputs, and
for every branch the event that creates it and the event that consumes it.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from methods import stan_data as sd
from methods.models import ModelSpec

STAN_FILE = HERE / "stan" / "mixed_model_Nsmooth_poisson_branch.stan"
# An absolute `file` makes ModelSpec.path point here, so methods' normaliser and
# Pathfinder wrapper apply unchanged.
SPEC = ModelSpec("mixed_Nsmooth_poisson_branch", str(STAN_FILE), "mixed", "smooth", "poisson")


def events_of(dem):
    """[(inputs, outputs)] per event, by branch name; outputs of an admixture are (source f, source 1 - f)."""
    return [((ev["children"], [ev["parent"]]) if ev["type"] == "MERGE" else ([ev["child"]], ev["parents"]))
            for ev in dem.ordered_events]


def topology_data(dem):
    idx = {name: i + 1 for i, name in enumerate(dem.nodes)}                 # 1-based branch index
    n, E = len(idx), len(dem.ordered_events)
    start, end = [0] * n, [0] * n
    event_admix, event_nodes, admix_event = [], [], []
    n_pred, pred, n_open, open_ = [], [], [], []
    live = set(dem.initial_leaves)
    for e, (ins, outs) in enumerate(events_of(dem), 1):
        if not set(ins) <= live:
            raise ValueError(f"event {e} consumes a branch that is not open: the listing breaks the partial order")
        before = sorted({start[idx[x] - 1] for x in ins} - {0})
        n_pred.append(len(before)); pred.append(before + [0] * (2 - len(before)))
        others = sorted(idx[x] for x in live - set(ins))
        n_open.append(len(others)); open_.append(others + [0] * (n - len(others)))
        if len(ins) == 2:
            event_admix.append(0); event_nodes.append([idx[ins[0]], idx[ins[1]], idx[outs[0]]])
        else:
            admix_event.append(e); event_admix.append(len(admix_event))
            event_nodes.append([idx[ins[0]], idx[outs[0]], idx[outs[1]]])
        for x in ins:
            end[idx[x] - 1] = e
        for x in outs:
            start[idx[x] - 1] = e
        live = (live - set(ins)) | set(outs)
    if len(live) != 1:
        raise ValueError("the graph does not end in a single root")
    smooth = sd.smooth_data(dem)
    return dict(n_leaves=len(dem.initial_leaves), n_nodes=n, n_events=E, n_admixture=len(admix_event),
                event_admix=event_admix, event_nodes=event_nodes, n_pred=n_pred, pred=pred, n_open=n_open,
                open=open_, admix_event=admix_event, start_event=start, end_event=end, root=idx[next(iter(live))],
                ne_parent=smooth["ne_parent"], ne_admix_idx=smooth["ne_admix_idx"])


def stan_data(dem, obs, edges, n_haploid, T_max):
    data = topology_data(dem)
    ibd = sd.ibd_data(edges, obs, n_haploid, T_max)
    data.update({k: ibd[k] for k in ("n_bins", "bin_length", "T_max", "cm", "n_samples", "ibd_count")})
    data.update(sd.snp_data(obs))
    return data


def heights_to_times(dem, height):
    """Waiting times (the model's `times`) that give the event heights `height` (one per event)."""
    data = topology_data(dem)
    return [height[e] - max([height[p - 1] for p in data["pred"][e][:data["n_pred"][e]]], default=0.0)
            for e in range(data["n_events"])]
