"""Exact log evidence from posterior draws: bridge sampling with a Gaussian-mixture proposal.

Used as an offline reference to calibrate the Pathfinder ELBO, not in the search itself.

    logZ = log p(data | graph) = log ∫ p(data | θ, graph) p(θ | graph) dθ

for a model normalised by `models.normalized_model` (every constant kept), i.e. the
marginal likelihood of the fitted summary statistics under the composite likelihood.

Why this estimator (tested 2026-09-28 on the admix_thru_b_plus_ancient models): the
posteriors are strongly skewed on Stan's unconstrained scale (weakly identified time
increments pile against their `times >= 1` bound: skew ~ -2, excess kurtosis 5-37), so
importance sampling from ANY single elliptical proposal fails (Pareto k 4-12 even for a
t fitted to the exact NUTS posterior), and plain / warp-3 bridge sampling reach only
0.16-0.3 nats.  Bridge sampling needs overlap rather than tail domination, and a mixture
proposal follows the skew: 6 components, 8,000 proposal draws -> 0.01-0.07 nats.

Procedure (Meng & Wong 1996; Overstall & Forster 2010 split): half the posterior draws
fit the mixture, the other half enter the bridge; the target is evaluated exactly with
CmdStan's log_prob (jacobian=True, the density Pathfinder's lp__ also uses).  The error
is the relative-MSE approximation of Fruhwirth-Schnatter (2004), ~ the sd of log Z.
"""
from __future__ import annotations

import csv
from pathlib import Path
import re

import numpy as np
from scipy.special import expit, logit, logsumexp

from .models import _block


def param_layout(spec, data):
    """[(name, size, lower, upper, scalar)] for every parameter, in declaration order."""
    block = re.sub(r"//[^\n]*", "", _block(spec.path.read_text(), "parameters"))
    ints = {k: v for k, v in data.items() if isinstance(v, (int, np.integer))}
    out = []
    for decl in (d.strip() for d in block.split(";")):
        if not decl:
            continue
        name = decl.split()[-1]
        lo = re.search(r"lower\s*=\s*([^,>]+)", decl)
        up = re.search(r"upper\s*=\s*([^,>]+)", decl)
        dims = re.findall(r"\[([^\]]+)\]", decl)
        size = 1
        for d in dims:
            for part in d.split(","):
                size *= int(eval(part.strip(), {"__builtins__": {}}, ints))
        out.append((name, size, float(lo.group(1)) if lo else None, float(up.group(1)) if up else None,
                    not dims))
    return [x for x in out if x[1] > 0]


def columns(layout):
    """CmdStan CSV column names of the parameters, in layout order."""
    return [name if scalar else f"{name}.{i + 1}" for name, size, _, _, scalar in layout for i in range(size)]


def unconstrain(layout, draws):
    """draws: mapping column -> values (e.g. a DataFrame of a Stan CSV) -> (n, dim)."""
    parts, j = [], 0
    names = columns(layout)
    for name, size, lo, up, _ in layout:
        x = np.column_stack([np.asarray(draws[c], float) for c in names[j:j + size]]); j += size
        parts.append(logit((x - lo) / (up - lo)) if lo is not None and up is not None
                     else np.log(x - lo) if lo is not None else x)
    return np.hstack(parts)


def constrain(layout, u):
    out, j = [], 0
    for name, size, lo, up, _ in layout:
        z = u[:, j:j + size]; j += size
        out.append(lo + (up - lo) * expit(z) if lo is not None and up is not None
                   else lo + np.exp(z) if lo is not None else z)
    return np.hstack(out)


def log_density(model, data, layout, u, path):
    """Exact log density on the unconstrained scale (with Jacobian) at draws u."""
    x = constrain(layout, u)
    path = Path(path)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["lp__"] + columns(layout))
        for row in x:
            writer.writerow([0.0] + [repr(float(v)) for v in row])
    return model.log_prob(params=str(path), data=data, jacobian=True, sig_figs=12)["lp__"].to_numpy()


def bridge(lp_post, lq_post, lp_prop, lq_prop, iters=1000, tol=1e-10):
    """Optimal bridge estimate of log Z and its relative-MSE error (~ sd of log Z)."""
    n1, n2 = len(lp_post), len(lp_prop)
    ls1, ls2 = np.log(n1 / (n1 + n2)), np.log(n2 / (n1 + n2))
    l1, l2 = lp_post - lq_post, lp_prop - lq_prop
    shift = np.median(l1)
    l1, l2 = l1 - shift, l2 - shift
    r = 0.0
    for _ in range(iters):
        num = logsumexp(l2 - np.logaddexp(ls1 + l2, ls2 + r)) - np.log(n2)
        den = logsumexp(-np.logaddexp(ls1 + l1, ls2 + r)) - np.log(n1)
        new = num - den
        done = abs(new - r) < tol
        r = new
        if done:
            break
    f2 = np.exp(l2 - np.logaddexp(ls1 + l2, ls2 + r))
    g1 = np.exp(-np.logaddexp(ls1 + l1, ls2 + r))
    error = np.sqrt(np.var(f2) / np.mean(f2)**2 / n2 + np.var(g1) / np.mean(g1)**2 / n1)
    return float(r + shift), float(error)


def log_evidence(model, data, spec, csv_files, folder, components=6, proposal_draws=8000, seed=0):
    """Bridge-sampling log Z from the NUTS output files of `model` fitted to `data`."""
    import pandas as pd
    from sklearn.mixture import GaussianMixture
    layout = param_layout(spec, data)
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    frames = [pd.read_csv(f, comment="#") for f in csv_files]
    u = unconstrain(layout, pd.concat(frames))
    lp = np.concatenate([model.log_prob(params=str(f), data=data, jacobian=True, sig_figs=12)["lp__"].to_numpy()
                         for f in csv_files])
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(u))
    fit_part, bridge_part = order[:len(u) // 2], order[len(u) // 2:]
    mixture = GaussianMixture(components, covariance_type="full", reg_covar=1e-6,
                              random_state=seed).fit(u[fit_part])
    x, _ = mixture.sample(proposal_draws)
    lp_x = log_density(model, data, layout, x, folder / "bridge_proposal.csv")
    logz, error = bridge(lp[bridge_part], mixture.score_samples(u[bridge_part]), lp_x, mixture.score_samples(x))
    return dict(logz=logz, error=error, posterior_draws=len(bridge_part), proposal_draws=proposal_draws,
                components=components, method="bridge_sampling_gmm")
