"""Posterior summaries of fit quality, and true-size references.

`fit_summary` turns weighted posterior draws into JSON-ready parameter
summaries, posterior-mean predictions and residuals.  `residual_tilt` and
`trajectory_reference` support the question of what a constant-Ne branch is
really estimating when the truth changed size.
"""
from __future__ import annotations

import numpy as np
from scipy.special import xlogy

from .fitting import PREDICTIONS, parameter_summary
from .ibd import pair_counts

PARAMETERS = ("cumulative_times", "admixture_fractions", "effective_N", "Ne", "Ne_ibd", "Ne_snp")
COMPONENTS = ("lp_ibd", "lp_snp", "chi2_ibd", "chi2_snp")


def residual_tilt(pearson, edges):
    """Least-squares slope of the Pearson residual against segment length, per
    pair, in residual units per cM.  A monotone tilt is the fingerprint of a
    branch whose Ne changed within the fitted interval: a constant Ne cannot
    make both the short (older) and long (recent) bins agree at once."""
    pearson = np.asarray(pearson)
    x = (np.asarray(edges[:-1]) + np.asarray(edges[1:])) / 2
    x = x - x.mean()
    out = {}
    for i, j in zip(*np.triu_indices(pearson.shape[-1])):
        y = pearson[:, i, j]
        out[f"{i},{j}"] = dict(slope=float(np.sum(x * (y - y.mean())) / np.sum(x * x)),
                               correlation=float(np.corrcoef(x, y)[0, 1]) if np.std(y) > 0 else 0.0)
    return out


def fit_summary(sv, weights, obs, n_haploid, edges):
    """Parameters, predictions, residuals and component likelihoods of one fit.

    IBD residuals (when the model predicts IBD): fraction residual, count
    residual, Pearson (k - lambda)/sqrt(lambda) and signed deviance residual of
    the Poisson counts, and their tilt across length.  SNP residuals: raw and
    standardised by w_se.
    """
    upper = np.triu_indices(len(n_haploid))
    result = {name: parameter_summary(sv[name], weights) for name in PARAMETERS if name in sv}
    pred = {name: np.tensordot(weights, sv[name], axes=1) for name in PREDICTIONS if name in sv}
    result["predictions"] = {name: value.tolist() for name, value in pred.items()}
    result["prediction_intervals"] = {name: parameter_summary(sv[name], weights)
                                      for name in PREDICTIONS if name in sv}
    result["likelihood_summaries"] = {name: parameter_summary(sv[name], weights)
                                      for name in COMPONENTS if name in sv}
    checks, residuals = {}, {}
    if "ibd_fraction" in pred:
        residual = obs["ibd_hat"] - pred["ibd_fraction"]
        checks["ibd_rmse_by_pair"] = np.sqrt(np.mean(residual[:, upper[0], upper[1]]**2, axis=0)).tolist()
        residuals["ibd_fraction"] = residual.tolist()
    if "ibd_number" in pred:
        expected = np.maximum(pred["ibd_number"] * obs["cm"] * pair_counts(n_haploid), 1e-12)
        counts = obs["ibd_count"]
        deviance = 2 * (xlogy(counts, counts / expected) - counts + expected)
        pearson = (counts - expected) / np.sqrt(expected)
        residuals.update(ibd_count=(counts - expected).tolist(), ibd_pearson=pearson.tolist(),
                         ibd_deviance=(np.sign(counts - expected) * np.sqrt(np.maximum(deviance, 0))).tolist())
        result["residual_tilt"] = residual_tilt(pearson, edges)
    if "W_centered" in pred:
        raw = obs["w_hat"] - pred["W_centered"]
        checks["snp_rmse"] = float(np.sqrt(np.mean(raw[upper]**2)))
        checks["snp_standardized_residual"] = (raw / obs["w_se"]).tolist()
        residuals.update(snp_raw=raw.tolist(), snp_standardized=(raw / obs["w_se"]).tolist())
    result["fit_checks"] = checks
    result["residuals"] = residuals
    result["observations"] = {k: np.asarray(v).tolist() for k, v in obs.items()}
    if "Ne_ibd" in sv and "Ne_snp" in sv:
        result["log_ne_ibd_over_snp"] = parameter_summary(np.log(sv["Ne_ibd"] / sv["Ne_snp"]), weights)
    return result


def trajectory_reference(size, t0, t1, points=2001):
    """What a constant-Ne branch could hope to recover from a size trajectory.

    size(t) -> haploid Ne over [t0, t1].  Returns the sizes at both ends, the
    harmonic mean (accumulated drift is duration / harmonic mean, so this is the
    SNP-side target) and the arithmetic mean.  An infinite t1 (root) is constant.
    """
    if not np.isfinite(t1):
        n = float(np.asarray(size(np.asarray(float(t0)))))
        return dict(t0=t0, t1=None, young=n, old=n, harmonic=n, arithmetic=n)
    grid = np.linspace(t0, t1, points)
    n = np.asarray(size(grid), float)
    return dict(t0=t0, t1=t1, young=float(n[0]), old=float(n[-1]),
                harmonic=float((t1 - t0) / np.trapz(1 / n, grid)),
                arithmetic=float(np.trapz(n, grid) / (t1 - t0)))


def true_parameters(dem):
    """(name, truth) for every event time then every admixture fraction of a fully
    parameterised topology, in the order of cumulative_times / admixture_fractions,
    for recovery and relative-error tables."""
    out = [(f"t_{ev['parent'] if ev['type'] == 'MERGE' else ev['child']}", t)
           for ev, t in zip(dem.ordered_events, dem.event_times())]
    out += [(f"f_{ev['child']}", f) for ev, f in
            zip([e for e in dem.ordered_events if e["type"] == "ADMIXTURE"], dem.admixture_fractions())]
    return out
