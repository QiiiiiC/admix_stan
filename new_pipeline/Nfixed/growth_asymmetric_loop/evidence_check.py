#!/usr/bin/env python3
"""Quick check of two evidence estimators on the winning backbone (g15_o1).

The full run's importance weights have Pareto k > 0.9 everywhere, so the
Pathfinder logZ is not a measurement.  Two cheap alternatives on the same fits:

  laplace   L-BFGS mode + Hessian (CmdStan `laplace`): logZ = lp(mode) +
            d/2 log 2pi + 1/2 log|Sigma|.  Also used as a Gaussian IS proposal.
  student-t heavy-tailed IS: multivariate t (nu=4, scale 1.5) centred either on
            the pooled Pathfinder draws or on the Laplace mode; log p evaluated
            with CmdStan `log_prob` in one batched call.

Runs on growth / true IBD / 10 subsamples / shared and separate.  Writes
runs/<out>/evidence_check/<variant>.json; `--report` collates to summary/.
"""
import argparse, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.special import logsumexp, gammaln
sys.path.insert(0, str(Path(__file__).parent))
from study import load_config, candidates, selections, aggregate, stan_data, save_json
from run_study import load_blocks, normalized_model, pool_folder, initial

HERE = Path(__file__).parent
NU, SCALE, DRAWS = 4.0, 1.5, 4000


def blocks_of(variant):
    return ("", ) if not variant.endswith("separate") else ("_ibd", "_snp")


def to_unconstrained(sv, variant, single=False):
    """Stan's transforms in declaration order: log(t-1), logit(f), mu, log sigma, log tau, Ne_raw.
    `single`: sv holds one point (scalars / 1-D); otherwise a batch with the draw axis first."""
    def col(a):
        a = np.asarray(a, dtype=float)
        return a.reshape(-1) if single else a.reshape(a.shape[0], -1)
    t = col(sv["times"]); f = col(sv["admixture_fractions"])
    parts = [np.log(t - 1), np.log(f / (1 - f))]
    for sfx in blocks_of(variant):
        parts += [col(sv["mu_log"+sfx]), np.log(col(sv["sigma_log"+sfx])), np.log(col(sv["tau"+sfx])), col(sv["Ne_raw"+sfx])]
    return np.concatenate(parts, axis=-1)


def _log_prob_call(exe, data_json, U, tmp):
    S, d = U.shape
    M = U.ravel().reshape(d, S).T
    ujson = tmp / "u.json"; out = tmp / "lp.csv"
    json.dump({"params_r": M.tolist()}, open(ujson, "w"))
    r = subprocess.run([str(exe), "log_prob", f"unconstrained_params={ujson}", "jacobian=1",
                        "data", f"file={data_json}", "output", f"file={out}"], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stdout[-300:] + r.stderr[-300:])
    lines = [ln.strip() for ln in open(out) if ln.strip() and not ln.startswith("#")]
    head = next(i for i, ln in enumerate(lines) if ln.startswith("lp__"))
    lp = [float(ln.split(",")[0]) for ln in lines[head+1:]]
    if len(lp) != S:
        raise RuntimeError(f"log_prob returned {len(lp)} rows for {S} points")
    return np.asarray(lp, dtype=float)


def batched_log_prob(exe, data_json, U, tmp, chunk=250):
    """CmdStan log_prob at many unconstrained points.  The JSON reader flattens
    2-D arrays column-major while log_prob slices the flat vector row-wise, so
    the matrix is permuted to compensate (verified against single-point calls).
    One point that throws inside the model (a far-tail draw overflowing an
    exp) aborts the whole call, so chunks that fail are redone point by point
    and such points get -inf: numerically zero posterior mass, counted as dropped."""
    lp = np.full(len(U), -np.inf)
    for a in range(0, len(U), chunk):
        try:
            lp[a:a+chunk] = _log_prob_call(exe, data_json, U[a:a+chunk], tmp)
        except RuntimeError:
            for i in range(a, min(a+chunk, len(U))):
                try:
                    lp[i] = _log_prob_call(exe, data_json, U[i:i+1], tmp)[0]
                except RuntimeError:
                    pass
    return lp


def mvn_logpdf(U, mu, cov):
    d = len(mu); L = np.linalg.cholesky(cov)
    z = np.linalg.solve(L, (U - mu).T)
    return -0.5*d*np.log(2*np.pi) - np.log(np.diag(L)).sum() - 0.5*(z**2).sum(0)


def mvt(rng, mu, cov, nu, n):
    d = len(mu); L = np.linalg.cholesky(cov)
    g = rng.chisquare(nu, n) / nu
    return mu + (rng.standard_normal((n, d)) @ L.T) / np.sqrt(g)[:, None]


def mvt_logpdf(U, mu, cov, nu):
    d = len(mu); L = np.linalg.cholesky(cov)
    z = np.linalg.solve(L, (U - mu).T); m = (z**2).sum(0)
    return (gammaln((nu+d)/2) - gammaln(nu/2) - 0.5*d*np.log(nu*np.pi) - np.log(np.diag(L)).sum()
            - 0.5*(nu+d)*np.log1p(m/nu))


def is_stats(lw):
    import arviz as az
    ok = np.isfinite(lw); lw = lw[ok]
    if ok.sum() < 2:
        return dict(logz=float("nan"), ess=float("nan"), pareto_k=float("nan"), dropped=int((~ok).sum()), max_weight=float("nan"))
    w = np.exp(lw - logsumexp(lw))
    try:
        _, k = az.psislw(lw)
    except Exception:          # degenerate weights (a handful of finite draws)
        k = float("nan")
    return dict(logz=float(logsumexp(lw) - np.log(len(lw))), ess=float(1/np.sum(w**2)),
                pareto_k=float(k), dropped=int((~ok).sum()), max_weight=float(w.max()))


def one_fit(model, data, variant, folder, tmp, rng):
    from cmdstanpy import write_stan_json
    res = json.loads((folder/"result.json").read_text())
    out = dict(pathfinder=dict(elbo=res["elbo"], logz=res["logz"], ess=res["ess"], pareto_k=res["pareto_k"]))
    z = np.load(folder/"posterior_draws.npz")
    U_pf = to_unconstrained(z, variant)
    lw_pf = z["log_weight"]
    data_json = tmp/"data.json"; write_stan_json(str(data_json), data)
    d = U_pf.shape[1]; out["d"] = int(d)

    # ---- Laplace: best of four L-BFGS starts (three dispersed + Pathfinder's best draw)
    best = None
    inits = [initial(data, variant, s) for s in (1, 7, 13)]
    top = int(np.argmax(lw_pf))
    inits.append({k: np.asarray(z[k])[top] for k in z.files if k not in ("log_weight", "weight")})
    for i, init in enumerate(inits):
        try:
            opt = model.optimize(data=data, inits=init, seed=100+i, algorithm="lbfgs", jacobian=True,
                                 iter=20000, output_dir=str(tmp/"opt"), show_console=False)
            lp = float(opt.optimized_params_dict["lp__"])
            if best is None or lp > best[0]:
                best = (lp, opt)
        except Exception as ex:
            out.setdefault("opt_errors", []).append(str(ex)[:200])
    lp_mode, opt = best
    lap = model.laplace_sample(data=data, mode=opt, draws=DRAWS, jacobian=True, seed=7,
                               output_dir=str(tmp/"lap"), show_console=False)
    sv = lap.stan_variables()
    U_lap = to_unconstrained(sv, variant)
    mu_lap = to_unconstrained(opt.stan_variables(), variant, single=True)
    cov_lap = np.cov(U_lap.T)
    sign, logdet = np.linalg.slogdet(cov_lap)
    out["laplace"] = dict(lp_mode=lp_mode, logz=float(lp_mode + 0.5*d*np.log(2*np.pi) + 0.5*logdet))
    log_p = np.asarray(lap.draws_pd()["log_p__"], dtype=float)
    out["laplace_is"] = is_stats(log_p - mvn_logpdf(U_lap, mu_lap, cov_lap))
    out["laplace_is"]["max_log_p_minus_mode"] = float(np.nanmax(log_p) - lp_mode) if np.isfinite(log_p).any() else float("nan")

    # ---- Student-t IS around the pooled Pathfinder draws and around the Laplace mode
    for tag, mu, cov in (("t_pathfinder", U_pf.mean(0), np.cov(U_pf.T)), ("t_laplace", mu_lap, cov_lap)):
        cov_t = SCALE**2 * cov + 1e-10*np.eye(d)
        U = mvt(rng, mu, cov_t, NU, DRAWS)
        lp = batched_log_prob(model.exe_file, data_json, U, tmp)
        out[tag] = is_stats(lp - mvt_logpdf(U, mu, cov_t, NU))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="poisson_shared")
    ap.add_argument("--out", type=Path, default=HERE/"runs"/"default")
    ap.add_argument("--scenario", default="growth"); ap.add_argument("--source", default="true")
    ap.add_argument("--candidate", default="g15_o1")
    a = ap.parse_args()
    c = load_config()
    cand = next(x for x in candidates() if x["id"] == a.candidate)
    model = normalized_model(a.variant, a.out)
    blocks = load_blocks(pool_folder(a.out, a.scenario, 0), c, [a.source])
    tmp = a.out/"evidence_check"/"tmp"/a.variant; tmp.mkdir(parents=True, exist_ok=True)
    dest = a.out/"evidence_check"/f"{a.scenario}_{a.source}_{a.candidate}_{a.variant}.json"
    results = json.loads(dest.read_text()) if dest.exists() else {}
    rng = np.random.default_rng(2026)
    for rep, selected in enumerate(selections(c, 0)):
        if str(rep) in results:
            continue
        t0 = time.monotonic()
        obs = aggregate(blocks, selected, a.source, c)
        data = stan_data(cand, obs, c)
        folder = pool_folder(a.out, a.scenario, 0)/"fits"/f"rep_{rep:03d}"/a.source/a.variant/a.candidate
        r = one_fit(model, data, a.variant, folder, tmp, rng)
        r["seconds"] = time.monotonic()-t0
        r = json.loads(json.dumps(r, default=str), parse_constant=lambda _: None)   # NaN -> null
        results[str(rep)] = r
        save_json(dest, results)
        f = lambda v, fmt: ("n/a" if v is None else format(v, fmt))
        print(f"rep {rep} {a.variant}: PF logz {f(r['pathfinder']['logz'],'.1f')} k {f(r['pathfinder']['pareto_k'],'.2f')} | "
              f"Laplace {f(r['laplace']['logz'],'.1f')} | " + " | ".join(
              f"{tag} {f(r[tag]['logz'],'.1f')} k {f(r[tag]['pareto_k'],'.2f')} ess {f(r[tag]['ess'],'.0f')}"
              for tag in ("laplace_is", "t_pathfinder", "t_laplace")) + f"  ({r['seconds']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
