"""Rank the 21 topologies by ELBO and turn that into a selection probability.

Two different "probabilities" are reported, because they answer two different
questions and people routinely conflate them:

  P_model   softmax over models of the ELBO, i.e. treat the ELBO as a stand-in
            for the log marginal likelihood and put a uniform prior on the 21
            models.  This is the posterior model probability that was asked for.
            It is an ANSWER TO A QUESTION ABOUT THE WORLD, and it inherits every
            assumption below.

  P_argmax  TWO-LEVEL bootstrap: resample the Pathfinder restarts with
            replacement, resample the draws within each, take the max ELBO over
            the resampled restarts (the same statistic the ranking uses), and
            count how often each model comes out on top.  This answers only "is
            the ranking stable against the noise in the ELBO estimates?"  It
            says nothing about whether the ELBO is the right criterion.

            Both levels are needed and the seed level dominates.  The draw-level
            MC error on one ELBO is ~0.2 nats; the spread ACROSS restarts of the
            same model reached 11 nats in the first run -- larger than the gap
            between the top two models.  A single-level bootstrap over draws
            reported P = 1.000 for a lead that restarting the optimiser erases.
            Reporting the seed range next to the ELBO is not decoration.

Both will usually read 1.000 / 0.000.  That is not confidence, it is arithmetic:
the composite likelihood has ~119 terms and the models differ by hundreds of
nats, so any softmax saturates.  The informative column is the RAW GAP, which is
why delta-ELBO is printed next to every probability.

Three standing caveats, repeated in the written report so they travel with it:

  1. The ELBO is a LOWER BOUND on log Z, and the gap KL(q||p) differs per model.
     A model whose posterior is more Gaussian scores better for that reason
     alone.  logZ (importance sampling) would correct this, but the ESS column
     shows the weights are degenerate here, so it cannot be relied on either.
  2. An admixture graph gives the admixed leaf two source branches with separate
     Ne, i.e. a piecewise Ne.  These data are known to want that.  A win by an
     admixture graph is not by itself evidence of admixture.
  3. Model probabilities are conditional on this list of 21 being exhaustive.
     Two admixture events, or a non-constant Ne without admixture, are not in it.
"""
import os, json, argparse
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tag", required=True, help="Run directory under 3pop/.")
ap.add_argument("--boot", type=int, default=4000)
a = ap.parse_args()

CMP = os.path.join(_HERE, a.tag, "comparison")

blob = json.load(open(os.path.join(CMP, "all_fits.json")))
POPS = blob["pops"]
R = sorted(blob["results"], key=lambda r: -r["elbo"])
LW = np.load(os.path.join(CMP, "logw.npz"))

elbo = np.array([r["elbo"] for r in R])
d = elbo - elbo.max()
p_model = np.exp(d - np.logaddexp.reduce(d))

# --- two-level bootstrap: restarts, then draws within a restart -------------
rng = np.random.default_rng(0)
B = a.boot
stat = np.empty((B, len(R)))
for k, r in enumerate(R):
    seeds = list(r["elbo_by_seed"].keys())
    W = [LW[f"{r['index']}_s{s}"] for s in seeds]
    C = r["elbo"] - r["elbo_raw"]                      # dropped-constant correction
    ns = len(seeds)
    # per-restart bootstrapped ELBOs, then max over a bootstrap resample of restarts
    per = np.empty((B, ns))
    for q, w in enumerate(W):
        per[:, q] = w[rng.integers(0, w.size, size=(B, w.size))].mean(axis=1) + C
    pick = rng.integers(0, ns, size=(B, ns))
    stat[:, k] = np.take_along_axis(per, pick, axis=1).max(axis=1)
p_argmax = np.bincount(stat.argmax(axis=1), minlength=len(R)) / B

for r, pm, pa in zip(R, p_model, p_argmax):
    r["p_model"], r["p_argmax"] = float(pm), float(pa)
    r["d_elbo"] = float(r["elbo"] - elbo.max())

# ------------------------------------------------------------------ console
def row(r):
    f = (f"{r['admixture_fractions'][0]:.3f}" if r.get("admixture_fractions") else "-")
    ci = r["chi2_ibd"] / max(r.get("n_ibd_obs", 1), 1)
    cs = r["chi2_snp"] / max(r.get("n_snp_obs", 1), 1)
    spread = r["elbo_seed_max"] - r["elbo_seed_min"]
    return (f"{r['index']:>3} {r['newick']:<32} {r['n_admix']:>3} "
            f"{r['elbo']:>10.1f} {r['d_elbo']:>9.1f} {spread:>7.1f} "
            f"{r['ess']:>6.1f} {ci:>7.2f} {cs:>7.2f} {f:>6} "
            f"{r['p_model']:>8.3f} {r['p_argmax']:>8.3f}")


print("=" * 132)
print(f"TOPOLOGY RANKING BY ELBO   leaves {POPS}   mixed_model_Nsmooth, Pathfinder")
print("=" * 132)
print(f"{'#':>3} {'topology':<32} {'adm':>3} {'ELBO':>10} {'dELBO':>9} {'seedrng':>7} "
      f"{'ESS':>6} {'X2/n_i':>7} {'X2/n_s':>7} {'f':>6} {'P_model':>8} {'P_amax':>8}")
for r in R:
    print(row(r))

best_tree = next((r for r in R if r["n_admix"] == 0), None)
best_adm = next((r for r in R if r["n_admix"] == 1), None)
print("\n" + "-" * 132)
print(f"best overall     : [{R[0]['index']:02d}] {R[0]['newick']}   ELBO {R[0]['elbo']:+.1f}")
if best_tree:
    print(f"best tree        : [{best_tree['index']:02d}] {best_tree['newick']}   "
          f"ELBO {best_tree['elbo']:+.1f}")
if best_adm:
    print(f"best 1-admixture : [{best_adm['index']:02d}] {best_adm['newick']}   "
          f"ELBO {best_adm['elbo']:+.1f}")
if best_tree and best_adm:
    print(f"admixture buys   : {best_adm['elbo'] - best_tree['elbo']:+.1f} nats "
          f"for {best_adm['n_nodes'] - best_tree['n_nodes']} extra Ne + "
          f"{best_adm['n_events'] - best_tree['n_events']} extra times + 1 fraction")

gap = R[0]["elbo"] - R[1]["elbo"]
spread0 = R[0]["elbo_seed_max"] - R[0]["elbo_seed_min"]
# "gap vs spread" is too crude: a single bad restart inflates the spread without
# saying anything about the ranking.  The statistic that answers the question is
# how many of the winner's INDIVIDUAL restarts beat the runner-up's BEST one --
# i.e. would a one-restart run of this pipeline still have picked the winner.
w_seeds = np.array(list(R[0]["elbo_by_seed"].values()))
n_beat = int((w_seeds > R[1]["elbo"]).sum())
print(f"\ntop-two gap      : {gap:.1f} nats; winner's restarts span {spread0:.1f} nats "
      f"over {R[0]['n_seeds']} seeds")
print(f"restart robustness: {n_beat}/{len(w_seeds)} of the winner's individual restarts "
      f"beat the runner-up's BEST restart"
      + ("   -> lead survives a single-restart run"
         if n_beat == len(w_seeds) else
         "   -> a single-restart run could have ranked them the other way"))
print(f"max ESS over all models: {max(r['ess'] for r in R):.1f} of "
      f"{max(r['n_draws'] for r in R)} draws"
      + ("   -> logZ is NOT usable; read the ELBO column"
         if max(r["ess"] for r in R) < 50 else ""))

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(1, 3, figsize=(19, 6.4),
                       gridspec_kw={"width_ratios": [2.0, 1.0, 1.0]})
y = np.arange(len(R))[::-1]
col = ["#1f77b4" if r["n_admix"] == 0 else "#d62728" for r in R]
top = R[0]["elbo"]
lo = [r["elbo"] - r["elbo_seed_min"] for r in R]     # bar reaches the MAX seed
ax[0].barh(y, [r["d_elbo"] for r in R], color=col,
           xerr=[lo, [0] * len(R)], error_kw=dict(ecolor="k", lw=1, capsize=2))
ax[0].set_yticks(y)
ax[0].set_yticklabels([f"[{r['index']:02d}] {r['newick']}" for r in R], fontsize=8)
ax[0].set_xlabel("ELBO - best  (nats);  whisker = worst restart")
ax[0].set_title("Topology ranking\nblue = tree, red = one admixture", fontsize=11)
ax[0].axvline(0, color="k", lw=.8)

ci = [r["chi2_ibd"] / max(r.get("n_ibd_obs", 1), 1) for r in R]
cs = [r["chi2_snp"] / max(r.get("n_snp_obs", 1), 1) for r in R]
ax[1].barh(y - .2, ci, .4, color="#1f77b4", label="IBD")
ax[1].barh(y + .2, cs, .4, color="#d62728", label="SNP")
ax[1].axvline(1, color="k", ls="--", lw=1)
ax[1].set_xscale("log"); ax[1].set_yticks([]); ax[1].legend(fontsize=8)
ax[1].set_xlabel("chi2 / n"); ax[1].set_title("fit quality per component\n"
                                              "1.0 = inside its own SEs", fontsize=11)

fr = [(r["admixture_fractions"][0] if r.get("admixture_fractions") else np.nan)
      for r in R]
ax[2].barh(y, fr, color="#7f7f7f")
for v in (0.05, 0.95):
    ax[2].axvline(v, color="firebrick", ls=":", lw=1)
ax[2].set_xlim(0, 1); ax[2].set_yticks([])
ax[2].set_xlabel("admixture fraction f")
ax[2].set_title("f (dotted = collapsed to a tree)", fontsize=11)
fig.suptitle(f"Three-leaf topology selection, {'/'.join(POPS)}", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(CMP, "elbo_ranking.png"), dpi=130)
print(f"\n[plot] {CMP}/elbo_ranking.png")

# ------------------------------------------------------------------ files
import csv
with open(os.path.join(CMP, "elbo_table.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["rank", "index", "newick", "n_admix", "admixed", "elbo", "elbo_se",
                "elbo_seed_min", "elbo_seed_max", "n_seeds",
                "d_elbo", "logz", "ess", "chi2_ibd_per_n", "chi2_snp_per_n",
                "admix_fraction", "tau", "p_model", "p_argmax", "elbo_raw",
                "lp_const", "seed", "secs"])
    for k, r in enumerate(R, 1):
        w.writerow([k, r["index"], r["newick"], r["n_admix"], r["admixed"],
                    f"{r['elbo']:.3f}", f"{r['elbo_se']:.3f}",
                    f"{r['elbo_seed_min']:.3f}", f"{r['elbo_seed_max']:.3f}",
                    r["n_seeds"], f"{r['d_elbo']:.3f}",
                    f"{r['logz']:.3f}", f"{r['ess']:.2f}",
                    f"{r['chi2_ibd']/max(r.get('n_ibd_obs',1),1):.3f}",
                    f"{r['chi2_snp']/max(r.get('n_snp_obs',1),1):.3f}",
                    (f"{r['admixture_fractions'][0]:.4f}"
                     if r.get("admixture_fractions") else ""),
                    f"{r['tau']:.3f}", f"{r['p_model']:.4g}", f"{r['p_argmax']:.4g}",
                    f"{r['elbo_raw']:.3f}", f"{r['lp_const']:.3f}",
                    r["seed"], f"{r['secs']:.0f}"])

L = [f"# Topology selection on {' / '.join(POPS)}", "",
     f"21 models (3 trees + 18 one-admixture graphs), `mixed_model_Nsmooth`, "
     f"Pathfinder with {len(set(r['seed'] for r in R))}+ seeds, best ELBO kept.", "",
     "## Ranking", "",
     "| rank | topology | adm | ELBO | dELBO | restart range | ESS | chi2/n IBD | chi2/n SNP | f | P_model | P_argmax |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for k, r in enumerate(R, 1):
    f = (f"{r['admixture_fractions'][0]:.3f}" if r.get("admixture_fractions") else "-")
    L.append(f"| {k} | [`{r['newick']}`](../{r['name']}/report.md) | {r['n_admix']} | "
             f"{r['elbo']:+.1f} | {r['d_elbo']:+.1f} | "
             f"{r['elbo_seed_max']-r['elbo_seed_min']:.1f} | {r['ess']:.1f} | "
             f"{r['chi2_ibd']/max(r.get('n_ibd_obs',1),1):.2f} | "
             f"{r['chi2_snp']/max(r.get('n_snp_obs',1),1):.2f} | {f} | "
             f"{r['p_model']:.3f} | {r['p_argmax']:.3f} |")
L += ["", "![ranking](elbo_ranking.png)", "",
      "## Selection", "",
      f"- **Best model: `{R[0]['newick']}`**, ELBO {R[0]['elbo']:+.1f}, "
      f"P_model = {R[0]['p_model']:.3f}, P_argmax = {R[0]['p_argmax']:.3f}.",
      f"- Runner-up `{R[1]['newick']}` is {R[1]['d_elbo']:+.1f} nats behind. The winner's "
      f"{R[0]['n_seeds']} restarts span {spread0:.1f} nats, and **{n_beat} of them beat the "
      f"runner-up's best restart** -- "
      + ("the lead does not depend on the restart budget."
         if n_beat == len(w_seeds) else
         f"so {len(w_seeds)-n_beat} in {len(w_seeds)} single-restart runs would have "
         "ranked the two the other way. Treat the top two as close.")]
if best_tree and best_adm:
    L.append(f"- Best tree `{best_tree['newick']}` vs best admixture graph "
             f"`{best_adm['newick']}`: {best_adm['elbo'] - best_tree['elbo']:+.1f} nats "
             f"for {best_adm['n_nodes']-best_tree['n_nodes']} extra Ne, "
             f"{best_adm['n_events']-best_tree['n_events']} extra times and 1 fraction.")
L += ["", "## How to read the two probabilities", "",
      "`P_model` is softmax over the 21 ELBOs with a uniform model prior -- the posterior "
      "model probability. `P_argmax` bootstraps the Pathfinder draws and counts how often "
      "each model wins, so it measures only whether the RANKING is stable against Monte "
      "Carlo error. Both saturate at 1.000/0.000 whenever the gaps exceed a few nats, "
      "which they do here; the informative number is dELBO.", "",
      "## Caveats that travel with these numbers", "",
      "1. The ELBO is a lower bound on log Z and the gap KL(q||p) differs per model, so a "
      "model whose posterior happens to be more Gaussian is rewarded for that alone. The "
      "importance-sampling `logZ` would correct it, but the ESS column shows the weights "
      f"are degenerate (max ESS {max(r['ess'] for r in R):.0f} of "
      f"{max(r['n_draws'] for r in R)}), so `logZ` cannot be trusted here either.",
      "2. An admixture graph splits the admixed leaf into two source branches with separate "
      "Ne -- it buys a piecewise Ne, which these data independently want. A win by an "
      "admixture graph is evidence for non-constant Ne at least as much as for admixture, "
      "and three leaves cannot separate the two. Check each folder's residual row for a "
      "monotone tilt before reading any result as admixture.",
      "3. These probabilities are conditional on the 21 listed models being the whole "
      "hypothesis space. Two admixture events, and non-constant Ne without admixture, are "
      "not in it."]
open(os.path.join(CMP, "report.md"), "w").write("\n".join(L) + "\n")
print(f"[save] {CMP}/elbo_table.csv  and  {CMP}/report.md")
