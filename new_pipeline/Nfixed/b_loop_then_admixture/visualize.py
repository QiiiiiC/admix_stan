"""Archive search/parameter/residual tables; render parameter and residual figures only."""
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from study import SCENARIOS, SOURCES, VARIANTS, candidates, edges, save_json

PAIRS = [(i, j) for i in range(3) for j in range(i, 3)]
PAIR_NAMES = ["a–a", "a–b", "a–c", "b–b", "b–c", "c–c"]
COLORS = ["#236FA6", "#D87921", "#269B79", "#925BB7"]
LABELS = ["Omit\nshared", "Omit\nseparate", "Explicit\nshared", "Explicit\nseparate"]
COMPARISONS = [("omitted_backbone", "poisson_shared"), ("omitted_backbone", "poisson_separate"),
               ("explicit_loop", "poisson_shared"), ("explicit_loop", "poisson_separate")]


def csv_write(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def build_summary(c, out):
    out = Path(out); dest = out / "summary"; dest.mkdir(parents=True, exist_ok=True)
    graphs = defaultdict(list)
    for cand in candidates(): graphs[cand["graph"]].append(cand)
    ranking_rows = []; selections = []; diagnostics = []; parameters = []
    residual_handle = (dest / "residuals.tmp.csv").open("w", newline="")
    residual_writer = csv.DictWriter(residual_handle, ["scenario","pool","replicate","source","variant",
        "candidate","graph","reliable","kind","bin","pair","residual"])
    residual_writer.writeheader()
    plotted = defaultdict(list); all_fit_index = []
    for scenario in SCENARIOS:
        for pool in range(c["pools"]):
            for rep in range(c["replicates"]):
                for source in SOURCES:
                    for variant in VARIANTS:
                        ident = dict(scenario=scenario, pool=pool, replicate=rep, source=source, variant=variant)
                        folder = out / scenario / f"pool_{pool:03d}" / "fits" / f"rep_{rep:03d}" / source / variant
                        fits = {p.parent.name: json.loads(p.read_text()) for p in folder.glob("*/result.json")}
                        scores = {}; best_orders = {}; reliability = {}
                        for gid, orders in graphs.items():
                            rows = [fits.get(order["id"]) for order in orders]
                            if any(row is None or row["status"] != "ok" for row in rows): continue
                            scores[gid] = dict(logz=float(logsumexp([r["logz"] for r in rows])-np.log(len(rows))),
                                               elbo=float(logsumexp([r["elbo"] for r in rows])-np.log(len(rows))))
                            best_orders[gid] = max(rows, key=lambda row: row["logz"])
                            reliability[gid] = all(row["reliable"] for row in rows)
                        for cid, result in fits.items():
                            base = dict(ident, candidate=cid, graph=result["candidate"]["graph"], reliable=result["reliable"])
                            diagnostics.append(dict(base, status=result["status"], ess=result.get("ess"),
                                                    pareto_k=result.get("pareto_k"), seconds=result["seconds"]))
                            all_fit_index.append(dict(base, result=str((folder / cid / "result.json").relative_to(out)),
                                                      posterior_draws=str((folder / cid / "posterior_draws.npz").relative_to(out))
                                                      if (folder / cid / "posterior_draws.npz").exists() else None))
                            if result["status"] != "ok": continue
                            for name, stat in result.get("semantic_parameters", {}).items():
                                truth = stat["truth"]
                                parameters.append(dict(base, parameter=name, component="demography", **stat,
                                    error=stat["mean"]-truth if truth is not None else None,
                                    covered=stat["lower"] <= truth <= stat["upper"] if truth is not None else None,
                                    width=stat["upper"]-stat["lower"]))
                            for key, component in (("Ne", "shared"), ("Ne_ibd", "IBD"), ("Ne_snp", "SNP")):
                                if key not in result: continue
                                for j, node in enumerate(result["candidate"]["nodes"]):
                                    stat = {s: result[key][s][j] for s in ("mean", "lower", "upper")}
                                    parameters.append(dict(base, parameter=node, component=component, **stat,
                                        truth=c["haploid_ne"], error=stat["mean"]-c["haploid_ne"],
                                        covered=stat["lower"] <= c["haploid_ne"] <= stat["upper"], width=stat["upper"]-stat["lower"]))
                            for kind, values in result["residuals"].items():
                                values = np.asarray(values)
                                for i, j in PAIRS:
                                    if values.ndim == 3:
                                        for b, val in enumerate(values[:, i, j]):
                                            residual_writer.writerow(dict(base, kind=kind, bin=b, pair=f"{i},{j}", residual=float(val)))
                                    else:
                                        residual_writer.writerow(dict(base, kind=kind, bin=None, pair=f"{i},{j}", residual=float(values[i,j])))
                        # Separate the two scientific questions, reusing the same fitted candidates.
                        for search, ids in (("loop_omitted", list(range(1,22))), ("including_loop", list(range(1,23)))):
                            complete = all(gid in scores for gid in ids)
                            row = dict(ident, search=search, complete=complete,
                                       reliable=complete and all(reliability[g] for g in ids),
                                       winner=None, best_order=None, backbone_correct=None, full_graph_correct=None,
                                       logz_gap=None, elbo_winner=None)
                            if complete:
                                ordered = sorted(ids, key=lambda gid: scores[gid]["logz"], reverse=True)
                                win = ordered[0]
                                row.update(winner=win, best_order=best_orders[win]["candidate"]["id"],
                                           backbone_correct=win in (15,22),
                                           full_graph_correct=win == (22 if scenario == "b_loop" else 15),
                                           logz_gap=scores[win]["logz"]-scores[ordered[1]]["logz"],
                                           elbo_winner=max(ids,key=lambda g:scores[g]["elbo"]))
                                plotted[(scenario, source, "winner_"+search, variant)].append(best_orders[win])
                                for rank, gid in enumerate(ordered, 1):
                                    ranking_rows.append(dict(ident, search=search, rank=rank, graph=gid,
                                        newick=graphs[gid][0]["newick"], best_order=best_orders[gid]["candidate"]["id"],
                                        reliable=reliability[gid], **scores[gid]))
                            selections.append(row)
                        for gid, label in ((15,"omitted_backbone"), (22,"explicit_loop")):
                            if gid in best_orders:
                                plotted[(scenario,source,label,variant)].append(best_orders[gid])
    common = ["scenario", "pool", "replicate", "source", "variant"]
    csv_write(dest/"search_rankings.csv", ranking_rows, common+["search","rank","graph","newick","best_order","reliable","logz","elbo"])
    csv_write(dest/"search_selection.csv", selections, common+["search","complete","reliable","winner","best_order","backbone_correct","full_graph_correct","logz_gap","elbo_winner"])
    csv_write(dest/"fit_diagnostics.csv", diagnostics, common+["candidate","graph","status","reliable","ess","pareto_k","seconds"])
    csv_write(dest/"parameter_estimates.csv", parameters, common+["candidate","graph","reliable","parameter","component","mean","lower","upper","truth","error","covered","width"])
    residual_handle.close()
    (dest / "residuals.tmp.csv").replace(dest / "residuals.csv")
    save_json(dest/"fit_index.json", all_fit_index)
    # Raw/reliable-only aggregates are separate, never turn failed diagnostics into missing raw estimates.
    grouped = defaultdict(list)
    for row in parameters:
        if row["truth"] is None: continue
        for scope in ("all", "reliable"):
            if scope == "reliable" and not row["reliable"]: continue
            grouped[tuple(row[k] for k in ("scenario","source","variant","candidate","parameter","component"))+(scope,)].append(row)
    metrics = []
    for key, rows in grouped.items():
        errors = np.array([r["error"] for r in rows])
        metrics.append(dict(zip(["scenario","source","variant","candidate","parameter","component","scope"],key)) |
                       dict(n=len(rows), bias=float(errors.mean()), rmse=float(np.sqrt(np.mean(errors**2))),
                            coverage=float(np.mean([r["covered"] for r in rows])),
                            mean_interval_width=float(np.mean([r["width"] for r in rows]))))
    csv_write(dest/"parameter_metrics.csv", metrics,
              ["scenario","source","variant","candidate","parameter","component","scope","n","bias","rmse","coverage","mean_interval_width"])
    save_json(dest/"summary_status.json", dict(expected_fits=29*len(VARIANTS)*len(SOURCES)*len(SCENARIOS)*c["replicates"]*c["pools"],
        attempted=len(diagnostics), failed=sum(r["status"] != "ok" for r in diagnostics),
        reliable=sum(r["reliable"] for r in diagnostics),
        note="Overlapping genome subsamples. Rankings are composite-likelihood importance estimates. Parameter plots use weighted posterior means, including unreliable fits."))
    figures(c, plotted, dest)
    print(f"Archived {len(diagnostics)} fits; parameter and residual visualizations in {dest / 'figures'}")


def box(ax, groups, truth, labels=LABELS):
    for j, group in enumerate(groups):
        if not group: continue
        bp = ax.boxplot([group], positions=[j], widths=.55, patch_artist=True, showfliers=False,
                        medianprops=dict(color="#222222"))
        bp["boxes"][0].set(facecolor=COLORS[j], alpha=.35)
        ax.scatter(j+np.linspace(-.14,.14,len(group)),group,s=14,color=COLORS[j],zorder=3)
    if truth is not None: ax.axhline(truth,color="#A42E32",ls="--",lw=1.1)
    ax.set_xticks(range(len(groups)),[f"{labels[i]}\nn={len(g)}" for i,g in enumerate(groups)],fontsize=7)
    ax.set_xlim(-.6,len(groups)-.4); ax.grid(axis="y",alpha=.15)
    ax.spines[["top","right"]].set_visible(False)


def save(fig, folder, name, subtitle):
    note = ("Residual summaries include fits that fail numerical diagnostics." if name.startswith("residuals_") else
            "Dots = subsample posterior means; boxes = across-subsample IQR. Estimates include fits that fail numerical diagnostics.")
    fig.text(.5,.015,subtitle+"\n"+note,ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.07,1,.95))
    fig.savefig(folder/(name+".png"),dpi=160); fig.savefig(folder/(name+".pdf")); plt.close(fig)


def figures(c, data, dest):
    folder=dest/"figures"; folder.mkdir(exist_ok=True)
    from plot_topologies import topology_figures
    topology_figures(c, folder)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
    for scenario in SCENARIOS:
        for mode in ("matched_backbone", "selected_topology"):
            comps = COMPARISONS if mode == "matched_backbone" else [
                ("winner_loop_omitted","poisson_shared"),("winner_loop_omitted","poisson_separate"),
                ("winner_including_loop","poisson_shared"),("winner_including_loop","poisson_separate")]
            params=["b_split","left","right","root","b_fraction"]
            fig,axes=plt.subplots(2,5,figsize=(18,8))
            for i,source in enumerate(SOURCES):
                for j,param in enumerate(params):
                    groups=[[r["semantic_parameters"][param]["mean"] for r in data[(scenario,source,label,variant)]
                             if param in r.get("semantic_parameters",{})] for label,variant in comps]
                    truth=c["fractions"]["b"] if param=="b_fraction" else c["times"][param]
                    box(axes[i,j],groups,truth,
                        LABELS if mode=="matched_backbone" else ["Restricted\nshared","Restricted\nseparate","Augmented\nshared","Augmented\nseparate"])
                    axes[i,j].set_title(param.replace("_"," ")+f" (truth {truth:g})")
                    axes[i,j].set_ylabel(source+(" / fraction" if param=="b_fraction" else " / generations"))
            fig.suptitle(f"{scenario}: demographic parameters — {mode.replace('_',' ')}",fontsize=16)
            save(fig,folder,f"parameters_{mode}_{scenario}",
                 "Best-scoring event order within each graph. Selected-topology panels omit estimates without a matching backbone parameter; counts are explicit.")
        # Loop estimates have no corresponding truth in the no-loop control.
        fig,axes=plt.subplots(2,3,figsize=(12,8))
        for i,source in enumerate(SOURCES):
            for j,param in enumerate(("loop_open","loop_close","loop_major_fraction")):
                groups=[[r["semantic_parameters"][param]["mean"] for r in data[(scenario,source,"explicit_loop",v)]] for v in VARIANTS]
                truth=(c["times"][param] if param!="loop_major_fraction" else max(c["fractions"]["loop"],1-c["fractions"]["loop"])) if scenario=="b_loop" else None
                box(axes[i,j],groups,truth,["Shared","Separate"]); axes[i,j].set_title(param.replace("_"," ")); axes[i,j].set_ylabel(source)
        fig.suptitle(f"{scenario}: explicit-loop parameters",fontsize=16)
        save(fig,folder,f"loop_parameters_{scenario}","Loop fraction is max(f,1-f), because loop source labels are interchangeable. No loop-parameter truth line in the no-loop control.")
        for source in SOURCES:
            for label in ("omitted_backbone","explicit_loop"):
                example=next((r for v in VARIANTS for r in data[(scenario,source,label,v)]),None)
                if example is None: continue
                nodes=example["candidate"]["nodes"]
                fig,axes=plt.subplots((len(nodes)+3)//4,4,figsize=(17,3.5*((len(nodes)+3)//4)),squeeze=False)
                for ax,node in zip(axes.flat,nodes):
                    groups=[]
                    for variant,component in (("poisson_shared","Ne"),("poisson_separate","Ne_ibd"),("poisson_separate","Ne_snp")):
                        groups.append([r[component]["mean"][r["candidate"]["nodes"].index(node)] for r in data[(scenario,source,label,variant)]])
                    box(ax,groups,c["haploid_ne"],["Shared","Separate\nIBD","Separate\nSNP"]); ax.set_title(node); ax.set_ylabel("Haploid Ne")
                for ax in list(axes.flat)[len(nodes):]: ax.set_visible(False)
                fig.suptitle(f"{scenario} / {source} / {label}: branch sizes",fontsize=15)
                save(fig,folder,f"Ne_{scenario}_{source}_{label}","Explicit-loop nodes describe different time intervals from the collapsed b branch. Individual loop-source Ne labels can exchange.")
            # Count and SNP residuals side by side; curves and ranges across all available subsamples.
            x=(edges(c)[:-1]+edges(c)[1:])/2
            fig,axes=plt.subplots(2,4,figsize=(18,8))
            for pidx,(i,j) in enumerate(PAIRS):
                ax=list(axes.flat)[pidx]
                for k,(label,variant) in enumerate(COMPARISONS):
                    arr=np.asarray([r["residuals"]["ibd_pearson"] for r in data[(scenario,source,label,variant)]])
                    if len(arr):
                        y=arr[:,:,i,j]; q=np.quantile(y,[.25,.5,.75],axis=0)
                        ax.plot(x,q[1],color=COLORS[k],lw=1.3)
                        ax.fill_between(x,q[0],q[2],color=COLORS[k],alpha=.12)
                ax.axhline(0,color="black",lw=.8); ax.set_title(PAIR_NAMES[pidx]); ax.set_xlabel("IBD length (cM)"); ax.set_ylabel("Count Pearson residual")
            for idx,kind in ((6,"snp_raw"),(7,"snp_standardized")):
                ax=list(axes.flat)[idx]
                for k,(label,variant) in enumerate(COMPARISONS):
                    arr=np.asarray([r["residuals"][kind] for r in data[(scenario,source,label,variant)]])
                    if len(arr):
                        y=np.asarray([arr[:,i,j] for i,j in PAIRS]).T
                        q=np.quantile(y,[.25,.5,.75],axis=0); xx=np.arange(6)+(k-1.5)*.12
                        ax.errorbar(xx,q[1],yerr=[q[1]-q[0],q[2]-q[1]],color=COLORS[k],fmt="o",ms=3,capsize=2)
                ax.set_xticks(range(6),PAIR_NAMES,rotation=30); ax.axhline(0,color="black",lw=.8); ax.set_title(kind.replace("_"," "))
            fig.suptitle(f"{scenario} / {source}: residuals (observed − fitted)",fontsize=16)
            fig.legend([Line2D([0],[0],color=color) for color in COLORS],[x.replace("\n"," ") for x in LABELS],loc="upper center",bbox_to_anchor=(.5,.955),ncol=4,frameon=False)
            save(fig,folder,f"residuals_{scenario}_{source}","Curves/points = median residual across subsamples; bands/bars = IQR. Descriptive fitted-data checks, not independent predictions.")
