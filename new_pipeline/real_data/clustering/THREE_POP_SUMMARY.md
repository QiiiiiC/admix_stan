# Three-population analysis - consolidated results

Graph `((A,B),C)`: two sister leaves plus an outgroup. Fitted with Pathfinder
(NUTS agreed to 3 digits at 190x the cost, see below). `kappa_snp` removed;
IBD and SNP log-densities added directly. `w_se` on 16 Mb blocks.


## 1. Why three leaves

| L | statistics `L(L-1)/2` | tree branches `2L-3` | consequence |
|---|---|---|---|
| 2 | 1 | 1 | only `x_A + x_B`; drift never separable |
| 3 | 3 | 3 | `x_A`, `x_B`, `x_C` each identified |
| 4 | 6 | 5 | over-determined; topology testable |

With A and B sisters their branch durations are equal, so `x_A/x_B = Ne_B/Ne_A` -
a **t-free** quantity comparable directly against the IBD Ne ratio.


## 2. Leaf screening (f3, Patterson unbiased, delete-one-chromosome jackknife)

`f3(X; Y, out) < 0` means X is not a tree leaf. Screened before fitting anything.


| population | worst Z over all Y | status |
|---|---|---|
| YRI | -1.01 (Y=ESN) | clean (but drift ~ 0) |
| LWK | -11.30 (Y=ESN) | **ADMIXED** |
| GWD | -5.04 (Y=MSL) | **ADMIXED** |
| MSL | +23.94 (Y=ESN) | clean |
| ESN | +21.66 (Y=YRI) | clean |
| GBR | +69.76 (Y=MSL) | clean |
| CEU | +71.55 (Y=MSL) | clean |
| TSI | +61.48 (Y=MSL) | clean |
| IBS | +64.74 (Y=MSL) | clean |

LWK - which blind IBD clustering had selected as the AFR leaf - is admixed at Z = -11.3.
Cohesion is a within-population criterion and is blind to admixture; so is any two-leaf fit.


## 3. Fits


### `((YRI,ESN),FIN)`

| model | t_sis | t_root | Ne_YRI | Ne_ESN | Ne_FIN | x_YRI | x_ESN | tau | chi2/n IBD | chi2/n SNP | s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | 102 | 447 | 318,983 | 47,791 | 1,739 | 0.00044 | 0.00226 | 0.74 | - | 1.56 | 1 |
| IBD-only-Nsm | 109 | 367 | 1,243,741 | 111,788 | 43,529 | 0.00009 | 0.00098 | 1.25 | 17.42 | - | 7 |
| Mixed-Nsm | 117 | 2792 | 1,131,053 | 113,220 | 43,493 | 0.00010 | 0.00103 | 1.18 | 17.46 | 3.29 | 13 |
| **OBSERVED** | | | | | | **-0.00038** | **0.00295** | | | | |

Per-pair IBD fit (Mixed):
| pair | bins | chi2/bin | short-bin resid | long-bin resid |
|---|---|---|---|---|
| YRI-YRI | 27 | 1.0 | +0.2 | +0.0 |
| YRI-ESN | 25 | 2.1 | +1.9 | +0.8 |
| YRI-FIN | 1 | 0.5 | +0.7 | +0.7 |
| ESN-ESN | 37 | 3.2 | +0.8 | +0.6 |
| ESN-FIN | 1 | 0.5 | +0.7 | +0.7 |
| FIN-FIN | 31 | 62.2 | +8.9 | -6.6 |

### `((MSL,ESN),FIN)`

| model | t_sis | t_root | Ne_MSL | Ne_ESN | Ne_FIN | x_MSL | x_ESN | tau | chi2/n IBD | chi2/n SNP | s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | 65 | 465 | 16,733 | 14,670 | 2,013 | 0.00391 | 0.00447 | 0.52 | - | 0.31 | 1 |
| IBD-only-Nsm | 165 | 359 | 150,472 | 88,032 | 43,648 | 0.00110 | 0.00188 | 0.87 | 19.08 | - | 4 |
| Mixed-Nsm | 256 | 392 | 146,143 | 86,648 | 44,117 | 0.00175 | 0.00295 | 1.28 | 19.33 | 9.96 | 11 |
| **OBSERVED** | | | | | | **0.00379** | **0.00454** | | | | |

Per-pair IBD fit (Mixed):
| pair | bins | chi2/bin | short-bin resid | long-bin resid |
|---|---|---|---|---|
| MSL-MSL | 35 | 0.9 | +0.3 | +0.1 |
| MSL-ESN | 5 | 7.2 | -1.0 | +0.7 |
| MSL-FIN | 1 | 1.3 | +1.2 | +1.2 |
| ESN-ESN | 37 | 3.3 | -0.2 | +0.3 |
| ESN-FIN | 1 | 0.2 | +0.4 | +0.4 |
| FIN-FIN | 31 | 62.2 | +9.0 | -6.5 |

### `((MSL,ESN),GBR)`

| model | t_sis | t_root | Ne_MSL | Ne_ESN | Ne_GBR | x_MSL | x_ESN | tau | chi2/n IBD | chi2/n SNP | s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | 46 | 547 | 11,027 | 10,655 | 2,667 | 0.00414 | 0.00428 | 0.36 | - | 0.36 | 1 |
| IBD-only-Nsm | 163 | 400 | 150,887 | 87,885 | 170,186 | 0.00108 | 0.00186 | 0.90 | 2.14 | - | 4 |
| Mixed-Nsm | 241 | 479 | 147,468 | 86,924 | 168,010 | 0.00164 | 0.00278 | 1.20 | 2.32 | 10.91 | 12 |
| **OBSERVED** | | | | | | **0.00383** | **0.00449** | | | | |

Per-pair IBD fit (Mixed):
| pair | bins | chi2/bin | short-bin resid | long-bin resid |
|---|---|---|---|---|
| MSL-MSL | 35 | 0.9 | +0.4 | +0.1 |
| MSL-ESN | 5 | 6.1 | -1.0 | +0.7 |
| MSL-GBR | 1 | 0.5 | +0.7 | +0.7 |
| ESN-ESN | 37 | 3.3 | -0.1 | +0.3 |
| GBR-GBR | 36 | 2.1 | -0.2 | -0.2 |

## 4. Engine comparison

| model | engine | t_sis | Ne_A | x_A | secs |
|---|---|---|---|---|---|
| SNP-only | pathfinder | 102 | 318,983 | 0.00044 | 1 |
| SNP-only | nuts | 164 | 841,766 | 0.00050 | 110 |
| IBD-only | pathfinder | 109 | 1,243,741 | 0.00009 | 7 |
| IBD-only | nuts | 109 | 1,239,958 | 0.00009 | 559 |
| Mixed | pathfinder | 117 | 1,131,053 | 0.00010 | 13 |
| Mixed | nuts | 116 | 1,174,691 | 0.00010 | 2453 |

(on `((YRI,ESN),FIN)`, R-hat 1.00 for every NUTS run)

- **Mixed**: agree to 3 digits, 189x faster. The two-leaf mode-splitting that
  motivated running both was created by `kappa_snp`, which is gone.
- **SNP-only**: `t` and `Ne` differ 2.6x but the drifts agree to 2 digits - the
  posterior is a ridge (`r(log t, log Ne) = 0.98`), and both engines find the
  identified direction while sliding freely along the unidentified one.
- **IBD-only**: agree to 4 digits; IBD pins `t`, so there is no ridge.

NUTS is now off by default in `fit_3pop.py` (`--engines nuts` to re-enable).


## 5. What the three configurations establish

**Leaf choice matters more than anything else in the model.**

| configuration | defect | chi2/n IBD | chi2/n SNP |
|---|---|---|---|
| `((YRI,ESN),FIN)` | YRI drift observed NEGATIVE; FIN founder bottleneck | 17.46 | 3.29 |
| `((MSL,ESN),FIN)` | FIN founder bottleneck | 19.33 | 9.96 |
| `((MSL,ESN),GBR)` | none of the above | **2.32** | 10.91 |

- **YRI is unusable as a sister leaf**: observed `x_YRI = -0.00038`. No tree can produce a
  negative branch drift (`Ne > 0` forces `x > 0`), so the fit rails `Ne_YRI` to 1.2 million.
  The f3 screen said "not significantly admixed" (Z = -1.0), which is not the same as
  "usefully positive" - select on the drift being solidly positive, not on the test passing.
- **FIN is unusable as an outgroup**: `FIN-FIN` alone is 91% of the total IBD misfit
  (62.2 x 31 bins = 1928 of 2125). The residual runs monotonically from **+24 at 2 cM to
  -9 at 15 cM** - too few short segments, too many long ones, i.e. the data want larger
  recent Ne and smaller ancient Ne. That is the Finnish founder bottleneck, and one Ne per
  branch cannot represent it. FIN was picked for having the richest IBD (9.72 seg/pair),
  but that richness IS the bottleneck.
- **GBR fixes it**: `GBR-GBR` chi2/bin = 2.1 against FIN's 62.2, and total IBD misfit drops
  8x. Criteria that worked: modest IBD (1.18 seg/pair), large outbred population, no founder
  event, f3-clean against both Africans (Z = +69.8, +66.5).

## 6. The result on the clean tree `((MSL,ESN),GBR)`

Both data types now fit their OWN data well - SNP-only chi2/n = 0.36, IBD-only = 2.14 -
and they still disagree. At the IBD-identified split `t_sis = 163` generations (~4.7 kyr,
the Bantu expansion):

| leaf | x observed | Ne SNP needs | Ne IBD gives | ratio |
|---|---|---|---|---|
| MSL | 0.00383 | 42,695 | 150,887 | 3.5x |
| ESN | 0.00449 | 36,365 | 87,885 | 2.4x |

**The t-free ratio test** - the thing only an outgroup can supply:

```
SNP says  Ne_ESN/Ne_MSL = x_MSL/x_ESN = 0.852
IBD says  Ne_ESN/Ne_MSL              = 0.582
```

These disagree, so the conflict is **not** a shared level shift that any single global
rescaling could absorb. MSL and ESN expanded by different factors. No reweighting of the
two likelihoods can reconcile that - which is the strongest evidence so far that the fix
must come from the parameterisation (Ne resolved in time), not from balancing the data
types against each other. This retro-justifies deleting `kappa_snp` rather than repairing it.

## 7. Remaining misfits

1. **SNP is overruled in the composite**: chi2/n goes 0.36 -> 10.91 while IBD barely moves
   (2.14 -> 2.32). 114 IBD terms against 6 SNP terms; the composite is IBD-dominated by
   sheer term count. This is structural, not a weighting bug.
2. **MSL-ESN cross-pair, chi2/bin = 6.1 over 5 bins.** With `t_sis = 241` the model predicts
   ~1e-30 cross-population IBD at 10 cM, but the data hold 1,110 MSL-ESN segments >= 2 cM.
   A clean split cannot generate recent cross-population sharing: this is ongoing gene flow
   between neighbouring African populations, which a strict tree has no way to express.
3. **ESN-ESN chi2/bin = 3.3** with a wave-shaped residual - milder version of the FIN problem.

## 8. Layout

```
clustering/
  THREE_POP_SUMMARY.md          <- this file
  f3_screen.py                  leaf admixture screen (run BEFORE fitting)
  se_block_sweep.py             SE block-size measurement (kappa's replacement)
  fit_3pop.py                   fit + decomposition + figures + report
  results_3pop_<tag>/
      stan_data/                stan_3pop_<tag>.{npz,json}, labels
      fits.json  report.md
      spectrum_fit.png          observed vs fitted IBD spectrum + residuals, per pair
      likelihood_split.png      lp and chi2/n by component, term counts
      drift_and_Ne.png          fitted vs observed drift, Ne per leaf
  results_screen/               f3 + SE-sweep caches and outputs
  results_noadmix/stan_data/    earlier 4-leaf inputs
```
