# Three-leaf fit  ((MSL,ESN),GBR)

engine: pathfinder   nodes: ['MSL', 'ESN', 'GBR', 'sister', 'root']


## Parameters

| model | engine | t_sis | t_root | Ne_MSL | Ne_ESN | Ne_GBR | x_MSL | x_ESN | tau | secs |
|---|---|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | pathfinder | 46 | 547 | 11,027 | 10,655 | 2,667 | 0.00414 | 0.00428 | 0.36 | 1 |
| IBD-only-Nsm | pathfinder | 163 | 400 | 150,887 | 87,885 | 170,186 | 0.00108 | 0.00186 | 0.90 | 4 |
| Mixed-Nsm | pathfinder | 241 | 479 | 147,468 | 86,924 | 168,010 | 0.00164 | 0.00278 | 1.20 | 12 |
| **OBSERVED** | (w_hat) | | | | | | **0.00383** | **0.00449** | | |

## Likelihood decomposition (per draw)

| model | engine | lp_ibd | lp_snp | n_ibd | n_snp | chi2/n ibd | chi2/n snp |
|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | pathfinder | - | +41.6 ± 0.9 | - | 6 | - | 0.36 |
| IBD-only-Nsm | pathfinder | +1294.2 ± 2.3 | - | 114 | - | 2.14 | - |
| Mixed-Nsm | pathfinder | +1285.0 ± 1.6 | +10.0 ± 1.7 | 114 | 6 | 2.32 | 10.91 |

## Per-pair IBD fit (Mixed, chi2 per observed bin)

| pair | bins | chi2/bin | short-bin resid | long-bin resid |
|---|---|---|---|---|
| MSL-MSL | 35 | 0.9 | +0.4 | +0.1 |
| MSL-ESN | 5 | 6.1 | -1.0 | +0.7 |
| MSL-GBR | 1 | 0.5 | +0.7 | +0.7 |
| ESN-ESN | 37 | 3.3 | -0.1 | +0.3 |
| GBR-GBR | 36 | 2.1 | -0.2 | -0.2 |

A monotone sign change from the short-bin to the long-bin column is the signature of Ne changing over time, which constant-per-branch Ne cannot fit.

## Figures

- `spectrum_fit.png` - observed vs fitted IBD spectrum + residuals, per pair
- `likelihood_split.png` - lp and chi2/n by component, and term counts
- `drift_and_Ne.png` - fitted vs observed SNP branch drift, and Ne per leaf
