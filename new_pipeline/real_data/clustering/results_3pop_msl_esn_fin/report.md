# Three-leaf fit  ((MSL,ESN),FIN)

engine: pathfinder   nodes: ['MSL', 'ESN', 'FIN', 'sister', 'root']


## Parameters

| model | engine | t_sis | t_root | Ne_MSL | Ne_ESN | Ne_FIN | x_MSL | x_ESN | tau | secs |
|---|---|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | pathfinder | 65 | 465 | 16,733 | 14,670 | 2,013 | 0.00391 | 0.00447 | 0.52 | 1 |
| IBD-only-Nsm | pathfinder | 165 | 359 | 150,472 | 88,032 | 43,648 | 0.00110 | 0.00188 | 0.87 | 4 |
| Mixed-Nsm | pathfinder | 256 | 392 | 146,143 | 86,648 | 44,117 | 0.00175 | 0.00295 | 1.28 | 11 |
| **OBSERVED** | (w_hat) | | | | | | **0.00379** | **0.00454** | | |

## Likelihood decomposition (per draw)

| model | engine | lp_ibd | lp_snp | n_ibd | n_snp | chi2/n ibd | chi2/n snp |
|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | pathfinder | - | +41.9 ± 0.8 | - | 6 | - | 0.31 |
| IBD-only-Nsm | pathfinder | +230.4 ± 1.2 | - | 110 | - | 19.08 | - |
| Mixed-Nsm | pathfinder | +218.0 ± 2.1 | +13.0 ± 1.6 | 110 | 6 | 19.33 | 9.96 |

## Per-pair IBD fit (Mixed, chi2 per observed bin)

| pair | bins | chi2/bin | short-bin resid | long-bin resid |
|---|---|---|---|---|
| MSL-MSL | 35 | 0.9 | +0.3 | +0.1 |
| MSL-ESN | 5 | 7.2 | -1.0 | +0.7 |
| MSL-FIN | 1 | 1.3 | +1.2 | +1.2 |
| ESN-ESN | 37 | 3.3 | -0.2 | +0.3 |
| ESN-FIN | 1 | 0.2 | +0.4 | +0.4 |
| FIN-FIN | 31 | 62.2 | +9.0 | -6.5 |

A monotone sign change from the short-bin to the long-bin column is the signature of Ne changing over time, which constant-per-branch Ne cannot fit.

## Figures

- `spectrum_fit.png` - observed vs fitted IBD spectrum + residuals, per pair
- `likelihood_split.png` - lp and chi2/n by component, and term counts
- `drift_and_Ne.png` - fitted vs observed SNP branch drift, and Ne per leaf
