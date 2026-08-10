# Three-leaf fit  ((YRI,ESN),FIN)

engine: pathfinder   nodes: ['YRI', 'ESN', 'FIN', 'sister', 'root']


## Parameters

| model | engine | t_sis | t_root | Ne_YRI | Ne_ESN | Ne_FIN | x_YRI | x_ESN | tau | secs |
|---|---|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | pathfinder | 102 | 447 | 318,983 | 47,791 | 1,739 | 0.00044 | 0.00226 | 0.74 | 1 |
| IBD-only-Nsm | pathfinder | 109 | 367 | 1,243,741 | 111,788 | 43,529 | 0.00009 | 0.00098 | 1.25 | 7 |
| Mixed-Nsm | pathfinder | 117 | 2792 | 1,131,053 | 113,220 | 43,493 | 0.00010 | 0.00103 | 1.18 | 13 |
| **OBSERVED** | (w_hat) | | | | | | **-0.00038** | **0.00295** | | |

## Likelihood decomposition (per draw)

| model | engine | lp_ibd | lp_snp | n_ibd | n_snp | chi2/n ibd | chi2/n snp |
|---|---|---|---|---|---|---|---|
| SNP-only-Nsm | pathfinder | - | +38.3 ± 2.0 | - | 6 | - | 1.56 |
| IBD-only-Nsm | pathfinder | +441.4 ± 1.7 | - | 122 | - | 17.42 | - |
| Mixed-Nsm | pathfinder | +438.1 ± 2.4 | +33.1 ± 0.6 | 122 | 6 | 17.46 | 3.29 |

## Per-pair IBD fit (Mixed, chi2 per observed bin)

| pair | bins | chi2/bin | short-bin resid | long-bin resid |
|---|---|---|---|---|
| YRI-YRI | 27 | 1.0 | +0.2 | +0.0 |
| YRI-ESN | 25 | 2.1 | +1.9 | +0.8 |
| YRI-FIN | 1 | 0.5 | +0.7 | +0.7 |
| ESN-ESN | 37 | 3.2 | +0.8 | +0.6 |
| ESN-FIN | 1 | 0.5 | +0.7 | +0.7 |
| FIN-FIN | 31 | 62.2 | +8.9 | -6.6 |

A monotone sign change from the short-bin to the long-bin column is the signature of Ne changing over time, which constant-per-branch Ne cannot fit.

## Figures

- `spectrum_fit.png` - observed vs fitted IBD spectrum + residuals, per pair
- `likelihood_split.png` - lp and chi2/n by component, and term counts
- `drift_and_Ne.png` - fitted vs observed SNP branch drift, and Ne per leaf
