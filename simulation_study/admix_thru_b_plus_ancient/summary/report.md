# Recent loop + ancient admixture (Nfixed): IBD-only vs SNP-only vs Mixed

Truth: loop on a 20→60→100 gen (α = 0.5), ancient admixture on c at 500 gen (70%/30%), merges 700/900/1000, haploid Ne 10,000 everywhere.  100 subsamples per genome length, blocks of 25 cM drawn without replacement from a pool of 100; subsamples overlap, so shares are sensitivities to genome selection, not independent-replicate rates.

ΔELBO = ELBO(T_true) − ELBO(alternative), each from three pooled Pathfinder starts on normalised models. 'wins' = share of subsamples with ΔELBO > 0; 'both' = share in which T_true beats both alternatives. SNP-only is fitted once per subsample and repeated under both IBD sources.

Fits: 12000 attempted, 12000 ok, 748 pass the ESS/Pareto-k screen.

## IBD source: true

| cM | model | median Δ recent | wins recent | median Δ ancient | wins ancient | both | n |
|---:|---|---:|---:|---:|---:|---:|---:|
| 50 | IBD-only | +2.8 | 84% | +0.8 | 82% | 68% | 100 |
| 50 | SNP-only | -4.7 | 1% | +75.2 | 100% | 1% | 100 |
| 50 | Mixed | +4.2 | 83% | +76.5 | 100% | 83% | 100 |
| 100 | IBD-only | +9.3 | 98% | +0.9 | 86% | 85% | 100 |
| 100 | SNP-only | -5.4 | 3% | +145.0 | 100% | 3% | 100 |
| 100 | Mixed | +11.1 | 99% | +146.0 | 100% | 99% | 100 |
| 150 | IBD-only | +15.3 | 100% | +0.8 | 88% | 88% | 100 |
| 150 | SNP-only | -5.7 | 1% | +217.0 | 100% | 1% | 100 |
| 150 | Mixed | +17.8 | 100% | +218.1 | 100% | 100% | 100 |
| 200 | IBD-only | +22.6 | 100% | +0.7 | 86% | 86% | 100 |
| 200 | SNP-only | -5.8 | 5% | +280.1 | 100% | 5% | 100 |
| 200 | Mixed | +26.9 | 100% | +278.6 | 100% | 100% | 100 |
| 300 | IBD-only | +34.4 | 100% | +0.7 | 87% | 87% | 100 |
| 300 | SNP-only | -5.3 | 5% | +415.7 | 100% | 5% | 100 |
| 300 | Mixed | +38.8 | 100% | +415.9 | 100% | 100% | 100 |
| 500 | IBD-only | +61.0 | 100% | +0.7 | 84% | 84% | 100 |
| 500 | SNP-only | -4.8 | 16% | +691.7 | 100% | 16% | 100 |
| 500 | Mixed | +68.9 | 100% | +693.3 | 100% | 100% | 100 |
| 750 | IBD-only | +96.6 | 100% | +0.8 | 88% | 88% | 100 |
| 750 | SNP-only | -3.7 | 15% | +1034.4 | 100% | 15% | 100 |
| 750 | Mixed | +107.6 | 100% | +1033.3 | 100% | 100% | 100 |
| 1000 | IBD-only | +133.8 | 100% | +0.6 | 78% | 78% | 100 |
| 1000 | SNP-only | -1.4 | 34% | +1369.7 | 100% | 34% | 100 |
| 1000 | Mixed | +148.6 | 100% | +1366.5 | 100% | 100% | 100 |

## IBD source: hapibd

| cM | model | median Δ recent | wins recent | median Δ ancient | wins ancient | both | n |
|---:|---|---:|---:|---:|---:|---:|---:|
| 50 | IBD-only | +3.6 | 78% | +0.6 | 89% | 70% | 100 |
| 50 | SNP-only | -4.7 | 1% | +75.2 | 100% | 1% | 100 |
| 50 | Mixed | +4.6 | 89% | +75.4 | 100% | 89% | 100 |
| 100 | IBD-only | +10.7 | 98% | +0.7 | 86% | 84% | 100 |
| 100 | SNP-only | -5.4 | 3% | +145.0 | 100% | 3% | 100 |
| 100 | Mixed | +12.1 | 99% | +144.3 | 100% | 99% | 100 |
| 150 | IBD-only | +15.5 | 100% | +0.6 | 82% | 82% | 100 |
| 150 | SNP-only | -5.7 | 1% | +217.0 | 100% | 1% | 100 |
| 150 | Mixed | +16.7 | 100% | +215.6 | 100% | 100% | 100 |
| 200 | IBD-only | +24.9 | 100% | +0.5 | 74% | 74% | 100 |
| 200 | SNP-only | -5.8 | 5% | +280.1 | 100% | 5% | 100 |
| 200 | Mixed | +27.5 | 100% | +277.4 | 100% | 100% | 100 |
| 300 | IBD-only | +36.0 | 100% | +0.5 | 72% | 72% | 100 |
| 300 | SNP-only | -5.3 | 5% | +415.7 | 100% | 5% | 100 |
| 300 | Mixed | +39.4 | 100% | +414.9 | 100% | 100% | 100 |
| 500 | IBD-only | +62.6 | 100% | +0.7 | 80% | 80% | 100 |
| 500 | SNP-only | -4.8 | 16% | +691.7 | 100% | 16% | 100 |
| 500 | Mixed | +65.8 | 100% | +691.6 | 100% | 100% | 100 |
| 750 | IBD-only | +95.7 | 100% | +0.9 | 80% | 80% | 100 |
| 750 | SNP-only | -3.7 | 15% | +1034.4 | 100% | 15% | 100 |
| 750 | Mixed | +103.1 | 100% | +1032.5 | 100% | 100% | 100 |
| 1000 | IBD-only | +132.2 | 100% | +1.1 | 81% | 81% | 100 |
| 1000 | SNP-only | -1.4 | 34% | +1369.7 | 100% | 34% | 100 |
| 1000 | Mixed | +142.7 | 100% | +1364.2 | 100% | 100% | 100 |

Figures: `figures/win_rate_*`, `figures/delbo_*`, `figures/parameters_*` (one per IBD source).
Tables: `contrasts.csv` (per subsample), `fits.csv`, `parameter_estimates.csv`.
