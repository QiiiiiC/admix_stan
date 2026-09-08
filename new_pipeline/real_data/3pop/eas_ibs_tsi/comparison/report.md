# Eight-model comparison on EAS / IBS / TSI

## Best topology within each variant

| variant | topology | ELBO | dELBO runner-up |
|---|---|---:|---:|
| Poisson, shared Ne | `(((EAS.1,IBS),TSI),EAS.2)` | -523.1 | -75.8 |
| Poisson, separate Ne | `((EAS,IBS.1),(IBS.2,TSI))` | -499.8 | -29.5 |
| Normal, shared Ne | `(((EAS.1,IBS),TSI),EAS.2)` | +4023.2 | -146.7 |
| Normal, separate Ne | `((EAS,IBS.1),(IBS.2,TSI))` | +4112.2 | -42.3 |
| Poisson, shared Ne, grid | `(((EAS.1,IBS),TSI),EAS.2)` | -540.7 | -196.8 |
| Poisson, separate Ne, grid | `((EAS,IBS.1),(IBS.2,TSI))` | -549.4 | -13.8 |
| Normal, shared Ne, grid | `(((EAS.1,IBS),TSI),EAS.2)` | +3994.5 | -88.0 |
| Normal, separate Ne, grid | `((EAS,IBS.1),(IBS.2,TSI))` | +4063.4 | -17.8 |

## Ne-model evidence gains

These differences are valid only within the same likelihood.

| likelihood | comparison | best-to-best gain | same winner? |
|---|---|---:|---:|
| Poisson | simple: separate - shared | +23.3 | no |
| Poisson | grid: separate - shared | -8.7 | no |
| Poisson | shared: grid - simple | -17.6 | yes |
| Poisson | separate: grid - simple | -49.6 | yes |
| Normal | simple: separate - shared | +89.1 | no |
| Normal | grid: separate - shared | +68.9 | no |
| Normal | shared: grid - simple | -28.6 | yes |
| Normal | separate: grid - simple | -48.9 | yes |

![comparison](eight_model_comparison.png)

## Rank consistency

ELBO values are comparable between shared- and separate-Ne models within one likelihood, but not between Poisson and Normal because they are masses/densities for different observed summaries.

| comparison | Spearman rank correlation | top-5 overlap |
|---|---:|---:|
| Poisson simple: shared vs separate | +0.956 | 5/5 |
| Poisson grid: shared vs separate | +0.955 | 5/5 |
| Poisson shared: simple vs grid | +0.974 | 5/5 |
| Poisson separate: simple vs grid | +0.995 | 5/5 |
| Normal simple: shared vs separate | +0.761 | 5/5 |
| Normal grid: shared vs separate | +0.878 | 5/5 |
| Normal shared: simple vs grid | +0.969 | 5/5 |
| Normal separate: simple vs grid | +0.957 | 5/5 |
| Simple shared: Poisson vs Normal | +0.713 | 5/5 |
| Simple separate: Poisson vs Normal | +0.913 | 5/5 |
| Grid shared: Poisson vs Normal | +0.704 | 5/5 |
| Grid separate: Poisson vs Normal | +0.855 | 5/5 |

## Winner diagnostics

| variant | IBD chi2/n | SNP chi2/n | admixture fraction | interpretation |
|---|---:|---:|---:|---|
| Poisson, shared Ne | 1.61 | 2.59 | 0.9929 | collapsed/near-tree |
| Poisson, separate Ne | 1.11 | 1.56 | 0.0055 | collapsed/near-tree |
| Normal, shared Ne | 1.48 | 2.44 | 0.9925 | collapsed/near-tree |
| Normal, separate Ne | 1.03 | 1.41 | 0.0000 | collapsed/near-tree |
| Poisson, shared Ne, grid | 1.53 | 2.94 | 0.9926 | collapsed/near-tree |
| Poisson, separate Ne, grid | 1.09 | 2.13 | 0.0056 | collapsed/near-tree |
| Normal, shared Ne, grid | 1.38 | 1.94 | 0.9922 | collapsed/near-tree |
| Normal, separate Ne, grid | 0.96 | 3.32 | 0.0008 | collapsed/near-tree |

## Recent Ne in grid-model winners

The ratio is Ne(0-5 generations) / Ne(5-10 generations); values above one indicate growth toward the present.

| variant | component | population | Ne 0-5 | Ne 5-10 | ratio | first event |
|---|---|---|---:|---:|---:|---:|
| Poisson, shared Ne, grid | shared | EAS | 4639564 | 3779480 | 1.23 | 129.0 |
| Poisson, shared Ne, grid | shared | IBS | 1659363 | 1201337 | 1.38 | 129.0 |
| Poisson, shared Ne, grid | shared | TSI | 985491 | 691514 | 1.43 | 129.0 |
| Poisson, separate Ne, grid | IBD | EAS | 4827530 | 3883861 | 1.24 | 72.1 |
| Poisson, separate Ne, grid | IBD | IBS | 1371458 | 1298804 | 1.06 | 72.1 |
| Poisson, separate Ne, grid | IBD | TSI | 1027199 | 807863 | 1.27 | 72.1 |
| Poisson, separate Ne, grid | SNP | EAS | 1190 | 1095 | 1.09 | 72.1 |
| Poisson, separate Ne, grid | SNP | IBS | 89132 | 89786 | 0.99 | 72.1 |
| Poisson, separate Ne, grid | SNP | TSI | 127689 | 135451 | 0.94 | 72.1 |
| Normal, shared Ne, grid | shared | EAS | 6665628 | 5451659 | 1.22 | 130.0 |
| Normal, shared Ne, grid | shared | IBS | 2993299 | 2036788 | 1.47 | 130.0 |
| Normal, shared Ne, grid | shared | TSI | 681492 | 577692 | 1.18 | 130.0 |
| Normal, separate Ne, grid | IBD | EAS | 4479374 | 3920338 | 1.14 | 76.2 |
| Normal, separate Ne, grid | IBD | IBS | 1444785 | 1229511 | 1.18 | 76.2 |
| Normal, separate Ne, grid | IBD | TSI | 961770 | 742528 | 1.30 | 76.2 |
| Normal, separate Ne, grid | SNP | EAS | 1134 | 1187 | 0.96 | 76.2 |
| Normal, separate Ne, grid | SNP | IBS | 110513 | 108422 | 1.02 | 76.2 |
| Normal, separate Ne, grid | SNP | TSI | 112389 | 112680 | 1.00 | 76.2 |

## Interpretation

Poisson and Normal ELBO levels are not subtracted from each other: the two likelihoods are densities/masses for different summaries and therefore use different base measures and units. Compare their topology ranks, residuals, and shared-versus-separate-Ne gains instead.

The data contain **116/222 empty unique pair-by-bin cells**. The winning grid separate-Ne Normal fit places 68/222 unique cells at its hard `1e-12` theory-SE floor. Consequently, its all-bin CLT likelihood is being used most aggressively exactly where the CLT is least justified.

The Poisson winner has admixture fraction **0.0056**. It is therefore best read as a near-tree model with an extra branch breakpoint if the fraction remains near a boundary, not as evidence for substantial admixture.

The grid comparison directly tests whether 0-10 generation growth explains the topology preference. Rank agreement and the winner-fraction diagnostics above should be considered together; a high ELBO for boundary admixture is evidence for remaining Ne misspecification rather than for gene flow.

## Conclusion

Poisson, shared Ne selects topology 5 (`(((EAS.1,IBS),TSI),EAS.2)`), f=0.9929, collapsing to tree 3 (`((IBS,TSI),EAS)`). Poisson, separate Ne selects topology 15 (`((EAS,IBS.1),(IBS.2,TSI))`), f=0.0055, collapsing to tree 3 (`((IBS,TSI),EAS)`). Normal, shared Ne selects topology 5 (`(((EAS.1,IBS),TSI),EAS.2)`), f=0.9925, collapsing to tree 3 (`((IBS,TSI),EAS)`). Normal, separate Ne selects topology 15 (`((EAS,IBS.1),(IBS.2,TSI))`), f=0.0000, collapsing to tree 3 (`((IBS,TSI),EAS)`). Poisson, shared Ne, grid selects topology 5 (`(((EAS.1,IBS),TSI),EAS.2)`), f=0.9926, collapsing to tree 3 (`((IBS,TSI),EAS)`). Poisson, separate Ne, grid selects topology 15 (`((EAS,IBS.1),(IBS.2,TSI))`), f=0.0056, collapsing to tree 3 (`((IBS,TSI),EAS)`). Normal, shared Ne, grid selects topology 5 (`(((EAS.1,IBS),TSI),EAS.2)`), f=0.9922, collapsing to tree 3 (`((IBS,TSI),EAS)`). Normal, separate Ne, grid selects topology 15 (`((EAS,IBS.1),(IBS.2,TSI))`), f=0.0008, collapsing to tree 3 (`((IBS,TSI),EAS)`).

All eight winners collapse to the same non-admixture topology **3**. The raw graph labels differ between shared- and separate-Ne parameterizations, but their boundary fractions make them equivalent at the population-tree level.

The grid separate-Ne Poisson winner has IBD chi2/n **1.09** but fraction **0.0056**. Its good count calibration does not turn that boundary edge into admixture evidence; the graph is acting as a tree with an additional branch-specific Ne change.

The grid separate-Ne Normal winner has fraction **0.0008**, but its IBD chi2/n is **0.96** and 68/222 cells use the variance floor. Its topology result is therefore sensitivity evidence, not a reliable resolution of the Poisson result.

The next misspecification test should use fixed absolute Ne breakpoints beyond generation 10 and truncate each branch trajectory at its event time. That lets events occur inside the grid, avoiding an artificial lower bound on the first event while testing whether the near-tree admixture edge disappears.
