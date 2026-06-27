# Comparison with GraphPad, Excel/manual workflows, and general dose-response tools

## The problem ARStat addresses

Anthelmintic resistance phenotyping assays often begin with raw counts or scores, but published results usually report normalized dose-response curves, IC50/EC50 values, and resistance ratios. The steps between raw measurements and final figures are frequently performed in Excel, GraphPad Prism, or custom scripts. This can create inconsistencies in response definitions, normalization, IC50 interpretation, multiple-testing treatment, and reporting.

ARStat is designed to standardize that middle step for parasitology assays.

## Comparison table

| Feature | Excel/manual workflow | GraphPad Prism | General R/Python packages | ARStat |
|---|---|---|---|---|
| No-code interface | Yes | Yes | Usually no | Yes |
| Starts from raw egg/larval/alive/dead counts | Manual setup required | Manual setup required | User-defined | Built in |
| Egg hatch-specific response calculation | Manual | Manual | User-defined | Built in |
| Larval development-specific response calculation | Manual | Manual | User-defined | Built in |
| Survival/mortality-specific response calculation | Manual | Manual | User-defined | Built in |
| Motility normalization to zero-dose control within isolate/drug group | Manual | Manual | User-defined | Built in |
| Four-parameter logistic IC50 fitting | Hard to do robustly | Yes | Yes | Yes |
| Bootstrap IC50 confidence intervals | Manual or difficult | Available in some workflows | Package-dependent | Built in |
| Fold-resistance versus susceptible reference | Manual | Manual | User-defined | Built in |
| Fold-resistance confidence intervals | Manual | Manual | User-defined | Built in when bootstrap/CI available |
| Pairwise dose-level tests | Manual | Possible | Possible | Built in, with BH and Bonferroni corrections |
| Transparent assay-specific warnings | No | Limited | User-defined | Built in |
| Downloadable methods text | No | No | No | Built in |
| Downloadable templates and example datasets | No | No | Package-dependent | Built in |
| Reproducible scriptable backend | No | Limited | Yes | Yes |
| Parasitology-specific focus | No | No | No | Yes |

## Positioning

ARStat is not intended to replace GraphPad Prism or general dose-response packages. Instead, it provides an assay-aware, parasitology-specific front end and reporting layer that standardizes how raw anthelmintic resistance phenotyping data are converted into fitted curves, IC50 estimates, fold-resistance ratios, plots, and methods text.
