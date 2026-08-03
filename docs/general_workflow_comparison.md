# Comparison with manual and general dose-response workflows

Anthelmintic resistance assays often begin with raw counts, continuous motility measurements, or normalized replicate responses. Published analyses usually report dose-response curves, IC50 values, confidence intervals, and resistance ratios. The steps between the raw assay and the final figure are frequently distributed across spreadsheets, tracking-software exports, general curve-fitting applications, or custom scripts.

ARStat provides a standardized, assay-aware workflow that accepts either:

1. raw assay measurements, including count endpoints or continuous motility/activity measurements; or
2. a normalized XY replicate table containing one dose column, individual replicate response columns, and optional experimental-group and drug columns.

For raw motility data, ARStat normalizes activity to the mean zero-dose control within each fitted group and calculates motility inhibition. For normalized XY data, ARStat calculates dose-level means, standard deviations, and sample sizes from the imported replicates. It does not require precomputed summary statistics.

| Feature | Manual spreadsheet workflow | General curve-fitting software | ARStat |
|---|---:|---:|---:|
| Raw assay templates | Partial | Partial | Yes |
| Assay-specific preprocessing | Partial | Partial | Yes |
| Egg hatch workflow | Partial | Partial | Yes |
| Larval development workflow | Partial | Partial | Yes |
| Motility workflow with control normalization | Partial | Partial | Yes |
| Survival/mortality workflow | Partial | Partial | Yes |
| Wide XY replicate import | Yes | Yes | Yes |
| Automatic IC50 and confidence intervals | Partial | Yes | Yes |
| Fold-resistance calculation | Partial | Partial | Yes |
| Adjusted dose-level tests | No | Partial | Yes |
| Reproducible export files and methods text | Partial | Partial | Yes |
