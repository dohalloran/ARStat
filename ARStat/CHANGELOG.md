# Changelog

## v1.1.1 — focused assay scope

- Focused ARStat on egg hatch, larval development, and survival/mortality assays.
- Removed the motility workflow, sample dataset, templates, benchmarks, tests, and documentation.

## v1.1.0 — multi-group normalized replicate input

- Added optional experimental-group and drug/compound column mapping for normalized XY replicate tables.
- Experimental groups can represent strains, isolates, genetic backgrounds, populations, or treatment groups.
- ARStat fits one curve per drug-by-group combination and retains reference-group fold-resistance calculations.
- Added single-dataset and multi-group normalized replicate templates.
- Normalized input continues to require individual replicate values; precomputed summary statistics are not accepted.

## 2026-07 normalized XY replicate import

- Added CSV/XLSX import for wide XY tables containing one dose column and individual replicate response columns.
- ARStat now calculates mean, standard deviation, and n from replicate values rather than requiring imported summary statistics.
- Added a normalized XY replicate template and vendor-neutral documentation.

## v1.0.0-rc1

Release-candidate package for public software testing and validation.

### Added

- Three supported assay workflows: egg hatch, larval development, and survival/mortality.
- Hookworm-style synthetic example datasets using *Ancylostoma caninum* WMD/KGR and assay-appropriate drugs.
- Simulated benchmark generator with known IC50 values.
- Example-validation script.
- Fold-resistance confidence interval columns.
- Benjamini-Hochberg and Bonferroni adjusted p-values for exploratory pairwise tests.
- Statistical notes and limitations page text.
- Documentation for validation, statistical methods, and workflow comparison.

### Fixed

- Zero-dose controls now display on log-scale plots at a symbolic tick labelled `0`; tick labels are applied after setting the log axis so Matplotlib does not overwrite them.
- Download buttons no longer trigger model refitting.
- Session-state handling protects against missing stored results.
- Deprecated Streamlit dataframe width arguments replaced.


## v1.0.0-rc3 plot display update

- Added a traditional raw-outcome IC50 plot display option, so hatch rate, development rate, and survival can be shown as descending curves while retaining the existing inhibition/mortality-based IC50 fitting.
- Added `survival_fraction` for survival assay display while preserving mortality/affected fraction as the fitted response.
- Included plot-display mode in result staleness detection, so toggling plot style prompts users to rerun instead of showing stale figures.
- Made plot captions conditional on the selected plot display mode.


## v1.0.0-rc2 audit fixes

- Corrected the larval development unit test so developed larvae are treated as the success count and undeveloped larvae as the failure count, matching the app presets and documentation.
- Added test assertions that larval-development inhibition increases at high dose and that fitted top values exceed fitted bottom values for the standard example.
- Improved fold-resistance bootstrap confidence intervals by using 10,000 ratio draws from stored IC50 bootstrap samples.
- Added a fit message when a converged 4PL model has a fitted top below the fitted bottom, which can indicate swapped columns, incorrect assay settings, or poor data quality.
