# Changelog

## v1.0.0-rc1

Release-candidate package for public software testing and validation.

### Added

- Four supported assay workflows: egg hatch, larval development, survival/mortality, and motility.
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

- Added a traditional raw-outcome IC50 plot display option, so hatch rate, development rate, survival, and motility can be shown as descending curves while retaining the existing inhibition/mortality-based IC50 fitting.
- Added `survival_fraction` for survival assay display while preserving mortality/affected fraction as the fitted response.
- Included plot-display mode in result staleness detection, so toggling plot style prompts users to rerun instead of showing stale figures.
- Made plot captions conditional on the selected plot display mode.


## v1.0.0-rc2 audit fixes

- Corrected the larval development unit test so developed larvae are treated as the success count and undeveloped larvae as the failure count, matching the app presets and documentation.
- Added test assertions that larval-development inhibition increases at high dose and that fitted top values exceed fitted bottom values for the standard example.
- Added warnings when motility scores exceed matched zero-dose controls; negative inhibition is retained in an unclipped helper column and clipped to 0 only for the fitted response.
- Improved fold-resistance bootstrap confidence intervals by using 10,000 ratio draws from stored IC50 bootstrap samples.
- Added a fit message when a converged 4PL model has a fitted top below the fitted bottom, which can indicate swapped columns, incorrect assay settings, or poor data quality.
