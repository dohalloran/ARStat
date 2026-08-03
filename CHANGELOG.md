# Changelog

## v1.2.1 — normalized input and column-mapping safeguards

- Retain normalized replicate values below 0% or above 100% instead of clipping them before curve fitting.
- Improved automatic percentage-versus-fraction detection so a single modest outlier does not rescale an entire fractional dataset.
- Display automatic scale detection as an informational note rather than a warning in the Streamlit interface.
- Added a warning when only one normalized replicate column is selected.
- Added conservative assay-column aliases, including an LDA suggestion of `L3` as developed and `L1` as undeveloped, while requiring user confirmation.
- Prevented the same count column from being mapped to both assay outcomes.
- Elevated decreasing inhibition/mortality fits to an error requiring review and explained why complemented responses can retain a similar IC50 midpoint.
- Clarified that fewer than four unique dose levels prevents 4PL IC50 calculation rather than merely making it unstable.
- Automatically remove entirely empty upload columns, including blank `Unnamed:` spreadsheet-export columns.
- Added seven regression tests for these safeguards; the suite now contains 23 tests.

## v1.2.0 — motility IC50 workflow restored

- Restored motility as a fully supported fourth assay workflow while retaining survival/mortality.
- Added raw continuous motility/activity input with normalization to the mean zero-dose control within each fitted group.
- Added support for pre-normalized motility values supplied as 0–100 percentages or 0–1 fractions.
- Added explicit retained-motility versus motility-inhibition response direction handling.
- Added descending relative-motility plotting while fitting the increasing motility-inhibition response.
- Added replicate-level Mann-Whitney dose comparisons for motility data.
- Added a motility template, synthetic example dataset, simulated benchmark, validation outputs, tests, and documentation.
- Updated application and citation version metadata to v1.2.0.

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
