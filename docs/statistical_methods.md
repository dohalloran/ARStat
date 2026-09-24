# Statistical methods implemented in ARStat

## Response calculations

### Egg hatch assay

Raw hatch fraction is calculated as:

```text
L1 / (L1 + eggs)
```

The modeled drug effect is hatch inhibition:

```text
1 - hatch_fraction
```

### Larval development assay

Raw development fraction is calculated as:

```text
developed / (developed + undeveloped)
```

The modeled drug effect is development inhibition:

```text
1 - development_fraction
```

### Motility assay

For raw continuous activity measurements, ARStat calculates a separate zero-dose mean for every fitted group, normally each drug-by-experimental-group combination:

```text
relative_motility = activity_measurement / mean_zero_dose_activity
motility_inhibition = 1 - relative_motility
```

The fitted response is motility inhibition. Relative motility is retained for the traditional descending plot and row-level exports.

Values can instead be imported as already normalized percentages (`0–100`) or fractions (`0–1`). The user specifies whether imported values represent retained motility or motility inhibition. Values outside the biological 0–1 interval are flagged but preserved because control normalization, background correction, hypermotility, or ordinary replicate variation can produce observations below 0% or above 100%. The fitted 4PL asymptotes remain constrained to 0–1.

## Dose-response model

ARStat fits an increasing four-parameter logistic model on the log10 dose scale:

```text
response = bottom + (top - bottom) / (1 + 10^((logIC50 - log10(dose)) * hill))
```

IC50 is returned on the original dose scale as `10^logIC50`.

Fitting details:

- Parameters are estimated by bounded nonlinear least squares with SciPy `curve_fit` (trust-region reflective algorithm). Bounds are 0–1 for `bottom` and `top`, 0.05–10 for `hill`, and log10 of the smallest positive dose − 3 to log10 of the largest positive dose + 3 for `logIC50`. An estimate that sits on a bound (for example, a Hill slope of 10 when no concentration falls within the response transition) should be interpreted cautiously.
- Zero-dose observations keep dose = 0 in all exported tables. Because log10(0) is undefined, the model evaluates them at one-tenth of the smallest positive dose in the fitted group.
- Starting values are taken from the data, and the fit is retried with Hill-slope starting values of 1, 0.5, 2, and 4 if the first attempt fails.

IC50 is the dose corresponding to the midpoint between the fitted lower and upper asymptotes. If the fitted top response is below 100%, the IC50 is not necessarily the dose producing an absolute 50% response.

For motility, the result is a motility-inhibition IC50 under the imported activity definition and normalization. It should not be interpreted automatically as a lethal-concentration endpoint.

## Confidence intervals

When enabled, ARStat estimates percentile 95% IC50 confidence intervals by nonparametric bootstrap resampling of rows within each fitted group. A confidence interval is reported when at least 25% of bootstrap fits converge and at least 20 bootstrap estimates are available.

These intervals are approximate. Resampling is not stratified by dose, and it does not represent plate-, day-, or experiment-level hierarchy. `scripts/simulate_bootstrap_coverage.py` simulates egg hatch data with a known IC50. With 120 simulated data sets per design, nominal 95% intervals contained the true IC50 in 91% of data sets for a design like the BCR validation experiments (7 concentrations × 3 wells, with well-to-well overdispersion) and in 92% for the bundled benchmark design (8 concentrations × 6 wells, binomial). The Monte Carlo standard error is about 2.5 percentage points, so coverage is slightly below nominal. Results are in `benchmarks/bootstrap_coverage_summary.csv`.

Count-based curve fitting can use total counts as fitting weights. Continuous motility and normalized XY responses are fitted without count weights. Normalized replicate observations are fitted as imported and are not clipped to 0–1.

## Fold-resistance

Fold-resistance versus the reference group is calculated as:

```text
IC50_test / IC50_reference
```

When bootstrap IC50 samples are available for both test and reference groups, ARStat estimates a bootstrap confidence interval for the ratio. Otherwise, if IC50 confidence limits are available, it provides an approximate log-scale propagated interval.

## Pairwise dose-level tests

For count-based assays, ARStat performs Fisher exact tests at each dose after pooling replicate counts within each group. Raw p-values are accompanied by Benjamini-Hochberg and Bonferroni adjusted p-values.

For motility and normalized replicate inputs, ARStat performs replicate-level Mann-Whitney U tests by default at each dose. These tests compare distributions at individual doses and are separate from the nonlinear curve fit. The output reports the number of replicates in each group and `exact_min_p`, the smallest two-sided P value an exact test can return for those sample sizes (2 / C(n1 + n2, n1); 0.10 for 3 versus 3 replicates). When values are tied, SciPy uses a normal approximation, which is unreliable at such small sample sizes.

## Limitations

- Pairwise dose-level tests are exploratory.
- Fisher exact tests pool replicate counts and do not model replicate-to-replicate overdispersion.
- Mann-Whitney tests do not account for repeated measures, plate effects, nested experiments, or other dependence structures.
- Raw motility normalization requires a representative zero-dose control within every fitted group.
- Normalized values below 0% or above 100% are preserved for all normalized XY workflows; complementary raw-outcome or inhibition values may therefore also fall outside 0–1. These flagged rows should be reviewed rather than assumed to be errors.
- A manual ordinal motility score is treated numerically and should only be used when scoring is consistent and scientifically justified.
- Zero-dose controls are included in fitting and summaries, but on log-scale plots they are displayed at a symbolic left-edge tick labelled `0`.
- Current ARStat models are not full mixed-effects, beta-binomial, repeated-measures, hormesis, or time-series models.
