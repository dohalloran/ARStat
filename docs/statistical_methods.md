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

Values can instead be imported as already normalized percentages (`0–100`) or fractions (`0–1`). The user specifies whether imported values represent retained motility or motility inhibition. Retained-motility values above 1 are flagged but preserved because they may represent hypermotility relative to the control mean or ordinary replicate variation. Imported inhibition values outside the biological 0–1 interval are flagged and clipped.

### Survival/mortality assay

Mortality fraction is calculated as:

```text
dead / (dead + alive)
```

The modeled drug effect is mortality/affected fraction. Survival fraction is retained for the traditional descending plot.

## Dose-response model

ARStat fits an increasing four-parameter logistic model on the log10 dose scale:

```text
response = bottom + (top - bottom) / (1 + 10^((logIC50 - log10(dose)) * hill))
```

IC50 is returned on the original dose scale as `10^logIC50`.

IC50 is the dose corresponding to the midpoint between the fitted lower and upper asymptotes. If the fitted top response is below 100%, the IC50 is not necessarily the dose producing an absolute 50% response.

For motility, the result is a motility-inhibition IC50 under the imported activity definition and normalization. It is not interchangeable with a survival or lethal-concentration endpoint.

## Confidence intervals

When enabled, ARStat estimates IC50 confidence intervals by nonparametric bootstrap resampling of rows within each fitted group. A confidence interval is reported when at least 25% of bootstrap fits converge and at least 20 bootstrap estimates are available.

Count-based curve fitting can use total counts as fitting weights. Continuous motility and normalized XY responses are fitted without count weights.

## Fold-resistance

Fold-resistance versus the reference group is calculated as:

```text
IC50_test / IC50_reference
```

When bootstrap IC50 samples are available for both test and reference groups, ARStat estimates a bootstrap confidence interval for the ratio. Otherwise, if IC50 confidence limits are available, it provides an approximate log-scale propagated interval.

## Pairwise dose-level tests

For count-based assays, ARStat performs Fisher exact tests at each dose after pooling replicate counts within each group. Raw p-values are accompanied by Benjamini-Hochberg and Bonferroni adjusted p-values.

For motility and normalized replicate inputs, ARStat performs replicate-level Mann-Whitney U tests by default at each dose. These tests compare distributions at individual doses and are separate from the nonlinear curve fit.

## Limitations

- Pairwise dose-level tests are exploratory.
- Fisher exact tests pool replicate counts and do not model replicate-to-replicate overdispersion.
- Mann-Whitney tests do not account for repeated measures, plate effects, nested experiments, or other dependence structures.
- Raw motility normalization requires a representative zero-dose control within every fitted group.
- Retained-motility values above 100% are preserved; negative fitted inhibition observations can therefore occur at low doses. These flagged rows should be reviewed rather than assumed to be errors.
- A manual ordinal motility score is treated numerically and should only be used when scoring is consistent and scientifically justified.
- Zero-dose controls are included in fitting and summaries, but on log-scale plots they are displayed at a symbolic left-edge tick labelled `0`.
- Current ARStat models are not full mixed-effects, beta-binomial, repeated-measures, hormesis, or time-series models.
