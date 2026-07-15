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

### Survival/mortality assay

Mortality fraction is calculated as:

```text
dead / (dead + alive)
```

The modeled drug effect is mortality/affected fraction.


## Dose-response model

ARStat fits an increasing four-parameter logistic model on the log10 dose scale:

```text
response = bottom + (top - bottom) / (1 + 10^((logIC50 - log10(dose)) * hill))
```

IC50 is returned on the original dose scale as `10^logIC50`.

Important: IC50 is the dose corresponding to the midpoint between the fitted lower and upper asymptotes. If the fitted top response is below 100%, the IC50 is not necessarily the dose producing an absolute 50% response.

## Confidence intervals

When enabled, ARStat estimates IC50 confidence intervals by nonparametric bootstrap resampling of rows within each fitted group. A confidence interval is reported when at least 25% of bootstrap fits converge and at least 20 bootstrap estimates are available.

## Fold-resistance

Fold-resistance versus the reference isolate is calculated as:

```text
IC50_test / IC50_reference
```

When bootstrap IC50 samples are available for both test and reference groups, ARStat estimates a bootstrap confidence interval for the ratio. Otherwise, if IC50 confidence limits are available, it provides an approximate log-scale propagated interval.

## Pairwise dose-level tests

For count-based assays, ARStat performs Fisher exact tests at each dose after pooling replicate counts within each group. Raw p-values are accompanied by Benjamini-Hochberg and Bonferroni adjusted p-values.

For normalized replicate inputs, ARStat performs Mann-Whitney U tests by default at each dose.

## Limitations

- Pairwise dose-level tests are exploratory.
- Fisher exact tests pool replicate counts and do not model replicate-to-replicate overdispersion.
- For highly variable wells, pooled Fisher tests may be anticonservative.
- Zero-dose controls are included in fitting and summaries, but on log-scale plots they are displayed at a symbolic left-edge tick labelled `0`.
- Current ARStat models are not full mixed-effects or beta-binomial models.
