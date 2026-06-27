# ARStat validation plan

## Goal

Validate that ARStat correctly transforms raw anthelmintic phenotyping data into standardized dose-response summaries, IC50 estimates, fold-resistance ratios, adjusted exploratory p-values, and publication-ready outputs.

## Validation layer 1: Unit tests

The included test suite checks:

1. Four-parameter logistic model behavior.
2. Count-response calculations for egg hatch, larval development, and survival/mortality assays.
3. Motility normalization within each isolate/drug group.
4. IC50 fitting and fold-resistance calculations.
5. Multiple-testing correction columns.
6. All four bundled hookworm-style example datasets.

Run:

```bash
pytest -q
```

## Validation layer 2: Simulated benchmark datasets with known IC50 values

The script `scripts/generate_simulated_benchmarks.py` creates synthetic datasets for all four assay workflows with known IC50 values. It then runs the same ARStat core functions used by the web app and reports the fitted IC50, expected IC50, absolute error, and percent error.

Run:

```bash
python scripts/generate_simulated_benchmarks.py
```

Outputs are written to `benchmarks/`:

- `simulated_egg_hatch.csv`
- `simulated_larval_development.csv`
- `simulated_survival.csv`
- `simulated_motility.csv`
- `benchmark_fit_accuracy_summary.csv`
- `benchmark_fold_resistance_summary.csv`

## Validation layer 3: Bundled example datasets

The script `scripts/run_all_examples.py` runs each bundled hookworm-style example through the full analysis workflow and exports summary tables.

Run:

```bash
python scripts/run_all_examples.py
```

Outputs are written to `validation_outputs/`.

## Validation layer 4: Real hookworm case studies

For external validation, at least one real dataset should be added if possible. The strongest case would include susceptible and resistant *Ancylostoma caninum* isolates, for example WMD and KGR, tested using one or more assay types.

Recommended real-data case studies:

1. Thiabendazole egg hatch assay.
2. Ivermectin larval development assay.
3. Ivermectin survival/mortality assay.
4. Ivermectin motility assay, if quantitative motility data are available.

If real data are not ready for public release, include synthetic demonstration data in the software repository and cite the real data as unavailable or available upon request only if journal policy permits. A public software paper is stronger if at least one real de-identified example dataset is included.
