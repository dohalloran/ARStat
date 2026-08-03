# ARStat validation plan

## Goal

Validate that ARStat correctly transforms raw anthelmintic phenotyping data into standardized dose-response summaries, IC50 estimates, fold-resistance ratios, adjusted exploratory p-values, and reproducible outputs.

## Validation layer 1: Unit tests

The included test suite checks:

1. Four-parameter logistic model behavior.
2. Count-response calculations for egg hatch, larval development, and survival/mortality assays.
3. Raw motility normalization, motility-inhibition calculation, and zero-dose control validation.
4. IC50 fitting and fold-resistance calculations.
5. Multiple-testing correction columns.
6. Normalized XY replicate import and multi-group handling.
7. All four bundled hookworm-style example workflows.

Run:

```bash
pytest -q
```

## Validation layer 2: Simulated benchmark datasets with known IC50 values

The script `scripts/generate_simulated_benchmarks.py` creates synthetic datasets for all four assay workflows with known IC50 values. Count endpoints are generated using binomial sampling. Motility endpoints are generated as continuous activity measurements that decline with the known inhibition curve and are then normalized through the same raw-motility function used by the app.

Run:

```bash
python scripts/generate_simulated_benchmarks.py
```

Outputs are written to `benchmarks/`, including:

- `simulated_egg_hatch.csv`
- `simulated_larval_development.csv`
- `simulated_motility.csv`
- `simulated_survival.csv`
- assay-specific analyzed-row, dose-summary, fit-summary, and fold-resistance files
- `benchmark_fit_accuracy_summary.csv`
- `benchmark_fold_resistance_summary.csv`
- `benchmark_warnings.csv`

## Validation layer 3: Bundled example datasets

The script `scripts/run_all_examples.py` runs each bundled hookworm-style example through the full analysis workflow and exports row-level data, dose summaries, IC50 fits, fold-resistance tables, and assay-appropriate pairwise tests.

Run:

```bash
python scripts/run_all_examples.py
```

Outputs are written to `validation_outputs/`.

## Validation layer 4: Real case studies

External validation should include susceptible and resistant parasite groups tested in independently generated experiments. Recommended case studies include:

1. Thiabendazole egg hatch assay.
2. Ivermectin larval development assay.
3. Ivermectin motility assay with archived raw activity measurements and zero-dose controls.
4. Ivermectin survival/mortality assay.

A motility validation dataset should document the parasite stage, exposure time, temperature, medium, stimulus conditions, tracking platform, activity metric, replicate structure, and whether values were normalized before import. A parameter-matched comparison with GraphPad Prism or another nonlinear-regression package should use identical replicate values, response direction, weighting, asymptote constraints, and confidence-interval settings.
