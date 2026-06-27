# ARStat user guide

## 1. Start the app

```bash
conda activate arstat
cd ARStat
streamlit run app.py
```

Open the local URL shown by Streamlit.

## 2. Choose an assay workflow

ARStat supports:

- Egg hatch
- Larval development
- Survival/mortality
- Motility

Each workflow has its own response calculation and expected raw data columns.

## 3. Load data

Use one of the bundled hookworm-style examples or upload your own CSV file. Blank templates are available in the app and in the `templates/` folder.

## 4. Confirm column settings

Check that ARStat has correctly identified:

- dose column
- strain/isolate column
- drug column
- response columns, such as `eggs` and `L1` for egg hatch assays
- reference isolate, usually the susceptible isolate

For bundled examples, WMD is used as the susceptible reference isolate and KGR as the resistant demonstration isolate.

## 5. Run ARStat

Click **Run ARStat**. The app will calculate assay-specific responses, fit dose-response models, estimate IC50 values, calculate fold-resistance, run exploratory pairwise tests, and create plots.

## 6. Interpret outputs

### IC50

The IC50 is the dose corresponding to the midpoint between the fitted lower and upper asymptotes of the dose-response curve. It is not always the dose producing an absolute 50% response if the fitted top response is less than 100%.

### Fold-resistance versus reference

Fold-resistance is calculated as:

```text
IC50_test / IC50_reference
```

A value of 4 means that the test isolate requires about four times more drug than the reference isolate to reach the same model-defined 50% response.

### Pairwise tests

Pairwise tests are exploratory dose-level comparisons. For count assays, ARStat pools replicate counts at each dose and uses Fisher exact tests. This does not model replicate-to-replicate overdispersion, so adjusted p-values should be interpreted cautiously.

## 7. Download outputs

ARStat can download:

- analyzed data
- dose summary table
- IC50 table
- fold-resistance table
- pairwise test table
- Excel workbook
- PNG plot
- methods text

## 8. Run validation scripts

From the repository root:

```bash
python scripts/generate_simulated_benchmarks.py
python scripts/run_all_examples.py
```

These scripts generate benchmark and validation tables for reproducibility checks.
