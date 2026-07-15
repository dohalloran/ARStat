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

Each workflow has its own response calculation and expected raw data columns.

## 3. Load data

Use one of the bundled hookworm-style examples or upload a CSV/XLSX file. Choose either **Raw assay measurements** or **Normalized XY replicate table**. The normalized layout uses one dose column and adjacent individual replicate columns; ARStat calculates mean, standard deviation, and n internally. Blank templates are available in the app and in the `templates/` folder.


### Normalized XY replicate layout

A normalized import should contain one dose/X column and one or more individual replicate/Y columns:

```text
dose,replicate_1,replicate_2,replicate_3
0,100,100,100
1,91,93,90
3,74,78,76
10,40,44,42
```

Do not import precomputed means, medians, standard deviations, or sample sizes. Select the replicate columns in the sidebar and indicate whether responses are percentages (`0–100`) or fractions (`0–1`).

For multiple experimental groups, use one row per group and dose:

```text
Group,Drug,Dose,Rep1,Rep2,Rep3
WMD,TBZ,0,100,99,101
WMD,TBZ,0.5,91,89,90
KGR,TBZ,0,100,100,99
KGR,TBZ,0.5,98,96,97
```

Map **Experimental group column** to any column representing strain, isolate, genetic background, population, or treatment. The drug/compound column is also optional. ARStat fits one curve per drug-by-group combination and enables fold-resistance calculations when a reference group is selected.

## 4. Confirm column settings

Check that ARStat has correctly identified:

- dose column
- experimental group column, when present
- drug/compound column, when present
- individual replicate/Y columns
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
