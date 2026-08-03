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
- Motility
- Survival/mortality

Each workflow has its own response calculation and expected raw data columns.

## 3. Load data

Use one of the bundled hookworm-style examples or upload a CSV/XLSX file. Choose either **Raw assay measurements** or **Normalized XY replicate table**. Blank templates are available in the app and in the `templates/` folder.

### Raw motility layout

A motility file should contain one row per replicate or well and one continuous activity column:

```text
strain,drug,dose,unit,replicate,motility
WMD,ivermectin,0,nM,1,1012
WMD,ivermectin,0,nM,2,988
WMD,ivermectin,10,nM,1,731
WMD,ivermectin,10,nM,2,705
```

The measurement may be tracking-derived activity units, velocity, distance moved, percentage motile, or a consistently defined manual motility score.

Choose the appropriate value scale:

- **Raw activity units:** ARStat divides every value by the mean zero-dose value within the same fitted group. A positive zero-dose control mean is required for every drug-by-group curve.
- **0–100 percent:** values are already normalized percentages.
- **0–1 fraction:** values are already normalized fractions.

For normalized values, specify whether higher values mean retained motility/activity or motility inhibition. ARStat always fits the increasing affected response, so retained motility is converted to `1 - motility` before fitting.

### Normalized XY replicate layout

A normalized import should contain one dose/X column and one or more individual replicate/Y columns:

```text
dose,replicate_1,replicate_2,replicate_3
0,100,100,100
1,91,93,90
3,74,78,76
10,40,44,42
```

Do not import precomputed means, medians, standard deviations, or sample sizes. Select the replicate columns and indicate whether responses are percentages (`0–100`) or fractions (`0–1`). Auto-detection is reported as an informational note. Normalized values below 0% or above 100% are retained, not clipped, and are flagged for review.

For multiple experimental groups, use one row per group and dose:

```text
Group,Drug,Dose,Rep1,Rep2,Rep3
WMD,TBZ,0,100,99,101
WMD,TBZ,0.5,91,89,90
KGR,TBZ,0,100,100,99
KGR,TBZ,0.5,98,96,97
```

Map **Experimental group column** to any column representing strain, isolate, genetic background, population, treatment, or time point. The drug/compound column is also optional. ARStat fits one curve per drug-by-group combination and enables fold-resistance calculations when a reference group is selected.

## 4. Confirm column settings

Check that ARStat has correctly identified:

- dose column
- experimental-group column, when present
- drug/compound column, when present
- individual replicate/Y columns for normalized XY input
- assay-specific count or motility columns
- for conventional LDA tables, ARStat may suggest `L3` as developed and `L1` as undeveloped; confirm that this matches the species, stage definitions, and protocol
- value scale and response direction for motility
- reference group, usually the susceptible isolate

For bundled examples, WMD is used as the susceptible reference isolate and KGR as the resistant demonstration isolate.

## 5. Run ARStat

Click **Run ARStat**. The app calculates assay-specific responses, fits dose-response models, estimates IC50 values, calculates fold-resistance, runs exploratory pairwise tests, and creates plots.

## 6. Interpret outputs

### IC50

The IC50 is the dose corresponding to the midpoint between the fitted lower and upper asymptotes of the dose-response curve. It is not always the dose producing an absolute 50% response if the fitted top response is less than 100%.

For motility, IC50 means the concentration producing the model-defined midpoint of **motility inhibition** relative to the selected control and activity definition. It does not establish that the worms are dead.

### Fold-resistance versus reference

Fold-resistance is calculated as:

```text
IC50_test / IC50_reference
```

A value of 4 means that the test group requires approximately four times more drug than the reference group to reach the same model-defined midpoint.

### Pairwise tests

Pairwise tests are exploratory dose-level comparisons:

- Count assays pool success/failure counts and use Fisher exact tests.
- Motility and normalized replicate inputs use Mann-Whitney U tests on replicate responses.

Both Benjamini-Hochberg and Bonferroni adjusted p-values are reported.

## 7. Review motility quality-control messages

Motility fits depend strongly on the measurement definition and controls. Review warnings for:

- missing zero-dose controls
- a zero or negative control mean
- values above 100% after normalization, which may indicate biological hypermotility or measurement noise
- normalized values below 0% or above 100%; these values are preserved for fitting and should be checked against the selected scale and normalization method
- fewer than four unique dose levels, which prevents a four-parameter IC50 fit
- failure to bracket the fitted IC50
- incomplete or decreasing response curves

Record parasite species, life stage, exposure duration, stimulus conditions, temperature, medium, tracking system, and the activity metric used when reporting a motility IC50.

## 8. Download outputs

ARStat can download:

- analyzed data
- dose summary table
- IC50 table
- fold-resistance table
- pairwise test table
- Excel workbook
- PNG plot
- methods text

## 9. Run validation scripts

From the repository root:

```bash
python scripts/generate_simulated_benchmarks.py
python scripts/run_all_examples.py
```

These scripts generate benchmark and validation tables for all four supported workflows.
