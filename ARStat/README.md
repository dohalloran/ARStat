# ARStat v1.1.1

[![Tests](https://github.com/dohalloran/ARStat/actions/workflows/tests.yml/badge.svg)](https://github.com/dohalloran/ARStat/actions/workflows/tests.yml)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/)

**Live app:** [https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/](https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/)

ARStat is a Streamlit web application and scriptable Python backend for reproducible statistical analysis and visualization of anthelmintic resistance dose-response assays.

ARStat supports three assay workflows:

- **Egg hatch assay**
- **Larval development assay**
- **Survival/mortality assay**

The app starts from raw assay measurements or normalized XY replicate tables, calculates assay-specific response variables, fits dose-response curves, estimates IC50 values, calculates fold-resistance versus a selected reference isolate, performs exploratory dose-level tests with multiple-testing correction, and exports analysis tables, plots, Excel workbooks, and methods text.

## Hookworm-style demonstration data

This release includes synthetic demonstration datasets using *Ancylostoma caninum* labels:

| Assay | Isolates | Drug |
|---|---|---|
| Egg hatch | WMD vs KGR | Thiabendazole |
| Larval development | WMD vs KGR | Ivermectin |
| Survival/mortality | WMD vs KGR | Ivermectin |

The bundled datasets are synthetic demonstration datasets for testing the interface and workflow. They are not primary experimental measurements and should not be cited as biological results.

## Quick start

### Use the live app

No installation required. Open the app directly in your browser:

[https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/](https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/)

### Run locally

```bash
conda create -n arstat python=3.11 -y
conda activate arstat
cd ARStat
pip install -r requirements.txt
streamlit run app.py
```

Then open the local URL shown by Streamlit, usually:

```text
http://localhost:8501
```


## Supported input layouts

### Raw assay measurements

Use assay-specific count or score columns. ARStat calculates the biological response before fitting the dose-response model.

### Normalized XY replicate table

Use the same wide XY layout commonly used by dose-response applications:

```text
dose,replicate_1,replicate_2,replicate_3
0,100,100,100
0.5,92,89,91
2.5,78,82,80
5,48,43,46
```

The dose/X column and individual replicate/Y columns are required. Optional experimental-group and drug columns allow several strains, isolates, genetic backgrounds, populations, or treatment groups to be analyzed in one file:

```text
Group,Drug,Dose,Rep1,Rep2,Rep3
WMD,TBZ,0,100,99,101
WMD,TBZ,0.5,91,89,90
KGR,TBZ,0,100,100,99
KGR,TBZ,0.5,98,96,97
```

ARStat maps the selected experimental-group column to its internal backward-compatible `strain` field and fits one curve per drug-by-group combination. Values may be entered as percentages (`0–100`) or fractions (`0–1`). ARStat reshapes the table internally and calculates the mean, standard deviation, and sample size. Do not precompute or import means, medians, standard deviations, or sample sizes. CSV and XLSX files are supported.

## Reproducibility and validation

Run the unit tests:

```bash
pip install pytest
pytest -q
```

Generate simulated benchmark datasets with known IC50 values:

```bash
python scripts/generate_simulated_benchmarks.py
```

Run all bundled example datasets through the full workflow:

```bash
python scripts/run_all_examples.py
```

Benchmark outputs are written to `benchmarks/`. Example validation outputs are written to `validation_outputs/`.

## What ARStat outputs

- analyzed row-level data
- dose summaries
- fitted IC50 table
- IC50 confidence intervals when bootstrap fitting is enabled
- fold-resistance versus a selected reference isolate
- fold-resistance confidence intervals when available
- per-dose pairwise tests with raw, Benjamini-Hochberg, and Bonferroni p-values
- downloadable Excel workbook
- export-ready PNG plot
- methods text

## Important statistical notes

- **Zero-dose controls:** log-scaled plots cannot display x=0, so ARStat shows zero-dose controls at a symbolic left-edge tick labelled `0`. These points are still included in calculations and summaries.
- **IC50 interpretation:** IC50 is the midpoint between the fitted lower and upper asymptotes of the four-parameter logistic model. If the fitted maximum response is below 100%, IC50 is not necessarily the dose giving 50% absolute response.
- **Pairwise tests:** pairwise dose-level p-values are exploratory. ARStat reports both raw and adjusted p-values.
- **Count assays:** Fisher exact tests pool replicate counts at each dose. This is useful as a simple exploratory test, but it does not model replicate-to-replicate overdispersion.

## Required columns

For raw assay measurements, files should include `strain`, `drug`, `dose`, `unit`, and `replicate` or equivalent columns that can be mapped in the interface.

For normalized XY replicate input, only a dose column and individual replicate response columns are required. Experimental-group and drug columns are optional. An experimental group can represent a strain, isolate, genetic background, population, treatment, or another category for which an independent curve should be fitted.

Assay-specific raw-count columns:

- Egg hatch: `eggs`, `L1`
- Larval development: `developed`, `undeveloped`
- Survival/mortality: `alive`, `dead`

## Repository structure

```text
app.py                         Streamlit app
arstat_core.py                 Reusable analysis backend
sample_data/                   Synthetic hookworm-style examples
templates/                     Blank CSV templates
tests/                         Unit tests
scripts/                       Benchmark and validation scripts
benchmarks/                    Simulated benchmark outputs
docs/                          User, validation, and statistical documentation
```

## License

MIT License.

## Citation

Please cite the archived software release listed in `CITATION.cff`. If no DOI has been minted yet, cite the GitHub repository and release tag.
