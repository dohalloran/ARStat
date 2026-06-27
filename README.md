# ARStat v1.0.0

ARStat is a Streamlit web application and scriptable Python backend for reproducible statistical analysis and visualization of anthelmintic resistance dose-response assays.

ARStat supports four assay workflows:

- **Egg hatch assay**
- **Larval development assay**
- **Survival/mortality assay**
- **Motility assay**

The app starts from raw assay measurements, calculates assay-specific response variables, fits dose-response curves, estimates IC50 values, calculates fold-resistance versus a selected reference isolate, performs exploratory dose-level tests with multiple-testing correction, and exports analysis tables, plots, Excel workbooks, and methods text.

## Hookworm-style demonstration data

This release includes synthetic demonstration datasets using *Ancylostoma caninum* labels:

| Assay | Isolates | Drug |
|---|---|---|
| Egg hatch | WMD vs KGR | Thiabendazole |
| Larval development | WMD vs KGR | Ivermectin |
| Survival/mortality | WMD vs KGR | Ivermectin |
| Motility | WMD vs KGR | Ivermectin |

The bundled datasets are synthetic demonstration datasets for testing the interface and workflow. They are not primary experimental measurements and should not be cited as biological results.

## Quick start

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

All assays should include `strain`, `drug`, `dose`, `unit`, and `replicate`.

Assay-specific columns:

- Egg hatch: `eggs`, `L1`
- Larval development: `developed`, `undeveloped`
- Survival/mortality: `alive`, `dead`
- Motility: `motility_score`

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
