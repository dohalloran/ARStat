# ARStat v1.3.1

[![Tests](https://github.com/dohalloran/ARStat/actions/workflows/tests.yml/badge.svg)](https://github.com/dohalloran/ARStat/actions/workflows/tests.yml)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/)

**Live app:** [https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/](https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/)

ARStat is a Streamlit web application and scriptable Python backend for reproducible statistical analysis and visualization of anthelmintic dose-response assays.

ARStat supports three assay workflows:

- **Egg hatch assay**
- **Larval development assay**
- **Motility assay**

The app starts from raw assay measurements or normalized XY replicate tables, calculates assay-specific response variables, fits four-parameter dose-response curves, estimates IC50 values, calculates fold-resistance versus a selected reference group, performs exploratory dose-level tests with multiple-testing correction, and exports analysis tables, plots, Excel workbooks, and methods text.

## Bundled sample data

The `sample_data` folder includes one **simulated/illustrative** example dataset for each supported workflow, using *Ancylostoma caninum* labels. These bundled examples demonstrate the input formats and analysis routes; they are not the empirical validation datasets used in the manuscript:

| Assay | Isolates | Drug | Raw endpoint | Provenance |
|---|---|---|---|---|
| Egg hatch | WMD vs KGR | Thiabendazole | Eggs and L1 larvae | Simulated/illustrative count data |
| Larval development | WMD vs KGR | Ivermectin | Developed and undeveloped larvae | Simulated/illustrative count data |
| Motility | WMD vs KGR | Ivermectin | Continuous activity units | Simulated from a known 4PL model (IC50 20 and 80 nM); not experimental data |

See `sample_data/DATA_DESCRIPTION.md` for the file-by-file description.

## Quick start

### Use the live app

No installation is required. Open the app directly in a browser:

[https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/](https://arstat-jm7varr6fck8uajj4lgs6t.streamlit.app/)

### Run locally

```bash
conda create -n arstat python=3.11 -y
conda activate arstat
cd ARStat
pip install -r requirements.txt
streamlit run app.py
```

Then open the local URL shown by Streamlit, usually `http://localhost:8501`.

## Supported input layouts

### Raw assay measurements

Use assay-specific count or measurement columns. ARStat calculates the biological response before fitting the dose-response model.

For motility assays, the selected column may contain continuous activity measurements from WormLab or another tracking system, mean velocity, distance moved, activity counts, or a consistently defined manual motility score. Raw activity values are normalized within each fitted drug-by-group combination to the mean zero-dose control. Each fitted group must therefore contain a positive zero-dose control mean. ARStat then models **motility inhibition** while displaying relative motility as a descending curve when the traditional plot mode is selected.

Motility values may also be supplied already normalized as percentages (`0–100`) or fractions (`0–1`). The user specifies whether those values represent retained motility or motility inhibition.

### Normalized XY replicate table

Use the wide XY layout commonly used by dose-response applications:

```text
dose,replicate_1,replicate_2,replicate_3
0,100,100,100
0.5,92,89,91
2.5,78,82,80
5,48,43,46
```

The dose/X column and individual replicate/Y columns are required. Optional experimental-group and drug columns allow several strains, isolates, genetic backgrounds, populations, treatments, time points, or other groups to be analyzed in one file:

```text
Group,Drug,Dose,Rep1,Rep2,Rep3
WMD,TBZ,0,100,99,101
WMD,TBZ,0.5,91,89,90
KGR,TBZ,0,100,100,99
KGR,TBZ,0.5,98,96,97
```

ARStat maps the selected experimental-group column to its backward-compatible internal `strain` field and fits one curve per drug-by-group combination. Values may be entered as percentages (`0–100`) or fractions (`0–1`). Values below 0% or above 100% are flagged but retained rather than clipped, because control normalization, background correction, or replicate variation can legitimately produce out-of-range observations. ARStat reshapes the table internally and calculates the mean, standard deviation, and sample size. Do not import precomputed means, medians, standard deviations, or sample sizes. CSV and XLSX files are supported.

## Reproducibility and validation

Run the unit tests from the repository root (`conftest.py` puts the repository on the import path, so no `PYTHONPATH` setting is needed):

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

Estimate the coverage of the bootstrap IC50 confidence intervals by simulation (about 10 minutes on 2 cores):

```bash
python scripts/simulate_bootstrap_coverage.py --n-sim 120
```

Benchmark outputs are written to `benchmarks/`. Example validation outputs are written to `validation_outputs/`.

## What ARStat outputs

- analyzed row-level data
- dose summaries
- fitted IC50 table
- bootstrap IC50 confidence intervals when enabled
- fold-resistance versus a selected reference group
- fold-resistance confidence intervals when available
- per-dose pairwise tests with raw, Benjamini-Hochberg, and Bonferroni p-values
- downloadable Excel workbook
- export-ready PNG plot
- methods text describing the selected workflow

## Important statistical notes

- **Zero-dose controls:** zero-dose observations keep dose = 0 in all exported tables. Because log10(0) is undefined, the 4PL model evaluates them at one-tenth of the smallest positive dose in the fitted group, and plots show them at a symbolic left-edge tick labelled `0`.
- **Model constraints:** the lower and upper asymptotes are constrained to 0–1 (on the inhibition scale), the Hill slope to 0.05–10, and log10(IC50) to within 3 decades of the tested positive doses. Check the IC50 table for estimates that sit on a bound.
- **Bootstrap confidence intervals:** percentile intervals from resampling replicate rows within each fitted curve. They ignore plate, day, and experiment structure. In simulations (`scripts/simulate_bootstrap_coverage.py`; 120 data sets per design), nominal 95% intervals contained the true IC50 in 91–92% of data sets with 3–6 wells per concentration, so treat them as approximate.
- **IC50 interpretation:** IC50 is the midpoint between the fitted lower and upper asymptotes of the four-parameter logistic model. If the fitted maximum response is below 100%, IC50 is not necessarily the dose giving 50% absolute response.
- **Motility interpretation:** the reported value is a motility-inhibition IC50 under the selected activity definition, exposure time, parasite stage, and normalization. It should not be interpreted automatically as a lethal concentration.
- **Pairwise tests:** pairwise dose-level p-values are exploratory. ARStat reports both raw and adjusted p-values.
- **Count assays:** Fisher exact tests pool replicate counts at each dose and do not model replicate-to-replicate overdispersion.
- **Continuous assays:** motility and normalized XY workflows use replicate-level Mann-Whitney U tests by default. With 3 replicates per group an exact two-sided test cannot reach P < 0.05; the `exact_min_p` column reports this limit for each comparison.

## Required columns

For raw assay measurements, the minimum required fields are a strain/isolate or other group column, a dose column, and the assay response measurements. Drug, dose-unit, and replicate/well identifiers are optional interface mappings or metadata.

For normalized XY input, only a dose column and individual replicate response columns are required. Experimental-group and drug columns are optional.

Assay-specific raw columns:

- Egg hatch: `eggs`, `L1`
- Larval development: `developed`, `undeveloped`; conventional `L3` and `L1` headings are suggested as developed and undeveloped, respectively, but must be confirmed against the protocol
- Motility: one continuous activity or motility column, such as `motility`

## Repository structure

```text
app.py                         Streamlit app
arstat_core.py                 Reusable analysis backend
sample_data/                   Bundled simulated/illustrative example datasets for all three workflows
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
