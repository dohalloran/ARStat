# Sample data description

The `sample_data` folder includes one **simulated/illustrative** example dataset for each supported ARStat workflow. These files demonstrate input structure and app behavior; they are **not** the biological validation datasets used in the ARStat manuscript. Empirical validation uses archived *Ancylostoma caninum* egg-hatch/larval-development experiments and the published IDEA-ms *Brugia* motility data supplied in Supplementary Data S1.

| File | Assay | Species label | Isolates | Drug | Data source / purpose |
|---|---|---|---|---|---|
| `egg_hatch_example.csv` | Egg hatch | *Ancylostoma caninum* (label only) | WMD, KGR | Thiabendazole | **Simulated/illustrative** egg and L1 count data used to demonstrate hatch-inhibition analysis. |
| `larval_development_example.csv` | Larval development | *Ancylostoma caninum* (label only) | WMD, KGR | Ivermectin | **Simulated/illustrative** developed/undeveloped count data used to demonstrate development-inhibition analysis. |
| `motility_example.csv` | Motility | *Ancylostoma caninum* (label only) | WMD, KGR | Ivermectin | **Simulated** continuous activity values (zero-dose mean ≈ 100 units, 6 replicates per dose) generated from a four-parameter logistic model with IC50 = 20 nM (WMD) and 80 nM (KGR), Hill slope 1.2, and asymptotes 0.02 and 0.95. Demonstrates control normalization and motility-inhibition IC50 analysis. |

WMD is used as the reference label and KGR as the comparison label in the bundled examples. The labels are retained for workflow realism and should not be interpreted as provenance claims for the synthetic values.
