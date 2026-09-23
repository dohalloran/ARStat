# Sample data description

The `sample_data` folder includes **real experimental data** for the egg hatch, larval development, and motility workflows. The survival/mortality dataset is an illustrative example included so that all four ARStat assay workflows can be exercised from the interface.

| File | Assay | Species label | Isolates | Drug | Data source / purpose |
|---|---|---|---|---|---|
| `egg_hatch_example.csv` | Egg hatch | *Ancylostoma caninum* | WMD, KGR | Thiabendazole | Real experimental data; demonstrates hatch-inhibition analysis from egg and L1 counts. |
| `larval_development_example.csv` | Larval development | *Ancylostoma caninum* | WMD, KGR | Ivermectin | Real experimental data; demonstrates development-inhibition analysis from developed and undeveloped counts. |
| `motility_example.csv` | Motility | *Ancylostoma caninum* | WMD, KGR | Ivermectin | Real experimental data; demonstrates control normalization and motility-inhibition IC50 analysis from continuous activity measurements. |
| `survival_example.csv` | Survival/mortality | *Ancylostoma caninum* | WMD, KGR | Ivermectin | Illustrative example; demonstrates mortality analysis from alive/dead counts. |

WMD is used as the susceptible reference isolate in the bundled examples. KGR is used as the comparison/resistant isolate.
