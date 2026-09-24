# Sample data description

The `sample_data` folder includes **real experimental data** for all three supported ARStat workflows: egg hatch, larval development, and motility.

| File | Assay | Species label | Isolates | Drug | Data source / purpose |
|---|---|---|---|---|---|
| `egg_hatch_example.csv` | Egg hatch | *Ancylostoma caninum* | WMD, KGR | Thiabendazole | Real experimental data; demonstrates hatch-inhibition analysis from egg and L1 counts. |
| `larval_development_example.csv` | Larval development | *Ancylostoma caninum* | WMD, KGR | Ivermectin | Real experimental data; demonstrates development-inhibition analysis from developed and undeveloped counts. |
| `motility_example.csv` | Motility | *Ancylostoma caninum* | WMD, KGR | Ivermectin | Real experimental data; demonstrates control normalization and motility-inhibition IC50 analysis from continuous activity measurements. |

WMD is used as the susceptible reference isolate in the bundled examples. KGR is used as the comparison/resistant isolate.
