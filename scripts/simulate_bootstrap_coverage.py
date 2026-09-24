"""Monte Carlo check of ARStat bootstrap IC50 confidence-interval coverage.

Simulates count-based egg hatch data with a known IC50, fits each data set with
``fit_dose_response(..., n_boot=200)``, and reports how often the nominal 95%
percentile interval contains the true IC50. Two designs are simulated:

* a design like the BCR validation experiments (7 concentrations x 3 wells,
  about 55 eggs per well, beta-binomial well-to-well overdispersion); and
* the bundled benchmark design (8 concentrations x 6 wells, about 85 eggs per
  well, binomial sampling).

Run from the repository root (about 10 minutes on 2 cores with the defaults):
    python scripts/simulate_bootstrap_coverage.py --n-sim 120
"""
from __future__ import annotations

import argparse
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from arstat_core import calculate_count_response, fit_dose_response, four_parameter_logistic  # noqa: E402

DESIGNS = {
    "BCR-like: 7 doses x 3 wells, overdispersed": dict(
        doses=[0, 0.5, 2.5, 5, 12.5, 25, 50.0], reps=3, true=5.0, hill=3.0, bottom=0.08, top=1.0, n=50, phi=49
    ),
    "Benchmark-like: 8 doses x 6 wells, binomial": dict(
        doses=[0, 3.125, 6.25, 12.5, 25, 50, 100, 200.0], reps=6, true=20.0, hill=1.2, bottom=0.02, top=0.95, n=80, phi=None
    ),
}


def simulate(seed: int, design: dict) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    floor = min(d for d in design["doses"] if d > 0) / 10
    rows = []
    for dose in design["doses"]:
        p = float(four_parameter_logistic(
            np.array([dose if dose > 0 else floor]), design["bottom"], design["top"], np.log10(design["true"]), design["hill"]
        )[0])
        for _ in range(design["reps"]):
            n = int(rng.poisson(design["n"])) + 5
            if design["phi"]:
                p_well = rng.beta(max(p * design["phi"], 1e-3), max((1 - p) * design["phi"], 1e-3))
            else:
                p_well = p
            unhatched = int(rng.binomial(n, p_well))
            rows.append({"group": "sim", "dose": dose, "eggs": unhatched, "L1": n - unhatched})
    return pd.DataFrame(rows)


def one(job: tuple[int, str, int]) -> dict:
    seed, name, n_boot = job
    warnings.filterwarnings("ignore")
    design = DESIGNS[name]
    prepared, _ = calculate_count_response(simulate(seed, design), "L1", "eggs", "Egg hatch")
    fit, _ = fit_dose_response(prepared, ["group"], n_boot=n_boot, random_seed=seed)
    row = fit.iloc[0]
    return {"design": name, "seed": seed, "true_IC50": design["true"], "IC50": row.IC50,
            "CI_low": row.IC50_CI_low, "CI_high": row.IC50_CI_high}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n-sim", type=int, default=120)
    parser.add_argument("--n-boot", type=int, default=200)
    parser.add_argument("--processes", type=int, default=2)
    args = parser.parse_args()

    jobs = [(seed, name, args.n_boot) for name in DESIGNS for seed in range(args.n_sim)]
    with Pool(args.processes) as pool:
        results = pd.DataFrame(pool.map(one, jobs))

    rows = []
    for name, g in results.groupby("design", sort=False):
        ok = g.CI_low.notna() & g.CI_high.notna()
        covered = (g.CI_low <= g.true_IC50) & (g.true_IC50 <= g.CI_high)
        coverage = covered[ok].mean()
        rows.append({
            "design": name,
            "n_simulated": len(g),
            "n_with_CI": int(ok.sum()),
            "nominal_level": 0.95,
            "coverage": coverage,
            "coverage_mc_se": np.sqrt(coverage * (1 - coverage) / ok.sum()),
            "median_log10_CI_width": np.median(np.log10(g.CI_high[ok] / g.CI_low[ok])),
            "median_percent_error_IC50": 100 * np.median(g.IC50 / g.true_IC50 - 1),
        })
    summary = pd.DataFrame(rows)
    out = ROOT / "benchmarks" / "bootstrap_coverage_summary.csv"
    summary.to_csv(out, index=False)
    print(summary.to_string(index=False))
    print("Wrote", out)


if __name__ == "__main__":
    main()
