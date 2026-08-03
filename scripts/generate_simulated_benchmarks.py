"""Generate simulated ARStat benchmark datasets with known IC50 values.

This script creates synthetic datasets for egg hatch, larval development,
motility, and survival/mortality workflows, then fits them using ARStat's core
analysis functions. The outputs can be used in validation reports and figures.

Run from the repository root:
    python scripts/generate_simulated_benchmarks.py
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from arstat_core import (  # noqa: E402
    calculate_count_response,
    calculate_resistance_ratios,
    fit_dose_response,
    four_parameter_logistic,
    prepare_motility_response,
    summarize_by_dose,
)

OUTDIR = ROOT / "benchmarks"
OUTDIR.mkdir(exist_ok=True)

RNG = np.random.default_rng(20260625)
DOSES = np.array([0, 3.125, 6.25, 12.5, 25, 50, 100, 200], dtype=float)
REPLICATES = 6
GROUPS = [
    {"strain": "WMD", "expected_ic50": 20.0, "bottom": 0.02, "top": 0.95, "hill": 1.2},
    {"strain": "KGR", "expected_ic50": 80.0, "bottom": 0.02, "top": 0.95, "hill": 1.2},
]

ASSAYS = {
    "Egg hatch": {
        "kind": "count",
        "file": "simulated_egg_hatch.csv",
        "assay_label": "egg_hatch",
        "drug": "thiabendazole",
        "unit": "uM",
        "success_col": "L1",
        "failure_col": "eggs",
        "total": 80,
    },
    "Larval development": {
        "kind": "count",
        "file": "simulated_larval_development.csv",
        "assay_label": "larval_development",
        "drug": "ivermectin",
        "unit": "nM",
        "success_col": "developed",
        "failure_col": "undeveloped",
        "total": 75,
    },
    "Survival": {
        "kind": "count",
        "file": "simulated_survival.csv",
        "assay_label": "survival",
        "drug": "ivermectin",
        "unit": "nM",
        "success_col": "dead",
        "failure_col": "alive",
        "total": 30,
    },
    "Motility": {
        "kind": "continuous",
        "file": "simulated_motility.csv",
        "assay_label": "motility",
        "drug": "ivermectin",
        "unit": "nM",
        "motility_col": "motility",
        "control_activity": 1000.0,
        "noise_sd": 35.0,
    },
}


def response_at_dose(dose, cfg):
    return float(
        four_parameter_logistic(
            np.array([dose]),
            cfg["bottom"],
            cfg["top"],
            np.log10(cfg["expected_ic50"]),
            cfg["hill"],
        )[0]
    )


def make_count_assay(assay_name, meta):
    rows = []
    for group in GROUPS:
        for dose in DOSES:
            p_effect = np.clip(response_at_dose(dose, group), 0.001, 0.999)
            for rep in range(1, REPLICATES + 1):
                n = int(max(5, RNG.poisson(meta["total"])))
                affected = int(RNG.binomial(n, p_effect))
                unaffected = n - affected

                if assay_name in ["Egg hatch", "Larval development"]:
                    success = unaffected
                    failure = affected
                else:  # Survival: success is dead/affected
                    success = affected
                    failure = unaffected

                rows.append(
                    {
                        "experiment_id": f"SIM_{meta['assay_label'].upper()}_001",
                        "assay": meta["assay_label"],
                        "species": "Ancylostoma_caninum",
                        "strain": group["strain"],
                        "expected_ic50": group["expected_ic50"],
                        "drug": meta["drug"],
                        "dose": dose,
                        "unit": meta["unit"],
                        "replicate": rep,
                        "well": f"R{rep}_D{dose:g}",
                        meta["success_col"]: success,
                        meta["failure_col"]: failure,
                    }
                )
    return pd.DataFrame(rows)


def make_motility_assay(meta):
    rows = []
    for group in GROUPS:
        for dose in DOSES:
            inhibition = np.clip(response_at_dose(dose, group), 0.0, 1.0)
            expected_motility = meta["control_activity"] * (1.0 - inhibition)
            for rep in range(1, REPLICATES + 1):
                motility = max(0.0, RNG.normal(expected_motility, meta["noise_sd"]))
                rows.append(
                    {
                        "experiment_id": f"SIM_{meta['assay_label'].upper()}_001",
                        "assay": meta["assay_label"],
                        "species": "Ancylostoma_caninum",
                        "strain": group["strain"],
                        "expected_ic50": group["expected_ic50"],
                        "drug": meta["drug"],
                        "dose": dose,
                        "unit": meta["unit"],
                        "replicate": rep,
                        "well": f"R{rep}_D{dose:g}",
                        meta["motility_col"]: round(motility, 3),
                    }
                )
    return pd.DataFrame(rows)


def analyze_dataset(assay_name, df, meta):
    group_cols = ["strain", "drug"]
    if meta["kind"] == "continuous":
        analyzed, warnings = prepare_motility_response(
            df,
            dose_col="dose",
            motility_col=meta["motility_col"],
            group_cols=group_cols,
            value_scale="raw",
            response_direction="raw_outcome",
        )
        total_col = None
    else:
        analyzed, warnings = calculate_count_response(
            df,
            meta["success_col"],
            meta["failure_col"],
            assay_name,
        )
        total_col = "total_count"

    fit_summary, fit_results = fit_dose_response(
        analyzed,
        group_cols=group_cols,
        dose_col="dose",
        response_col="response_fraction",
        total_col=total_col,
        n_boot=100,
        random_seed=20260625,
    )
    fit_summary["assay_workflow"] = assay_name
    expected = df[["strain", "drug", "expected_ic50"]].drop_duplicates()
    fit_summary = fit_summary.merge(expected, on=["strain", "drug"], how="left")
    fit_summary["absolute_error"] = fit_summary["IC50"] - fit_summary["expected_ic50"]
    fit_summary["percent_error"] = 100 * fit_summary["absolute_error"] / fit_summary["expected_ic50"]

    rr = calculate_resistance_ratios(
        fit_summary,
        group_col="strain",
        reference_group="WMD",
        fit_results=fit_results,
        group_cols=group_cols,
    )
    dose_summary = summarize_by_dose(analyzed, group_cols=group_cols, dose_col="dose")

    return analyzed, fit_summary, rr, dose_summary, warnings


def main():
    all_fits = []
    all_ratios = []
    all_warnings = []

    for assay_name, meta in ASSAYS.items():
        if meta["kind"] == "continuous":
            df = make_motility_assay(meta)
        else:
            df = make_count_assay(assay_name, meta)

        data_path = OUTDIR / meta["file"]
        df.to_csv(data_path, index=False)

        analyzed, fit_summary, rr, dose_summary, warnings = analyze_dataset(assay_name, df, meta)
        slug = meta["assay_label"]
        analyzed.to_csv(OUTDIR / f"{slug}_analyzed_rows.csv", index=False)
        fit_summary.to_csv(OUTDIR / f"{slug}_fit_summary.csv", index=False)
        rr.to_csv(OUTDIR / f"{slug}_fold_resistance.csv", index=False)
        dose_summary.to_csv(OUTDIR / f"{slug}_dose_summary.csv", index=False)
        all_fits.append(fit_summary)
        all_ratios.append(rr.assign(assay_workflow=assay_name))
        for warning in warnings:
            all_warnings.append({"assay_workflow": assay_name, "warning": warning})

    combined_fits = pd.concat(all_fits, ignore_index=True)
    combined_ratios = pd.concat(all_ratios, ignore_index=True)
    combined_fits.to_csv(OUTDIR / "benchmark_fit_accuracy_summary.csv", index=False)
    combined_ratios.to_csv(OUTDIR / "benchmark_fold_resistance_summary.csv", index=False)
    pd.DataFrame(all_warnings, columns=["assay_workflow", "warning"]).to_csv(
        OUTDIR / "benchmark_warnings.csv", index=False
    )

    print("Wrote simulated benchmark datasets and summary tables to:", OUTDIR)
    print(
        combined_fits[
            ["assay_workflow", "strain", "drug", "expected_ic50", "IC50", "percent_error"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
