"""Run all bundled ARStat example datasets and export validation summaries.

Run from repository root:
    python scripts/run_all_examples.py
"""
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from arstat_core import (  # noqa: E402
    calculate_count_response,
    calculate_resistance_ratios,
    fit_dose_response,
    pairwise_count_tests,
    summarize_by_dose,
)

OUTDIR = ROOT / "validation_outputs"
OUTDIR.mkdir(exist_ok=True)

EXAMPLES = {
    "Egg hatch": {
        "path": ROOT / "sample_data" / "egg_hatch_example.csv",
        "success_col": "L1",
        "failure_col": "eggs",
        "reference": "WMD",
    },
    "Larval development": {
        "path": ROOT / "sample_data" / "larval_development_example.csv",
        "success_col": "developed",
        "failure_col": "undeveloped",
        "reference": "WMD",
    },
    "Survival": {
        "path": ROOT / "sample_data" / "survival_example.csv",
        "success_col": "dead",
        "failure_col": "alive",
        "reference": "WMD",
    },
}

def main():
    rows = []
    for assay_name, cfg in EXAMPLES.items():
        df = pd.read_csv(cfg["path"])
        group_cols = ["strain", "drug"]
        analyzed, warnings = calculate_count_response(df, cfg["success_col"], cfg["failure_col"], assay_name)
        total_col = "total_count"
        tests = pairwise_count_tests(analyzed, comparison_col="strain", dose_col="dose", stratify_cols=["drug"])
        fit_summary, fit_results = fit_dose_response(analyzed, group_cols=group_cols, total_col=total_col, n_boot=100)
        rr = calculate_resistance_ratios(fit_summary, group_col="strain", reference_group=cfg["reference"], fit_results=fit_results, group_cols=group_cols)
        dose_summary = summarize_by_dose(analyzed, group_cols=group_cols)
        slug = assay_name.lower().replace(" ", "_")
        analyzed.to_csv(OUTDIR / f"{slug}_analyzed_rows.csv", index=False)
        dose_summary.to_csv(OUTDIR / f"{slug}_dose_summary.csv", index=False)
        fit_summary.to_csv(OUTDIR / f"{slug}_fit_summary.csv", index=False)
        rr.to_csv(OUTDIR / f"{slug}_fold_resistance.csv", index=False)
        tests.to_csv(OUTDIR / f"{slug}_pairwise_tests.csv", index=False)
        for _, r in fit_summary.iterrows():
            rows.append({"assay_workflow": assay_name, "strain": r.get("strain"), "drug": r.get("drug"), "IC50": r.get("IC50"), "converged": r.get("converged"), "message": r.get("message"), "warnings": " | ".join(warnings)})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUTDIR / "all_examples_fit_summary.csv", index=False)
    print("Wrote validation outputs to:", OUTDIR)
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
