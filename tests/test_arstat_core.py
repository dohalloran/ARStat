import numpy as np
import pandas as pd

from arstat_core import (
    calculate_count_response,
    calculate_motility_response,
    calculate_resistance_ratios,
    fit_dose_response,
    four_parameter_logistic,
    pairwise_count_tests,
    pairwise_continuous_tests,
    summarize_by_dose,
)


def _assert_basic_fit(prepared, group_cols, dose_col="dose"):
    summary, results = fit_dose_response(prepared, group_cols=group_cols, dose_col=dose_col)
    assert len(summary) == 2
    assert "IC50" in summary.columns
    assert summary["IC50"].notna().all()
    assert (summary["IC50"] > 0).all()
    dose_summary = summarize_by_dose(prepared, group_cols=group_cols, dose_col=dose_col)
    assert not dose_summary.empty
    rr = calculate_resistance_ratios(summary, group_col="strain", reference_group="WMD", fit_results=results, group_cols=group_cols)
    assert "fold_resistance_vs_reference" in rr.columns
    assert "fold_resistance_CI_low" in rr.columns
    assert "fold_resistance_CI_high" in rr.columns
    assert rr["fold_resistance_vs_reference"].notna().all()
    return summary, rr


def test_egg_hatch_workflow():
    df = pd.read_csv("sample_data/egg_hatch_example.csv")
    prepared, warnings = calculate_count_response(
        df, success_col="L1", failure_col="eggs", assay_name="Egg hatch"
    )
    assert warnings == []
    assert prepared["response_fraction"].between(0, 1).all()
    _assert_basic_fit(prepared, group_cols=["drug", "strain"])


def test_larval_development_workflow():
    df = pd.read_csv("sample_data/larval_development_example.csv")
    prepared, warnings = calculate_count_response(
        df,
        success_col="developed",
        failure_col="undeveloped",
        assay_name="Larval development",
    )
    assert warnings == []
    assert prepared["response_fraction"].between(0, 1).all()
    zero_mean = prepared.loc[prepared["dose"] == 0, "response_fraction"].mean()
    high_mean = prepared.loc[prepared["dose"] == prepared["dose"].max(), "response_fraction"].mean()
    assert high_mean > zero_mean
    summary, _ = _assert_basic_fit(prepared, group_cols=["drug", "strain"])
    assert (summary["top"] > summary["bottom"]).all()


def test_survival_workflow():
    df = pd.read_csv("sample_data/survival_example.csv")
    prepared, warnings = calculate_count_response(
        df, success_col="dead", failure_col="alive", assay_name="Survival"
    )
    assert warnings == []
    assert prepared["response_fraction"].between(0, 1).all()
    _assert_basic_fit(prepared, group_cols=["drug", "strain"])


def test_motility_workflow():
    df = pd.read_csv("sample_data/motility_example.csv")
    prepared, warnings = calculate_motility_response(
        df,
        score_col="motility_score",
        group_cols=["drug", "strain"],
        dose_col="dose",
    )
    assert all("missing" not in warning.lower() for warning in warnings)
    assert prepared["response_fraction"].between(0, 1).all()
    _assert_basic_fit(prepared, group_cols=["drug", "strain"])


def test_count_pairwise_tests_include_multiple_testing_adjustment():
    df = pd.read_csv("sample_data/egg_hatch_example.csv")
    prepared, _ = calculate_count_response(
        df, success_col="L1", failure_col="eggs", assay_name="Egg hatch"
    )
    tests = pairwise_count_tests(prepared, comparison_col="strain", dose_col="dose", stratify_cols=["drug"])
    assert not tests.empty
    assert "p_value_bh" in tests.columns
    assert "p_value_bonferroni" in tests.columns


def test_motility_pairwise_tests_include_multiple_testing_adjustment():
    df = pd.read_csv("sample_data/motility_example.csv")
    prepared, _ = calculate_motility_response(
        df, score_col="motility_score", group_cols=["drug", "strain"], dose_col="dose"
    )
    tests = pairwise_continuous_tests(prepared, comparison_col="strain", dose_col="dose", stratify_cols=["drug"])
    assert not tests.empty
    assert "p_value_bh" in tests.columns
    assert "p_value_bonferroni" in tests.columns


def test_bootstrap_resistance_ratio_confidence_interval_columns():
    df = pd.read_csv("sample_data/survival_example.csv")
    prepared, _ = calculate_count_response(
        df, success_col="dead", failure_col="alive", assay_name="Survival"
    )
    summary, results = fit_dose_response(
        prepared, group_cols=["drug", "strain"], dose_col="dose", n_boot=30
    )
    rr = calculate_resistance_ratios(
        summary,
        group_col="strain",
        reference_group="WMD",
        fit_results=results,
        group_cols=["drug", "strain"],
    )
    assert "fold_resistance_CI_method" in rr.columns
    assert "fold_resistance_CI_low" in rr.columns
    assert "fold_resistance_CI_high" in rr.columns


def test_motility_above_control_emits_warning_and_clips_to_zero():
    df = pd.DataFrame(
        {
            "drug": ["IVM", "IVM", "IVM"],
            "strain": ["WMD", "WMD", "WMD"],
            "dose": [0, 0, 1],
            "motility_score": [10, 10, 12],
        }
    )
    prepared, warnings = calculate_motility_response(
        df, score_col="motility_score", group_cols=["drug", "strain"], dose_col="dose"
    )
    assert any("above the matched zero-dose control" in warning for warning in warnings)
    stimulated = prepared.loc[prepared["dose"] == 1].iloc[0]
    assert np.isclose(stimulated["motility_inhibition_fraction_unclipped"], -0.2)
    assert np.isclose(stimulated["response_fraction"], 0.0)


def test_decreasing_fit_is_flagged_in_message():
    df = pd.DataFrame(
        {
            "drug": ["Drug"] * 12,
            "strain": ["Test"] * 12,
            "dose": [0, 0, 1, 1, 3, 3, 10, 10, 30, 30, 100, 100],
            "response_fraction": [0.95, 0.92, 0.85, 0.82, 0.65, 0.62, 0.42, 0.40, 0.20, 0.18, 0.05, 0.04],
            "total_count": [100] * 12,
        }
    )
    summary, _ = fit_dose_response(df, group_cols=["drug", "strain"], dose_col="dose")
    assert summary.loc[0, "converged"]
    assert summary.loc[0, "top"] < summary.loc[0, "bottom"]
    assert "fitted top is below fitted bottom" in summary.loc[0, "message"]


def test_raw_outcome_columns_exist_for_traditional_plot_mode():
    egg, _ = calculate_count_response(
        pd.read_csv("sample_data/egg_hatch_example.csv"),
        success_col="L1",
        failure_col="eggs",
        assay_name="Egg hatch",
    )
    larval, _ = calculate_count_response(
        pd.read_csv("sample_data/larval_development_example.csv"),
        success_col="developed",
        failure_col="undeveloped",
        assay_name="Larval development",
    )
    survival, _ = calculate_count_response(
        pd.read_csv("sample_data/survival_example.csv"),
        success_col="dead",
        failure_col="alive",
        assay_name="Survival",
    )
    motility, _ = calculate_motility_response(
        pd.read_csv("sample_data/motility_example.csv"),
        score_col="motility_score",
        group_cols=["drug", "strain"],
        dose_col="dose",
    )
    assert "hatch_fraction" in egg.columns
    assert "development_fraction" in larval.columns
    assert "survival_fraction" in survival.columns
    assert "normalized_motility" in motility.columns


def test_traditional_curve_is_complement_of_inhibition_fit():
    df = pd.read_csv("sample_data/egg_hatch_example.csv")
    prepared, _ = calculate_count_response(
        df, success_col="L1", failure_col="eggs", assay_name="Egg hatch"
    )
    summary, _ = fit_dose_response(prepared, group_cols=["drug", "strain"], dose_col="dose")
    row = summary.loc[summary["strain"] == "WMD"].iloc[0]
    positive = prepared[(prepared["strain"] == "WMD") & (prepared["dose"] > 0)].copy()
    predicted_inhibition = four_parameter_logistic(
        positive["dose"].to_numpy(dtype=float),
        row["bottom"],
        row["top"],
        np.log10(row["IC50"]),
        row["hill_slope"],
    )
    predicted_hatch = 1 - predicted_inhibition
    residual = np.abs(predicted_hatch - positive["hatch_fraction"].to_numpy(dtype=float))
    assert residual.mean() < 0.20
