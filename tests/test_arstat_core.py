import numpy as np
import pandas as pd

from arstat_core import (
    INFO_PREFIX,
    assay_warnings,
    calculate_count_response,
    calculate_resistance_ratios,
    drop_entirely_empty_columns,
    fit_dose_response,
    four_parameter_logistic,
    pairwise_count_tests,
    pairwise_continuous_tests,
    prepare_motility_response,
    prepare_normalized_xy_response,
    suggest_count_columns,
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


def test_motility_workflow_normalizes_to_group_controls_and_fits_ic50():
    df = pd.read_csv("sample_data/motility_example.csv")
    prepared, warnings = prepare_motility_response(
        df,
        dose_col="dose",
        motility_col="motility",
        group_cols=["drug", "strain"],
        value_scale="raw",
        response_direction="raw_outcome",
    )
    assert any("normalized to the mean zero-dose control" in warning for warning in warnings)
    assert "motility_fraction" in prepared.columns
    assert "motility_inhibition_fraction" in prepared.columns
    zero = prepared.loc[prepared["dose"] == 0].groupby("strain")["motility_fraction"].mean()
    assert np.allclose(zero.to_numpy(), 1.0)
    high = prepared.loc[prepared["dose"] == prepared["dose"].max()].groupby("strain")["response_fraction"].mean()
    assert (high > 0.65).all()
    summary, _ = _assert_basic_fit(prepared, group_cols=["drug", "strain"])
    assert (summary["top"] > summary["bottom"]).all()



def test_motility_hypermotility_is_retained_not_clipped():
    df = pd.DataFrame(
        {
            "strain": ["WT"] * 8,
            "drug": ["Drug"] * 8,
            "dose": [0, 0, 1, 1, 10, 10, 100, 100],
            "motility": [90, 110, 120, 105, 55, 50, 5, 8],
        }
    )
    prepared, warnings = prepare_motility_response(
        df,
        dose_col="dose",
        motility_col="motility",
        group_cols=["drug", "strain"],
        value_scale="raw",
        response_direction="raw_outcome",
    )
    assert prepared["motility_fraction"].max() > 1.0
    assert prepared["response_fraction"].min() < 0.0
    assert any("were retained" in warning for warning in warnings)

def test_motility_raw_values_require_zero_dose_control_per_group():
    df = pd.DataFrame(
        {
            "strain": ["WT"] * 4,
            "drug": ["Drug"] * 4,
            "dose": [1, 3, 10, 30],
            "motility": [100, 80, 40, 10],
        }
    )
    try:
        prepare_motility_response(
            df,
            dose_col="dose",
            motility_col="motility",
            group_cols=["drug", "strain"],
            value_scale="raw",
        )
    except ValueError as exc:
        assert "zero-dose control" in str(exc)
    else:
        raise AssertionError("Expected raw motility preparation to require a zero-dose control.")


def test_count_pairwise_tests_include_multiple_testing_adjustment():
    df = pd.read_csv("sample_data/egg_hatch_example.csv")
    prepared, _ = calculate_count_response(
        df, success_col="L1", failure_col="eggs", assay_name="Egg hatch"
    )
    tests = pairwise_count_tests(prepared, comparison_col="strain", dose_col="dose", stratify_cols=["drug"])
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
    assert "hatch_fraction" in egg.columns
    assert "development_fraction" in larval.columns
    assert "survival_fraction" in survival.columns


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


def test_normalized_xy_replicate_import_percent_scale():
    wide = pd.DataFrame({
        "dose": [0, 1, 3, 10, 30],
        "replicate_1": [100, 92, 75, 40, 8],
        "replicate_2": [100, 90, 78, 43, 5],
        "replicate_3": [100, 91, 76, 41, 6],
    })
    prepared, warnings = prepare_normalized_xy_response(
        wide,
        dose_col="dose",
        replicate_cols=["replicate_1", "replicate_2", "replicate_3"],
        assay_name="Egg hatch",
        dataset_label="BCR",
        drug_label="TBZ",
        unit="uM",
        value_scale="percent",
        response_direction="raw_outcome",
    )
    assert len(prepared) == 15
    assert prepared["hatch_fraction"].between(0, 1).all()
    assert prepared["response_fraction"].between(0, 1).all()
    assert np.isclose(prepared.loc[prepared["dose"] == 0, "hatch_fraction"].mean(), 1.0)
    assert prepared.loc[prepared["dose"] == 30, "response_fraction"].mean() > 0.9
    summary = summarize_by_dose(prepared, group_cols=["drug", "strain"], dose_col="dose")
    assert (summary["n"] == 3).all()
    assert summary["sd_response"].notna().all()
    fit_summary, _ = fit_dose_response(
        prepared, group_cols=["drug", "strain"], dose_col="dose", total_col=None
    )
    assert fit_summary.loc[0, "converged"]
    assert fit_summary.loc[0, "IC50"] > 0


def test_normalized_xy_import_does_not_require_summary_columns():
    wide = pd.DataFrame({
        "concentration": [0, 2, 20, 50, 100],
        "well_A": [1.0, 0.95, 0.80, 0.20, 0.02],
        "well_B": [1.0, 0.94, 0.82, 0.18, 0.01],
    })
    prepared, _ = prepare_normalized_xy_response(
        wide,
        dose_col="concentration",
        replicate_cols=["well_A", "well_B"],
        assay_name="Larval development",
        value_scale="fraction",
        response_direction="raw_outcome",
    )
    assert set(prepared["replicate"]) == {"well_A", "well_B"}
    assert "std" not in wide.columns
    assert "n" not in wide.columns
    assert "development_fraction" in prepared.columns


def test_normalized_xy_multigroup_import_preserves_groups_and_drugs():
    wide = pd.DataFrame({
        "Group": ["WMD"] * 6 + ["KGR"] * 6,
        "Drug": ["TBZ"] * 12,
        "Dose": [0, 0.5, 2.5, 5, 12.5, 25] * 2,
        "Rep1": [100, 92, 78, 48, 8, 0, 100, 98, 91, 78, 45, 12],
        "Rep2": [99, 90, 80, 45, 5, 0, 100, 96, 89, 75, 42, 10],
        "Rep3": [100, 91, 79, 47, 6, 0, 99, 97, 90, 77, 44, 11],
    })
    prepared, warnings = prepare_normalized_xy_response(
        wide,
        dose_col="Dose",
        replicate_cols=["Rep1", "Rep2", "Rep3"],
        assay_name="Egg hatch",
        group_col="Group",
        drug_col="Drug",
        unit="uM",
        value_scale="percent",
        response_direction="raw_outcome",
    )
    assert len(prepared) == 36
    assert set(prepared["strain"]) == {"WMD", "KGR"}
    assert set(prepared["drug"]) == {"TBZ"}
    assert prepared.groupby(["drug", "strain", "dose"]).size().eq(3).all()
    summary, results = fit_dose_response(
        prepared,
        group_cols=["drug", "strain"],
        dose_col="dose",
        total_col=None,
    )
    assert len(summary) == 2
    assert summary["converged"].all()
    rr = calculate_resistance_ratios(
        summary,
        group_col="strain",
        reference_group="WMD",
        fit_results=results,
        group_cols=["drug", "strain"],
    )
    assert set(rr["strain"]) == {"WMD", "KGR"}
    assert rr["fold_resistance_vs_reference"].notna().all()


def test_normalized_xy_multigroup_allows_group_without_drug_column():
    wide = pd.DataFrame({
        "Background": ["WT", "WT", "Mutant", "Mutant"],
        "Dose": [0, 10, 0, 10],
        "Y1": [100, 20, 100, 60],
        "Y2": [99, 22, 98, 62],
    })
    prepared, _ = prepare_normalized_xy_response(
        wide,
        dose_col="Dose",
        replicate_cols=["Y1", "Y2"],
        assay_name="Larval development",
        group_col="Background",
        drug_label="IVM",
        value_scale="percent",
        response_direction="raw_outcome",
    )
    assert set(prepared["strain"]) == {"WT", "Mutant"}
    assert set(prepared["drug"]) == {"IVM"}
    assert prepared["development_fraction"].between(0, 1).all()


def test_normalized_xy_motility_and_continuous_pairwise_tests():
    wide = pd.DataFrame(
        {
            "Group": ["WMD"] * 5 + ["KGR"] * 5,
            "Drug": ["IVM"] * 10,
            "Dose": [0, 3, 10, 30, 100] * 2,
            "Rep1": [105, 92, 70, 25, 3, 102, 99, 92, 67, 18],
            "Rep2": [98, 90, 68, 22, 2, 99, 98, 90, 65, 15],
            "Rep3": [100, 91, 69, 24, 4, 101, 97, 91, 66, 16],
        }
    )
    prepared, warnings = prepare_normalized_xy_response(
        wide,
        dose_col="Dose",
        replicate_cols=["Rep1", "Rep2", "Rep3"],
        assay_name="Motility",
        group_col="Group",
        drug_col="Drug",
        unit="nM",
        value_scale="percent",
        response_direction="raw_outcome",
    )
    assert prepared["motility_fraction"].max() > 1.0
    assert prepared["response_fraction"].min() < 0.0
    assert any("were retained" in warning for warning in warnings)
    summary, _ = fit_dose_response(
        prepared,
        group_cols=["drug", "strain"],
        dose_col="dose",
        total_col=None,
    )
    assert summary["converged"].all()
    tests = pairwise_continuous_tests(
        prepared,
        comparison_col="strain",
        dose_col="dose",
        response_col="response_fraction",
        stratify_cols=["drug"],
    )
    assert not tests.empty
    assert {"p_value", "p_value_bh", "p_value_bonferroni"}.issubset(tests.columns)


def test_normalized_xy_out_of_range_values_are_retained_not_clipped():
    wide = pd.DataFrame(
        {
            "dose": [0, 1, 10, 100],
            "rep1": [105, 92, 35, -2],
            "rep2": [100, 88, 30, 1],
        }
    )
    prepared, warnings = prepare_normalized_xy_response(
        wide,
        dose_col="dose",
        replicate_cols=["rep1", "rep2"],
        assay_name="Larval development",
        value_scale="percent",
        response_direction="raw_outcome",
    )
    assert prepared["development_fraction"].max() == 1.05
    assert prepared["development_fraction"].min() == -0.02
    assert np.isclose(prepared["response_fraction"].min(), -0.05)
    assert np.isclose(prepared["response_fraction"].max(), 1.02)
    assert any("were retained" in warning for warning in warnings)


def test_auto_scale_detection_is_informational_and_robust_to_small_fraction_outlier():
    wide = pd.DataFrame(
        {
            "dose": [0, 1, 10, 100],
            "rep1": [1.0, 0.8, 0.3, 1.6],
        }
    )
    prepared, messages = prepare_normalized_xy_response(
        wide,
        dose_col="dose",
        replicate_cols=["rep1"],
        assay_name="Egg hatch",
        value_scale="auto",
        response_direction="raw_outcome",
    )
    assert np.isclose(prepared.loc[prepared["dose"] == 0, "hatch_fraction"].iloc[0], 1.0)
    assert any(message.startswith(INFO_PREFIX) for message in messages)
    assert any("conservatively assumed fractions" in message for message in messages)


def test_single_normalized_replicate_generates_variance_warning():
    wide = pd.DataFrame({"dose": [0, 1, 10, 100], "rep1": [100, 80, 30, 5]})
    _, warnings = prepare_normalized_xy_response(
        wide,
        dose_col="dose",
        replicate_cols=["rep1"],
        assay_name="Egg hatch",
        value_scale="percent",
        response_direction="raw_outcome",
    )
    assert any("Only one replicate response column" in warning for warning in warnings)


def test_lda_l1_l3_aliases_are_suggested_conservatively():
    developed, undeveloped = suggest_count_columns(
        ["strain", "dose", "L1", "L3"],
        "Larval development",
    )
    assert developed == "L3"
    assert undeveloped == "L1"


def test_same_count_column_cannot_be_used_twice():
    df = pd.DataFrame({"L1": [10, 5]})
    try:
        calculate_count_response(df, "L1", "L1", "Egg hatch")
    except ValueError as exc:
        assert "must use different columns" in str(exc)
    else:
        raise AssertionError("Expected duplicate response-column mapping to fail.")


def test_entirely_empty_columns_are_removed():
    df = pd.DataFrame({"dose": [0, 1], "Rep1": [100, 50], "Unnamed: 4": [np.nan, np.nan]})
    cleaned, notes = drop_entirely_empty_columns(df)
    assert list(cleaned.columns) == ["dose", "Rep1"]
    assert any("Unnamed: 4" in note for note in notes)


def test_insufficient_dose_warning_says_no_ic50_will_be_calculated():
    df = pd.DataFrame(
        {
            "drug": ["IVM"] * 6,
            "strain": ["BCR"] * 6,
            "dose": [0, 0, 10, 10, 100, 100],
        }
    )
    warnings = assay_warnings(df, ["drug", "strain"], dose_col="dose")
    assert any("no IC50 will be calculated" in warning for warning in warnings)
    assert any("three positive concentrations plus a zero-dose control" in warning for warning in warnings)
