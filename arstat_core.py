"""Core analysis functions for ARStat.

ARStat is designed for antinematodal / anthelmintic dose-response assays.
The Streamlit app imports these functions, but they can also be used directly
from Python scripts or notebooks.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import fisher_exact, mannwhitneyu, ttest_ind


ASSAY_PRESETS = {
    "Egg hatch": {
        "success_default": "L1",
        "failure_default": "eggs",
        "success_label": "Hatched larvae / L1",
        "failure_label": "Unhatched eggs",
        "raw_fraction": "hatch_fraction",
        "effect_fraction": "inhibition_fraction",
        "effect_label": "Hatch inhibition",
    },
    "Larval development": {
        "success_default": "developed",
        "failure_default": "undeveloped",
        "success_label": "Developed larvae",
        "failure_label": "Undeveloped larvae",
        "raw_fraction": "development_fraction",
        "effect_fraction": "inhibition_fraction",
        "effect_label": "Development inhibition",
    },
    "Survival": {
        "success_default": "dead",
        "failure_default": "alive",
        "success_label": "Dead / affected",
        "failure_label": "Alive / unaffected",
        "raw_fraction": "mortality_fraction",
        "effect_fraction": "mortality_fraction",
        "effect_label": "Mortality / affected fraction",
    },
    "Motility": {
        "score_default": "motility_score",
        "raw_fraction": "normalized_motility",
        "effect_fraction": "motility_inhibition_fraction",
        "effect_label": "Motility inhibition",
    },
}


@dataclass
class FitResult:
    group_key: tuple
    bottom: float
    top: float
    log_ic50: float
    ic50: float
    hill: float
    converged: bool
    message: str
    n: int
    unique_doses: int
    ic50_ci_low: Optional[float] = None
    ic50_ci_high: Optional[float] = None
    ic50_bootstrap_samples: Optional[list[float]] = None

    def as_dict(self, group_names: Iterable[str]) -> dict:
        out = {name: value for name, value in zip(group_names, self.group_key)}
        out.update(
            {
                "IC50": self.ic50,
                "IC50_CI_low": self.ic50_ci_low,
                "IC50_CI_high": self.ic50_ci_high,
                "bottom": self.bottom,
                "top": self.top,
                "hill_slope": self.hill,
                "converged": self.converged,
                "message": self.message,
                "n_observations": self.n,
                "unique_doses": self.unique_doses,
            }
        )
        return out


def coerce_numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def validate_common_columns(df: pd.DataFrame, columns: Iterable[str]) -> list[str]:
    warnings: list[str] = []
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    for col in columns:
        if df[col].isna().all():
            warnings.append(f"Column '{col}' is entirely missing.")
    return warnings


def calculate_count_response(
    df: pd.DataFrame,
    success_col: str,
    failure_col: str,
    assay_name: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Calculate count-based fractions for egg hatch, LDA, or survival assays."""
    warnings = validate_common_columns(df, [success_col, failure_col])
    out = coerce_numeric(df, [success_col, failure_col])

    if out[[success_col, failure_col]].isna().any().any():
        warnings.append("Some count values could not be converted to numbers and were set to missing.")

    if (out[[success_col, failure_col]] < 0).any().any():
        raise ValueError("Count columns cannot contain negative values.")

    out["success_count"] = out[success_col]
    out["failure_count"] = out[failure_col]
    out["total_count"] = out["success_count"] + out["failure_count"]

    zero_total = out["total_count"] <= 0
    if zero_total.any():
        warnings.append(f"{zero_total.sum()} rows have zero total count and were removed.")
        out = out.loc[~zero_total].copy()

    assay = ASSAY_PRESETS[assay_name]
    raw_fraction = assay["raw_fraction"]
    effect_fraction = assay["effect_fraction"]

    out[raw_fraction] = out["success_count"] / out["total_count"]

    if assay_name in ["Egg hatch", "Larval development"]:
        out[effect_fraction] = 1 - out[raw_fraction]
    else:
        # Survival preset treats success as dead/affected, so raw fraction is already effect.
        out[effect_fraction] = out[raw_fraction]
        out["survival_fraction"] = 1 - out[raw_fraction]

    out["response_fraction"] = out[effect_fraction].clip(0, 1)
    out["response_percent"] = out["response_fraction"] * 100
    return out, warnings


def calculate_motility_response(
    df: pd.DataFrame,
    score_col: str,
    group_cols: list[str],
    dose_col: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Normalize motility score to the zero-dose control within each group."""
    warnings = validate_common_columns(df, [score_col, dose_col])
    out = coerce_numeric(df, [score_col, dose_col])

    if out[[score_col, dose_col]].isna().any().any():
        warnings.append("Some motility or dose values could not be converted to numbers and were set to missing.")

    if (out[score_col] < 0).any():
        warnings.append("Negative motility scores were detected. Please confirm these are valid.")

    # Build control mean by group. Usually group_cols = [strain, drug].
    control = out.loc[out[dose_col] == 0].groupby(group_cols, dropna=False)[score_col].mean()

    def get_control(row):
        key = tuple(row[c] for c in group_cols)
        try:
            return control.loc[key]
        except KeyError:
            return np.nan

    out["control_mean"] = out.apply(get_control, axis=1)
    missing_control = out["control_mean"].isna()
    if missing_control.any():
        warnings.append(
            f"{missing_control.sum()} rows lack a matching zero-dose control; motility normalization may be incomplete."
        )

    out["normalized_motility"] = out[score_col] / out["control_mean"]
    above_control = (out[dose_col] > 0) & (out["normalized_motility"] > 1)
    if above_control.any():
        warnings.append(
            f"{int(above_control.sum())} nonzero-dose motility rows had scores above the matched zero-dose control mean; "
            "negative inhibition was clipped to 0. Consider inspecting these rows for stimulation/hormesis."
        )
    out["motility_inhibition_fraction_unclipped"] = 1 - out["normalized_motility"]
    out["motility_inhibition_fraction"] = out["motility_inhibition_fraction_unclipped"].clip(0, 1)
    out["response_fraction"] = out["motility_inhibition_fraction"]
    out["response_percent"] = out["response_fraction"] * 100
    return out, warnings


def four_parameter_logistic(dose: np.ndarray, bottom: float, top: float, log_ic50: float, hill: float) -> np.ndarray:
    """Increasing four-parameter logistic function on log10 dose scale.

    Response approaches bottom at low doses and top at high doses.
    IC50 is returned on the original dose scale as 10 ** log_ic50.
    """
    dose = np.asarray(dose, dtype=float)
    positive = dose[dose > 0]
    if positive.size == 0:
        positive_floor = 1e-12
    else:
        positive_floor = positive.min() / 10.0
    dose_safe = np.where(dose <= 0, positive_floor, dose)
    x = np.log10(dose_safe)
    return bottom + (top - bottom) / (1 + 10 ** ((log_ic50 - x) * hill))


def _initial_parameters(x: np.ndarray, y: np.ndarray) -> list[float]:
    positive = x[x > 0]
    if positive.size == 0:
        log_ic50 = 0.0
    else:
        y_min = float(np.nanmin(y))
        y_max = float(np.nanmax(y))
        mid = y_min + (y_max - y_min) / 2
        idx = int(np.nanargmin(np.abs(y - mid)))
        guess_dose = x[idx] if x[idx] > 0 else np.nanmedian(positive)
        log_ic50 = float(np.log10(guess_dose))

    bottom = max(0.0, min(0.25, float(np.nanmin(y))))
    top = min(1.0, max(0.75, float(np.nanmax(y))))
    if top <= bottom:
        top = min(1.0, bottom + 0.2)
    return [bottom, top, log_ic50, 1.0]


def _fit_one_group(data: pd.DataFrame, dose_col: str, response_col: str, total_col: Optional[str] = None) -> tuple[np.ndarray, np.ndarray, str]:
    clean = data[[dose_col, response_col] + ([total_col] if total_col and total_col in data.columns else [])].dropna()
    x = clean[dose_col].astype(float).to_numpy()
    y = clean[response_col].astype(float).to_numpy()
    if len(x) < 4:
        raise RuntimeError("Fewer than four usable observations.")
    if len(np.unique(x)) < 4:
        raise RuntimeError("Fewer than four unique dose levels; four-parameter logistic fitting is unreliable.")
    if np.all(x <= 0):
        raise RuntimeError("No positive dose values available for IC50 fitting.")

    p0 = _initial_parameters(x, y)
    positive = x[x > 0]
    log_min = np.log10(positive.min()) - 3
    log_max = np.log10(positive.max()) + 3
    bounds = ([0.0, 0.0, log_min, 0.05], [1.0, 1.0, log_max, 10.0])

    sigma = None
    if total_col and total_col in clean.columns:
        # Approximate binomial standard error with a floor to avoid overweighting p=0 or p=1.
        n = clean[total_col].astype(float).clip(lower=1).to_numpy()
        p = np.clip(y, 0.01, 0.99)
        sigma = np.sqrt(p * (1 - p) / n)
        sigma = np.clip(sigma, 0.02, None)

    last_error = ""
    for hill_guess in [1.0, 0.5, 2.0, 4.0]:
        p0_try = p0.copy()
        p0_try[3] = hill_guess
        try:
            popt, pcov = curve_fit(
                four_parameter_logistic,
                x,
                y,
                p0=p0_try,
                bounds=bounds,
                sigma=sigma,
                absolute_sigma=False,
                maxfev=30000,
            )
            return popt, pcov, "OK"
        except Exception as exc:  # pragma: no cover - useful runtime fallback
            last_error = str(exc)
            continue
    raise RuntimeError(last_error or "Curve fit failed.")


def fit_dose_response(
    df: pd.DataFrame,
    group_cols: list[str],
    dose_col: str = "dose",
    response_col: str = "response_fraction",
    total_col: Optional[str] = "total_count",
    n_boot: int = 0,
    random_seed: int = 42,
) -> tuple[pd.DataFrame, dict[tuple, FitResult]]:
    """Fit a four-parameter logistic model within each group."""
    results: dict[tuple, FitResult] = {}
    rng = np.random.default_rng(random_seed)

    if not group_cols:
        group_cols = ["__all__"]
        df = df.copy()
        df["__all__"] = "All data"

    for group_key, group_data in df.groupby(group_cols, dropna=False):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)

        n = int(group_data[[dose_col, response_col]].dropna().shape[0])
        unique_doses = int(group_data[dose_col].dropna().nunique())
        try:
            popt, pcov, message = _fit_one_group(group_data, dose_col, response_col, total_col)
            bottom, top, log_ic50, hill = [float(v) for v in popt]
            ic50 = float(10 ** log_ic50)
            result = FitResult(group_key, bottom, top, log_ic50, ic50, hill, True, message, n, unique_doses)
            if np.isfinite(bottom) and np.isfinite(top) and top < bottom:
                result.message += (
                    "; fitted top is below fitted bottom, indicating a decreasing dose-response curve. "
                    "For inhibition/mortality endpoints this may indicate swapped response columns, wrong assay settings, or poor data quality."
                )

            if n_boot > 0:
                boot_ic50 = []
                for _ in range(n_boot):
                    idx = rng.integers(0, len(group_data), len(group_data))
                    boot = group_data.iloc[idx].copy()
                    try:
                        boot_popt, _, _ = _fit_one_group(boot, dose_col, response_col, total_col)
                        boot_ic50.append(float(10 ** boot_popt[2]))
                    except Exception:
                        continue
                if len(boot_ic50) >= max(20, n_boot * 0.25):
                    result.ic50_ci_low = float(np.percentile(boot_ic50, 2.5))
                    result.ic50_ci_high = float(np.percentile(boot_ic50, 97.5))
                    result.ic50_bootstrap_samples = [float(v) for v in boot_ic50]
                else:
                    result.message += "; bootstrap CI unavailable because too few bootstrap fits converged"

        except Exception as exc:
            result = FitResult(
                group_key=group_key,
                bottom=np.nan,
                top=np.nan,
                log_ic50=np.nan,
                ic50=np.nan,
                hill=np.nan,
                converged=False,
                message=str(exc),
                n=n,
                unique_doses=unique_doses,
            )
        results[group_key] = result

    summary = pd.DataFrame([r.as_dict(group_cols) for r in results.values()])
    return summary, results


def summarize_by_dose(
    df: pd.DataFrame,
    group_cols: list[str],
    dose_col: str = "dose",
    response_col: str = "response_fraction",
) -> pd.DataFrame:
    summary = (
        df.groupby(group_cols + [dose_col], dropna=False)
        .agg(
            mean_response=(response_col, "mean"),
            sd_response=(response_col, "std"),
            n=(response_col, "count"),
            mean_response_percent=("response_percent", "mean"),
            sd_response_percent=("response_percent", "std"),
        )
        .reset_index()
    )
    summary["se_response"] = summary["sd_response"] / np.sqrt(summary["n"].clip(lower=1))
    summary["se_response_percent"] = summary["sd_response_percent"] / np.sqrt(summary["n"].clip(lower=1))
    return summary


def _bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, preserving NaNs."""
    pvalues = np.asarray(pvalues, dtype=float)
    adjusted = np.full_like(pvalues, np.nan, dtype=float)
    mask = np.isfinite(pvalues)
    p = pvalues[mask]
    m = len(p)
    if m == 0:
        return adjusted
    order = np.argsort(p)
    ranked = p[order]
    raw_adj = ranked * m / np.arange(1, m + 1)
    monotone = np.minimum.accumulate(raw_adj[::-1])[::-1]
    monotone = np.clip(monotone, 0, 1)
    out = np.empty_like(monotone)
    out[order] = monotone
    adjusted[mask] = out
    return adjusted


def add_p_value_adjustments(test_table: pd.DataFrame, p_col: str = "p_value") -> pd.DataFrame:
    """Add Bonferroni and Benjamini-Hochberg adjusted p-values to a test table."""
    out = test_table.copy()
    if out.empty or p_col not in out.columns:
        return out
    p = pd.to_numeric(out[p_col], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(p)
    m = int(finite.sum())
    out["p_value_bonferroni"] = np.nan
    out.loc[finite, "p_value_bonferroni"] = np.clip(p[finite] * m, 0, 1)
    out["p_value_bh"] = _bh_adjust(p)
    out["p_value_adjustment_note"] = (
        "Adjusted over all finite p-values in this table; use raw p-values only as exploratory dose-level tests."
    )
    return out


def _approx_rr_ci_from_ic50_ci(row: pd.Series, ref_row: pd.Series, rr: float) -> tuple[float, float, str]:
    """Approximate fold-resistance CI from independent log-scale IC50 CIs."""
    try:
        vals = [
            float(row["IC50"]),
            float(row["IC50_CI_low"]),
            float(row["IC50_CI_high"]),
            float(ref_row["IC50"]),
            float(ref_row["IC50_CI_low"]),
            float(ref_row["IC50_CI_high"]),
        ]
    except Exception:
        return np.nan, np.nan, "not available"
    if not all(np.isfinite(v) and v > 0 for v in vals) or not np.isfinite(rr) or rr <= 0:
        return np.nan, np.nan, "not available"
    _, low, high, _, ref_low, ref_high = vals
    se_log = (np.log(high) - np.log(low)) / (2 * 1.96)
    ref_se_log = (np.log(ref_high) - np.log(ref_low)) / (2 * 1.96)
    se_log_rr = float(np.sqrt(se_log ** 2 + ref_se_log ** 2))
    low_rr = float(np.exp(np.log(rr) - 1.96 * se_log_rr))
    high_rr = float(np.exp(np.log(rr) + 1.96 * se_log_rr))
    return low_rr, high_rr, "approximate log-scale propagation from IC50 CIs"


def _bootstrap_rr_ci(test_result: FitResult, ref_result: FitResult, random_seed: int = 42) -> tuple[float, float, str]:
    """Estimate a fold-resistance CI from stored bootstrap IC50 samples."""
    test_samples = getattr(test_result, "ic50_bootstrap_samples", None)
    ref_samples = getattr(ref_result, "ic50_bootstrap_samples", None)
    if not test_samples or not ref_samples or len(test_samples) < 20 or len(ref_samples) < 20:
        return np.nan, np.nan, "not available"
    test = np.asarray(test_samples, dtype=float)
    ref = np.asarray(ref_samples, dtype=float)
    test = test[np.isfinite(test) & (test > 0)]
    ref = ref[np.isfinite(ref) & (ref > 0)]
    if len(test) < 20 or len(ref) < 20:
        return np.nan, np.nan, "not available"
    rng = np.random.default_rng(random_seed)
    n_ratio_draws = 10000
    ratios = rng.choice(test, size=n_ratio_draws, replace=True) / rng.choice(ref, size=n_ratio_draws, replace=True)
    return float(np.percentile(ratios, 2.5)), float(np.percentile(ratios, 97.5)), "bootstrap ratio from stored IC50 bootstrap samples"


def calculate_resistance_ratios(
    fit_summary: pd.DataFrame,
    group_col: str,
    reference_group: str,
    fit_results: Optional[dict[tuple, FitResult]] = None,
    group_cols: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Calculate fold-resistance versus a reference group, with optional uncertainty intervals.

    The point estimate is IC50_test / IC50_reference. If bootstrap IC50 samples are
    available, ARStat reports an empirical bootstrap CI for the ratio. Otherwise, if
    IC50 confidence limits are available, it reports an approximate log-scale CI.
    """
    if group_col not in fit_summary.columns:
        raise ValueError(f"Reference group column '{group_col}' is not in the fit table.")
    if group_cols is None:
        group_cols = [group_col]
    other_cols = [c for c in group_cols if c != group_col and c in fit_summary.columns]

    out_parts = []
    grouped = fit_summary.groupby(other_cols, dropna=False) if other_cols else [((), fit_summary)]
    for stratum, sub in grouped:
        sub = sub.copy()
        ref_rows = sub.loc[sub[group_col].astype(str) == str(reference_group)]
        if ref_rows.empty:
            sub["reference_group"] = reference_group
            sub["resistance_ratio_message"] = f"Reference group '{reference_group}' was not found in this stratum."
            sub["resistance_ratio"] = np.nan
            sub["fold_resistance_vs_reference"] = np.nan
            sub["fold_resistance_CI_low"] = np.nan
            sub["fold_resistance_CI_high"] = np.nan
            sub["fold_resistance_CI_method"] = "not available"
            out_parts.append(sub)
            continue
        ref_row = ref_rows.iloc[0]
        ref_ic50 = float(ref_row["IC50"]) if pd.notna(ref_row["IC50"]) else np.nan
        if not np.isfinite(ref_ic50) or ref_ic50 <= 0:
            sub["reference_group"] = reference_group
            sub["resistance_ratio_message"] = f"Reference group '{reference_group}' does not have a valid positive IC50."
            sub["resistance_ratio"] = np.nan
            sub["fold_resistance_vs_reference"] = np.nan
            sub["fold_resistance_CI_low"] = np.nan
            sub["fold_resistance_CI_high"] = np.nan
            sub["fold_resistance_CI_method"] = "not available"
            out_parts.append(sub)
            continue

        sub["reference_group"] = reference_group
        sub["resistance_ratio"] = sub["IC50"] / ref_ic50
        sub["fold_resistance_vs_reference"] = sub["resistance_ratio"]
        sub["fold_resistance_CI_low"] = np.nan
        sub["fold_resistance_CI_high"] = np.nan
        sub["fold_resistance_CI_method"] = "not available"
        sub["resistance_ratio_message"] = ""

        # Reference key for stored FitResult lookup.
        ref_key = tuple(ref_row[c] for c in group_cols if c in ref_row.index)
        ref_result = fit_results.get(ref_key) if fit_results else None

        for idx, row in sub.iterrows():
            rr = float(row["resistance_ratio"]) if pd.notna(row["resistance_ratio"]) else np.nan
            row_key = tuple(row[c] for c in group_cols if c in row.index)
            test_result = fit_results.get(row_key) if fit_results else None
            if test_result is not None and ref_result is not None:
                low, high, method = _bootstrap_rr_ci(test_result, ref_result)
            else:
                low, high, method = np.nan, np.nan, "not available"
            if not np.isfinite(low) or not np.isfinite(high):
                low, high, method = _approx_rr_ci_from_ic50_ci(row, ref_row, rr)
            sub.at[idx, "fold_resistance_CI_low"] = low
            sub.at[idx, "fold_resistance_CI_high"] = high
            sub.at[idx, "fold_resistance_CI_method"] = method

        out_parts.append(sub)

    return pd.concat(out_parts, ignore_index=True)


def assay_warnings(df: pd.DataFrame, group_cols: list[str], dose_col: str = "dose") -> list[str]:
    warnings: list[str] = []
    if dose_col not in df.columns:
        return [f"Dose column '{dose_col}' is missing."]
    dose_numeric = pd.to_numeric(df[dose_col], errors="coerce")
    if dose_numeric.isna().any():
        warnings.append("Some dose values are non-numeric or missing.")
    for group_key, group_data in df.assign(__dose=dose_numeric).groupby(group_cols, dropna=False):
        label = ", ".join(map(str, group_key if isinstance(group_key, tuple) else (group_key,)))
        doses = sorted(group_data["__dose"].dropna().unique())
        if 0 not in doses:
            warnings.append(f"{label}: missing zero-dose control.")
        if len(doses) < 4:
            warnings.append(f"{label}: fewer than four unique dose levels; IC50 may be unstable.")
        positive = [d for d in doses if d > 0]
        if len(positive) < 3:
            warnings.append(f"{label}: fewer than three positive dose levels.")
    return warnings


def pairwise_count_tests(
    df: pd.DataFrame,
    comparison_col: str,
    dose_col: str,
    stratify_cols: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Fisher exact tests for count-based assays at each dose.

    This compares success/failure counts between all pairs of comparison groups within
    each drug/dose stratum.
    """
    stratify_cols = stratify_cols or []
    required = [comparison_col, dose_col, "success_count", "failure_count"] + stratify_cols
    validate_common_columns(df, required)
    rows = []
    for stratum, sub in df.groupby(stratify_cols + [dose_col], dropna=False):
        if not isinstance(stratum, tuple):
            stratum = (stratum,)
        groups = sorted(sub[comparison_col].dropna().unique())
        for g1, g2 in combinations(groups, 2):
            a = sub.loc[sub[comparison_col] == g1, ["success_count", "failure_count"]].sum()
            b = sub.loc[sub[comparison_col] == g2, ["success_count", "failure_count"]].sum()
            table = [[a["success_count"], a["failure_count"]], [b["success_count"], b["failure_count"]]]
            try:
                odds_ratio, pvalue = fisher_exact(table)
            except Exception:
                odds_ratio, pvalue = np.nan, np.nan
            row = {name: value for name, value in zip(stratify_cols + [dose_col], stratum)}
            row.update({"group_1": g1, "group_2": g2, "odds_ratio": odds_ratio, "p_value": pvalue})
            rows.append(row)
    return add_p_value_adjustments(pd.DataFrame(rows))


def pairwise_continuous_tests(
    df: pd.DataFrame,
    comparison_col: str,
    dose_col: str,
    response_col: str = "response_fraction",
    stratify_cols: Optional[list[str]] = None,
    test: str = "mannwhitney",
) -> pd.DataFrame:
    """Pairwise per-dose tests for continuous responses, useful for motility scores."""
    stratify_cols = stratify_cols or []
    validate_common_columns(df, [comparison_col, dose_col, response_col] + stratify_cols)
    rows = []
    for stratum, sub in df.groupby(stratify_cols + [dose_col], dropna=False):
        if not isinstance(stratum, tuple):
            stratum = (stratum,)
        groups = sorted(sub[comparison_col].dropna().unique())
        for g1, g2 in combinations(groups, 2):
            y1 = sub.loc[sub[comparison_col] == g1, response_col].dropna().to_numpy()
            y2 = sub.loc[sub[comparison_col] == g2, response_col].dropna().to_numpy()
            if len(y1) < 2 or len(y2) < 2:
                statistic, pvalue = np.nan, np.nan
            elif test == "t-test":
                statistic, pvalue = ttest_ind(y1, y2, equal_var=False)
            else:
                statistic, pvalue = mannwhitneyu(y1, y2, alternative="two-sided")
            row = {name: value for name, value in zip(stratify_cols + [dose_col], stratum)}
            row.update({"group_1": g1, "group_2": g2, "statistic": statistic, "p_value": pvalue})
            rows.append(row)
    return add_p_value_adjustments(pd.DataFrame(rows))
