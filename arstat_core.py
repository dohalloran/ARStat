"""Core analysis functions for ARStat.

ARStat is designed for antinematodal / anthelmintic dose-response assays.
The Streamlit app imports these functions, but they can also be used directly
from Python scripts or notebooks.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import fisher_exact, mannwhitneyu, ttest_ind


INFO_PREFIX = "INFO: "


def dataframe_signature(df: pd.DataFrame) -> str:
    """Return a content-aware signature for Streamlit state/result invalidation.

    The signature includes column names and dtypes as well as row values/index, so
    a corrected re-upload with the same filename and dimensions cannot silently
    inherit stale auto-detection or stored results.
    """
    digest = hashlib.sha256()
    schema = [(str(col), str(dtype)) for col, dtype in zip(df.columns, df.dtypes)]
    digest.update(repr(schema).encode("utf-8"))
    try:
        row_hashes = pd.util.hash_pandas_object(df, index=True).to_numpy(dtype="uint64", copy=False)
        digest.update(row_hashes.tobytes())
    except Exception:
        digest.update(df.to_csv(index=True).encode("utf-8"))
    return f"{digest.hexdigest()}_{df.shape[0]}x{df.shape[1]}"


def read_table_bytes(raw: bytes, filename: str = "") -> pd.DataFrame:
    """Read an uploaded CSV or .xlsx file from bytes.

    CSV files are decoded as UTF-8 first (a byte-order mark is handled). CSVs
    saved by Excel on Windows are usually Windows-1252, where ``µ`` is byte
    0xB5 and is invalid UTF-8, so those files fall back to Windows-1252 rather
    than being rejected. Files that decode as neither are reported as not text.
    """
    from io import BytesIO

    if str(filename).lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(BytesIO(raw))
    try:
        return pd.read_csv(BytesIO(raw), encoding="utf-8")
    except UnicodeDecodeError:
        pass
    try:
        return pd.read_csv(BytesIO(raw), encoding="cp1252")
    except UnicodeDecodeError as exc:
        raise ValueError(
            "the file is not a readable text CSV; it could not be decoded as UTF-8 or Windows-1252"
        ) from exc


def normalize_dose_unit(value: object) -> str:
    """Normalize common molar dose-unit spellings for display and comparison.

    This treats ``uM``, ``μM`` and ``µM`` as equivalent without changing unknown
    units, whose capitalization may be meaningful.
    """
    text = str(value).strip().replace("μ", "µ")
    aliases = {
        "um": "µM",
        "µm": "µM",
        "nm": "nM",
        "mm": "mM",
        "pm": "pM",
    }
    return aliases.get(text.lower(), text)


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
    "Motility": {
        "measurement_default": "motility",
        "measurement_label": "Motility / activity measurement",
        "raw_fraction": "motility_fraction",
        "effect_fraction": "motility_inhibition_fraction",
        "effect_label": "Motility inhibition",
    },
}


ASSAY_COLUMN_ALIASES = {
    "Egg hatch": {
        "success": ["l1", "hatched", "hatched larvae", "larvae", "hatch"],
        "failure": ["eggs", "unhatched eggs", "unhatched", "egg"],
    },
    "Larval development": {
        "success": ["developed", "developed larvae", "l3", "l3 larvae", "infective larvae"],
        "failure": ["undeveloped", "undeveloped larvae", "l1", "l2", "l1 l2", "early larvae"],
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


def _canonical_column_name(value: object) -> str:
    """Normalize a column label for conservative alias matching."""
    text = str(value).strip().lower()
    return " ".join("".join(ch if ch.isalnum() else " " for ch in text).split())


def suggest_count_columns(columns: Iterable[str], assay_name: str) -> tuple[Optional[str], Optional[str]]:
    """Suggest assay count columns without overriding user confirmation.

    Matching is intentionally conservative and based on exact normalized aliases.
    In particular, ``L3`` is suggested as developed and ``L1`` as undeveloped for
    conventional larval-development assays, but the app still exposes both choices
    so investigators can confirm protocol-specific definitions.
    """
    aliases = ASSAY_COLUMN_ALIASES.get(assay_name)
    if not aliases:
        return None, None

    column_list = list(columns)
    normalized = {_canonical_column_name(col): col for col in column_list}

    def first_match(candidates: Iterable[str], excluded: Optional[str] = None) -> Optional[str]:
        for candidate in candidates:
            match = normalized.get(_canonical_column_name(candidate))
            if match is not None and match != excluded:
                return match
        return None

    success = first_match(aliases["success"])
    failure = first_match(aliases["failure"], excluded=success)
    return success, failure


ASSAY_VALUE_ALIASES = {
    "egg hatch": "Egg hatch",
    "egg_hatch": "Egg hatch",
    "egg-hatch": "Egg hatch",
    "eha": "Egg hatch",
    "larval development": "Larval development",
    "larval_development": "Larval development",
    "larval-development": "Larval development",
    "lda": "Larval development",
    "motility": "Motility",
}


def detect_declared_assays(df: pd.DataFrame) -> list[str]:
    """Return recognized assay names found in an ``assay`` metadata column."""
    assay_col = next((c for c in df.columns if _canonical_column_name(c) == "assay"), None)
    if assay_col is None:
        return []

    recognized: set[str] = set()
    for value in df[assay_col].dropna():
        key = str(value).strip().lower()
        normalized_key = _canonical_column_name(key)
        match = ASSAY_VALUE_ALIASES.get(key) or ASSAY_VALUE_ALIASES.get(normalized_key)
        if match:
            recognized.add(match)
    return sorted(recognized)


def detect_unrecognized_assay_labels(df: pd.DataFrame) -> list[str]:
    """Return values in an ``assay`` metadata column that ARStat does not support.

    For example, a file may declare a retired or otherwise unsupported assay;
    such files should not be run silently through another
    assay's preset.
    """
    assay_col = next((c for c in df.columns if _canonical_column_name(c) == "assay"), None)
    if assay_col is None:
        return []
    unknown: set[str] = set()
    for value in df[assay_col].dropna():
        key = str(value).strip().lower()
        if not key:
            continue
        if not (ASSAY_VALUE_ALIASES.get(key) or ASSAY_VALUE_ALIASES.get(_canonical_column_name(key))):
            unknown.add(str(value).strip())
    return sorted(unknown)


def infer_declared_assay(df: pd.DataFrame) -> Optional[str]:
    """Return one recognized assay declared by metadata, or ``None`` if ambiguous."""
    recognized = detect_declared_assays(df)
    return recognized[0] if len(recognized) == 1 else None


def detect_raw_assay_types(df: pd.DataFrame) -> list[str]:
    """Detect strong raw-assay signatures from assay-specific measurement columns.

    This intentionally ignores a generic ``assay`` metadata column by itself. The
    function is used to keep raw count/activity tables from being accidentally
    analyzed as normalized XY replicate tables.
    """
    detected: list[str] = []
    columns = list(df.columns)

    for assay_name in ("Egg hatch", "Larval development"):
        success, failure = suggest_count_columns(columns, assay_name)
        if success is not None and failure is not None:
            detected.append(assay_name)

    normalized_columns = {_canonical_column_name(c) for c in columns}
    # A single wide XY response column could legitimately be named "motility" or
    # "activity". Treat motility as a strong raw signature only when long-form
    # replicate/well metadata or an explicit motility assay declaration is also present.
    has_motility_measurement = any(name in normalized_columns for name in {"motility", "motility score", "activity", "activity score"})
    has_long_form_id = any(name in normalized_columns for name in {"replicate", "rep", "well", "well id"})
    if has_motility_measurement and (has_long_form_id or infer_declared_assay(df) == "Motility"):
        detected.append("Motility")

    return detected


# Column names that identify one row per replicate or well. A normalized XY
# table stores replicates as separate columns, so a file containing one of these
# identifier columns is long-form data and must not be analyzed as XY.
LONG_FORM_ID_NAMES = (
    "replicate", "rep", "replicate id", "replicate number", "rep id",
    "well", "well id", "wells", "well position",
)

# Additional identifier-like names that should never be offered as numeric
# response columns even when their values happen to be numbers.
ID_LIKE_NAMES = set(LONG_FORM_ID_NAMES) | {
    "experiment id", "experiment", "id", "sample id", "sample", "plate", "plate id",
    "row", "column", "index", "run", "batch", "date",
}

# Wide XY replicate columns such as Rep1, rep_2, Replicate 3, Y1.
XY_REPLICATE_PATTERN = re.compile(r"^(rep|replicate|y)\s?\d+$")


def find_column(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    """Return the first column matching a candidate name, ignoring case and punctuation."""
    lookup: dict[str, str] = {}
    for col in columns:
        lookup.setdefault(_canonical_column_name(col), col)
    for candidate in candidates:
        match = lookup.get(_canonical_column_name(candidate))
        if match is not None:
            return match
    return None


def is_id_like_column(column: object) -> bool:
    """Return True for replicate, well, and other identifier columns."""
    return _canonical_column_name(column) in ID_LIKE_NAMES


def is_xy_replicate_column(column: object) -> bool:
    """Return True for wide replicate column names such as Rep1 or Y2."""
    return bool(XY_REPLICATE_PATTERN.match(_canonical_column_name(column)))


def mostly_numeric(series: pd.Series, threshold: float = 0.8) -> bool:
    """Return True when at least ``threshold`` of supplied cells are numeric."""
    supplied = series.notna() & series.astype(str).str.strip().ne("")
    n_supplied = int(supplied.sum())
    if n_supplied == 0:
        return False
    n_numeric = int(pd.to_numeric(series.loc[supplied], errors="coerce").notna().sum())
    return n_numeric / n_supplied >= threshold


def detect_long_form_layout(df: pd.DataFrame) -> Optional[str]:
    """Explain why a table is long-form (one row per replicate/well), or return None.

    Raw assay files are recognised by their measurement columns in
    :func:`detect_raw_assay_types`. This check catches long-form files whose
    measurement column has a name ARStat does not recognise (for example a
    motility file with a ``thrashes`` column). Such files contain a replicate or
    well identifier column, which a normalized XY table never has.
    """
    id_col = find_column(df.columns, LONG_FORM_ID_NAMES)
    if id_col is not None:
        return (
            f"it has a per-row replicate/well identifier column ('{id_col}'), so each row is one replicate "
            "rather than one dose"
        )
    return None


def looks_like_normalized_xy(df: pd.DataFrame) -> bool:
    """Return True for wide tables with two or more RepN/YN replicate columns and no raw-assay signature."""
    rep_cols = [c for c in df.columns if is_xy_replicate_column(c)]
    if len(rep_cols) < 2:
        return False
    return not detect_raw_assay_types(df) and detect_long_form_layout(df) is None


def xy_repeated_dose_rows(
    df: pd.DataFrame,
    dose_col: str,
    group_col: Optional[str] = None,
    drug_col: Optional[str] = None,
) -> int:
    """Count rows whose dose repeats within the same group/drug in a wide XY table.

    A normalized XY table has one row per dose for each group/drug combination,
    with replicates in separate columns. Repeated doses mean either that the
    group/drug column has not been selected (several curves would be pooled) or
    that the file is long-form replicate data.
    """
    keys = [c for c in (group_col, drug_col) if c] + [dose_col]
    work = df[keys].copy()
    work[dose_col] = pd.to_numeric(work[dose_col], errors="coerce")
    work = work.dropna(subset=[dose_col])
    if work.empty:
        return 0
    return int(work.duplicated(subset=keys, keep="first").sum())


def find_duplicate_headers(columns: Iterable[str]) -> list[str]:
    """Return header names that were duplicated in the original file.

    ``pandas.read_csv`` and ``read_excel`` silently rename a repeated header
    ``L1`` to ``L1.1``, so the duplicates never reach the app under the same
    name. Detect that renaming so a second, different ``L1`` column cannot be
    ignored without warning.
    """
    column_list = [str(c) for c in columns]
    present = set(column_list)
    duplicated: list[str] = []
    for col in column_list:
        match = re.match(r"^(.*)\.(\d+)$", col)
        if match and match.group(1) in present and match.group(1) not in duplicated:
            duplicated.append(match.group(1))
    return duplicated


def count_columns_look_like_proportions(df: pd.DataFrame, success_col: str, failure_col: str) -> bool:
    """Return True when two 'count' columns are really proportions (every row totals at most 1)."""
    success = pd.to_numeric(df[success_col], errors="coerce")
    failure = pd.to_numeric(df[failure_col], errors="coerce")
    total = (success + failure).dropna()
    if total.empty:
        return False
    values = pd.concat([success, failure]).dropna()
    has_fraction = bool(((values % 1).abs() > 1e-9).any())
    return has_fraction and bool((total <= 1 + 1e-9).all())


def drop_entirely_empty_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Remove columns containing no usable values and report what was removed."""
    empty_cols = []
    for col in df.columns:
        series = df[col]
        nonempty = series.notna() & series.astype(str).str.strip().ne("")
        if not nonempty.any():
            empty_cols.append(col)
    if not empty_cols:
        return df.copy(), []
    cleaned = df.drop(columns=empty_cols).copy()
    note = f"Removed {len(empty_cols)} entirely empty column(s): {', '.join(map(str, empty_cols))}."
    return cleaned, [note]


def calculate_count_response(
    df: pd.DataFrame,
    success_col: str,
    failure_col: str,
    assay_name: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Calculate count-based fractions for egg hatch or larval-development assays."""
    if assay_name == "Motility":
        raise ValueError("Motility is a continuous-response assay; use prepare_motility_response().")
    if success_col == failure_col:
        raise ValueError(
            "The affected/success and unaffected/failure measurements must use different columns. "
            "Confirm the assay column mapping before analysis."
        )
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

    out[effect_fraction] = 1 - out[raw_fraction]

    out["response_fraction"] = out[effect_fraction].clip(0, 1)
    out["response_percent"] = out["response_fraction"] * 100
    return out, warnings


def prepare_motility_response(
    df: pd.DataFrame,
    dose_col: str,
    motility_col: str,
    group_cols: list[str],
    value_scale: str = "raw",
    response_direction: str = "raw_outcome",
) -> tuple[pd.DataFrame, list[str]]:
    """Prepare long-format motility measurements for IC50 fitting.

    Parameters
    ----------
    value_scale
        ``"raw"`` normalizes each observation to the mean zero-dose motility
        within its fitted group. ``"percent"`` treats values as 0--100 and
        ``"fraction"`` treats values as 0--1.
    response_direction
        ``"raw_outcome"`` means higher values represent greater motility and
        motility declines with dose. ``"effect"`` means imported values already
        represent motility inhibition or affected response.
    """
    required = [dose_col, motility_col] + list(group_cols)
    warnings = validate_common_columns(df, required)
    out = coerce_numeric(df, [dose_col, motility_col])

    before = len(out)
    out = out.dropna(subset=[dose_col, motility_col]).copy()
    removed = before - len(out)
    if removed:
        warnings.append(f"{removed} rows with missing or non-numeric dose/motility values were omitted.")
    if out.empty:
        raise ValueError("No usable motility measurements were found.")
    if (out[dose_col] < 0).any():
        raise ValueError("Dose values cannot be negative.")

    out["motility_measurement"] = out[motility_col].astype(float)

    if value_scale == "raw":
        if (out["motility_measurement"] < 0).any():
            raise ValueError("Raw motility/activity measurements cannot be negative.")

        control = out.loc[out[dose_col] == 0].copy()
        if control.empty:
            raise ValueError(
                "Raw motility measurements require at least one zero-dose control in every fitted group."
            )
        control_means = (
            control.groupby(group_cols, dropna=False)["motility_measurement"]
            .mean()
            .rename("control_mean_motility")
            .reset_index()
        )
        invalid_control = ~np.isfinite(control_means["control_mean_motility"]) | (
            control_means["control_mean_motility"] <= 0
        )
        if invalid_control.any():
            bad = control_means.loc[invalid_control, group_cols].apply(lambda col: col.map(str)).agg(" | ".join, axis=1)
            raise ValueError(
                "Zero-dose mean motility must be positive for every fitted group. Invalid groups: "
                + ", ".join(bad.tolist())
            )

        out = out.merge(control_means, on=group_cols, how="left", validate="many_to_one")
        missing_control = out["control_mean_motility"].isna()
        if missing_control.any():
            missing_groups = (
                out.loc[missing_control, group_cols]
                .drop_duplicates()
                .apply(lambda col: col.map(str))
                .agg(" | ".join, axis=1)
                .tolist()
            )
            raise ValueError(
                "Raw motility measurements require a zero-dose control for every fitted group. Missing groups: "
                + ", ".join(missing_groups)
            )
        imported_fraction = out["motility_measurement"] / out["control_mean_motility"]
        warnings.append("Raw motility measurements were normalized to the mean zero-dose control within each fitted group.")
    elif value_scale == "percent":
        out["control_mean_motility"] = np.nan
        imported_fraction = out["motility_measurement"] / 100.0
    elif value_scale == "fraction":
        out["control_mean_motility"] = np.nan
        imported_fraction = out["motility_measurement"]
    else:
        raise ValueError("value_scale must be 'raw', 'percent', or 'fraction'.")

    if response_direction not in {"raw_outcome", "effect"}:
        raise ValueError("response_direction must be 'raw_outcome' or 'effect'.")

    out["normalized_input_fraction"] = imported_fraction.astype(float)
    if (out["normalized_input_fraction"] < 0).any():
        raise ValueError("Motility/activity values cannot be negative after scaling or normalization.")

    if response_direction == "raw_outcome":
        raw = out["normalized_input_fraction"]
        above_control = raw > 1
        if above_control.any():
            warnings.append(
                f"{int(above_control.sum())} relative motility values exceeded the zero-dose mean (100%). "
                "They were retained because they may reflect hypermotility or ordinary replicate variation."
            )
        effect = 1 - raw
    else:
        outside = ~out["normalized_input_fraction"].between(0, 1)
        if outside.any():
            warnings.append(
                f"{int(outside.sum())} imported motility-inhibition values fell outside 0-1 and were retained. "
                "This can occur after control normalization or background correction; inspect these cells."
            )
        effect = out["normalized_input_fraction"]
        raw = 1 - effect

    out["motility_fraction"] = raw
    out["motility_inhibition_fraction"] = effect
    out["response_fraction"] = out["motility_inhibition_fraction"]
    out["response_percent"] = out["response_fraction"] * 100
    out["raw_outcome_percent"] = out["motility_fraction"] * 100
    return out, warnings



def _label_column_with_fallback(values: pd.Series, fallback: str) -> tuple[pd.Series, int]:
    """Return stripped text labels with blank or missing cells replaced by ``fallback``.

    ``Series.astype(str)`` turns missing cells into the literal strings ``"nan"`` or
    ``"None"`` under pandas 2.x, which would silently create a spurious group.
    Missing values are therefore detected before conversion to text.
    """
    text = values.astype("string").str.strip()
    missing = (text.isna() | text.eq("")).fillna(True).to_numpy(dtype=bool)
    labels = text.astype(object).to_numpy(copy=True)
    labels[missing] = fallback
    return pd.Series(labels, index=values.index).astype(str), int(missing.sum())


def prepare_normalized_xy_response(
    df: pd.DataFrame,
    dose_col: str,
    replicate_cols: list[str],
    assay_name: str,
    group_col: str | None = None,
    drug_col: str | None = None,
    dataset_label: str = "Dataset 1",
    drug_label: str = "Drug",
    unit: str = "",
    value_scale: str = "auto",
    response_direction: str = "raw_outcome",
) -> tuple[pd.DataFrame, list[str]]:
    """Convert a wide XY replicate table into ARStat's long analysis format.

    The input contains one X/dose column and one or more Y columns containing
    individual replicate responses. Optional group and drug columns permit
    several strains, isolates, genetic backgrounds, populations, or treatment
    groups to be analyzed in the same table. ARStat calculates dose-level mean,
    standard deviation, and sample size from the replicate values.

    Parameters
    ----------
    group_col
        Optional column identifying the experimental group. Values are mapped
        internally to ARStat's backward-compatible ``strain`` field.
    drug_col
        Optional column identifying the drug. When omitted, ``drug_label`` is
        applied to all rows.
    value_scale
        ``"auto"`` uses the overall magnitude of the imported values to
        distinguish percentages from fractions; ``"percent"`` divides values
        by 100; ``"fraction"`` leaves them on the 0--1 scale.
    response_direction
        ``"raw_outcome"`` means hatch, development, or motility decreases as dose increases. ``"effect"`` means the imported
        values already represent inhibition or affected response.
    """
    id_cols = [dose_col]
    for optional_col in (group_col, drug_col):
        if optional_col and optional_col not in id_cols:
            id_cols.append(optional_col)
    required = id_cols + list(replicate_cols)
    warnings = validate_common_columns(df, required)
    if not replicate_cols:
        raise ValueError("Select at least one replicate response column.")
    if len(replicate_cols) == 1:
        warnings.append(
            "Only one replicate response column was selected. ARStat can fit the curve, but dose-level "
            "variation and replicate-based comparisons cannot be estimated reliably."
        )
    overlap = set(replicate_cols).intersection(id_cols)
    if overlap:
        raise ValueError(f"Columns cannot be used as both identifiers and replicates: {sorted(overlap)}")

    wide = df.copy()
    wide = coerce_numeric(wide, [dose_col] + list(replicate_cols))
    if wide[dose_col].isna().any():
        warnings.append("Some dose values could not be converted to numbers and were removed.")

    long = wide[id_cols + list(replicate_cols)].melt(
        id_vars=id_cols,
        value_vars=list(replicate_cols),
        var_name="replicate",
        value_name="imported_response",
    )
    before = len(long)
    long = long.dropna(subset=[dose_col, "imported_response"]).copy()
    removed = before - len(long)
    if removed:
        warnings.append(f"{removed} missing or non-numeric dose/response cells were omitted.")
    if long.empty:
        raise ValueError("No usable normalized replicate values were found.")

    finite = long["imported_response"].to_numpy(dtype=float)
    scale = value_scale
    if scale == "auto":
        max_abs = float(np.nanmax(np.abs(finite)))
        median_abs = float(np.nanmedian(np.abs(finite)))
        if max_abs <= 1.5:
            scale = "fraction"
            detection_detail = "all values were compatible with a 0-1 fraction scale"
        elif max_abs >= 5 or median_abs > 1.5:
            scale = "percent"
            detection_detail = "the overall response magnitude was consistent with percentages"
        else:
            scale = "fraction"
            detection_detail = (
                "values between 1.5 and 5 are ambiguous, so ARStat conservatively assumed fractions; "
                "select the scale manually if these are percentages"
            )
        warnings.append(
            INFO_PREFIX
            + f"Normalized response scale was detected as {'0-100 percent' if scale == 'percent' else '0-1 fraction'}; "
            + detection_detail
            + "."
        )
    if scale == "percent":
        long["normalized_input_fraction"] = long["imported_response"] / 100.0
    elif scale == "fraction":
        long["normalized_input_fraction"] = long["imported_response"]
    else:
        raise ValueError("value_scale must be 'auto', 'percent', or 'fraction'.")

    if response_direction not in {"raw_outcome", "effect"}:
        raise ValueError("response_direction must be 'raw_outcome' or 'effect'.")

    outside = ~long["normalized_input_fraction"].between(0, 1)
    if outside.any():
        response_description = "raw-outcome" if response_direction == "raw_outcome" else "inhibition/affected-response"
        warnings.append(
            f"{int(outside.sum())} normalized {response_description} values fell outside 0-1 after scaling and were retained. "
            "ARStat fits the imported replicate values without clipping because control normalization, background "
            "correction, or ordinary replicate variation can produce values below 0% or above 100%."
        )

    if assay_name == "Motility" and response_direction == "raw_outcome":
        if (long["normalized_input_fraction"] < 0).any():
            warnings.append(
                "Negative relative-motility values were retained because the normalized table may include "
                "background-corrected measurements; confirm that the selected response scale is correct."
            )
        above_control = long["normalized_input_fraction"] > 1
        if above_control.any():
            warnings.append(
                f"{int(above_control.sum())} relative motility values exceeded 100% and were retained because "
                "they may reflect hypermotility or replicate variation."
            )

    if response_direction == "raw_outcome":
        raw = long["normalized_input_fraction"]
        effect = 1 - raw
    else:
        effect = long["normalized_input_fraction"]
        raw = 1 - effect

    raw_col = {
        "Egg hatch": "hatch_fraction",
        "Larval development": "development_fraction",
        "Motility": "motility_fraction",
    }[assay_name]
    effect_col = ASSAY_PRESETS[assay_name]["effect_fraction"]

    long = long.rename(columns={dose_col: "dose"})
    fallback = str(dataset_label).strip() or "Dataset 1"
    if group_col:
        long["strain"], n_missing_group = _label_column_with_fallback(long[group_col], fallback)
        if n_missing_group:
            warnings.append(f"{n_missing_group} rows had a missing group value and were assigned '{fallback}'.")
    else:
        long["strain"] = fallback

    fallback_drug = str(drug_label).strip() or "Drug"
    if drug_col:
        long["drug"], n_missing_drug = _label_column_with_fallback(long[drug_col], fallback_drug)
        if n_missing_drug:
            warnings.append(f"{n_missing_drug} rows had a missing drug value and were assigned '{fallback_drug}'.")
    else:
        long["drug"] = fallback_drug

    long["unit"] = unit
    long["assay"] = assay_name
    long[raw_col] = raw
    long[effect_col] = effect
    if assay_name == "Motility":
        long["motility_inhibition_fraction"] = effect
    long["response_fraction"] = effect
    long["response_percent"] = long["response_fraction"] * 100
    long["raw_outcome_percent"] = long[raw_col] * 100
    return long, warnings

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
                    "For inhibition/affected-response endpoints this may indicate swapped response columns, wrong assay settings, or poor data quality."
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
            warnings.append(
                f"{label}: only {len(doses)} unique dose level(s) were found. At least four are required for "
                "the four-parameter logistic model, so no IC50 will be calculated."
            )
        positive = [d for d in doses if d > 0]
        if len(positive) < 3:
            warnings.append(
                f"{label}: only {len(positive)} positive dose level(s) were found. Include at least three positive "
                "concentrations plus a zero-dose control for IC50 fitting."
            )
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
    """Pairwise per-dose tests for continuous responses, useful for normalized replicate responses."""
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
            row.update({"n_group_1": len(y1), "n_group_2": len(y2)})
            if test != "t-test":
                row["exact_min_p"] = mann_whitney_min_two_sided_p(len(y1), len(y2))
            rows.append(row)
    return add_p_value_adjustments(pd.DataFrame(rows))


def mann_whitney_min_two_sided_p(n1: int, n2: int) -> float:
    """Smallest two-sided P value an exact Mann-Whitney U test can return.

    With complete separation of the two samples, the exact two-sided P value is
    2 / C(n1 + n2, n1). For 3 versus 3 replicates this is 0.10, so an exact test
    can never reach P < 0.05 regardless of the effect size. (When values are
    tied, SciPy switches to a normal approximation, which can return smaller
    P values but is unreliable at such small sample sizes.)
    """
    if n1 < 1 or n2 < 1:
        return np.nan
    return float(min(1.0, 2.0 / comb(n1 + n2, n1)))
