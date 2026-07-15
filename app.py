"""ARStat web app.

Run locally:
    streamlit run app.py
"""

from __future__ import annotations

from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

from arstat_core import (
    ASSAY_PRESETS,
    assay_warnings,
    calculate_count_response,
    calculate_resistance_ratios,
    fit_dose_response,
    four_parameter_logistic,
    pairwise_continuous_tests,
    pairwise_count_tests,
    prepare_normalized_xy_response,
    summarize_by_dose,
)

APP_DIR = Path(__file__).parent
SAMPLE_DIR = APP_DIR / "sample_data"
TEMPLATE_DIR = APP_DIR / "templates"

st.set_page_config(
    page_title="ARStat",
    page_icon="🪱",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("ARStat")
st.caption("Reproducible statistical analysis and visualization of anthelmintic resistance dose-response assays")

with st.expander("What ARStat does", expanded=False):
    st.markdown(
        """
        ARStat turns raw assay counts or normalized replicate responses into standardized dose-response outputs.

        The bundled examples use illustrative hookworm datasets: Ancylostoma caninum WMD as the susceptible reference and KGR as a resistant isolate. Thiabendazole is used for the egg hatch example; ivermectin is used for larval development and survival/mortality examples.

        - assay-specific response calculations
        - input validation and warnings
        - four-parameter logistic IC50 estimation
        - optional bootstrap confidence intervals
        - fold-resistance estimates relative to a reference strain or isolate
        - fold-resistance confidence intervals when uncertainty estimates are available
        - per-dose pairwise tests with multiple-testing-adjusted p-values
        - downloadable CSV tables, Excel workbooks, and publication-ready PNG figures

        This version supports egg hatch, larval development, and survival/mortality assays.
        """
    )


def read_example(name: str) -> pd.DataFrame:
    return pd.read_csv(SAMPLE_DIR / name)


def read_uploaded_table(uploaded) -> pd.DataFrame:
    """Read a CSV or Excel worksheet uploaded through Streamlit."""
    name = getattr(uploaded, "name", "").lower()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(uploaded)
    return pd.read_csv(uploaded)


def dataframe_to_excel_bytes(tables: dict[str, pd.DataFrame]) -> bytes:
    output = BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        for sheet_name, table in tables.items():
            safe_name = sheet_name[:31]
            table.to_excel(writer, sheet_name=safe_name, index=False)
    return output.getvalue()


def csv_download(dataframe: pd.DataFrame) -> bytes:
    return dataframe.to_csv(index=False).encode("utf-8")


def png_plot_download(fig) -> bytes:
    output = BytesIO()
    fig.savefig(output, format="png", dpi=300, bbox_inches="tight")
    return output.getvalue()


def make_dose_response_plot(
    data: pd.DataFrame,
    fit_results: dict,
    group_cols: list[str],
    dose_col: str,
    response_col: str,
    y_label: str,
    dose_unit: str = "",
    plot_mode: str = "effect",
):
    """Create a publication-style dose-response plot using matplotlib.

    Zero-dose controls cannot be placed at x=0 on a logarithmic x-axis, so ARStat
    plots them at a symbolic left-edge position labelled "0". The fitted 4PL
    curve is still fit using the true dose values.
    """
    fig, ax = plt.subplots(figsize=(8, 5.5))

    label_col = "__plot_group"
    plot_data = data.copy()
    plot_data[dose_col] = pd.to_numeric(plot_data[dose_col], errors="coerce")
    plot_data[label_col] = plot_data[group_cols].astype(str).agg(" | ".join, axis=1)

    positive_all = plot_data.loc[plot_data[dose_col] > 0, dose_col].dropna().astype(float)
    if positive_all.empty:
        zero_plot_dose = 1e-6
        x_max_global = 1.0
    else:
        zero_plot_dose = float(positive_all.min()) / 10.0
        x_max_global = float(positive_all.max()) * 2.0

    plot_data["__plot_dose"] = plot_data[dose_col].where(plot_data[dose_col] > 0, zero_plot_dose)

    # Raw points. Zero-dose controls are shown at a symbolic left-edge tick.
    for group_key, group_data in plot_data.groupby(label_col, dropna=False):
        ax.scatter(
            group_data["__plot_dose"].astype(float),
            group_data[response_col].astype(float) * 100,
            alpha=0.75,
            label=f"{group_key} raw",
        )

    # Fitted 4PL curves. The model is always fit to effect response
    # (inhibition/mortality). In traditional raw-outcome mode, display the
    # complementary raw outcome so the curves decline with increasing dose.
    for result_key, result in fit_results.items():
        if not result.converged:
            continue
        group_label = " | ".join(map(str, result_key))
        sub = data.copy()
        for col, value in zip(group_cols, result_key):
            sub = sub.loc[sub[col].astype(str) == str(value)]
        positive_doses = sub.loc[sub[dose_col] > 0, dose_col].dropna().astype(float)
        if positive_doses.empty:
            continue
        x_min = max(float(positive_doses.min()) / 10.0, zero_plot_dose)
        x_max = float(positive_doses.max()) * 2.0
        x_curve = np.logspace(np.log10(x_min), np.log10(x_max), 250)
        y_curve = four_parameter_logistic(
            x_curve,
            result.bottom,
            result.top,
            result.log_ic50,
            result.hill,
        )
        if plot_mode == "raw_outcome":
            y_curve = 1 - y_curve
        ax.plot(x_curve, y_curve * 100, label=f"{group_label} 4PL fit")
        ax.axvline(result.ic50, linestyle="--", alpha=0.35)

    # Show a true zero-dose label despite using a log-scaled axis.
    positive_ticks = sorted(pd.unique(positive_all))
    if len(positive_ticks) > 8:
        keep_idx = np.linspace(0, len(positive_ticks) - 1, 8).round().astype(int)
        positive_ticks = [positive_ticks[i] for i in sorted(set(keep_idx))]
    xticks = [zero_plot_dose] + [float(v) for v in positive_ticks]
    xticklabels = ["0"] + [f"{float(v):g}" for v in positive_ticks]
    # Set the log scale before custom ticks; otherwise Matplotlib
    # can reset the ticker/formatter and overwrite the symbolic "0" label.
    ax.set_xscale("log")
    ax.set_xlim(zero_plot_dose / 1.8, x_max_global)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    unit_text = f" ({dose_unit})" if dose_unit else ""
    ax.set_xlabel(f"Drug concentration{unit_text}; zero-dose controls shown at symbolic left tick")
    ax.set_ylabel(y_label)
    ax.set_ylim(0, 105)
    if plot_mode == "raw_outcome":
        ax.set_title("Traditional IC50 curve: raw outcome declines with dose")
    else:
        ax.set_title("Dose-response curve: inhibition / affected response")
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return fig


def raw_plot_settings(assay_name: str) -> tuple[str, str]:
    """Return the column and y-axis label for traditional raw-outcome IC50 plots."""
    if assay_name == "Egg hatch":
        return "hatch_fraction", "Hatch rate (%)"
    if assay_name == "Larval development":
        return "development_fraction", "Development rate (%)"
    if assay_name == "Survival":
        return "survival_fraction", "Survival (%)"
    return "response_fraction", "Response (%)"


def make_method_text(
    assay_name: str,
    dose_col: str,
    group_cols: list[str],
    response_label: str,
    n_boot: int,
    dose_unit: str = "",
) -> str:
    bootstrap_sentence = (
        f"Bootstrap 95% confidence intervals for IC50 were estimated using {n_boot} within-group resamples. "
        "When available, fold-resistance confidence intervals were estimated from stored bootstrap IC50 samples; otherwise, approximate log-scale confidence intervals were calculated from IC50 confidence limits. "
        if n_boot > 0
        else "Bootstrap confidence intervals were not requested, so fold-resistance confidence intervals may be unavailable. "
    )
    unit_sentence = f" Doses were interpreted in {dose_unit}." if dose_unit else ""
    return (
        f"Dose-response data were analyzed using ARStat. For the {assay_name.lower()} assay, "
        f"raw measurements were converted to {response_label.lower()} and modeled as a function of "
        f"the concentration column '{dose_col}'."
        f"{unit_sentence} Curves were fit separately for each combination of "
        f"{', '.join(group_cols)} using a four-parameter logistic model. IC50 values were estimated "
        "as the midpoint between the fitted lower and upper asymptotes of the model; therefore, IC50 is relative to the fitted response range and is not necessarily the dose giving 50% absolute response when the fitted top or bottom differs from 100% or 0%. "
        f"{bootstrap_sentence}"
        "Fold resistance was calculated by dividing each IC50 estimate by the IC50 of the user-selected susceptible/reference group. "
        "Dose-level pairwise tests were reported as exploratory analyses with Benjamini-Hochberg and Bonferroni adjusted p-values. For count-based assays, Fisher exact tests pooled counts across replicate wells at each dose and therefore do not model replicate-to-replicate overdispersion."
    )

sample_options = {
    "Egg hatch example": "egg_hatch_example.csv",
    "Larval development example": "larval_development_example.csv",
    "Survival example": "survival_example.csv",
}

example_presets = {
    "Egg hatch example": {"assay": "Egg hatch", "reference": "WMD"},
    "Larval development example": {"assay": "Larval development", "reference": "WMD"},
    "Survival example": {"assay": "Survival", "reference": "WMD"},
}

template_files = {
    "Egg hatch template": "egg_hatch_template.csv",
    "Larval development template": "larval_development_template.csv",
    "Survival template": "survival_template.csv",
    "Normalized XY single-dataset template": "normalized_xy_replicates_template.csv",
    "Normalized XY multi-group template": "normalized_xy_multigroup_template.csv",
}


def choose_reference(values, preferred=None):
    """Pick a sensible reference group from the current dataset."""
    clean = [str(v) for v in values if pd.notna(v)]
    if preferred and str(preferred) in clean:
        return str(preferred)
    priority = ["WMD", "Susceptible", "susceptible", "Control", "control", "WT", "Wildtype", "N2"]
    for candidate in priority:
        if candidate in clean:
            return candidate
    return clean[0] if clean else ""


def show_download_library(location=st.sidebar):
    """Show download buttons for every sample and template CSV."""
    location.markdown("**Sample data**")
    for label, filename in sample_options.items():
        location.download_button(
            f"Download {label.lower()}",
            data=(SAMPLE_DIR / filename).read_bytes(),
            file_name=filename,
            mime="text/csv",
            key=f"sample-download-{filename}",
        )
    location.markdown("**Blank templates**")
    for label, filename in template_files.items():
        location.download_button(
            f"Download {label.lower()}",
            data=(TEMPLATE_DIR / filename).read_bytes(),
            file_name=filename,
            mime="text/csv",
            key=f"template-download-{filename}",
        )


st.sidebar.header("Navigation")
page = st.sidebar.radio("Page", ["Analyze data", "How to / user guide", "Why ARStat / comparison"], index=0)

if page == "Why ARStat / comparison":
    st.header("Why ARStat?")
    st.markdown(
        """
        ARStat is being built for a specific gap: parasitology labs often collect raw count data
        from egg hatch assays, larval development assays, and mortality/survival assays,
        but the analysis is commonly performed with a mixture of spreadsheets, general dose-response tools,
        and manually edited figures.

        **The problem ARStat solves**

        ARStat converts raw assay measurements into a consistent, reproducible analysis workflow:

        - assay-specific response calculations, such as hatch inhibition, development inhibition, or mortality
        - automatic data checks for missing controls, too few dose levels, zero counts, and non-numeric dose entries
        - IC50 estimation with a four-parameter logistic model
        - resistance ratios relative to a selected susceptible/control isolate
        - dose-level pairwise tests
        - downloadable prepared data, IC50 tables, Excel workbooks, figures, and methods text

        **How ARStat differs from general dose-response software**

        General tools can fit curves, but they usually expect users to calculate the biological response first.
        ARStat starts from raw parasitology assay counts or scores and applies assay-specific rules before fitting curves.

        | Tool type | Strength | Limitation that ARStat addresses |
        |---|---|---|
        | Commercial curve-fitting software | Excellent general curve fitting and figures | Usually requires users to structure and preprocess assay responses manually |
        | R drc / dr4pl | Powerful programmable dose-response modeling | Requires R coding and user-defined preprocessing |
        | General IC50 calculators | Simple and fast | Usually handle concentration/response pairs, not raw assay counts |
        | eggCounts / FECRT tools | Strong for fecal egg count reduction tests | Focused on FECRT rather than in vitro dose-response assays |
        | ARStat | Assay-aware, no-code, reproducible reports for anthelmintic phenotyping | Current version uses 4PL fitting; future versions may add replicate-aware beta-binomial or mixed-effects models |

        The main contribution is not that ARStat invents a new logistic equation. The main contribution is that
        it standardizes a messy, under-supported analysis workflow used in anthelmintic resistance phenotyping.
        """
    )
    st.stop()

if page == "How to / user guide":
    st.header("How to use ARStat")
    st.markdown(
        """
        ARStat is designed to take raw assay measurements and return a reproducible dose-response analysis.

        **Basic workflow**

        1. Choose an example dataset or upload your own CSV.
        2. Confirm the assay type. Example datasets automatically select the correct assay settings.
        3. Check that ARStat detected the correct strain, drug, dose, replicate, and count/score columns.
        4. Choose the reference strain or isolate for resistance-ratio calculations.
        5. Click **Run ARStat**.
        6. Inspect the fitted curve, IC50 table, dose summary, pairwise tests, and warnings.
        7. Download the Excel workbook and figure for your records or reports.

        **Important interpretation note**

        ARStat fits a four-parameter logistic curve to the assay effect. A right-shifted curve usually indicates reduced drug sensitivity.
        The fit should always be inspected visually, especially when there are few dose levels, no clear lower/upper plateau, or high replicate variability.
        """
    )

    st.subheader("Required columns")
    required = pd.DataFrame(
        [
            {"Assay": "Egg hatch", "Required measurements": "L1 and eggs", "Default response": "Hatch inhibition = 1 - L1/(L1 + eggs)"},
            {"Assay": "Larval development", "Required measurements": "developed and undeveloped", "Default response": "Development inhibition = 1 - developed/(developed + undeveloped)"},
            {"Assay": "Survival", "Required measurements": "dead and alive", "Default response": "Mortality / affected fraction = dead/(dead + alive)"},
        ]
    )
    st.dataframe(required, width='stretch')

    st.subheader("Download sample data and blank templates")
    show_download_library(location=st)

    st.subheader("Recommended experimental design")
    st.markdown(
        """
        - Include a zero-dose or vehicle control for each strain/isolate and drug.
        - Use at least 4 dose levels; 6 to 8 is better for stable IC50 estimation.
        - Include at least 3 to 4 biological or technical replicates per dose.
        - Make sure the tested dose range spans both weak and strong response.
        - Use the same dose units within one uploaded file.
        """
    )
    st.stop()

st.sidebar.header("1. Data")
input_layout = st.sidebar.radio(
    "Input layout",
    ["Raw assay measurements", "Normalized XY replicate table"],
    index=0,
    help=(
        "Raw measurements use assay-specific counts or scores. The normalized XY layout uses one dose column "
        "and adjacent columns containing individual replicate responses; ARStat calculates mean, SD, and n internally."
    ),
)

if input_layout == "Raw assay measurements":
    source = st.sidebar.radio("Data source", ["Use example data", "Upload CSV"], index=0)
    if source == "Use example data":
        sample_label = st.sidebar.selectbox("Example", list(sample_options.keys()), key="sample_label")
        if st.session_state.get("last_sample_label") != sample_label:
            preset = example_presets.get(sample_label, {})
            st.session_state["assay_type"] = preset.get("assay", "Egg hatch")
            st.session_state["reference_group_select"] = preset.get("reference", "None") or "None"
            st.session_state["last_sample_label"] = sample_label
        df = read_example(sample_options[sample_label])
    else:
        uploaded = st.sidebar.file_uploader("Upload CSV or Excel file", type=["csv", "xlsx"], key="raw_upload")
        if uploaded is None:
            st.info("Upload a CSV file or switch to an example dataset.")
            st.stop()
        df = read_uploaded_table(uploaded)
        sample_label = "uploaded"
        st.session_state["last_sample_label"] = None
else:
    source = "Upload CSV"
    sample_label = "normalized_xy"
    uploaded = st.sidebar.file_uploader(
        "Upload normalized XY replicate table",
        type=["csv", "xlsx"],
        key="normalized_xy_upload",
        help="Use one dose/X column and one or more individual replicate/Y columns. Do not precompute mean, SD, or n.",
    )
    if uploaded is None:
        st.info("Upload a normalized XY replicate table or switch to raw assay measurements.")
        st.stop()
    df = read_uploaded_table(uploaded)
    st.session_state["last_sample_label"] = None

with st.sidebar.expander("Download sample data/templates", expanded=False):
    show_download_library(location=st.sidebar)

cols = list(df.columns)
if df.empty:
    st.error("The selected dataset is empty. Upload a CSV with at least one data row.")
    st.stop()
if len(cols) != len(set(cols)):
    duplicated = sorted({c for c in cols if cols.count(c) > 1})
    st.error(f"Duplicate column names detected: {duplicated}. Please rename duplicate columns before analysis.")
    st.stop()

st.sidebar.header("2. Assay settings")
assay_names = list(ASSAY_PRESETS.keys())
default_assay = st.session_state.get("assay_type", "Egg hatch")
if default_assay not in assay_names:
    st.session_state["assay_type"] = "Egg hatch"
assay_name = st.sidebar.selectbox("Assay type", assay_names, key="assay_type")
assay = ASSAY_PRESETS[assay_name]

def default_col(name: str, fallback_index: int = 0) -> str:
    if name in cols:
        return name
    return cols[fallback_index] if cols else ""

if input_layout == "Raw assay measurements":
    strain_col = st.sidebar.selectbox(
        "Strain / isolate column", cols, index=cols.index(default_col("strain")) if "strain" in cols else 0
    )
    dose_col = st.sidebar.selectbox(
        "Dose column", cols, index=cols.index(default_col("dose")) if "dose" in cols else min(1, len(cols) - 1)
    )
    replicate_col = st.sidebar.selectbox(
        "Replicate / well column", ["None"] + cols, index=(["None"] + cols).index("replicate") if "replicate" in cols else 0
    )
    drug_col = st.sidebar.selectbox(
        "Drug column", ["None"] + cols, index=(["None"] + cols).index("drug") if "drug" in cols else 0
    )
    include_drug_in_groups = drug_col != "None"
    group_cols = [drug_col, strain_col] if include_drug_in_groups else [strain_col]
    dose_unit = st.sidebar.text_input(
        "Dose unit", value="µM", help="Used for plot labels and generated methods text."
    )
    success_default = assay.get("success_default", cols[0])
    failure_default = assay.get("failure_default", cols[0])
    success_col = st.sidebar.selectbox(
        assay.get("success_label", "Success count"), cols,
        index=cols.index(default_col(success_default)) if success_default in cols else 0,
    )
    failure_col = st.sidebar.selectbox(
        assay.get("failure_label", "Failure count"), cols,
        index=cols.index(default_col(failure_default)) if failure_default in cols else 0,
    )
    strain_values = sorted([str(v) for v in df[strain_col].dropna().unique()]) if strain_col in df.columns else []
else:
    dose_col = st.sidebar.selectbox(
        "Dose / X column",
        cols,
        index=cols.index("dose") if "dose" in cols else 0,
    )
    identifier_options = ["None"] + [c for c in cols if c != dose_col]
    group_candidates = ["group", "Group", "strain", "Strain", "isolate", "Isolate", "genetic_background", "population", "treatment"]
    detected_group = next((c for c in group_candidates if c in cols and c != dose_col), "None")
    group_col_input = st.sidebar.selectbox(
        "Experimental group column (optional)",
        identifier_options,
        index=identifier_options.index(detected_group),
        help="Examples include strain, isolate, genetic background, population, or treatment group.",
    )
    drug_candidates = ["drug", "Drug", "compound", "Compound"]
    detected_drug = next((c for c in drug_candidates if c in cols and c not in {dose_col, group_col_input}), "None")
    drug_options = ["None"] + [c for c in cols if c not in {dose_col, group_col_input}]
    drug_col_input = st.sidebar.selectbox(
        "Drug / compound column (optional)",
        drug_options,
        index=drug_options.index(detected_drug) if detected_drug in drug_options else 0,
    )
    excluded = {dose_col}
    if group_col_input != "None":
        excluded.add(group_col_input)
    if drug_col_input != "None":
        excluded.add(drug_col_input)
    available_y = [c for c in cols if c not in excluded]
    likely_replicates = [c for c in available_y if str(c).lower().startswith(("rep", "y"))]
    replicate_cols = st.sidebar.multiselect(
        "Replicate / Y columns",
        available_y,
        default=likely_replicates or available_y,
        help="Select individual replicate response columns. ARStat calculates mean, SD, and n internally.",
    )
    dataset_label = st.sidebar.text_input(
        "Single-dataset group label",
        value="Dataset 1",
        help="Used only when no experimental group column is selected.",
    )
    drug_label = st.sidebar.text_input(
        "Single-dataset drug label",
        value="Drug",
        help="Used only when no drug/compound column is selected.",
    )
    normalized_group_col = None if group_col_input == "None" else group_col_input
    normalized_drug_col = None if drug_col_input == "None" else drug_col_input
    dose_unit = st.sidebar.text_input("Dose unit", value="µM")
    normalized_scale_label = st.sidebar.selectbox(
        "Response scale", ["Auto-detect", "0–100 percent", "0–1 fraction"], index=0
    )
    normalized_scale = {"Auto-detect": "auto", "0–100 percent": "percent", "0–1 fraction": "fraction"}[normalized_scale_label]
    normalized_direction_label = st.sidebar.radio(
        "Imported response represents",
        ["Raw outcome (e.g., % hatched; decreases with dose)", "Inhibition / affected response (increases with dose)"],
        index=0,
    )
    normalized_direction = "raw_outcome" if normalized_direction_label.startswith("Raw outcome") else "effect"
    strain_col = "strain"
    drug_col = "drug"
    replicate_col = "replicate"
    include_drug_in_groups = True
    group_cols = [drug_col, strain_col]
    if normalized_group_col:
        strain_values = sorted(str(v) for v in df[normalized_group_col].dropna().unique())
    else:
        strain_values = [dataset_label] if dataset_label else ["Dataset 1"]

st.sidebar.header("3. Model settings")
n_boot = st.sidebar.select_slider(
    "Bootstrap IC50 confidence intervals",
    options=[0, 50, 100, 250, 500],
    value=0,
    help="More bootstrap resamples give more stable CIs but run more slowly.",
)

plot_display = st.sidebar.radio(
    "Dose-response plot display",
    [
        "Traditional raw outcome (descending IC50 curve)",
        "Inhibition / affected response (ascending curve)",
    ],
    index=0,
    help="The IC50 model is fit to inhibition/mortality. The traditional view displays the complementary raw assay outcome, such as hatch rate or survival, so the curve declines with increasing dose.",
)

preferred_ref = st.session_state.get("reference_group_select", None)
default_ref = "" if input_layout == "Normalized XY replicate table" else choose_reference(strain_values, preferred=preferred_ref)
reference_options = ["None"] + strain_values
if default_ref and default_ref not in reference_options:
    reference_options.append(default_ref)

# Keep the widget state valid when switching examples or uploading new files, but
# do not also pass an explicit default index to the selectbox. Passing both a
# widget default and a Session State value causes Streamlit's warning:
# "widget ... was created with a default value but also had its value set".
if st.session_state.get("reference_group_select") not in reference_options:
    st.session_state["reference_group_select"] = default_ref if default_ref in reference_options else "None"

reference_group = st.sidebar.selectbox(
    "Reference strain / isolate for resistance ratios",
    reference_options,
    key="reference_group_select",
    help="Choose the susceptible/control strain when calculating resistance ratios.",
)
if reference_group == "None":
    reference_group = ""

if not reference_group and len(strain_values) > 1:
    st.sidebar.info("Select a susceptible/control reference strain to calculate fold resistance.")

run_button = st.sidebar.button("Run ARStat", key="run_arstat_btn", type="primary")
if st.sidebar.button("Clear stored results", key="clear_results_btn"):
    st.session_state.pop("arstat_results", None)
    st.session_state.pop("arstat_config", None)
    st.rerun()

# A stable signature lets ARStat know whether the displayed results match the
# current inputs. Streamlit reruns the script after downloads and table clicks;
# this prevents expensive model fitting from running again unless the user
# explicitly clicks Run ARStat or changes settings.
def dataframe_signature(dataframe: pd.DataFrame) -> str:
    try:
        hashed = pd.util.hash_pandas_object(dataframe, index=True).astype("uint64")
        return str(int(hashed.sum())) + f"_{dataframe.shape[0]}x{dataframe.shape[1]}"
    except Exception:
        return f"{dataframe.shape}_{list(dataframe.columns)}"

current_config = {
    "input_layout": input_layout,
    "data_source": source,
    "sample_label": sample_label if source == "Use example data" else "uploaded",
    "data_signature": dataframe_signature(df),
    "assay_name": assay_name,
    "dose_col": dose_col,
    "dose_unit": dose_unit,
    "reference_group": reference_group,
    "n_boot": n_boot,
    "plot_display": plot_display,
}
if input_layout == "Raw assay measurements":
    current_config.update({
        "strain_col": strain_col,
        "replicate_col": replicate_col,
        "drug_col": drug_col,
        "include_drug_in_groups": include_drug_in_groups,
        "group_cols": group_cols,
    })
    current_config["success_col"] = success_col
    current_config["failure_col"] = failure_col
else:
    current_config.update({
        "replicate_cols": tuple(replicate_cols),
        "normalized_group_col": normalized_group_col,
        "normalized_drug_col": normalized_drug_col,
        "dataset_label": dataset_label,
        "drug_label": drug_label,
        "normalized_scale": normalized_scale,
        "normalized_direction": normalized_direction,
    })

st.subheader("Input preview")
st.caption(f"Loaded {len(df):,} rows and {len(df.columns):,} columns. Confirm column mappings in the sidebar before running the analysis.")
if source == "Use example data":
    st.success(f"Loaded {sample_label}. Assay settings were set to **{assay_name}** automatically. The bundled examples are illustrative hookworm-style sample data, not primary experimental measurements.")
st.dataframe(df.head(20), width='stretch')


def compute_arstat_results():
    """Run the analysis once and return display-ready results.

    The returned dict is stored in Streamlit Session State. Download buttons
    still trigger Streamlit reruns, but reruns reuse this dict rather than
    refitting the model.
    """
    working = df.copy()
    working[dose_col] = pd.to_numeric(working[dose_col], errors="coerce")

    if input_layout == "Normalized XY replicate table":
        analysis_df, response_warnings = prepare_normalized_xy_response(
            working,
            dose_col=dose_col,
            replicate_cols=replicate_cols,
            assay_name=assay_name,
            group_col=normalized_group_col,
            drug_col=normalized_drug_col,
            dataset_label=dataset_label,
            drug_label=drug_label,
            unit=dose_unit,
            value_scale=normalized_scale,
            response_direction=normalized_direction,
        )
        analysis_dose_col = "dose"
    else:
        analysis_df, response_warnings = calculate_count_response(
            working,
            success_col=success_col,
            failure_col=failure_col,
            assay_name=assay_name,
        )
        analysis_dose_col = dose_col

    validation_warnings = assay_warnings(analysis_df, group_cols=group_cols, dose_col=analysis_dose_col)
    all_warnings = response_warnings + validation_warnings

    statistical_notes = [
        "Zero-dose controls are included in model fitting and summary tables. Because log-scale plots cannot display a true dose of 0, zero-dose controls are shown at a symbolic left-edge tick labelled 0.",
        "Pairwise dose-level tests are exploratory. Raw p-values are reported together with Benjamini-Hochberg and Bonferroni adjusted p-values.",
    ]
    if input_layout == "Raw assay measurements":
        statistical_notes.append(
            "For count-based assays, Fisher exact tests pool replicate counts at each dose. This approach does not model replicate-to-replicate overdispersion and may be anticonservative when wells are highly variable."
        )
    statistical_notes.append(
        "IC50 values are estimated from the fitted dose-response curve and represent the dose corresponding to the midpoint between the fitted lower and upper asymptotes. Therefore, if the maximum fitted response is below 100%, the IC50 is not necessarily the dose producing an absolute 50% response."
    )

    fit_summary, fit_results = fit_dose_response(
        analysis_df,
        group_cols=group_cols,
        dose_col=analysis_dose_col,
        response_col="response_fraction",
        total_col=("total_count" if input_layout == "Raw assay measurements" else None),
        n_boot=n_boot,
    )

    if "message" in fit_summary.columns:
        decreasing_fit_messages = fit_summary.loc[
            fit_summary["message"].astype(str).str.contains("fitted top is below fitted bottom", na=False),
            group_cols + ["message"],
        ]
        for _, row in decreasing_fit_messages.iterrows():
            label = ", ".join(str(row[col]) for col in group_cols)
            all_warnings.append(f"{label}: {row['message']}")

    dose_summary = summarize_by_dose(
        analysis_df,
        group_cols=group_cols,
        dose_col=analysis_dose_col,
        response_col="response_fraction",
    )

    rr_table = pd.DataFrame()
    rr_message = ""
    if reference_group:
        try:
            rr_table = calculate_resistance_ratios(
                fit_summary,
                group_col=strain_col,
                reference_group=reference_group,
                fit_results=fit_results,
                group_cols=group_cols,
            )
        except Exception as exc:
            rr_message = f"Resistance ratios could not be calculated: {exc}"

    try:
        if input_layout == "Normalized XY replicate table":
            test_table = pairwise_continuous_tests(
                analysis_df,
                comparison_col=strain_col,
                dose_col=analysis_dose_col,
                response_col="response_fraction",
                stratify_cols=[drug_col] if include_drug_in_groups else [],
            )
        else:
            test_table = pairwise_count_tests(
                analysis_df,
                comparison_col=strain_col,
                dose_col=analysis_dose_col,
                stratify_cols=[drug_col] if include_drug_in_groups else [],
            )
        test_message = ""
    except Exception as exc:
        test_table = pd.DataFrame({"message": [f"Pairwise tests unavailable: {exc}"]})
        test_message = test_table.loc[0, "message"]

    if plot_display.startswith("Traditional raw outcome"):
        plot_response_col, response_label = raw_plot_settings(assay_name)
        plot_mode = "raw_outcome"
    else:
        plot_response_col = "response_fraction"
        response_label = assay.get("effect_label", "Response") + " (%)"
        plot_mode = "effect"

    fig = make_dose_response_plot(
        analysis_df,
        fit_results=fit_results,
        group_cols=group_cols,
        dose_col=analysis_dose_col,
        response_col=plot_response_col,
        y_label=response_label,
        dose_unit=dose_unit,
        plot_mode=plot_mode,
    )
    plot_png = png_plot_download(fig)
    plt.close(fig)

    methods_text = make_method_text(
        assay_name=assay_name,
        dose_col=analysis_dose_col,
        group_cols=group_cols,
        response_label=assay.get("effect_label", "response"),
        n_boot=n_boot,
        dose_unit=dose_unit,
    )
    if input_layout == "Normalized XY replicate table":
        methods_text += (
            " Normalized responses were imported from a wide XY replicate table containing one dose column, "
            "optional experimental-group and drug columns, and individual replicate response columns. ARStat "
            "reshaped the table to long format, fitted one curve per drug-by-group combination, and calculated "
            "dose-level means, standard deviations, and sample sizes from the replicate values."
        )

    all_tables = {
        "prepared_data": analysis_df,
        "dose_summary": dose_summary,
        "ic50_summary": fit_summary,
        "resistance_ratios": rr_table if not rr_table.empty else pd.DataFrame(),
        "pairwise_tests": test_table,
    }
    excel_bytes = dataframe_to_excel_bytes(all_tables)

    converged_count = int(fit_summary["converged"].sum()) if "converged" in fit_summary else 0

    return {
        "analysis_df": analysis_df,
        "fit_summary": fit_summary,
        "dose_summary": dose_summary,
        "rr_table": rr_table,
        "rr_message": rr_message,
        "test_table": test_table,
        "test_message": test_message,
        "all_warnings": all_warnings,
        "statistical_notes": statistical_notes,
        "methods_text": methods_text,
        "plot_png": plot_png,
        "excel_bytes": excel_bytes,
        "converged_count": converged_count,
        "fit_group_count": len(fit_summary),
        "rows_analyzed": len(analysis_df),
        "dose_levels": analysis_df[analysis_dose_col].nunique(),
        "response_label": response_label,
        "plot_mode": plot_mode,
    }


def render_arstat_results(results: dict):
    if results["all_warnings"]:
        with st.expander("Data checks and warnings", expanded=True):
            for warning in results["all_warnings"]:
                st.warning(warning)

    if results.get("statistical_notes"):
        with st.expander("Statistical notes and limitations", expanded=True):
            for note in results["statistical_notes"]:
                st.info(note)

    st.subheader("Results")
    c1, c2, c3 = st.columns(3)
    c1.metric("Groups fitted", f"{results['converged_count']}/{results['fit_group_count']}")
    c2.metric("Rows analyzed", results["rows_analyzed"])
    c3.metric("Dose levels", results["dose_levels"])

    if results.get("plot_mode") == "raw_outcome":
        plot_caption = (
            "Dose-response curve. In traditional raw-outcome mode, the fitted "
            "inhibition/mortality model is displayed as the complementary raw "
            "assay outcome, so curves decline with increasing dose. Zero-dose "
            "controls are shown at the symbolic left tick labelled 0 because "
            "log-scaled axes cannot display x=0."
        )
    else:
        plot_caption = (
            "Dose-response curve. In inhibition/affected-response mode, the "
            "fitted response increases with dose. Zero-dose controls are shown "
            "at the symbolic left tick labelled 0 because log-scaled axes cannot "
            "display x=0."
        )
    st.image(results["plot_png"], caption=plot_caption)

    st.markdown("### IC50 estimates")
    st.dataframe(results["fit_summary"], width='stretch')

    if results["rr_message"]:
        st.warning(results["rr_message"])
    elif not results["rr_table"].empty:
        st.markdown("### Fold resistance vs reference")
        st.info("Fold resistance = IC50 of each group divided by the IC50 of the selected reference group. Confidence intervals are reported when bootstrap IC50 samples or IC50 confidence limits are available.")
        st.dataframe(results["rr_table"], width='stretch')

    st.markdown("### Dose-level summary")
    st.dataframe(results["dose_summary"], width='stretch')

    st.markdown("### Pairwise per-dose tests")
    st.info(
        "These tests are exploratory. Use the adjusted p-values for multiple comparisons. "
        "For count-based assays, Fisher exact tests pool replicate wells at each dose and do not model overdispersion."
    )
    if results["test_message"]:
        st.info(results["test_message"])
    st.dataframe(results["test_table"], width='stretch')

    st.markdown("### Suggested methods text")
    st.text_area("Copy/edit this for a methods section", value=results["methods_text"], height=130)

    st.download_button(
        "Download methods text",
        data=results["methods_text"].encode("utf-8"),
        file_name="arstat_methods_text.txt",
        mime="text/plain",
        key="download-methods-text",
        on_click="ignore",
    )

    st.subheader("Downloads")
    c1, c2, c3, c4 = st.columns(4)
    c1.download_button(
        "Prepared data CSV",
        data=csv_download(results["analysis_df"]),
        file_name="arstat_prepared_data.csv",
        mime="text/csv",
        key="download-prepared-data-csv",
        on_click="ignore",
    )
    c2.download_button(
        "IC50 table CSV",
        data=csv_download(results["fit_summary"]),
        file_name="arstat_ic50_summary.csv",
        mime="text/csv",
        key="download-ic50-table-csv",
        on_click="ignore",
    )
    c3.download_button(
        "All tables Excel",
        data=results["excel_bytes"],
        file_name="arstat_results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="download-all-tables-excel",
        on_click="ignore",
    )
    c4.download_button(
        "Plot PNG",
        data=results["plot_png"],
        file_name="arstat_dose_response_plot.png",
        mime="image/png",
        key="download-plot-png",
        on_click="ignore",
    )


if run_button:
    try:
        with st.spinner("Running ARStat..."):
            st.session_state["arstat_results"] = compute_arstat_results()
            st.session_state["arstat_config"] = current_config
    except Exception as exc:
        st.error(f"Could not run ARStat: {exc}")
        st.stop()

stored_results = st.session_state.get("arstat_results")
stored_config = st.session_state.get("arstat_config")

if stored_results is None:
    st.info("Adjust the settings in the sidebar and click **Run ARStat**.")
    st.stop()

if stored_config != current_config:
    st.warning(
        "The input data or settings have changed since the displayed results were generated. "
        "Click **Run ARStat** again to update the results."
    )
    st.stop()

if stored_results is not None:
    render_arstat_results(stored_results)

st.caption("ARStat v1.0.0-rc1. Example datasets use illustrative hookworm assay data; downloads reuse stored results.")
