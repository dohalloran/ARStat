"""UI regression tests for the Streamlit app using Streamlit's AppTest harness.

These run the real app.py script headlessly, so they cover widget-state logic
(example locking, unit locking) that the core-function tests cannot reach.
"""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP_PATH = Path(__file__).resolve().parents[1] / "app.py"

EXAMPLES = {
    "Egg hatch example": ("Egg hatch", "µM"),
    "Larval development example": ("Larval development", "nM"),
    "Motility example": ("Motility", "nM"),
}


def _app():
    at = AppTest.from_file(str(APP_PATH), default_timeout=60)
    at.run()
    assert not at.exception
    return at


def test_example_selection_forces_and_locks_matching_assay_type():
    at = _app()
    # Visit examples in a scrambled order, including revisits, to catch stale widget state.
    order = ["Motility example", "Egg hatch example",
             "Larval development example", "Egg hatch example", "Motility example"]
    for label in order:
        at.selectbox(key="sample_label").select(label).run()
        assert not at.exception
        assay_box = at.selectbox(key="assay_type")
        expected_assay, expected_unit = EXAMPLES[label]
        assert assay_box.value == expected_assay, label
        assert assay_box.disabled, label
        unit_box = at.text_input(key="dose_unit_input")
        assert unit_box.value == expected_unit, label
        assert unit_box.disabled, label


def test_every_example_runs_and_fits_both_isolates():
    at = _app()
    for label in EXAMPLES:
        at.selectbox(key="sample_label").select(label).run()
        at.button(key="run_arstat_btn").click().run()
        assert not at.exception, label
        assert not at.error, [e.value for e in at.error]
        fitted = next(m for m in at.metric if m.label == "Groups fitted")
        assert fitted.value == "2/2", label


def test_example_provenance_is_reported_accurately():
    at = _app()
    for label, expected in [
        ("Egg hatch example", "simulated/illustrative egg-hatch"),
        ("Larval development example", "simulated/illustrative larval-development"),
        ("Motility example", "simulated motility data"),
    ]:
        at.selectbox(key="sample_label").select(label).run()
        messages = " ".join(s.value for s in at.success)
        assert expected in messages, (label, messages)
        assert "real experimental" not in messages
