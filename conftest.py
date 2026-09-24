"""Pytest configuration shared by all ARStat tests.

Puts the repository root on ``sys.path`` so ``import arstat_core`` works with a
plain ``pytest`` call (no ``PYTHONPATH=.`` needed), and runs every test from the
repository root so relative paths such as ``sample_data/...`` resolve even when
pytest is launched from another directory.
"""
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture(autouse=True)
def _run_from_repo_root(monkeypatch):
    monkeypatch.chdir(ROOT)
