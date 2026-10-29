"""Shared fixtures and helpers for csdid tests."""
import os

import pandas as pd
import pytest

DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: bootstrap-heavy tests (R skip_on_cran); deselect with -m 'not slow'",
    )


@pytest.fixture(scope="module")
def sim_data():
    """Default simulated dataset (R equivalent: build_sim_dataset(reset.sim()))."""
    return pd.read_csv(os.path.join(DATA_DIR, "sim_data.csv"))


@pytest.fixture(scope="module")
def mpdta():
    """Minimum wage panel data (mpdta)."""
    return pd.read_csv(os.path.join(DATA_DIR, "mpdta.csv"))
