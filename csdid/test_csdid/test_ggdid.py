"""
Tests for plotting functions (ggdid equivalents).
Translated from R package 'did' tests/testthat/test-ggdid.R (11 test_that blocks).

R tests check ggdid() returns ggplot objects. Python equivalents:
- plot_attgt() returns matplotlib Figure
- plot_aggte() returns matplotlib Axes
"""
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from helpers import build_sim_data

# Use non-interactive backend so plots don't pop up during tests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


@pytest.fixture(scope="module")
def mp_plot():
    """Fitted ATTgt object for plotting tests."""
    data = build_sim_data(n=1000, time_periods=4, te=1.0, seed=20260401)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method="dr", bstrap=False)
    return obj


@pytest.fixture(scope="module")
def agg_dynamic(mp_plot):
    """Dynamic aggte result."""
    import copy
    obj = copy.copy(mp_plot)
    obj.MP = dict(mp_plot.MP)
    obj.did_object = dict(mp_plot.did_object)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec="dynamic", bstrap=False)
    return obj


@pytest.fixture(scope="module")
def agg_group(mp_plot):
    """Group aggte result."""
    import copy
    obj = copy.copy(mp_plot)
    obj.MP = dict(mp_plot.MP)
    obj.did_object = dict(mp_plot.did_object)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec="group", bstrap=False)
    return obj


@pytest.fixture(scope="module")
def agg_calendar(mp_plot):
    """Calendar aggte result."""
    import copy
    obj = copy.copy(mp_plot)
    obj.MP = dict(mp_plot.MP)
    obj.did_object = dict(mp_plot.did_object)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec="calendar", bstrap=False)
    return obj


# ─────────────────────────────────────────────────────────────
# plot_attgt (equivalent to ggdid.MP)
# ─────────────────────────────────────────────────────────────

def test_plot_attgt_returns_figure(mp_plot):
    """plot_attgt() returns a matplotlib Figure."""
    fig = mp_plot.plot_attgt()
    assert isinstance(fig, Figure)
    plt.close("all")


def test_plot_attgt_with_single_group(mp_plot):
    """plot_attgt() works with a single group specified."""
    groups = np.unique(mp_plot.did_object["group"]).astype(int)
    fig = mp_plot.plot_attgt(group=[groups[0]])
    assert isinstance(fig, Figure)
    plt.close("all")


def test_plot_attgt_with_custom_labels(mp_plot):
    """plot_attgt() accepts custom labels."""
    fig = mp_plot.plot_attgt(xlab="Time", ylab="Effect", title="Test")
    assert isinstance(fig, Figure)
    plt.close("all")


def test_plot_attgt_invalid_group_raises(mp_plot):
    """plot_attgt() raises for non-existent group values."""
    with pytest.raises((ValueError, IndexError)):
        mp_plot.plot_attgt(group=[9999])
    plt.close("all")


# ─────────────────────────────────────────────────────────────
# plot_aggte (equivalent to ggdid.AGGTEobj)
# ─────────────────────────────────────────────────────────────

def test_plot_aggte_dynamic(agg_dynamic):
    """plot_aggte() works for dynamic type."""
    ax = agg_dynamic.plot_aggte()
    assert isinstance(ax, Axes)
    plt.close("all")


def test_plot_aggte_group(agg_group):
    """plot_aggte() works for group type."""
    ax = agg_group.plot_aggte()
    assert isinstance(ax, Axes)
    plt.close("all")


def test_plot_aggte_calendar(agg_calendar):
    """plot_aggte() works for calendar type."""
    ax = agg_calendar.plot_aggte()
    assert isinstance(ax, Axes)
    plt.close("all")


def test_plot_aggte_custom_labels(agg_dynamic):
    """plot_aggte() accepts custom labels and theme settings."""
    ax = agg_dynamic.plot_aggte(
        xlab="Event Time", ylab="ATT", title="Event Study",
        theming=True, legend=True,
    )
    assert isinstance(ax, Axes)
    plt.close("all")


def test_plot_aggte_no_theming(agg_dynamic):
    """plot_aggte() works with theming=False."""
    ax = agg_dynamic.plot_aggte(theming=False)
    assert isinstance(ax, Axes)
    plt.close("all")


def test_plot_aggte_no_ref_line(agg_dynamic):
    """plot_aggte() works with ref_line=None."""
    ax = agg_dynamic.plot_aggte(ref_line=None)
    assert isinstance(ax, Axes)
    plt.close("all")


# ─────────────────────────────────────────────────────────────
# plot_aggte for all types × theming settings
# ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("typec", ["dynamic", "group", "calendar"])
@pytest.mark.parametrize("theming", [True, False])
def test_plot_aggte_all_types_theming(mp_plot, typec, theming):
    """plot_aggte() works for all types with theming on/off."""
    import copy
    obj = copy.copy(mp_plot)
    obj.MP = dict(mp_plot.MP)
    obj.did_object = dict(mp_plot.did_object)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec=typec, bstrap=False)
    ax = obj.plot_aggte(theming=theming)
    assert isinstance(ax, Axes)
    plt.close("all")
