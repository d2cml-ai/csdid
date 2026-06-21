"""
Tests for user-reported bug fixes.
Translated from R package 'did' tests/testthat/test-user_bug_fixes.R
"""
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from helpers import build_sim_data


DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")


@pytest.fixture(scope="module")
def mpdta():
    return pd.read_csv(os.path.join(DATA_DIR, "mpdta.csv"))


# ─────────────────────────────────────────────────────────────
# Bug: column named t1 crashes code
# ─────────────────────────────────────────────────────────────

def test_column_named_t1_does_not_crash(mpdta):
    """R: having column named t1 causes code to crash."""
    data = mpdta.copy()

    # First run without t1 column
    out1 = ATTgt(
        yname="lemp", gname="first.treat", idname="countyreal",
        tname="year", data=data, control_group="notyettreated",
    ).fit(est_method="reg", bstrap=False)
    assert out1.MP is not None

    # Add t1 column and run again
    data["t1"] = 1
    out2 = ATTgt(
        yname="lemp", gname="first.treat", idname="countyreal",
        tname="year", data=data, control_group="notyettreated",
    ).fit(est_method="reg", bstrap=False)
    assert out2.MP is not None


# ─────────────────────────────────────────────────────────────
# Bug: missing covariates should warn but proceed
# ─────────────────────────────────────────────────────────────

def test_missing_covariates(mpdta, capsys):
    """R: missing covariates — should warn about missing data but run.
    csdid now routes the resulting unbalanced panel through RC estimators (matches R)."""
    data = mpdta.copy()
    data.loc[data.index[0], "lpop"] = np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = ATTgt(
            yname="lemp", gname="first.treat", idname="countyreal",
            tname="year", xformla="lemp~lpop", data=data,
            control_group="notyettreated",
        ).fit(est_method="reg", bstrap=False)

    assert out.MP is not None


# ─────────────────────────────────────────────────────────────
# Bug: fewer time periods than groups
# ─────────────────────────────────────────────────────────────

def test_fewer_periods_than_groups():
    """R: fewer time periods than groups — from GitHub issue #56."""
    data = build_sim_data(
        n=1000, time_periods=6, te=0, te_e=[1, 2, 3, 4, 5, 6], seed=56
    )
    # Drop periods 2 and 5 — now 4 periods but up to 6 group values
    data = data[~data["period"].isin([2, 5])].reset_index(drop=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method="dr", bstrap=False)

    groups = np.asarray(res.MP["group"])
    t_arr = np.asarray(res.MP["t"])
    att = np.asarray(res.MP["att"], dtype=float)

    # Group 2 at t=3 should have ATT ≈ 2 (exposure=1→te_e[1]=2)
    idx = np.where((groups == 2) & (t_arr == 3))[0]
    if len(idx) > 0:
        assert abs(att[idx[0]] - 2) < 1.0

    # Dynamic aggregation should work
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dyn = res.aggte(typec="dynamic")
    assert dyn.atte is not None
    dyn_att = np.asarray(dyn.atte["att_egt"], dtype=float)
    assert not np.all(np.isnan(dyn_att))

    # Group and calendar aggregations
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grp = res.aggte(typec="group")
    assert grp.atte is not None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal = res.aggte(typec="calendar")
    cal_att = np.asarray(cal.atte["att_egt"], dtype=float)
    assert not np.isnan(cal_att[0])


# ─────────────────────────────────────────────────────────────
# Bug: 0 pre-treatment estimates when outcomes are 0
# ─────────────────────────────────────────────────────────────

def test_zero_pretreatment_outcomes():
    """R: 0 pre-treatment estimates when outcomes are 0 — issue #126."""
    data = build_sim_data(n=1000, time_periods=10, te=1.0, seed=126)
    data = data[data["G"] > 0].copy()       # drop never-treated
    data = data[data["G"] > 6].copy()        # keep only late groups
    data = data[data["period"] > 5].copy()   # keep only later periods
    data.loc[data["period"] < data["G"], "Y"] = 0  # zero pre-treatment

    for bp in ["universal", "varying"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=data, control_group="notyettreated",
            ).fit(base_period=bp, bstrap=False)

        groups = np.asarray(res.MP["group"])
        t_arr = np.asarray(res.MP["t"])
        att = np.asarray(res.MP["att"], dtype=float)

        # Pre-treatment ATT for G=9 at t=7 should be 0
        idx = np.where((groups == 9) & (t_arr == 7))[0]
        if len(idx) > 0:
            assert abs(att[idx[0]]) < 0.01, f"base_period={bp}"


# ─────────────────────────────────────────────────────────────
# Bug: variables not in dataset should error
# ─────────────────────────────────────────────────────────────

def test_variables_not_in_dataset():
    """R: referencing variables not in dataset should error.

    csdid now raises a ValueError for unknown formula variables (matches R),
    instead of silently falling back to an intercept-only model.
    """
    data = build_sim_data(n=500, time_periods=3, seed=99)
    with pytest.raises(Exception):
        ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X2",  # X2 not in data
            control_group="notyettreated",
        ).fit(est_method="dr", bstrap=False)


# ─────────────────────────────────────────────────────────────
# Bug: groups treated after max(t) within anticipation window
# ─────────────────────────────────────────────────────────────

def test_anticipation_window_coercion():
    """R: groups treated after max(t) but within anticipation window
    should NOT be coerced to never-treated.

    csdid now coerces groups treated after max(t)+anticipation to never-treated
    and respects the anticipation window (matches R did v2.5.1).
    """
    np.random.seed(20250228)
    n = 600
    ids = np.repeat(np.arange(1, n + 1), 5)
    times = np.tile(np.arange(1, 6), n)
    group_vals = np.concatenate([np.full(200, 0), np.full(200, 4), np.full(200, 6)])
    G = np.repeat(group_vals, 5)
    Y = np.random.randn(len(ids))
    dt = pd.DataFrame({"id": ids, "time": times, "group": G, "y": Y})

    # anticipation=0: group 6 is beyond max(t)=5, coerced to never-treated
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res0 = ATTgt(
            yname="y", gname="group", idname="id", tname="time",
            data=dt, anticipation=0,
            control_group="nevertreated",
        ).fit(bstrap=False)

    # Only group 4 should appear (R coerces group 6 to never-treated)
    assert all(np.asarray(res0.MP["group"]) == 4)

    # anticipation=2: group 6 anticipates at period 4, should remain treated
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res2 = ATTgt(
            yname="y", gname="group", idname="id", tname="time",
            data=dt, anticipation=2,
            control_group="nevertreated",
        ).fit(bstrap=False)

    assert 6 in np.asarray(res2.MP["group"])
    assert len(np.unique(res2.MP["group"])) > len(np.unique(res0.MP["group"]))


# ─────────────────────────────────────────────────────────────
# Cross est_method: key bug fixes should pass across all methods
# ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("em", ["dr", "reg", "ipw"])
def test_fewer_periods_than_groups_all_methods(em):
    """Fewer time periods than groups works across all methods."""
    data = build_sim_data(
        n=1000, time_periods=6, te=0, te_e=[1, 2, 3, 4, 5, 6], seed=56
    )
    data = data[~data["period"].isin([2, 5])].reset_index(drop=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method=em, bstrap=False)

    att = np.asarray(res.MP["att"], dtype=float)
    assert np.any(~np.isnan(att))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic")
    assert not np.isnan(res.atte["overall_att"])
