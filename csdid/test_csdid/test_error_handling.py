"""
Tests for error handling and validation.
Translated from R package 'did' tests/testthat/test-error-handling.R

R tests that require faster_mode, fix_weights, or R-specific validation
(e.g. dreamerr checks) are skipped — those features are not in csdid Python.
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
def sim_data():
    return pd.read_csv(os.path.join(DATA_DIR, "sim_data.csv"))


# ─────────────────────────────────────────────────────────────
# att_gt() validation errors
# ─────────────────────────────────────────────────────────────

def test_att_gt_errors_on_missing_column(sim_data):
    """R: att_gt errors on missing column name."""
    with pytest.raises(Exception):
        ATTgt(
            yname="nonexistent", tname="period", idname="id",
            gname="G", data=sim_data,
        )


def test_att_gt_errors_on_bad_idname(sim_data):
    """R: att_gt errors on incorrectly specified id."""
    with pytest.raises(Exception):
        ATTgt(
            yname="Y", tname="period", idname="brant",
            gname="G", data=sim_data,
        )


def test_att_gt_coerces_no_nevertreated(sim_data):
    """R did v2.5.1: no never-treated group + nevertreated control warns and
    coerces the last cohort to never-treated (it no longer raises)."""
    data = sim_data[sim_data["G"] > 0].copy()
    with pytest.warns(UserWarning, match="never-treated"):
        ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, control_group="nevertreated",
        )


def test_att_gt_warns_on_small_groups(sim_data):
    """R: att_gt warns on small groups."""
    data = sim_data.copy()
    treated_ids = data.loc[data["G"] > 0, "id"].unique()
    control_ids = data.loc[data["G"] == 0, "id"].unique()[:50]
    keep_ids = np.concatenate([treated_ids[:2], control_ids])
    small = data[data["id"].isin(keep_ids)].copy()

    with pytest.warns(UserWarning, match="small groups"):
        ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=small,
        ).fit(bstrap=False)


def test_att_gt_drops_first_period_units(sim_data):
    """R: att_gt warns when units are treated in the first period."""
    data = sim_data.copy()
    first_per = int(data["period"].min())
    treated_g = np.sort(data.loc[data["G"] > 0, "G"].unique())
    extra = data.loc[data["G"] == treated_g[0]].copy()
    extra["G"] = first_per
    extra["id"] = extra["id"] + data["id"].max()
    data = pd.concat([data, extra], ignore_index=True)

    with pytest.warns(UserWarning, match="already treated|[Dd]ropped"):
        ATTgt(
            yname="Y", tname="period", idname="id", gname="G", data=data,
        ).fit(bstrap=False)


def test_att_gt_handles_na_in_outcome(sim_data, capsys):
    """R: att_gt warns/prints when dropping rows with missing data.
    csdid drops NA rows (printing a message) and routes the resulting
    unbalanced panel through RC estimators (matches R)."""
    data = sim_data.copy()
    data.loc[data.index[:5], "Y"] = np.nan

    ATTgt(
        yname="Y", tname="period", idname="id", gname="G", data=data,
    ).fit(bstrap=False)

    out = capsys.readouterr().out.lower()
    assert "dropped" in out or "missing" in out


# ─────────────────────────────────────────────────────────────
# aggte() validation — tests that apply to Python
# ─────────────────────────────────────────────────────────────

def test_aggte_handles_na_with_na_rm(sim_data):
    """R: aggte with na.rm=TRUE drops NA ATTs and proceeds."""
    obj = ATTgt(
        yname="Y", tname="period", idname="id", gname="G", data=sim_data,
    ).fit(bstrap=False)

    # Inject NA
    obj.MP["att"] = list(obj.MP["att"])
    obj.MP["att"][0] = np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec="dynamic", na_rm=True)

    assert obj.atte is not None
    assert not np.isnan(obj.atte["overall_att"])


# ─────────────────────────────────────────────────────────────
# Skipped tests (R-only features)
# ─────────────────────────────────────────────────────────────

def test_fix_weights_validation(sim_data):
    """fix_weights validation matches R: bad value errors; base/first need a panel."""
    with pytest.raises(ValueError, match="fix_weights"):
        ATTgt(yname="Y", tname="period", idname="id", gname="G",
              data=sim_data, fix_weights="nope")
    with pytest.raises(ValueError, match="repeated cross sections"):
        ATTgt(yname="Y", tname="period", idname="id", gname="G",
              data=sim_data, panel=False, fix_weights="first_period")


def test_faster_mode_missing_column():
    """faster_mode still validates inputs: a missing column raises (matches standard)."""
    data = build_sim_data(n=500, time_periods=3, seed=99)
    with pytest.raises(Exception):
        ATTgt(yname="Y", tname="period", idname="id", gname="nonexistent_col",
              data=data, faster_mode=True)


def test_invalid_est_method(sim_data):
    """An unknown est_method string raises a clear error (matches R)."""
    with pytest.raises(ValueError, match="est_method"):
        ATTgt(
            yname="Y", tname="period", idname="id", gname="G", data=sim_data,
        ).fit(est_method="bogus", bstrap=False)


def test_treatment_reversals(sim_data):
    """Time-varying gname (reversible treatment) raises ValueError (matches R)."""
    data = sim_data.copy()
    uid = data["id"].iloc[0]
    idx = data.index[data["id"] == uid]
    # Flip this unit's cohort across periods -> treatment is no longer irreversible.
    data.loc[idx[: len(idx) // 2], "G"] = 0
    data.loc[idx[len(idx) // 2:], "G"] = 2
    with pytest.raises(ValueError, match="irreversible"):
        ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data)
