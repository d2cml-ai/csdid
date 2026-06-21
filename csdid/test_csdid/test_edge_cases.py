"""
Tests for edge cases and boundary conditions.
Translated from R package 'did' tests/testthat/test-edge-cases.R (26 assertions).
"""

import os
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt

EST_METHODS = ["dr", "reg", "ipw"]


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def sim_data():
    """Load simulated dataset (equivalent to did::build_sim_dataset(did::reset.sim()))."""
    data_path = os.path.join(
        os.path.dirname(__file__), os.pardir, "data", "sim_data.csv"
    )
    return pd.read_csv(data_path)


# =============================================================================
# Helpers
# =============================================================================

def _make_single_group_data(seed=20260401):
    """Single treated group (G=4, 50 units) + never-treated (G=0, 150 units),
    5 periods, treatment effect of +1 starting at period 4."""
    np.random.seed(seed)
    n_ids, n_periods = 200, 5
    ids = np.repeat(np.arange(1, n_ids + 1), n_periods)
    periods = np.tile(np.arange(1, n_periods + 1), n_ids)
    groups = np.repeat(
        np.concatenate([np.full(50, 4), np.full(150, 0)]), n_periods
    )
    y = np.random.randn(n_ids * n_periods)

    df = pd.DataFrame({"id": ids, "period": periods, "G": groups, "Y": y})
    df.loc[(df["G"] == 4) & (df["period"] >= 4), "Y"] += 1
    return df


def _make_no_nevertreated_data(seed=20260401):
    """All units eventually treated: groups 3 (100 units) and 5 (100 units),
    6 periods, no never-treated (G=0)."""
    np.random.seed(seed)
    n_ids, n_periods = 200, 6
    ids = np.repeat(np.arange(1, n_ids + 1), n_periods)
    periods = np.tile(np.arange(1, n_periods + 1), n_ids)
    groups = np.repeat(
        np.concatenate([np.full(100, 3), np.full(100, 5)]), n_periods
    )
    y = np.random.randn(n_ids * n_periods)
    return pd.DataFrame({"id": ids, "period": periods, "G": groups, "Y": y})


# =============================================================================
# Tests — single treated group
# =============================================================================

@pytest.mark.parametrize("em", EST_METHODS)
def test_single_treated_group_produces_valid_att_gt(em):
    """R: single treated group produces valid att_gt"""
    data = _make_single_group_data()

    result = ATTgt(
        yname="Y", tname="period", idname="id", gname="G", data=data
    ).fit(est_method=em, bstrap=False)

    assert isinstance(result.MP, dict)
    att = np.asarray(result.MP["att"], dtype=float)
    assert np.any(~np.isnan(att))


@pytest.mark.parametrize("agg_type", ["simple", "dynamic", "group", "calendar"])
def test_single_treated_group_aggte_types(agg_type):
    """R: single treated group works with all 4 aggte types"""
    data = _make_single_group_data()

    mp = ATTgt(
        yname="Y", tname="period", idname="id", gname="G", data=data
    ).fit(bstrap=False)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        agg = mp.aggte(typec=agg_type)

    assert agg.atte is not None
    assert not np.isnan(agg.atte["overall_att"]), (
        f"overall_att is NaN for type={agg_type}"
    )


# =============================================================================
# Tests — two-period data
# =============================================================================

@pytest.mark.parametrize("em", EST_METHODS)
def test_two_period_data_with_universal_base_period(em):
    """R: two-period data with universal base period"""
    np.random.seed(20260401)
    n_ids = 200
    ids = np.repeat(np.arange(1, n_ids + 1), 2)
    periods = np.tile([1, 2], n_ids)
    groups = np.repeat(
        np.concatenate([np.full(50, 2), np.full(150, 0)]), 2
    )
    y = np.random.randn(n_ids * 2)

    df = pd.DataFrame({"id": ids, "period": periods, "G": groups, "Y": y})
    df.loc[(df["G"] == 2) & (df["period"] == 2), "Y"] += 1

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = ATTgt(
            yname="Y", tname="period", idname="id", gname="G", data=df
        ).fit(est_method=em, base_period="universal", bstrap=False)

    assert isinstance(result.MP, dict)
    att = np.asarray(result.MP["att"], dtype=float)
    assert np.any(~np.isnan(att))


# =============================================================================
# Tests — no never-treated group
# =============================================================================

def test_no_nevertreated_works_with_notyettreated():
    """R: data with no never-treated group works with notyettreated"""
    df = _make_no_nevertreated_data()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=df, control_group="notyettreated",
        ).fit(bstrap=False)

    assert isinstance(result.MP, dict)


def test_no_nevertreated_coerces_with_nevertreated():
    """R did v2.5.1: data with no never-treated group + control_group='nevertreated'
    warns and coerces the last cohort to never-treated (it does not error).
    csdid now matches this behavior.
    """
    df = _make_no_nevertreated_data()

    with pytest.warns(UserWarning, match="never-treated"):
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=df, control_group="nevertreated",
        )
    # Coercion creates a never-treated (0) group in the processed data.
    assert (obj.dp["data"]["G"] == 0).any()
    # And the run produces finite ATT estimates.
    res = obj.fit(bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert np.isfinite(att).any()


# =============================================================================
# Tests — groups treated in first period
# =============================================================================

def test_groups_treated_in_first_period_are_dropped(sim_data):
    """R: groups treated in first period are dropped (with a warning).
    csdid issues a UserWarning naming the dropped units (matches R's warning())."""
    data = sim_data.copy()
    first_per = int(data["period"].min())

    treated_groups = np.sort(data.loc[data["G"] > 0, "G"].unique())
    first_group = treated_groups[0]

    extra = data.loc[data["G"] == first_group].copy()
    extra["G"] = first_per
    extra["id"] = extra["id"] + data["id"].max()
    data = pd.concat([data, extra], ignore_index=True)

    with pytest.warns(UserWarning, match="already treated|[Dd]ropped"):
        ATTgt(
            yname="Y", tname="period", idname="id", gname="G", data=data
        ).fit(bstrap=False)


# =============================================================================
# Tests — non-consecutive periods / groups
# =============================================================================

@pytest.mark.parametrize("em", EST_METHODS)
def test_non_consecutive_time_periods(em):
    """R: non-consecutive time periods work"""
    np.random.seed(20260401)
    n_ids = 200
    time_periods = [2000, 2003, 2007, 2010]
    n_t = len(time_periods)
    ids = np.repeat(np.arange(1, n_ids + 1), n_t)
    periods = np.tile(time_periods, n_ids)
    groups = np.repeat(
        np.concatenate([np.full(50, 2007), np.full(150, 0)]), n_t
    )
    y = np.random.randn(n_ids * n_t)

    df = pd.DataFrame({"id": ids, "period": periods, "G": groups, "Y": y})
    df.loc[(df["G"] == 2007) & (df["period"] >= 2007), "Y"] += 1

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = ATTgt(
            yname="Y", tname="period", idname="id", gname="G", data=df
        ).fit(est_method=em, bstrap=False)

    assert isinstance(result.MP, dict)
    att = np.asarray(result.MP["att"], dtype=float)
    assert np.any(~np.isnan(att))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        agg = result.aggte(typec="dynamic")
    assert agg.atte is not None


@pytest.mark.parametrize("em", EST_METHODS)
def test_non_consecutive_group_values(em):
    """R: non-consecutive group values work"""
    np.random.seed(20260401)
    n_ids, n_periods = 200, 5
    ids = np.repeat(np.arange(1, n_ids + 1), n_periods)
    periods = np.tile(np.arange(1, n_periods + 1), n_ids)
    groups = np.repeat(
        np.concatenate([np.full(50, 3), np.full(50, 5), np.full(100, 0)]),
        n_periods,
    )
    y = np.random.randn(n_ids * n_periods)

    df = pd.DataFrame({"id": ids, "period": periods, "G": groups, "Y": y})

    result = ATTgt(
        yname="Y", tname="period", idname="id", gname="G", data=df
    ).fit(est_method=em, bstrap=False)

    assert isinstance(result.MP, dict)
    assert len(np.unique(result.MP["group"])) >= 2


# =============================================================================
# Tests — unbalanced panel
# =============================================================================

def test_allow_unbalanced_panel_with_balanced_data(sim_data):
    """R: allow_unbalanced_panel=TRUE with balanced data proceeds normally"""
    result = ATTgt(
        yname="Y", tname="period", idname="id", gname="G",
        data=sim_data, allow_unbalanced_panel=True,
    ).fit(bstrap=False)

    assert isinstance(result.MP, dict)
    att = np.asarray(result.MP["att"], dtype=float)
    assert np.any(~np.isnan(att))


def test_allow_unbalanced_panel_with_truly_unbalanced_data(sim_data):
    """R: allow_unbalanced_panel=TRUE with truly unbalanced data.
    csdid now routes truly unbalanced panels through the RC estimators (matches R)."""
    data = sim_data.copy()
    data = data.drop([0, 4, 9]).reset_index(drop=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, allow_unbalanced_panel=True,
        ).fit(bstrap=False)

    assert isinstance(result.MP, dict)
    att = np.asarray(result.MP["att"], dtype=float)
    assert np.any(~np.isnan(att))


# =============================================================================
# Tests — single post-treatment period
# =============================================================================

@pytest.mark.parametrize("em", EST_METHODS)
def test_single_post_treatment_period(em):
    """R: single post-treatment period produces valid results"""
    np.random.seed(20260401)
    n_ids, n_periods = 200, 3
    ids = np.repeat(np.arange(1, n_ids + 1), n_periods)
    periods = np.tile(np.arange(1, n_periods + 1), n_ids)
    groups = np.repeat(
        np.concatenate([np.full(50, 3), np.full(150, 0)]), n_periods
    )
    y = np.random.randn(n_ids * n_periods)

    df = pd.DataFrame({"id": ids, "period": periods, "G": groups, "Y": y})
    df.loc[(df["G"] == 3) & (df["period"] == 3), "Y"] += 1

    result = ATTgt(
        yname="Y", tname="period", idname="id", gname="G", data=df
    ).fit(est_method=em, bstrap=False)

    assert isinstance(result.MP, dict)
    groups_arr = np.asarray(result.MP["group"])
    t_arr = np.asarray(result.MP["t"])
    att_arr = np.asarray(result.MP["att"], dtype=float)
    post_atts = att_arr[groups_arr <= t_arr]
    assert np.any(~np.isnan(post_atts))


# =============================================================================
# Tests — RCS with unbalanced data (works unlike panel=True)
# =============================================================================

@pytest.mark.parametrize("em", EST_METHODS)
def test_unbalanced_data_rcs_works(sim_data, em):
    """Unbalanced data works with panel=False (repeated cross-sections)."""
    data = sim_data.copy()
    data = data.drop([0, 4, 9]).reset_index(drop=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, panel=False,
        ).fit(est_method=em, bstrap=False)

    assert isinstance(result.MP, dict)
    att = np.asarray(result.MP["att"], dtype=float)
    assert np.any(~np.isnan(att))
