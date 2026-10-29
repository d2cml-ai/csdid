"""
Tests for tidy/summary representations of results.
Translated from R package 'did' tests/testthat/test-tidy.R (9 test_that blocks).

R tests check broom::tidy() and nobs() S3 methods. Python equivalents:
- tidy.MP → summ_attgt() returns self with .summary2 DataFrame
- MP result dict has expected keys (group, att, t, etc.)
- atte dict has expected keys (overall_att, att_egt, se_egt, etc.)
- nobs → DIDparams['n'] gives number of units
"""
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt


DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")


@pytest.fixture(scope="module")
def mp_tidy():
    """Fitted ATTgt object on mpdta for tidy tests."""
    data = pd.read_csv(os.path.join(DATA_DIR, "mpdta.csv"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = ATTgt(
            yname="lemp", gname="first.treat", idname="countyreal",
            tname="year", data=data,
        ).fit(est_method="reg", bstrap=False)
    return obj


@pytest.fixture(scope="module")
def aggte_results(mp_tidy):
    """Dict of aggte results keyed by type."""
    import copy
    results = {}
    for typec in ["simple", "dynamic", "group", "calendar"]:
        obj = copy.copy(mp_tidy)
        obj.MP = dict(mp_tidy.MP)
        obj.did_object = dict(mp_tidy.did_object)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obj.aggte(typec=typec, bstrap=False)
        results[typec] = obj.atte
    return results


# ─────────────────────────────────────────────────────────────
# summ_attgt (equivalent to tidy.MP)
# ─────────────────────────────────────────────────────────────

def test_summ_attgt_returns_dataframe(mp_tidy):
    """summ_attgt() returns a DataFrame with expected columns."""
    mp_tidy.summ_attgt()
    df = mp_tidy.summary2
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert "Group" in df.columns
    assert "Time" in df.columns
    assert "ATT(g, t)" in df.columns


def test_mp_result_has_expected_keys(mp_tidy):
    """MP dict has group, att, t, DIDparams, inffunc, n."""
    expected = {"group", "att", "t", "DIDparams", "inffunc", "n"}
    assert expected.issubset(set(mp_tidy.MP.keys()))


def test_att_values_are_finite(mp_tidy):
    """ATT(g,t) estimates are all finite."""
    att = np.asarray(mp_tidy.MP["att"], dtype=float)
    assert np.all(np.isfinite(att))


def test_se_values_are_positive(mp_tidy):
    """SE values are positive where ATT is finite."""
    se = np.asarray(mp_tidy.results["se"], dtype=float)
    att = np.asarray(mp_tidy.results["att"], dtype=float)
    finite = np.isfinite(att) & np.isfinite(se)
    assert np.all(se[finite] > 0)


# ─────────────────────────────────────────────────────────────
# aggte results (equivalent to tidy.AGGTEobj)
# ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("typec", ["simple", "dynamic", "group", "calendar"])
def test_aggte_has_expected_keys(aggte_results, typec):
    """aggte atte dict has expected keys."""
    atte = aggte_results[typec]
    assert "overall_att" in atte
    assert "overall_se" in atte
    assert "DIDparams" in atte


def test_aggte_dynamic_has_event_times(aggte_results):
    """Dynamic aggte has egt, att_egt, se_egt."""
    atte = aggte_results["dynamic"]
    assert atte["egt"] is not None
    assert atte["att_egt"] is not None
    assert atte["se_egt"] is not None


def test_aggte_group_has_group_ids(aggte_results):
    """Group aggte egt values are group identifiers."""
    atte = aggte_results["group"]
    assert atte["egt"] is not None
    egt = np.asarray(atte["egt"])
    assert len(egt) > 0


def test_aggte_calendar_has_periods(aggte_results):
    """Calendar aggte egt values are time periods."""
    atte = aggte_results["calendar"]
    assert atte["egt"] is not None
    egt = np.asarray(atte["egt"])
    assert len(egt) > 0


# ─────────────────────────────────────────────────────────────
# nobs (equivalent to nobs.MP / nobs.AGGTEobj)
# ─────────────────────────────────────────────────────────────

def test_nobs_from_didparams(mp_tidy):
    """DIDparams['n'] gives number of unique units."""
    n = mp_tidy.MP["DIDparams"]["n"]
    assert isinstance(n, (int, np.integer))
    assert n > 0


def test_nobs_matches_data(mp_tidy):
    """DIDparams n matches unique ids in the dataset."""
    data = mp_tidy.MP["DIDparams"]["data"]
    idname = mp_tidy.MP["DIDparams"]["idname"]
    expected_n = data[idname].nunique()
    assert mp_tidy.MP["DIDparams"]["n"] == expected_n


# ─────────────────────────────────────────────────────────────
# Cross est_method: verify tidy output structure is consistent
# ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("em", ["dr", "reg", "ipw"])
def test_summ_attgt_across_methods(em):
    """summ_attgt() produces valid DataFrame across all est_methods."""
    from helpers import build_sim_data
    data = build_sim_data(n=1000, time_periods=4, te=1.0, seed=20260401)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method=em, bstrap=False)
    obj.summ_attgt()
    df = obj.summary2
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert "ATT(g, t)" in df.columns
