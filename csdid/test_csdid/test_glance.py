"""
Tests for glance-style summary information.
Translated from R package 'did' tests/testthat/test-glance.R (11 test_that blocks).

R tests check broom::glance() S3 methods returning one-row DataFrames with
nobs, ngroup, ntime, control.group, est.method. Python equivalent: verify
that DIDparams contains these values and they are reasonable.

Tests requiring faster_mode are skipped — not in Python csdid.
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
def mp_glance():
    """Fitted ATTgt on simulated data for glance tests."""
    data = build_sim_data(n=1000, time_periods=4, te=1.0, seed=20260401)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method="dr", bstrap=False)
    return obj


@pytest.fixture(scope="module")
def agg_results(mp_glance):
    """Dict of aggte results for all 4 types."""
    import copy
    results = {}
    for typec in ["simple", "dynamic", "group", "calendar"]:
        obj = copy.copy(mp_glance)
        obj.MP = dict(mp_glance.MP)
        obj.did_object = dict(mp_glance.did_object)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obj.aggte(typec=typec, bstrap=False)
        results[typec] = obj.atte
    return results


# ─────────────────────────────────────────────────────────────
# glance.MP equivalents
# ─────────────────────────────────────────────────────────────

def test_didparams_has_expected_keys(mp_glance):
    """DIDparams has nobs, ngroup, ntime, control_group, est_method."""
    dp = mp_glance.MP["DIDparams"]
    assert "n" in dp          # nobs
    assert "nG" in dp         # ngroup
    assert "nT" in dp         # ntime
    assert "control_group" in dp
    assert "est_method" in dp


def test_didparams_values_are_reasonable(mp_glance):
    """DIDparams values are positive and correct types."""
    dp = mp_glance.MP["DIDparams"]
    assert dp["n"] > 0
    assert dp["nG"] > 0
    assert dp["nT"] > 0
    assert dp["control_group"] in ["nevertreated", "notyettreated"]
    assert dp["est_method"] in ["dr", "reg", "ipw"]


def test_nobs_matches_unique_ids(mp_glance):
    """DIDparams n matches number of unique units in data."""
    dp = mp_glance.MP["DIDparams"]
    data = dp["data"]
    expected = data[dp["idname"]].nunique()
    assert dp["n"] == expected


def test_ngroup_matches_treatment_groups(mp_glance):
    """DIDparams nG matches number of treatment groups."""
    dp = mp_glance.MP["DIDparams"]
    assert dp["nG"] == len(dp["glist"])


def test_ntime_matches_time_periods(mp_glance):
    """DIDparams nT matches number of time periods."""
    dp = mp_glance.MP["DIDparams"]
    assert dp["nT"] == len(dp["tlist"])


# ─────────────────────────────────────────────────────────────
# glance.AGGTEobj equivalents
# ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("typec", ["simple", "dynamic", "group", "calendar"])
def test_aggte_has_didparams(agg_results, typec):
    """aggte result has DIDparams with expected keys."""
    dp = agg_results[typec]["DIDparams"]
    assert dp is not None
    assert "n" in dp
    assert "nG" in dp
    assert "nT" in dp


@pytest.mark.parametrize("typec", ["simple", "dynamic", "group", "calendar"])
def test_aggte_didparams_not_null(agg_results, typec):
    """aggte DIDparams values are not None or NaN."""
    dp = agg_results[typec]["DIDparams"]
    assert dp["n"] is not None and not np.isnan(dp["n"])
    assert dp["nG"] is not None and not np.isnan(dp["nG"])
    assert dp["nT"] is not None and not np.isnan(dp["nT"])


def test_aggte_and_mp_nobs_agree(mp_glance, agg_results):
    """nobs from MP and all aggte types should agree."""
    n_mp = mp_glance.MP["DIDparams"]["n"]
    for typec in ["simple", "dynamic", "group", "calendar"]:
        n_agg = agg_results[typec]["DIDparams"]["n"]
        assert n_mp == n_agg, f"nobs mismatch for typec={typec}"


# ─────────────────────────────────────────────────────────────
# Skipped: faster_mode comparison
# ─────────────────────────────────────────────────────────────

def test_glance_faster_mode():
    """glance-style summary works with faster_mode."""
    data = build_sim_data(n=1000, time_periods=4, te=1.0, seed=20260401)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = ATTgt(yname="Y", tname="period", idname="id", gname="G",
                    data=data, xformla="Y~X", faster_mode=True).fit(est_method="dr", bstrap=False)
    assert obj.dp["n"] == data["id"].nunique()
    assert np.isfinite(np.asarray(obj.MP["att"], dtype=float)).any()


def test_glance_faster_mode_agreement():
    """faster_mode and standard agree on glance values and point estimates."""
    data = build_sim_data(n=1000, time_periods=4, te=1.0, seed=20260401)
    common = dict(yname="Y", tname="period", idname="id", gname="G", data=data, xformla="Y~X")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        std = ATTgt(**common, faster_mode=False).fit(est_method="dr", bstrap=False)
        fast = ATTgt(**common, faster_mode=True).fit(est_method="dr", bstrap=False)
    assert std.dp["n"] == fast.dp["n"]
    assert std.dp["nG"] == fast.dp["nG"] and std.dp["nT"] == fast.dp["nT"]
    np.testing.assert_allclose(
        np.asarray(std.MP["att"], dtype=float),
        np.asarray(fast.MP["att"], dtype=float), atol=1e-9)


# ─────────────────────────────────────────────────────────────
# Cross est_method: verify DIDparams consistency
# ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("em", ["dr", "reg", "ipw"])
def test_didparams_consistent_across_methods(em):
    """DIDparams has correct est_method across all methods."""
    data = build_sim_data(n=1000, time_periods=4, te=1.0, seed=20260401)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method=em, bstrap=False)
    dp = obj.MP["DIDparams"]
    assert dp["est_method"] == em
    assert dp["n"] > 0
    assert dp["nG"] > 0
    assert dp["nT"] > 0
