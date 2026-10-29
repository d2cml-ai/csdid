"""
Comprehensive tests for aggte() aggregation types.
Translated from R package 'did' tests/testthat/test-aggte-comprehensive.R (45 assertions).
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from helpers import build_sim_data


# ─────────────────────────────────────────────────────────────
# Shared fixture: compute att_gt once, reuse across tests
# ─────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def mp_agg():
    """ATTgt object fitted on simulated data with constant te=1 (R reset.sim())."""
    data = build_sim_data(n=1000, time_periods=4, te=1.0, seed=9142024)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method="dr", bstrap=False)
    return obj


def _aggte(mp, typec, **kwargs):
    """Run aggte and return atte dict, suppressing warnings."""
    import copy
    obj = copy.copy(mp)
    obj.MP = dict(mp.MP)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec=typec, **kwargs)
    return obj.atte


# ═════════════════════════════════════════════════════════════
# type = "simple"
# ═════════════════════════════════════════════════════════════

def test_simple_valid_overall_att(mp_agg):
    """aggte simple returns valid overall ATT."""
    agg = _aggte(mp_agg, "simple")
    assert agg is not None
    assert not np.isnan(agg["overall_att"])
    assert abs(agg["overall_att"] - 1) < 0.5


def test_simple_valid_se(mp_agg):
    """aggte simple returns valid SE."""
    agg = _aggte(mp_agg, "simple")
    assert agg["overall_se"] > 0
    assert not np.isnan(agg["overall_se"])


def test_simple_no_egt(mp_agg):
    """aggte simple has no egt component."""
    agg = _aggte(mp_agg, "simple")
    assert agg.get("egt") is None


# ═════════════════════════════════════════════════════════════
# type = "dynamic"
# ═════════════════════════════════════════════════════════════

def test_dynamic_returns_event_times(mp_agg):
    """aggte dynamic returns event-time specific ATTs."""
    agg = _aggte(mp_agg, "dynamic")
    egt = np.asarray(agg["egt"])
    assert len(egt) > 0
    assert np.any(egt < 0)     # pre-treatment
    assert np.any(egt >= 0)    # post-treatment


def test_dynamic_event_times_sorted(mp_agg):
    """aggte dynamic event times are sorted."""
    agg = _aggte(mp_agg, "dynamic")
    egt = np.asarray(agg["egt"], dtype=float)
    assert np.all(np.diff(egt) >= 0)


def test_dynamic_overall_att(mp_agg):
    """aggte dynamic overall.att averages post-treatment event times."""
    agg = _aggte(mp_agg, "dynamic")
    assert not np.isnan(agg["overall_att"])
    assert abs(agg["overall_att"] - 1) < 0.5


def test_dynamic_min_e_filters(mp_agg):
    """aggte dynamic min_e filters event times."""
    full = _aggte(mp_agg, "dynamic")
    filt = _aggte(mp_agg, "dynamic", min_e=-1)
    egt_filt = np.asarray(filt["egt"], dtype=float)
    assert np.min(egt_filt) >= -1
    assert len(egt_filt) <= len(np.asarray(full["egt"]))


def test_dynamic_max_e_filters(mp_agg):
    """aggte dynamic max_e filters event times."""
    full = _aggte(mp_agg, "dynamic")
    filt = _aggte(mp_agg, "dynamic", max_e=1)
    egt_filt = np.asarray(filt["egt"], dtype=float)
    assert np.max(egt_filt) <= 1
    assert len(egt_filt) <= len(np.asarray(full["egt"]))


def test_dynamic_min_max_together(mp_agg):
    """aggte dynamic min_e and max_e together."""
    agg = _aggte(mp_agg, "dynamic", min_e=-1, max_e=1)
    egt = np.asarray(agg["egt"], dtype=float)
    assert np.min(egt) >= -1
    assert np.max(egt) <= 1


def test_dynamic_balance_e(mp_agg):
    """aggte dynamic balance_e filters groups."""
    unbal = _aggte(mp_agg, "dynamic")
    bal = _aggte(mp_agg, "dynamic", balance_e=1)
    assert len(np.asarray(bal["egt"])) <= len(np.asarray(unbal["egt"]))


def test_dynamic_positive_se(mp_agg):
    """aggte dynamic SEs are positive where ATT is not NA."""
    agg = _aggte(mp_agg, "dynamic")
    att = np.asarray(agg["att_egt"], dtype=float)
    se = np.asarray(agg["se_egt"]).flatten()
    non_na = ~np.isnan(att)
    assert np.all(se[non_na] > 0)


# ═════════════════════════════════════════════════════════════
# type = "group"
# ═════════════════════════════════════════════════════════════

def test_group_returns_per_group(mp_agg):
    """aggte group returns per-group ATTs matching MP groups."""
    agg = _aggte(mp_agg, "group")
    egt = np.asarray(agg["egt"])
    mp_groups = np.unique(mp_agg.MP["group"])
    assert all(g in mp_groups for g in egt)


def test_group_overall_att(mp_agg):
    """aggte group overall.att is reasonable."""
    agg = _aggte(mp_agg, "group")
    assert not np.isnan(agg["overall_att"])
    assert abs(agg["overall_att"] - 1) < 0.5


def test_group_positive_se(mp_agg):
    """aggte group SEs are positive for each group."""
    agg = _aggte(mp_agg, "group")
    att = np.asarray(agg["att_egt"], dtype=float)
    se = np.asarray(agg["se_egt"]).flatten()
    non_na = ~np.isnan(att)
    assert np.all(se[non_na] > 0)


# ═════════════════════════════════════════════════════════════
# type = "calendar"
# ═════════════════════════════════════════════════════════════

def test_calendar_returns_periods(mp_agg):
    """aggte calendar returns per-period ATTs."""
    agg = _aggte(mp_agg, "calendar")
    assert agg["egt"] is not None
    assert len(np.asarray(agg["egt"])) > 0


def test_calendar_overall_att(mp_agg):
    """aggte calendar overall.att is reasonable."""
    agg = _aggte(mp_agg, "calendar")
    assert not np.isnan(agg["overall_att"])
    assert abs(agg["overall_att"] - 1) < 0.5


def test_calendar_positive_se(mp_agg):
    """aggte calendar SEs are positive for each period."""
    agg = _aggte(mp_agg, "calendar")
    att = np.asarray(agg["att_egt"], dtype=float)
    se = np.asarray(agg["se_egt"]).flatten()
    non_na = ~np.isnan(att)
    assert np.all(se[non_na] > 0)


def test_calendar_post_treatment_only(mp_agg):
    """aggte calendar only includes post-treatment periods."""
    agg = _aggte(mp_agg, "calendar")
    egt = np.asarray(agg["egt"], dtype=float)
    min_group = np.min(mp_agg.MP["group"])
    assert np.all(egt >= min_group)


# ═════════════════════════════════════════════════════════════
# na.rm behavior
# ═════════════════════════════════════════════════════════════

def test_na_rm_drops_na_and_proceeds(mp_agg):
    """aggte with na_rm=True drops NA ATTs and proceeds."""
    import copy
    obj = copy.deepcopy(mp_agg)
    obj.MP["att"] = list(obj.MP["att"])
    obj.MP["att"][0] = np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec="dynamic", na_rm=True)

    assert obj.atte is not None
    assert not np.isnan(obj.atte["overall_att"])


# ═════════════════════════════════════════════════════════════
# Cross-type consistency
# ═════════════════════════════════════════════════════════════

@pytest.mark.parametrize("typec", ["simple", "dynamic", "group", "calendar"])
def test_all_types_return_dict(mp_agg, typec):
    """All aggte types return a valid atte dict."""
    agg = _aggte(mp_agg, typec)
    assert isinstance(agg, dict)


@pytest.mark.parametrize("typec", ["simple", "dynamic", "group", "calendar"])
def test_all_types_have_didparams(mp_agg, typec):
    """All aggte types have non-None DIDparams."""
    agg = _aggte(mp_agg, typec)
    assert agg.get("DIDparams") is not None


def test_preserves_overridden_settings(mp_agg):
    """aggte preserves overridden bstrap/alp settings."""
    agg = _aggte(mp_agg, "dynamic", bstrap=False, alp=0.01)
    assert agg["DIDparams"]["bstrap"] is False
    assert agg["DIDparams"]["alp"] == 0.01


# ═════════════════════════════════════════════════════════════
# Cross-method tests: parametrize key checks across est_method
# ═════════════════════════════════════════════════════════════

EST_METHODS = ["dr", "reg", "ipw"]


@pytest.fixture(scope="module", params=EST_METHODS)
def mp_by_method(request):
    """ATTgt fitted with each est_method."""
    data = build_sim_data(n=1000, time_periods=4, te=1.0, seed=9142024)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method=request.param, bstrap=False)
    return obj


def test_simple_valid_across_methods(mp_by_method):
    """aggte simple returns valid overall ATT across est_methods."""
    agg = _aggte(mp_by_method, "simple")
    assert not np.isnan(agg["overall_att"])
    assert abs(agg["overall_att"] - 1) < 0.5


def test_dynamic_valid_across_methods(mp_by_method):
    """aggte dynamic returns event times across est_methods."""
    agg = _aggte(mp_by_method, "dynamic")
    egt = np.asarray(agg["egt"])
    assert len(egt) > 0
    assert not np.isnan(agg["overall_att"])


def test_group_valid_across_methods(mp_by_method):
    """aggte group returns per-group ATTs across est_methods."""
    agg = _aggte(mp_by_method, "group")
    assert not np.isnan(agg["overall_att"])
    assert abs(agg["overall_att"] - 1) < 0.5


def test_calendar_valid_across_methods(mp_by_method):
    """aggte calendar returns per-period ATTs across est_methods."""
    agg = _aggte(mp_by_method, "calendar")
    assert not np.isnan(agg["overall_att"])


def test_balance_e_across_methods(mp_by_method):
    """aggte dynamic balance_e filters groups across est_methods."""
    unbal = _aggte(mp_by_method, "dynamic")
    bal = _aggte(mp_by_method, "dynamic", balance_e=1)
    assert len(np.asarray(bal["egt"])) <= len(np.asarray(unbal["egt"]))


def test_na_rm_across_methods(mp_by_method):
    """aggte with na_rm=True across est_methods."""
    import copy
    obj = copy.deepcopy(mp_by_method)
    obj.MP["att"] = list(obj.MP["att"])
    obj.MP["att"][0] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec="dynamic", na_rm=True)
    assert obj.atte is not None
    assert not np.isnan(obj.atte["overall_att"])
