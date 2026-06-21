"""
Tests for att_gt function.
Translated from R package 'did' tests/testthat/test-att_gt.R (82 assertions).

Tests that require faster_mode, fix_weights, or custom est_method (R function)
are skipped — those features are not implemented in csdid Python.
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


@pytest.fixture(scope="module")
def mpdta():
    return pd.read_csv(os.path.join(DATA_DIR, "mpdta.csv"))


# ═════════════════════════════════════════════════════════════
# Basic estimation methods (Tests 1-2)
# ═════════════════════════════════════════════════════════════

def test_att_gt_dr_and_reg(sim_data):
    """R: att_gt works w/o dynamics, time effects, or group effects."""
    for em in ["dr", "reg"]:
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=sim_data, xformla="Y~X",
        ).fit(est_method=em, bstrap=False)
        att = np.asarray(res.MP["att"], dtype=float)
        assert abs(att[0] - 1) < 0.5, f"est_method={em}"


def test_att_gt_ipw(sim_data):
    """R: att_gt works using ipw."""
    for em in ["dr", "ipw"]:
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=sim_data,
        ).fit(est_method=em, bstrap=False)
        att = np.asarray(res.MP["att"], dtype=float)
        assert abs(att[0] - 1) < 0.5, f"est_method={em}"


# ═════════════════════════════════════════════════════════════
# Two-period case (Test 3)
# ═════════════════════════════════════════════════════════════

def test_two_period_case():
    """R: two period case — all 4 aggte types."""
    data = build_sim_data(n=5000, time_periods=2, te=1.0, seed=3)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method="reg", bstrap=False)

    for tp in ["simple", "dynamic", "group", "calendar"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res.aggte(typec=tp)
        assert abs(res.atte["overall_att"] - 1) < 0.5, f"type={tp}"


# ═════════════════════════════════════════════════════════════
# No covariates (Test 4)
# ═════════════════════════════════════════════════════════════

def test_no_covariates(sim_data):
    """R: no covariates case."""
    for em in ["dr", "reg"]:
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=sim_data,  # xformla=None → intercept only
        ).fit(est_method=em, bstrap=False)
        att = np.asarray(res.MP["att"], dtype=float)
        assert abs(att[0] - 1) < 0.5, f"est_method={em}"


# ═════════════════════════════════════════════════════════════
# Repeated cross sections (Tests 5-7)
# ═════════════════════════════════════════════════════════════

def test_repeated_cross_section(sim_data):
    """R: repeated cross section — dr and reg."""
    for em in ["dr", "reg"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=sim_data, xformla="Y~X", panel=False,
            ).fit(est_method=em, bstrap=False)
        att = np.asarray(res.MP["att"], dtype=float)
        assert abs(att[0] - 1) < 0.5, f"est_method={em}"


def test_ipw_repeated_cross_sections(sim_data):
    """R: ipw repeated cross sections."""
    for em in ["dr", "ipw"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=sim_data, xformla="Y~X", panel=False,
            ).fit(est_method=em, bstrap=False)
        att = np.asarray(res.MP["att"], dtype=float)
        assert abs(att[0] - 1) < 0.5, f"est_method={em}"


def test_rc_dynamic_effects():
    """R: repeated cross sections dynamic effects."""
    data = build_sim_data(n=2000, time_periods=4, te_e=[1, 2, 3, 4], seed=7)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", panel=False,
        ).fit(est_method="dr", bstrap=False)
        res.aggte(typec="dynamic")

    egt = np.asarray(res.atte["egt"], dtype=float)
    att_egt = np.asarray(res.atte["att_egt"], dtype=float)
    idx = np.where(egt == 2)[0]
    if len(idx) > 0:
        assert abs(att_egt[idx[0]] - 3) < 1.0  # te_e[2] = 3


# ═════════════════════════════════════════════════════════════
# Unbalanced panel (Test 8)
# ═════════════════════════════════════════════════════════════

def test_unbalanced_panel(sim_data):
    """R: unbalanced panel — dr with allow_unbalanced_panel.
    csdid now routes truly unbalanced panels through RC estimators (matches R)."""
    data = sim_data.copy()
    data = data.drop(data.index[1]).reset_index(drop=True)

    res = ATTgt(
        yname="Y", tname="period", idname="id", gname="G",
        data=data, xformla="Y~X", allow_unbalanced_panel=True,
    ).fit(est_method="dr", bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert abs(att[0] - 1) < 0.5


# ═════════════════════════════════════════════════════════════
# Not-yet-treated comparison group (Test 9)
# ═════════════════════════════════════════════════════════════

def test_notyettreated_rc(sim_data):
    """R: not yet treated comparison group — RC."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=sim_data, xformla="Y~X", panel=False,
            control_group="notyettreated",
        ).fit(est_method="dr", bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert abs(att[0] - 1) < 0.5


def test_notyettreated_no_nevertreated(sim_data):
    """R: no never-treated group, notyettreated control."""
    data = sim_data[sim_data["G"] > 0].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", panel=False,
            control_group="notyettreated",
        ).fit(est_method="dr", bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert abs(att[0] - 1) < 0.5


def test_nevertreated_coerces_no_nevertreated(sim_data):
    """R did v2.5.1: nevertreated control with no never-treated group warns and
    coerces the last cohort to never-treated (it does not error)."""
    data = sim_data[sim_data["G"] > 0].copy()
    with pytest.warns(UserWarning, match="never-treated"):
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, control_group="nevertreated",
        )
    assert (obj.dp["data"]["G"] == 0).any()


# ═════════════════════════════════════════════════════════════
# Aggregations (Test 10)
# ═════════════════════════════════════════════════════════════

def test_aggregation_dynamic():
    """R: dynamic aggregation with exposure-varying effects."""
    data = build_sim_data(n=2000, time_periods=4, te=0, te_e=[1, 2, 3, 4], seed=10)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", panel=False,
        ).fit(est_method="reg", bstrap=False)
        res.aggte(typec="dynamic")

    egt = np.asarray(res.atte["egt"], dtype=float)
    att_egt = np.asarray(res.atte["att_egt"], dtype=float)
    idx = np.where(egt == 2)[0]
    if len(idx) > 0:
        # At exposure 2: te_e[2] = 3
        assert abs(att_egt[idx[0]] - 3) < 1.0


def test_aggregation_balance_e():
    """R: dynamic aggregation with balance_e."""
    data = build_sim_data(n=2000, time_periods=4, te=0, te_e=[1, 2, 3, 4], seed=10)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", panel=False,
        ).fit(est_method="dr", bstrap=False)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic")
        egt_unbal = np.asarray(res.atte["egt"])

        res.aggte(typec="dynamic", balance_e=1)
        egt_bal = np.asarray(res.atte["egt"])

    assert len(egt_bal) <= len(egt_unbal)


# ═════════════════════════════════════════════════════════════
# Unequally spaced groups (Test 11)
# ═════════════════════════════════════════════════════════════

def test_unequally_spaced_groups():
    """R: unequally spaced groups."""
    data = build_sim_data(n=2000, time_periods=8, te=0, te_e=[1, 2, 3, 4, 5, 6, 7, 8], seed=11)
    keep_periods = [1, 2, 5, 7]
    data = data[data["period"].isin(keep_periods)].copy()
    data = data[(data["G"].isin([0] + keep_periods))].copy()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", panel=False,
        ).fit(est_method="reg", bstrap=False)
        res.aggte(typec="dynamic")

    egt = np.asarray(res.atte["egt"], dtype=float)
    att_egt = np.asarray(res.atte["att_egt"], dtype=float)
    idx = np.where(egt == 2)[0]
    if len(idx) > 0:
        assert abs(att_egt[idx[0]] - 3) < 1.5


# ═════════════════════════════════════════════════════════════
# Units treated in first period (Test 12)
# ═════════════════════════════════════════════════════════════

def test_first_period_treatment(sim_data):
    """R: some units treated in first period (dropped with a warning)."""
    data = sim_data.copy()
    data = data[data["period"] >= 2].reset_index(drop=True)

    with pytest.warns(UserWarning, match="already treated|[Dd]ropped"):
        ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", panel=False,
        ).fit(est_method="reg", bstrap=False)


# ═════════════════════════════════════════════════════════════
# Min / max length of exposures (Test 13)
# ═════════════════════════════════════════════════════════════

def test_min_max_exposures():
    """R: min and max length of exposures."""
    data = build_sim_data(n=2000, time_periods=4, te=0, te_e=[1, 2, 3, 4], seed=13)

    res = ATTgt(
        yname="Y", tname="period", idname="id", gname="G",
        data=data,
    ).fit(est_method="reg", bstrap=False, base_period="varying")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic", min_e=-1, max_e=1)

    egt = np.asarray(res.atte["egt"], dtype=float)
    att_egt = np.asarray(res.atte["att_egt"], dtype=float)
    idx = np.where(egt == 1)[0]
    if len(idx) > 0:
        assert abs(att_egt[idx[0]] - 2) < 0.5


# ═════════════════════════════════════════════════════════════
# Anticipation (Test 14)
# ═════════════════════════════════════════════════════════════

def test_anticipation():
    """R: anticipation parameter shifts base period correctly."""
    data = build_sim_data(n=2000, time_periods=5, te=0, te_e=[-1, 0, 1, 2, 3], seed=14)
    # Shift G forward by 1 to simulate anticipation
    data["G"] = data["G"].apply(lambda g: 0 if g == 0 else g + 1)
    data = data[data["G"] <= 5].copy()  # drop last period group

    # With correct anticipation=1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res1 = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", anticipation=1,
        ).fit(est_method="dr", bstrap=False)
        res1.aggte(typec="dynamic")

    egt1 = np.asarray(res1.atte["egt"], dtype=float)
    att1 = np.asarray(res1.atte["att_egt"], dtype=float)
    idx1 = np.where(egt1 == 2)[0]
    if len(idx1) > 0:
        assert abs(att1[idx1[0]] - 2) < 1.0

    # Without anticipation (anticipation=0) — overstates effects
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res0 = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", anticipation=0,
        ).fit(est_method="dr", bstrap=False)
        res0.aggte(typec="dynamic")

    egt0 = np.asarray(res0.atte["egt"], dtype=float)
    att0 = np.asarray(res0.atte["att_egt"], dtype=float)
    idx0 = np.where(egt0 == 2)[0]
    if len(idx0) > 0:
        assert abs(att0[idx0[0]] - 3) < 1.0


# ═════════════════════════════════════════════════════════════
# Significance level and confidence bands (Test 15)
# ═════════════════════════════════════════════════════════════

def test_significance_level(sim_data):
    """R: significance level and uniform confidence bands."""
    np.random.seed(1234)
    res05 = ATTgt(
        yname="Y", tname="period", idname="id", gname="G",
        data=sim_data, xformla="Y~X", alp=0.05,
    ).fit(est_method="dr", bstrap=False)

    np.random.seed(1234)
    res_pw = ATTgt(
        yname="Y", tname="period", idname="id", gname="G",
        data=sim_data, xformla="Y~X", alp=0.05, cband=False,
    ).fit(est_method="dr", bstrap=False)

    att05 = np.asarray(res05.results["att"], dtype=float)
    se05 = np.asarray(res05.results["se"], dtype=float)
    c05 = res05.results["c"]

    att_pw = np.asarray(res_pw.results["att"], dtype=float)
    se_pw = np.asarray(res_pw.results["se"], dtype=float)
    c_pw = res_pw.results["c"]

    # Uniform band (cband=True by default) should be at least as wide as pointwise
    upper05 = att05[0] + c05 * se05[0]
    upper_pw = att_pw[0] + c_pw * se_pw[0]
    assert upper05 >= upper_pw - 1e-10


# ═════════════════════════════════════════════════════════════
# Malformed data (Test 16)
# ═════════════════════════════════════════════════════════════

def test_malformed_data_bad_idname(sim_data):
    """R: incorrectly specified id."""
    with pytest.raises(Exception):
        ATTgt(
            yname="Y", tname="period", idname="brant", gname="G",
            data=sim_data, xformla="Y~X",
        )


# ═════════════════════════════════════════════════════════════
# Varying vs universal base period (Test 17)
# ═════════════════════════════════════════════════════════════

def test_varying_vs_universal_base():
    """R: varying or universal base period."""
    data = build_sim_data(n=2000, time_periods=8, te=0, te_e=[1, 2, 3, 4, 5, 6, 7, 8], seed=17)
    # Keep only early groups and add pre-treatment effects
    data = data[(data["G"] <= 5) | (data["G"] == 0)].copy()
    data["G"] = data["G"].apply(lambda g: 0 if g == 0 else g + 3)

    for bp in ["varying", "universal"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=data, xformla="Y~X",
            ).fit(est_method="dr", bstrap=False, base_period=bp)
            res.aggte(typec="dynamic")

        assert res.atte is not None
        assert not np.isnan(res.atte["overall_att"])


# ═════════════════════════════════════════════════════════════
# Small groups (Test 18)
# ═════════════════════════════════════════════════════════════

def test_small_groups(sim_data):
    """R: small groups — code should still compute but warn."""
    data = sim_data.copy()
    g2_ids = data.loc[data["G"] == 2, "id"].unique()
    keep_id = g2_ids[0]
    data = data[(data["G"] != 2) | (data["id"] == keep_id)].copy()

    with pytest.warns(UserWarning, match="small groups"):
        for em in ["dr", "reg"]:
            res = ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=data, xformla="Y~X",
            ).fit(est_method=em, bstrap=False)

            groups = np.asarray(res.MP["group"])
            t_arr = np.asarray(res.MP["t"])
            att = np.asarray(res.MP["att"], dtype=float)

            # Group 3 at t=3 should still be approximately 1
            idx = np.where((groups == 3) & (t_arr == 3))[0]
            if len(idx) > 0:
                assert abs(att[idx[0]] - 1) < 0.5, f"est_method={em}"


# ═════════════════════════════════════════════════════════════
# Sampling weights (Test 21)
# ═════════════════════════════════════════════════════════════

def test_sampling_weights(sim_data):
    """R: sampling weights — re-weighting should give same result as subset."""
    data = sim_data.copy()
    keep_ids = data["id"].unique()  # keep all
    data["w"] = 1 * data["id"].isin(keep_ids)
    data2 = data[data["id"].isin(keep_ids)].copy()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res_w = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", weights_name="w",
            control_group="notyettreated",
        ).fit(est_method="reg", bstrap=False)
        res_s = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data2, xformla="Y~X",
            control_group="notyettreated",
        ).fit(est_method="reg", bstrap=False)

    att_w = np.asarray(res_w.MP["att"], dtype=float)
    att_s = np.asarray(res_s.MP["att"], dtype=float)
    np.testing.assert_allclose(att_w[0], att_s[0], atol=1e-6)


# ═════════════════════════════════════════════════════════════
# Column naming (Test 22)
# ═════════════════════════════════════════════════════════════

def test_column_naming_gname(sim_data):
    """R: works when user column is literally named 'gname'."""
    data = sim_data.copy()
    data = data.rename(columns={"G": "gname", "period": "tname", "id": "idname"})

    res = ATTgt(
        yname="Y", tname="tname", idname="idname", gname="gname",
        data=data, xformla="Y~X",
    ).fit(est_method="reg", bstrap=False)

    att = np.asarray(res.MP["att"], dtype=float)
    assert not np.all(np.isnan(att))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="simple")
    assert not np.isnan(res.atte["overall_att"])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic")
    assert not np.isnan(res.atte["overall_att"])


# ═════════════════════════════════════════════════════════════
# Clustered SEs (Test 38 — partial, skip faster_mode parts)
# ═════════════════════════════════════════════════════════════

def test_clustered_se(sim_data):
    """R: clustered standard errors — numeric cluster variable."""
    data = sim_data.copy()
    data["cluster"] = data["cluster"].astype(float)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", clustervar="cluster",
        ).fit(est_method="dr", bstrap=False)

    att = np.asarray(res.MP["att"], dtype=float)
    assert abs(att[0] - 1) < 0.5


# ═════════════════════════════════════════════════════════════
# Skipped tests (R-only features)
# ═════════════════════════════════════════════════════════════

def test_custom_est_method(mpdta):
    """A custom Python callable est_method reproduces the built-in equivalent.

    csdid accepts a callable with the panel signature (y1, y0, D, i_weights,
    covariates) -> (att, inf_func).
    """
    from drdid import reg_did

    def my_reg(y1, y0, D, i_weights, covariates):
        return reg_did.reg_did_panel(y1, y0, D, i_weights=i_weights, covariates=covariates)

    common = dict(yname="lemp", tname="year", idname="countyreal",
                  gname="first.treat", data=mpdta, control_group="nevertreated",
                  xformla="lemp~1")
    res_custom = ATTgt(**common).fit(est_method=my_reg, bstrap=False)
    res_builtin = ATTgt(**common).fit(est_method="reg", bstrap=False)
    np.testing.assert_allclose(
        np.asarray(res_custom.MP["att"], dtype=float),
        np.asarray(res_builtin.MP["att"], dtype=float), atol=1e-9)


class TestFasterMode:
    """faster_mode must reproduce the standard path exactly (bit-identical)."""

    @staticmethod
    def _both(data, control_group="nevertreated", panel=True, xformla="Y~X",
              est_method="dr", base_period="varying", allow_unbalanced_panel=True):
        common = dict(yname="Y", tname="period", idname="id", gname="G", data=data,
                      control_group=control_group, panel=panel, xformla=xformla,
                      allow_unbalanced_panel=allow_unbalanced_panel)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            std = ATTgt(**common, faster_mode=False).fit(est_method, base_period=base_period, bstrap=False)
            fast = ATTgt(**common, faster_mode=True).fit(est_method, base_period=base_period, bstrap=False)
        a = np.asarray(std.MP["att"], dtype=float)
        b = np.asarray(fast.MP["att"], dtype=float)
        np.testing.assert_allclose(a, b, atol=1e-9, equal_nan=True)
        return std, fast

    def test_panel(self):
        self._both(build_sim_data(n=1000, time_periods=4, te=1.0, seed=1))

    def test_rcs(self):
        self._both(build_sim_data(n=1000, time_periods=4, te=1.0, seed=2), panel=False)

    def test_unbalanced(self):
        d = build_sim_data(n=1000, time_periods=4, te=1.0, seed=3)
        d = d.drop(d.index[[1, 7, 20]]).reset_index(drop=True)
        self._both(d, allow_unbalanced_panel=True)

    def test_filtered(self):
        self._both(build_sim_data(n=1200, time_periods=5, te=0, te_e=[1, 2, 3, 4, 5], seed=4),
                   control_group="notyettreated")

    def test_time_indexing_rc(self):
        self._both(build_sim_data(n=1000, time_periods=4, te=1.0, seed=5),
                   panel=False, base_period="universal")

    def test_time_indexing_panel(self):
        self._both(build_sim_data(n=1000, time_periods=4, te=1.0, seed=6), base_period="universal")

    def test_time_indexing_nonconsec(self):
        d = build_sim_data(n=1200, time_periods=6, te=0, te_e=[1, 2, 3, 4, 5, 6], seed=7)
        d = d[~d["period"].isin([3])].reset_index(drop=True)  # non-consecutive periods
        self._both(d)

    def test_time_indexing_universal(self):
        self._both(build_sim_data(n=1000, time_periods=5, te=0, te_e=[-1, 0, 1, 2, 3], seed=8),
                   base_period="universal")


class TestFixWeights:
    """fix_weights parameter (R did v2.5.x): None/varying/base_period/first_period."""

    @staticmethod
    def _tvw_data(seed=11):
        data = build_sim_data(n=1500, time_periods=4, te=1.0, seed=seed)
        rng = np.random.default_rng(seed)
        data["wt"] = 1.0 + 0.2 * data["period"] + rng.uniform(0, 0.5, len(data))
        return data

    def _fit(self, data, fix_weights=None, **kw):
        kw.setdefault("xformla", "Y~X")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return ATTgt(
                yname="Y", tname="period", idname="id", gname="G", data=data,
                weights_name="wt", fix_weights=fix_weights, **kw,
            ).fit(est_method="reg", bstrap=False)

    def test_fix_weights_options(self):
        data = self._tvw_data()
        for fw in [None, "varying", "base_period", "first_period"]:
            att = np.asarray(self._fit(data, fw).MP["att"], dtype=float)
            assert np.isfinite(att).any(), f"fix_weights={fw}"

    def test_time_invariant(self):
        # Time-invariant weights: every mode must give the same point estimates.
        data = build_sim_data(n=1200, time_periods=4, te=1.0, seed=7)
        data["wt"] = 1.0 + (data["id"] % 5)  # constant within unit
        a0 = np.asarray(self._fit(data, None).MP["att"], dtype=float)
        for fw, tol in [("base_period", 1e-10), ("first_period", 1e-10), ("varying", 1e-7)]:
            a = np.asarray(self._fit(data, fw).MP["att"], dtype=float)
            np.testing.assert_allclose(a, a0, atol=tol, err_msg=f"fix_weights={fw}")

    def test_tv_weights_match(self):
        # Time-varying weights: base_period and first_period should differ.
        data = self._tvw_data()
        a_base = np.asarray(self._fit(data, "base_period").MP["att"], dtype=float)
        a_first = np.asarray(self._fit(data, "first_period").MP["att"], dtype=float)
        assert not np.allclose(a_base, a_first)

    def test_nyt_tv_weights(self):
        data = self._tvw_data()
        att = np.asarray(self._fit(data, "base_period", control_group="notyettreated").MP["att"], dtype=float)
        assert np.isfinite(att).any()

    def test_rc_tv_weights(self):
        data = self._tvw_data()
        att = np.asarray(self._fit(data, "varying", panel=False).MP["att"], dtype=float)
        assert np.isfinite(att).any()

    def test_validation(self):
        data = self._tvw_data()
        with pytest.raises(ValueError, match="fix_weights"):
            ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                  weights_name="wt", fix_weights="bogus")
        # base_period / first_period not allowed for true repeated cross sections
        with pytest.raises(ValueError, match="repeated cross sections"):
            ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                  weights_name="wt", panel=False, fix_weights="base_period")

    def test_unbalanced(self):
        data = self._tvw_data()
        data = data.drop(data.index[1]).reset_index(drop=True)  # make unbalanced
        att = np.asarray(self._fit(data, "first_period", allow_unbalanced_panel=True).MP["att"], dtype=float)
        assert np.isfinite(att).any()


class TestIFConsistency:
    """Influence functions are identical between faster_mode and the standard path."""

    @staticmethod
    def _if_equal(data, control_group="nevertreated", panel=True, xformla="Y~X",
                  allow_unbalanced_panel=True):
        common = dict(yname="Y", tname="period", idname="id", gname="G", data=data,
                      control_group=control_group, panel=panel, xformla=xformla,
                      allow_unbalanced_panel=allow_unbalanced_panel)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            std = ATTgt(**common, faster_mode=False).fit("dr", bstrap=False)
            fast = ATTgt(**common, faster_mode=True).fit("dr", bstrap=False)
        ia = np.asarray(std.MP["inffunc"]["inffunc"], dtype=float)
        ib = np.asarray(fast.MP["inffunc"]["inffunc"], dtype=float)
        np.testing.assert_allclose(ia, ib, atol=1e-9)

    def test_balanced(self):
        self._if_equal(build_sim_data(n=1000, time_periods=4, te=1.0, seed=11))

    def test_notyettreated(self):
        self._if_equal(build_sim_data(n=1200, time_periods=5, te=0, te_e=[1, 2, 3, 4, 5], seed=12),
                       control_group="notyettreated")

    def test_rc(self):
        self._if_equal(build_sim_data(n=1000, time_periods=4, te=1.0, seed=13), panel=False)

    def test_unbalanced(self):
        d = build_sim_data(n=1000, time_periods=4, te=1.0, seed=14)
        d = d.drop(d.index[[2, 9]]).reset_index(drop=True)
        self._if_equal(d, allow_unbalanced_panel=True)

    def test_no_covariates(self):
        self._if_equal(build_sim_data(n=1000, time_periods=4, te=1.0, seed=15), xformla="Y~1")


# ═════════════════════════════════════════════════════════════
# Parametrized tests across all three est_methods
# ═════════════════════════════════════════════════════════════

EST_METHODS_ALL = ["dr", "reg", "ipw"]


@pytest.mark.parametrize("em", EST_METHODS_ALL)
def test_basic_estimation_all_methods(sim_data, em):
    """All three estimation methods produce valid ATTs with covariates."""
    res = ATTgt(
        yname="Y", tname="period", idname="id", gname="G",
        data=sim_data, xformla="Y~X",
    ).fit(est_method=em, bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert abs(att[0] - 1) < 0.5, f"est_method={em}"


@pytest.mark.parametrize("em", EST_METHODS_ALL)
def test_no_covariates_all_methods(sim_data, em):
    """All three methods work without covariates."""
    res = ATTgt(
        yname="Y", tname="period", idname="id", gname="G",
        data=sim_data,
    ).fit(est_method=em, bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert abs(att[0] - 1) < 0.5, f"est_method={em}"


@pytest.mark.parametrize("em", EST_METHODS_ALL)
def test_repeated_cross_section_all_methods(sim_data, em):
    """All three methods work with panel=False."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=sim_data, xformla="Y~X", panel=False,
        ).fit(est_method=em, bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert abs(att[0] - 1) < 0.5, f"est_method={em}"


@pytest.mark.parametrize("em", EST_METHODS_ALL)
def test_notyettreated_all_methods(sim_data, em):
    """All three methods work with notyettreated control group."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=sim_data, xformla="Y~X",
            control_group="notyettreated",
        ).fit(est_method=em, bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert abs(att[0] - 1) < 0.5, f"est_method={em}"


@pytest.mark.parametrize("em", EST_METHODS_ALL)
def test_anticipation_all_methods(em):
    """Anticipation parameter works with all three methods."""
    data = build_sim_data(n=2000, time_periods=5, te=0, te_e=[-1, 0, 1, 2, 3], seed=14)
    data["G"] = data["G"].apply(lambda g: 0 if g == 0 else g + 1)
    data = data[data["G"] <= 5].copy()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X", anticipation=1,
        ).fit(est_method=em, bstrap=False)
        res.aggte(typec="dynamic")

    assert res.atte is not None
    assert not np.isnan(res.atte["overall_att"])


@pytest.mark.parametrize("em", EST_METHODS_ALL)
def test_column_naming_variations(sim_data, em):
    """Works with user columns literally named 'group', 'time', etc."""
    data = sim_data.copy()
    data = data.rename(columns={"G": "group", "period": "time", "id": "unit"})
    res = ATTgt(
        yname="Y", tname="time", idname="unit", gname="group",
        data=data,
    ).fit(est_method=em, bstrap=False)
    att = np.asarray(res.MP["att"], dtype=float)
    assert not np.all(np.isnan(att))


@pytest.mark.parametrize("em", EST_METHODS_ALL)
def test_two_period_all_methods(em):
    """Two-period case works with all methods and all aggte types."""
    data = build_sim_data(n=5000, time_periods=2, te=1.0, seed=3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y~X",
        ).fit(est_method=em, bstrap=False)
    for tp in ["simple", "dynamic", "group", "calendar"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res.aggte(typec=tp)
        assert abs(res.atte["overall_att"] - 1) < 0.5, f"em={em}, type={tp}"
