"""
Parametric combination tests for csdid.

Systematically tests across major axes: est_method, control_group, base_period,
panel, anticipation, bstrap, and aggte type. Uses a single data fixture to
avoid redundant generation.
"""
import copy
import warnings

import numpy as np
import pytest

from csdid.att_gt import ATTgt
from helpers import build_sim_data


# ═════════════════════════════════════════════════════════════
# Shared data fixtures (module-scoped for performance)
# ═════════════════════════════════════════════════════════════

@pytest.fixture(scope="module")
def panel_data():
    """Standard balanced panel with constant te=1."""
    return build_sim_data(n=1000, time_periods=4, te=1.0, seed=9142024)


@pytest.fixture(scope="module")
def dynamic_data():
    """Panel with exposure-varying treatment effects."""
    return build_sim_data(n=2000, time_periods=4, te=0, te_e=[1, 2, 3, 4], seed=10)


# ═════════════════════════════════════════════════════════════
# Combo A: est_method × control_group × base_period (3×2×2 = 12)
# Each combo gets 3 checks → 36 tests
# ═════════════════════════════════════════════════════════════

EST_METHODS = ["dr", "reg", "ipw"]
CONTROL_GROUPS = ["nevertreated", "notyettreated"]
BASE_PERIODS = ["varying", "universal"]

_combo_a_params = [
    (em, cg, bp)
    for em in EST_METHODS
    for cg in CONTROL_GROUPS
    for bp in BASE_PERIODS
]


@pytest.mark.parametrize("em,cg,bp", _combo_a_params,
                         ids=[f"{em}-{cg}-{bp}" for em, cg, bp in _combo_a_params])
class TestComboA:
    """est_method × control_group × base_period combinations."""

    def _fit(self, panel_data, em, cg, bp):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=panel_data, xformla="Y~X", control_group=cg,
            ).fit(est_method=em, bstrap=False, base_period=bp)

    def test_att_gt_valid(self, panel_data, em, cg, bp):
        """ATT(g,t) produces valid (non-all-NaN) results."""
        res = self._fit(panel_data, em, cg, bp)
        att = np.asarray(res.MP["att"], dtype=float)
        assert np.any(~np.isnan(att))

    def test_se_positive_finite(self, panel_data, em, cg, bp):
        """SEs are positive and finite where ATT is non-trivially defined."""
        res = self._fit(panel_data, em, cg, bp)
        se = np.asarray(res.results["se"], dtype=float)
        att = np.asarray(res.results["att"], dtype=float)
        # universal base period produces SE=0 for reference periods (expected)
        valid = np.isfinite(att) & np.isfinite(se) & (se > 0)
        assert np.sum(valid) > 0

    def test_aggte_simple(self, panel_data, em, cg, bp):
        """aggte('simple') produces valid overall ATT."""
        res = self._fit(panel_data, em, cg, bp)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res.aggte(typec="simple", bstrap=False, na_rm=True)
        assert not np.isnan(res.atte["overall_att"])


# ═════════════════════════════════════════════════════════════
# Combo B: est_method × panel (3×2 = 6)
# Each combo gets 2 checks → 12 tests
# ═════════════════════════════════════════════════════════════

_combo_b_params = [
    (em, panel)
    for em in EST_METHODS
    for panel in [True, False]
]


@pytest.mark.parametrize("em,panel", _combo_b_params,
                         ids=[f"{em}-panel{panel}" for em, panel in _combo_b_params])
class TestComboB:
    """est_method × panel combinations."""

    def _fit(self, panel_data, em, panel):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=panel_data, xformla="Y~X", panel=panel,
            ).fit(est_method=em, bstrap=False)

    def test_produces_results(self, panel_data, em, panel):
        """Both panel=True and panel=False produce ATT results."""
        res = self._fit(panel_data, em, panel)
        att = np.asarray(res.MP["att"], dtype=float)
        assert np.any(~np.isnan(att))

    def test_aggte_dynamic(self, panel_data, em, panel):
        """aggte('dynamic') works for both panel settings."""
        res = self._fit(panel_data, em, panel)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res.aggte(typec="dynamic", bstrap=False)
        assert res.atte is not None
        assert not np.isnan(res.atte["overall_att"])


# ═════════════════════════════════════════════════════════════
# Combo C: anticipation × est_method (3×3 = 9)
# Each combo gets 2 checks → 18 tests
# ═════════════════════════════════════════════════════════════

ANTICIPATIONS = [0, 1, 2]

_combo_c_params = [
    (ant, em)
    for ant in ANTICIPATIONS
    for em in EST_METHODS
]


@pytest.fixture(scope="module")
def anticipation_data():
    """Data with enough periods for anticipation=2."""
    return build_sim_data(n=1000, time_periods=6, te=1.0, seed=42)


@pytest.mark.parametrize("ant,em", _combo_c_params,
                         ids=[f"ant{ant}-{em}" for ant, em in _combo_c_params])
class TestComboC:
    """anticipation × est_method combinations."""

    def _fit(self, anticipation_data, ant, em):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=anticipation_data, xformla="Y~X", anticipation=ant,
            ).fit(est_method=em, bstrap=False)

    def test_produces_results(self, anticipation_data, ant, em):
        """ATTgt with anticipation parameter produces valid results."""
        res = self._fit(anticipation_data, ant, em)
        att = np.asarray(res.MP["att"], dtype=float)
        assert np.any(~np.isnan(att))

    def test_se_finite(self, anticipation_data, ant, em):
        """SEs are finite with anticipation parameter."""
        res = self._fit(anticipation_data, ant, em)
        se = np.asarray(res.results["se"], dtype=float)
        assert np.any(np.isfinite(se) & (se > 0))


# ═════════════════════════════════════════════════════════════
# Combo D: bstrap × est_method (2×3 = 6)
# Each combo gets 2 checks → 12 tests
# ═════════════════════════════════════════════════════════════

_combo_d_params = [
    (bstrap, em)
    for bstrap in [True, False]
    for em in EST_METHODS
]


@pytest.mark.parametrize("bstrap,em", _combo_d_params,
                         ids=[f"bstrap{bstrap}-{em}" for bstrap, em in _combo_d_params])
class TestComboD:
    """bstrap × est_method combinations."""

    def test_se_finite(self, panel_data, bstrap, em):
        """SEs are finite with both bstrap=True and bstrap=False."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=panel_data, xformla="Y~X", biters=100,
            ).fit(est_method=em, bstrap=bstrap)
        se = np.asarray(res.results["se"], dtype=float)
        assert np.any(np.isfinite(se) & (se > 0))

    def test_att_reasonable(self, panel_data, bstrap, em):
        """ATT estimates are reasonable regardless of bstrap setting."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(
                yname="Y", tname="period", idname="id", gname="G",
                data=panel_data, xformla="Y~X", biters=100,
            ).fit(est_method=em, bstrap=bstrap)
        att = np.asarray(res.MP["att"], dtype=float)
        non_na = att[~np.isnan(att)]
        assert len(non_na) > 0
        assert np.all(np.isfinite(non_na))


# ═════════════════════════════════════════════════════════════
# Combo E: Full integration combos (select 24 representative)
# est_method × control_group × base_period × panel
# ═════════════════════════════════════════════════════════════

_combo_e_params = [
    (em, cg, bp, panel)
    for em in EST_METHODS
    for cg in CONTROL_GROUPS
    for bp in BASE_PERIODS
    for panel in [True, False]
]


@pytest.mark.parametrize("em,cg,bp,panel", _combo_e_params,
                         ids=[f"{em}-{cg}-{bp}-panel{panel}"
                              for em, cg, bp, panel in _combo_e_params])
def test_full_integration(panel_data, em, cg, bp, panel):
    """Full pipeline: fit → aggte → check for each combo."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=panel_data, xformla="Y~X", panel=panel,
            control_group=cg,
        ).fit(est_method=em, bstrap=False, base_period=bp)
        res.aggte(typec="simple", bstrap=False, na_rm=True)
    assert not np.isnan(res.atte["overall_att"])


# ═════════════════════════════════════════════════════════════
# Combo F: aggte types across methods (4×3 = 12)
# ═════════════════════════════════════════════════════════════

AGGTE_TYPES = ["simple", "dynamic", "group", "calendar"]

_combo_f_params = [
    (typec, em)
    for typec in AGGTE_TYPES
    for em in EST_METHODS
]


@pytest.mark.parametrize("typec,em", _combo_f_params,
                         ids=[f"{typec}-{em}" for typec, em in _combo_f_params])
def test_aggte_type_by_method(panel_data, typec, em):
    """Each aggte type works with each est_method."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=panel_data, xformla="Y~X",
        ).fit(est_method=em, bstrap=False)
        res.aggte(typec=typec, bstrap=False)
    assert res.atte is not None
    assert not np.isnan(res.atte["overall_att"])
    se = res.atte["overall_se"]
    assert np.isfinite(se) and se > 0


# ═════════════════════════════════════════════════════════════
# Combo G: RCS dynamic effects across methods (3 methods × 2 checks = 6)
# ═════════════════════════════════════════════════════════════

@pytest.mark.parametrize("em", EST_METHODS)
def test_rcs_dynamic_effects_by_method(dynamic_data, em):
    """RCS dynamic effects are detected across all methods."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=dynamic_data, xformla="Y~X", panel=False,
        ).fit(est_method=em, bstrap=False)
        res.aggte(typec="dynamic", bstrap=False)
    egt = np.asarray(res.atte["egt"], dtype=float)
    att_egt = np.asarray(res.atte["att_egt"], dtype=float)
    # Post-treatment effects should be positive
    post = att_egt[egt >= 0]
    assert np.any(post[~np.isnan(post)] > 0)


@pytest.mark.parametrize("em", EST_METHODS)
def test_rcs_vs_panel_both_work(panel_data, em):
    """Both panel and RCS produce finite results for same data."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res_p = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=panel_data, xformla="Y~X", panel=True,
        ).fit(est_method=em, bstrap=False)
        res_r = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=panel_data, xformla="Y~X", panel=False,
        ).fit(est_method=em, bstrap=False)
    att_p = np.asarray(res_p.MP["att"], dtype=float)
    att_r = np.asarray(res_r.MP["att"], dtype=float)
    assert np.any(~np.isnan(att_p))
    assert np.any(~np.isnan(att_r))
