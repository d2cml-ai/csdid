"""R did PR279: outcome units, observed balance support, and bootstrap edges.

Scale identities and explicit cohort/cell oracles check statistical behavior
without depending on the R interpreter or its random-number generator.
"""
import copy
import importlib
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from csdid.att_gt import ATTgt, _wald_pretest
from csdid.aggte_fnc.compute_aggte import compute_aggte
from csdid.utils.mboot import mboot


def _data():
    rng = np.random.default_rng(730)
    n, periods = 120, 4
    ids = np.repeat(np.arange(n), periods)
    t = np.tile(np.arange(1, periods + 1), n)
    g = np.repeat(np.repeat([0, 3, 4], n // 3), periods)
    x = np.repeat(rng.normal(size=n), periods)
    cluster = ids % 20
    shocks = rng.normal(size=(20, periods))
    y = (np.repeat(rng.normal(size=n), periods) + .3 * t + .2 * x * t +
         shocks[cluster, t - 1] + rng.normal(size=len(ids)) + ((g > 0) & (t >= g)))
    return pd.DataFrame(dict(id=ids, period=t, G=g, X=x, Y=y,
                             weight=np.exp(.2 * x), cluster=cluster))


def _fit(data, faster, bootstrap, clustered, method="reg", **kwargs):
    np.random.seed(731)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ATTgt("Y", "period", "id", "G", data=data, xformla="Y ~ X",
                     weights_name="weight", faster_mode=faster, cband=bootstrap,
                     biters=199, clustervar="cluster" if clustered else None,
                     **kwargs).fit(est_method=method, bstrap=bootstrap)


def _aggregate(fit, kind):
    np.random.seed(732)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fit.aggte(typec=kind).atte


@pytest.mark.parametrize("faster", [False, True])
@pytest.mark.parametrize("bootstrap", [False, True])
@pytest.mark.parametrize("clustered", [False, True])
@pytest.mark.parametrize("scale", [1e-12, 1e6])
@pytest.mark.parametrize("method", ["dr", "reg", "ipw"])
def test_outcome_units_preserve_inference(faster, bootstrap, clustered, scale, method):
    data = _data()
    reference = _fit(data, faster, bootstrap, clustered, method)
    scaled = _fit(data.assign(Y=data.Y * scale), faster, bootstrap, clustered, method)
    assert np.all(np.isfinite(scaled.results["se"]) & (scaled.results["se"] > 0))
    for field in ("att", "se", "l_se", "u_se"):
        np.testing.assert_allclose(np.asarray(scaled.results[field]) / scale,
                                   reference.results[field], atol=1e-10, rtol=1e-9)
    np.testing.assert_allclose(scaled.results["c"], reference.results["c"], atol=1e-10)
    if bootstrap and clustered:
        assert scaled.MP["W"] is reference.MP["W"] is None
    else:
        np.testing.assert_allclose(scaled.MP["W"], reference.MP["W"], atol=1e-9)
        assert scaled.MP["Wpval"] == reference.MP["Wpval"]
    for kind in ("simple", "dynamic", "group", "calendar"):
        expected, actual = _aggregate(reference, kind), _aggregate(scaled, kind)
        for field in ("overall_att", "overall_se"):
            np.testing.assert_allclose(actual[field] / scale, expected[field],
                                       atol=1e-10, rtol=1e-9)
        assert np.all(np.asarray(actual["overall_se"]) > 0)
        if kind != "simple":
            for field in ("att_egt", "se_egt"):
                np.testing.assert_allclose(np.asarray(actual[field]) / scale, expected[field],
                                           atol=1e-10, rtol=1e-9)
            np.testing.assert_allclose(actual["crit_val_egt"], expected["crit_val_egt"], atol=1e-10)


def _bootstrap_params(n=120, clustered=False, biters=199):
    return dict(data=pd.DataFrame(dict(id=np.arange(n), period=1,
                                      cluster=np.arange(n) % 20)),
                idname="id", tname="period", clustervars="cluster" if clustered else None,
                biters=biters, alp=.05, panel=True, true_repeated_cross_sections=False)


@pytest.mark.parametrize("clustered", [False, True])
@pytest.mark.parametrize("scale", [1e-6, 1e-100, 1e-200, 1e6])
def test_bootstrap_dimension_support_does_not_square_tiny_draws(clustered, scale):
    inf = np.column_stack([np.random.default_rng(501).normal(size=(120, 3)),
                           np.zeros(120), np.full(120, np.nan), np.full(120, np.inf)])
    params = _bootstrap_params(clustered=clustered)
    with np.errstate(invalid="ignore"):
        np.random.seed(502)
        reference = mboot(inf, params)
        np.random.seed(502)
        actual = mboot(inf * scale, params)
    assert actual["bres"].shape == reference["bres"].shape == (199, 3)
    np.testing.assert_allclose(actual["se"] / scale, reference["se"], rtol=1e-12)
    np.testing.assert_allclose(actual["bres"] / scale, reference["bres"], atol=1e-12)
    np.testing.assert_allclose(actual["crit_val"], reference["crit_val"], atol=1e-12)


def test_single_draw_retains_live_columns_and_marks_unavailable_scales():
    inf = np.column_stack([np.arange(1, 121), 2 * np.arange(1, 121),
                           3 * np.arange(1, 121), np.zeros(120), np.full(120, np.inf)])
    np.random.seed(834)
    with np.errstate(invalid="ignore"):
        actual = mboot(inf, _bootstrap_params(biters=1))
    assert actual["bres"].shape == (1, 3)
    assert actual["se"].shape == (5,)
    assert np.isnan(actual["se"]).all()
    assert np.isnan(actual["crit_val"])


@pytest.mark.parametrize("faster", [False, True])
def test_single_draw_fit_falls_back_to_pointwise_but_retains_analytic_wald(faster):
    np.random.seed(835)
    with pytest.warns(UserWarning, match="Falling back to pointwise"):
        fit = ATTgt("Y", "period", "id", "G", data=_data(), xformla="Y ~ X",
                    faster_mode=faster, cband=True, biters=1).fit(bstrap=True)
    assert np.isnan(fit.results["se"]).all()
    assert not fit.dp["cband"]
    assert fit.results["c"] == norm.ppf(.975)
    assert np.isfinite(fit.MP["W"])


@pytest.mark.parametrize("faster", [False, True])
@pytest.mark.parametrize("bootstrap", [False, True])
@pytest.mark.parametrize("nonfinite", [0., np.inf])
@pytest.mark.parametrize("panel", [False, True])
def test_zero_and_nonfinite_custom_influence_never_produce_inference(faster, bootstrap, nonfinite, panel):
    def estimator(y1=None, y0=None, D=None, i_weights=None, covariates=None,
                  y=None, post=None):
        return 0., np.full(len(D), nonfinite)
    data = _data()
    if not panel:
        data["id"] = np.arange(len(data))
    fit = _fit(data, faster, bootstrap, False, method=estimator, panel=panel)
    assert np.isnan(fit.results["se"]).all()
    assert fit.MP["W"] is None
    assert np.isfinite(fit.results["c"])
    if bootstrap:
        assert not fit.dp["cband"]
    for kind in ("simple", "dynamic", "group", "calendar"):
        assert np.isnan(_aggregate(fit, kind)["overall_se"]).all()


@pytest.mark.parametrize("faster", [False, True])
@pytest.mark.parametrize("bootstrap", [False, True])
def test_low_variance_cohort_remains_available_beside_normal_variance_cohort(faster, bootstrap):
    data = _data()
    rng = np.random.default_rng(733)
    noise_scale = np.where(data.G == 4, 1., 1e-11)
    data["Y"] = (rng.normal(size=120)[data.id] + noise_scale * rng.normal(size=len(data)) +
                 ((data.G > 0) & (data.period >= data.G)))
    fit = _fit(data, faster, bootstrap, False, method="dr")
    small = np.asarray(fit.results["se"])[np.asarray(fit.MP["group"]) == 3]
    large = np.asarray(fit.results["se"])[np.asarray(fit.MP["group"]) == 4]
    assert np.all(np.isfinite(small) & (small > 0))
    assert small.max() < 1e-10
    assert large.min() > 1e-2


def _balance_mp(periods=(1, 2, 3, 4, 5), cohorts=(2, 3)):
    groups = np.repeat(cohorts, len(periods) - 1)
    times = np.tile(periods[1:], len(cohorts))
    unit_groups = np.repeat([0, *cohorts], 10)
    data = pd.DataFrame(dict(id=np.repeat(np.arange(len(unit_groups)), len(periods)),
                             period=np.tile(periods, len(unit_groups)),
                             G=np.repeat(unit_groups, len(periods)), w1=1.))
    inf = np.random.default_rng(711).normal(size=(len(unit_groups), len(groups)))
    inf -= inf.mean(axis=0)
    return dict(group=groups, t=times, att=groups + .1 * (times - groups), n=len(unit_groups),
                inffunc={"inffunc": inf}, DIDparams=dict(
                    data=data, idname="id", tname="period", gname="G", panel=True,
                    tlist=np.array(periods), glist=np.array(cohorts), bstrap=False,
                    cband=False, biters=99, clustervars=None, alp=.05,
                    control_group="nevertreated", est_method="reg", anticipation=0))


@pytest.mark.parametrize("max_e", [0, np.inf])
def test_balancing_precedes_display_limits_and_drops_incomplete_cohort(max_e):
    mp = _balance_mp()
    mp["att"][(mp["group"] == 2) & (mp["t"] == 3)] = np.nan
    actual = compute_aggte(copy.deepcopy(mp), typec="dynamic", balance_e=1,
                           min_e=0, max_e=max_e, na_rm=True)
    selected = (mp["group"] == 3) & (mp["t"] >= 3) & (mp["t"] <= min(4, 3 + max_e))
    np.testing.assert_equal(actual["egt"], mp["t"][selected] - 3)
    np.testing.assert_allclose(actual["att_egt"], mp["att"][selected])
    np.testing.assert_allclose(actual["inf_function"]["dynamic_inf_func_e"],
                               mp["inffunc"]["inffunc"][:, selected], atol=1e-14)


def test_irregular_grid_balance_requires_common_observed_support():
    mp = _balance_mp(periods=(1, 2, 4, 5), cohorts=(2, 4))
    actual = compute_aggte(copy.deepcopy(mp), typec="dynamic", balance_e=1, min_e=0)
    selected = (mp["group"] == 4) & (mp["t"] >= 4)
    np.testing.assert_equal(actual["egt"], [0, 1])
    np.testing.assert_allclose(actual["att_egt"], mp["att"][selected])


def test_globally_unobserved_event_time_is_not_required():
    mp = _balance_mp(periods=(1, 3, 5), cohorts=(3, 5))
    actual = compute_aggte(copy.deepcopy(mp), typec="dynamic", balance_e=2, min_e=0)
    selected = mp["group"] == 3
    np.testing.assert_equal(actual["egt"], [0, 2])
    np.testing.assert_allclose(actual["att_egt"], mp["att"][selected])


@pytest.mark.parametrize("balance_e", [1, 10])
def test_all_missing_event_support_cannot_disappear_during_na_removal(balance_e):
    mp = _balance_mp()
    mp["att"][mp["t"] - mp["group"] == 1] = np.nan
    with pytest.raises(ValueError, match="No groups have non-missing ATT"):
        compute_aggte(mp, typec="dynamic", balance_e=balance_e, na_rm=True)


def test_anticipation_cells_are_excluded_from_wald_restrictions():
    # Orthogonal columns make the exact Wald sum transparent: only the first
    # cell precedes g-a; the large anticipation ATT must not enter the test.
    influence = np.eye(3)
    result = _wald_pretest(influence, [2., 100., 200.], [5, 5, 5], [2, 3, 4],
                           dict(anticipation=2, clustervars=None, bstrap=False))
    assert result[0] == 36.


@pytest.mark.parametrize("cohort", [4, 5])
def test_universal_base_warning_counts_actual_irregular_time_grid(cohort):
    # Last observed period before g-a is 1, rather than the nonexistent g-a-1.
    with pytest.warns(UserWarning, match="No pre-treatment ATT.*exist"):
        result = _wald_pretest(np.zeros((1, 3)), [0.], [cohort], [1],
                               dict(anticipation=2, base_period="universal",
                                    tlist=[1, 4, 8], clustervars=None, bstrap=False))
    assert result == (None, None)


def test_zero_iqr_column_does_not_remove_other_studentized_columns(monkeypatch):
    module = importlib.import_module("csdid.utils.mboot")
    draws = np.column_stack([np.arange(1., 200.), np.ones(199)])
    monkeypatch.setattr(module, "run_multiplier_bootstrap", lambda *args: draws)
    actual = mboot(np.ones((120, 2)), _bootstrap_params())
    assert np.isfinite(actual["se"][0])
    assert np.isnan(actual["se"][1])
    assert np.isfinite(actual["crit_val"])


def test_balanced_dynamic_stress_against_explicit_cell_ratio_oracle():
    rng = np.random.default_rng(20261026)
    for _ in range(200):
        periods = tuple(np.r_[1, np.sort(rng.choice(np.arange(2, 11), 4, replace=False))])
        cohorts = tuple(rng.choice(periods[1:4], 2, replace=False))
        mp = _balance_mp(periods=periods, cohorts=cohorts)
        full_g, full_t = mp["group"], mp["t"]
        mp["att"][rng.random(len(full_t)) < .25] = np.nan
        horizon = int(rng.integers(0, 12))
        max_t = max(full_t)
        eligible = {g for g in cohorts if max_t - g >= horizon}
        required = {t - g for g, t in zip(full_g, full_t)
                    if g in eligible and 0 <= t - g <= horizon}
        retained = [g for g in eligible if required and all(
            any(gg == g and tt - gg == event and np.isfinite(value)
                for gg, tt, value in zip(full_g, full_t, mp["att"]))
            for event in required)]
        if not retained:
            with pytest.raises(ValueError, match="No groups have non-missing ATT"):
                compute_aggte(copy.deepcopy(mp), typec="dynamic", balance_e=horizon,
                               min_e=0, na_rm=True)
            continue
        actual = compute_aggte(copy.deepcopy(mp), typec="dynamic", balance_e=horizon,
                               min_e=0, na_rm=True)
        assert set(actual["egt"]) == required
        first = mp["DIDparams"]["data"].drop_duplicates("id")
        unit_groups = first.G.to_numpy()
        for index, event in enumerate(actual["egt"]):
            selected = np.array([g in retained and t - g == event
                                 for g, t in zip(full_g, full_t)])
            values = mp["att"][selected]
            mass_by_unit = (unit_groups[:, None] == full_g[selected]).astype(float)
            mass = mass_by_unit.mean(axis=0)
            mean = mass @ values / mass.sum()
            influence = mp["inffunc"]["inffunc"][:, selected] @ (mass / mass.sum())
            influence += mass_by_unit @ (values - mean) / mass.sum()
            np.testing.assert_allclose(actual["att_egt"][index], mean, atol=1e-14)
            np.testing.assert_allclose(actual["inf_function"]["dynamic_inf_func_e"][:, index],
                                       influence, atol=1e-14)
            np.testing.assert_allclose(np.asarray(actual["se_egt"]).ravel()[index],
                                       np.linalg.norm(influence) / mp["n"], atol=1e-14)
