"""Weighted feasibility regressions from R did PRs #278/#279.

The public-fit oracles use deletion of zero-mass observations and direct 2x2
estimators; count/normalization boundaries are checked without unsafe fits.
"""
import importlib
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from csdid.attgt_fnc import drdid_trim
from csdid.attgt_fnc.compute_att_gt_shared import (
    overlap_check_fail, rcond_check_fail, rc_outcome_has_too_few_obs,
)


def _panel():
    rng = np.random.default_rng(50)
    units = pd.DataFrame({"id": np.arange(100), "G": np.repeat([0, 3], 50),
                          "X": rng.normal(size=100), "weight": 1.,
                          "a": rng.normal(size=100)})
    data = units.merge(pd.DataFrame({"period": np.arange(1, 5)}), how="cross")
    data["Y"] = (data.a + data.period * (.2 + .3 * data.X)
                 + rng.normal(size=len(data)) + 2 * ((data.G == 3) & (data.period >= 3)))
    return data


def _fit(data, mode, method="dr", route="panel", formula="~X"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ATTgt(
            yname="Y", tname="period", idname="id", gname="G", data=data,
            xformla=formula, weights_name="weight", faster_mode=mode,
            panel=route != "rc", allow_unbalanced_panel=route == "unbalanced",
            fix_weights="varying" if route == "varying" else None,
        ).fit(est_method=method, bstrap=False)


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("method", ["reg", "dr", "ipw"])
@pytest.mark.parametrize("route", ["panel", "varying", "rc", "unbalanced"])
def test_zero_mass_extreme_covariate_matches_deleting_observation(mode, method, route):
    data = _panel()
    if route == "unbalanced":
        data = data.loc[(data.id * 7 + data.period) % 17 != 0]
    extra = pd.DataFrame({"id": 100, "G": 0 if method == "reg" else 3,
                          "X": 1e12 if method == "reg" else 1000.,
                          "weight": 0., "a": 0., "period": np.arange(1, 5), "Y": 0.})
    augmented = pd.concat([data, extra], ignore_index=True)
    if route == "rc":
        data = data.copy()
        data["id"] = np.arange(len(data))
        augmented["id"] = np.arange(len(augmented))
    expected = _fit(data, mode, method, route)
    actual = _fit(augmented, mode, method, route)
    assert np.isfinite(actual.MP["att"]).all()
    np.testing.assert_allclose(actual.MP["att"], expected.MP["att"], atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(actual.results["se"], expected.results["se"], atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("mode", [False, True])
def test_weighted_control_gram_keeps_tiny_weight_extreme_covariate(mode):
    data = _panel()
    extra = pd.DataFrame({"id": 100, "G": 0, "X": 1e12,
                          "weight": 1e-24, "a": 0., "period": np.arange(1, 5), "Y": 0.})
    data = pd.concat([data, extra], ignore_index=True)
    fit = _fit(data, mode, "reg")
    pre = data[data.period == 2].sort_values("id")
    post = data[data.period == 3].sort_values("id")
    from drdid.reg_did import reg_did_panel
    att, _ = reg_did_panel(post.Y.to_numpy(), pre.Y.to_numpy(), (pre.G == 3).to_numpy(),
                          covariates=np.column_stack([np.ones(len(pre)), pre.X]),
                          i_weights=pre.weight.to_numpy())
    assert np.isfinite(fit.MP["att"]).all()
    np.testing.assert_allclose(np.asarray(fit.MP["att"])[np.asarray(fit.MP["t"]) == 3], att,
                               atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("scale", [1e-100, 1., 1e100])
def test_weighted_control_gram_rejects_lost_direction_independent_of_units(scale):
    x = np.column_stack([np.ones(6), [-1, 1, -1, 1, 0, 0]])
    d = np.array([0, 0, 0, 0, 1, 1])
    w = np.array([1, 1e-20, 1, 1e-20, 1, 1]) * scale
    assert not rcond_check_fail(x, d)
    assert rcond_check_fail(x, d, i_weights=w)
    assert not rcond_check_fail(x, d, i_weights=np.full(6, scale))


def test_overlap_uses_positive_normalized_support_and_empty_support_fails():
    x = np.column_stack([np.ones(4), [-1, 0, 1, 1000]])
    d = np.array([0, 1, 0, 1])
    w = np.array([1e200, 1e200, 1e200, 1e-200])
    assert overlap_check_fail(x, d, i_weights=w) == overlap_check_fail(x[:3], d[:3])
    assert overlap_check_fail(np.ones((3, 1)), [0, 1, 1], i_weights=[0, 1, 1])
    assert rcond_check_fail(np.ones((3, 1)), [0, 1, 1], i_weights=[0, 1, 1])
    assert overlap_check_fail(np.ones((3, 1)), [0, 1, 1], i_weights=[0, 0, 0])


@pytest.mark.parametrize("scale", [1., 2.**-900])
def test_rc_support_counts_weights_after_local_normalization(scale):
    d, post = np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1])
    w = np.array([2.**1000, 2.**1000, 2.**-100, 2.**-100]) * scale
    assert np.all(np.isfinite(w) & (w > 0))
    np.testing.assert_array_equal(w / np.mean(w), [2, 2, 0, 0])
    assert rc_outcome_has_too_few_obs(d, post, w, 1)
    assert not rc_outcome_has_too_few_obs(d, post, w, 1, include_treated=False)
    assert rc_outcome_has_too_few_obs(d, post, w[::-1], 1, include_treated=False)


def test_rc_support_accepts_smallest_positive_subnormal_and_square_design():
    d, post = np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1])
    retained = np.array([2.**1000, 2.**1000, 2.**-75, 2.**-75])
    lost = np.array([2.**1000, 2.**1000, 2.**-76, 2.**-76])
    assert (retained / np.mean(retained))[2] == np.nextafter(0., 1.)
    assert not rc_outcome_has_too_few_obs(d, post, retained, 1)
    assert rc_outcome_has_too_few_obs(d, post, lost, 1)
    for level in (2.**-1000, 1., 2.**1000):
        assert not rc_outcome_has_too_few_obs(d, post, np.full(4, level), 1)


def _tiny_rc(n_treated=3):
    cube = np.array(np.meshgrid([-1., 1.], [-1., 1.], [-1., 1.])).reshape(3, -1).T
    regular = np.tile(cube, (2, 1))
    tiny = (np.array([[-1, -1, -1], [1, 1, -1], [0, 0, 1]]) if n_treated == 3
            else np.vstack([np.zeros(3), np.eye(3)]))
    frames = []
    for period in range(1, 4):
        for g in (0, 2, 3):
            x = tiny if g == 2 else regular
            frame = pd.DataFrame(x, columns=["X", "X2", "X3"])
            frame["period"], frame["G"], frame["weight"] = period, g, 1.
            frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    data["id"] = np.arange(len(data))
    data["Y"] = (.2 * data.period + .3 * data.X - .1 * data.X2 + .2 * data.X3
                 + ((data.G > 0) & (data.period >= data.G)) + .1 * np.sin(data.id))
    return data


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("method", ["dr", "reg", "ipw"])
@pytest.mark.parametrize("n_treated", [3, 4])
def test_rc_tiny_and_square_treated_designs_gate_only_required_outcome_fits(
        monkeypatch, mode, method, n_treated):
    data = _tiny_rc(n_treated)
    original = drdid_trim.drdid_rc
    entered = []

    def backend(y, post, D, **kwargs):
        entered.append(np.sum(D))
        assert np.sum(D) != 6, "undersized DR outcome design reached backend"
        return original(y=y, post=post, D=D, **kwargs)

    monkeypatch.setattr(drdid_trim, "drdid_rc", backend)
    fit = _fit(data, mode, method, "rc", "~X + X2 + X3")
    failed = (np.asarray(fit.MP["group"]) == 2) & (method == "dr") & (n_treated == 3)
    np.testing.assert_array_equal(np.isnan(fit.MP["att"]), failed)
    np.testing.assert_array_equal(np.isnan(fit.results["se"]), failed)
    np.testing.assert_array_equal(np.isnan(fit.MP["inffunc"]["inffunc"]).all(axis=0), failed)
    if method == "dr":
        assert len(entered) == (2 if n_treated == 3 else 4)


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("route", ["rc", "varying", "unbalanced"])
def test_second_local_normalization_is_checked_before_rc_backend(monkeypatch, mode, route):
    rows = []
    for d in (0, 1):
        for period in range(1, 5):
            for j in range(16):
                rows.append((1000 * (d + 1) + j, 2 * d, period,
                             6. if d == 0 and period <= 2 else np.nextafter(0., 1.)))
    data = pd.DataFrame(rows, columns=["id", "G", "period", "weight"])
    data["Y"] = data.period.astype(float)
    if route == "rc":
        data["id"] = np.arange(len(data))
    elif route == "unbalanced":
        extra = data.iloc[[-1]].copy()
        extra["id"] = 3000
        data = pd.concat([data, extra], ignore_index=True)
    seen = []

    def backend(y, post, **kwargs):
        seen.extend(np.unique(y[np.asarray(post) == 1]))
        raise ValueError("safe non-target backend probe")

    monkeypatch.setattr(drdid_trim, "drdid_rc", backend)
    # Isolate the physical-size preflight from the separate downstream Gram
    # conditioning guard, whose division can also underflow at this boundary.
    engine = importlib.import_module(
        "csdid.attgt_fnc.compute_att_gt2" if mode else "csdid.attgt_fnc.compute_att_gt")
    monkeypatch.setattr(engine, "drdid_design_singular", lambda *a, **k: False)
    fit = _fit(data, mode, "dr", route, "~1")
    assert seen == [3., 4.]
    assert np.isnan(np.asarray(fit.MP["att"])[np.asarray(fit.MP["t"]) == 2]).all()


@pytest.mark.parametrize("method", ["dr", "ipw"])
def test_fast_overlap_cache_tracks_time_varying_positive_support(method):
    data = _panel()
    # X is invariant, while the zero-mass outlier becomes supported at period2.
    extra = pd.DataFrame({"id": 100, "G": 3, "X": 1000., "a": 0.,
                          "period": np.arange(1, 5), "Y": 0.,
                          "weight": [0., 1., 1., 1.]})
    data = pd.concat([data, extra], ignore_index=True)
    expected = _fit(data, False, method)
    actual = _fit(data, True, method)
    np.testing.assert_array_equal(actual.MP["att"], expected.MP["att"])
    np.testing.assert_array_equal(actual.MP["inffunc"]["inffunc"], expected.MP["inffunc"]["inffunc"])
    assert np.isfinite(actual.MP["att"][0])
    assert np.isnan(actual.MP["att"][1:]).all()


def test_intercept_weighted_rank_uses_positive_support_without_underflowing_gram():
    # Dividing the scalar Gram by the cell size would round it to zero, while
    # its reciprocal condition is exactly1 whenever supported controls exist.
    assert not rcond_check_fail(np.ones((3, 1)), [0, 1, 1],
                                i_weights=[np.nextafter(0., 1.), 1., 1.])


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("route", ["panel", "rc"])
def test_no_user_weights_preserve_unweighted_preliminary_guards(monkeypatch, mode, route):
    engine = importlib.import_module(
        "csdid.attgt_fnc.compute_att_gt2" if mode else "csdid.attgt_fnc.compute_att_gt")
    entered = {"overlap_check_fail": 0, "rcond_check_fail": 0}
    for name in entered:
        original = getattr(engine, name)

        def guard(*args, _name=name, _original=original, **kwargs):
            assert kwargs.get("i_weights") is None
            entered[_name] += 1
            return _original(*args, **kwargs)

        monkeypatch.setattr(engine, name, guard)
    data = _panel()
    if route == "rc":
        data["id"] = np.arange(len(data))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                    xformla="~X", faster_mode=mode, panel=route != "rc").fit(bstrap=False)
    assert np.isfinite(fit.MP["att"]).all()
    assert all(count == 3 for count in entered.values())


def test_randomized_rc_support_decisions_match_explicit_group_period_counts():
    rng = np.random.default_rng(278279)
    for _ in range(1000):
        n = rng.integers(4, 80)
        d, post = rng.integers(0, 2, size=(2, n))
        weights = np.exp2(rng.integers(-1074, 1001, size=n).astype(float))
        weights[rng.random(n) < .15] = 0
        n_covariates = int(rng.integers(1, 8))
        normalized = weights / np.mean(weights)
        for include_treated in (False, True):
            expected = any(np.count_nonzero((d == group) & (post == period)
                           & np.isfinite(normalized) & (normalized > 0)) < n_covariates
                           for group in ([0, 1] if include_treated else [0])
                           for period in (0, 1))
            assert rc_outcome_has_too_few_obs(d, post, weights, n_covariates,
                                             include_treated) == expected
