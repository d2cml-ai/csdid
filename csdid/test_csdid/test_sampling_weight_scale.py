"""Sampling weights change mass, not the units of estimates or inference.

Public comparisons use the same seeded data as the independent live R audit.
The cohort-weight oracle differentiates a ratio of empirical totals directly;
it does not copy either implementation of ``wif``.
"""
import contextlib
import io
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from csdid.aggte_fnc.utils import wif, _wif_legacy


SCALES = (1e-200, 1e-160, 1e155, 1e307)
AGGREGATIONS = ("simple", "dynamic", "group", "calendar")


def _two_period_data(panel):
    rng = np.random.default_rng(8201)
    n = 240
    x, z = rng.normal(size=(2, n))
    probability = 1 / (1 + np.exp(.3 - .5 * x + .2 * z - .4 * np.sin(2 * x)))
    treated = (rng.random(n) < probability).astype(float)
    y0 = .8 + .3 * x - .2 * z + rng.normal(size=n)
    y1 = y0 + .6 + .4 * x + 1.1 * treated + rng.normal(size=n)
    weights = np.exp(.2 * x - .15 * z)
    weights[::11] = 0
    data = pd.DataFrame({
        "id": np.tile(np.arange(n), 2), "period": np.repeat([1, 2], n),
        "G": np.tile(2 * treated, 2), "Y": np.r_[y0, y1],
        "X": np.tile(x, 2), "wt": np.tile(weights, 2),
    })
    if not panel:
        data["id"] = np.arange(len(data))
    return data


def _fit(data, panel, method, mode):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return ATTgt(
            yname="Y", tname="period", idname="id", gname="G", data=data,
            panel=panel, weights_name="wt", xformla="~X",
            faster_mode=mode, cband=False,
        ).fit(method, bstrap=False)


def _aggregate(model, kind):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        model.aggte(kind, bstrap=False, cband=False)
    return model.atte


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("scale", SCALES)
def test_public_weight_units_preserve_cell_and_all_aggregation_inference(
        panel, method, mode, scale):
    baseline = _two_period_data(panel)
    scaled = baseline.copy(deep=True)
    scaled["wt"] *= scale
    original = scaled.copy(deep=True)
    expected = _fit(baseline, panel, method, mode)
    actual = _fit(scaled, panel, method, mode)

    for model in (expected, actual):
        assert np.isfinite(model.MP["att"]).all()
        assert np.isfinite(model.results["se"]).all()
        assert np.isfinite(model.MP["inffunc"]["inffunc"]).all()
    np.testing.assert_allclose(actual.MP["att"], expected.MP["att"],
                               atol=2e-10, rtol=2e-10)
    np.testing.assert_allclose(actual.results["se"], expected.results["se"],
                               atol=2e-10, rtol=2e-10)
    np.testing.assert_allclose(actual.MP["inffunc"]["inffunc"],
                               expected.MP["inffunc"]["inffunc"],
                               atol=2e-9, rtol=2e-9)

    # With one estimable cohort-time cell, every aggregation nests that cell.
    # This checks inference independently of the baseline implementation.
    for kind in AGGREGATIONS:
        reference, result = _aggregate(expected, kind), _aggregate(actual, kind)
        for field in ("overall_att", "overall_se"):
            assert np.isfinite(result[field])
            np.testing.assert_allclose(result[field], reference[field],
                                       atol=2e-10, rtol=2e-10)
        assert result["inf_function"].keys() == reference["inf_function"].keys()
        for key, influence in result["inf_function"].items():
            influence = np.asarray(influence).ravel()
            assert np.isfinite(influence).all()
            np.testing.assert_allclose(
                influence, np.asarray(reference["inf_function"][key]).ravel(),
                atol=2e-9, rtol=2e-9)
            np.testing.assert_allclose(
                influence, np.asarray(actual.MP["inffunc"]["inffunc"]).ravel(),
                atol=2e-9, rtol=2e-9)
        np.testing.assert_allclose(result["overall_att"], actual.MP["att"][0],
                                   atol=2e-10, rtol=2e-10)
        np.testing.assert_allclose(result["overall_se"], actual.results["se"][0],
                                   atol=2e-10, rtol=2e-10)

    pd.testing.assert_frame_equal(scaled, original)
    # Internal numeric safeguards must not rewrite the user's original column.
    np.testing.assert_array_equal(np.sort(actual.dp["data"]["wt"]),
                                  np.sort(original["wt"]))


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("weight_units", ["unit", "subnormal"])
def test_sparse_float32_weight_units_preserve_all_public_inference(
        panel, method, mode, weight_units):
    baseline = _two_period_data(panel)
    mask = np.arange(len(baseline)) % 5 == 0
    baseline["wt"] = mask.astype(np.float64)
    scaled = baseline.copy(deep=True)
    scale = (np.nextafter(np.float32(0), np.float32(1))
             if weight_units == "subnormal" else np.float32(1))
    scaled["wt"] = mask.astype(np.float32) * scale
    original = scaled.copy(deep=True)
    assert scaled["wt"].dtype == np.float32
    if weight_units == "subnormal":
        assert np.mean(scaled["wt"].to_numpy()) == 0
    assert np.mean(scaled["wt"].to_numpy(dtype=np.float64)) > 0

    expected = _fit(baseline, panel, method, mode)
    actual = _fit(scaled, panel, method, mode)
    for model in (expected, actual):
        assert np.isfinite(model.MP["att"]).all()
        assert np.isfinite(model.results["se"]).all()
        assert np.isfinite(model.MP["inffunc"]["inffunc"]).all()
    np.testing.assert_allclose(actual.MP["att"], expected.MP["att"],
                               atol=2e-10, rtol=2e-10)
    np.testing.assert_allclose(actual.results["se"], expected.results["se"],
                               atol=2e-10, rtol=2e-10)
    np.testing.assert_allclose(actual.MP["inffunc"]["inffunc"],
                               expected.MP["inffunc"]["inffunc"],
                               atol=2e-9, rtol=2e-9)
    for kind in AGGREGATIONS:
        reference, result = _aggregate(expected, kind), _aggregate(actual, kind)
        for field in ("overall_att", "overall_se", "att_egt", "se_egt"):
            if reference[field] is None:
                assert result[field] is None
                continue
            assert np.isfinite(result[field]).all()
            np.testing.assert_allclose(result[field], reference[field],
                                       atol=2e-10, rtol=2e-10)
        assert result["inf_function"].keys() == reference["inf_function"].keys()
        for key, influence in result["inf_function"].items():
            assert np.isfinite(influence).all()
            np.testing.assert_allclose(
                influence, reference["inf_function"][key], atol=2e-9, rtol=2e-9)
    pd.testing.assert_frame_equal(scaled, original)
    assert actual.dp["data"]["wt"].dtype == np.float32
    np.testing.assert_array_equal(np.sort(actual.dp["data"]["wt"]),
                                  np.sort(original["wt"]))


def _cohort_sample(layout):
    unit_groups = np.array([0, 2, 2, 4, 4, 6, 6, 0])
    weights = np.array([0., 1., 2., 0., 3., 4., 1., 2.])
    if layout == "duplicates":
        groups, keepers = np.array([2, 4, 2, 6]), np.array([3, 0, 2])
    else:
        groups = np.array([2, 4, 6])
        keepers = np.array([1]) if layout == "single" else np.array([2, 0, 1])
    return unit_groups, weights, groups, keepers


def _ratio_of_totals(weights, membership):
    totals = weights @ membership
    return totals / totals.sum()


def _ratio_influence(weights, membership):
    probability = _ratio_of_totals(weights, membership)
    denominator = np.mean(weights * membership.sum(axis=1))
    return (weights[:, None]
            * (membership - membership.sum(axis=1)[:, None] * probability)
            / denominator)


@pytest.mark.parametrize("layout", ["distinct", "duplicates", "single"])
@pytest.mark.parametrize("scale", (1.,) + SCALES)
@pytest.mark.parametrize("function", [wif, _wif_legacy], ids=["optimized", "fallback"])
@pytest.mark.parametrize("pandas_groups", [False, True], ids=["array", "pandas-index"])
def test_cohort_weight_influence_matches_independent_ratio_derivative_in_any_units(
        layout, scale, function, pandas_groups):
    unit_groups, weights, groups, keepers = _cohort_sample(layout)
    membership = (unit_groups[:, None] == groups[keepers]).astype(float)
    expected = _ratio_influence(weights, membership)
    pg = np.array([np.mean(weights * (unit_groups == group)) for group in groups])
    scaled_pg, scaled_weights = pg * scale, weights * scale
    before_pg, before_weights = scaled_pg.copy(), scaled_weights.copy()
    labels = pd.Index(unit_groups) if pandas_groups else unit_groups
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = function(keepers, scaled_pg, scaled_weights, labels, groups)
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=2e-12)
    np.testing.assert_allclose(actual.sum(axis=1), 0., atol=2e-12)
    np.testing.assert_allclose(actual.mean(axis=0), 0., atol=2e-12)
    np.testing.assert_allclose(actual[weights == 0], 0., atol=2e-12)
    np.testing.assert_array_equal(scaled_pg, before_pg)
    np.testing.assert_array_equal(scaled_weights, before_weights)


@pytest.mark.parametrize("layout", ["distinct", "duplicates", "single"])
def test_cohort_ratio_oracle_agrees_with_refitted_empirical_mass_derivatives(layout):
    unit_groups, weights, groups, keepers = _cohort_sample(layout)
    membership = (unit_groups[:, None] == groups[keepers]).astype(float)
    analytic = _ratio_influence(weights, membership)
    h = 1e-4
    derivatives = []
    for row in range(len(weights)):
        plus, minus = weights.copy(), weights.copy()
        plus[row] *= 1 + h
        minus[row] *= 1 - h
        derivatives.append(len(weights) * (
            _ratio_of_totals(plus, membership)
            - _ratio_of_totals(minus, membership)) / (2 * h))
    np.testing.assert_allclose(analytic, derivatives, atol=2e-8, rtol=2e-8)


@pytest.mark.parametrize("layout", ["distinct", "duplicates", "single"])
@pytest.mark.parametrize("scale", [1., 1e-30, 1e30])
@pytest.mark.parametrize("function", [wif, _wif_legacy], ids=["optimized", "fallback"])
@pytest.mark.parametrize("pandas_groups", [False, True], ids=["array", "pandas-index"])
def test_float32_cohort_weight_units_preserve_the_ratio_influence(
        layout, scale, function, pandas_groups):
    unit_groups, weights, groups, keepers = _cohort_sample(layout)
    membership = (unit_groups[:, None] == groups[keepers]).astype(float)
    expected = _ratio_influence(weights, membership)
    weights = (weights * scale).astype(np.float32)
    pg = np.asarray([np.mean(weights * (unit_groups == group))
                     for group in groups], dtype=np.float32)
    labels = pd.Index(unit_groups) if pandas_groups else unit_groups
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = function(keepers, pg, weights, labels, groups)
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, expected, atol=2e-6, rtol=2e-6)


@pytest.mark.parametrize("function", [wif, _wif_legacy], ids=["optimized", "fallback"])
@pytest.mark.parametrize("mass_dtype", [np.float32, np.float64], ids=["native-mass", "double-mass"])
def test_float32_subnormal_cohort_means_preserve_the_ratio_influence(
        function, mass_dtype):
    unit_groups, weights, groups, keepers = _cohort_sample("distinct")
    membership = (unit_groups[:, None] == groups[keepers]).astype(float)
    expected = _ratio_influence(weights, membership)
    smallest = np.nextafter(np.float32(0), np.float32(1))
    tiny_weights = weights.astype(np.float32) * smallest
    pg = np.asarray([np.mean(tiny_weights.astype(mass_dtype)
                             * (unit_groups == group)) for group in groups],
                    dtype=mass_dtype)
    actual = function(keepers, pg, tiny_weights, unit_groups, groups)
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=2e-12)


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("mode", [False, True])
def test_unrelated_large_weight_cohort_does_not_erase_local_positive_support(
        panel, method, mode):
    rng = np.random.default_rng(5172)
    n_group = 60
    units = pd.DataFrame({
        "id": np.arange(3 * n_group), "G": np.repeat([0, 2, 3], n_group),
        "X": rng.normal(size=3 * n_group), "a": rng.normal(size=3 * n_group),
    })
    data = units.merge(pd.DataFrame({"period": [1, 2, 3]}), how="cross")
    data["Y"] = (data.a + .3 * data.period + .1 * data.X * data.period
                 + ((data.G > 0) & (data.period >= data.G))
                 + rng.normal(size=len(data)))
    data["wt"] = np.where(data.G == 3, 1e200, 1e-200)
    if not panel:
        data["id"] = np.arange(len(data))
    # G=3 is excluded from the G=2 vs never-treated comparison. Its units
    # cannot change that comparison's locally positive sampling-weight support.
    reference = _fit(data.loc[data.G != 3].copy(), panel, method, mode)
    actual = _fit(data, panel, method, mode)
    actual_cell = np.flatnonzero((np.asarray(actual.MP["group"]) == 2)
                                & (np.asarray(actual.MP["t"]) == 2))[0]
    reference_cell = np.flatnonzero((np.asarray(reference.MP["group"]) == 2)
                                   & (np.asarray(reference.MP["t"]) == 2))[0]
    assert np.isfinite(actual.MP["att"][actual_cell])
    np.testing.assert_allclose(actual.MP["att"][actual_cell],
                               reference.MP["att"][reference_cell],
                               atol=2e-10, rtol=2e-10)
    full_if = np.asarray(actual.MP["inffunc"]["inffunc"])[:, actual_cell]
    ref_if = np.asarray(reference.MP["inffunc"]["inffunc"])[:, reference_cell]
    np.testing.assert_allclose(full_if[:len(ref_if)] / len(full_if),
                               ref_if / len(ref_if), atol=2e-10, rtol=2e-10)
    np.testing.assert_array_equal(full_if[len(ref_if):], 0.)
