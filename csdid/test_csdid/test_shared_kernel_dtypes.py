"""Shared estimators use the represented numeric values in double precision.

The large-integer oracle is a weighted-free difference of group means, computed
after conversion to doubles. It does not invoke another estimator's equations.
"""
import contextlib
import io
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from csdid.attgt_fnc import drdid_trim


KERNELS = ("drdid_panel", "std_ipw_did_panel", "drdid_rc", "std_ipw_did_rc")
AGGREGATIONS = ("simple", "dynamic", "group", "calendar")


def _large_integer_outcomes():
    rng = np.random.default_rng(74129)
    n = 120
    treated = np.repeat([0, 1], n // 2)
    y0 = (np.where(treated == 0, 10**18, -9 * 10**18).astype(np.int64)
          + rng.integers(-9, 9, n, dtype=np.int64) * 10**14)
    y1 = (np.full(n, 9 * 10**18, dtype=np.int64)
          + rng.integers(-9, 9, n, dtype=np.int64) * 10**14)
    # Every input is representable; subtraction in the input dtype is not.
    double_change = y1.astype(float) - y0.astype(float)
    assert np.any((y1 - y0).astype(float) != double_change)
    return y1, y0, treated


def _kernel_inputs(name, dtype):
    rng = np.random.default_rng(8201)
    if dtype == np.int64:
        y1, y0, treated = _large_integer_outcomes()
        n = len(treated)
        cov = np.column_stack((np.ones(n), rng.integers(-3, 4, n)))
        weights = rng.integers(1, 5, n)
    else:
        n = 240
        x, z = rng.normal(size=(2, n))
        cov = np.column_stack((np.ones(n), x, z))
        probability = 1 / (1 + np.exp(.3 - .5 * x + .2 * z))
        treated = (rng.random(n) < probability).astype(float)
        y0 = .8 + .3 * x - .2 * z + rng.normal(size=n)
        y1 = y0 + .6 + .4 * x + 1.1 * treated + rng.normal(size=n)
        weights = np.exp(.2 * x - .15 * z)
    post = (rng.random(n) < .55).astype(int)
    outcomes = (np.where(post, y1, y0), post) if name.endswith("_rc") else (y1, y0)
    return (tuple(np.asarray(v, dtype=dtype) for v in (*outcomes, treated)),
            np.asarray(cov, dtype=dtype), np.asarray(weights, dtype=dtype))


@pytest.mark.parametrize("name", KERNELS)
@pytest.mark.parametrize("dtype", [np.float32, np.int64], ids=["float32", "large-int64"])
def test_local_kernels_use_identical_represented_values_in_double_precision(name, dtype):
    args, cov, weights = _kernel_inputs(name, dtype)
    original = tuple(values.copy() for values in (*args, cov, weights))
    kernel = getattr(drdid_trim, name)
    expected = kernel(*(values.astype(float) for values in args),
                      covariates=cov.astype(float), i_weights=weights.astype(float))
    actual = kernel(*args, covariates=cov, i_weights=weights)
    assert np.isfinite(actual[0]) and np.isfinite(actual[1]).all()
    np.testing.assert_allclose(actual[0], expected[0], atol=1e-7, rtol=1e-6)
    np.testing.assert_allclose(actual[1], expected[1], atol=1e-7, rtol=1e-6)
    for values, before in zip((*args, cov, weights), original):
        np.testing.assert_array_equal(values, before)


def _integer_frame(panel):
    y1, y0, treated = _large_integer_outcomes()
    n = len(treated)
    data = pd.DataFrame({"id": np.tile(np.arange(n), 2),
                         "period": np.repeat([1, 2], n),
                         "G": np.tile(2 * treated, 2), "Y": np.r_[y0, y1]})
    if not panel:
        data["id"] = np.arange(len(data))
    return data


def _closed_form(data, panel):
    outcome = data.Y.to_numpy(dtype=float)
    treated = data.G.to_numpy() == 2
    if panel:
        n = len(data) // 2
        change, treated = outcome[n:] - outcome[:n], treated[:n]
        mean_treated, mean_control = np.mean(change[treated]), np.mean(change[~treated])
        influence = (treated * (change - mean_treated) / np.mean(treated)
                     - (~treated) * (change - mean_control) / np.mean(~treated))
        return mean_treated - mean_control, influence
    post = data.period.to_numpy() == 2
    influence = np.zeros(len(data))
    att = 0.
    for d, t, sign in ((True, True, 1), (True, False, -1),
                       (False, True, -1), (False, False, 1)):
        cell = (treated == d) & (post == t)
        mean = np.mean(outcome[cell])
        att += sign * mean
        influence += sign * cell * (outcome - mean) / np.mean(cell)
    return att, influence


def _fit(data, panel, method, mode, covariates=False):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                     panel=panel, faster_mode=mode, cband=False,
                     weights_name="wt" if "wt" in data else None,
                     xformla="~X" if covariates else "~1").fit(method, bstrap=False)


def _aggregate(model, kind):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        model.aggte(kind, bstrap=False, cband=False)
    return model.atte


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("mode", [False, True])
def test_public_large_integer_outcomes_match_independent_difference_of_means(
        panel, method, mode):
    data = _integer_frame(panel)
    original = data.copy(deep=True)
    expected_att, expected_if = _closed_form(data, panel)
    model = _fit(data, panel, method, mode)
    actual_if = np.asarray(model.MP["inffunc"]["inffunc"]).ravel()
    expected_se = np.sqrt(np.mean(expected_if**2) / len(expected_if))
    assert np.isfinite(model.MP["att"]).all() and np.isfinite(actual_if).all()
    np.testing.assert_allclose(model.MP["att"], [expected_att], atol=1e-7, rtol=1e-6)
    np.testing.assert_allclose(actual_if, expected_if, atol=1e-7, rtol=1e-6)
    np.testing.assert_allclose(model.results["se"], [expected_se], atol=1e-7, rtol=1e-6)
    for kind in AGGREGATIONS:
        result = _aggregate(model, kind)
        np.testing.assert_allclose(result["overall_att"], expected_att, atol=1e-7, rtol=1e-6)
        np.testing.assert_allclose(result["overall_se"], expected_se, atol=1e-7, rtol=1e-6)
        for influence in result["inf_function"].values():
            np.testing.assert_allclose(np.asarray(influence).ravel(), expected_if,
                                       atol=1e-7, rtol=1e-6)
    pd.testing.assert_frame_equal(data, original)


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("mode", [False, True])
def test_public_float32_inputs_match_their_identical_double_representations(panel, method, mode):
    args, cov, weights = _kernel_inputs("drdid_panel", np.float32)
    y1, y0, treated = args
    n = len(treated)
    data = pd.DataFrame({"id": np.tile(np.arange(n), 2),
                         "period": np.repeat([1, 2], n),
                         "G": np.tile(2 * treated, 2), "Y": np.r_[y0, y1],
                         "X": np.tile(cov[:, 1], 2), "wt": np.tile(weights, 2)})
    if not panel:
        data["id"] = np.arange(len(data))
    original = data.copy(deep=True)
    reference = data.astype({"Y": float, "X": float, "wt": float})
    expected, actual = _fit(reference, panel, method, mode, True), _fit(data, panel, method, mode, True)
    for field in ("att", "se", "l_se", "u_se"):
        np.testing.assert_allclose(actual.results[field], expected.results[field], atol=1e-7, rtol=1e-6)
    np.testing.assert_allclose(actual.MP["inffunc"]["inffunc"],
                               expected.MP["inffunc"]["inffunc"], atol=1e-7, rtol=1e-6)
    for kind in AGGREGATIONS:
        result, baseline = _aggregate(actual, kind), _aggregate(expected, kind)
        for field in ("overall_att", "overall_se"):
            np.testing.assert_allclose(result[field], baseline[field], atol=1e-7, rtol=1e-6)
        for key, influence in result["inf_function"].items():
            np.testing.assert_allclose(influence, baseline["inf_function"][key], atol=1e-7, rtol=1e-6)
    pd.testing.assert_frame_equal(data, original)


@pytest.mark.parametrize("panel,mode", [(True, True), (False, False), (False, True)])
def test_custom_estimators_keep_original_outcome_and_weight_dtypes(panel, mode):
    data = _integer_frame(panel)
    data["wt"] = (np.arange(len(data)) % 7).astype(np.float32)
    original = data.copy(deep=True)
    observed = []

    def custom(*args, **kwargs):
        outcome = args[0] if args else kwargs["y"]
        weights = kwargs["i_weights"]
        observed.append((np.asarray(outcome).dtype, np.asarray(weights).dtype))
        return 1., np.linspace(-1, 1, len(outcome))

    _fit(data, panel, custom, mode)
    assert observed and all(types == (np.dtype(np.int64), np.dtype(np.float32))
                            for types in observed)
    pd.testing.assert_frame_equal(data, original)
