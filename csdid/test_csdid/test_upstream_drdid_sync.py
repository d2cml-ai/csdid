"""Regression tests for current DRDID's influence functions and fitted designs.

DRDID #20's old sign survives intercept-only and mean-zero IF tests. These use
a misspecified propensity model and covariates, and differentiate the actual
refitted ATT rather than reproducing its influence-function algebra.
The design checks also retain feasible regression-only RC estimates when a
treated cell is too small to support the extra outcome regression used by DR.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from csdid.attgt_fnc.compute_att_gt_shared import drdid_design_singular
from csdid.attgt_fnc.drdid_trim import drdid_rc


def _rc_data(weight_kind):
    rng = np.random.default_rng(20261021)
    n = 320
    x = rng.normal(size=n)
    z = rng.normal(size=n)
    # The fitted logit includes x and z, but omits the nonlinear term.
    p = 1 / (1 + np.exp(0.5 - 1.5 * np.maximum(x, 0) ** 2 + 0.3 * z))
    d = (rng.uniform(size=n) < p).astype(float)
    post = (rng.uniform(size=n) < 0.5).astype(float)
    y = (1 + 2 * x - z + 0.5 * d + post * (0.5 + x + d)
         + rng.normal(size=n) * np.where((d == 0) & (post == 0), 3, 1))
    cov = np.column_stack((np.ones(n), x, z))
    weights = np.ones(n)
    if weight_kind != "unit":
        weights = rng.uniform(0.2, 2.5, size=n)
    if weight_kind == "zeros":
        weights[::13] = 0
    return y, post, d, cov, weights


@pytest.mark.parametrize("weight_kind", ["unit", "unequal", "zeros"])
@pytest.mark.parametrize("trim_level", [0.995, 0.75])
def test_rc_if_matches_refitted_empirical_derivatives(weight_kind, trim_level):
    y, post, d, cov, weights = _rc_data(weight_kind)
    _, influence = drdid_rc(y, post, d, cov, weights, trim_level)
    n = len(y)
    # Cover all four cells, including high-leverage observations.  The defect
    # affects the pre-period controls through their estimated OR coefficients.
    indices = []
    for dd in (0, 1):
        for tt in (0, 1):
            cell = np.flatnonzero((d == dd) & (post == tt))
            indices.extend(cell[np.argsort(np.abs(influence[cell]))[-3:]])
    indices.extend(np.flatnonzero(weights == 0)[:2])
    h = 0.001
    numerical = []
    for i in indices:
        plus, minus = weights.copy(), weights.copy()
        plus[i] *= 1 + h
        minus[i] *= 1 - h
        att_plus = drdid_rc(y, post, d, cov, plus, trim_level)[0]
        att_minus = drdid_rc(y, post, d, cov, minus, trim_level)[0]
        # Multiplying row i's mass by 1+h moves the empirical law by h/(n+h).
        eps_plus, eps_minus = h / (n + h), -h / (n - h)
        numerical.append((att_plus - att_minus) / (eps_plus - eps_minus))
    np.testing.assert_allclose(influence[indices], numerical, atol=2e-5, rtol=2e-6)


@pytest.mark.parametrize("weight_kind", ["unit", "unequal", "zeros"])
def test_rc_period_reversal_negates_att_and_influence(weight_kind):
    y, post, d, cov, weights = _rc_data(weight_kind)
    att, influence = drdid_rc(y, post, d, cov, weights)
    reversed_att, reversed_influence = drdid_rc(y, 1 - post, d, cov, weights)
    np.testing.assert_allclose(reversed_att, -att, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(reversed_influence, -influence, atol=1e-10, rtol=1e-10)


def _tiny_treated_rc():
    rng = np.random.default_rng(107)
    post = np.r_[np.zeros(40), np.ones(40), 0, 1]
    d = np.r_[np.zeros(80), 1, 1]
    x = np.r_[rng.normal(size=80), 0.4, 0.4]
    y = (1 + 2 * x + post * (0.5 + x) + 2 * d * post
         + rng.normal(scale=0.1, size=82))
    return pd.DataFrame({
        "id": np.arange(82), "period": post + 1, "G": 2 * d, "X": x,
        "Y": y, "wt": rng.uniform(0.5, 2, size=82),
    })


@pytest.mark.parametrize("method", ["reg", "dr"])
@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("singular_controls", [False, True])
def test_rc_guard_checks_only_the_outcome_designs_actually_fitted(
        method, mode, singular_controls):
    data = _tiny_treated_rc()
    if singular_controls:
        # The pooled control design is full rank, but its pre-period cell is not.
        data.loc[(data.G == 0) & (data.period == 1), "X"] = 0
    cov = np.column_stack((np.ones(len(data)), data.X))
    expected_failure = singular_controls or method == "dr"
    assert drdid_design_singular(
        cov, data.wt, data.G == 2, method, post=data.period == 2,
    ) == expected_failure
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = ATTgt(
            yname="Y", tname="period", idname="id", gname="G", data=data,
            xformla="Y ~ X", weights_name="wt", panel=False, faster_mode=mode,
        ).fit(est_method=method, bstrap=False)
    if expected_failure:
        assert np.isnan(fit.MP["att"]).all()
        assert np.isnan(fit.results["se"]).all()
        assert any("singular" in str(w.message) for w in caught)
    else:
        # Executed with source DRDID 1.4.0 eaac3b4 and did 74b88eb, both modes:
        # did::att_gt(..., panel=FALSE, est_method="reg", xformla=~X,
        #             weightsname="wt", bstrap=FALSE). One treated row per cell
        # cannot support a treated OR, which this estimator never fits.
        np.testing.assert_allclose(fit.MP["att"], [2.04760499329232], atol=1e-10, rtol=0)
        np.testing.assert_allclose(fit.results["se"], [0.0209179458588539], atol=1e-10, rtol=0)
        assert not any("singular" in str(w.message) for w in caught)
