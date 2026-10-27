"""Final-score inference contracts for the four owned DR/IPW estimators.

Current R DRDID (eaac3b49) rebuilds propensity information at its final capped
scores. Statsmodels IRLS's stored covariance can retain preceding-iteration
weights. Full current-R numerical fixtures independently gate all IF entries;
these tests also differentiate the fitted estimator itself and guard against
reintroducing dependence on that incidental covariance state.
"""
import numpy as np
import pytest
import statsmodels.api as sm

from csdid.attgt_fnc import drdid_trim


METHODS = ("drdid_panel", "std_ipw_did_panel", "drdid_rc", "std_ipw_did_rc")


def _sample(name):
    rng = np.random.default_rng(20261024)
    n = 240
    x, z = rng.normal(size=(2, n))
    covariates = np.column_stack((np.ones(n), x, z))
    probability = 1 / (1 + np.exp(.2 - .35*x + .15*z))
    treated = (rng.random(n) < probability).astype(float)
    post = (rng.random(n) < .5).astype(float)
    y0 = .3*x - .2*z + .4*x*z + rng.normal(size=n)
    y1 = y0 + .2 + .3*x + .2*z*z + treated*(1 + .1*x) + rng.normal(size=n)
    weights = np.exp(.2*x - .1*z)
    outcomes = (post*y1 + (1-post)*y0, post) if name.endswith("_rc") else (y1, y0)
    return (*outcomes, treated), covariates, weights


@pytest.mark.parametrize("name", METHODS)
def test_full_influence_differentiates_the_refitted_estimator(name):
    """A numerical derivative supplies an independent influence-function oracle.

    Perturb original sampling mass, refit all nuisance models, and compare the
    derivative of ATT to the empirical IF projection. This exercises every DR
    nuisance contribution as well as normalized IPW, without coding a Hessian
    formula in the test.
    """
    args, covariates, weights = _sample(name)
    estimator = getattr(drdid_trim, name)
    _, influence = estimator(*args, covariates=covariates, i_weights=weights)
    rng = np.random.default_rng(1037)
    directions = [rng.normal(size=len(weights)), np.sin(covariates[:, 1]),
                  np.cos(covariates[:, 2])]
    step = 1e-4
    for direction in directions:
        plus = estimator(*args, covariates=covariates,
                         i_weights=weights*(1 + step*direction))[0]
        minus = estimator(*args, covariates=covariates,
                          i_weights=weights*(1 - step*direction))[0]
        derivative = (plus - minus) / (2*step)
        np.testing.assert_allclose(np.mean(influence*direction), derivative,
                                   atol=1e-7, rtol=1e-6)


@pytest.mark.parametrize("name", METHODS)
def test_inference_ignores_stale_solver_covariance_without_changing_fit(name, monkeypatch):
    """Changing only GLM's stored covariance cannot alter this estimator's IF."""
    args, covariates, weights = _sample(name)
    estimator = getattr(drdid_trim, name)
    expected_att, expected_if = estimator(*args, covariates=covariates, i_weights=weights)
    original_fit = sm.GLM.fit
    fits = []

    def fit(model, *fit_args, **fit_kwargs):
        assert fit_args == () and fit_kwargs == {}, "GLM defaults were changed"
        result = original_fit(model, *fit_args, **fit_kwargs)
        fits.append(result)

        def stale_covariance(*args, **kwargs):
            raise AssertionError("Inference must use final scores, not stored IRLS covariance")

        result.cov_params = stale_covariance
        return result

    monkeypatch.setattr(sm.GLM, "fit", fit)
    actual_att, actual_if = estimator(*args, covariates=covariates, i_weights=weights)
    assert len(fits) == 1
    assert actual_att == expected_att
    np.testing.assert_array_equal(actual_if, expected_if)


@pytest.mark.parametrize("name", ["drdid_panel", "std_ipw_did_panel"])
def test_final_score_inference_preserves_clean_propensity_reuse(name, monkeypatch):
    args, covariates, weights = _sample(name)
    estimator = getattr(drdid_trim, name)
    original_fit = sm.GLM.fit
    calls = []

    def fit(model, *fit_args, **fit_kwargs):
        calls.append(1)
        return original_fit(model, *fit_args, **fit_kwargs)

    monkeypatch.setattr(sm.GLM, "fit", fit)
    cache = {}
    first = estimator(*args, covariates=covariates, i_weights=weights, _pscore_cache=cache)
    second = estimator(*args, covariates=covariates, i_weights=weights, _pscore_cache=cache)
    assert len(calls) == 1 and "fit" in cache
    assert first[0] == second[0]
    np.testing.assert_array_equal(first[1], second[1])


@pytest.mark.parametrize("name", ["std_ipw_did_panel", "std_ipw_did_rc"])
def test_singular_propensity_information_is_not_a_pseudoinverse_estimate(name):
    args, covariates, weights = _sample(name)
    duplicated = np.column_stack((covariates, covariates[:, 1]))
    with pytest.raises(ValueError, match="propensity score design matrix is singular"):
        getattr(drdid_trim, name)(*args, covariates=duplicated, i_weights=weights)
