"""Exact parity and behavioral boundaries for request-local panel fit reuse."""
import cProfile
import importlib
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from csdid.attgt_fnc import drdid_trim
from csdid.attgt_fnc.preprocess_did import pre_process_did
from helpers import build_sim_data


fast = importlib.import_module("csdid.attgt_fnc.compute_att_gt2")


def _data(control="nevertreated", fix_weights=None, varying=False,
          formula="~X", anticipation=0):
    data = build_sim_data(n=240, time_periods=5, seed=816)
    if varying:
        data["X"] += .13 * data.period * np.sin(data.id)
    data["weight"] = .3 + (data.id % 13) / 7
    if varying:
        data["weight"] += .02 * data.period * (1 + data.id % 5)
    data["id"] = 31 * data.id - 40
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pre_process_did(
            yname="Y", tname="period", idname="id", gname="G",
            data=data.sample(frac=1, random_state=74), xformla=formula,
            weights_name="weight", control_group=control,
            fix_weights=fix_weights, anticipation=anticipation)


def _run(dp, monkeypatch, method="dr", base="varying", cached=True,
         warning_filter="always"):
    np.random.seed(691)
    profiler = cProfile.Profile()
    with monkeypatch.context() as patch, warnings.catch_warnings(record=True) as caught:
        if not cached:
            patch.setattr(fast, "_CACHEABLE_PANEL_ESTIMATORS", ())
        warnings.simplefilter(warning_filter)
        output, influence = profiler.runcall(
            fast.compute_att_gt2, dp, est_method=method, base_period=base)
    fits = sum(entry.callcount for entry in profiler.getstats()
               if getattr(entry.code, "co_name", None) == "_fit_irls")
    return (output, influence,
            [(warning.category, str(warning.message)) for warning in caught],
            np.random.get_state(), fits)


def _assert_same(expected, actual):
    for key in ("group", "year", "post"):
        np.testing.assert_array_equal(actual[0][key], expected[0][key])
    np.testing.assert_array_equal(actual[0]["att"], expected[0]["att"])
    np.testing.assert_array_equal(actual[1], expected[1])
    assert actual[2] == expected[2]
    for left, right in zip(actual[3], expected[3]):
        np.testing.assert_equal(left, right)


@pytest.mark.parametrize("method", ["dr", "ipw"])
@pytest.mark.parametrize("control", ["nevertreated", "notyettreated"])
@pytest.mark.parametrize("base", ["varying", "universal"])
@pytest.mark.parametrize("fix_weights", [None, "first_period", "base_period"])
@pytest.mark.parametrize("varying", [False, True])
def test_cached_estimates_influence_warnings_and_rng_are_exact(
        monkeypatch, method, control, base, fix_weights, varying):
    dp = _data(control, fix_weights, varying)
    expected = _run(dp, monkeypatch, method, base, cached=False)
    actual = _run(dp, monkeypatch, method, base)
    _assert_same(expected, actual)
    assert actual[4] <= expected[4]


@pytest.mark.parametrize("method", ["dr", "ipw"])
def test_invariant_design_reduces_glm_fits_to_two_per_cohort(monkeypatch, method):
    dp = _data()
    expected = _run(dp, monkeypatch, method, cached=False)
    actual = _run(dp, monkeypatch, method)
    _assert_same(expected, actual)
    assert actual[4] == 2 * len(dp["glist"])
    assert expected[4] == 2 * len(expected[0]["att"])
    assert actual[4] < expected[4]


@pytest.mark.parametrize("anticipation", [0, 1])
def test_time_varying_designs_and_changing_controls_keep_their_own_fits(
        monkeypatch, anticipation):
    dp = _data("notyettreated", varying=True, anticipation=anticipation)
    expected = _run(dp, monkeypatch, base="universal", cached=False)
    actual = _run(dp, monkeypatch, base="universal")
    _assert_same(expected, actual)


@pytest.mark.parametrize("failure", ["overlap", "singular", "both"])
@pytest.mark.parametrize("warning_filter", ["always", "default", "once", "module"])
def test_invalid_cells_preserve_guard_order_and_per_cell_warnings(
        monkeypatch, failure, warning_filter):
    formula = "~X" if failure == "overlap" else "~X+I(2*X)"
    dp = _data(formula=formula)
    if failure in ("overlap", "both"):
        # Every treated cohort separates from never-treated controls. With a
        # duplicate column too, overlap must still be checked before singularity.
        data = dp["data"]
        data[dp["xcov_cols"][1]] = (data.G > 0).astype(float)
        if failure == "both":
            data[dp["xcov_cols"][2]] = 2 * data[dp["xcov_cols"][1]]
    expected = _run(dp, monkeypatch, cached=False, warning_filter=warning_filter)
    actual = _run(dp, monkeypatch, warning_filter=warning_filter)
    _assert_same(expected, actual)
    assert np.isnan(actual[0]["att"]).all()
    assert np.isnan(actual[1]).all()
    assert len(actual[2]) == len(actual[0]["att"])
    message = "singular" if failure == "singular" else "overlap condition"
    assert all(message in text for _, text in actual[2])


@pytest.mark.parametrize("kind", ["callable", "patched_estimator", "patched_glm",
                                  "patched_solver", "patched_overlap", "patched_rcond",
                                  "patched_design"])
def test_replaced_functions_retain_per_cell_execution(monkeypatch, kind):
    dp = _data()
    expected = _run(dp, monkeypatch, cached=False)
    calls = []
    method = "dr"
    if kind in ("callable", "patched_estimator"):
        original = drdid_trim.drdid_panel

        # No private keyword is accepted: existing estimators keep their ABI.
        def estimator(y1, y0, D, covariates=None, i_weights=None):
            calls.append(1)
            return original(y1, y0, D, covariates=covariates, i_weights=i_weights)

        if kind == "callable":
            method = estimator
        else:
            monkeypatch.setattr(drdid_trim, "drdid_panel", estimator)
    elif kind in ("patched_overlap", "patched_rcond", "patched_design"):
        name = {"patched_overlap": "overlap_check_fail", "patched_rcond": "rcond_check_fail",
                "patched_design": "drdid_design_singular"}[kind]
        original = getattr(fast, name)

        def guard(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)

        monkeypatch.setattr(fast, name, guard)
    else:
        field = "fit" if kind == "patched_glm" else "_fit_irls"
        original = getattr(drdid_trim.sm.GLM, field)

        def fit(model, *args, **kwargs):
            calls.append(1)
            return original(model, *args, **kwargs)

        monkeypatch.setattr(drdid_trim.sm.GLM, field, fit)
    actual = _run(dp, monkeypatch, method=method)
    _assert_same(expected, actual)
    cells = len(actual[0]["att"])
    assert len(calls) == (2 * cells if kind in ("patched_glm", "patched_solver") else cells)


@pytest.mark.parametrize("method,name", [("dr", "drdid_panel"), ("ipw", "std_ipw_did_panel")])
def test_estimator_replaced_before_fast_module_import_retains_public_signature(
        monkeypatch, method, name):
    dp = _data()
    expected = _run(dp, monkeypatch, method=method, cached=False)
    original = getattr(drdid_trim, name)
    calls = []

    def estimator(y1, y0, D, covariates=None, i_weights=None):
        calls.append(1)
        return original(y1, y0, D, covariates=covariates, i_weights=i_weights)

    monkeypatch.setattr(drdid_trim, name, estimator)
    spec = importlib.util.spec_from_file_location("_fresh_fast_module", fast.__file__)
    fresh_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fresh_module)
    output, influence = fresh_module.compute_att_gt2(dp, est_method=method)
    np.testing.assert_array_equal(output["att"], expected[0]["att"])
    np.testing.assert_array_equal(influence, expected[1])
    assert len(calls) == len(output["att"])


@pytest.mark.parametrize("name", ["overlap_check_fail", "rcond_check_fail", "drdid_design_singular"])
def test_guard_replaced_before_fast_module_import_runs_for_each_cell(monkeypatch, name):
    dp = _data()
    expected = _run(dp, monkeypatch, cached=False)
    shared = importlib.import_module("csdid.attgt_fnc.compute_att_gt_shared")
    original = getattr(shared, name)
    calls = []

    def guard(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(shared, name, guard)
    spec = importlib.util.spec_from_file_location("_fresh_fast_guard_module", fast.__file__)
    fresh_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fresh_module)
    output, influence = fresh_module.compute_att_gt2(dp)
    np.testing.assert_array_equal(output["att"], expected[0]["att"])
    np.testing.assert_array_equal(influence, expected[1])
    assert len(calls) == len(output["att"])


@pytest.mark.parametrize("warning_filter", ["always", "default", "once", "module", "ignore", "error"])
@pytest.mark.parametrize("fail", [False, True])
def test_weighted_fit_warning_filters_and_exception_order(monkeypatch, warning_filter, fail):
    """Force the cache hook to exercise diagnostics, even with a patched GLM.

    Multivariate overlap guards enter the suppressing context before every fit;
    this also pins the warning-registry behavior of an overlap cache hit.
    """
    calls = []

    class Model:
        def __init__(self, *args, **kwargs):
            pass

        def fit(self):
            calls.append(1)
            warnings.warn("solver diagnostic", RuntimeWarning)
            if fail:
                raise ValueError("fit failure after warning")
            return SimpleNamespace(converged=True, params=np.array([0., 1.]))

    monkeypatch.setattr(drdid_trim.sm, "GLM", Model)

    def run(cached):
        calls.clear()
        cache = {} if cached else None
        errors = []
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter(warning_filter)
            for _ in range(3):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                try:
                    drdid_trim._fit_panel_pscore(
                        np.array([0., 1.]), np.ones((2, 2)), np.ones(2), cache)
                except Exception as error:
                    errors.append((type(error), str(error)))
        return [(w.category, str(w.message), w.filename, w.lineno) for w in caught], errors, len(calls)

    expected, actual = run(False), run(True)
    assert actual[:2] == expected[:2]
    assert actual[2] == (1 if warning_filter == "ignore" and not fail else 3)


@pytest.mark.parametrize("warning_filter", ["always", "default", "once", "module", "ignore", "error"])
def test_intercept_only_retains_original_warning_behavior(monkeypatch, warning_filter):
    dp = _data(formula="~1")
    assert len(dp["xcov_cols"]) == 1
    # Force a warning at the fit helper boundary without changing GLM, so cache
    # eligibility cannot accidentally be disabled by the patched-GLM safeguard.
    original = drdid_trim._fit_panel_pscore
    seen_cache = []

    def fit(D, int_cov, weights, cache):
        seen_cache.append(cache)
        warnings.warn("intercept fit diagnostic", RuntimeWarning)
        return original(D, int_cov, weights, cache)

    monkeypatch.setattr(drdid_trim, "_fit_panel_pscore", fit)
    if warning_filter == "error":
        errors = []
        for cached in (False, True):
            with pytest.raises(UserWarning) as caught:
                _run(dp, monkeypatch, cached=cached, warning_filter=warning_filter)
            errors.append(str(caught.value))
        assert errors[0] == errors[1]
    else:
        expected = _run(dp, monkeypatch, cached=False, warning_filter=warning_filter)
        actual = _run(dp, monkeypatch, warning_filter=warning_filter)
        _assert_same(expected, actual)
    assert all(entry is None for entry in seen_cache)


@pytest.mark.parametrize("converged,params", [(False, [0.]), (True, [np.nan]), (True, [np.inf])])
def test_failed_diagnostics_do_not_enter_propensity_cache(monkeypatch, converged, params):
    calls = []

    class Model:
        def __init__(self, *args, **kwargs):
            pass

        def fit(self):
            calls.append(1)
            return SimpleNamespace(converged=converged, params=np.array(params))

    monkeypatch.setattr(drdid_trim.sm, "GLM", Model)
    cache = {}
    for _ in range(2):
        drdid_trim._fit_panel_pscore(np.array([0., 1.]), np.ones((2, 2)), np.ones(2), cache)
    assert len(calls) == 2
    assert cache == {}


def test_cache_is_request_local_and_does_not_grow_with_periods(monkeypatch):
    dp = _data(varying=True)
    original = fast._att_gt_panel_cell
    cache_ids, observed = [], []

    def cell(*args, **kwargs):
        cache = args[-1]
        result = original(*args, **kwargs)
        if cache is not None:
            observed.append(cache)
            cache_ids.append(id(cache))
            assert set(cache) <= {"cov_keys", "weight_keys", "overlap_key", "overlap_failed",
                                  "pscore_key", "pscore"}
            assert len(cache.get("pscore", {})) <= 1
        return result

    monkeypatch.setattr(fast, "_att_gt_panel_cell", cell)
    first = _run(dp, monkeypatch)
    first_ids = set(cache_ids)
    cache_ids.clear()
    second = _run(dp, monkeypatch)
    _assert_same(first, second)
    assert len(first_ids) == len(dp["glist"])
    assert not first_ids.intersection(cache_ids)
