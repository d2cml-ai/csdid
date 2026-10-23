"""Adversarial equivalence checks for aggregation allocation and metadata reuse."""
import contextlib
import copy
import importlib
import io
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.aggte_fnc.utils import get_se, wif, _make_cached_get_se
from csdid.att_gt import ATTgt
from helpers import build_sim_data


aggregation = importlib.import_module("csdid.aggte_fnc.compute_aggte")
aggregation_utils = importlib.import_module("csdid.aggte_fnc.utils")


def _reference_wif(keepers, pg, weights_ind, G, group):
    # Original two-matrix implementation, including its casting/reduction order.
    pg, group = np.array(pg), np.array(group)
    first = np.empty((len(weights_ind), len(keepers)))
    for i, k in enumerate(keepers):
        numerator = weights_ind * 1 * (G == group[k]) - pg[k]
        first[:, i] = (numerator[:, None] / np.sum(pg[keepers])).squeeze()
    second = np.empty((len(weights_ind), len(keepers)))
    for i, k in enumerate(keepers):
        numerator = weights_ind * 1 * (G == group[k]) - pg[k]
        second[:, i] = numerator.squeeze()
    second = np.sum(second, axis=1)
    multiplier = pg[keepers] / sum(pg[keepers]) ** 2
    return first - np.outer(second, multiplier)


def _invoke(mp, options):
    np.random.seed(1926)
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = aggregation.compute_aggte(copy.deepcopy(mp), **options)
        except Exception as exc:
            result = (type(exc), str(exc))
    return result, [(type(w.message), str(w.message)) for w in caught], np.random.get_state()


def _assert_same(expected, actual):
    a, wa, sa = expected
    b, wb, sb = actual
    assert wa == wb
    for old, new in zip(sa, sb):
        np.testing.assert_array_equal(old, new)
    if isinstance(a, tuple):
        assert a == b
        return
    assert not isinstance(b, tuple), b
    assert a.keys() == b.keys()
    for key in a:
        if key == "DIDparams":
            assert a[key].keys() == b[key].keys()
        elif isinstance(a[key], dict):
            for child in a[key]:
                np.testing.assert_array_equal(a[key][child], b[key][child])
        elif a[key] is None or isinstance(a[key], str):
            assert a[key] == b[key]
        else:
            np.testing.assert_array_equal(a[key], b[key])


@pytest.fixture(scope="module", params=[(False, False), (False, True), (True, False), (True, True)])
def fitted_mp(request):
    faster_mode, panel = request.param
    data = build_sim_data(n=180, time_periods=5, seed=813)
    data["weight"] = 0.2 + (data.id % 17) / 4
    if not panel:
        data["id"] = np.arange(len(data))
    fit = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                panel=panel, faster_mode=faster_mode, xformla="Y~X",
                weights_name="weight", clustervar="cluster", cband=False, biters=17)
    with contextlib.redirect_stdout(io.StringIO()):
        fit.fit(bstrap=False)
    return fit.MP


@pytest.mark.parametrize("kind", ["simple", "group", "dynamic", "calendar"])
@pytest.mark.parametrize("cluster", [None, "cluster"])
@pytest.mark.parametrize("inference", ["analytic", "bootstrap", "analytic_cband"])
@pytest.mark.parametrize("case", ["full", "window", "missing"])
def test_aggregation_exact_equivalence(fitted_mp, monkeypatch, kind, cluster, inference, case):
    mp = copy.deepcopy(fitted_mp)
    mp["DIDparams"]["clustervars"] = cluster
    options = dict(typec=kind, bstrap=inference == "bootstrap",
                   cband=inference != "analytic", biters=17)
    if case == "window":
        options.update(min_e=0, max_e=2, balance_e=2)
    elif case == "missing":
        mp["att"][0] = np.nan
        options["na_rm"] = True
    actual = _invoke(mp, options)
    with monkeypatch.context() as patch:
        patch.setattr(aggregation, "wif", _reference_wif)
        patch.setattr(aggregation, "_make_cached_get_se", lambda helper: helper)
        expected = _invoke(mp, options)
    _assert_same(expected, actual)


@pytest.mark.parametrize("dtype", [np.float64, np.float32, np.int64])
@pytest.mark.parametrize("n,k", [(1, 1), (9, 5), (50, 20)])
@pytest.mark.parametrize("selection", ["empty", "reverse", "duplicates"])
def test_wif_preserves_dtype_empty_and_duplicate_semantics(dtype, n, k, selection):
    rng = np.random.default_rng(188)
    unit_groups = rng.choice([None, 1, 2, 3], n).tolist()
    groups = rng.choice([1, 2, 3], k)
    weights = rng.integers(0, 5, n).astype(dtype)
    pg = np.asarray([np.mean(weights * (np.asarray(unit_groups) == g))
                     for g in groups], dtype=dtype)
    keepers = {"empty": [], "reverse": list(range(k))[::-1],
               "duplicates": rng.choice(k, k).tolist()}[selection]
    with warnings.catch_warnings(record=True) as old_warnings:
        warnings.simplefilter("always")
        expected = _reference_wif(keepers, pg, weights, unit_groups, groups)
    with warnings.catch_warnings(record=True) as new_warnings:
        warnings.simplefilter("always")
        actual = wif(keepers, pg, weights, unit_groups, groups)
    np.testing.assert_array_equal(actual, expected)
    assert [str(w.message) for w in old_warnings] == [str(w.message) for w in new_warnings]


@pytest.mark.parametrize("pg,weights", [
    ([0., 0.], [1., 2., 3.]),
    ([np.inf, -np.inf], [1., 2., 3.]),
    ([.2, .3], [np.inf, 2., np.nan]),
    ([1e308, 1e308], [1., 2., 3.]),
    ([1e-300, 1e-300], [1., 2., 3.]),
    ([-1e308, 2e307], [1e308, 2., 3.]),
])
def test_wif_retains_exceptional_arithmetic_warnings(pg, weights):
    args = ([1, 0], np.array(pg), np.array(weights), [1, 2, None], [1, 2])
    with warnings.catch_warnings(record=True) as old_warnings:
        warnings.simplefilter("always")
        expected = _reference_wif(*args)
    with warnings.catch_warnings(record=True) as new_warnings:
        warnings.simplefilter("always")
        actual = wif(*args)
    np.testing.assert_array_equal(actual, expected)
    assert [str(w.message) for w in old_warnings] == [str(w.message) for w in new_warnings]


@pytest.mark.parametrize("unit_groups,groups", [(["a", "b"], ["a"]),
                                              ([None, None], [None])])
def test_wif_preserves_legacy_non_numeric_cohort_validation(unit_groups, groups):
    with pytest.raises(ValueError, match="cond should be a logical vector"):
        wif([0], [.5], np.ones(2), unit_groups, groups)


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_cached_cluster_se_matches_public_helper_including_order(panel, dtype):
    data = pd.DataFrame({"id": np.repeat([8, 2, 4, 1, 3, 7], 2),
                         "time": np.tile([2004, 2000], 6),
                         "cluster": np.repeat([2., 1., np.nan, 2., 3., 1.], 2)})
    dp = {"data": data, "idname": "id", "tname": "time", "panel": panel,
          "clustervars": "cluster", "bstrap": False, "alp": .05, "cband": False}
    cached = _make_cached_get_se(get_se)
    for influence in ([1e12, 3., -1e12, 8., 1., -2.], [0., 0., 0., 0., 0., 0.],
                      [np.nan, 0., 4., 8., 2., 1.]):
        influence = np.asarray(influence, dtype=dtype)[:, None]
        np.testing.assert_array_equal(cached(influence, dp), get_se(influence, dp))
    # A mismatched IF length must retain the public helper's refusal.
    with pytest.raises(IndexError):
        cached(np.zeros(7), dp)


def test_cluster_metadata_rebuilt_after_in_place_data_mutation(fitted_mp, monkeypatch):
    mp = copy.deepcopy(fitted_mp)
    dp = mp["DIDparams"]
    original_keys = set(dp)
    calls = []
    set_index = pd.DataFrame.set_index

    def spy(frame, *args, **kwargs):
        if list(frame.columns) == [dp["idname"], "cluster"]:
            calls.append(1)
        return set_index(frame, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "set_index", spy)
    with contextlib.redirect_stdout(io.StringIO()):
        aggregation.compute_aggte(mp, typec="dynamic", bstrap=False, cband=False)
        assert len(calls) == 1
        dp["data"].loc[:, "cluster"] = 1
        actual = aggregation.compute_aggte(mp, typec="dynamic", bstrap=False, cband=False)
        assert len(calls) == 2
    assert set(dp) == original_keys
    with monkeypatch.context() as patch:
        patch.setattr(aggregation, "_make_cached_get_se", lambda helper: helper)
        with contextlib.redirect_stdout(io.StringIO()):
            expected = aggregation.compute_aggte(mp, typec="dynamic", bstrap=False, cband=False)
    for key in ("se_egt", "overall_se"):
        np.testing.assert_array_equal(actual[key], expected[key])


@pytest.mark.parametrize("override_both_modules", [False, True])
def test_overridden_public_se_helper_still_receives_every_call(fitted_mp, monkeypatch, override_both_modules):
    calls = []

    def overridden(influence, DIDparams=None):
        calls.append(influence.shape)
        return np.nan

    monkeypatch.setattr(aggregation, "get_se", overridden)
    if override_both_modules:
        monkeypatch.setattr(aggregation_utils, "get_se", overridden)
    with contextlib.redirect_stdout(io.StringIO()):
        result = aggregation.compute_aggte(copy.deepcopy(fitted_mp), typec="dynamic",
                                           bstrap=False, cband=False)
    assert len(calls) == len(result["egt"]) + 1
    assert np.isnan(result["overall_se"])
    assert np.all(np.isnan(result["se_egt"]))


def test_mixed_integer_float_time_mapping_retains_numpy_equality():
    start = 2**53
    cohorts = np.array([0, 0, start + 1, start + 1, start + 5, start + 5])
    years = np.array([start, start + 4, start + 8], dtype=float)
    data = pd.DataFrame({"id": np.repeat(np.arange(6), 3), "time": np.tile(years, 6),
                         "cohort": np.repeat(cohorts, 3), "w1": np.ones(18)})
    groups = np.repeat(np.array([start + 1, start + 5]), 2)
    times = np.tile(years[1:], 2)
    influence = np.random.default_rng(83).normal(size=(6, 4))
    att = np.array([2., 8., 3., 5.])
    mp = {"group": groups, "t": times, "att": att, "n": 6,
          "inffunc": {"inffunc": influence}, "DIDparams": {
              "data": data, "idname": "id", "tname": "time", "gname": "cohort",
              "panel": True, "tlist": years, "glist": np.unique(groups),
              "bstrap": False, "cband": False, "clustervars": None, "biters": 17,
              "alp": .05, "control_group": "nevertreated", "anticipation": 0,
              "est_method": "dr"}}
    mass = (cohorts[:, None] == groups).astype(float)
    pg = mass.mean(axis=0)
    estimate = np.sum(att * pg) / np.sum(pg)
    expected_if = influence @ (pg / pg.sum()) + mass @ (att - estimate) / pg.sum()
    with contextlib.redirect_stdout(io.StringIO()):
        result = aggregation.compute_aggte(mp, typec="simple")
    np.testing.assert_allclose(result["overall_att"], estimate)
    np.testing.assert_allclose(result["inf_function"]["simple_att"].ravel(), expected_if)
    np.testing.assert_allclose(result["overall_se"], np.linalg.norm(expected_if) / 6)
