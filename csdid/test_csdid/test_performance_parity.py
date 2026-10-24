"""Full-output parity and public-default checks for the 0.4.3 optimizations."""
import contextlib
import inspect
import io
import warnings

import numpy as np
import pytest

from csdid.att_gt import ATTgt
from benchmark_v043 import AGG_FIELDS, FIT_FIELDS, compare
from helpers import build_sim_data


def _data(design):
    data = build_sim_data(n=400, time_periods=5, seed=9813)
    # Shuffled, non-dense identifiers, factor covariates, unequal time-varying
    # weights and covariates, and original calendar units exercise the caches.
    data["factor"] = (data.id % 3).astype(str)
    data["X"] += .04 * data.period * np.sin(data.id)
    data["wt"] = np.exp(.2 * data.X + .015 * data.period * (data.id % 2))
    data["cluster"] = data.id % 23
    data["id"] = 100 + 13 * data.id
    if design == "unbalanced":
        rng = np.random.default_rng(741)
        data = data.loc[rng.random(len(data)) > .12].copy()
    elif design == "rc":
        data["id"] = 7000 + 17 * np.arange(len(data))
    clock = {0: 0, 1: 2000, 2: 2002, 3: 2005, 4: 2009, 5: 2014}
    data["period"] = data.period.map(clock)
    data["G"] = data.G.map(clock)
    return data.sample(frac=1, random_state=14).reset_index(drop=True)


def _run(design, control, method="dr", base="varying", mode=None,
         bootstrap=False, clustered=True):
    data = _data(design)
    options = dict(yname="Y", tname="period", idname="id", gname="G", data=data,
                   xformla="Y~X+C(factor)", weights_name="wt",
                   control_group=control, panel=design != "rc",
                   allow_unbalanced_panel=design == "unbalanced",
                   clustervar="cluster" if clustered else None, cband=bootstrap, biters=79)
    if mode is not None:
        options["faster_mode"] = mode
    np.random.seed(7982)
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        obj = ATTgt(**options).fit(method, base_period=base, bstrap=bootstrap)
        fit = {field: obj.results[field] for field in FIT_FIELDS}
        fit.update(group=obj.MP["group"], t=obj.MP["t"], W=obj.MP["W"],
                   Wpval=obj.MP["Wpval"], inffunc=obj.MP["inffunc"]["inffunc"],
                   cband=obj.dp["cband"])
        aggs = {}
        for kind in ("simple", "group", "dynamic", "calendar"):
            obj.aggte(typec=kind, bstrap=bootstrap, cband=bootstrap)
            aggs[kind] = {field: obj.atte.get(field) for field in AGG_FIELDS}
            aggs[kind]["inf_function"] = obj.atte["inf_function"]
    return dict(fit=fit, aggs=aggs,
                warnings=[(w.category.__name__, str(w.message)) for w in caught],
                rng_state=np.random.get_state())


def _assert_parity(expected, actual, exact=False):
    result = compare(expected, actual, atol=0. if exact else 1e-9,
                     rtol=0. if exact else 1e-9)
    assert not result["issues"], result["issues"][:10]
    # RNG state must match exactly even for the numerically tolerant mode check.
    for left, right in zip(expected["rng_state"], actual["rng_state"]):
        np.testing.assert_equal(left, right)


@pytest.mark.parametrize("design", ["panel", "unbalanced", "rc"])
@pytest.mark.parametrize("control", ["nevertreated", "notyettreated"])
def test_default_uses_fast_mode_with_complete_outputs(design, control):
    assert inspect.signature(ATTgt).parameters["faster_mode"].default is True
    default = _run(design, control)
    _assert_parity(_run(design, control, mode=True), default, exact=True)
    _assert_parity(_run(design, control, mode=False), default)


@pytest.mark.parametrize("design", ["panel", "unbalanced", "rc"])
@pytest.mark.parametrize("control", ["nevertreated", "notyettreated"])
@pytest.mark.parametrize("method", ["dr", "reg", "ipw"])
@pytest.mark.parametrize("base", ["varying", "universal"])
def test_both_modes_preserve_full_influence_and_inference(design, control, method, base):
    _assert_parity(_run(design, control, method, base, mode=False),
                   _run(design, control, method, base, mode=True))


@pytest.mark.parametrize("design", ["panel", "unbalanced", "rc"])
@pytest.mark.parametrize("clustered", [False, True])
def test_bootstrap_bands_and_rng_match_in_both_modes(design, clustered):
    options = dict(bootstrap=True, clustered=clustered, base="universal")
    standard = _run(design, "nevertreated", mode=False, **options)
    fast = _run(design, "nevertreated", mode=True, **options)
    _assert_parity(standard, fast)
    _assert_parity(fast, _run(design, "nevertreated", **options), exact=True)


@pytest.mark.parametrize("bad", [
    {"x": np.array([1., np.nan])},
    {"x": np.array([1., np.inf])},
    {"x": np.array([1., -np.inf])},
    {"x": np.array([[1., 2.]])},
    {"x": np.array([1., 2.]), "extra": 3.},
    {},
])
def test_benchmark_comparator_rejects_missing_support_and_nonfinite_masks(bad):
    good = {"x": np.array([1., 2.])}
    assert compare(good, bad)["issues"]
    assert compare(bad, good)["issues"]
