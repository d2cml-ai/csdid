"""Adversarial regressions for R did master 74b88eb (PRs #273 and #275).

Balanced-sample identities test the actual public estimator and inference,
including missing covariates, missing rows and loss of the fallback cohort.
"""
import importlib
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from csdid.att_gt import ATTgt


def _panel(cohorts=(0, 3, 4), periods=range(1, 6), n_per=24):
    rng = np.random.default_rng(864)
    ids = np.repeat(np.arange(n_per * len(cohorts)), len(periods))
    t = np.tile(list(periods), n_per * len(cohorts))
    g = np.repeat(np.repeat(cohorts, n_per), len(periods))
    x = rng.normal(size=len(ids))
    y = 1 + x + 0.3 * t + ((g > 0) & (t >= g)) + rng.normal(size=len(ids))
    return pd.DataFrame({"id": ids, "period": t, "G": g, "X": x, "Y": y})


def _fit(data, method="dr", base="varying", **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y ~ X", **kwargs,
        ).fit(est_method=method, base_period=base, bstrap=False)


def _same_fit(fit, ref):
    for key in ("group", "t", "n"):
        np.testing.assert_array_equal(fit.MP[key], ref.MP[key])
    for key in ("att", "se"):
        np.testing.assert_allclose(fit.results[key], ref.results[key], atol=1e-11, rtol=1e-11)
    np.testing.assert_allclose(
        fit.MP["inffunc"]["inffunc"], ref.MP["inffunc"]["inffunc"],
        atol=1e-10, rtol=1e-10,
    )


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("control", ["nevertreated", "notyettreated"])
@pytest.mark.parametrize("method", ["dr", "reg", "ipw"])
@pytest.mark.parametrize("base", ["varying", "universal"])
@pytest.mark.parametrize("varying", [None, "varying"])
def test_never_treated_removed_by_balancing_matches_settled_sample(
        mode, control, method, base, varying):
    data = _panel()
    data.loc[(data.G == 0) & (data.period == 2), "X"] = np.nan
    opts = dict(faster_mode=mode, control_group=control, fix_weights=varying)
    fit = _fit(data, method, base, **opts)
    ref = _fit(data[data.G != 0], method, base, **opts)
    _same_fit(fit, ref)
    assert np.isfinite(fit.MP["att"]).all()
    assert 4 not in fit.MP["group"]


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("control", ["nevertreated", "notyettreated"])
@pytest.mark.parametrize("cohorts,kill", [
    ((0, 3, 4, 5), 4), ((0, 3, 4, 5), 5), ((3, 4, 5), 5),
])
def test_missing_entire_cohort_drops_only_that_cohort(mode, control, cohorts, kill):
    data = _panel(cohorts=cohorts)
    damaged = data[~((data.G == kill) & (data.period == 2))]
    opts = dict(faster_mode=mode, control_group=control)
    fit = _fit(damaged, **opts)
    _same_fit(fit, _fit(data[data.G != kill], **opts))
    assert kill not in fit.MP["group"]
    assert np.isfinite(fit.MP["att"]).all()


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("control", ["nevertreated", "notyettreated"])
def test_balancing_fallback_honors_anticipation_and_warns(mode, control):
    data = _panel()
    data.loc[(data.G == 0) & (data.period == 2), "X"] = np.nan
    with pytest.warns(UserWarning) as caught:
        obj = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=data, xformla="Y ~ X", faster_mode=mode,
            control_group=control, anticipation=1,
        )
    assert sum("after converting to balanced panel" in str(w.message) for w in caught) == 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = obj.fit(bstrap=False)
    np.testing.assert_array_equal(fit.MP["t"], [2])
    _same_fit(fit, _fit(data[data.G != 0], faster_mode=mode,
                       control_group=control, anticipation=1))


@pytest.mark.parametrize("mode", [False, True])
def test_unbalanced_allowed_retains_never_treated_and_latest_cohort(mode):
    data = _panel()
    data.loc[(data.G == 0) & (data.period == 2), "X"] = np.nan
    fit = _fit(data, faster_mode=mode, allow_unbalanced_panel=True)
    assert 4 in fit.MP["group"]
    assert fit.MP["n"] == 72
    np.testing.assert_array_equal(fit.dp["tlist"], [1, 2, 3, 4, 5])


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("failure", ["no_cohort", "no_period", "one_period"])
def test_empty_settled_sample_errors_name_the_cause(mode, failure):
    if failure == "no_cohort":
        data = _panel(cohorts=(0, 3))
        data.loc[(data.G == 0) & (data.period == 2), "X"] = np.nan
        opts, match = {}, "no other treated cohort remains"
    elif failure == "no_period":
        data, opts, match = _panel(cohorts=(1,)), {}, "No periods before"
    else:
        data = _panel(cohorts=(0, 3, 5), periods=(1, 4, 7))
        data.loc[(data.G == 0) & (data.period == 4), "X"] = np.nan
        opts, match = {"anticipation": 1}, "Only one time period remains"
    with pytest.raises(ValueError, match=match):
        _fit(data, faster_mode=mode, **opts)


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("same_period", [False, True])
def test_true_rc_rejects_repeated_ids_before_panel_specific_checks(mode, same_period):
    data = _panel()
    if same_period:
        data.loc[1, "period"] = data.loc[0, "period"]
    with pytest.raises(ValueError, match="idname must be unique when panel = FALSE"):
        _fit(data, panel=False, faster_mode=mode)


@pytest.mark.parametrize("mode", [False, True])
def test_true_rc_unique_ids_and_omitted_ids_are_identical(mode):
    data = _panel()
    data["id"] = np.arange(len(data))
    opts = dict(yname="Y", tname="period", gname="G", data=data,
                xformla="Y ~ X", panel=False, faster_mode=mode)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        supplied = ATTgt(idname="id", **opts).fit(bstrap=False)
        omitted = ATTgt(idname=None, **opts).fit(bstrap=False)
    _same_fit(supplied, omitted)


@pytest.mark.parametrize("mode", [False, True])
def test_missing_rc_ids_are_dropped_before_duplicate_id_check(mode):
    data = _panel()
    data["id"] = np.arange(len(data), dtype=float)
    data.loc[[0, 7], "id"] = np.nan
    _same_fit(_fit(data, panel=False, faster_mode=mode),
              _fit(data.dropna(subset=["id"]), panel=False, faster_mode=mode))


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("missing", ["Y", "X"])
def test_positive_weight_mean_is_checked_on_finite_sample(mode, missing):
    data = _panel()
    data["weight"] = 0.0
    data.loc[0, "weight"] = 1.0
    data.loc[0, missing] = np.nan
    with pytest.raises(ValueError, match="non-negative with a positive mean"):
        _fit(data, faster_mode=mode, weights_name="weight")


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("missing", ["Y", "X"])
def test_negative_weight_on_excluded_row_matches_explicit_deletion(mode, missing):
    data = _panel()
    data["weight"] = 1.0
    data.loc[0, "weight"] = -1.0
    data.loc[0, missing] = np.nan
    _same_fit(_fit(data, faster_mode=mode, weights_name="weight"),
              _fit(data.drop(index=0), faster_mode=mode, weights_name="weight"))


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("periods,cband", [(2, True), (4, False), (4, True)])
def test_bootstrap_se_and_band_critical_value_respect_resolved_cband(
        monkeypatch, mode, periods, cband):
    module = importlib.import_module("csdid.att_gt")
    bootstrap_se = 0.123

    def bootstrap(inffunc, dp, **kwargs):
        return {"crit_val": 4.5, "se": np.full(inffunc.shape[1], bootstrap_se),
                "V": np.eye(inffunc.shape[1])}

    monkeypatch.setattr(module, "mboot", bootstrap)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = ATTgt(
            yname="Y", tname="period", idname="id", gname="G",
            data=_panel(cohorts=(0, 2), periods=range(1, periods + 1)),
            xformla="Y ~ X", faster_mode=mode, cband=cband, alp=0.1,
        ).fit(est_method="reg", bstrap=True)
    expected_band = cband and periods > 2
    assert fit.dp["cband"] == expected_band
    np.testing.assert_allclose(fit.results["se"], bootstrap_se)
    np.testing.assert_allclose(fit.results["c"],
                               4.5 if expected_band else norm.ppf(0.95))


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("route", ["panel", "rc", "unbalanced", "varying"])
@pytest.mark.parametrize("compute_inffunc", [False, True])
def test_nan_estimate_warns_and_marks_entire_if_missing(
        monkeypatch, mode, route, compute_inffunc):
    # Deliberately return a NaN ATT with a finite IF. Only period 2 fails;
    # healthy later cells must survive with finite inference.
    def backend(*args, **kwargs):
        if "post" in kwargs:
            y, post = kwargs["y"], kwargs["post"]
            period = np.asarray(y)[np.asarray(post) == 1].mean()
        else:
            y = args[0]
            period = np.asarray(y).mean()
        return (np.nan if period == 2 else 1.0), np.linspace(-1, 1, len(y))

    module = importlib.import_module(
        "csdid.attgt_fnc.compute_att_gt2" if mode else
        "csdid.attgt_fnc.compute_att_gt")
    monkeypatch.setattr(module, "select_estimators", lambda method: (backend, backend))
    data = _panel(cohorts=(0, 2), periods=range(1, 5))
    data["Y"] = data.period.astype(float)
    if route == "rc":
        data["id"] = np.arange(len(data))
    elif route == "unbalanced":
        data = data.drop(data.index[-1])
    opts = dict(panel=route != "rc", allow_unbalanced_panel=route == "unbalanced",
                fix_weights="varying" if route == "varying" else None,
                compute_inffunc=compute_inffunc, faster_mode=mode)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = ATTgt(yname="Y", tname="period", idname="id", gname="G",
                    data=data, xformla="Y ~ X", **opts).fit(
                        est_method="reg", bstrap=False)
    assert sum("(g,t) = (2,2) is NaN" in str(w.message) for w in caught) == 1
    failed = np.asarray(fit.MP["t"]) == 2
    np.testing.assert_array_equal(np.isnan(fit.MP["att"]), failed)
    if compute_inffunc:
        np.testing.assert_array_equal(np.isnan(fit.results["se"]), failed)
        assert np.isnan(fit.MP["inffunc"]["inffunc"][:, failed]).all()
        assert np.isfinite(fit.MP["inffunc"]["inffunc"][:, ~failed]).all()
    else:
        assert fit.MP["inffunc"] is None
