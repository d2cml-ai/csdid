"""Anticipation and column-name regressions for R did PRs #276/#279.

Intercept-only cells have a separate weighted difference-in-means oracle,
including the complete influence vector on the original unit sample.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt


def _data():
    rng = np.random.default_rng(261026)
    n, periods = 90, np.arange(1, 9)
    ids = np.repeat(np.arange(n), len(periods))
    time = np.tile(periods, n)
    cohort = np.repeat(np.repeat([0, 5, 7], n // 3), len(periods))
    alpha = np.repeat(rng.normal(size=n), len(periods))
    anticipation = (cohort > 0) & (time >= cohort - 2) & (time < cohort)
    treated = (cohort > 0) & (time >= cohort)
    return pd.DataFrame({
        "id": ids, "period": time, "G": cohort,
        "Y": alpha + .2 * time + anticipation + 2 * treated + rng.normal(size=len(ids)),
        "weight": np.exp(.15 * np.sin(ids + time)),
    })


def _fit(data, mode, method="dr", base="varying", **options):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ATTgt(yname="Y", tname="period", idname="id", gname="G",
                     data=data, faster_mode=mode, **options).fit(
                         est_method=method, base_period=base, bstrap=False)


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("method", ["dr", "reg", "ipw"])
@pytest.mark.parametrize("base", ["varying", "universal"])
@pytest.mark.parametrize("weighted", [False, True])
def test_anticipation_cells_match_independent_weighted_did(mode, method, base, weighted):
    data = _data()
    fit = _fit(data, mode, method, base, anticipation=2,
               weights_name="weight" if weighted else None)
    n = data.id.nunique()
    for j, (g, t) in enumerate(zip(fit.MP["group"], fit.MP["t"])):
        if t < g - 2:
            continue
        before = data[data.period == g - 3].set_index("id")
        after = data[data.period == t].set_index("id")
        keep = before.G.isin([0, g])
        d = (before.G == g).to_numpy()
        weight = before.weight.to_numpy() if weighted else np.ones(n)
        delta = (after.Y - before.Y).to_numpy()
        mean_t = np.average(delta[keep & d], weights=weight[keep & d])
        mean_c = np.average(delta[keep & ~d], weights=weight[keep & ~d])
        influence = np.zeros(n)
        for mask, mean, sign in ((keep & d, mean_t, 1), (keep & ~d, mean_c, -1)):
            influence[mask] = sign * n * weight[mask] * (delta[mask] - mean) / weight[mask].sum()
        np.testing.assert_allclose(fit.MP["att"][j], mean_t - mean_c, atol=1e-11, rtol=1e-11)
        np.testing.assert_allclose(fit.MP["inffunc"]["inffunc"][:, j], influence,
                                   atol=1e-10, rtol=1e-10)
        np.testing.assert_allclose(fit.results["se"][j], np.sqrt(np.mean(influence ** 2) / n),
                                   atol=1e-11, rtol=1e-11)


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("control", ["nevertreated", "notyettreated"])
@pytest.mark.parametrize("time_name", ["calendar period", "time-period", "cutoff_t", "args"])
@pytest.mark.parametrize("never", [False, True])
def test_time_column_names_do_not_change_sample_or_results(mode, control, time_name, never):
    data = _data()
    if not never:
        data = data[data.G != 0]
    expected = _fit(data, mode, control_group=control)
    renamed = data.rename(columns={"period": time_name})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = ATTgt(yname="Y", tname=time_name, idname="id", gname="G",
                       data=renamed, faster_mode=mode, control_group=control).fit(bstrap=False)
    for key in ("group", "t", "att"):
        np.testing.assert_array_equal(actual.MP[key], expected.MP[key])
    np.testing.assert_array_equal(actual.results["se"], expected.results["se"])
    np.testing.assert_array_equal(actual.MP["inffunc"]["inffunc"],
                                  expected.MP["inffunc"]["inffunc"])


@pytest.mark.parametrize("mode", [False, True])
@pytest.mark.parametrize("route", ["panel", "unbalanced", "rc", "no_never", "balance_drop", "missing"])
@pytest.mark.parametrize("column", ["Y", "period", "id", "G", "X", "weight", "cluster"])
def test_argument_named_columns_preserve_preprocessing_across_routes(mode, route, column):
    data = _data()
    data["X"] = np.sin(data.id)
    data["cluster"] = data.id % 13
    options = dict(panel=route != "rc", allow_unbalanced_panel=route == "unbalanced")
    if route == "unbalanced":
        data = data.iloc[1:]
    elif route == "rc":
        data["id"] = np.arange(len(data))
    elif route == "no_never":
        data = data[data.G != 0]
    elif route == "balance_drop":
        data = data[~((data.G == 0) & (data.period == 1))]
    elif route == "missing":
        data.loc[0, "Y"] = np.inf
    names = dict(yname="Y", tname="period", idname="id", gname="G",
                 weights_name="weight", clustervar="cluster")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reference = ATTgt(data=data, xformla="~ X", faster_mode=mode,
                          **names, **options).fit(bstrap=False)
        renamed = {key: "args" if value == column else value for key, value in names.items()}
        actual = ATTgt(data=data.rename(columns={column: "args"}),
                       xformla="~ args" if column == "X" else "~ X", faster_mode=mode,
                       **renamed, **options).fit(bstrap=False)
    for key in ("group", "t", "att"):
        np.testing.assert_array_equal(actual.MP[key], reference.MP[key])
    np.testing.assert_array_equal(actual.results["se"], reference.results["se"])
    np.testing.assert_array_equal(actual.MP["inffunc"]["inffunc"],
                                  reference.MP["inffunc"]["inffunc"])

