"""Applicable numerical contracts from current R DRDID's test suite.

These tests exercise csdid's dispatched traditional DR and normalized IPW
estimators. Public ATTgt tests also cover regression, whose two-by-two kernel is
provided by the separate Python drdid dependency. IPT, TWFE, non-normalized IPW,
and R's exponential-weight refit bootstrap are not substituted for these methods.
The file-level audit is recorded in r_test_alignment_drdid.json.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from csdid.attgt_fnc import drdid_trim
from csdid.attgt_fnc.compute_att_gt_shared import select_estimators


METHODS = ("drdid_panel", "std_ipw_did_panel", "drdid_rc", "std_ipw_did_rc")


def _sample(zero_weights=False):
    rng = np.random.default_rng(8201)
    n = 240
    x, z = rng.normal(size=(2, n))
    cov = np.column_stack((np.ones(n), x, z))
    probability = 1 / (1 + np.exp(.3 - .5 * x + .2 * z - .4 * np.sin(2 * x)))
    treated = (rng.random(n) < probability).astype(float)
    post = (rng.random(n) < .55).astype(float)
    y0 = .8 + .3 * x - .2 * z + rng.normal(size=n)
    y1 = y0 + .6 + .4 * x + 1.1 * treated + rng.normal(size=n)
    weights = np.exp(.2 * x - .15 * z)
    if zero_weights:
        weights[::11] = 0
    return y1, y0, treated, post, cov, weights


def _arguments(name, sample):
    y1, y0, treated, post, cov, weights = sample
    outcomes = (post * y1 + (1 - post) * y0, post) if name.endswith("_rc") else (y1, y0)
    return (*outcomes, treated), cov, weights


@pytest.mark.parametrize("name", METHODS)
@pytest.mark.parametrize("zero_weights", [False, True])
def test_common_weight_scale_preserves_att_and_full_influence(name, zero_weights):
    """R sampling-weight tests compare the IF too, not only the ATT."""
    args, cov, weights = _arguments(name, _sample(zero_weights))
    estimator = getattr(drdid_trim, name)
    att, influence = estimator(*args, covariates=cov, i_weights=weights)
    assert np.isfinite(att) and np.isfinite(influence).all()
    for multiplier in (.001, 17.3):
        scaled_att, scaled_influence = estimator(
            *args, covariates=cov, i_weights=multiplier * weights)
        np.testing.assert_allclose(scaled_att, att, atol=1e-11, rtol=1e-11)
        np.testing.assert_allclose(scaled_influence, influence, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("name", METHODS)
def test_zero_weight_rows_preserve_observation_contributions_and_se(name):
    """Adding zero mass changes n and IF, but leaves IF/n and its SE unchanged."""
    args, cov, weights = _arguments(name, _sample(zero_weights=True))
    estimator = getattr(drdid_trim, name)
    att, influence = estimator(*args, covariates=cov, i_weights=weights)
    n = len(weights)
    extended_args = tuple(np.r_[values, [0., 1., 0., 1.]] for values in args)
    extended_cov = np.vstack((cov, [[1., -4., 2.], [1., 3., -1.],
                                   [1., .2, .1], [1., -1., 3.]]))
    extended_att, extended_if = estimator(
        *extended_args, covariates=extended_cov, i_weights=np.r_[weights, np.zeros(4)])
    np.testing.assert_allclose(extended_att, att, atol=1e-11, rtol=1e-11)
    np.testing.assert_array_equal(extended_if[n:], np.zeros(4))
    np.testing.assert_allclose(extended_if[:n] / (n + 4), influence / n,
                               atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(np.linalg.norm(extended_if) / (n + 4),
                               np.linalg.norm(influence) / n, atol=1e-11, rtol=1e-11)


@pytest.mark.parametrize("name", ["std_ipw_did_panel", "std_ipw_did_rc"])
@pytest.mark.parametrize("zero_weights", [False, True])
def test_normalized_ipw_influence_differentiates_refitted_observation_weights(name, zero_weights):
    """Independent derivative oracle for the normalized estimator actually used.

    R's ipw_did_rc normalization tests concern a different, non-normalized
    functional. Only their empirical-derivative test strategy is reused here.
    Every perturbed call renormalizes weights and refits the propensity model.
    """
    sample = _sample(zero_weights)
    args, cov, weights = _arguments(name, sample)
    estimator = getattr(drdid_trim, name)
    _, influence = estimator(*args, covariates=cov, i_weights=weights)
    treated, post = sample[2:4]
    indices = []
    for d in (0, 1):
        for t in (0, 1):
            cell = np.flatnonzero((treated == d) & (post == t))
            indices.extend(cell[np.argsort(np.abs(influence[cell]))[-2:]])
    indices.extend(np.flatnonzero(weights == 0)[:2])
    h, n = .001, len(weights)
    derivatives = []
    for i in indices:
        plus, minus = weights.copy(), weights.copy()
        plus[i] *= 1 + h
        minus[i] *= 1 - h
        att_plus = estimator(*args, covariates=cov, i_weights=plus)[0]
        att_minus = estimator(*args, covariates=cov, i_weights=minus)[0]
        derivatives.append((att_plus - att_minus) / (h / (n + h) + h / (n - h)))
    np.testing.assert_allclose(influence[indices], derivatives, atol=2e-5, rtol=2e-6)


@pytest.mark.parametrize("name", METHODS)
def test_column_treatment_indicator_preserves_vector_result(name):
    args, cov, weights = _arguments(name, _sample())
    estimator = getattr(drdid_trim, name)
    expected = estimator(*args, covariates=cov, i_weights=weights)
    actual = estimator(*args[:2], args[2][:, None], covariates=cov, i_weights=weights)
    np.testing.assert_equal(actual[0], expected[0])
    np.testing.assert_array_equal(actual[1], expected[1])


def _public_frame(panel):
    y1, y0, treated, _, cov, weights = _sample(zero_weights=True)
    n = len(y1)
    data = pd.DataFrame(dict(id=np.tile(np.arange(n), 2),
                             period=np.repeat([1, 2], n), G=np.tile(2 * treated, 2),
                             Y=np.r_[y0, y1], X=np.tile(cov[:, 1], 2),
                             wt=np.tile(weights, 2)))
    if not panel:
        data["id"] = np.arange(len(data))
    return data


def _public_fit(data, panel, method, *, bootstrap=False, covariates=False, mode=True):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                     panel=panel, weights_name="wt", xformla="~X" if covariates else "~1",
                     faster_mode=mode, cband=False, biters=1999).fit(method, bstrap=bootstrap)


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("mode", [False, True])
def test_no_covariate_public_methods_share_closed_form_att_and_influence(panel, mode):
    """Covers the supported OR dispatch too, without assuming panel and RC SEs match."""
    data = _public_frame(panel)
    pre, post = data.loc[data.period == 1], data.loc[data.period == 2]

    def contrast(frame):
        return (np.average(frame.loc[frame.G == 2, "Y"], weights=frame.loc[frame.G == 2, "wt"])
                - np.average(frame.loc[frame.G == 0, "Y"], weights=frame.loc[frame.G == 0, "wt"]))

    expected = contrast(post) - contrast(pre)
    results = [_public_fit(data, panel, method, mode=mode) for method in ("dr", "ipw", "reg")]
    for result in results:
        np.testing.assert_allclose(result.MP["att"], [expected], atol=1e-11, rtol=1e-11)
        np.testing.assert_allclose(result.results["se"], results[0].results["se"],
                                   atol=1e-10, rtol=1e-10)
        np.testing.assert_allclose(result.MP["inffunc"]["inffunc"],
                                   results[0].MP["inffunc"]["inffunc"], atol=1e-9, rtol=1e-9)


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
def test_public_multiplier_bootstrap_preserves_point_estimate_and_influence(panel, method):
    data = _public_frame(panel)
    analytic = _public_fit(data, panel, method, covariates=True)
    before = np.random.get_state()
    try:
        np.random.seed(9012)
        bootstrap = _public_fit(data, panel, method, covariates=True, bootstrap=True)
    finally:
        np.random.set_state(before)
    np.testing.assert_array_equal(bootstrap.MP["att"], analytic.MP["att"])
    np.testing.assert_array_equal(bootstrap.MP["inffunc"]["inffunc"],
                                  analytic.MP["inffunc"]["inffunc"])
    assert np.isfinite(bootstrap.results["se"]).all()
    np.testing.assert_allclose(bootstrap.results["se"], analytic.results["se"], rtol=.2)


@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("mode", [False, True])
def test_weighted_public_two_by_two_matches_explicit_estimator_inputs(panel, method, mode):
    """Check long-form wrapper alignment against explicit paired/RC arrays."""
    data = _public_frame(panel)
    result = _public_fit(data.sample(frac=1, random_state=94), panel, method,
                         covariates=True, mode=mode)
    estimator = select_estimators(method)[0 if panel else 1]
    if panel:
        pre = data.loc[data.period == 1].set_index("id").sort_index()
        post = data.loc[data.period == 2].set_index("id").reindex(pre.index)
        args = (post.Y.to_numpy(), pre.Y.to_numpy(), (pre.G == 2).to_numpy(dtype=float))
        rows, ids = pre, pre.index
    else:
        rows = data.sort_values("id")
        args = (rows.Y.to_numpy(), (rows.period == 2).to_numpy(dtype=float),
                (rows.G == 2).to_numpy(dtype=float))
        ids = rows.id
    cov = np.column_stack((np.ones(len(rows)), rows.X.to_numpy()))
    expected_att, expected_if = estimator(*args, covariates=cov, i_weights=rows.wt.to_numpy())
    np.testing.assert_allclose(result.MP["att"], [expected_att], atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(result.results["se"], [np.linalg.norm(expected_if) / len(rows)],
                               atol=1e-10, rtol=1e-10)
    actual_rows = result.dp["data"]
    if panel:
        actual_rows = actual_rows.loc[actual_rows.period == 1]
    actual_if = pd.Series(np.asarray(result.MP["inffunc"]["inffunc"]).ravel(),
                          index=actual_rows.id).reindex(ids).to_numpy()
    np.testing.assert_allclose(actual_if, expected_if, atol=1e-9, rtol=1e-9)


@pytest.mark.parametrize("cell", range(4))
@pytest.mark.parametrize("mode", [False, True])
def test_each_undersized_rc_outcome_cell_is_guarded_before_estimation(cell, mode):
    """Port R's four-cell support fixture to csdid's warning-and-NA contract."""
    from itertools import product

    # The R test calls its low-level estimator directly. Enlarge unaffected
    # cells so csdid's earlier public cohort-size check allows this design to
    # reach the per-outcome-cell guard, keeping the affected cell at 3 rows.
    support = np.tile(np.asarray(list(product((-1., 1.), repeat=3))), (3, 1))
    cov = np.tile(support, (4, 1))
    treated = np.repeat([0., 0., 1., 1.], 24)
    post = np.repeat([0., 1., 0., 1.], 24)
    data = pd.DataFrame(dict(id=np.arange(96), period=post + 1, G=2 * treated,
                             Y=.3 * cov[:, 0] + .2 * cov[:, 1] + treated * post
                             + np.arange(96) / 100, X=cov[:, 0], Z=cov[:, 1], Q=cov[:, 2]))
    data = data.drop(index=np.arange(cell * 24 + 3, (cell + 1) * 24))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                       panel=False, xformla="~X+Z+Q", faster_mode=mode).fit("dr", bstrap=False)
    assert any("singular" in str(w.message) for w in caught)
    assert np.isnan(result.MP["att"]).all()
    assert np.isnan(result.MP["inffunc"]["inffunc"]).all()


@pytest.mark.parametrize("missing", ["pre_Y", "post_Y", "pre_X", "post_X", "pre_wt", "post_wt", "mixed"])
@pytest.mark.parametrize("panel", [False, True])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("mode", [False, True])
def test_missing_inputs_match_independently_retained_rows_and_units(missing, panel, method, mode):
    """Missingness must preserve outcomes, designs, weights, and every IF row.

    Panel comparisons explicitly request complete pairs. RC comparisons retain
    complete rows independently. The reference sample is constructed without
    calling package preprocessing, and IF comparisons use original user IDs.
    """
    data = _public_frame(True)
    # Retain the original, deliberately sparse identifiers through both routes.
    data["id"] = 100 + 13 * data.id
    if not panel:
        data["id"] = 500 + 17 * np.arange(len(data))
    data["unused"] = np.nan
    if missing == "mixed":
        data.loc[3, "X"] = np.nan
        data.loc[245, "wt"] = np.nan
        data.loc[68, "Y"] = np.nan
    else:
        period, column = missing.split("_")
        data.loc[3 if period == "pre" else 243, column] = np.nan
    # Independent imbalance and shuffled input prevent positional pairing from
    # accidentally satisfying the numerical comparison.
    data = data.drop(index=9).sample(frac=1, random_state=623)
    complete = np.isfinite(data[["Y", "X", "wt"]]).all(axis=1)
    if panel:
        observed_twice = data.groupby("id").size().eq(2)
        retained = set(observed_twice[observed_twice].index) - set(data.loc[~complete, "id"])
        reference = data.loc[data.id.isin(retained)].copy()
    else:
        reference = data.loc[complete].copy()
        retained = set(reference.id)

    def fit(frame):
        return ATTgt(yname="Y", tname="period", idname="id", gname="G", data=frame,
                     panel=panel, allow_unbalanced_panel=False, weights_name="wt", xformla="~X",
                     faster_mode=mode, cband=False).fit(method, bstrap=False)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = fit(data)
    with warnings.catch_warnings(record=True) as expected_warnings:
        warnings.simplefilter("always")
        expected = fit(reference)
    assert not any("missing" in str(w.message).lower() or "dropped" in str(w.message).lower()
                   for w in expected_warnings)
    assert any("missing" in str(w.message).lower() or "dropped" in str(w.message).lower()
               for w in caught)
    assert actual.dp["panel"] == expected.dp["panel"] == panel
    assert actual.dp["n"] == expected.dp["n"] == len(retained)
    assert set(actual.dp["data"].id) == set(expected.dp["data"].id) == retained
    for field in ("att", "se", "c", "l_se", "u_se"):
        np.testing.assert_allclose(actual.results[field], expected.results[field],
                                   atol=1e-10, rtol=1e-10)

    def aligned_influence(result):
        frame = result.dp["data"]
        if panel:
            frame = frame.loc[frame.period == result.dp["tlist"][0]]
        # Data are in the unit/rowid order used by the influence-function matrix.
        influence = np.asarray(result.MP["inffunc"]["inffunc"])
        assert len(frame) == influence.shape[0]
        assert frame.id.is_unique
        return pd.DataFrame(influence, index=frame.id).sort_index()

    pd.testing.assert_frame_equal(aligned_influence(actual), aligned_influence(expected),
                                   check_exact=False, atol=1e-10, rtol=1e-10)


def _manual_weighted_ratio(effects, influence, cohorts, unit_cohorts, weights):
    """Differentiate the cohort-weighted ratio directly, without aggte helpers."""
    mass_by_unit = weights[:, None] * (unit_cohorts[:, None] == cohorts)
    masses = mass_by_unit.mean(axis=0)
    denominator = masses.sum()
    estimate = effects @ masses / denominator
    derivative = influence @ (masses / denominator)
    derivative += mass_by_unit @ (effects - estimate) / denominator
    return estimate, derivative


@pytest.mark.parametrize("design", ["balanced", "unbalanced", "rc"])
@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("mode", [False, True])
def test_varying_weight_aggregations_use_sampling_design_cohort_masses(design, method, mode):
    """An RC estimator inside a balanced panel retains first-period group masses.

    Genuine unbalanced panels use each unit's mean weight; true repeated cross
    sections use each observation's weight. All IF rows are aligned by user ID.
    """
    rng = np.random.default_rng(7620)
    count = 180
    unit_ids = 100 + 17 * np.arange(count)
    cohorts = np.repeat([0, 3, 4], count // 3)
    x = rng.normal(size=count)
    frames = []
    for period in range(1, 5):
        effect = (cohorts > 0) * (period >= cohorts) * (1 + .4 * cohorts + .2 * period)
        frames.append(pd.DataFrame(dict(
            id=unit_ids, period=period, G=cohorts, X=x + .05 * period,
            Y=.4 * x + .2 * period + effect + rng.normal(size=count),
            wt=np.exp(.1 * x + .15 * period * (cohorts == 3)
                      - .12 * period * (cohorts == 4)))))
    data = pd.concat(frames, ignore_index=True)
    if design == "unbalanced":
        data = data.drop(index=[2, count + 67, 2 * count + 130])
    elif design == "rc":
        data["id"] = 500 + 19 * np.arange(len(data))
    data = data.sample(frac=1, random_state=201)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                       panel=design != "rc", allow_unbalanced_panel=True,
                       weights_name="wt", fix_weights="varying", xformla="~X",
                       faster_mode=mode, cband=False).fit(method, bstrap=False)

    retained = result.dp["data"]
    internal_id = result.dp["idname"]
    identity = retained[[internal_id, "id"]].iloc[:, [0, -1]].copy()
    identity.columns = ["internal_id", "user_id"]
    identity = identity.drop_duplicates().set_index("internal_id")["user_id"]
    if_ids = identity.reindex(retained[internal_id].unique()).to_numpy()
    assert len(if_ids) == len(set(if_ids)) == result.MP["n"]
    if design == "balanced":
        units = data.loc[data.period == 1].set_index("id")[["G", "wt"]]
        # Ensure this fixture distinguishes first-period and mean-period masses.
        first_mass = units.groupby("G").wt.sum()
        mean_mass = data.groupby("id")[["G", "wt"]].mean().groupby("G").wt.sum()
        assert not np.allclose(first_mass / first_mass.sum(), mean_mass / mean_mass.sum())
    elif design == "unbalanced":
        units = data.groupby("id")[["G", "wt"]].mean()
    else:
        units = data.set_index("id")[["G", "wt"]]
    units = units.reindex(if_ids)
    assert units.notna().all().all()
    unit_cohorts, weights = units.G.to_numpy(), units.wt.to_numpy()
    group, times, effects = (np.asarray(result.MP[key]) for key in ("group", "t", "att"))
    influence = np.asarray(result.MP["inffunc"]["inffunc"])

    def ratio(selected):
        return _manual_weighted_ratio(effects[selected], influence[:, selected],
                                      group[selected], unit_cohorts, weights)

    for kind in ("simple", "group", "dynamic", "calendar"):
        result.aggte(typec=kind, bstrap=False, cband=False)
        actual = result.atte
        if kind == "simple":
            expected_att, expected_if = ratio(times >= group)
            overall_key = "simple_att"
        else:
            if kind == "group":
                support = np.unique(group)
                selected = [(group == g) & (times >= g) for g in support]
                event_att = np.array([effects[keep].mean() for keep in selected])
                event_if = np.column_stack([influence[:, keep].mean(axis=1) for keep in selected])
                expected_att, expected_if = _manual_weighted_ratio(
                    event_att, event_if, support, unit_cohorts, weights)
                event_key, overall_key = "selective_inf_func_g", "selective_inf_func"
            else:
                if kind == "dynamic":
                    support = np.unique(times - group)
                    selected = [times - group == e for e in support]
                    positive = support >= 0
                    event_key, overall_key = "dynamic_inf_func_e", "dynamic_inf_func"
                else:
                    support = np.unique(times[times >= group])
                    selected = [(times == t) & (times >= group) for t in support]
                    positive = np.ones(len(support), dtype=bool)
                    event_key, overall_key = "calendar_inf_func_t", "calendar_inf_func"
                pairs = [ratio(keep) for keep in selected]
                event_att = np.array([pair[0] for pair in pairs])
                event_if = np.column_stack([pair[1] for pair in pairs])
                expected_att = event_att[positive].mean()
                expected_if = event_if[:, positive].mean(axis=1)
            np.testing.assert_array_equal(actual["egt"], support)
            np.testing.assert_allclose(actual["att_egt"], event_att, atol=1e-12, rtol=1e-12)
            np.testing.assert_allclose(actual["inf_function"][event_key], event_if,
                                       atol=1e-12, rtol=1e-12)
            np.testing.assert_allclose(np.asarray(actual["se_egt"]).reshape(-1),
                                       np.sqrt(np.sum(event_if ** 2, axis=0)) / len(if_ids),
                                       atol=1e-12, rtol=1e-12)
        np.testing.assert_allclose(actual["overall_att"], expected_att, atol=1e-12, rtol=1e-12)
        np.testing.assert_allclose(actual["inf_function"][overall_key].reshape(-1),
                                   expected_if, atol=1e-12, rtol=1e-12)
        np.testing.assert_allclose(actual["overall_se"], np.linalg.norm(expected_if) / len(if_ids),
                                   atol=1e-12, rtol=1e-12)
