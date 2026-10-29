"""Independent before/after statistical audit; no R or project test helpers.

Run from any directory with the dependency-complete Python interpreter::

    python inference_adversarial.py --baseline /clean/old-checkout \
        --candidate /clean/new-checkout --output inference_adversarial.json.gz \
        --summary inference_adversarial_summary.json

Each source runs in a separate process, with imports pinned and verified before
any case runs. The sample-mean, ratio-derivative, covariance, and order-statistic
oracles below do not import csdid implementations. Failed imports abort; they
cannot become evidence of an estimator defect. Inputs, seeds, source identities,
per-case numerical errors and all failed-case definitions remain in the JSON.

Scope: intercept-only panel reductions and synthetic valid group-time objects;
bootstrap calculations conditional on explicit multipliers. This establishes
finite-sample implementation identities, not causal identification or nominal
repeated-sample rejection rates. Limits of floating-point exponent range remain.
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import hashlib
import gzip
import importlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

POLICY = {"outcome_normalized_atol": 1e-10, "outcome_normalized_rtol": 1e-9,
          "wald_atol": 1e-8, "wald_rtol": 1e-9,
          "support_and_nonfinite_masks": "exact",
          "dimensionful_comparison": "Divide by the known absolute outcome multiplier before comparing. This independent equation gate does not replace the permanent raw-scale R parity gate or clear its retained large-scale IF cancellation failures."}
ORACLE_CORRECTIONS = [
    {
        "stage": "First exploratory 1,000-balance-case run",
        "observed": "23 apparent current-source failures; no product conclusion was drawn.",
        "diagnosis": "The first oracle retained missing negative-event cells after na_rm, although those cells do not enter the post-treatment balance requirement and must be omitted from the displayed event-specific averages. It also assigned statistical meaning to the overall-IF placeholder when the display window contained no post-treatment event.",
        "correction": "Build displayed events and their cell averages from finite estimates, while retaining original support solely for the balance requirement. With no displayed post-treatment event, require unavailable overall ATT/SE and mark its internal IF placeholder not applicable.",
        "unchanged": "Data seeds, all 1,000 case definitions, source checkouts, numerical tolerances, and production code.",
        "disposition": "Correct oracle; rerun both isolated sources. All 1,000 current-source cases then pass. This is not counted as a product fix."
    }
]
THREADS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def git(source, *args):
    return subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metrics(actual, expected, units=1., wald=False):
    a = np.asarray(np.nan if actual is None else actual, dtype=float)
    b = np.asarray(np.nan if expected is None else expected, dtype=float)
    same_shape = a.shape == b.shape
    if not same_shape:
        return dict(passed=False, actual_shape=list(a.shape), expected_shape=list(b.shape))
    masks = bool(np.array_equal(np.isnan(a), np.isnan(b)) and
                 np.array_equal(np.isposinf(a), np.isposinf(b)) and
                 np.array_equal(np.isneginf(a), np.isneginf(b)))
    finite = np.isfinite(a) & np.isfinite(b)
    error = np.abs(a[finite] / units - b[finite] / units)
    bound = ((POLICY["wald_atol"] if wald else POLICY["outcome_normalized_atol"]) +
             (POLICY["wald_rtol"] if wald else POLICY["outcome_normalized_rtol"]) *
             np.abs(b[finite] / units))
    return dict(passed=masks and bool(np.all(error <= bound)), same_nonfinite_masks=masks,
                max_abs_error=float(np.max(np.abs(a[finite] - b[finite]), initial=0)),
                max_unit_normalized_abs_error=float(np.max(error, initial=0)),
                finite_comparisons=int(finite.sum()))


def record(name, parameters, fields, **extra):
    return dict(name=name, parameters=parameters, passed=all(v["passed"] for v in fields.values()),
                fields=fields, **extra)


def make_panel(clock):
    rng = np.random.default_rng(180226)
    if clock == "regular":
        periods, cohorts = np.arange(1, 9), np.array([0, 5, 7])
    elif clock == "irregular_negative":
        periods, cohorts = np.array([-4, -1, 1, 2, 5, 6, 8, 11]), np.array([0, 6, 11])
    else:
        periods = np.array([.5, 1.5, 3.5, 4.5, 5.5, 7.5, 8.5, 9.5])
        cohorts = np.array([0., 5.5, 8.5])
    n = 72
    ids = np.repeat(np.arange(n), len(periods))
    t = np.tile(periods, n)
    g = np.repeat(np.repeat(cohorts, [21, 24, 27]), len(periods))
    a = np.repeat(rng.normal(size=n), len(periods))
    noise = rng.normal(scale=.3, size=len(t))
    y = a + .2 * t + noise + ((g > 0) & (g - 2 <= t) & (t < g)) + 2 * ((g > 0) & (t >= g))
    return pd.DataFrame(dict(id=ids, period=t, G=g, Y=y,
                             weight=np.exp(.3 * np.sin(ids + t)), cluster=ids % 17))


def panel_oracle(data, groups, times, anticipation, base, control, clustered):
    periods = np.sort(data.period.unique())
    ids = np.sort(data.id.unique())
    n = len(ids)
    cohort = data[data.period == periods[0]].sort_values("id").G.to_numpy()
    cluster = ids % 17
    estimates, influence = [], []
    for g, t in zip(groups, times):
        unaffected = max(tt for tt in periods if tt < g - anticipation)
        before = unaffected if base == "universal" or t >= g - anticipation else max(tt for tt in periods if tt < t)
        if t == before:
            estimates.append(0.)
            influence.append(np.zeros(n))
            continue
        p = data[data.period == before].sort_values("id")
        q = data[data.period == t].sort_values("id")
        earlier = data[data.period == min(t, before)].sort_values("id")
        dy, w = q.Y.to_numpy() - p.Y.to_numpy(), earlier.weight.to_numpy()
        treated = cohort == g
        controls = cohort == 0
        if control == "notyettreated":
            controls |= (cohort > max(t, before) + anticipation) & (cohort != g)
        means, inf = [], np.zeros(n)
        for mask, sign in ((treated, 1), (controls, -1)):
            average = sum(w[mask] * dy[mask]) / sum(w[mask])
            means.append(average)
            inf[mask] = sign * n * w[mask] * (dy[mask] - average) / sum(w[mask])
        estimates.append(means[0] - means[1])
        influence.append(inf)
    estimates, influence = np.asarray(estimates), np.column_stack(influence)
    if clustered:
        summed = np.vstack([influence[cluster == c].sum(axis=0) for c in sorted(set(cluster))])
    else:
        summed = influence
    covariance = summed.T @ summed / n
    se = np.sqrt(np.diag(covariance) / n)
    se[(se <= 0) | ~np.isfinite(se)] = np.nan
    pre = (np.asarray(times) < np.asarray(groups) - anticipation) & np.isfinite(se)
    if pre.any():
        v = covariance[np.ix_(pre, pre)]
        if 1 / np.linalg.cond(v, p=1) > np.finfo(float).eps:
            wald = n * estimates[pre] @ np.linalg.solve(v, estimates[pre])
            pvalue = round(float(chi2.sf(wald, sum(pre))), 5)
        else:
            wald = pvalue = None
    else:
        wald = pvalue = None
    return estimates, influence, se, wald, pvalue


def audit_panels(ATTgt):
    rows = []
    for clock in ("regular", "irregular_negative", "fractional"):
        template = make_panel(clock)
        for faster in (False, True):
            for method in ("dr", "reg", "ipw"):
                for base in ("varying", "universal"):
                    for control in ("nevertreated", "notyettreated"):
                        for anticipation in (0, 2):
                            for scale, clustered in ((1., False), (-1e-12, True), (1e6, True)):
                                params = dict(clock=clock, faster_mode=faster, method=method,
                                              base=base, control=control, anticipation=anticipation,
                                              scale=scale, clustered=clustered)
                                data = template.assign(Y=template.Y * scale)
                                try:
                                    fit = ATTgt("Y", "period", "id", "G", data=data, xformla="~1",
                                                weights_name="weight", anticipation=anticipation,
                                                control_group=control, faster_mode=faster,
                                                clustervar="cluster" if clustered else None).fit(
                                                    est_method=method, base_period=base, bstrap=False)
                                    a, inf, se, wald, pvalue = panel_oracle(
                                        data, fit.MP["group"], fit.MP["t"], anticipation, base, control, clustered)
                                    expected_groups = np.sort(data.G.unique()[data.G.unique() > 0])
                                    expected_times = np.sort(data.period.unique())
                                    if base == "varying":
                                        expected_times = expected_times[1:]
                                    fields = dict(
                                                  group_support=dict(passed=bool(np.array_equal(fit.MP["group"], np.repeat(expected_groups, len(expected_times))))),
                                                  time_support=dict(passed=bool(np.array_equal(fit.MP["t"], np.tile(expected_times, len(expected_groups))))),
                                                  sample_size=dict(passed=int(fit.MP["n"]) == data.id.nunique()),
                                                  att=metrics(fit.MP["att"], a, abs(scale)),
                                                  inffunc=metrics(fit.MP["inffunc"]["inffunc"], inf, abs(scale)),
                                                  se=metrics(fit.results["se"], se, abs(scale)),
                                                  l_se=metrics(fit.results["l_se"], a - norm.ppf(.975) * se, abs(scale)),
                                                  u_se=metrics(fit.results["u_se"], a + norm.ppf(.975) * se, abs(scale)),
                                                  critical=metrics(fit.results["c"], norm.ppf(.975)),
                                                  W=metrics(fit.MP["W"], wald, wald=True),
                                                  Wpval=metrics(fit.MP["Wpval"], pvalue, wald=True))
                                    rows.append(record("public_panel_equations", params, fields,
                                                       number_of_cells=len(a),
                                                       expected_W=None if wald is None else float(wald),
                                                       observed_W=fit.MP["W"]))
                                except Exception as exc:
                                    rows.append(record("public_panel_equations", params, {},
                                                       passed_exception=False, exception=type(exc).__name__, error=str(exc)))
                                    rows[-1]["passed"] = False
    return rows


def make_mp(rng, case):
    # Include shifted/fractional/large clocks and cohorts off the observed grid.
    time_units = (.5, 1., 3.)[case % 3]
    origin = (-10., 0., 1e8)[(case // 3) % 3]
    grid = origin + time_units * np.r_[0, np.cumsum(rng.integers(1, 4, size=5))]
    candidates = grid[1:4] + (0. if case % 4 else .25 * time_units)
    if candidates.min() <= 0:
        grid += 20
        candidates += 20
    cohorts = candidates[:int(rng.integers(2, 4))]
    groups = np.repeat(cohorts, len(grid) - 1)
    times = np.tile(grid[1:], len(cohorts))
    sizes = rng.integers(4, 16, len(cohorts) + 1)
    unit_group = np.repeat(np.r_[0., cohorts], sizes)
    n = len(unit_group)
    unit_weights = np.exp(rng.normal(scale=.4, size=n))
    data = pd.DataFrame(dict(id=np.repeat(np.arange(n), len(grid)),
                             period=np.tile(grid, n), G=np.repeat(unit_group, len(grid)),
                             w1=np.repeat(unit_weights, len(grid))))
    influence = rng.normal(size=(n, len(groups)))
    influence -= influence.mean(axis=0)
    att = 2 + .001 * (times - groups) + rng.normal(size=len(groups))
    missing = rng.random(len(att)) < .2
    att[missing] = np.nan
    influence[:, missing] = np.nan
    scale = (1., -1e-12, 1e6)[case % 3]
    horizon = int(rng.integers(0, 13))
    maximum = (0., math.inf, float(rng.integers(1, 9)))[case % 3]
    minimum = (0., -math.inf, -2.)[case % 3]
    return dict(group=groups, t=times, att=att * scale, n=n,
                inffunc={"inffunc": influence * scale}, DIDparams=dict(
                    data=data, idname="id", tname="period", gname="G", panel=True,
                    tlist=grid, glist=cohorts, bstrap=False, cband=False, biters=199,
                    clustervars=None, alp=.05, control_group="nevertreated",
                    est_method="reg", anticipation=0)), horizon, minimum, maximum, scale


def ratio_oracle(mp, horizon, minimum, maximum):
    g, t, a = (mp[k] for k in ("group", "t", "att"))
    full = list(zip(g, t, a))
    max_t, first_t = max(t), min(mp["DIDparams"]["tlist"])
    eligible = {group for group in g if max_t - group >= horizon}
    required = {time - group for group, time, _ in full
                if group in eligible and 0 <= time - group <= horizon}
    complete = {group for group in eligible if required and all(
        any(gg == group and tt - gg == event and np.isfinite(value)
            for gg, tt, value in full) for event in required)}
    if not complete:
        return "no_complete_cohort", None
    events = sorted({tt - gg for gg, tt, value in full if gg in complete and np.isfinite(value) and
                     horizon - max_t + first_t <= tt - gg <= horizon and minimum <= tt - gg <= maximum})
    if not events:
        return "empty_window", None
    units = mp["DIDparams"]["data"].drop_duplicates("id")
    estimates, influences = [], []
    for e in events:
        columns = [i for i, (gg, tt, _) in enumerate(full) if gg in complete and tt - gg == e and np.isfinite(a[i])]
        values = a[columns]
        mass = units.w1.to_numpy()[:, None] * (units.G.to_numpy()[:, None] == g[columns])
        probabilities = mass.mean(axis=0)
        mean = float(probabilities @ values / probabilities.sum())
        influence = mp["inffunc"]["inffunc"][:, columns] @ (probabilities / probabilities.sum())
        influence += mass @ (values - mean) / probabilities.sum()
        estimates.append(mean)
        influences.append(influence)
    estimates, influences = np.asarray(estimates), np.column_stack(influences)
    se = np.sqrt((influences ** 2).sum(axis=0)) / mp["n"]
    se[(se <= 0) | ~np.isfinite(se)] = np.nan
    post = np.asarray(events) >= 0
    overall_att = estimates[post].mean() if post.any() else np.nan
    overall_if = influences[:, post].mean(axis=1) if post.any() else np.full(mp["n"], np.nan)
    overall_se = np.sqrt(sum(overall_if ** 2)) / mp["n"]
    if not np.isfinite(overall_se) or overall_se <= 0:
        overall_se = np.nan
    return None, dict(egt=events, att=estimates, inffunc=influences, se=se,
                      overall_att=overall_att, overall_if=overall_if, overall_se=overall_se, has_post=bool(post.any()),
                      cohorts=sorted(complete), required_event_times=sorted(required))


def audit_balance(compute_aggte, count):
    rng = np.random.default_rng(202610261)
    rows = []
    for case in range(count):
        mp, horizon, minimum, maximum, scale = make_mp(rng, case)
        expected_error, expected = ratio_oracle(mp, horizon, minimum, maximum)
        params = dict(case=case, horizon=horizon, min_e=minimum if np.isfinite(minimum) else "-Inf",
                      max_e=maximum if np.isfinite(maximum) else "Inf", scale=scale,
                      periods=mp["DIDparams"]["tlist"].tolist(), cohorts=mp["DIDparams"]["glist"].tolist(),
                      missing_cells=np.flatnonzero(~np.isfinite(mp["att"])).tolist())
        try:
            out = compute_aggte(copy.deepcopy(mp), typec="dynamic", balance_e=horizon,
                                min_e=minimum, max_e=maximum, na_rm=True)
            if expected_error:
                rows.append(record("dynamic_ratio_equations", params,
                                   dict(error_contract=dict(passed=False)), expected_error=expected_error,
                                   observed_events=list(out["egt"])))
                continue
            fields = dict(egt=metrics(out["egt"], expected["egt"]),
                          att=metrics(out["att_egt"], expected["att"], abs(scale)),
                          inffunc=metrics(out["inf_function"]["dynamic_inf_func_e"], expected["inffunc"], abs(scale)),
                          se=metrics(np.asarray(out["se_egt"]).ravel(), expected["se"], abs(scale)),
                          overall_att=metrics(out["overall_att"], expected["overall_att"], abs(scale)),
                          overall_if=(metrics(np.asarray(out["inf_function"]["dynamic_inf_func"]).ravel(),
                                              expected["overall_if"], abs(scale)) if expected["has_post"] else
                                      dict(passed=True, applicable=False, reason="No displayed post-treatment event; overall target and SE must be unavailable. Its internal IF placeholder has no statistical interpretation.")),
                          overall_se=metrics(out["overall_se"], expected["overall_se"], abs(scale)))
            rows.append(record("dynamic_ratio_equations", params, fields,
                               expected_cohorts=expected["cohorts"], required_event_times=expected["required_event_times"]))
        except Exception as exc:
            correct_error = (isinstance(exc, ValueError) and
                             ((expected_error == "no_complete_cohort" and "No groups have non-missing ATT" in str(exc)) or
                              (expected_error == "empty_window" and "No event times" in str(exc))))
            rows.append(record("dynamic_ratio_equations", params,
                               dict(error_contract=dict(passed=correct_error)), expected_error=expected_error,
                               exception=type(exc).__name__, error=str(exc)))
    return rows


def order_quantile(array, probability):
    ordered = np.sort(array, axis=0)
    return ordered[math.ceil(probability * len(ordered)) - 1]


def bootstrap_oracle(influence, multipliers, clusters):
    n = len(influence)
    if clusters is None:
        grouped = influence
    else:
        grouped = np.vstack([influence[clusters == c].sum(axis=0) for c in sorted(set(clusters))])
    nc = len(grouped)
    bres = np.sqrt(nc) * (multipliers @ grouped / nc)
    live = np.array([all(np.isfinite(bres[:, j])) and any(bres[:, j] != 0)
                     for j in range(bres.shape[1])])
    selected = bres[:, live]
    se = np.full(influence.shape[1], np.nan)
    if not live.any():
        return selected, se, np.nan
    scales = (order_quantile(selected, .75) - order_quantile(selected, .25)) / (norm.ppf(.75) - norm.ppf(.25))
    positive = np.isfinite(scales) & (scales > 0)
    se[np.flatnonzero(live)[positive]] = scales[positive] * np.sqrt(nc) / n
    if positive.any():
        maxima = np.max(np.abs(selected[:, positive] / scales[positive]), axis=1)
        crit = order_quantile(maxima[:, None], .95)[0]
    else:
        crit = np.nan
    return selected, se, crit


def audit_bootstrap(mboot_module, ATTgt):
    rng = np.random.default_rng(202610262)
    rows = []
    original = mboot_module.run_multiplier_bootstrap
    try:
        for clustered in (False, True):
            for draws in (1, 2, 5, 199):
                for scale in (1., -1e-12, 1e-100, -1e-200, 1e6):
                    n = 64
                    inf = np.c_[rng.normal(size=(n, 3)), np.zeros(n), np.full(n, np.nan), np.full(n, np.inf)] * scale
                    clusters = np.arange(n) % 13 if clustered else None
                    multipliers = rng.normal(size=(draws, 13 if clustered else n))
                    mboot_module.run_multiplier_bootstrap = lambda x, b, pl, cores: multipliers @ x / len(x)
                    params = dict(data=pd.DataFrame(dict(id=np.arange(n), period=1, cluster=np.arange(n) % 13)),
                                  idname="id", tname="period", clustervars="cluster" if clustered else None,
                                  biters=draws, alp=.05, panel=True, true_repeated_cross_sections=False)
                    with np.errstate(all="ignore"):
                        expected_bres, expected_se, expected_crit = bootstrap_oracle(inf, multipliers, clusters)
                        out = mboot_module.mboot(inf, params)
                    fields = dict(bres=metrics(out["bres"], expected_bres, abs(scale)),
                                  se=metrics(out["se"], expected_se, abs(scale)),
                                  crit=metrics(out["crit_val"], expected_crit))
                    rows.append(record("frozen_multiplier_order_statistics",
                                       dict(clustered=clustered, biters=draws, scale=scale), fields))
    finally:
        mboot_module.run_multiplier_bootstrap = original
    # Exact degeneracy of built-in methods, using public fit and all aggregations.
    template = make_panel("regular")
    for faster in (False, True):
        for method in ("dr", "reg", "ipw"):
            for panel in (False, True):
                data = template.assign(Y=0.)
                if not panel:
                    data["id"] = np.arange(len(data))
                fit = ATTgt("Y", "period", "id", "G", data=data, xformla="~1", panel=panel,
                            faster_mode=faster, cband=True, biters=99).fit(est_method=method, bstrap=True)
                fields = dict(att=metrics(fit.MP["att"], np.zeros(len(fit.MP["att"]))),
                              inffunc=metrics(fit.MP["inffunc"]["inffunc"],
                                              np.zeros_like(fit.MP["inffunc"]["inffunc"])),
                              se=metrics(fit.results["se"], np.full(len(fit.MP["att"]), np.nan)),
                              crit=metrics(fit.results["c"], norm.ppf(.975)),
                              cband_flag=dict(passed=fit.dp["cband"] is False),
                              W=metrics(fit.MP["W"], None))
                for kind in ("simple", "dynamic", "group", "calendar"):
                    out = fit.aggte(typec=kind, cband=True).atte
                    fields[kind+"_overall_se"] = metrics(out["overall_se"], np.full_like(out["overall_se"], np.nan))
                    fields[kind+"_att"] = metrics(out["overall_att"], 0.)
                    for component, values in out["inf_function"].items():
                        fields[kind+"_IF_"+component] = metrics(values, np.zeros_like(values))
                    if kind != "simple":
                        fields[kind+"_component_att"] = metrics(out["att_egt"], np.zeros_like(out["att_egt"]))
                        fields[kind+"_component_se"] = metrics(out["se_egt"], np.full_like(out["se_egt"], np.nan))
                        fields[kind+"_critical"] = metrics(out["crit_val_egt"], norm.ppf(.975))
                        fields[kind+"_cband"] = dict(passed=out["DIDparams"]["cband"] is False)
                rows.append(record("public_zero_outcome_exact_degeneracy",
                                   dict(faster_mode=faster, method=method, panel=panel), fields))
    return rows


def worker(source, output, balance_count):
    source = source.resolve()
    sys.path.insert(0, str(source))
    att_module = importlib.import_module("csdid.att_gt")
    aggregate_module = importlib.import_module("csdid.aggte_fnc.compute_aggte")
    mboot_module = importlib.import_module("csdid.utils.mboot")
    modules = (att_module, aggregate_module, mboot_module)
    for module in modules:
        assert Path(module.__file__).resolve().is_relative_to(source), module.__file__
    rows = []
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rows += audit_panels(att_module.ATTgt)
        rows += audit_balance(aggregate_module.compute_aggte, balance_count)
        rows += audit_bootstrap(mboot_module, att_module.ATTgt)
    counts = {}
    for row in rows:
        counts.setdefault(row["name"], {"cases": 0, "passed": 0, "failed": 0})
        counts[row["name"]]["cases"] += 1
        counts[row["name"]]["passed" if row["passed"] else "failed"] += 1
    payload = dict(source=str(source), commit=git(source, "rev-parse", "HEAD"),
                   tracked_status=git(source, "status", "--porcelain", "--untracked-files=no"),
                   verified_imports={m.__name__:str(Path(m.__file__).resolve()) for m in modules},
                   source_sha256={name:digest(source / name) for name in (
                       "csdid/att_gt.py", "csdid/aggte_fnc/compute_aggte.py", "csdid/aggte_fnc/utils.py",
                       "csdid/utils/mboot.py", "csdid/attgt_fnc/compute_att_gt_shared.py",
                       "csdid/attgt_fnc/compute_att_gt.py", "csdid/attgt_fnc/compute_att_gt2.py",
                       "csdid/attgt_fnc/preprocess_did.py")},
                   input_panel_sha256={clock:hashlib.sha256(make_panel(clock).to_csv(index=False,
                       float_format="%.17g").encode()).hexdigest() for clock in
                       ("regular", "irregular_negative", "fractional")},
                   input_seeds=[180226, 202610261, 202610262], counts=counts,
                   warning_counts={name:sum(type(w.message).__name__ == name for w in caught)
                                   for name in sorted({type(w.message).__name__ for w in caught})},
                   rows=rows)
    output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--balance-count", type=int, default=1000)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.output, args.balance_count)
        return 0
    if args.baseline is None or args.candidate is None:
        parser.error("--baseline and --candidate required")
    results = {}
    env = dict(os.environ, **{key:"1" for key in THREADS})
    env.pop("PYTHONPATH", None)
    for label, source in (("before", args.baseline), ("after", args.candidate)):
        target = args.output.with_name(args.output.stem + f"_{label}.json")
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(source.resolve()),
                        "--output", str(target.resolve()), "--balance-count", str(args.balance_count)],
                       cwd="/tmp", env=env, check=True)
        results[label] = json.loads(target.read_text())
        target.unlink()
    before, after = results["before"], results["after"]
    assert len(before["rows"]) == len(after["rows"])
    transitions = {"corrected": [], "preserved_correct": [], "regressed": [], "still_failing": []}
    for old, new in zip(before["rows"], after["rows"]):
        assert old["name"] == new["name"] and old["parameters"] == new["parameters"]
        category = ("preserved_correct" if old["passed"] else "corrected") if new["passed"] else (
            "regressed" if old["passed"] else "still_failing")
        transitions[category].append(dict(name=new["name"], parameters=new["parameters"]))
    failure_kinds = {}
    for old in before["rows"]:
        if old["passed"]:
            continue
        if old["name"] == "dynamic_ratio_equations" and old.get("expected_error"):
            label = ("diagnostic_only_refusal_message" if old.get("exception") == "ValueError"
                     else "returned_aggregate_when_refusal_required")
        elif old["name"] == "dynamic_ratio_equations":
            label = "aggregate_selection_numerical_or_scale_error"
        elif old["name"] == "public_panel_equations":
            label = "public_ATT_IF_SE_interval_or_Wald_equation_error"
        elif old["name"] == "frozen_multiplier_order_statistics":
            label = "bootstrap_dimension_scale_or_critical_error"
        else:
            label = "unavailable_scale_critical_value_or_fallback_error"
        failure_kinds[label] = failure_kinds.get(label, 0) + 1
    artifact = dict(schema_version=1, audit_date="2026-10-26", policy=POLICY,
                    oracle_corrections_after_results=ORACLE_CORRECTIONS,
                    baseline_failure_classification=failure_kinds,
                    oracle_scope=__doc__, harness_sha256=digest(__file__), python=sys.version,
                    versions=dict(numpy=np.__version__, pandas=pd.__version__),
                    counts={k:len(v) for k,v in transitions.items()},
                    numerical_and_decision_counts=dict(
                        total_cases=len(after["rows"]),
                        baseline_correct=len(transitions["preserved_correct"]) + failure_kinds.get("diagnostic_only_refusal_message", 0),
                        candidate_correct=sum(row["passed"] for row in after["rows"]),
                        corrected_excluding_diagnostic_only=len(transitions["corrected"]) - failure_kinds.get("diagnostic_only_refusal_message", 0),
                        diagnostic_only_message_changes=failure_kinds.get("diagnostic_only_refusal_message", 0),
                        regressions=len(transitions["regressed"])),
                    transitions=transitions, sources=results)
    serialized = (json.dumps(artifact, indent=2, allow_nan=False) + "\n").encode()
    if args.output.suffix == ".gz":
        args.output.write_bytes(gzip.compress(serialized, mtime=0))
    else:
        args.output.write_bytes(serialized)
    if args.summary:
        summary = {key:artifact[key] for key in ("schema_version", "audit_date", "policy",
            "oracle_corrections_after_results", "baseline_failure_classification", "harness_sha256", "counts", "numerical_and_decision_counts")}
        summary["source_identity"] = {key:{field:value[field] for field in (
            "commit", "tracked_status", "verified_imports", "source_sha256", "input_panel_sha256", "counts")}
            for key,value in results.items()}
        summary["details"] = dict(path=str(args.output), sha256=digest(args.output),
                                  uncompressed_sha256=hashlib.sha256(serialized).hexdigest())
        summary["limits"] = [
            "The baseline diagnostic-only refusal messages are separate from numerical or support failures.",
            "These seeded analytical identities are not Monte Carlo coverage or proof of identifying assumptions.",
            "An unavailable overall target has no statistical IF; its internal placeholder is not compared.",
            "The corrected oracle and unchanged numerical thresholds are documented after the exploratory run.",
            "Outcome-unit-normalized analytic identities are separate from the permanent raw-scale R parity gate; the retained large-scale cancellation misses are not cleared."]
        args.summary.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(counts=artifact["counts"], before=before["counts"], after=after["counts"]), indent=2))
    return 1 if transitions["regressed"] or transitions["still_failing"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
