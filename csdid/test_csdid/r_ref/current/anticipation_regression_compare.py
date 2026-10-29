"""Independent anticipation oracle and before/after regression audit.

Use the same Python environment for clean baseline/candidate checkouts::

    python anticipation_regression_compare.py --baseline /clean/before \
      --candidate /clean/after --output /tmp/csdid-anticipation-proof

This compares complete public and influence-function output, exact support and
nonfinite masks, caller-data preservation, and NumPy RNG state. Anticipation
estimates are checked separately against direct weighted differences of means
and their analytic influence functions; this oracle does not call R, DRDID, or
the package's cell-selection/estimation helpers. Agreement with old code alone
is regression evidence, not proof that either implementation is correct.
"""
import argparse
import contextlib
import gzip
import hashlib
import importlib.metadata
import io
import itertools
import json
import os
from pathlib import Path
import pickle
import platform
import subprocess
import sys
import time
import warnings

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")

import numpy as np
import pandas as pd
from scipy.stats import chi2


ATOL, RTOL = 1e-10, 1e-9
SHA = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def cases():
    rows = []
    for design, method, fast, base, control, cov, sample in itertools.product(
        ("panel", "rc", "unbalanced"), ("dr", "reg", "ipw"), (False, True),
        ("varying", "universal"), ("nevertreated", "notyettreated"),
        ("invariant", "varying"), ("clean", "missing_y", "missing_row", "no_never")
    ):
        rows.append(dict(kind="regression", design=design, method=method, fast=fast,
            base=base, control=control, cov=cov, sample=sample,
            weighted=bool(len(rows) % 2), clustered=bool((len(rows) // 2) % 2), anticipation=0))
    assert len(rows) == 576
    for method, fast, base, control, weighted, anticipation in itertools.product(
        ("dr", "reg", "ipw"), (False, True), ("varying", "universal"),
        ("nevertreated", "notyettreated"), (False, True), (1, 2)
    ):
        rows.append(dict(kind="anticipation_oracle", design="panel", method=method,
            fast=fast, base=base, control=control, cov="intercept", sample="clean",
            weighted=weighted, clustered=False, anticipation=anticipation))
    assert len(rows) == 672
    for design, fast in itertools.product(("panel", "rc", "unbalanced"), (False, True)):
        rows.append(dict(kind="bootstrap_rng", design=design, method="dr", fast=fast,
            base="varying", control="notyettreated", cov="varying", sample="missing_y",
            weighted=True, clustered=True, anticipation=0))
    for i, row in enumerate(rows):
        row["name"] = f"{i:03d}_{row['kind']}_{row['design']}_{row['method']}_{int(row['fast'])}_{row['base']}_{row['control']}_{row['cov']}_{row['sample']}_a{row['anticipation']}"
    return rows


def data_for(z):
    rng = np.random.default_rng(260106)
    n, periods = 120, np.arange(1, 11)
    uid = np.arange(101, 101 + n)
    cohort = np.repeat([0, 5, 7, 9], n // 4)
    alpha, x = rng.normal(size=n), rng.normal(size=n)
    frames = []
    for t in periods:
        xt = x + (.1 * np.cos(uid / 7 + t) if z["cov"] == "varying" else 0)
        y = alpha + .2 * t + .35 * xt + .07 * t * x + rng.normal(scale=.6, size=n)
        y += .8 * ((cohort > 0) & (t >= cohort - 2) & (t < cohort))
        y += np.where((cohort > 0) & (t >= cohort), 1.5 + .1 * (t - cohort), 0)
        frames.append(pd.DataFrame(dict(id=uid, period=t, G=cohort, Y=y, X=xt,
            weight=np.exp(.13 * x + .15 * np.sin(uid / 9 + t)), cluster=uid % 19)))
    data = pd.concat(frames, ignore_index=True)
    if z["design"] == "unbalanced":
        data = data[(data.id * 7 + data.period * 3) % 23 != 0]
    if z["sample"] == "missing_y":
        data.loc[(data.id == 104) & (data.period == 3), "Y"] = np.nan
    elif z["sample"] == "missing_row":
        data = data[~((data.id == 136) & (data.period == 4))]
    elif z["sample"] == "no_never":
        data = data[data.G != 0]
    data = data.sample(frac=1, random_state=1109).reset_index(drop=True)
    if z["design"] == "rc":
        data["id"] = np.arange(100001, 100001 + len(data))
    return data


def pack(value):
    if value is None:
        return None
    values = np.asarray(value, dtype=float)
    flat = values.reshape(-1).copy()
    masks = {name: np.flatnonzero(fn(flat)).tolist() for name, fn in (
        ("nan", np.isnan), ("posinf", np.isposinf), ("neginf", np.isneginf))}
    flat[~np.isfinite(flat)] = 0
    return dict(shape=list(values.shape), values=flat.tolist(), **masks)


def unpack(value):
    if value is None:
        return None
    out = np.asarray(value["values"], dtype=float)
    for key, val in (("nan", np.nan), ("posinf", np.inf), ("neginf", -np.inf)):
        out[value[key]] = val
    return out.reshape(value["shape"])


def rng_hash():
    return hashlib.sha256(pickle.dumps(np.random.get_state())).hexdigest()


def unit_ids(obj):
    dp = obj.MP["DIDparams"]
    internal, data = dp["idname"], dp["data"]
    ids = data[internal].unique()
    if internal != "id":
        mapping = data[[internal, "id"]].drop_duplicates()
        assert mapping[internal].is_unique
        ids = mapping.set_index(internal)["id"].reindex(ids).to_numpy()
    assert len(ids) == len(set(ids))
    return ids


def evaluate(ATTgt, z):
    data = data_for(z)
    untouched = data.copy(deep=True)
    np.random.seed(998124)
    initial_rng = rng_hash()
    bootstrap = z["kind"] == "bootstrap_rng"
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        obj = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data,
            xformla="~ 1" if z["cov"] == "intercept" else "~ X", faster_mode=z["fast"],
            panel=z["design"] != "rc", allow_unbalanced_panel=z["design"] == "unbalanced",
            weights_name="weight" if z["weighted"] else None,
            clustervar="cluster" if z["clustered"] else None,
            control_group=z["control"], anticipation=z["anticipation"],
            cband=bootstrap, biters=63).fit(est_method=z["method"],
                base_period=z["base"], bstrap=bootstrap)
        ids = unit_ids(obj)
        order = np.argsort(ids)
        fit = {key: pack(obj.results[key]) for key in ("att", "se", "c", "l_se", "u_se")}
        fit.update({key: pack(obj.MP[key]) for key in ("group", "t", "n", "W", "Wpval")})
        fit["inffunc"] = pack(np.asarray(obj.MP["inffunc"]["inffunc"])[order])
        fit["unit_ids"] = pack(ids[order])
        fit["cband"] = pack(obj.dp["cband"])
        after_fit_rng = rng_hash()
        aggs = {}
        for kind in ("simple", "group", "dynamic", "calendar"):
            try:
                obj.aggte(typec=kind, bstrap=bootstrap, cband=bootstrap, na_rm=True)
                aggs[kind] = {key: pack(obj.atte[key]) for key in (
                    "overall_att", "overall_se", "egt", "att_egt", "se_egt", "crit_val_egt")}
                aggs[kind]["inf_function"] = {
                    key: pack(np.asarray(value)[order])
                    for key, value in obj.atte["inf_function"].items()}
            except Exception as exc:
                aggs[kind] = {"error": type(exc).__name__ + ": " + str(exc)}
        after_aggregation_rng = rng_hash()
    pd.testing.assert_frame_equal(data, untouched, check_exact=True)
    return dict(fit=fit, aggs=aggs, caller_data_unchanged=True,
        rng=dict(initial=initial_rng, after_fit=after_fit_rng, after_aggregation=after_aggregation_rng),
        warnings=[str(w.message) for w in caught])


def worker(repo, output, manifest):
    sys.path.insert(0, str(repo))
    from csdid.att_gt import ATTgt
    loaded = Path(sys.modules[ATTgt.__module__].__file__).resolve()
    assert loaded.is_relative_to(repo.resolve()), loaded
    output.mkdir(parents=True, exist_ok=True)
    rows = json.loads(manifest.read_text())
    for i, z in enumerate(rows):
        try:
            result = evaluate(ATTgt, z)
        except Exception as exc:
            result = {"error": type(exc).__name__ + ": " + str(exc)}
        with gzip.open(output / (z["name"] + ".json.gz"), "wt") as stream:
            json.dump(result, stream)
        if i % 48 == 0:
            print(f"{repo.name}: {i}/{len(rows)}", flush=True)


def known_effect_worker(repo, output):
    """No-noise counterexample: known anticipation=.8 and treatment=1.5."""
    sys.path.insert(0, str(repo))
    from csdid.att_gt import ATTgt
    ids = np.arange(60)
    frames = []
    for t in range(1, 9):
        g = np.where(ids < 30, 0, 5)
        y = np.sin(ids / 3) + .2 * t + .8 * ((g == 5) & (t >= 3) & (t < 5))
        y += 1.5 * ((g == 5) & (t >= 5))
        frames.append(pd.DataFrame(dict(id=ids, period=t, G=g, Y=y,
                                       weight=np.exp(.1 * np.cos(ids + t)))))
    data = pd.concat(frames, ignore_index=True)
    records = []
    for method, fast, base in itertools.product(("dr", "reg", "ipw"),
                                               (False, True), ("varying", "universal")):
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            fit = ATTgt(yname="Y", tname="period", idname="id", gname="G",
                data=data, faster_mode=fast, weights_name="weight", anticipation=2).fit(
                    est_method=method, base_period=base, bstrap=False)
        times, att = np.asarray(fit.MP["t"]), np.asarray(fit.MP["att"])
        expected = np.where(times >= 5, 1.5, np.where(times >= 3, .8, 0))
        records.append(dict(method=method, fast=fast, base=base, t=times.tolist(),
            att=att.tolist(), expected=expected.tolist(),
            max_abs_error=float(np.max(np.abs(att-expected))),
            passed=bool(np.allclose(att, expected, atol=1e-12, rtol=1e-12))))
    output.write_text(json.dumps(records, indent=2) + "\n")


def compare_numeric(a, b, field, issues, metrics):
    if a is None or b is None:
        if not (a is None and b is None):
            issues.append(field + ": absent field differs")
        return
    aa, bb = unpack(a), unpack(b)
    if aa.shape != bb.shape:
        issues.append(field + ": shape differs")
        return
    for label, mask in (("nan", np.isnan), ("posinf", np.isposinf), ("neginf", np.isneginf)):
        if not np.array_equal(mask(aa), mask(bb)):
            issues.append(field + ": " + label + " mask differs")
            return
    finite = np.isfinite(aa)
    metrics["elements"] += int(aa.size)
    difference = float(np.max(np.abs(aa[finite] - bb[finite]), initial=0))
    metrics["max_abs"] = max(metrics["max_abs"], difference)
    if not np.allclose(aa[finite], bb[finite], atol=ATOL, rtol=RTOL):
        issues.append(f"{field}: numeric difference max={difference:.12g}")


def compare_tree(a, b, path, issues, metrics):
    if isinstance(a, dict) and "shape" in a:
        compare_numeric(a, b, path, issues, metrics)
    elif isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            issues.append(path + ": fields differ")
            return
        for key in a:
            if key != "warnings":
                compare_tree(a[key], b[key], path + "." + key, issues, metrics)
    elif a is None or b is None:
        if a != b:
            issues.append(path + ": absent field differs")
    elif a != b:
        issues.append(path + ": value differs")


def oracle(z, result):
    data = data_for(z)
    ids = unpack(result["fit"]["unit_ids"])
    groups, times = (unpack(result["fit"][k]) for k in ("group", "t"))
    periods = np.sort(data.period.unique())
    estimates, influence = [], []
    for g, t in zip(groups, times):
        base = (max(p for p in periods if p < g - z["anticipation"])
                if z["base"] == "universal" or t >= g - z["anticipation"]
                else max(p for p in periods if p < t))
        before = data[data.period == base].set_index("id").reindex(ids)
        after = data[data.period == t].set_index("id").reindex(ids)
        reference = data[data.period == min(t, base)].set_index("id").reindex(ids)
        weights = reference.weight.to_numpy() if z["weighted"] else np.ones(len(ids))
        cohort = before.G.to_numpy()
        treated = cohort == g
        control = cohort == 0
        if z["control"] == "notyettreated":
            control |= (cohort > max(t, base) + z["anticipation"]) & ~treated
        delta = (after.Y - before.Y).to_numpy()
        mean_t, mean_c = (np.average(delta[m], weights=weights[m]) for m in (treated, control))
        effect = mean_t - mean_c
        inf = np.zeros(len(ids))
        for mask, mean, sign in ((treated, mean_t, 1), (control, mean_c, -1)):
            inf[mask] = sign * len(ids) * weights[mask] * (delta[mask] - mean) / weights[mask].sum()
        estimates.append(effect)
        influence.append(inf)
    estimates, influence = np.asarray(estimates), np.asarray(influence).T
    se = np.sqrt(np.mean(influence ** 2, axis=0) / len(ids))
    se[se <= 0] = np.nan
    pre = (times < groups - z["anticipation"]) & np.isfinite(se)
    pre_att, pre_if = estimates[pre], influence[:, pre]
    cov = pre_if.T @ pre_if / len(ids) ** 2
    wald = float(pre_att @ np.linalg.solve(cov, pre_att))
    # The public API reports five decimal places; preserve that presentation
    # contract rather than comparing its rounded value to an unrounded oracle.
    wp = round(float(chi2.sf(wald, int(pre.sum()))), 5)
    return dict(att=pack(estimates), se=pack(se), inffunc=pack(influence),
                W=pack(wald), Wpval=pack(wp))


def summarize(output, manifest, baseline, candidate):
    checks, expected_old_oracle_failures = [], []
    for z in json.loads(manifest.read_text()):
        def read(label):
            with gzip.open(output / label / (z["name"] + ".json.gz"), "rt") as stream:
                return json.load(stream)
        old, new = read("baseline"), read("candidate")
        issues, metrics = [], dict(elements=0, max_abs=0.)
        row = dict(case=z["name"], kind=z["kind"])
        if "error" in old or "error" in new:
            issues.append("A fit failed: " + repr((old.get("error"), new.get("error"))))
        elif z["kind"] != "anticipation_oracle":
            compare_tree(old, new, "output", issues, metrics)
            if z["kind"] == "regression" and len(set(new["rng"].values())) != 1:
                issues.append("Analytic fitting/aggregation changed NumPy RNG state")
        else:
            expected = oracle(z, new)
            for field in expected:
                compare_numeric(new["fit"][field], expected[field], "direct_oracle." + field, issues, metrics)
            group, t = (unpack(new["fit"][k]) for k in ("group", "t"))
            unaffected = (t < group - z["anticipation"]) | (t >= group)
            compare_numeric(old["fit"]["unit_ids"], new["fit"]["unit_ids"], "unaffected.unit_ids", issues, metrics)
            for field in ("group", "t", "att", "se", "l_se", "u_se"):
                compare_numeric(pack(unpack(old["fit"][field])[unaffected]),
                                pack(unpack(new["fit"][field])[unaffected]),
                                "unaffected." + field, issues, metrics)
            compare_numeric(pack(unpack(old["fit"]["inffunc"])[:, unaffected]),
                            pack(unpack(new["fit"]["inffunc"])[:, unaffected]),
                            "unaffected.inffunc", issues, metrics)
            for kind in ("simple", "group", "calendar"):
                compare_tree(old["aggs"][kind], new["aggs"][kind], "unaffected." + kind, issues, metrics)
            old_issues, old_metrics = [], dict(elements=0, max_abs=0.)
            for field in expected:
                compare_numeric(old["fit"][field], expected[field], "old_direct_oracle." + field, old_issues, old_metrics)
            expected_old_oracle_failures.append(dict(case=z["name"], issues=old_issues,
                changed_window_att=float(np.max(np.abs(unpack(old["fit"]["att"]) - unpack(new["fit"]["att"])))),
                old_wald=unpack(old["fit"]["W"]).tolist(), new_wald=unpack(new["fit"]["W"]).tolist()))
        row.update(passed=not issues, issues=issues, **metrics)
        checks.append(row)
    def source(path):
        commit = subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
        status = subprocess.check_output(["git", "-C", str(path), "status", "--porcelain", "--untracked-files=no"], text=True)
        assert not status, str(path) + ": dirty source"
        files = subprocess.check_output(["git", "-C", str(path), "ls-files", "csdid"], text=True).splitlines()
        return dict(commit=commit, tracked_source_clean=True,
            runtime_files_sha256={name: SHA(path/name) for name in files if name.endswith(".py") and "/test_csdid/" not in name})
    result = dict(audit_date="2026-10-26", baseline=source(baseline), candidate=source(candidate),
        tolerance=dict(atol=ATOL, rtol=RTOL, support_and_nonfinite_masks="exact"),
        runtime=dict(python=sys.version, platform=platform.platform(), packages={name:importlib.metadata.version(name)
            for name in ("numpy", "pandas", "scipy", "statsmodels", "drdid")}),
        script_sha256=SHA(__file__), case_manifest_sha256=SHA(manifest), cases=json.loads(manifest.read_text()),
        oracle_presentation=dict(Wpval_decimal_places=5),
        audit_corrections=["The first independent oracle emitted unrounded chi-square p-values and missed 48 anticipation=2 public Wpval comparisons by at most 5e-6. ATT/SE/full IF/Wald statistics already passed. The oracle now applies the existing five-decimal public presentation contract; no estimator, input or tolerance changed."],
        checks=checks, previous_version_direct_oracle_failures=expected_old_oracle_failures,
        compared_elements=sum(r["elements"] for r in checks),
        passed=sum(r["passed"] for r in checks), total=len(checks),
        limitations=["A finite differential grid cannot prove absence of all regressions.",
            "The independent oracle covers intercept-only balanced panels; covariate-rich, RC and unbalanced cases use before/after regression comparisons.",
            "Warning wording/counts are recorded in raw outputs but excluded from equality because intentional diagnostic changes differ.",
            "Nonzero anticipation changes the window and Wald restrictions; dynamic aggregation is excluded from before/after equality there."])
    result["known_effect_counterexample"] = {
        label: json.loads((output / (label + "_known_effect.json")).read_text())
        for label in ("baseline", "candidate")}
    assert all(row["passed"] for row in result["known_effect_counterexample"]["candidate"])
    assert sum(not row["passed"] for row in result["known_effect_counterexample"]["baseline"]) == 6
    with gzip.open(output / "summary.json.gz", "wt") as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps({k:result[k] for k in ("passed", "total", "compared_elements")}, indent=2))
    for row in checks:
        if not row["passed"]:
            print(row["case"], row["issues"])
    return all(r["passed"] for r in checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--known-effect-worker", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    if args.worker:
        if args.known_effect_worker:
            known_effect_worker(args.worker.resolve(), output)
        else:
            worker(args.worker.resolve(), output, args.manifest.resolve())
        return 0
    if args.baseline is None or args.candidate is None:
        parser.error("--baseline and --candidate are required")
    output.mkdir(parents=True, exist_ok=True)
    manifest = output / "cases.json"
    if not args.summarize_only:
        manifest.write_text(json.dumps(cases(), indent=2) + "\n")
        for label, repo in (("baseline", args.baseline), ("candidate", args.candidate)):
            with (output / (label + ".log")).open("w") as log:
                subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(repo.resolve()),
                    "--manifest", str(manifest), "--output", str(output/label)], stdout=log,
                    stderr=subprocess.STDOUT, check=True)
    for label, repo in (("baseline", args.baseline), ("candidate", args.candidate)):
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(repo.resolve()),
            "--known-effect-worker", "--output", str(output / (label + "_known_effect.json"))], check=True)
    return 0 if summarize(output, manifest, args.baseline.resolve(), args.candidate.resolve()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
