"""Independent and before/after evidence for the October 2026 support fixes.

Run with the project's test interpreter (uses pure base-R estimator oracles)::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python \\
      csdid/test_csdid/validate_support_adversarial.py \\
      --baseline /path/to/2a36b48 --candidate /path/to/81fa180 --output /tmp/support-audit

No current-only helper is imported into the old process. The healthy grid
compares complete public estimates, inference, all aggregate influence functions,
warnings and RNG state by exact canonical hashes. Targeted records retain full
arrays locally. Oracle tolerances are declared here before comparison.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import warnings

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(variable, "1")

import numpy as np
import pandas as pd

PINS = {"baseline": "2a36b48e2c435d87c6dfaecda36734a4336c78ec",
        "candidate": "81fa18073a239f40173c9a9a60e5c5d4934073be"}
ATOL, RTOL = 1e-7, 1e-6


def canonical(value):
    if isinstance(value, dict):
        return {k: canonical(v) for k, v in sorted(value.items())}
    if isinstance(value, np.ndarray):
        return {"shape": list(value.shape), "values": canonical(value.tolist())}
    if isinstance(value, (list, tuple)):
        return [canonical(v) for v in value]
    if isinstance(value, np.generic):
        return canonical(value.item())
    return value


def fingerprint(value):
    return hashlib.sha256(json.dumps(canonical(value), sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def source_info(repo):
    sha = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    paths = sorted(p for p in (repo / "csdid").rglob("*.py")
                   if "test_csdid" not in p.parts)
    changed = subprocess.check_output(["git", "-C", str(repo), "diff", "--name-only", "HEAD", "--", "csdid"], text=True).splitlines()
    runtime_changed = [name for name in changed if not name.startswith("csdid/test_csdid/")]
    assert not runtime_changed, "Uncommitted runtime changes: " + str(runtime_changed)
    return {"commit": sha, "tracked_runtime_clean": True, "runtime_sha256": {
        str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def panel_data():
    rng = np.random.default_rng(27827931)
    n, periods = 160, 5
    ids = np.repeat(np.arange(n), periods)
    t = np.tile(np.arange(1, periods + 1), n)
    groups = np.repeat(np.repeat([0, 3, 4, 5], n // 4), periods)
    x = np.repeat(rng.uniform(-1., 1., n), periods) + .06 * np.sin(ids + t)
    y = (.2 * t + .3 * x + .1 * t * x + np.repeat(rng.normal(size=n), periods)
         + rng.normal(scale=.4, size=len(ids)) + ((groups > 0) & (t >= groups)) * (1. + .1 * x))
    weight = np.exp(.2 * x + .1 * np.sin(ids + .4 * t))
    return pd.DataFrame(dict(id=ids, period=t, G=groups, X=x, Y=y, weight=weight))


def tiny_rc(n_treated):
    cube = np.array(list(itertools.product([-1., 1.], repeat=3)))
    regular = np.tile(cube, (2, 1))
    tiny = (np.array([[-1., -1., -1.], [1., 1., -1.], [0., 0., 1.]])
            if n_treated == 3 else np.vstack([np.zeros(3), np.eye(3)]))
    frames = []
    for period, group in itertools.product(range(1, 4), (0, 2, 3)):
        x = tiny if group == 2 else regular
        frame = pd.DataFrame(x, columns=["X", "X2", "X3"])
        frame["period"], frame["G"], frame["weight"] = period, group, 1.
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    data["id"] = np.arange(len(data))
    data["Y"] = (.2 * data.period + .3 * data.X - .1 * data.X2 + .2 * data.X3
                 + ((data.G > 0) & (data.period >= data.G)) + .1 * np.sin(data.id))
    return data


def prepare(output):
    data = panel_data()
    data.to_csv(output / "panel.csv", index=False, float_format="%.17g")
    rc = data.copy()
    rc["id"] = np.arange(len(rc))
    rc.to_csv(output / "rc.csv", index=False, float_format="%.17g")
    unbalanced = data[(data.id * 7 + data.period * 3) % 13 != 0]
    unbalanced.to_csv(output / "unbalanced.csv", index=False, float_format="%.17g")
    healthy = []
    for route in ("panel", "rc", "unbalanced"):
        anchors = (None,) if route == "rc" else (None, "first_period", "base_period", "varying")
        for anchor, method, control, base, profile, fast in itertools.product(
                anchors, ("dr", "reg", "ipw"), ("nevertreated", "notyettreated"),
                ("varying", "universal"), ("none", "unit", "positive", "timezero"), (False, True)):
            healthy.append(dict(name=f"healthy_{len(healthy):04d}", route=route,
                data=route + ".csv", anchor=anchor, method=method, control=control,
                base=base, profile=profile, fast=fast, formula="~X"))
    assert len(healthy) == 864
    targets = []
    def add(name, filename, method, route="panel", formula="~X", oracle=True, **extra):
        for fast in (False, True):
            targets.append(dict(name=f"{name}_{str(fast).lower()}", route=route,
                data=filename, anchor=None, method=method, control="nevertreated",
                base="varying", profile="positive", fast=fast, formula=formula,
                oracle=oracle, **extra))
    extra = pd.DataFrame(dict(id=160, period=np.arange(1, 6), G=3, X=1000., Y=0., weight=0.))
    pd.concat([data, extra], ignore_index=True).to_csv(output / "zero_treated.csv", index=False, float_format="%.17g")
    for method in ("dr", "ipw"):
        add("zero_treated_" + method, "zero_treated.csv", method)
    for mass, name in ((0., "zero_control"), (1e-24, "tiny_control")):
        extra = extra.assign(G=0, X=1e12, weight=mass)
        pd.concat([data, extra], ignore_index=True).to_csv(output / (name + ".csv"), index=False, float_format="%.17g")
        add(name + "_reg", name + ".csv", "reg")
    # Intercept-only closed-form means/IF do not fit nuisance models.
    for route, zero in itertools.product(("panel", "rc"), (False, True)):
        sample = data.copy() if route == "panel" else rc.copy()
        if zero:
            sample.loc[(sample.id % 9 == 0), "weight"] = 0.
        filename = f"mean_{route}_{zero}.csv"
        sample.to_csv(output / filename, index=False, float_format="%.17g")
        for method in ("dr", "reg", "ipw"):
            add(f"mean_{route}_{zero}_{method}", filename, method, route, "~1")
    # Positive raw weights disappear only in the 2x2 estimator's normalization.
    loss_panel = data[data.G.isin([0, 3, 4])].copy()
    loss_panel["G"] = loss_panel.G.replace({3: 2, 4: 3})
    loss_panel["weight"] = np.where(loss_panel.G == 2, 2. ** -100, 2. ** 1000)
    for route in ("rc", "varying", "unbalanced"):
        loss = loss_panel.copy()
        if route == "rc":
            loss["id"] = np.arange(len(loss))
        elif route == "unbalanced":
            loss = loss.iloc[:-1].copy()
        filename = f"normalization_loss_{route}.csv"
        loss.to_csv(output / filename, index=False, float_format="%.17g")
        for do_inf in (False, True):
            add(f"normalization_loss_{route}_inf{do_inf}", filename, "dr", route, "~1",
                oracle=False, safe_probe=True, do_inf=do_inf)
    for n_treated, method in itertools.product((3, 4), ("dr", "reg", "ipw")):
        filename = f"tiny_{n_treated}.csv"
        tiny_rc(n_treated).to_csv(output / filename, index=False, float_format="%.17g")
        add(f"tiny_{n_treated}_{method}", filename, method, "rc", "~X+X2+X3", oracle=False)
    extra = pd.DataFrame(dict(id=160, period=np.arange(1, 6), G=3, X=1000., Y=0.,
                              weight=[0., 1., 1., 1., 1.]))
    cache_base = data.copy()
    cache_base["X"] = cache_base.groupby("id")["X"].transform("first")
    pd.concat([cache_base, extra], ignore_index=True).to_csv(output / "cache_support.csv", index=False, float_format="%.17g")
    for method in ("dr", "ipw"):
        add("cache_support_" + method, "cache_support.csv", method, oracle=False)
        for original in list(targets[-2:]):
            targets.append(dict(original,
                name=f"cache_support_{method}_notyet_{str(original['fast']).lower()}", control="notyettreated"))
    (output / "cases.json").write_text(json.dumps(dict(healthy=healthy, targeted=targets), indent=2))
    return healthy, targets


def weight_profile(data, profile):
    data = data.copy()
    if profile == "unit":
        data["weight"] = 1.
    elif profile == "timezero":
        data.loc[(data.id + 3 * data.period) % 11 == 0, "weight"] = 0.
    return data


def snapshot_fit(fit):
    params = fit.MP["DIDparams"]
    frame, internal_id = params["data"], params["idname"]
    ids = frame[internal_id].unique()
    if internal_id != "id":
        mapping = frame[[internal_id, "id"]].drop_duplicates()
        ids = mapping.set_index(internal_id)["id"].reindex(ids).to_numpy()
    influence = fit.MP["inffunc"]
    return dict(group=fit.MP["group"], t=fit.MP["t"],
                results=fit.results, W=fit.MP["W"], Wpval=fit.MP["Wpval"],
                unit_ids=ids, influence=None if influence is None else influence["inffunc"])


def fit_case(ATTgt, data, case, aggregate=True):
    np.random.seed(640278)
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        fit = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=data.copy(),
            xformla=case["formula"], faster_mode=case["fast"], panel=case["route"] != "rc",
            allow_unbalanced_panel=case["route"] == "unbalanced",
            fix_weights="varying" if case["route"] == "varying" else case["anchor"],
            weights_name=None if case["profile"] == "none" else "weight",
            control_group=case["control"], compute_inffunc=case.get("do_inf", True)
        ).fit(est_method=case["method"], base_period=case["base"], bstrap=False)
        result = snapshot_fit(fit)
        if aggregate:
            result["aggregations"] = {}
            for kind in ("simple", "group", "dynamic", "calendar"):
                fit.aggte(typec=kind, bstrap=False)
                result["aggregations"][kind] = {key: value for key, value in fit.atte.items()
                                                  if key not in ("DIDparams", "call")}
    result["warnings"] = [(warning.category.__name__, str(warning.message)) for warning in caught]
    result["rng_state"] = np.random.get_state()
    return result


def pack_arrays(result):
    if isinstance(result, dict):
        return {k: pack_arrays(v) for k, v in result.items()}
    if isinstance(result, np.ndarray):
        return result.tolist()
    if isinstance(result, (tuple, list)):
        return [pack_arrays(v) for v in result]
    if isinstance(result, np.generic):
        return result.item()
    return result


def worker(args):
    sys.path.insert(0, str(args.repo.resolve()))
    from csdid.att_gt import ATTgt
    from csdid.attgt_fnc import drdid_trim
    cases = json.loads((args.output / "cases.json").read_text())
    saved_healthy = []
    if args.targets_only:
        saved = json.loads((args.output / f"{args.role}.json").read_text())
        assert saved["source"]["commit"] == PINS[args.role]
        saved_healthy = saved["healthy"]
    report = {"healthy": saved_healthy, "targeted": {}, "source": source_info(args.repo)}
    data_cache = {}
    for case in ([] if args.targets_only else cases["healthy"]):
        key = case["data"], case["profile"]
        if key not in data_cache:
            data_cache[key] = weight_profile(pd.read_csv(args.output / case["data"]), case["profile"])
        result = fit_case(ATTgt, data_cache[key], case)
        report["healthy"].append(dict(name=case["name"], fingerprint=fingerprint(result)))
    for case in cases["targeted"]:
        data = pd.read_csv(args.output / case["data"])
        entries = []
        original = drdid_trim.drdid_rc
        if case.get("safe_probe"):
            def safe_backend(y, post, D, covariates=None, i_weights=None):
                with np.errstate(under="ignore"):
                    normalized = i_weights / np.mean(i_weights)
                lost = bool(np.any((normalized <= 0) | ~np.isfinite(normalized)))
                entries.append(dict(lost=lost, treated_rows=int(np.sum(D))))
                if lost:
                    raise ValueError("safe probe: normalized support loss reached backend")
                return original(y=y, post=post, D=D, covariates=covariates, i_weights=i_weights)
            drdid_trim.drdid_rc = safe_backend
        try:
            result = fit_case(ATTgt, data, case, aggregate=False)
        finally:
            drdid_trim.drdid_rc = original
        if case.get("safe_probe"):
            result["backend_entries"] = entries
        report["targeted"][case["name"]] = pack_arrays(result)
    if args.role == "candidate":
        # Seed a narrowly scoped defect without changing repository files.
        import csdid.att_gt as att_module
        source = (args.repo / "csdid/attgt_fnc/compute_att_gt2.py").read_text()
        original_key = "nuisance_cache['weight_keys'][weight_idx])"
        assert source.count(original_key) == 2
        # Only overlap_key loses the weight key; the propensity-fit key stays correct.
        mutant_source = source.replace(original_key, "None)", 1)
        module = type(sys)("_support_cache_mutant")
        exec(compile(mutant_source, "_support_cache_mutant.py", "exec"), module.__dict__)
        original_fast = att_module.compute_att_gt2
        att_module.compute_att_gt2 = module.compute_att_gt2
        try:
            for case in cases["targeted"]:
                if case["name"].startswith("cache_support") and case["fast"]:
                    result = fit_case(ATTgt, pd.read_csv(args.output / case["data"]), case, aggregate=False)
                    report["targeted"][case["name"] + "_weight_key_mutant"] = pack_arrays(result)
        finally:
            att_module.compute_att_gt2 = original_fast
    (args.output / f"{args.role}.json").write_text(json.dumps(report))
    print(args.role, "healthy fits", len(report["healthy"]), "targeted", len(report["targeted"]), flush=True)


def compare_numeric(actual, expected, label):
    actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    assert actual.shape == expected.shape, label + ": shape"
    np.testing.assert_array_equal(np.isnan(actual), np.isnan(expected), err_msg=label + ": NaN mask")
    np.testing.assert_array_equal(np.isposinf(actual), np.isposinf(expected), err_msg=label + ": +Inf mask")
    np.testing.assert_array_equal(np.isneginf(actual), np.isneginf(expected), err_msg=label + ": -Inf mask")
    np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=RTOL, err_msg=label)
    finite = np.isfinite(actual)
    return float(np.max(np.abs(actual[finite] - expected[finite]))) if finite.any() else 0.


def compare(args, healthy, targets):
    before, after = [json.loads((args.output / f"{role}.json").read_text()) for role in ("baseline", "candidate")]
    assert len(before["healthy"]) == len(after["healthy"]) == len(healthy)
    assert [z["name"] for z in before["healthy"]] == [z["name"] for z in healthy]
    assert [z["name"] for z in after["healthy"]] == [z["name"] for z in healthy]
    for role, report in (("baseline", before), ("candidate", after)):
        assert report["source"]["commit"] == PINS[role]
    mismatches = [z["name"] for z, a, b in zip(healthy, before["healthy"], after["healthy"])
                  if a["fingerprint"] != b["fingerprint"]]
    assert not mismatches, "Healthy full-output differences: " + str(mismatches)
    references = json.loads((args.output / "oracles.json").read_text())
    checks = []
    for case in targets:
        name = case["name"]
        old, new = before["targeted"][name], after["targeted"][name]
        old_att, new_att = np.asarray(old["results"]["att"]), np.asarray(new["results"]["att"])
        record = dict(case=name, old_finite_att=int(np.isfinite(old_att).sum()),
                      new_finite_att=int(np.isfinite(new_att).sum()))
        if case["oracle"]:
            ref = references[name]
            order = pd.Index(new["unit_ids"]).get_indexer(ref["unit_ids"])
            assert np.all(order >= 0) and len(order) == len(new["unit_ids"])
            assert list(new["group"]) == list(ref["group"]) and list(new["t"]) == list(ref["t"])
            record["oracle_att_max_abs"] = compare_numeric(new_att, ref["att"], name + " ATT")
            record["oracle_if_max_abs"] = compare_numeric(np.asarray(new["influence"])[order], ref["influence"], name + " IF")
            assert np.isfinite(new_att).all()
            if name.startswith(("zero_control", "tiny_control", "zero_treated")):
                assert new_att.size > np.isfinite(old_att).sum(), "Missing before/after failure: " + name
                record["disposition"] = "corrected false feasibility rejection; current matches independent base-R ATT and full IF"
                unaffected = np.asarray(new["group"]) != 3 if name.startswith("zero_treated") else None
                if unaffected is not None:
                    for field in ("att", "se", "l_se", "u_se"):
                        np.testing.assert_array_equal(np.asarray(new["results"][field])[unaffected],
                                                      np.asarray(old["results"][field])[unaffected])
                    np.testing.assert_array_equal(np.asarray(new["influence"])[:, unaffected], np.asarray(old["influence"])[:, unaffected])
                    record["unaffected_cells_exact"] = int(unaffected.sum())
            else:
                assert fingerprint(old) == fingerprint(new), "Ordinary analytic oracle before/after differs"
                record["disposition"] = "unchanged and agrees with closed-form weighted means and influence function"
        elif case.get("safe_probe"):
            old_entries, new_entries = old["backend_entries"], new["backend_entries"]
            assert any(e["lost"] for e in old_entries), "Old process did not enter unsupported backend"
            assert new_entries and not any(e["lost"] for e in new_entries)
            affected = np.asarray(new["group"]) == 2
            for field in ("att", "se", "l_se", "u_se"):
                np.testing.assert_array_equal(np.asarray(new["results"][field])[~affected],
                                              np.asarray(old["results"][field])[~affected])
            assert np.isfinite(new_att[~affected]).all() and np.isnan(new_att[affected]).all()
            if case.get("do_inf"):
                np.testing.assert_array_equal(np.asarray(new["influence"])[:, ~affected], np.asarray(old["influence"])[:, ~affected])
            record.update(old_unsupported_backend_entries=sum(e["lost"] for e in old_entries),
                          new_unsupported_backend_entries=sum(e["lost"] for e in new_entries),
                          unchanged_healthy_cells=int((~affected).sum()),
                          disposition="normalized-support loss rejected before backend; healthy estimates/IF exact")
        elif name.startswith("tiny_"):
            expected_bad = (np.asarray(new["group"]) == 2) & (case["data"] == "tiny_3.csv") & (case["method"] == "dr")
            np.testing.assert_array_equal(np.isnan(new_att), expected_bad)
            for field in ("att", "se", "l_se", "u_se"):
                np.testing.assert_array_equal(new["results"][field], old["results"][field])
            np.testing.assert_array_equal(new["influence"], old["influence"])
            record["disposition"] = "existing numerical rejection preserved; explicit physical preflight added" if expected_bad.any() else "square/full-rank or non-treated-outcome method preserved exactly"
        elif name.startswith("cache_support"):
            counterpart = name[:-4] + "false" if case["fast"] else name[:-5] + "true"
            other = after["targeted"][counterpart]
            np.testing.assert_array_equal(new_att, other["results"]["att"])
            np.testing.assert_array_equal(new["influence"], other["influence"])
            record["disposition"] = "weight-aware overlap correct and fast/standard identical"
            if case["fast"]:
                mutant = after["targeted"][name + "_weight_key_mutant"]
                mask_changed = np.isnan(np.asarray(mutant["results"]["att"])) != np.isnan(new_att)
                assert mask_changed.any(), "Weight cache-key mutant escaped test"
                record["seeded_weight_cache_defect_detected_cells"] = int(mask_changed.sum())
        checks.append(record)
    summary = dict(audit_date="2026-10-26", source_pins=PINS,
        sources={"baseline": before["source"], "candidate": after["source"]},
        harness_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in (Path(__file__), Path(__file__).with_suffix(".R"))},
        oracle_session_info=(args.output / "oracle-sessionInfo.txt").read_text(),
        oracle_policy=dict(atol=ATOL, rtol=RTOL, support_and_nonfinite_masks="exact"),
        harness_corrections=["The initial seeded-cache fixture had time-varying X, so it did not reuse a cache entry and could not detect deletion of the weight key. The corrected fixture fixes X by unit while varying weight support; production code was unchanged."],
        oracle_warning_messages=sorted({message for ref in references.values() for message in ref.get("warnings", [])}),
        independent_oracle="Base R only: explicit normalized four-cell means, weighted linear normal equations, logistic score/Hessian and estimator influence equations; no did/DRDID code called.",
        healthy_comparisons=len(healthy), healthy_full_output_exact=len(healthy)-len(mismatches),
        healthy_discrepancies=mismatches,
        healthy_combined_sha256=fingerprint(before["healthy"]),
        healthy_scope="All methods, both controls/bases/modes; panel/unbalanced anchors default/first/base/varying and RC default; absent/unit/nonconstant/time-varying-zero weights; all four aggregations, complete IFs, warnings and RNG states.",
        targeted_totals=dict(
            checks=len(checks), independent_base_r_oracles=sum("oracle_att_max_abs" in z for z in checks),
            independent_max_att=max(z.get("oracle_att_max_abs", 0.) for z in checks),
            independent_max_if=max(z.get("oracle_if_max_abs", 0.) for z in checks),
            restored_group_time_results=sum(z["new_finite_att"] - z["old_finite_att"] for z in checks
                                           if z["disposition"].startswith("corrected false")),
            old_unsupported_backend_entries=sum(z.get("old_unsupported_backend_entries", 0) for z in checks),
            new_unsupported_backend_entries=sum(z.get("new_unsupported_backend_entries", 0) for z in checks),
            seeded_weight_cache_defects_detected=sum("seeded_weight_cache_defect_detected_cells" in z for z in checks)),
        targeted_checks=checks,
        input_sha256={name: hashlib.sha256((args.output / name).read_bytes()).hexdigest()
                      for name in sorted({case["data"] for case in healthy + targets})},
        limitations=["Finite targeted and grid evidence cannot prove absence of all regressions.",
                    "Physical counts are necessary, not full-rank certificates; extreme floating-point overflow and solver-threshold boundaries remain numerical limits.",
                    "Old tiny-design cases already rejected by their numerical Gram check are reported as preserved rejections, not newly fixed estimates."])
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"healthy_exact": summary["healthy_full_output_exact"], "targeted": len(checks),
                      "oracle_max_if": max(z.get("oracle_if_max_abs", 0.) for z in checks)}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--role", choices=["baseline", "candidate"])
    parser.add_argument("--compare-only", action="store_true")
    parser.add_argument("--targets-only", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.worker:
        worker(args)
        return
    if args.compare_only:
        cases = json.loads((args.output / "cases.json").read_text())
        compare(args, cases["healthy"], cases["targeted"])
        return
    for role, repo in (("baseline", args.baseline), ("candidate", args.candidate)):
        assert repo is not None and source_info(repo)["commit"] == PINS[role], role + " source pin"
    healthy, targeted = prepare(args.output)
    for role, repo in (("baseline", args.baseline), ("candidate", args.candidate)):
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", "--repo", str(repo),
                        "--role", role, "--output", str(args.output)], check=True)
    subprocess.run(["Rscript", str(Path(__file__).with_suffix(".R")), str(args.output)], check=True)
    compare(args, healthy, targeted)


if __name__ == "__main__":
    main()
