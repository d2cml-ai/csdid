#!/usr/bin/env python3
"""Replay full-output R/Python weight-scale comparisons without changing fixtures."""
import argparse
import copy
import gzip
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd


SCALES = (1., 1e-200, 1e-160, 1e155, 1e307)
DATASETS = ("two_panel", "two_rc", "multi_panel", "multi_rc")
DID_COMMIT = "75168454d6aa7faa2f9167b05ae7ad4f38f9e27b"
DRDID_COMMIT = "eaac3b49cd225d0d7706c633858555a5f204e21e"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(path, *args):
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


def source_record(path, package=None, expected=None):
    path = Path(path).resolve()
    head = git(path, "rev-parse", "HEAD")
    if expected is not None and head != expected:
        raise ValueError("Unexpected source commit for {}: {}".format(package, head))
    if package:
        changed = git(path, "status", "--porcelain", "--untracked-files=no")
        if changed:
            raise ValueError("R sources must be clean: " + str(path))
        names = git(path, "ls-files", "DESCRIPTION", "NAMESPACE", "R", "src").splitlines()
        files = {name: sha256(path / name) for name in names if (path / name).is_file()}
    else:
        files = {str(file.relative_to(path)): sha256(file)
                 for file in sorted((path / "csdid").rglob("*.py"))
                 if "test_csdid" not in file.parts}
    return {"git_head": head, "source_files_sha256": files,
            "source_tree_sha256": hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()}


def frozen_hashes(root):
    folder = root / "csdid" / "test_csdid" / "r_ref" / "current"
    return {str(file.relative_to(folder)): sha256(file)
            for file in sorted(folder.rglob("*")) if file.is_file()}


def generate_two_period_inputs(output):
    rng = np.random.default_rng(8201)
    n = 240
    x, z = rng.normal(size=(2, n))
    probability = 1 / (1 + np.exp(.3 - .5 * x + .2 * z - .4 * np.sin(2 * x)))
    treated = (rng.random(n) < probability).astype(float)
    y0 = .8 + .3 * x - .2 * z + rng.normal(size=n)
    y1 = y0 + .6 + .4 * x + 1.1 * treated + rng.normal(size=n)
    weights = np.exp(.2 * x - .15 * z)
    weights[::11] = 0
    frame = pd.DataFrame(dict(id=np.tile(np.arange(n), 2), period=np.repeat([1, 2], n),
                              G=np.tile(2 * treated, 2), Y=np.r_[y0, y1],
                              X=np.tile(x, 2), wt=np.tile(weights, 2)))
    frame.to_csv(output / "two_panel.csv", index=False)
    frame["id"] = np.arange(len(frame))
    frame.to_csv(output / "two_rc.csv", index=False)


def case_spec(dataset, method, scale, filename):
    panel = dataset.endswith("panel")
    ctor = dict(yname="Y", tname="period", idname="id", gname="G", xformla="~X",
                panel=panel, weights_name="wt" if dataset.startswith("two") else "weight",
                cband=False)
    fit = dict(est_method=method, base_period="varying", bstrap=False)
    aggregation = dict(bstrap=False, cband=False, na_rm=False)
    r_ctor = dict(ctor)
    r_ctor["weightsname"] = r_ctor.pop("weights_name")
    return dict(name="{}_{}_{}".format(dataset, method, format(scale, ".0e")),
                dataset=dataset, scale=scale, data=filename, constructor=ctor,
                fit=fit, aggregation=aggregation,
                expected_r_kwargs=dict(r_ctor, **fit),
                expected_r_aggregation_kwargs=dict(bstrap=False, cband=False, **{"na.rm": False}))


def numeric_payload(value):
    if value is None:
        return None
    array = np.atleast_1d(np.asarray(value, dtype=float))
    flat = array.ravel().copy()
    masks = {key: np.flatnonzero(function(flat)).tolist() for key, function in (
        ("nan", np.isnan), ("posinf", np.isposinf), ("neginf", np.isneginf))}
    flat[~np.isfinite(flat)] = 0
    return dict(shape=list(array.shape), values=flat.tolist(), **masks)


def encoded_python(actual):
    result = dict(unit_ids=np.asarray(actual["unit_ids"]).tolist(), fit={}, aggs={})
    for key, value in actual["fit"].items():
        result["fit"][key] = bool(value) if key == "cband" else numeric_payload(value)
    for kind, aggregation in actual["aggs"].items():
        result["aggs"][kind] = {key: numeric_payload(value) for key, value in aggregation.items()
                                if key != "inf_function"}
        result["aggs"][kind]["inf_function"] = {
            key: numeric_payload(np.asarray(value).reshape(len(actual["unit_ids"]), -1))
            for key, value in aggregation["inf_function"].items()}
    return result


def decoded_reference(reference, module):
    return dict(unit_ids=reference["unit_ids"],
                fit={key: value if key == "cband" else module.decode_numeric(value)
                     for key, value in reference["fit"].items()},
                aggs={kind: dict({key: module.decode_numeric(value) for key, value in agg.items()
                                 if key != "inf_function"}, inf_function={
                                     key: module.decode_numeric(value)
                                     for key, value in agg["inf_function"].items()})
                      for kind, agg in reference["aggs"].items()})


def compare_diagnostics(actual, reference, module, maximum):
    cells = module.support_order(zip(actual["fit"]["group"], actual["fit"]["t"]),
                                zip(module.decode_numeric(reference["fit"]["group"]),
                                    module.decode_numeric(reference["fit"]["t"])), "cells")
    units = module.support_order(actual["unit_ids"], reference["unit_ids"], "units")
    pairs = [("att", np.asarray(actual["fit"]["att"])[cells],
              module.decode_numeric(reference["fit"]["att"])),
             ("cell_if", np.asarray(actual["fit"]["inffunc"])[np.ix_(units, cells)],
              module.decode_numeric(reference["fit"]["inffunc"])),
             ("cell_se", np.asarray(actual["fit"]["se"])[cells],
              module.decode_numeric(reference["fit"]["se"]))]
    for kind in module.KINDS:
        agg, ref = actual["aggs"][kind], reference["aggs"][kind]
        pairs.append(("aggregate_se", agg["overall_se"], module.decode_numeric(ref["overall_se"])))
        for field, value in agg["inf_function"].items():
            pairs.append(("aggregate_if", np.asarray(value).reshape(len(units), -1)[units],
                          module.decode_numeric(ref["inf_function"][field])))
    for field, left, right in pairs:
        left, right = np.asarray(left).ravel(), np.asarray(right).ravel()
        finite = np.isfinite(left) & np.isfinite(right)
        if finite.any():
            maximum[field] = max(maximum.get(field, 0.), float(np.max(np.abs(left[finite] - right[finite]))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csdid-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--did", type=Path, required=True)
    parser.add_argument("--drdid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path,
                        help="Use existing panel.csv/rc.csv and multi_panel.csv/multi_rc.csv authorities")
    options = parser.parse_args()
    root, output = options.csdid_root.resolve(), options.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("--output must be new or empty")
    if output == root or root in output.parents:
        parser.error("--output must be outside the repository")
    output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    before = frozen_hashes(root)
    sources = dict(csdid=source_record(root), did=source_record(options.did, "did", DID_COMMIT),
                   DRDID=source_record(options.drdid, "DRDID", DRDID_COMMIT))
    script = Path(__file__).with_name("weight_scale_references.R")
    r_command = ["Rscript", str(script), "--output", str(output),
                 "--did", str(options.did.resolve()), "--drdid", str(options.drdid.resolve())]
    if options.input_dir:
        for dataset in DATASETS:
            filename = dataset.replace("two_", "") + ".csv"
            shutil.copyfile(options.input_dir / filename, output / (dataset + ".csv"))
    else:
        generate_two_period_inputs(output)
        subprocess.run(r_command + ["--generate-inputs"], check=True)
    cases = []
    for dataset in DATASETS:
        original = pd.read_csv(output / (dataset + ".csv"),
                               converters={key: float for key in ("X", "Y", "weight", "wt")})
        weight = "wt" if dataset.startswith("two") else "weight"
        for scale in SCALES:
            data = original.copy(deep=True)
            data[weight] *= scale
            filename = "{}_{}.csv".format(dataset, format(scale, ".0e"))
            data.to_csv(output / filename, index=False)
            for method in ("dr", "ipw", "reg"):
                cases.append(case_spec(dataset, method, scale, filename))
    (output / "cases.json").write_text(json.dumps(cases, indent=2) + "\n")
    with (output / "r-generation.log").open("w") as stream:
        subprocess.run(r_command, check=True, stdout=stream, stderr=subprocess.STDOUT)
    with gzip.open(output / "r_references.json.gz", "rt") as stream:
        references = json.load(stream)
    os.environ["CSDID_R_REFERENCE_DIR"] = str(output)
    os.environ.pop("CSDID_R_REFERENCE_SOURCES", None)
    sys.path.insert(0, str(root))
    sys.path.insert(1, str(root / "csdid" / "test_csdid"))
    spec = importlib.util.spec_from_file_location("weight_scale_strict_gate",
        root / "csdid" / "test_csdid" / "test_r_current_parity.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    categories = {name: dict(passed=0, failed=0) for name in (
        "cross_language", "python_modes", "r_modes", "python_weight_units", "r_weight_units")}
    failures, maximum = [], {}
    actuals = {}
    def check(category, name, actual, reference):
        try:
            module._assert_case_output(actual, reference)
        except Exception as error:
            categories[category]["failed"] += 1
            failures.append(dict(category=category, case=name,
                                 error="{}: {}".format(type(error).__name__, error)))
        else:
            categories[category]["passed"] += 1
    for case in cases:
        for fast in (False, True):
            label = case["name"] + ("__fast" if fast else "__standard")
            try:
                actuals[label] = module._run_case(case, fast)
            except Exception as error:
                failures.append(dict(category="python_execution", case=label,
                                     error="{}: {}".format(type(error).__name__, error)))
                continue
            for r_fast in (False, True):
                r_label = case["name"] + ("__fast" if r_fast else "__standard")
                check("cross_language", label + " vs " + r_label,
                      actuals[label], references[r_label])
                compare_diagnostics(actuals[label], references[r_label], module, maximum)
        standard, fast = case["name"] + "__standard", case["name"] + "__fast"
        if standard in actuals and fast in actuals:
            check("python_modes", case["name"], actuals[fast], encoded_python(actuals[standard]))
        check("r_modes", case["name"], decoded_reference(references[fast], module), references[standard])
        baseline = "{}_{}_{}".format(case["dataset"], case["fit"]["est_method"], format(1., ".0e"))
        if case["scale"] != 1.:
            for mode in ("__fast", "__standard"):
                if case["name"] + mode in actuals and baseline + mode in actuals:
                    check("python_weight_units", case["name"] + mode,
                          actuals[case["name"] + mode], encoded_python(actuals[baseline + mode]))
                check("r_weight_units", case["name"] + mode,
                      decoded_reference(references[case["name"] + mode], module), references[baseline + mode])
    with gzip.open(output / "python_payloads.json.gz", "wt") as stream:
        json.dump({name: encoded_python(value) for name, value in actuals.items()}, stream)
    after = frozen_hashes(root)
    final_sources = dict(csdid=source_record(root), did=source_record(options.did, "did", DID_COMMIT),
                         DRDID=source_record(options.drdid, "DRDID", DRDID_COMMIT))
    unchanged = sources == final_sources and before == after
    record = dict(schema_version=1, status="passed" if not failures and unchanged else "failed",
        elapsed_seconds=round(time.time() - started, 3), source_identity=sources,
        sources_and_frozen_fixtures_unchanged=unchanged,
        frozen_fixture_hashes_sha256=hashlib.sha256(json.dumps(before, sort_keys=True).encode()).hexdigest(),
        input_sha256={dataset + ".csv": sha256(output / (dataset + ".csv")) for dataset in DATASETS},
        scripts_sha256={Path(__file__).name: sha256(__file__), script.name: sha256(script)},
        comparator_sha256={name: sha256(root / "csdid" / "test_csdid" / name)
                           for name in ("test_r_current_parity.py", "r_reference_checks.py")},
        cases=len(cases), r_fits=2 * len(cases), python_fits=len(actuals), checks=categories,
        tolerance_policy=dict(public_and_full_influence=dict(atol=1e-7, rtol=1e-6),
                              wald=dict(atol=1e-5, rtol=1e-5), support_shape_nonfinite_masks="exact"),
        maximum_cross_language_absolute_errors=maximum,
        packages={name: importlib.metadata.version(name) for name in (
            "csdid", "numpy", "pandas", "scipy", "statsmodels", "drdid", "pytest")},
        payload_sha256={name: sha256(output / name) for name in (
            "r_references.json.gz", "python_payloads.json.gz", "cases.json", "sessionInfo.txt")},
        failures=failures,
        scope="Analytical inference; identical data and estimators/settings. No random bootstrap comparison.")
    (output / "result.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({key: value for key, value in record.items()
                      if key in ("status", "cases", "r_fits", "python_fits", "checks",
                                 "maximum_cross_language_absolute_errors", "elapsed_seconds")}, indent=2))
    return 0 if record["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
