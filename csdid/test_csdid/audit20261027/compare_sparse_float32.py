#!/usr/bin/env python3
"""Compare sparse smallest-positive float32 weights with their R double values."""
import argparse
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pytest

import compare_weight_scale as audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csdid-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--did", type=Path, required=True)
    parser.add_argument("--drdid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root, output = args.csdid_root.resolve(), args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("--output must be new or empty")
    if output == root or root in output.parents:
        parser.error("--output must be outside the repository")
    output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    before = audit.frozen_hashes(root)
    sources = dict(csdid=audit.source_record(root),
                   did=audit.source_record(args.did, "did", audit.DID_COMMIT),
                   DRDID=audit.source_record(args.drdid, "DRDID", audit.DRDID_COMMIT))
    sys.path.insert(0, str(root))
    sys.path.insert(1, str(root / "csdid" / "test_csdid"))
    from test_sampling_weight_scale import _two_period_data
    tiny = np.nextafter(np.float32(0), np.float32(1))
    cases, frames, inputs = [], {}, {}
    for panel in (False, True):
        dataset = "two_panel" if panel else "two_rc"
        frame = _two_period_data(panel)
        sparse = (np.arange(len(frame)) % 5 == 0).astype(np.float32)
        tiny_weights = sparse * tiny
        inputs[dataset] = dict(rows=len(frame), positive_rows=int(np.count_nonzero(sparse)),
            positive_tiny_value=float(tiny), tiny_mean_float32=float(tiny_weights.mean()),
            tiny_mean_float64=float(tiny_weights.mean(dtype=np.float64)))
        for scale in (np.float32(1), tiny):
            data = frame.copy(deep=True)
            data["wt"] = sparse * scale
            assert data["wt"].dtype == np.float32
            filename = "{}_sparse_{}.csv".format(dataset, format(float(scale), ".0e"))
            disk = data.copy(deep=True)
            # R has double numeric vectors. Preserve exactly the already
            # represented float32 weights, rather than rounding their mean.
            disk["wt"] = disk["wt"].astype(np.float64)
            disk.to_csv(output / filename, index=False, float_format="%.17g")
            for method in ("dr", "ipw", "reg"):
                case = audit.case_spec(dataset, method, float(scale), filename)
                cases.append(case)
                frames[case["name"]] = data
    (output / "cases.json").write_text(json.dumps(cases, indent=2) + "\n")
    r_script = Path(__file__).with_name("weight_scale_references.R")
    with (output / "r-generation.log").open("w") as stream:
        subprocess.run(["Rscript", str(r_script), "--did", str(args.did.resolve()),
                        "--drdid", str(args.drdid.resolve()), "--output", str(output)],
                       check=True, stdout=stream, stderr=subprocess.STDOUT)
    with gzip.open(output / "r_references.json.gz", "rt") as stream:
        references = json.load(stream)
    os.environ["CSDID_R_REFERENCE_DIR"] = str(output)
    os.environ.pop("CSDID_R_REFERENCE_SOURCES", None)
    spec = importlib.util.spec_from_file_location("float32_strict_gate",
        root / "csdid" / "test_csdid" / "test_r_current_parity.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    categories = {name: dict(passed=0, failed=0) for name in (
        "cross_language", "python_modes", "r_modes", "python_weight_units", "r_weight_units")}
    failures, maximum, actuals = [], {}, {}
    def check(category, name, actual, reference):
        try:
            gate._assert_case_output(actual, reference)
        except Exception as error:
            categories[category]["failed"] += 1
            failures.append(dict(category=category, case=name, error=str(error)))
        else:
            categories[category]["passed"] += 1
    for case in cases:
        original_read_csv = gate.pd.read_csv
        def read_preserving_float32(path, *positional, **keywords):
            if Path(path).resolve() == (output / case["data"]).resolve():
                return frames[case["name"]].copy(deep=True)
            return original_read_csv(path, *positional, **keywords)
        for fast in (False, True):
            label = case["name"] + ("__fast" if fast else "__standard")
            try:
                # Reuse the exact permanent output collector while ensuring
                # Python really receives float32, which CSV parsing would lose.
                with pytest.MonkeyPatch.context() as patch:
                    patch.setattr(gate.pd, "read_csv", read_preserving_float32)
                    actuals[label] = gate._run_case(case, fast)
            except Exception as error:
                failures.append(dict(category="python_execution", case=label,
                                     error="{}: {}".format(type(error).__name__, error)))
                continue
            for r_fast in (False, True):
                r_label = case["name"] + ("__fast" if r_fast else "__standard")
                check("cross_language", label + " vs " + r_label, actuals[label], references[r_label])
                audit.compare_diagnostics(actuals[label], references[r_label], gate, maximum)
        standard, fast = case["name"] + "__standard", case["name"] + "__fast"
        if standard in actuals and fast in actuals:
            check("python_modes", case["name"], actuals[fast], audit.encoded_python(actuals[standard]))
        check("r_modes", case["name"], audit.decoded_reference(references[fast], gate), references[standard])
        if case["scale"] != 1.:
            baseline = "{}_{}_1e+00".format(case["dataset"], case["fit"]["est_method"])
            for mode in ("__fast", "__standard"):
                if case["name"] + mode in actuals and baseline + mode in actuals:
                    check("python_weight_units", case["name"] + mode,
                          actuals[case["name"] + mode], audit.encoded_python(actuals[baseline + mode]))
                check("r_weight_units", case["name"] + mode,
                      audit.decoded_reference(references[case["name"] + mode], gate), references[baseline + mode])
    with gzip.open(output / "python_payloads.json.gz", "wt") as stream:
        json.dump({name: audit.encoded_python(value) for name, value in actuals.items()}, stream)
    final_sources = dict(csdid=audit.source_record(root),
                         did=audit.source_record(args.did, "did", audit.DID_COMMIT),
                         DRDID=audit.source_record(args.drdid, "DRDID", audit.DRDID_COMMIT))
    unchanged = sources == final_sources and before == audit.frozen_hashes(root)
    record = dict(schema_version=1, status="passed" if not failures and unchanged else "failed",
        elapsed_seconds=round(time.time() - started, 3), source_identity=sources,
        sources_and_frozen_fixtures_unchanged=unchanged, input_description=inputs,
        cases=len(cases), r_fits=2 * len(cases), python_fits=len(actuals), checks=categories,
        maximum_cross_language_absolute_errors=maximum,
        tolerance_policy=dict(public_and_full_influence=dict(atol=1e-7, rtol=1e-6),
                              wald=dict(atol=1e-5, rtol=1e-5), support_shape_nonfinite_masks="exact"),
        scripts_sha256={path.name: audit.sha256(path) for path in (
            Path(__file__), Path(audit.__file__), r_script)},
        dgp_test_sha256=audit.sha256(root / "csdid" / "test_csdid" / "test_sampling_weight_scale.py"),
        payload_sha256={name: audit.sha256(output / name) for name in (
            "r_references.json.gz", "python_payloads.json.gz", "cases.json", "sessionInfo.txt")},
        shared_csv_sha256={name: audit.sha256(output / name) for name in sorted({z["data"] for z in cases})},
        failures=failures,
        scope="Sparse float32 in-memory Python weights; R receives their exact double representations. Analytical inference.")
    (output / "result.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({key: record[key] for key in ("status", "r_fits", "python_fits", "checks",
                      "input_description", "maximum_cross_language_absolute_errors")}, indent=2))
    return 0 if record["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
