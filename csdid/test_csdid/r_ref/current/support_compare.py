"""Compare weighted guard regressions against clean current R sources.

Example (run from the repository root with its test environment)::

    python csdid/test_csdid/r_ref/current/support_compare.py \\
      --did /path/to/did --drdid /path/to/DRDID --output /tmp/support-compare

The deterministic inputs come from the portable support regressions. The R
export preserves original unit identities; masks are exact and public outputs
and every IF component use the existing atol=1e-7, rtol=1e-6 policy.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(variable, "1")

import numpy as np
import pandas as pd


def make_inputs(fixtures, output):
    cases = []
    for route in ("panel", "varying", "rc", "unbalanced"):
        for method in ("reg", "dr", "ipw"):
            for outlier in (False, True):
                data = fixtures._panel()
                if route == "unbalanced":
                    data = data[(data.id * 7 + data.period) % 17 != 0]
                if outlier:
                    extra = pd.DataFrame({
                        "id": 100, "G": 0 if method == "reg" else 3,
                        "X": 1e12 if method == "reg" else 1000.,
                        "weight": 0., "a": 0., "period": np.arange(1, 5), "Y": 0.,
                    })
                    data = pd.concat([data, extra], ignore_index=True)
                if route == "rc":
                    data["id"] = np.arange(len(data))
                name = f"{route}_{method}_{int(outlier)}"
                data.to_csv(output / f"{name}.csv", index=False)
                cases.append(dict(name=name, route=route, method=method))
    (output / "cases.json").write_text(json.dumps(cases))
    return cases


def compare_case(fixtures, data, case, mode, reference):
    fit = fixtures._fit(data, mode, case["method"], case["route"])
    record = {"case": f"{case['name']}_{str(mode).upper()}"}
    for field, actual in (("att", fit.MP["att"]), ("se", fit.results["se"])):
        expected = np.array([np.nan if x is None else x for x in reference[field]])
        actual = np.asarray(actual)
        np.testing.assert_array_equal(np.isnan(actual), np.isnan(expected))
        np.testing.assert_allclose(actual, expected, atol=1e-7, rtol=1e-6)
        record[field + "_max_abs"] = float(np.nanmax(np.abs(actual - expected)))

    # Match authoritative row identities, never coincidentally equal positions.
    expected = np.array(reference["if_values"]).reshape(reference["if_shape"])
    actual = np.asarray(fit.MP["inffunc"]["inffunc"])
    internal_id = fit.dp["idname"]
    ids = fit.dp["data"][internal_id].unique()
    if internal_id != "id":
        mapping = fit.dp["data"][[internal_id, "id"]].drop_duplicates()
        assert mapping[internal_id].is_unique
        ids = mapping.set_index(internal_id)["id"].reindex(ids).to_numpy()
    assert pd.Index(ids).is_unique
    order = pd.Index(ids).get_indexer(reference["if_ids"])
    assert len(order) == len(ids) and len(set(reference["if_ids"])) == len(ids)
    assert np.all(order >= 0)
    actual = actual[order]
    np.testing.assert_array_equal(np.isnan(actual), np.isnan(expected))
    np.testing.assert_allclose(actual, expected, atol=1e-7, rtol=1e-6)
    record["if_max_abs"] = float(np.nanmax(np.abs(actual - expected)))
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--did", type=Path, required=True)
    parser.add_argument("--drdid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[4]
    sys.path.insert(0, str(repo))
    fixture_path = repo / "csdid/test_csdid/test_upstream_20261026_support.py"
    spec = importlib.util.spec_from_file_location("support_tests", fixture_path)
    fixtures = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixtures)
    cases = make_inputs(fixtures, output)
    subprocess.run([
        "Rscript", str(Path(__file__).with_suffix(".R")), str(output),
        str(args.did.resolve()), str(args.drdid.resolve()),
    ], check=True)
    references = json.loads((output / "r-results.json").read_text())
    checks = []
    for case in cases:
        data = pd.read_csv(output / f"{case['name']}.csv")
        for mode in (False, True):
            key = f"{case['name']}_{str(mode).upper()}"
            checks.append(compare_case(fixtures, data, case, mode, references[key]))
    summary = {
        "audit_date": "2026-10-26",
        "sources": json.loads((output / "r-sources.json").read_text()),
        "fit_comparisons": len(checks), "atol": 1e-7, "rtol": 1e-6,
        "nonfinite_masks": "exact",
        "input_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(output.glob("*.csv"))
        },
        "checks": checks,
    }
    (output / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({
        "fits": len(checks),
        "max_att": max(x["att_max_abs"] for x in checks),
        "max_se": max(x["se_max_abs"] for x in checks),
        "max_if": max(x["if_max_abs"] for x in checks),
    }))


if __name__ == "__main__":
    main()
