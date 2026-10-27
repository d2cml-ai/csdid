"""Fresh differential audit against the pinned R sources and frozen inputs.

Run with the test environment's Python from the repository root::

    python csdid/test_csdid/r_ref/current/stress_compare.py \
      --did /clean/did --drdid /clean/DRDID --output /tmp/csdid-stress

This does not overwrite the permanent 25-case references. The output retains
the expanded case manifest, exact R harness, source provenance, R outputs,
per-case Python comparisons and any failures. Comparison tolerances and
nonfinite/support checks are the permanent numerical gate's unchanged ones.
"""
import argparse
import copy
import gzip
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


HERE = Path(__file__).resolve().parent
TESTS = HERE.parents[1]
REPO = HERE.parents[3]


def prepare(output, balance_grid=False):
    output.mkdir(parents=True, exist_ok=True)
    for name in ("panel.csv", "rc.csv", "unbalanced.csv",
                 "multipliers_unit.csv", "multipliers_cluster.csv"):
        shutil.copyfile(HERE / name, output / name)
    original = json.loads((HERE / "cases.json").read_text())
    cases = []
    for design in ("panel", "rc", "unbalanced"):
        for method in ("dr", "reg", "ipw"):
            prototype = next(z for z in original if z["name"] ==
                             f"{design}_{method}_never_varying")
            for base in ("varying", "universal"):
                for control in ("nevertreated", "notyettreated"):
                    for anticipation in (0, 1, 2):
                        z = copy.deepcopy(prototype)
                        z["name"] = f"{design}_{method}_{base}_{control}_a{anticipation}"
                        z["constructor"].update(control_group=control,
                                                anticipation=anticipation)
                        z["fit"]["base_period"] = base
                        z["expected_r_kwargs"].update(control_group=control,
                            anticipation=anticipation, base_period=base)
                        cases.append(z)
                for scale in (1e-10, 1e10):
                    z = copy.deepcopy(prototype)
                    z["name"] = f"{design}_{method}_{base}_scale{scale:g}"
                    z["outcome_scale"] = scale
                    z["constructor"]["anticipation"] = 1
                    z["fit"]["base_period"] = base
                    z["expected_r_kwargs"].update(anticipation=1, base_period=base)
                    z["data"] = f"{design}_scale{scale:g}.csv"
                    d = pd.read_csv(HERE / prototype["data"],
                                    converters={k: float for k in ("X", "Y", "weight")})
                    d["Y"] *= scale
                    d.to_csv(output / z["data"], index=False, float_format="%.17g")
                    cases.append(z)
    assert len(cases) == 144
    if balance_grid:
        cases = []
        # These clocks intentionally lack event-time integers between observed
        # periods. Cohort support is determined before display limits apply.
        for periods, cohorts, balance in (((1, 2, 4, 5), (0, 4, 5), 1),
                                          ((1, 3, 5), (0, 3, 5), 2)):
            d = pd.read_csv(HERE / "panel.csv",
                            converters={k: float for k in ("X", "Y", "weight")})
            d = d[d["period"].isin(periods) & d["G"].isin(cohorts)]
            filename = f"panel_irregular_balance{balance}.csv"
            d.to_csv(output / filename, index=False, float_format="%.17g")
            for method in ("dr", "reg", "ipw"):
                prototype = next(z for z in original if z["name"] ==
                                 f"panel_{method}_never_varying")
                for base in ("varying", "universal"):
                    for display_max in (0, 99):
                        z = copy.deepcopy(prototype)
                        z["name"] = f"irregular_b{balance}_{method}_{base}_max{display_max}"
                        z["data"] = filename
                        z["fit"]["base_period"] = base
                        z["expected_r_kwargs"]["base_period"] = base
                        options = dict(balance_e=balance, min_e=0, max_e=display_max)
                        z["aggregation"].update(options)
                        z["expected_r_aggregation_kwargs"].update(options)
                        cases.append(z)
        assert len(cases) == 24
    (output / "stress_cases.json").write_text(json.dumps(cases, indent=2) + "\n")
    source = (HERE / "generate.R").read_text()
    sentinel = 'jsonlite::write_json(cases, file.path(output, "cases.json"), pretty = TRUE'
    assert source.count(sentinel) == 1
    source = source.replace(sentinel,
        'cases <- jsonlite::fromJSON(file.path(output, "stress_cases.json"), simplifyVector = FALSE)\n' + sentinel)
    old = 'return(list(list(field = path, same_support = same_support, max_abs = max_abs)))'
    assert source.count(old) == 1
    source = source.replace(old, '''wald <- path %in% c("fit.W", "fit.Wpval")
    atol <- if (wald) 1e-5 else 1e-7
    rtol <- if (wald) 1e-5 else 1e-6
    within <- same_support && all(abs(av - bv) <= atol + rtol * abs(bv))
    return(list(list(field = path, same_support = same_support, max_abs = max_abs,
                     within_tolerance = within)))''')
    (output / "generate_stress.R").write_text(source)
    return cases


def normalized_output(actual, reference, scale):
    """Supplemental unit diagnosis; never changes the primary comparison."""
    actual, reference = copy.deepcopy(actual), copy.deepcopy(reference)
    for field in ("att", "se", "l_se", "u_se", "inffunc"):
        actual["fit"][field] = np.asarray(actual["fit"][field]) / scale
        payload = reference["fit"][field]
        payload["values"] = [value / scale for value in payload["values"]]
    for kind, aggregation in actual["aggs"].items():
        ref = reference["aggs"][kind]
        for field in ("overall_att", "overall_se", "att_egt", "se_egt"):
            if aggregation[field] is not None:
                aggregation[field] = np.asarray(aggregation[field]) / scale
                payload = ref[field]
                payload["values"] = [value / scale for value in payload["values"]]
        for field in aggregation["inf_function"]:
            aggregation["inf_function"][field] = np.asarray(aggregation["inf_function"][field]) / scale
            payload = ref["inf_function"][field]
            payload["values"] = [value / scale for value in payload["values"]]
    return actual, reference


def compare(output):
    os.environ["CSDID_R_REFERENCE_DIR"] = str(output)
    sys.path[:0] = [str(TESTS), str(REPO)]
    spec = importlib.util.spec_from_file_location("current_parity", TESTS / "test_r_current_parity.py")
    parity = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parity)
    with gzip.open(output / "references.json.gz", "rt") as stream:
        refs = json.load(stream)
    cases = json.loads((output / "cases.json").read_text())
    assert set(refs) == {z["name"] for z in cases}
    result = {"case_count": len(cases), "comparisons": [], "failures": [],
              "tolerance_policy": json.loads((output / "manifest.json").read_text())["tolerance_policy"]}
    checks = json.loads((output / "r_mode_checks.json").read_text())
    result["r_mode_failures"] = [{"case": name, **field} for name, fields in checks.items()
        for field in fields if not field["same_support"] or not field["within_tolerance"]]
    for i, z in enumerate(cases):
        for fast in (False, True):
            started = time.monotonic()
            row = {"case": z["name"], "faster_mode": fast}
            actual = None
            try:
                actual = parity._run_case(z, fast)
                parity._assert_case_output(actual, refs[z["name"]])
                row["passed"] = True
            except Exception as exc:
                row.update(passed=False, error=str(exc), exception=type(exc).__name__)
            scale = z.get("outcome_scale")
            if scale is None and "_scale" in z["name"]:
                scale = float(z["name"].rsplit("_scale", 1)[1])
            if scale is not None and actual is not None:
                # Dividing outcome-unit fields back by the declared scale tests
                # whether raw-scale failures are amplified numerical residuals.
                # Record this separately: a normalized pass cannot erase a raw
                # comparison failure or alter the script's nonzero exit status.
                try:
                    aa, rr = normalized_output(actual, refs[z["name"]], scale)
                    parity._assert_case_output(aa, rr)
                    row["scale_normalized_passed"] = True
                except Exception as exc:
                    row.update(scale_normalized_passed=False, normalized_error=str(exc))
            if not row["passed"]:
                result["failures"].append(row.copy())
            row["seconds"] = time.monotonic() - started
            result["comparisons"].append(row)
        if i % 12 == 0:
            print(f"Compared {i + 1}/{len(cases)} cases; failures={len(result['failures'])}", flush=True)
    result["passed_comparisons"] = sum(r["passed"] for r in result["comparisons"])
    (output / "comparison_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "comparisons"}, indent=2))
    return not result["failures"] and not result["r_mode_failures"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--did", type=Path)
    parser.add_argument("--drdid", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--compare-only", action="store_true")
    parser.add_argument("--balance-grid", action="store_true",
                        help="run the separate 24-case irregular-clock balance_e grid")
    args = parser.parse_args()
    output = args.output.resolve()
    if not args.compare_only:
        if args.did is None or args.drdid is None:
            parser.error("--did and --drdid are required unless --compare-only is set")
        prepare(output, balance_grid=args.balance_grid)
        env = dict(os.environ, **{k: "1" for k in ("OPENBLAS_NUM_THREADS",
                   "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")})
        with (output / "r_generation.log").open("w") as log:
            subprocess.run(["Rscript", str(output / "generate_stress.R"),
                "--did", str(args.did.resolve()), "--drdid", str(args.drdid.resolve()),
                "--output", str(output)], env=env, stdout=log,
                stderr=subprocess.STDOUT, check=True)
    return 0 if compare(output) else 1


if __name__ == "__main__":
    raise SystemExit(main())
