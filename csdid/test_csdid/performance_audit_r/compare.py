"""Compare R exports with one benchmark_v043 Python worker snapshot.

Public analytic outputs are a required gate (atol=1e-7, rtol=1e-6; Wald 1e-5).
Full IF discrepancies are reported separately at the same componentwise
tolerance because the packages retain different nuisance solvers. Bootstrap
SEs/bands are not a parity gate: this benchmark uses independent native RNGs.
"""
import argparse
import gzip
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_v043 import load_snapshot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-snapshot", type=Path, required=True,
                        help="Worker snapshot stem, without .json/.npz")
    parser.add_argument("--r-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-influence", action="store_true")
    args = parser.parse_args()
    py = load_snapshot(args.python_snapshot)
    with gzip.open(args.r_output / "snapshot.json.gz", "rt") as stream:
        rr = json.load(stream)
    issues, influence_issues, maxima, counts = [], [], {}, {}

    def array(value):
        if value is None:
            return np.array([], dtype=float)
        a = np.asarray(value)
        if a.dtype.kind in "USO":
            def decode(v):
                if v in ("NA", "NaN", None):
                    return float("nan")
                return float(v)
            a = np.vectorize(decode, otypes=[float])(a)
        return np.asarray(a, dtype=float)

    def check(case, field, left, right, influence=False, matrix=False):
        bucket = influence_issues if influence else issues
        a, b = array(left), array(right)
        if matrix:
            a = a.reshape(-1, 1) if a.ndim == 1 else a
            b = b.reshape(-1, 1) if b.ndim == 1 else b
        else:
            a, b = a.ravel(), b.ravel()
        counts[field] = counts.get(field, 0) + 1
        if a.shape != b.shape:
            bucket.append(dict(case=case, field=field, kind="shape",
                               python=list(a.shape), R=list(b.shape)))
            return
        for label, mask in (("finite", np.isfinite), ("nan", np.isnan),
                            ("posinf", np.isposinf), ("neginf", np.isneginf)):
            if not np.array_equal(mask(a), mask(b)):
                bucket.append(dict(case=case, field=field, kind=label + "_mask"))
                return
        finite = np.isfinite(a)
        if not finite.any():
            return
        difference = float(np.max(np.abs(a[finite] - b[finite])))
        maxima[field] = max(maxima.get(field, 0.), difference)
        atol, rtol = (1e-5, 1e-5) if field in ("W", "Wpval") else (1e-7, 1e-6)
        if not np.allclose(a[finite], b[finite], atol=atol, rtol=rtol):
            bucket.append(dict(case=case, field=field, kind="value", max_abs=difference))

    if py.keys() != rr.keys():
        issues.append(dict(field="cases", kind="support",
                           python_only=sorted(py.keys() - rr.keys()),
                           R_only=sorted(rr.keys() - py.keys())))
    for name in sorted(py.keys() & rr.keys()):
        p, r = py[name], rr[name]
        boot = name.startswith("bootstrap/")
        ids = array(r["unit_ids"]).ravel()
        # benchmark_v043.dataset() uses exactly these IDs; R has explicitly
        # aligned all IF rows to these IDs before export.
        if not np.array_equal(ids, np.arange(1, len(ids) + 1)):
            issues.append(dict(case=name, field="unit_ids", kind="support"))
            continue
        pk = list(zip(array(p["fit"]["group"]), array(p["fit"]["t"])))
        rk = list(zip(array(r["fit"]["group"]), array(r["fit"]["t"])))
        if len(pk) != len(set(pk)) or len(rk) != len(set(rk)) or set(pk) != set(rk):
            issues.append(dict(case=name, field="cells", kind="support"))
            continue
        order = [rk.index(cell) for cell in pk]
        for field in (["att"] if boot else ["att", "se", "l_se", "u_se"]):
            check(name, field, p["fit"][field], array(r["fit"][field])[order])
        if not boot:
            for field in ("c", "W", "Wpval"):
                check(name, field, p["fit"][field], r["fit"][field])
        check(name, "inffunc", p["fit"]["inffunc"],
              array(r["fit"]["inffunc"])[:, order], influence=True, matrix=True)
        if p["aggs"].keys() != r["aggs"].keys():
            issues.append(dict(case=name, field="aggs", kind="support"))
            continue
        for kind in p["aggs"]:
            pa, ra = p["aggs"][kind], r["aggs"][kind]
            fields = ["overall_att", "egt", "att_egt"]
            if not boot:
                fields += ["overall_se", "se_egt", "crit_val_egt"]
            for field in fields:
                check(name, kind + "." + field, pa[field], ra[field])
            if pa["inf_function"].keys() != ra["inf_function"].keys():
                influence_issues.append(dict(case=name, field=kind + ".inf_function",
                                             kind="support"))
                continue
            for field in pa["inf_function"]:
                check(name, kind + "." + field, pa["inf_function"][field],
                      ra["inf_function"][field], influence=True, matrix=True)
    result = dict(public_issues=issues, influence_issues=influence_issues,
                  max_abs=maxima, counts=counts, atol=1e-7, rtol=1e-6,
                  Wald_atol=1e-5, Wald_rtol=1e-5,
                  bootstrap="Only ATT and full IF compared; independent native draws make SE/band equality invalid",
                  warning_scope="Warnings are preserved in each snapshot; cross-language text is not a parity gate")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2))
    print(json.dumps(dict(cases=len(py.keys() & rr.keys()), public_issues=len(issues),
                          influence_issues=len(influence_issues), max_abs=maxima), indent=2))
    if issues or (args.require_influence and influence_issues):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
