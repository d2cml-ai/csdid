"""Restricted 2x2 ATT/SE parity and warm timing; not a staggered-DiD ranking.

Run from the csdid checkout with PYTHONPATH set to that checkout, using the
separate PyFixest comparison environment. Dataset construction is outside timing;
the first-difference PyFixest variant includes its pivot inside timing.
"""
import argparse
import gc
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time

# Keep thread configuration fixed before loading any numerical library.
for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
import pandas as pd
import pyfixest as pf
from threadpoolctl import threadpool_info, threadpool_limits

from csdid.att_gt import ATTgt


def dataset(n):
    rng = np.random.default_rng(730 + n)
    treated = rng.random(n) < 0.37
    unit_effect = rng.normal(size=n)
    baseline = unit_effect + rng.normal(size=n)
    followup = unit_effect + 0.4 + rng.normal(size=n) + treated * (
        1.0 + 0.3 * rng.normal(size=n)
    )
    frame = pd.DataFrame({
        "id": np.repeat(np.arange(n), 2),
        "t": np.tile([1, 2], n),
        "g": np.repeat(np.where(treated, 2, 0), 2),
        "y": np.column_stack([baseline, followup]).ravel(),
    })
    frame["treatment"] = ((frame.g == 2) & (frame.t == 2)).astype(int)
    delta = followup - baseline
    att = delta[treated].mean() - delta[~treated].mean()
    # Direct analytical influence-function variance, with no finite-sample
    # adjustment, is the common reference for both regression formulations.
    se = np.sqrt(
        np.sum((delta[treated] - delta[treated].mean()) ** 2) / treated.sum() ** 2
        + np.sum((delta[~treated] - delta[~treated].mean()) ** 2) / (~treated).sum() ** 2
    )
    return frame, np.array([att, se])


def csdid_fit(frame, faster):
    fit = ATTgt(
        yname="y", tname="t", idname="id", gname="g", data=frame,
        xformla="~1", cband=False, faster_mode=faster,
    ).fit(est_method="dr", bstrap=False)
    return np.array([fit.results["att"][0], fit.results["se"][0]])


def pyfixest_twfe(frame):
    fit = pf.feols(
        "y ~ treatment | id + t", data=frame, vcov={"CRV1": "id"},
        ssc=pf.ssc(k_adj=False, G_adj=False),
        demeaner=pf.MapDemeaner(backend="rust"),
    )
    return np.array([fit.coef()["treatment"], fit.se()["treatment"]])


def pyfixest_first_difference(frame):
    wide = frame.pivot(index="id", columns="t", values="y")
    treated = frame.loc[frame.t == 1].set_index("id").g.reindex(wide.index).eq(2)
    differences = pd.DataFrame({"dy": wide[2] - wide[1], "D": treated.astype(int)})
    fit = pf.feols(
        "dy ~ D", data=differences, vcov="hetero",
        ssc=pf.ssc(k_adj=False, G_adj=False),
    )
    return np.array([fit.coef()["D"], fit.se()["D"]])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--units", default="2000,10000,50000")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--parity-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    methods = {
        "csdid_standard_dr": lambda frame: csdid_fit(frame, False),
        "csdid_faster_dr": lambda frame: csdid_fit(frame, True),
        "pyfixest_twfe_crv0": pyfixest_twfe,
        "pyfixest_first_difference_hc0": pyfixest_first_difference,
    }
    packages = ["pyfixest", "numpy", "pandas", "scipy", "statsmodels", "drdid"]
    payload = {
        "scope": "Balanced two-period panel, one treatment cohort, never-treated controls, no covariates, no weights; equal ATT and unadjusted analytical SE only. Not a staggered-DiD ranking.",
        "python": platform.python_version(), "platform": platform.platform(),
        "packages": {name: importlib.metadata.version(name) for name in packages},
        "audited_pyfixest_master": "4245747fad805ec32c013ea4c0170a020e9ceded",
        "csdid_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "source_sha256": {
            path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
            for path in sorted(str(path) for path in Path("csdid").rglob("*.py")
                               if "test_csdid" not in path.parts)
        },
        "thread_environment": {name: os.environ.get(name) for name in (
            "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "RAYON_NUM_THREADS")},
        "rows": [],
    }
    with threadpool_limits(limits=1):
        payload["threadpools"] = threadpool_info()
        for n in map(int, args.units.split(",")):
            frame, reference = dataset(n)
            original = frame.copy(deep=True)
            # One untimed warm-up/parity check for every method.
            observed = {name: method(frame) for name, method in methods.items()}
            for values in observed.values():
                np.testing.assert_allclose(values, reference, rtol=1e-9, atol=1e-11)
            pd.testing.assert_frame_equal(frame, original)
            samples = {name: [] for name in methods}
            if not args.parity_only:
                # Rotate method order to distribute drift across implementations.
                names = list(methods)
                for repeat in range(args.repeats):
                    offset = repeat % len(names)
                    for name in names[offset:] + names[:offset]:
                        gc.collect()
                        start = time.perf_counter()
                        values = methods[name](frame)
                        samples[name].append(time.perf_counter() - start)
                        np.testing.assert_allclose(values, reference, rtol=1e-9, atol=1e-11)
            for name in methods:
                row = {
                    "units": n, "rows": len(frame), "method": name,
                    "att": float(observed[name][0]), "se": float(observed[name][1]),
                    "max_abs_reference_error": float(np.max(np.abs(observed[name] - reference))),
                    "seconds": samples[name],
                    "median_seconds": statistics.median(samples[name]) if samples[name] else None,
                }
                payload["rows"].append(row)
                print(json.dumps(row), flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
