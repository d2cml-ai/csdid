"""Broader serial performance audit against a fixed source snapshot.

Run with the same interpreter/dependencies and one BLAS thread, for example::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python \
      csdid/test_csdid/benchmark_speed_audit.py --baseline /tmp/csdid-before \
      --output /tmp/csdid-speed-audit --reps 5 --profiles

Each source/case runs in its own process, with untimed data generation, warmup,
and optional profiling. Timers cover construction, fit, and all four public
aggregations. Source order alternates. Complete influence functions, public
results, warnings and RNG states must agree exactly with the baseline; timing
alone is never a correctness check. Peak RSS includes the interpreter, inputs,
warmup and result arrays, but excludes snapshot serialization. It is a process
high-water mark, not an allocation measurement for the estimator alone.
"""
import argparse
import contextlib
import cProfile
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import pstats
import resource
import subprocess
import sys
import time
import warnings

import numpy as np

from benchmark_v043 import (
    AGG_FIELDS, FIT_FIELDS, compare, dataset, load_snapshot, save_snapshot,
    source_hashes,
)

REPO = Path(__file__).resolve().parents[2]
CASES = {
    "panel_dr": {},
    "panel_ipw": {"method": "ipw"},
    "panel_reg": {"method": "reg"},
    "intercept": {"covariates": 0, "weighted": False},
    "varying_inputs": {"varying": True},
    "notyet": {"control": "notyettreated"},
    "universal": {"base": "universal"},
    "wide": {"covariates": 12},
    "long": {"n": 4000, "periods": 24},
    "large": {"n": 50000},
    "rc": {"n": 2000, "periods": 8, "design": "rc"},
    "unbalanced": {"n": 2000, "periods": 8, "design": "unbalanced"},
    "bootstrap": {"n": 2000, "periods": 8, "boot": True},
    "large_bootstrap": {"n": 20000, "periods": 8, "boot": True, "biters": 999},
    "few_clusters": {"n": 6000, "periods": 8, "boot": True,
                     "clusters": 5, "biters": 199},
    "many_clusters": {"n": 6000, "periods": 8, "boot": True,
                      "clusters": 3000, "biters": 199},
}


def specification(name):
    return dict(n=8000, periods=10, method="dr", covariates=1, weighted=True,
                varying=False, control="nevertreated", base="varying",
                design="panel", boot=False, clusters=None, biters=399) | CASES[name]


def make_data(spec):
    data = dataset(spec["n"], spec["periods"])
    rng = np.random.default_rng(190421)
    if spec["varying"]:
        data["X"] += rng.normal(scale=.2, size=len(data))
        data["wt"] *= np.exp(rng.normal(scale=.1, size=len(data)))
    for j in range(1, spec["covariates"]):
        data[f"X{j}"] = np.tile(rng.normal(size=spec["n"]), spec["periods"])
    if spec["clusters"] is not None:
        data["cluster"] = data["id"] % spec["clusters"]
    if spec["design"] == "rc":
        data["id"] = np.arange(1, len(data) + 1)
        data["X"] += rng.normal(scale=.2, size=len(data))
        data["Y"] += rng.normal(scale=.2, size=len(data))
    elif spec["design"] == "unbalanced":
        data = data.loc[rng.random(len(data)) > .12].copy()
    return data


def execute(attgt, data, spec, fast):
    covs = ["X"] + [f"X{j}" for j in range(1, spec["covariates"])]
    formula = "Y~" + ("+".join(covs) if spec["covariates"] else "1")
    np.random.seed(1234)
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        start = time.perf_counter()
        obj = attgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                    control_group=spec["control"], xformla=formula, faster_mode=fast,
                    weights_name="wt" if spec["weighted"] else None,
                    clustervar="cluster" if spec["clusters"] else None,
                    panel=spec["design"] != "rc",
                    allow_unbalanced_panel=spec["design"] == "unbalanced",
                    biters=spec["biters"], cband=spec["boot"])
        constructed = time.perf_counter()
        obj.fit(spec["method"], base_period=spec["base"], bstrap=spec["boot"])
        fitted = time.perf_counter()
        fit = {key: obj.results[key] for key in FIT_FIELDS}
        fit.update(group=obj.MP["group"], t=obj.MP["t"], W=obj.MP["W"],
                   Wpval=obj.MP["Wpval"], inffunc=obj.MP["inffunc"]["inffunc"])
        aggs = {}
        for kind in ("simple", "group", "dynamic", "calendar"):
            obj.aggte(typec=kind, bstrap=spec["boot"], cband=spec["boot"])
            aggs[kind] = {key: obj.atte.get(key) for key in AGG_FIELDS}
            aggs[kind]["inf_function"] = obj.atte["inf_function"]
        end = time.perf_counter()
    return (dict(fit=fit, aggs=aggs,
                 warnings=[(w.category.__name__, str(w.message)) for w in caught],
                 rng_state=np.random.get_state()),
            dict(construct_s=constructed-start, fit_s=fitted-constructed,
                 aggregation_s=end-fitted, total_s=end-start))


def worker(args):
    sys.path.insert(0, str(args.repo.resolve()))
    from csdid.att_gt import ATTgt

    spec = specification(args.case)
    data = make_data(spec)
    execute(ATTgt, data, spec, args.fast)  # untimed warmup
    if args.profiles:
        profile = cProfile.Profile()
        profile.runcall(execute, ATTgt, data, spec, args.fast)
        profile.dump_stats(str(args.output / "profile.prof"))
        with (args.output / "profile.txt").open("w") as stream:
            pstats.Stats(profile, stream=stream).strip_dirs().sort_stats("cumulative").print_stats(60)
    result, timings = execute(ATTgt, data, spec, args.fast)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    timings["peak_rss_bytes"] = int(rss if sys.platform == "darwin" else rss * 1024)
    (args.output / "timings.json").write_text(json.dumps(timings, indent=2))
    save_snapshot(result, args.output / "snapshot")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--repo", type=Path, default=REPO)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--cases", default=",".join(CASES))
    parser.add_argument("--modes", default="true", help="true, false, or true,false")
    parser.add_argument("--profiles", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--case", choices=CASES, help=argparse.SUPPRESS)
    parser.add_argument("--fast", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.worker:
        worker(args)
        return
    cases = args.cases.split(",")
    modes = args.modes.split(",")
    if args.baseline is None or args.reps < 1 or any(c not in CASES for c in cases):
        parser.error("Supply --baseline, positive --reps and valid --cases")
    if any(m not in ("true", "false") for m in modes):
        parser.error("--modes must contain true and/or false")
    sources = dict(baseline=args.baseline.resolve(), candidate=args.repo.resolve())
    numpy_config = io.StringIO()
    with contextlib.redirect_stdout(numpy_config):
        np.show_config()
    environment = dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
                       cpu_count=os.cpu_count(), numpy_config=numpy_config.getvalue(),
                       packages={p: importlib.metadata.version(p) for p in ("numpy", "pandas", "scipy", "statsmodels", "drdid")},
                       threads={p: os.environ.get(p) for p in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
                       sources={k: str(v) for k, v in sources.items()},
                       source_hashes={k: source_hashes(v) for k, v in sources.items()},
                       cases={c: specification(c) for c in cases}, modes=modes,
                       benchmark_sha256=__import__("hashlib").sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output / "environment.json").write_text(json.dumps(environment, indent=2))
    records, checks = [], []
    for case in cases:
        for mode in modes:
            reference = None
            for rep in range(args.reps):
                for source in (("baseline", "candidate") if rep % 2 == 0 else ("candidate", "baseline")):
                    run_dir = args.output / f"{case}-fast{mode}-{rep}-{source}"
                    command = [sys.executable, str(Path(__file__).resolve()), "--worker",
                               "--repo", str(sources[source]), "--output", str(run_dir), "--case", case]
                    if mode == "true":
                        command.append("--fast")
                    if args.profiles and rep == 0:
                        command.append("--profiles")
                    subprocess.run(command, check=True)
                    snapshot = load_snapshot(run_dir / "snapshot")
                    if reference is None:
                        reference = snapshot
                    check = compare(reference, snapshot)
                    checks.append(dict(case=case, fast=mode, source=source, rep=rep,
                                       issues=check["issues"], max_abs=max(check["max_abs"].values(), default=0.)))
                    (args.output / "checks.json").write_text(json.dumps(checks, indent=2))
                    if check["issues"]:
                        raise AssertionError(f"Full-output parity failure: {case}, {mode}, {source}; see checks.json")
                    records.append(dict(case=case, fast=mode, source=source, rep=rep,
                                        **json.loads((run_dir / "timings.json").read_text())))
                    (args.output / "raw.json").write_text(json.dumps(records, indent=2))
                    if rep > 0:
                        for suffix in (".json", ".npz"):
                            (run_dir / "snapshot").with_suffix(suffix).unlink()
            medians = {source: float(np.median([r["total_s"] for r in records
                       if r["case"] == case and r["fast"] == mode and r["source"] == source]))
                       for source in sources}
            print(f"{case}/fast={mode}: {medians['baseline']:.4f} -> {medians['candidate']:.4f}s "
                  f"({medians['baseline']/medians['candidate']:.2f}x); exact full-output parity", flush=True)
    summary = []
    for case in cases:
        for mode in modes:
            row = dict(case=case, fast=mode)
            for source in sources:
                matching = [r for r in records if (r["case"], r["fast"], r["source"]) == (case, mode, source)]
                row[source] = {key: dict(median=float(np.median([r[key] for r in matching])),
                                         min=min(r[key] for r in matching), max=max(r[key] for r in matching))
                               for key in ("construct_s", "fit_s", "aggregation_s", "total_s", "peak_rss_bytes")}
            row["speedup"] = row["baseline"]["total_s"]["median"] / row["candidate"]["total_s"]["median"]
            summary.append(row)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
