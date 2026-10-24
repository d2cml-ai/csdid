r"""Serial before/after public-API benchmark, with complete numerical comparisons.

Run from any directory, with an aligned pre-performance source snapshot::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
      csdid/test_csdid/benchmark_v043.py --baseline /path/to/pre-speed-tree \
      --output /tmp/csdid-benchmark --reps 5 --profiles

Both faster_mode values are explicit. Dataset generation, imports, warmups,
serialization and comparisons are outside timers. Construction/preprocessing,
fit and four public aggte calls are included. Sources run in separate serial
processes; source order and mode order alternate across rounds. Every measured
result is compared with its baseline, including complete influence functions,
finite masks, warnings and the final NumPy RNG state. There is no speed threshold.
Output contains raw runs, median/range summaries, source hashes, environment,
reference arrays, comparison diagnostics and optional cProfile files. Use the
same interpreter and BLAS settings for both sources, and avoid competing jobs.
"""

import argparse
import contextlib
import cProfile
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import pstats
import subprocess
import sys
import time
import warnings

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
CASES = {
    "analytic": (8000, 10, False, False),
    "clustered": (8000, 10, False, True),
    "bootstrap": (2000, 8, True, False),
}
FIT_FIELDS = ("att", "se", "c", "l_se", "u_se")
AGG_FIELDS = ("overall_att", "overall_se", "egt", "att_egt", "se_egt", "crit_val_egt")


def dataset(n, periods):
    """Fixed weighted panel matching the original performance audit workload."""
    rng = np.random.default_rng(71828)
    groups = [0] + list(range(2, periods + 1))
    cohorts = np.repeat(groups, n // len(groups))
    cohorts = np.pad(cohorts, (0, n - len(cohorts)))
    ids = np.arange(1, n + 1)
    x, alpha = rng.standard_normal(n), rng.standard_normal(n) * .3
    frames = []
    for period in range(1, periods + 1):
        y = alpha + .3 * x + .1 * period + rng.standard_normal(n)
        y += ((cohorts > 0) & (period >= cohorts))
        frames.append(pd.DataFrame(dict(id=ids, period=period, G=cohorts, X=x,
                                        Y=y, cluster=ids % 400, wt=np.exp(.15 * x))))
    return pd.concat(frames, ignore_index=True)


def execute(attgt, data, fast, boot, clustered):
    np.random.seed(1234)
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        start = time.perf_counter()
        obj = attgt(yname="Y", tname="period", idname="id", gname="G", data=data,
                    control_group="nevertreated", xformla="Y~X", faster_mode=fast,
                    weights_name="wt", clustervar="cluster" if clustered else None,
                    biters=399, cband=boot)
        constructed = time.perf_counter()
        obj.fit("dr", bstrap=boot)
        fitted = time.perf_counter()
        fit = {key: obj.results[key] for key in FIT_FIELDS}
        fit.update(group=obj.MP["group"], t=obj.MP["t"], W=obj.MP["W"],
                   Wpval=obj.MP["Wpval"], inffunc=obj.MP["inffunc"]["inffunc"])
        aggs = {}
        for kind in ("simple", "group", "dynamic", "calendar"):
            obj.aggte(typec=kind, bstrap=boot, cband=boot)
            aggs[kind] = {key: obj.atte.get(key) for key in AGG_FIELDS}
            aggs[kind]["inf_function"] = obj.atte["inf_function"]
        end = time.perf_counter()
    timing = dict(construct_s=constructed-start, fit_s=fitted-constructed,
                  aggregation_s=end-fitted, total_s=end-start)
    result = dict(fit=fit, aggs=aggs,
                  warnings=[(w.category.__name__, str(w.message)) for w in caught],
                  rng_state=np.random.get_state())
    return result, timing


def save_snapshot(value, path):
    arrays = {}

    def encode(obj):
        if isinstance(obj, dict):
            return {str(key): encode(val) for key, val in obj.items()}
        if isinstance(obj, np.ndarray):
            key = f"array_{len(arrays)}"
            arrays[key] = obj
            return {"__array__": key}
        if isinstance(obj, (list, tuple)):
            return [encode(val) for val in obj]
        if isinstance(obj, np.generic):
            return obj.item()
        return obj

    metadata = encode(value)
    path.with_suffix(".json").write_text(json.dumps(metadata, allow_nan=True))
    np.savez_compressed(path.with_suffix(".npz"), **arrays)


def load_snapshot(path):
    with np.load(path.with_suffix(".npz"), allow_pickle=False) as arrays:
        def decode(obj):
            if isinstance(obj, dict):
                if set(obj) == {"__array__"}:
                    return arrays[obj["__array__"]].copy()
                return {key: decode(val) for key, val in obj.items()}
            if isinstance(obj, list):
                return [decode(val) for val in obj]
            return obj
        return decode(json.loads(path.with_suffix(".json").read_text()))


def compare(left, right, atol=0., rtol=0.):
    """Require symmetric keys and matching shapes/nonfinite masks before values."""
    issues, maxima = [], {}

    def visit(a, b, key):
        if isinstance(a, dict) or isinstance(b, dict):
            if not isinstance(a, dict) or not isinstance(b, dict) or a.keys() != b.keys():
                issues.append(dict(field=key, kind="keys", expected=list(a) if isinstance(a, dict) else type(a).__name__,
                                   actual=list(b) if isinstance(b, dict) else type(b).__name__))
                return
            for name in a:
                visit(a[name], b[name], f"{key}.{name}")
            return
        if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
            if not isinstance(a, (list, tuple)) or not isinstance(b, (list, tuple)) or len(a) != len(b):
                issues.append(dict(field=key, kind="sequence"))
                return
            for i, (aa, bb) in enumerate(zip(a, b)):
                visit(aa, bb, f"{key}[{i}]")
            return
        if a is None or b is None or isinstance(a, str) or isinstance(b, str):
            if a != b:
                issues.append(dict(field=key, kind="value", expected=a, actual=b))
            return
        a, b = np.asarray(a), np.asarray(b)
        if a.shape != b.shape:
            issues.append(dict(field=key, kind="shape", expected=a.shape, actual=b.shape))
            return
        for mask_name, mask_fn in (("finite", np.isfinite), ("nan", np.isnan),
                                   ("posinf", np.isposinf), ("neginf", np.isneginf)):
            if not np.array_equal(mask_fn(a), mask_fn(b)):
                issues.append(dict(field=key, kind=mask_name+"_mask"))
                return
        finite = np.isfinite(a)
        if not np.any(finite):
            return
        aa, bb = a[finite].astype(float), b[finite].astype(float)
        maximum = float(np.max(np.abs(aa-bb)))
        maxima[key] = maximum
        if not np.allclose(aa, bb, atol=atol, rtol=rtol):
            issues.append(dict(field=key, kind="numeric", max_abs=maximum,
                               max_rel=float(np.max(np.abs(aa-bb)/np.maximum(np.abs(aa), 1e-300)))))

    visit(left, right, "result")
    return dict(issues=issues, max_abs=maxima, atol=atol, rtol=rtol)


def source_hashes(repo):
    return {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted((repo / "csdid").rglob("*.py"))
            if "test_csdid" not in path.parts}


def worker(args):
    sys.path.insert(0, str(args.repo.resolve()))
    from csdid.att_gt import ATTgt

    result, timings = {}, {}
    modes = (True, False) if args.reverse else (False, True)
    for name in args.cases.split(","):
        n, periods, boot, clustered = CASES[name]
        data = dataset(n, periods)
        for fast in modes:
            execute(ATTgt, data, fast, boot, clustered)  # untimed warmup
        for fast in modes:
            key = f"{name}/fast={fast}"
            if args.profiles:
                profile = cProfile.Profile()
                profile.enable()
                execute(ATTgt, data, fast, boot, clustered)
                profile.disable()
                profile_path = args.output / f"{name}-fast{int(fast)}"
                profile.dump_stats(str(profile_path) + ".prof")
                with profile_path.with_suffix(".txt").open("w") as stream:
                    pstats.Stats(profile, stream=stream).strip_dirs().sort_stats("cumulative").print_stats(50)
            result[key], timings[key] = execute(ATTgt, data, fast, boot, clustered)
    save_snapshot(result, args.output / "snapshot")
    (args.output / "timings.json").write_text(json.dumps(timings, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--repo", type=Path, default=REPO)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--cases", default=",".join(CASES))
    parser.add_argument("--profiles", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--reverse", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.worker:
        worker(args)
        return
    if args.baseline is None:
        parser.error("--baseline must name an aligned pre-performance source tree")
    if args.reps < 1 or any(case not in CASES for case in args.cases.split(",")):
        parser.error("--reps must be positive and --cases must select known workloads")
    sources = dict(baseline=args.baseline.resolve(), candidate=args.repo.resolve())
    numpy_config = io.StringIO()
    with contextlib.redirect_stdout(numpy_config):
        np.show_config()
    environment = dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
                       cpu_count=os.cpu_count(), numpy_config=numpy_config.getvalue(),
                       packages={name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy", "statsmodels")},
                       thread_env={name: os.environ.get(name) for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
                       source_paths={name: str(path) for name, path in sources.items()},
                       source_hashes={name: source_hashes(path) for name, path in sources.items()})
    (args.output / "environment.json").write_text(json.dumps(environment, indent=2))
    records, comparisons, reference = {name: [] for name in sources}, [], None
    for rep in range(args.reps):
        for name in (("baseline", "candidate") if rep % 2 == 0 else ("candidate", "baseline")):
            run_dir = args.output / f"{rep:02d}-{name}"
            command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--repo", str(sources[name]),
                       "--output", str(run_dir), "--cases", args.cases]
            if rep % 2:
                command.append("--reverse")
            if args.profiles and rep == 0:
                command.append("--profiles")
            subprocess.run(command, check=True)
            snapshot = load_snapshot(run_dir / "snapshot")
            if reference is None:
                reference = snapshot
            check = compare(reference, snapshot)
            comparisons.append(dict(source=name, round=rep, **check))
            (args.output / "comparisons.json").write_text(json.dumps(comparisons, indent=2))
            if check["issues"]:
                raise AssertionError(f"{name}, round {rep}: {len(check['issues'])} parity issues; see comparisons.json")
            records[name].append(json.loads((run_dir / "timings.json").read_text()))
            if rep > 0:
                # The first pair preserves the full reference arrays. Later
                # rounds retain their comparisons and raw timings without
                # duplicating hundreds of MB of identical influence functions.
                for suffix in (".json", ".npz"):
                    (run_dir / "snapshot").with_suffix(suffix).unlink()
            print(f"{name}, round {rep+1}: exact full-output parity", flush=True)
    summary = {}
    for name, runs in records.items():
        summary[name] = {}
        for case in runs[0]:
            summary[name][case] = {}
            for phase in runs[0][case]:
                values = [run[case][phase] for run in runs]
                summary[name][case][phase] = dict(median_s=float(np.median(values)), min_s=min(values), max_s=max(values))
    (args.output / "results.json").write_text(json.dumps(dict(raw=records, summary=summary), indent=2))
    for case in summary["baseline"]:
        old = summary["baseline"][case]["total_s"]["median_s"]
        new = summary["candidate"][case]["total_s"]["median_s"]
        print(f"{case}: {old*1000:.1f} -> {new*1000:.1f} ms ({old/new:.2f}x)")


if __name__ == "__main__":
    main()
