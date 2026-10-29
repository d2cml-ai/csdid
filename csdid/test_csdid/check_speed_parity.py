"""Compare complete outputs on 84 configurations against an isolated source tree.

Uses the adversarial datasets and public-API workload in test_performance_parity:
shuffled large IDs, factors, varying weights/covariates, irregular dates, panel,
unbalanced panel and repeated cross sections, both controls/bases/modes, all
three estimators, and seeded bootstrap with and without extra clustering.
Run with --baseline /path/to/source --output /path/to/evidence. No timing claims.
"""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import sys

from benchmark_v043 import compare, load_snapshot, save_snapshot, source_hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.worker:
        sys.path.insert(0, str(args.repo.resolve()))
        from test_performance_parity import _run
        results = {}
        for design, control, method, base, fast in itertools.product(
                ("panel", "unbalanced", "rc"), ("nevertreated", "notyettreated"),
                ("dr", "reg", "ipw"), ("varying", "universal"), (False, True)):
            key = f"{design}/{control}/{method}/{base}/fast={fast}"
            results[key] = _run(design, control, method, base, mode=fast)
        for design, clustered, fast in itertools.product(
                ("panel", "unbalanced", "rc"), (False, True), (False, True)):
            key = f"bootstrap/{design}/clustered={clustered}/fast={fast}"
            results[key] = _run(design, "nevertreated", base="universal", mode=fast,
                                bootstrap=True, clustered=clustered)
        save_snapshot(results, args.output / "snapshot")
        return
    if args.baseline is None:
        parser.error("--baseline is required")
    sources = dict(baseline=args.baseline, candidate=args.repo)
    for name, source in sources.items():
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker",
                        "--repo", str(source), "--output", str(args.output / name)], check=True)
    baseline = load_snapshot(args.output / "baseline" / "snapshot")
    candidate = load_snapshot(args.output / "candidate" / "snapshot")
    check = compare(baseline, candidate)
    check.update(configurations=len(baseline), aggregation_requests=4*len(baseline),
                 source_hashes={name: source_hashes(path) for name, path in sources.items()})
    (args.output / "comparison.json").write_text(json.dumps(check, indent=2))
    if check["issues"]:
        raise AssertionError(f"{len(check['issues'])} differences; see comparison.json")
    print(f"{len(baseline)} configurations, {4*len(baseline)} aggregations: exact full-output parity")


if __name__ == "__main__":
    main()
