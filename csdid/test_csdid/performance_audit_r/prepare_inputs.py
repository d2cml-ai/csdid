"""Export exactly benchmark_v043.dataset(), without running an estimator."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_v043 import CASES, dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, (n, periods, boot, clustered) in CASES.items():
        path = args.output / f"{name}.csv"
        data = dataset(n, periods)
        data.to_csv(path, index=False, float_format="%.17g")
        payload = path.read_bytes()
        manifest[name] = dict(n=n, periods=periods, rows=len(data),
                              bootstrap=boot, clustered=clustered,
                              filename=path.name,
                              sha256=hashlib.sha256(payload).hexdigest(),
                              md5=hashlib.md5(payload).hexdigest())
    (args.output / "inputs.json").write_text(json.dumps(manifest, indent=2))
    print(f"Exported three benchmark datasets to {args.output}")


if __name__ == "__main__":
    main()
