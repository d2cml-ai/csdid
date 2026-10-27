"""Replay permanent R references with an isolated installed wheel.

Build a wheel, install it with pip --no-deps --target /tmp/installed-csdid, then
run this script from outside the source checkout with --installed pointing to
that target. The current test helpers supply reference checks only; every
csdid runtime import must resolve under the installed target.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--installed", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--tests", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    target, tests = args.installed.resolve(), args.tests.resolve()
    sys.path.insert(0, str(target))
    import csdid
    assert Path(csdid.__file__).resolve().parent.parent == target
    reference = tests / "r_ref/current"
    os.environ["CSDID_R_REFERENCE_DIR"] = str(reference)
    sys.path.insert(0, str(tests))
    spec = importlib.util.spec_from_file_location("wheel_reference_checks", tests / "test_r_current_parity.py")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    with gzip.open(reference / "references.json.gz", "rt") as stream:
        refs = json.load(stream)
    results = []
    for case in checks.CASES:
        for fast in (False, True):
            actual = checks._run_case(case, fast)
            checks._assert_case_output(actual, refs[case["name"]])
            results.append(dict(case=case["name"], faster_mode=fast, passed=True))
    loaded = {name: str(Path(module.__file__).resolve()) for name, module in sys.modules.items()
              if (name == "csdid" or name.startswith("csdid.")) and getattr(module, "__file__", None)}
    assert loaded and all(target in Path(path).parents for path in loaded.values())
    with zipfile.ZipFile(args.wheel) as archive:
        members = archive.namelist()
        assert not any("test_csdid" in path or path.endswith((".R", ".csv", ".json.gz")) for path in members)
        source_hashes = {name: hashlib.sha256(archive.read(name)).hexdigest()
                         for name in members if name.startswith("csdid/") and name.endswith(".py")}
    args.output.write_text(json.dumps(dict(
        audit_date="2026-10-26", comparisons=results, passed=len(results),
        loaded_runtime_paths=loaded, wheel_sha256=hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
        wheel_size_bytes=args.wheel.stat().st_size, wheel_runtime_sha256=source_hashes,
        test_and_R_reference_artifacts_excluded=True,
        scope="Same pinned R numerical gate, including full influence functions and shared-draw bootstrap, on installed wheel imports outside checkout.",
    ), indent=2) + "\n")
    print(f"Installed wheel: {len(results)} complete numerical workflows passed; runtime imports isolated.")


if __name__ == "__main__":
    main()
