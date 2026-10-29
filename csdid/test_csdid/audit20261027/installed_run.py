#!/usr/bin/env python3
"""Run an R comparison using only the wheels from a passed joint-package gate.

Invoke with that gate's fresh Python and ``-I``. Audit helpers may read source
fixtures and insert their helper directory on sys.path; preloading both runtime
packages and checking every resulting module origin prevents a checkout from
silently replacing the implementation under review.
"""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import runpy
import sys


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--joint", type=Path, required=True)
    parser.add_argument("--audit", choices=("compare_weight_scale", "compare_sparse_float32"), required=True)
    parser.add_argument("--guard-receipt", type=Path, required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    joint = args.joint.resolve()
    receipt_file = joint / "receipt.json"
    receipt = json.loads(receipt_file.read_text())
    if receipt["status"] != "passed":
        parser.error("The joint installed-wheel gate must have passed")
    guard_file = args.guard_receipt.resolve()
    if guard_file.exists():
        parser.error("The guard receipt must be new")
    helper = joint / "_joint_gate.py"
    if digest(helper) != receipt["driver_sha256"]:
        parser.error("Installed-wheel guard differs from the reviewed joint driver")
    spec = importlib.util.spec_from_file_location("installed_joint_origin_guard", helper)
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    expected = {}
    for identity in receipt["wheels"].values():
        expected.update(identity["runtime_files_sha256"])
    site_packages = receipt["installed"]["site_packages"]
    result = dict(schema_version=1, status="running", joint_receipt_sha256=digest(receipt_file),
                  wrapper_sha256=digest(__file__), executable=sys.executable,
                  isolated_mode=bool(sys.flags.isolated), audit=args.audit)
    returncode = 2
    try:
        if not sys.flags.isolated:
            raise ValueError("Use the joint environment's Python with -I")
        # These imports occur before the audit can add a source helper path.
        import csdid
        import drdid
        gate.verify_installed_files(site_packages, expected)
        result["origins_before"] = gate.collect_and_validate_origins(site_packages, expected)
        result["packages"] = {name: importlib.metadata.version(name) for name in ("csdid", "drdid")}
        for name, version in result["packages"].items():
            if version != receipt["wheels"][name]["version"]:
                raise ValueError("Installed package version differs from joint gate: " + name)
        audit = Path(__file__).with_name(args.audit + ".py")
        arguments = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
        result["audit_sha256"] = digest(audit)
        result["audit_arguments"] = arguments
        sys.path.insert(0, str(audit.parent))
        sys.argv = [str(audit)] + arguments
        try:
            runpy.run_path(str(audit), run_name="__main__")
        except SystemExit as exit_status:
            returncode = 0 if exit_status.code is None else int(exit_status.code)
        else:
            returncode = 0
        gate.verify_installed_files(site_packages, expected)
        result["origins_after"] = gate.collect_and_validate_origins(site_packages, expected)
        if returncode == 0:
            output_at = arguments.index("--output")
            audit_result = Path(arguments[output_at + 1]) / "result.json"
            numeric_receipt = json.loads(audit_result.read_text())
            if numeric_receipt["status"] != "passed":
                raise ValueError("Numerical audit did not produce a passing receipt")
            if (numeric_receipt["source_identity"]["csdid"]["source_files_sha256"]
                    != receipt["sources"]["csdid"]["runtime_files_sha256"]):
                raise ValueError("Audited CSDID source differs from the jointly installed wheel")
            result["numeric_receipt_sha256"] = digest(audit_result)
        result["status"] = "passed" if returncode == 0 else "failed"
    except Exception as error:
        result["status"] = "failed"
        result["error"] = "{}: {}".format(type(error).__name__, error)
        returncode = 2
    result["returncode"] = returncode
    guard_file.parent.mkdir(parents=True, exist_ok=True)
    guard_file.write_text(json.dumps(result, indent=2) + "\n")
    return returncode


if __name__ == "__main__":
    sys.exit(main())
