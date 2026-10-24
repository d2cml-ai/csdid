"""Compare csdid to supplied R source checkouts without updating reviewed fixtures.

The GitHub workflow checks out the latest default branches before invoking this
script. Local runs use exactly the supplied revisions and never fetch or modify
the source repositories. Reports identify those revisions, not an unverified
claim that a local checkout is current.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone


REPO = Path(__file__).resolve().parents[3]
TESTS = REPO / "csdid" / "test_csdid"
PINNED = TESTS / "r_ref" / "current"
INPUTS = ("panel.csv", "rc.csv", "unbalanced.csv",
          "multipliers_unit.csv", "multipliers_cluster.csv")
THREADS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS")


def _git(path, *args):
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


def source_identity(path, package):
    """Record a clean source tree and its package version before execution."""
    commit = _git(path, "rev-parse", "HEAD")
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("Invalid Git commit for " + package)
    dirty = _git(path, "status", "--porcelain", "--untracked-files=no")
    untracked = _git(path, "ls-files", "--others", "--exclude-standard", "--",
                     "DESCRIPTION", "NAMESPACE", "R", "src", "tests")
    if dirty or untracked:
        raise ValueError(package + ": source checkout must be clean")
    description = (path / "DESCRIPTION").read_text()
    names = re.findall(r"^Package:\s*(\S+)\s*$", description, re.MULTILINE)
    versions = re.findall(r"^Version:\s*(\S+)\s*$", description, re.MULTILINE)
    if names != [package] or len(versions) != 1:
        raise ValueError(package + ": invalid DESCRIPTION identity")
    return {"commit": commit, "version": versions[0]}


def _hashes(directory):
    return {str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(directory.rglob("*")) if path.is_file()}


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def _command(command, output, environment, timeout=1200):
    """Retain full diagnostics even when a subprocess fails or times out."""
    with output.open("w") as log:
        try:
            result = subprocess.run(command, cwd=REPO, env=environment,
                                    stdout=log, stderr=subprocess.STDOUT,
                                    timeout=timeout, check=False)
            return result.returncode
        except (OSError, subprocess.TimeoutExpired) as error:
            log.write("\n" + str(error) + "\n")
            return -1


def _prepare_output(output, did, drdid):
    for protected in (PINNED.resolve(), did, drdid):
        if output == protected or protected in output.parents or output in protected.parents:
            raise ValueError("Output must be separate from references and R sources")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Output directory must be new or empty: " + str(output))
    output.mkdir(parents=True, exist_ok=True)


def run(did, drdid, output):
    did, drdid, output = (Path(path).resolve() for path in (did, drdid, output))
    _prepare_output(output, did, drdid)
    report = {"schema_version": 1, "started_utc": datetime.now(timezone.utc).isoformat(),
              "status": "infrastructure_error", "exit_code": 2, "stages": {},
              "reference_policy": "Reviewed fixtures, inputs and tolerances are never updated",
              "python": sys.version, "source_checkouts": {"did": str(did), "DRDID": str(drdid)}}
    before = _hashes(PINNED)
    try:
        sources = {"did": source_identity(did, "did"),
                   "DRDID": source_identity(drdid, "DRDID")}
        lock = {"schema_version": 1, "sources": sources}
        report["sources"] = sources
        _write_json(output / "source-lock.json", lock)
        environment = os.environ.copy()
        # Do not inherit an override from a caller's previous comparison.
        environment.pop("CSDID_R_REFERENCE_DIR", None)
        environment.pop("CSDID_R_REFERENCE_SOURCES", None)
        environment.update({key: "1" for key in THREADS})
        environment["PYTHONPATH"] = str(REPO) + os.pathsep + environment.get("PYTHONPATH", "")
        audit_code = _command(
            [sys.executable, str(Path(__file__).with_name("audit.py")),
             "--did", str(did), "--drdid", str(drdid),
             "--output", str(output / "audit.json"), "--summary", str(output / "audit.md")],
            output / "audit.log", environment)
        report["stages"]["test_coverage_audit"] = audit_code
        if audit_code not in (0, 3):
            raise RuntimeError("Upstream test audit failed; see audit.log")
        references = output / "references"
        references.mkdir()
        for name in INPUTS + ("generate.R",):
            shutil.copy2(PINNED / name, references / name)
        generation_code = _command(
            ["Rscript", str(references / "generate.R"), "--did", str(did),
             "--drdid", str(drdid), "--source-lock", str(output / "source-lock.json")],
            output / "r-generation.log", environment)
        report["stages"]["r_generation"] = generation_code
        if generation_code != 0:
            raise RuntimeError("R reference generation failed; see r-generation.log")
        for name in INPUTS:
            if hashlib.sha256((references / name).read_bytes()).hexdigest() != before[name]:
                raise RuntimeError("Frozen comparison input changed: " + name)
        environment["CSDID_R_REFERENCE_DIR"] = str(references)
        environment["CSDID_R_REFERENCE_SOURCES"] = str(output / "source-lock.json")
        pytest_code = _command(
            [sys.executable, "-m", "pytest", str(TESTS / "test_r_current_parity.py"),
             "-q", "--tb=short", "--junitxml=" + str(output / "junit.xml")],
            output / "pytest.log", environment)
        report["stages"]["numerical_comparison"] = pytest_code
        if pytest_code == 1:
            report.update(status="numerical_mismatch", exit_code=1)
        elif pytest_code != 0:
            raise RuntimeError("Python comparison did not complete; see pytest.log")
        elif audit_code == 3:
            report.update(status="coverage_review_required", exit_code=3)
        else:
            report.update(status="passed", exit_code=0)
        # Freeze identity across the complete run, including external R code.
        if sources != {"did": source_identity(did, "did"),
                       "DRDID": source_identity(drdid, "DRDID")}:
            raise RuntimeError("R source revision changed during the run")
    except Exception as error:
        report.update(status="infrastructure_error", exit_code=2,
                      error="{}: {}".format(type(error).__name__, error))
    finally:
        report["reviewed_references_unchanged"] = before == _hashes(PINNED)
        if not report["reviewed_references_unchanged"]:
            report.update(status="infrastructure_error", exit_code=2,
                          error="Reviewed reference files changed during the run")
        report["finished_utc"] = datetime.now(timezone.utc).isoformat()
        _write_json(output / "report.json", report)
        lines = ["# R upstream comparison", "", "Status: **{}**".format(report["status"]), "",
                 "This run uses the recorded source revisions; reviewed Python references were not refreshed.", ""]
        for package, identity in report.get("sources", {}).items():
            lines.append("- {} {}: `{}`".format(package, identity["version"], identity["commit"]))
        lines.extend(["", "## Checks", ""])
        for stage, code in report["stages"].items():
            lines.append("- {}: exit {}".format(stage, code))
        if "error" in report:
            lines.extend(["", "Error: " + report["error"]])
        if (output / "audit.md").exists():
            lines.extend(["", (output / "audit.md").read_text()])
        lines.extend(["", "See the uploaded artifact for logs, JUnit results, source lock and generated R references.",
                      "Changed R tests require a coverage review; passing numerical cases do not port new assertions."])
        (output / "summary.md").write_text("\n".join(lines) + "\n")
    print("R upstream comparison: {} (exit {}); {}".format(
        report["status"], report["exit_code"], output / "summary.md"))
    return report["exit_code"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--did", type=Path, required=True)
    parser.add_argument("--drdid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args.did, args.drdid, args.output)
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
