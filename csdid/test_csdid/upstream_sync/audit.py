"""Compare local R Git trees with the committed numerical/coverage baseline.

This stdlib-only command never fetches sources or translates tests. Every
tracked Git-tree change requires human review, including runtime sources,
tests, helpers, fixtures, documentation and build metadata. Test and other
source changes remain separate evidence categories.
Exit codes: 0 = no tracked-tree drift, 3 = review required, 2 = invalid evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys


TEST_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = TEST_ROOT / "r_ref" / "current" / "manifest.json"
DEFAULT_CROSSWALKS = {
    "did": TEST_ROOT / "r_test_alignment_did.json",
    "DRDID": TEST_ROOT / "r_test_alignment_drdid.json",
}
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_COMMIT = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_DECLARATION = re.compile(r"^\s*(?:testthat::)?test_that\s*\(")


class AuditError(ValueError):
    """The available evidence does not support a drift determination."""


def _unique_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise AuditError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _read_json(path):
    try:
        result = json.loads(Path(path).read_text(encoding="utf-8"),
                            object_pairs_hook=_unique_keys)
    except (OSError, UnicodeError, ValueError) as exc:
        raise AuditError(f"Cannot read valid JSON from {path}: {exc}") from exc
    if not isinstance(result, dict) or type(result.get("schema_version")) is not int \
            or result["schema_version"] != 1:
        raise AuditError(f"Expected schema_version 1 object in {path}")
    return result


def _valid_path(path):
    return (isinstance(path, str) and bool(path)
            and not path.startswith("/")
            and all(part not in ("", ".", "..") for part in path.split("/")))


def _direct_test(path):
    path = PurePosixPath(path)
    return (path.parent == PurePosixPath("tests/testthat")
            and path.name.startswith("test") and path.name.endswith(".R"))


def _git(repo, *args):
    try:
        result = subprocess.run(["git", "-C", str(repo), *args],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                check=False)
    except OSError as exc:
        raise AuditError(f"Cannot run Git: {exc}") from exc
    if result.returncode:
        detail = result.stderr.decode("utf-8", errors="replace").strip()
        raise AuditError(f"Git {' '.join(args)} failed in {repo}: {detail}")
    return result.stdout


def _tree(repo, commit):
    result = {}
    for record in _git(repo, "ls-tree", "-rz", "--full-tree", commit).split(b"\0"):
        if not record:
            continue
        metadata, path = record.split(b"\t", 1)
        mode, kind, oid = metadata.decode("ascii").split()
        result[path.decode("utf-8", errors="surrogateescape")] = {
            "mode": mode, "kind": kind, "oid": oid,
        }
    return result


def _hash(content):
    return hashlib.sha256(content).hexdigest()


def _declarations(path, content):
    """Lexical declaration locations, not a count of executed R assertions."""
    if content is None or not path.endswith(".R"):
        return []
    return [{"line": number, "expression": line.strip()}
            for number, line in enumerate(content.decode("utf-8", errors="replace").splitlines(), 1)
            if _DECLARATION.match(line)]


def _baseline(path):
    manifest = _read_json(path)
    sources = manifest.get("sources")
    if not isinstance(sources, dict) or set(sources) != {"did", "DRDID"}:
        raise AuditError("Baseline sources must contain exactly did and DRDID")
    for package, entry in sources.items():
        if not isinstance(entry, dict) or not isinstance(entry.get("commit"), str) \
                or not _COMMIT.fullmatch(entry["commit"]):
            raise AuditError(f"Invalid baseline commit for {package}")
        hashes = entry.get("files_sha256")
        if not isinstance(hashes, dict) or not hashes:
            raise AuditError(f"Missing/nonempty files_sha256 required for {package}")
        for name, digest in hashes.items():
            if not _valid_path(name) or not isinstance(digest, str) or not _SHA256.fullmatch(digest):
                raise AuditError(f"Invalid baseline file hash for {package}: {name}")
    return sources


def _crosswalk(path, package, commit):
    crosswalk = _read_json(path)
    commit_key, entries_key, hash_key = (
        ("r_commit", "entries", "r_file_sha256") if package == "did"
        else ("r_revision", "files", "r_sha256")
    )
    if crosswalk.get(commit_key) != commit:
        raise AuditError(f"{package} crosswalk commit does not match numerical baseline")
    entries = crosswalk.get(entries_key)
    if not isinstance(entries, list):
        raise AuditError(f"Invalid {package} crosswalk {entries_key}")
    result = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise AuditError(f"Invalid {package} crosswalk entry")
        path = entry.get("r_file")
        digest = entry.get(hash_key)
        if not _valid_path(path) or not _direct_test(path) \
                or not isinstance(digest, str) or not _SHA256.fullmatch(digest) \
                or not isinstance(entry.get("status"), str) \
                or entry.get("status") not in {"covered", "partial", "gap", "not_applicable", "not applicable"}:
            raise AuditError(f"Invalid {package} crosswalk test/hash/status: {path}")
        if path in result:
            raise AuditError(f"Duplicate {package} crosswalk test: {path}")
        result[path] = {"sha256": digest, "coverage": entry}
    return result


def _audit_package(package, repo, baseline, crosswalk_path):
    commit = baseline["commit"]
    # A shallow checkout without the pinned tree must fail, never look unchanged.
    _git(repo, "cat-file", "-e", commit + "^{commit}")
    head = _git(repo, "rev-parse", "--verify", "HEAD^{commit}").decode("ascii").strip()
    if _git(repo, "status", "--porcelain", "--untracked-files=no"):
        raise AuditError(f"{package} checkout has tracked modifications; commit or remove them before auditing")
    if _git(repo, "ls-files", "--others", "--exclude-standard", "--", "tests/"):
        raise AuditError(f"{package} checkout has untracked tests; commit or remove them before auditing")
    old_tree, new_tree = _tree(repo, commit), _tree(repo, head)
    cache = {}

    def content(entry):
        if entry is None or entry["kind"] != "blob":
            return None
        oid = entry["oid"]
        if oid not in cache:
            cache[oid] = _git(repo, "cat-file", "blob", oid)
        return cache[oid]

    for path, digest in baseline["files_sha256"].items():
        original = content(old_tree.get(path))
        if original is None or _hash(original) != digest:
            raise AuditError(f"{package} baseline source hash mismatch: {path}")
    coverage = _crosswalk(crosswalk_path, package, commit)
    baseline_tests = {path for path in old_tree if _direct_test(path)}
    if set(coverage) != baseline_tests:
        missing = sorted(baseline_tests - set(coverage))
        extra = sorted(set(coverage) - baseline_tests)
        raise AuditError(f"{package} baseline crosswalk file support differs: missing={missing}, extra={extra}")
    for path, entry in coverage.items():
        original = content(old_tree[path])
        if original is None or _hash(original) != entry["sha256"]:
            raise AuditError(f"{package} baseline crosswalk hash mismatch: {path}")

    test_changes, source_changes = [], []
    for path in sorted(set(old_tree) | set(new_tree)):
        old, new = old_tree.get(path), new_tree.get(path)
        if old == new:
            continue
        old_content, new_content = content(old), content(new)
        change = {
            "path": path,
            "status": "added" if old is None else "removed" if new is None else "modified",
            "old": None if old is None else {**old, "sha256": None if old_content is None else _hash(old_content)},
            "new": None if new is None else {**new, "sha256": None if new_content is None else _hash(new_content)},
        }
        if path.startswith("tests/"):
            change.update({
                "is_testthat_file": _direct_test(path),
                "old_test_that": _declarations(path, old_content),
                "new_test_that": _declarations(path, new_content),
                "previous_coverage": coverage[path]["coverage"] if path in coverage else None,
                "coverage_review": "required",
            })
            test_changes.append(change)
        else:
            source_changes.append(change)
    status = "manual_review_required" if test_changes or source_changes else "unchanged"
    return {
        "baseline_commit": commit, "current_commit": head, "status": status,
        "baseline_crosswalk_test_count": len(coverage),
        "test_changes": test_changes, "source_changes": source_changes,
        "test_change_count": len(test_changes), "source_change_count": len(source_changes),
    }


def audit(did, drdid, baseline=DEFAULT_BASELINE, did_crosswalk=DEFAULT_CROSSWALKS["did"],
          drdid_crosswalk=DEFAULT_CROSSWALKS["DRDID"]):
    """Return deterministic evidence for two clean, already-fetched checkouts."""
    sources = _baseline(baseline)
    packages = {}
    for package, repo, crosswalk in (("did", did, did_crosswalk), ("DRDID", drdid, drdid_crosswalk)):
        packages[package] = _audit_package(package, repo, sources[package], crosswalk)
    test_drift = any(item["test_changes"] for item in packages.values())
    source_drift = any(item["source_changes"] for item in packages.values())
    return {
        "schema_version": 1,
        "status": "manual_review_required" if test_drift or source_drift else "unchanged",
        "exit_code": 3 if test_drift or source_drift else 0,
        "comparison": "pinned baseline Git tree versus checkout HEAD; no network access",
        "review_drift_scope": "all tracked Git-tree paths, including runtime sources, tests, documentation and build metadata",
        "test_drift_scope": "all tracked paths under tests/, including helpers and fixtures",
        "coverage_limit": "Declarations and prior crosswalk entries are evidence for manual review; automation does not port tests or establish assertion coverage.",
        "packages": packages,
    }


def _markdown(report):
    lines = ["# R upstream change audit", "", f"Status: **{report['status']}** (exit {report['exit_code']}).", ""]
    if report["status"] == "structural_failure":
        return "\n".join(lines + [report["error"], ""])
    lines.extend([report["coverage_limit"], "", "Review scope: " + report["review_drift_scope"] + ".",
                  "Test inventory: " + report["test_drift_scope"] + ".", ""])
    for package, result in report["packages"].items():
        lines.extend([
            f"## {package}", "",
            f"Baseline `{result['baseline_commit']}` → current `{result['current_commit']}`.", "",
            f"{result['test_change_count']} test-tree changes; {result['source_change_count']} other tracked-file changes.", "",
        ])
        if result["test_changes"]:
            lines.extend(["| Test-tree file | Change | Prior coverage |", "| --- | --- | --- |"])
            for change in result["test_changes"]:
                coverage = change["previous_coverage"]
                status = coverage["status"] if coverage else "no baseline test-file entry"
                path = change["path"].replace("|", "\\|").replace("\n", " ")
                lines.append(f"| `{path}` | {change['status']} | {status} |")
            lines.extend(["", "Every listed test-tree change requires manual coverage review; prior coverage is not a new-tree verdict.", ""])
        if result["source_changes"]:
            lines.extend(["Other tracked-file changes (manual source review required):", ""])
            for change in result["source_changes"]:
                lines.append(f"- {change['status']}: `{change['path']}`")
            lines.append("")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--did", type=Path, required=True)
    parser.add_argument("--drdid", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--did-crosswalk", type=Path, default=DEFAULT_CROSSWALKS["did"])
    parser.add_argument("--drdid-crosswalk", type=Path, default=DEFAULT_CROSSWALKS["DRDID"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        report = audit(args.did, args.drdid, args.baseline, args.did_crosswalk, args.drdid_crosswalk)
    except AuditError as exc:
        report = {"schema_version": 1, "status": "structural_failure", "exit_code": 2, "error": str(exc)}
    try:
        for destination in (args.output, args.summary):
            destination.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        args.summary.write_text(_markdown(report), encoding="utf-8")
    except OSError as exc:
        print(f"Cannot write audit artifacts: {exc}", file=sys.stderr)
        return 2
    print(f"R upstream audit: {report['status']} (exit {report['exit_code']})")
    return report["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())
