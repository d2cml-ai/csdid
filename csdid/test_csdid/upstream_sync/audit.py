"""Compare local R Git trees with the committed numerical/coverage baseline.

This stdlib-only command never fetches sources or translates tests. Every
tracked Git-tree change requires human review, including runtime sources,
tests, helpers, fixtures, documentation and build metadata. Test and other
source changes remain separate evidence categories. Exact reviewed automation
additions and a bounded DESCRIPTION version bump are retained in the report
without requiring a numerical rebaseline.
Exit codes: 0 = no unreviewed drift, 3 = review required, 2 = invalid evidence.
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
DEFAULT_REVIEWED_AUTOMATION = Path(__file__).with_name("reviewed_automation.json")
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


def _reviewed_automation(path, sources):
    """Read a narrow, baseline-bound review record; never ignore a directory."""
    record = _read_json(path)
    if set(record) != {"schema_version", "packages"} \
            or not isinstance(record["packages"], dict) \
            or set(record["packages"]) != set(sources):
        raise AuditError("Reviewed automation must identify exactly did and DRDID")
    reviewed, versions = {}, {}
    for package, entry in record["packages"].items():
        if not isinstance(entry, dict) \
                or set(entry) - {"version_metadata"} != {"baseline_commit", "rationale", "additions"} \
                or entry["baseline_commit"] != sources[package]["commit"]:
            raise AuditError(f"{package} reviewed automation baseline does not match numerical baseline")
        if not isinstance(entry["rationale"], str) or not entry["rationale"].strip() \
                or not isinstance(entry["additions"], list):
            raise AuditError(f"Invalid {package} reviewed automation rationale/additions")
        additions = {}
        for addition in entry["additions"]:
            if not isinstance(addition, dict) \
                    or set(addition) != {"path", "status", "kind", "mode", "sha256"}:
                raise AuditError(f"Invalid {package} reviewed automation addition")
            name = addition["path"]
            if not _valid_path(name) or "\\" in name \
                    or PurePosixPath(name).parent != PurePosixPath(".github/workflows") \
                    or PurePosixPath(name).suffix not in {".yml", ".yaml"} \
                    or addition["status"] != "added" \
                    or addition["kind"] != "blob" or addition["mode"] != "100644" \
                    or not isinstance(addition["sha256"], str) \
                    or not _SHA256.fullmatch(addition["sha256"]):
                raise AuditError(f"Invalid {package} reviewed automation path/status/blob/mode/hash")
            if name in additions:
                raise AuditError(f"Duplicate {package} reviewed automation path: {name}")
            additions[name] = addition
        reviewed[package] = additions
        metadata = entry.get("version_metadata", [])
        if not isinstance(metadata, list) or len(metadata) > 1:
            raise AuditError(f"Invalid {package} reviewed version metadata records")
        versions[package] = metadata[0] if metadata else None
        if not metadata:
            continue
        version = metadata[0]
        if not isinstance(version, dict) or set(version) != {
                "path", "status", "kind", "mode", "before_sha256", "after_sha256",
                "before_version", "after_version", "rationale"}:
            raise AuditError(f"Invalid {package} reviewed version metadata record")
        before, after = version["before_version"], version["after_version"]
        if package != "did" or version["path"] != "DESCRIPTION" \
                or version["status"] != "modified" or version["kind"] != "blob" \
                or version["mode"] != "100644" \
                or not isinstance(version["rationale"], str) or not version["rationale"].strip() \
                or any(not isinstance(version[key], str) or not _SHA256.fullmatch(version[key])
                       for key in ("before_sha256", "after_sha256")) \
                or any(not isinstance(value, str) or not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+", value)
                       for value in (before, after)):
            raise AuditError(f"Invalid {package} reviewed DESCRIPTION identity/hash/version")
        old_parts, new_parts = before.split("."), after.split(".")
        if old_parts[:3] != new_parts[:3] or int(new_parts[3]) != int(old_parts[3]) + 1:
            raise AuditError(f"{package} reviewed DESCRIPTION must increment only the fourth version component")
    return reviewed, versions


def _description_version(content):
    """Locate exactly one complete Version field without normalizing bytes."""
    if content is None:
        return None
    fields, offset = [], 0
    for line in content.splitlines(keepends=True):
        if line.startswith(b"Version:"):
            fields.append((offset, line))
        offset += len(line)
    if len(fields) != 1:
        return None
    offset, line = fields[0]
    match = re.fullmatch(rb"Version:[ \t]*([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)[ \t]*(?:\r?\n)?", line)
    if match is None:
        return None
    return (match.group(1).decode("ascii"), offset + match.start(1), offset + match.end(1))


def _reviewed_description_content(package, record, tree, content):
    """Validate both recorded hashes against a version-only baseline edit."""
    if record is None:
        return None
    old = tree.get("DESCRIPTION")
    original = content(old)
    field = _description_version(original)
    if old is None or old["kind"] != record["kind"] or old["mode"] != record["mode"] \
            or original is None or _hash(original) != record["before_sha256"] \
            or field is None or field[0] != record["before_version"]:
        raise AuditError(f"{package} reviewed DESCRIPTION baseline content/mode/version mismatch")
    updated = original[:field[1]] + record["after_version"].encode("ascii") + original[field[2]:]
    if _hash(updated) != record["after_sha256"]:
        raise AuditError(f"{package} reviewed DESCRIPTION after hash is not the exact version-only edit")
    return updated


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


def _audit_package(package, repo, baseline, crosswalk_path, reviewed_automation, version_metadata):
    commit = baseline["commit"]
    # A shallow checkout without the pinned tree must fail, never look unchanged.
    _git(repo, "cat-file", "-e", commit + "^{commit}")
    head = _git(repo, "rev-parse", "--verify", "HEAD^{commit}").decode("ascii").strip()
    if _git(repo, "status", "--porcelain", "--untracked-files=no"):
        raise AuditError(f"{package} checkout has tracked modifications; commit or remove them before auditing")
    if _git(repo, "ls-files", "--others", "--exclude-standard", "--", "tests/"):
        raise AuditError(f"{package} checkout has untracked tests; commit or remove them before auditing")
    old_tree, new_tree = _tree(repo, commit), _tree(repo, head)
    if any(path in old_tree for path in reviewed_automation):
        raise AuditError(f"{package} reviewed automation must be absent from the baseline tree")
    cache = {}

    def content(entry):
        if entry is None or entry["kind"] != "blob":
            return None
        oid = entry["oid"]
        if oid not in cache:
            cache[oid] = _git(repo, "cat-file", "blob", oid)
        return cache[oid]

    reviewed_description = _reviewed_description_content(package, version_metadata, old_tree, content)
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
            reviewed = reviewed_automation.get(path)
            accepted = (reviewed is not None and old is None and new is not None
                        and all(change[key] == reviewed[key] for key in ("path", "status"))
                        and all(change["new"][key] == reviewed[key]
                                for key in ("kind", "mode", "sha256")))
            change["source_review"] = "reviewed_automation" if accepted else "required"
            if (version_metadata is not None and path == "DESCRIPTION"
                    and change["status"] == version_metadata["status"]
                    and old is not None and new is not None
                    and all(old[key] == new[key] == version_metadata[key] for key in ("kind", "mode"))
                    and change["old"]["sha256"] == version_metadata["before_sha256"]
                    and change["new"]["sha256"] == version_metadata["after_sha256"]
                    and _description_version(new_content) is not None
                    and _description_version(new_content)[0] == version_metadata["after_version"]
                    and new_content == reviewed_description):
                change["source_review"] = "reviewed_version_metadata"
            source_changes.append(change)
    reviewed_paths = [item["path"] for item in source_changes
                      if item["source_review"] == "reviewed_automation"]
    unreviewed_paths = [item["path"] for item in source_changes
                        if item["source_review"] == "required"]
    version_paths = [item["path"] for item in source_changes
                     if item["source_review"] == "reviewed_version_metadata"]
    review_required = bool(test_changes or unreviewed_paths)
    status = ("manual_review_required" if review_required else
              "reviewed_changes_only" if version_paths else
              "reviewed_automation_only" if reviewed_paths else "unchanged")
    return {
        "baseline_commit": commit, "current_commit": head, "status": status,
        "baseline_crosswalk_test_count": len(coverage),
        "test_changes": test_changes, "source_changes": source_changes,
        "test_change_count": len(test_changes), "source_change_count": len(source_changes),
        "source_review_required": review_required,
        "reviewed_automation_paths": reviewed_paths,
        "reviewed_automation_change_count": len(reviewed_paths),
        "reviewed_version_metadata_paths": version_paths,
        "reviewed_version_metadata_change_count": len(version_paths),
        "unreviewed_source_paths": unreviewed_paths,
        "unreviewed_source_change_count": len(unreviewed_paths),
    }


def audit(did, drdid, baseline=DEFAULT_BASELINE, did_crosswalk=DEFAULT_CROSSWALKS["did"],
          drdid_crosswalk=DEFAULT_CROSSWALKS["DRDID"],
          reviewed_automation=DEFAULT_REVIEWED_AUTOMATION):
    """Return deterministic evidence for two clean, already-fetched checkouts."""
    sources = _baseline(baseline)
    reviewed, versions = _reviewed_automation(reviewed_automation, sources)
    packages = {}
    for package, repo, crosswalk in (("did", did, did_crosswalk), ("DRDID", drdid, drdid_crosswalk)):
        packages[package] = _audit_package(package, repo, sources[package], crosswalk,
                                           reviewed[package], versions[package])
    review_required = any(item["source_review_required"] for item in packages.values())
    accepted = any(item["reviewed_automation_paths"] for item in packages.values())
    accepted_metadata = any(item["reviewed_version_metadata_paths"] for item in packages.values())
    return {
        "schema_version": 1,
        "status": ("manual_review_required" if review_required else
                   "reviewed_changes_only" if accepted_metadata else
                   "reviewed_automation_only" if accepted else "unchanged"),
        "exit_code": 3 if review_required else 0,
        "source_review_required": review_required,
        "reviewed_automation_record_sha256": _hash(Path(reviewed_automation).read_bytes()),
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
            label = ("manual source review required" if result["unreviewed_source_paths"]
                     else "reviewed changes only" if result["reviewed_version_metadata_paths"]
                     else "reviewed automation only")
            lines.extend([f"Other tracked-file changes ({label}):", ""])
            for change in result["source_changes"]:
                review = {"reviewed_automation": "reviewed automation",
                          "reviewed_version_metadata": "reviewed version metadata",
                          "required": "manual review required"}[change["source_review"]]
                lines.append(f"- {change['status']}: `{change['path']}` ({review})")
            lines.append(f"\n{result['reviewed_automation_change_count']} reviewed automation additions; "
                         f"{result['unreviewed_source_change_count']} unreviewed other tracked-file changes.")
            lines.append(f"{result['reviewed_version_metadata_change_count']} reviewed DESCRIPTION version changes.")
            lines.append("")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--did", type=Path, required=True)
    parser.add_argument("--drdid", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--did-crosswalk", type=Path, default=DEFAULT_CROSSWALKS["did"])
    parser.add_argument("--drdid-crosswalk", type=Path, default=DEFAULT_CROSSWALKS["DRDID"])
    parser.add_argument("--reviewed-automation", type=Path, default=DEFAULT_REVIEWED_AUTOMATION)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        report = audit(args.did, args.drdid, args.baseline, args.did_crosswalk,
                       args.drdid_crosswalk, args.reviewed_automation)
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
