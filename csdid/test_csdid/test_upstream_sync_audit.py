"""Exercise the offline drift gate against real, temporary Git histories."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


AUDIT_PATH = Path(__file__).parent / "upstream_sync" / "audit.py"
_SPEC = importlib.util.spec_from_file_location("upstream_sync_audit", AUDIT_PATH)
audit_module = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(audit_module)


def _git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.decode().strip()


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _json(path, value):
    _write(path, json.dumps(value))


def _sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def _commit(repo):
    _git(repo, "add", "--all")
    _git(repo, "commit", "-m", "test checkpoint")
    return _git(repo, "rev-parse", "HEAD")


@pytest.fixture
def histories(tmp_path):
    repos, crosswalks, sources = {}, {}, {}
    texts = {
        "tests/testthat/test-existing.R": 'test_that("existing assertion", {\n  expect_equal(1, 1)\n})\n',
        "tests/testthat/test-removed.R": 'testthat::test_that("removed assertion", { expect_true(TRUE) })\n',
    }
    for package in ("did", "DRDID"):
        repo = tmp_path / package
        repo.mkdir()
        _git(repo, "init")
        _git(repo, "config", "user.name", "Audit Test")
        _git(repo, "config", "user.email", "audit-test@example.invalid")
        _git(repo, "config", "commit.gpgsign", "false")
        for path, text in texts.items():
            _write(repo / path, text)
        _write(repo / "tests/testthat/helper-fixture.R", "make_fixture <- function() 1\n")
        _write(repo / "tests/testthat/_snaps/output.md", "baseline snapshot\n")
        _write(repo / "DESCRIPTION", f"Package: {package}\nVersion: 1.0.0\n")
        _write(repo / "R/estimator.R", "estimate <- function() 1\n")
        commit = _commit(repo)
        sources[package] = {"commit": commit, "files_sha256": {
            "R/estimator.R": _sha("estimate <- function() 1\n"),
        }}
        is_did = package == "did"
        crosswalk = {
            "schema_version": 1,
            "r_commit" if is_did else "r_revision": commit,
            "entries" if is_did else "files": [
                {"r_file": path, "status": "partial", "python_tests": ["test_original_intent"],
                 "remaining_gaps": ["another assertion still requires review"],
                 "r_file_sha256" if is_did else "r_sha256": _sha(text)}
                for path, text in texts.items()
            ],
        }
        crosswalk_path = tmp_path / f"{package}-crosswalk.json"
        _json(crosswalk_path, crosswalk)
        repos[package] = repo
        crosswalks[package] = crosswalk_path
    manifest = tmp_path / "manifest.json"
    _json(manifest, {"schema_version": 1, "sources": sources})
    reviewed = tmp_path / "reviewed-automation.json"
    _json(reviewed, {"schema_version": 1, "packages": {
        package: {"baseline_commit": source["commit"], "rationale": "Temporary test review", "additions": []}
        for package, source in sources.items()}})
    return {"repos": repos, "crosswalks": crosswalks, "manifest": manifest,
            "sources": sources, "texts": texts, "root": tmp_path, "reviewed": reviewed}


def _audit(history):
    return audit_module.audit(history["repos"]["did"], history["repos"]["DRDID"],
                              history["manifest"], history["crosswalks"]["did"],
                              history["crosswalks"]["DRDID"], history["reviewed"])


def _args(history):
    return ["--did", str(history["repos"]["did"]), "--drdid", str(history["repos"]["DRDID"]),
            "--baseline", str(history["manifest"]),
            "--did-crosswalk", str(history["crosswalks"]["did"]),
            "--drdid-crosswalk", str(history["crosswalks"]["DRDID"]),
            "--reviewed-automation", str(history["reviewed"]),
            "--output", str(history["root"] / "report.json"),
            "--summary", str(history["root"] / "summary.md")]


def test_unchanged_histories_and_cli_are_deterministic(histories):
    report = _audit(histories)
    assert report["status"] == "unchanged"
    assert report["exit_code"] == 0
    for package, result in report["packages"].items():
        assert result["baseline_commit"] == result["current_commit"] == histories["sources"][package]["commit"]
        assert result["baseline_crosswalk_test_count"] == 2
        assert result["test_changes"] == result["source_changes"] == []
    # Exercise the standalone entry point from an unrelated working directory.
    command = [sys.executable, str(AUDIT_PATH.resolve()), *_args(histories)]
    first = subprocess.run(command, cwd=histories["root"], capture_output=True, text=True)
    assert first.returncode == 0, first.stderr
    json_bytes = (histories["root"] / "report.json").read_bytes()
    markdown_bytes = (histories["root"] / "summary.md").read_bytes()
    assert json.loads(json_bytes) == report
    assert audit_module.main(_args(histories)) == 0
    assert (histories["root"] / "report.json").read_bytes() == json_bytes
    assert (histories["root"] / "summary.md").read_bytes() == markdown_bytes
    assert audit_module.DEFAULT_BASELINE == AUDIT_PATH.resolve().parents[1] / "r_ref/current/manifest.json"


@pytest.mark.parametrize("package", ["did", "DRDID"])
def test_source_only_changes_require_review_without_claiming_test_drift(histories, package):
    repo = histories["repos"][package]
    updated = "estimate <- function() 2\n"
    _write(repo / "R/estimator.R", updated)
    _write(repo / "NEWS.md", "A code change.\n")
    current = _commit(repo)
    report = _audit(histories)
    assert report["status"] == "manual_review_required"
    assert report["exit_code"] == 3
    result = report["packages"][package]
    assert result["status"] == "manual_review_required"
    assert result["current_commit"] == current
    assert result["test_changes"] == []
    assert [(item["path"], item["status"]) for item in result["source_changes"]] == [
        ("NEWS.md", "added"), ("R/estimator.R", "modified"),
    ]
    change = result["source_changes"][1]
    assert change["old"]["sha256"] == histories["sources"][package]["files_sha256"]["R/estimator.R"]
    assert change["new"]["sha256"] == _sha(updated)
    assert audit_module.main(_args(histories)) == 3
    summary = (histories["root"] / "summary.md").read_text()
    assert "manual source review required" in summary
    assert "0 test-tree changes; 2 other tracked-file changes" in summary


@pytest.mark.parametrize("path", ["src/estimator.cpp", "DESCRIPTION", "README.md",
                                  ".github/workflows/check.yml"])
def test_every_tracked_tree_category_requires_review(histories, path):
    repo = histories["repos"]["did"]
    _write(repo / path, "changed tracked source or metadata\n")
    _commit(repo)
    report = _audit(histories)
    assert report["exit_code"] == 3
    result = report["packages"]["did"]
    assert result["test_changes"] == []
    assert [change["path"] for change in result["source_changes"]] == [path]


def _review_workflow(history, package, text="name: reviewed automation\n"):
    path = ".github/workflows/python-downstream.yml"
    record = json.loads(history["reviewed"].read_text())
    record["packages"][package]["additions"] = [{
        "path": path, "status": "added", "kind": "blob", "mode": "100644", "sha256": _sha(text),
    }]
    _json(history["reviewed"], record)
    _write(history["repos"][package] / path, text)
    return path


@pytest.mark.parametrize("package", ["did", "DRDID"])
def test_exact_reviewed_workflow_addition_retains_raw_drift_and_succeeds(histories, package):
    path = _review_workflow(histories, package)
    _commit(histories["repos"][package])
    report = _audit(histories)
    assert report["status"] == "reviewed_automation_only"
    assert report["exit_code"] == 0
    assert report["source_review_required"] is False
    assert report["reviewed_automation_record_sha256"] == hashlib.sha256(histories["reviewed"].read_bytes()).hexdigest()
    result = report["packages"][package]
    assert result["source_change_count"] == result["reviewed_automation_change_count"] == 1
    assert result["unreviewed_source_change_count"] == 0
    assert result["reviewed_automation_paths"] == [path]
    assert result["test_changes"] == []
    change, = result["source_changes"]
    assert change["path"] == path and change["status"] == "added" and change["old"] is None
    assert change["new"]["mode"] == "100644" and change["new"]["kind"] == "blob"
    assert change["source_review"] == "reviewed_automation"
    assert audit_module.main(_args(histories)) == 0
    summary = (histories["root"] / "summary.md").read_text()
    assert "1 reviewed automation additions; 0 unreviewed other tracked-file changes" in summary
    assert path in summary and "reviewed_automation_only" in summary


@pytest.mark.parametrize("tamper", ["content", "executable", "symlink", "other_source", "test"])
def test_reviewed_automation_cannot_hide_modified_bytes_modes_or_other_drift(histories, tamper):
    path = _review_workflow(histories, "did")
    repo = histories["repos"]["did"]
    if tamper == "content":
        _write(repo / path, "name: different automation\n")
    elif tamper == "executable":
        (repo / path).chmod(0o755)
    elif tamper == "symlink":
        (repo / path).unlink()
        (repo / path).symlink_to("../../R/estimator.R")
    else:
        _write(repo / ("NEWS.md" if tamper == "other_source" else "tests/fixtures/extra.csv"), "unreviewed\n")
    _commit(repo)
    report = _audit(histories)
    assert report["status"] == "manual_review_required" and report["exit_code"] == 3
    assert report["source_review_required"] is True
    result = report["packages"]["did"]
    assert result["source_review_required"] is True
    if tamper in {"content", "executable", "symlink"}:
        assert result["reviewed_automation_paths"] == []
        assert result["unreviewed_source_paths"] == [path]
    else:
        assert result["reviewed_automation_paths"] == [path]
        assert result["test_change_count"] == (1 if tamper == "test" else 0)
        assert result["unreviewed_source_paths"] == ([] if tamper == "test" else ["NEWS.md"])


@pytest.mark.parametrize("field,value", [
    ("path", "R/estimator.R"), ("path", "../workflow.yml"),
    ("path", ".github/workflows/nested/extra.yml"), ("status", "modified"),
    ("kind", "commit"), ("mode", "100755"), ("sha256", "not-a-hash"),
])
def test_invalid_review_record_fails_closed(histories, field, value):
    _review_workflow(histories, "did")
    _commit(histories["repos"]["did"])
    record = json.loads(histories["reviewed"].read_text())
    record["packages"]["did"]["additions"][0][field] = value
    _json(histories["reviewed"], record)
    assert audit_module.main(_args(histories)) == 2
    assert json.loads((histories["root"] / "report.json").read_text())["status"] == "structural_failure"


@pytest.mark.parametrize("defect", ["baseline", "duplicate", "missing", "schema"])
def test_review_records_are_required_unique_and_bound_to_numerical_baseline(histories, defect):
    _review_workflow(histories, "did")
    _commit(histories["repos"]["did"])
    record = json.loads(histories["reviewed"].read_text())
    if defect == "baseline":
        record["packages"]["did"]["baseline_commit"] = "a" * 40
    elif defect == "duplicate":
        record["packages"]["did"]["additions"] *= 2
    elif defect == "schema":
        record["ignore_directory"] = ".github"
    if defect == "missing":
        histories["reviewed"].unlink()
    else:
        _json(histories["reviewed"], record)
    assert audit_module.main(_args(histories)) == 2


def test_review_record_cannot_approve_a_preexisting_workflow(histories):
    path = _review_workflow(histories, "did")
    repo = histories["repos"]["did"]
    baseline = _commit(repo)
    manifest = json.loads(histories["manifest"].read_text())
    manifest["sources"]["did"]["commit"] = baseline
    _json(histories["manifest"], manifest)
    crosswalk = json.loads(histories["crosswalks"]["did"].read_text())
    crosswalk["r_commit"] = baseline
    _json(histories["crosswalks"]["did"], crosswalk)
    record = json.loads(histories["reviewed"].read_text())
    record["packages"]["did"]["baseline_commit"] = baseline
    _json(histories["reviewed"], record)
    _write(repo / path, "name: changed existing workflow\n")
    _commit(repo)
    with pytest.raises(audit_module.AuditError, match="absent from the baseline"):
        _audit(histories)


@pytest.mark.parametrize("reviewed,comparison_code,expected_status,expected_exit", [
    (False, 0, "coverage_review_required", 3),
    (False, 1, "numerical_mismatch", 1),
    (True, 0, "passed", 0),
    (True, 1, "numerical_mismatch", 1),
])
def test_real_git_audit_keeps_runner_numerical_stages_and_failure_precedence(
        histories, monkeypatch, reviewed, comparison_code, expected_status, expected_exit):
    # The Git audit is real; R generation and estimator execution are stubbed
    # here. The separate live source-only probe checks actual R/Python numbers.
    spec = importlib.util.spec_from_file_location(
        "source_only_runner", AUDIT_PATH.with_name("run.py"))
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    repo = histories["repos"]["did"]
    if reviewed:
        changed_path = _review_workflow(histories, "did")
    else:
        changed_path = "R/estimator.R"
        _write(repo / changed_path, "estimate <- function() 2\n")
    _commit(repo)
    pinned = histories["root"] / "reviewed"
    pinned.mkdir()
    for name in runner.INPUTS + ("generate.R", "manifest.json"):
        (pinned / name).write_text("reviewed " + name)
    monkeypatch.setattr(runner, "PINNED", pinned)
    before = runner._hashes(pinned)
    output = histories["root"] / "runner-output"
    numerical_stages = []

    def execute(command, log, environment, **kwargs):
        if log.stem == "audit":
            args = _args(histories)
            args[args.index("--output") + 1] = str(output / "audit.json")
            args[args.index("--summary") + 1] = str(output / "audit.md")
            with log.open("w") as stream:
                return subprocess.run([sys.executable, str(AUDIT_PATH), *args],
                                      env=environment, stdout=stream,
                                      stderr=subprocess.STDOUT, check=False).returncode
        numerical_stages.append(log.stem)
        log.write_text("stubbed numerical subprocess\n")
        return comparison_code if log.stem == "pytest" else 0

    monkeypatch.setattr(runner, "_command", execute)
    assert runner.run(repo, histories["repos"]["DRDID"], output) == expected_exit
    report = json.loads((output / "report.json").read_text())
    audit = json.loads((output / "audit.json").read_text())
    audit_code = 0 if reviewed else 3
    assert audit["exit_code"] == audit_code
    assert audit["status"] == ("reviewed_automation_only" if reviewed else "manual_review_required")
    assert audit["packages"]["did"]["test_changes"] == []
    assert [change["path"] for change in audit["packages"]["did"]["source_changes"]] == [changed_path]
    assert numerical_stages == ["r-generation", "pytest"]
    assert report["stages"] == {"test_coverage_audit": audit_code, "r_generation": 0,
                                "numerical_comparison": comparison_code}
    assert report["status"] == expected_status
    assert runner._hashes(pinned) == before


@pytest.mark.parametrize("package", ["did", "DRDID"])
def test_added_modified_removed_tests_require_review_with_exact_evidence(histories, package):
    repo = histories["repos"][package]
    added = '# setup\n\ntest_that("new assertion", { expect_true(TRUE) })\n'
    modified = '# changed setup\n  testthat::test_that("changed assertion", { expect_equal(2, 2) })\n'
    _write(repo / "tests/testthat/test-added.R", added)
    _write(repo / "tests/testthat/test-existing.R", modified)
    (repo / "tests/testthat/test-removed.R").unlink()
    _write(repo / "R/estimator.R", "estimate <- function() 3\n")
    _commit(repo)
    report = _audit(histories)
    assert report["status"] == "manual_review_required"
    assert report["exit_code"] == 3
    result = report["packages"][package]
    assert result["source_change_count"] == 1
    assert result["test_change_count"] == 3
    changes = result["test_changes"]
    assert [(item["path"], item["status"]) for item in changes] == [
        ("tests/testthat/test-added.R", "added"),
        ("tests/testthat/test-existing.R", "modified"),
        ("tests/testthat/test-removed.R", "removed"),
    ]
    assert all(item["coverage_review"] == "required" and item["is_testthat_file"] for item in changes)
    assert changes[0]["old"] is None
    assert changes[0]["new"]["sha256"] == _sha(added)
    assert changes[0]["previous_coverage"] is None
    assert changes[0]["new_test_that"] == [{"line": 3, "expression": added.splitlines()[2]}]
    assert changes[1]["old"]["sha256"] == _sha(histories["texts"]["tests/testthat/test-existing.R"])
    assert changes[1]["new"]["sha256"] == _sha(modified)
    assert changes[1]["old_test_that"][0]["line"] == 1
    assert changes[1]["new_test_that"] == [{"line": 2, "expression": modified.splitlines()[1].strip()}]
    assert changes[1]["previous_coverage"]["status"] == "partial"
    assert changes[1]["previous_coverage"]["remaining_gaps"] == ["another assertion still requires review"]
    assert changes[2]["new"] is None
    assert changes[2]["new_test_that"] == []
    assert changes[2]["old"]["sha256"] == _sha(histories["texts"]["tests/testthat/test-removed.R"])
    assert audit_module.main(_args(histories)) == 3
    assert "Every listed test-tree change requires manual coverage review" in (histories["root"] / "summary.md").read_text()


@pytest.mark.parametrize("path", ["tests/testthat/helper-fixture.R", "tests/testthat/_snaps/output.md",
                                  "tests/fixtures/input.csv", "tests/testthat.R"])
def test_helpers_fixtures_snapshots_and_runner_trigger_review_without_crosswalk_entry(histories, path):
    repo = histories["repos"]["did"]
    _write(repo / path, "changed test infrastructure\n")
    _commit(repo)
    report = _audit(histories)
    assert report["exit_code"] == 3
    result = report["packages"]["did"]
    assert result["source_changes"] == []
    assert len(result["test_changes"]) == 1
    change = result["test_changes"][0]
    assert change["path"] == path
    assert change["previous_coverage"] is None
    assert not change["is_testthat_file"]
    assert change["new"]["sha256"] == _sha("changed test infrastructure\n")


def test_missing_and_malformed_baselines_fail_closed_and_write_error_artifacts(histories):
    manifest = histories["manifest"]
    valid = manifest.read_text()
    malformed = [None, "{", "[]", '{"schema_version":1,"schema_version":1,"sources":{}}',
                 '{"schema_version":1,"sources":{}}',
                 valid.replace('"schema_version": 1', '"schema_version": true'),
                 valid.replace(histories["sources"]["did"]["commit"], "main"),
                 valid.replace('"files_sha256"', '"missing_hashes"')]
    for payload in malformed:
        if payload is None:
            manifest.unlink()
        else:
            manifest.write_text(payload)
        assert audit_module.main(_args(histories)) == 2
        report = json.loads((histories["root"] / "report.json").read_text())
        assert report["status"] == "structural_failure"
        assert report["exit_code"] == 2
        assert report["error"]
        assert "packages" not in report


def test_unavailable_baseline_git_object_is_failure_not_empty_diff(histories):
    manifest = json.loads(histories["manifest"].read_text())
    manifest["sources"]["did"]["commit"] = "a" * 40
    _json(histories["manifest"], manifest)
    record = json.loads(histories["reviewed"].read_text())
    record["packages"]["did"]["baseline_commit"] = "a" * 40
    _json(histories["reviewed"], record)
    with pytest.raises(audit_module.AuditError, match="cat-file.*failed"):
        _audit(histories)


def test_source_manifest_hash_mismatch_fails_even_without_changes(histories):
    manifest = json.loads(histories["manifest"].read_text())
    manifest["sources"]["did"]["files_sha256"]["R/estimator.R"] = "0" * 64
    _json(histories["manifest"], manifest)
    with pytest.raises(audit_module.AuditError, match="baseline source hash mismatch"):
        _audit(histories)


@pytest.mark.parametrize("package,entries_key,commit_key,hash_key", [
    ("did", "entries", "r_commit", "r_file_sha256"),
    ("DRDID", "files", "r_revision", "r_sha256"),
])
def test_crosswalk_is_required_pinned_complete_and_content_verified(histories, package, entries_key, commit_key, hash_key):
    path = histories["crosswalks"][package]
    valid = json.loads(path.read_text())
    path.unlink()
    with pytest.raises(audit_module.AuditError, match="Cannot read valid JSON"):
        _audit(histories)
    corruptions = [
        (lambda item: item.update({commit_key: "b" * 40}), "commit does not match"),
        (lambda item: item[entries_key][0].update({hash_key: "0" * 64}), "crosswalk hash mismatch"),
        (lambda item: item[entries_key].pop(), "file support differs"),
        (lambda item: item[entries_key].append(item[entries_key][0]), "Duplicate.*crosswalk"),
        (lambda item: item[entries_key][0].update({"r_file": "../escape.R"}), "Invalid.*crosswalk"),
        (lambda item: item[entries_key][0].update({"status": []}), "Invalid.*crosswalk"),
    ]
    for mutate, message in corruptions:
        item = json.loads(json.dumps(valid))
        mutate(item)
        _json(path, item)
        with pytest.raises(audit_module.AuditError, match=message):
            _audit(histories)
        assert audit_module.main(_args(histories)) == 2
        assert json.loads((histories["root"] / "report.json").read_text())["status"] == "structural_failure"


@pytest.mark.parametrize("path,message", [
    ("R/estimator.R", "tracked modifications"),
    ("tests/testthat/test-untracked.R", "untracked tests"),
])
def test_uncommitted_evidence_cannot_silently_disappear(histories, path, message):
    _write(histories["repos"]["did"] / path, "uncommitted change\n")
    with pytest.raises(audit_module.AuditError, match=message):
        _audit(histories)
