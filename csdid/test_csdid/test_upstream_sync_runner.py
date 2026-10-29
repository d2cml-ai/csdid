"""Failure-path tests for the non-mutating live-R comparison runner."""
import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


SCRIPT = Path(__file__).parent / "upstream_sync" / "run.py"
SPEC = importlib.util.spec_from_file_location("upstream_sync_runner", SCRIPT)
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    pinned = tmp_path / "reviewed"
    pinned.mkdir()
    for name in runner.INPUTS + ("generate.R", "manifest.json"):
        (pinned / name).write_text("reviewed " + name)
    did, drdid = tmp_path / "did", tmp_path / "DRDID"
    did.mkdir()
    drdid.mkdir()
    monkeypatch.setattr(runner, "PINNED", pinned)
    monkeypatch.setattr(runner, "source_identity", lambda path, package: {
        "commit": "a" * 40 if package == "did" else "b" * 40, "version": "1.0.0"})
    return pinned, did, drdid, tmp_path / "output"


def commands(monkeypatch, codes=None, mutate=None):
    calls = []
    codes = codes or {}

    def execute(command, log, environment, **kwargs):
        calls.append((command, environment.copy()))
        log.write_text("retained child diagnostics\n")
        stage = log.stem
        if mutate is not None:
            mutate(stage, log.parent)
        return codes.get(stage, 0)

    monkeypatch.setattr(runner, "_command", execute)
    return calls


@pytest.mark.parametrize("codes,status,exit_code,stages", [
    ({}, "passed", 0, 3),
    ({"audit": 3}, "coverage_review_required", 3, 3),
    ({"pytest": 1}, "numerical_mismatch", 1, 3),
    ({"audit": 3, "pytest": 1}, "numerical_mismatch", 1, 3),
    ({"audit": 2}, "infrastructure_error", 2, 1),
    ({"r-generation": 1}, "infrastructure_error", 2, 2),
    ({"r-generation": -1}, "infrastructure_error", 2, 2),
    ({"pytest": 5}, "infrastructure_error", 2, 3),
])
def test_run_preserves_failure_and_coverage_signals(sandbox, monkeypatch, codes, status, exit_code, stages):
    pinned, did, drdid, output = sandbox
    before = runner._hashes(pinned)
    calls = commands(monkeypatch, codes)
    monkeypatch.setenv("CSDID_R_REFERENCE_DIR", "stale-external-directory")
    monkeypatch.setenv("CSDID_R_REFERENCE_SOURCES", "stale-source-lock")
    assert runner.run(did, drdid, output) == exit_code
    report = json.loads((output / "report.json").read_text())
    assert report["status"] == status
    assert len(calls) == stages
    assert report["reviewed_references_unchanged"]
    assert runner._hashes(pinned) == before
    assert status in (output / "summary.md").read_text()
    assert "CSDID_R_REFERENCE_DIR" not in calls[0][1]
    if stages == 3:
        command, environment = calls[-1]
        assert environment["CSDID_R_REFERENCE_DIR"] == str(output / "references")
        assert environment["CSDID_R_REFERENCE_SOURCES"] == str(output / "source-lock.json")
        assert "--junitxml=" + str(output / "junit.xml") in command
        assert all(environment[key] == "1" for key in runner.THREADS)


def test_generator_cannot_replace_frozen_inputs(sandbox, monkeypatch):
    _, did, drdid, output = sandbox

    def mutate(stage, directory):
        if stage == "r-generation":
            (directory / "references" / "panel.csv").write_text("replacement")

    calls = commands(monkeypatch, mutate=mutate)
    assert runner.run(did, drdid, output) == 2
    report = json.loads((output / "report.json").read_text())
    assert "Frozen comparison input changed" in report["error"]
    assert len(calls) == 2


def test_reviewed_reference_mutation_fails_even_after_passing_comparison(sandbox, monkeypatch):
    pinned, did, drdid, output = sandbox

    def mutate(stage, directory):
        if stage == "pytest":
            (pinned / "manifest.json").write_text("unexpected mutation")

    commands(monkeypatch, mutate=mutate)
    assert runner.run(did, drdid, output) == 2
    report = json.loads((output / "report.json").read_text())
    assert report["reviewed_references_unchanged"] is False


@pytest.mark.parametrize("target", ["reviewed", "reviewed_child", "ancestor", "did_child", "drdid_child"])
def test_protected_output_locations_are_rejected_without_writes(sandbox, target):
    pinned, did, drdid, _ = sandbox
    output = {"reviewed": pinned, "reviewed_child": pinned / "results",
              "ancestor": pinned.parent, "did_child": did / "results",
              "drdid_child": drdid / "results"}[target]
    before = runner._hashes(pinned.parent)
    with pytest.raises(ValueError, match="separate"):
        runner.run(did, drdid, output)
    assert runner._hashes(pinned.parent) == before


def test_previous_artifacts_are_never_overwritten(sandbox):
    _, did, drdid, output = sandbox
    output.mkdir()
    previous = output / "previous.log"
    previous.write_text("keep evidence")
    with pytest.raises(ValueError, match="new or empty"):
        runner.run(did, drdid, output)
    assert previous.read_text() == "keep evidence"


def test_source_revision_moving_during_run_fails(sandbox, monkeypatch):
    _, did, drdid, output = sandbox
    original = runner.source_identity
    calls = []

    def identity(path, package):
        result = original(path, package)
        calls.append(package)
        if len(calls) > 2:
            result["commit"] = "c" * 40
        return result

    monkeypatch.setattr(runner, "source_identity", identity)
    commands(monkeypatch)
    assert runner.run(did, drdid, output) == 2
    assert "revision changed" in json.loads((output / "report.json").read_text())["error"]


@pytest.mark.parametrize("defect", [None, "dirty", "untracked", "package", "version"])
def test_source_identity_uses_clean_git_and_description(tmp_path, defect):
    repo = tmp_path / "source"
    repo.mkdir()
    def git(*args):
        return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()
    git("init", "-q")
    git("config", "user.email", "test@example.invalid")
    git("config", "user.name", "Test")
    name = "wrong" if defect == "package" else "did"
    (repo / "DESCRIPTION").write_text("Package: " + name + "\n" + (
        "Version: 1.0.0\n" if defect != "version" else ""))
    git("add", "DESCRIPTION")
    git("commit", "-qm", "fixture")
    if defect == "dirty":
        (repo / "DESCRIPTION").write_text("Package: did\nVersion: 2.0.0\n")
    if defect == "untracked":
        (repo / "R").mkdir()
        (repo / "R" / "extra.R").write_text("print('untracked source')\n")
    if defect is None:
        # Exercise the unpatched source_identity in a real temporary Git repo.
        assert runner.source_identity(repo, "did") == {"commit": git("rev-parse", "HEAD"), "version": "1.0.0"}
    else:
        with pytest.raises(ValueError):
            runner.source_identity(repo, "did")
