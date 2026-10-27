"""Failure probes for the installed-wheel gate, independent of estimator code."""
import importlib.util
import json
from pathlib import Path
import types
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest


SPEC = importlib.util.spec_from_file_location("joint_wheel_gate_under_test", Path(__file__).with_name("run.py"))
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


def test_peer_revision_mismatch_fails_before_any_build_or_package_install(tmp_path, monkeypatch):
    csdid, peer = tmp_path / "csdid-source", tmp_path / "peer-source"
    csdid.mkdir()
    peer.mkdir()
    lock = tmp_path / "reviewed.json"
    lock.write_text(json.dumps(dict(schema_version=1, repository="d2cml-ai/DRDIDpy",
                                   commit="a" * 40, version="1.1.6")))
    def forbidden(*args, **kwargs):
        raise AssertionError("A rejected peer revision must not execute a child process")
    monkeypatch.setattr(gate.subprocess, "run", forbidden)
    monkeypatch.setattr(gate.subprocess, "check_output", forbidden)
    output = tmp_path / "result"
    assert gate.run(csdid, peer, output, lock, drdid_ref="b" * 40) == 2
    receipt = json.loads((output / "receipt.json").read_text())
    assert receipt["status"] == "failed"
    assert receipt["error"]["stage"] == "peer-lock"
    assert "differs from reviewed" in receipt["error"]["message"]
    assert receipt["commands"] == []
    assert not (output / "venv").exists()


def _installed_modules(tmp_path):
    site = tmp_path / "venv/lib/site-packages"
    expected, modules = {}, {}
    for package in ("csdid", "drdid"):
        filename = site / package / "__init__.py"
        filename.parent.mkdir(parents=True)
        filename.write_text("# exact reviewed candidate\n")
        expected[str(filename.relative_to(site))] = gate.sha256(filename)
        modules[package] = types.SimpleNamespace(__file__=str(filename))
    return site, expected, modules


def test_origin_guard_rejects_a_checkout_even_with_the_correct_package_name(tmp_path):
    site, expected, modules = _installed_modules(tmp_path)
    source = tmp_path / "checkout/csdid/att_gt.py"
    source.parent.mkdir(parents=True)
    source.write_text("# source checkout must not supply this test\n")
    modules["csdid.att_gt"] = types.SimpleNamespace(__file__=str(source))
    with pytest.raises(ValueError, match="outside installed wheels"):
        gate.collect_and_validate_origins(site, expected, modules)


def test_origin_guard_rejects_a_modified_same_version_wheel(tmp_path):
    site, expected, modules = _installed_modules(tmp_path)
    assert set(gate.collect_and_validate_origins(site, expected, modules)) == {"csdid", "drdid"}
    (site / "drdid/__init__.py").write_text("# a stale or modified same-version package\n")
    with pytest.raises(ValueError, match="differs from reviewed wheel/source"):
        gate.collect_and_validate_origins(site, expected, modules)


def test_inherited_pytest_selector_cannot_turn_a_real_failure_into_a_partial_pass(tmp_path):
    # Without the cleaned environment, '-k keep' executes only the passing test;
    # a nonexistent external plugin would also prevent canonical collection.
    (tmp_path / "test_keep.py").write_text("def test_keep():\n    assert True\n")
    (tmp_path / "test_drop.py").write_text("def test_must_also_run():\n    assert False, 'must not be silently deselected'\n")
    config = tmp_path / "pytest.ini"
    config.write_text("[pytest]\naddopts =\n")
    injected = dict(os.environ, PYTEST_ADDOPTS="-k keep", PYTEST_PLUGINS="missing_external_plugin",
                    CSDID_R_REFERENCE_DIR="/unreviewed/csdid/fixtures",
                    CSDID_R_REFERENCE_SOURCES="/unreviewed/lock.json",
                    DRDID_R_REFERENCE="/unreviewed/drdid.json", DRDID_R_REVISION="unreviewed")
    clean = gate.clean_environment(injected)
    for name in ("CSDID_R_REFERENCE_DIR", "CSDID_R_REFERENCE_SOURCES",
                 "DRDID_R_REFERENCE", "DRDID_R_REVISION"):
        assert name not in clean
    junit = tmp_path / "complete.xml"
    process = subprocess.run([sys.executable, "-m", "pytest", "-c", str(config),
                              str(tmp_path), "-q", "--tb=no", "--junitxml=" + str(junit)],
                             env=clean, cwd=str(tmp_path), capture_output=True, text=True)
    assert process.returncode == 1, process.stdout + process.stderr
    cases = list(ET.parse(str(junit)).getroot().iter("testcase"))
    assert len(cases) == 2
    assert sum(case.find("failure") is not None for case in cases) == 1


def test_dirty_working_lock_cannot_override_a_committed_release_pair(tmp_path):
    source = tmp_path / "archive"
    path = source / "csdid/test_csdid/joint_sync/python_drdid.json"
    path.parent.mkdir(parents=True)
    committed = dict(schema_version=1, repository="d2cml-ai/DRDIDpy",
                     commit="a" * 40, version="1.1.7")
    path.write_text(json.dumps(committed))
    dirty = dict(committed, commit="b" * 40)
    with pytest.raises(ValueError, match="differs from the committed"):
        gate.bind_archived_peer_lock(source, dirty, "HEAD")
    assert gate.bind_archived_peer_lock(source, committed, "HEAD")["mode"] == "committed_peer_lock"
    path.unlink()
    with pytest.raises(ValueError, match="HEAD must contain"):
        gate.bind_archived_peer_lock(source, committed, "HEAD")
    assert "historical_explicit_ref" in gate.bind_archived_peer_lock(
        source, committed, "c" * 40)["mode"]


def test_empty_junit_cannot_be_reported_as_a_successful_joint_suite(tmp_path):
    empty = tmp_path / "empty.xml"
    empty.write_text('<testsuites><testsuite tests="0" failures="0" errors="0" skipped="0"/></testsuites>')
    with pytest.raises(ValueError, match="Empty JUnit"):
        gate._junit_counts(empty)
