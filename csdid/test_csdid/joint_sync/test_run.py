"""Failure probes for the installed-wheel gate, independent of estimator code."""
import importlib.util
import json
from pathlib import Path
import types

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
