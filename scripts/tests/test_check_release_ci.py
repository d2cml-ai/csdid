"""Offline release-gate counterexamples; no package imports or uploads."""
import copy
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import HTTPError
import zipfile


SPEC = importlib.util.spec_from_file_location(
    "check_release_ci", Path(__file__).parents[1] / "check_release_ci.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)
SHA = "a" * 40
REPO = "d2cml-ai/csdid"
IDENTITY = {"full_name": REPO, "id": 123}


class ReleaseChecks(unittest.TestCase):
    def setUp(self):
        self.runs = {}
        self.checks = []
        for index, (name, (workflow, path)) in enumerate(gate.REQUIRED.items(), 1):
            run_id = index * 10
            self.runs[str(run_id)] = dict(
                name=workflow, path=path, event="push", head_branch="main", head_sha=SHA,
                status="completed", conclusion="success", repository=IDENTITY,
                head_repository=IDENTITY, check_suite_id=run_id * 10)
            self.checks.append(dict(
                id=index, name=name, head_sha=SHA, status="completed", conclusion="success",
                app={"id": 15368, "slug": "github-actions"},
                check_suite={"id": run_id * 10},
                details_url="https://github.com/%s/actions/runs/%s/job/%s" % (REPO, run_id, index)))

    def api(self, path):
        if path == "/branches/main":
            return {"commit": {"sha": SHA}}
        if path.startswith("/commits/"):
            page = int(path.rsplit("=", 1)[1])
            return {"check_runs": self.checks[(page - 1) * 100:page * 100],
                    "total_count": len(self.checks)}
        return self.runs[path.rsplit("/", 1)[1]]

    def ready(self):
        return gate.ci_ready(REPO, "123", SHA, self.api)[0]

    def test_all_three_trusted_exact_sha_main_push_checks_are_required(self):
        self.assertTrue(self.ready())
        self.checks.pop()
        self.assertFalse(self.ready())

    def test_pending_failure_and_newer_failed_attempt_cannot_use_old_success(self):
        for status, conclusion in [("queued", None), ("in_progress", None), ("completed", "failure")]:
            with self.subTest(status=status):
                newer = dict(self.checks[0], id=1000, status=status, conclusion=conclusion)
                self.checks.append(newer)
                self.assertFalse(self.ready())
                self.checks.pop()
        self.runs["10"]["status"] = "in_progress"
        self.assertFalse(self.ready())

    def test_other_apps_and_fork_or_pr_checks_cannot_supply_required_names(self):
        for mutation in ("app", "fork", "pull_request", "branch", "workflow", "path", "sha", "suite"):
            with self.subTest(mutation=mutation):
                checks, runs = copy.deepcopy(self.checks), copy.deepcopy(self.runs)
                if mutation == "app":
                    self.checks[0]["app"] = {"id": 7, "slug": "github-actions"}
                elif mutation == "fork":
                    self.runs["10"]["head_repository"] = {"full_name": "fork/csdid", "id": 9}
                elif mutation == "pull_request":
                    self.runs["10"]["event"] = "pull_request"
                elif mutation == "branch":
                    self.runs["10"]["head_branch"] = "feature"
                elif mutation == "workflow":
                    self.runs["10"]["name"] = "Unrelated workflow"
                elif mutation == "path":
                    self.runs["10"]["path"] = ".github/workflows/forged.yml"
                elif mutation == "sha":
                    self.checks[0]["head_sha"] = "b" * 40
                else:
                    self.checks[0]["check_suite"]["id"] = 9
                self.assertFalse(self.ready())
                self.checks, self.runs = checks, runs

    def test_untrusted_duplicate_name_does_not_replace_real_github_actions_result(self):
        self.checks.append(dict(self.checks[0], id=9999, conclusion="failure", app={"id": 7}))
        self.assertTrue(self.ready())

    def test_unverifiable_details_url_fails_closed(self):
        self.checks[0]["details_url"] = "https://attacker.invalid/actions/runs/10/job/1"
        with self.assertRaises(gate.GateError):
            self.ready()

    def test_pagination_and_superseded_main_head(self):
        self.checks = [{"name": "irrelevant"}] * 100 + self.checks
        self.assertTrue(self.ready())
        original = self.api
        def moved(path):
            return {"commit": {"sha": "b" * 40}} if path == "/branches/main" else original(path)
        self.assertFalse(gate.ci_ready(REPO, 123, SHA, moved)[0])

    def test_candidate_requires_main_push_same_repository_or_main_manual_dispatch(self):
        run = dict(self.runs["10"], name="CI")
        event = {"repository": IDENTITY, "workflow_run": run}
        self.assertEqual(gate.candidate("workflow_run", event, REPO, 123, "refs/heads/main", "b" * 40), SHA)
        for field, value in [("event", "pull_request"), ("head_branch", "feature"),
                             ("conclusion", "failure"), ("head_sha", "main"),
                             ("head_repository", {"full_name": REPO, "id": 999})]:
            bad = copy.deepcopy(event)
            bad["workflow_run"][field] = value
            with self.subTest(field=field), self.assertRaises(gate.GateError):
                gate.candidate("workflow_run", bad, REPO, 123, "refs/heads/main", SHA)
        manual = {"repository": IDENTITY, "ref": "main"}
        self.assertEqual(gate.candidate("workflow_dispatch", manual, REPO, 123, "refs/heads/main", SHA), SHA)
        with self.assertRaises(gate.GateError):
            gate.candidate("workflow_dispatch", manual, REPO, 123, "refs/heads/feature", SHA)


class BackendWheelChecks(unittest.TestCase):
    def setUp(self):
        self.files = {"drdid/__init__.py": b"", "drdid/_version_.py": b"version = '1.1.7'\n"}

    def wheel(self, files=None, version="1.1.7", pure="true"):
        output = io.BytesIO()
        with zipfile.ZipFile(output, "w") as archive:
            for name, content in (self.files if files is None else files).items():
                archive.writestr(name, content)
            archive.writestr("drdid-1.1.7.dist-info/METADATA", "Name: drdid\nVersion: " + version + "\n")
            archive.writestr("drdid-1.1.7.dist-info/WHEEL", "Root-Is-Purelib: " + pure + "\n")
        value = output.getvalue()
        artifact = {"packagetype": "bdist_wheel", "filename": "drdid-1.1.7-py3-none-any.whl",
                    "url": "https://files.pythonhosted.org/packages/reviewed.whl",
                    "digests": {"sha256": hashlib.sha256(value).hexdigest()}}
        return value, {"info": {"version": "1.1.7"}, "wheel": artifact, "urls": [artifact]}

    def test_wheel_must_match_pypi_digest_all_python_bytes_version_and_purity(self):
        value, release = self.wheel()
        gate.verify_backend_wheel(value, release, "1.1.7", self.files)
        bad = copy.deepcopy(release)
        bad["wheel"]["digests"]["sha256"] = "0" * 64
        with self.assertRaises(gate.GateError):
            gate.verify_backend_wheel(value, bad, "1.1.7", self.files)
        for changes in (dict(self.files, **{"drdid/extra.py": b"unexpected"}),
                        {"drdid/__init__.py": b"changed"},
                        dict(self.files, **{"sitecustomize.py": b"unexpected"}),
                        dict(self.files, **{"injected.pth": b"import unexpected"}),
                        dict(self.files, **{"drdid/_estimators.so": b"unexpected"})):
            value, release = self.wheel(changes)
            with self.assertRaises(gate.GateError):
                gate.verify_backend_wheel(value, release, "1.1.7", self.files)
        for options in ({"version": "1.1.8"}, {"pure": "false"}):
            value, release = self.wheel(**options)
            with self.assertRaises(gate.GateError):
                gate.verify_backend_wheel(value, release, "1.1.7", self.files)

    def test_real_reviewed_git_checkout_pypi_missing_retry_and_success(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            peer = root / "peer"
            peer.mkdir()
            def git(*args):
                return subprocess.check_output(["git", "-C", str(peer), *args], stderr=subprocess.DEVNULL).decode().strip()
            git("init", "-q")
            git("config", "user.email", "release-test@example.invalid")
            git("config", "user.name", "Release test")
            git("config", "commit.gpgsign", "false")
            for name, content in self.files.items():
                target = peer / name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(content)
            git("add", ".")
            git("commit", "-qm", "reviewed source")
            lock = root / "lock.json"
            payload = dict(schema_version=1, repository="d2cml-ai/DRDIDpy", commit=git("rev-parse", "HEAD"), version="1.1.7")
            lock.write_text(json.dumps(payload))
            with patch.object(gate, "fetch_json", side_effect=HTTPError("url", 404, "Not found", {}, None)):
                with self.assertRaisesRegex(gate.GateError, "Publish reviewed DRDID 1.1.7 first"):
                    gate.backend_ready(lock, peer)
            value, release = self.wheel()
            with patch.object(gate, "fetch_json", return_value=release), patch.object(gate, "urlopen", return_value=io.BytesIO(value)):
                self.assertIn("matches every reviewed Python source byte", gate.backend_ready(lock, peer))
            payload["commit"] = "b" * 40
            lock.write_text(json.dumps(payload))
            with self.assertRaisesRegex(gate.GateError, "exact reviewed commit"):
                gate.backend_ready(lock, peer)

    def test_lock_cannot_select_another_repository_or_short_commit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lock.json"
            valid = dict(schema_version=1, repository="d2cml-ai/DRDIDpy", commit=SHA, version="1.1.7")
            for field, value in [("repository", "fork/DRDIDpy"), ("commit", "main"), ("schema_version", True)]:
                payload = dict(valid, **{field: value})
                path.write_text(json.dumps(payload))
                with self.subTest(field=field), self.assertRaises(gate.GateError):
                    gate.read_lock(path)


if __name__ == "__main__":
    unittest.main()
