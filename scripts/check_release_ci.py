"""Read-only release checks. This script never installs or uploads a package."""
import argparse
import ast
import email.parser
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from urllib.error import HTTPError
from urllib.parse import urlparse
from urllib.request import Request, urlopen
import zipfile


ACTION_APP_ID = 15368
REQUIRED = {
    "Tests (Python 3.8)": ("CI", ".github/workflows/ci-tests.yml"),
    "Tests (Python 3.12)": ("CI", ".github/workflows/ci-tests.yml"),
    "Python self-tests and current R upstream parity":
        ("R upstream parity", ".github/workflows/r-upstream-parity.yml"),
}


class GateError(ValueError):
    pass


def valid_sha(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{40}", value) is not None


def same_repo(value, repository, repository_id):
    return (isinstance(value, dict) and value.get("full_name", "").lower() == repository.lower()
            and str(value.get("id")) == str(repository_id))


def candidate(event_name, event, repository, repository_id, ref, sha):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repository):
        raise GateError("Invalid repository identity")
    if not same_repo(event.get("repository"), repository, repository_id):
        raise GateError("The event belongs to a different repository")
    if event_name == "workflow_dispatch":
        if ref != "refs/heads/main" or event.get("ref") not in ("main", "refs/heads/main"):
            raise GateError("Manual publication is allowed only on main")
    elif event_name == "workflow_run":
        run = event.get("workflow_run", {})
        if (run.get("event") != "push" or run.get("head_branch") != "main"
                or run.get("status") != "completed" or run.get("conclusion") != "success"
                or run.get("name") not in ("CI", "R upstream parity")
                or not same_repo(run.get("head_repository"), repository, repository_id)):
            raise GateError("Publication requires a successful same-repository main push workflow")
        sha = run.get("head_sha")
    else:
        raise GateError("Unsupported publication event")
    if not valid_sha(sha):
        raise GateError("A complete tested commit SHA is required")
    return sha


def fetch_json(url, token=None):
    headers = {"Accept": "application/json", "User-Agent": "csdid-release-check"}
    if token:
        headers["Accept"] = "application/vnd.github+json"
        headers["Authorization"] = "Bearer " + token
        headers["X-GitHub-Api-Version"] = "2022-11-28"
    with urlopen(Request(url, headers=headers), timeout=30) as response:
        return json.load(response)


def ci_ready(repository, repository_id, sha, api):
    """Require the latest trusted main-push attempt for every required context."""
    if api("/branches/main")["commit"]["sha"] != sha:
        return False, "Skipped: the tested commit is no longer the main head."
    checks = []
    for page in range(1, 101):
        payload = api("/commits/%s/check-runs?filter=all&per_page=100&page=%d" % (sha, page))
        batch = payload["check_runs"]
        checks.extend(batch)
        if len(checks) >= payload["total_count"]:
            break
        if not batch:
            raise GateError("Incomplete check-run pagination")
    else:
        raise GateError("Too many check runs to establish release readiness")
    runs, latest = {}, {}
    for check in checks:
        name = check.get("name")
        app = check.get("app", {})
        if (name not in REQUIRED or app.get("id") != ACTION_APP_ID
                or app.get("slug") != "github-actions" or check.get("head_sha") != sha):
            continue
        details = urlparse(check.get("details_url", ""))
        match = re.fullmatch(r"/" + re.escape(repository) + r"/actions/runs/([0-9]+)/job/[0-9]+",
                             details.path, flags=re.IGNORECASE)
        if details.scheme != "https" or details.netloc != "github.com" or not match:
            raise GateError("Cannot verify the workflow source of " + name)
        run_id = match.group(1)
        if run_id not in runs:
            runs[run_id] = api("/actions/runs/" + run_id)
        run = runs[run_id]
        workflow, path = REQUIRED[name]
        if (run.get("event") != "push" or run.get("head_branch") != "main"
                or run.get("head_sha") != sha or run.get("name") != workflow
                or run.get("path", "").split("@", 1)[0] != path
                or check.get("check_suite", {}).get("id") != run.get("check_suite_id")
                or not same_repo(run.get("repository"), repository, repository_id)
                or not same_repo(run.get("head_repository"), repository, repository_id)):
            continue
        if name not in latest or check["id"] > latest[name][0]["id"]:
            latest[name] = (check, run)
    missing = sorted(set(REQUIRED) - set(latest))
    if missing:
        return False, "Skipped: required trusted main-push checks are missing or pending: " + ", ".join(missing)
    blocked = [name for name, (check, run) in latest.items()
               if check.get("status") != "completed" or check.get("conclusion") != "success"
               or run.get("status") != "completed" or run.get("conclusion") != "success"]
    if blocked:
        return False, "Skipped: required checks have not all succeeded: " + ", ".join(sorted(blocked))
    return True, "All three required GitHub Actions main-push checks succeeded for " + sha


def verify_backend_wheel(wheel_bytes, release, expected_version, expected_files):
    if release.get("info", {}).get("version") != expected_version:
        raise GateError("PyPI backend version does not match the reviewed lock")
    expected_hash = release["wheel"]["digests"]["sha256"]
    if hashlib.sha256(wheel_bytes).hexdigest() != expected_hash:
        raise GateError("Backend wheel does not match its PyPI SHA256")
    with zipfile.ZipFile(io.BytesIO(wheel_bytes)) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise GateError("Duplicate backend wheel entries")
        if sum(item.file_size for item in archive.infolist()) > 64 * 1024 * 1024:
            raise GateError("Backend wheel expanded size exceeds the expected limit")
        code = {name: archive.read(name) for name in names if name.endswith(".py")}
        unexpected_package_files = [name for name in names if name.startswith("drdid/")
                                    and not name.endswith("/") and name not in expected_files]
        if (code != expected_files or unexpected_package_files
                or any(name.endswith((".pth", ".so", ".pyd", ".dll", ".dylib")) for name in names)):
            raise GateError("Published drdid Python files differ from the reviewed source")
        metadata_paths = [name for name in names if name.endswith(".dist-info/METADATA")]
        wheel_paths = [name for name in names if name.endswith(".dist-info/WHEEL")]
        if len(metadata_paths) != 1 or len(wheel_paths) != 1:
            raise GateError("Invalid backend wheel metadata")
        metadata = email.parser.BytesParser().parsebytes(archive.read(metadata_paths[0]))
        wheel = email.parser.BytesParser().parsebytes(archive.read(wheel_paths[0]))
        if metadata.get("Name", "").lower() != "drdid" or metadata.get("Version") != expected_version:
            raise GateError("Backend wheel package/version does not match the reviewed lock")
        if wheel.get("Root-Is-Purelib", "").lower() != "true":
            raise GateError("Expected a pure Python drdid wheel")


def read_lock(lock_path):
    lock = json.loads(Path(lock_path).read_text())
    # The release lock is authored and reviewed alongside the package version.
    if (set(lock) != {"schema_version", "repository", "commit", "version"}
            or type(lock["schema_version"]) is not int or lock["schema_version"] != 1):
        raise GateError("Invalid reviewed backend release lock schema")
    repository, sha, version = lock["repository"], lock["commit"], lock["version"]
    if repository != "d2cml-ai/DRDIDpy" or not valid_sha(sha) or not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        raise GateError("Invalid reviewed backend release lock")
    return lock


def backend_ready(lock_path, peer_path):
    peer = read_lock(lock_path)
    sha, version = peer["commit"], peer["version"]
    def git(*args):
        return subprocess.check_output(["git", "-C", str(peer_path), *args])
    if git("rev-parse", "HEAD").decode().strip() != sha or git("status", "--porcelain").strip():
        raise GateError("Backend source must be clean at the exact reviewed commit")
    paths = git("ls-files", "-z", "--", "drdid").decode().split("\0")
    expected = {path: git("show", sha + ":" + path) for path in paths if path.endswith(".py")}
    if not expected or "drdid/_version_.py" not in expected:
        raise GateError("Reviewed backend Python source is incomplete")
    assignments = [node for node in ast.parse(expected["drdid/_version_.py"]).body
                   if isinstance(node, ast.Assign) and any(
                       isinstance(target, ast.Name) and target.id == "version" for target in node.targets)]
    if len(assignments) != 1 or ast.literal_eval(assignments[0].value) != version:
        raise GateError("Reviewed backend source version does not match the release lock")
    try:
        release = fetch_json("https://pypi.org/pypi/drdid/" + version + "/json")
    except HTTPError as error:
        if error.code == 404:
            raise GateError("Publish reviewed DRDID " + version + " first, then retry this workflow manually on main") from error
        raise
    wheels = [item for item in release["urls"] if item.get("packagetype") == "bdist_wheel"
              and item.get("filename", "").endswith("-none-any.whl") and not item.get("yanked")]
    if len(wheels) != 1:
        raise GateError("Expected one published pure Python backend wheel; publish DRDID first")
    artifact = wheels[0]
    url = urlparse(artifact["url"])
    if url.scheme != "https" or url.netloc != "files.pythonhosted.org":
        raise GateError("Unexpected backend wheel download host")
    with urlopen(artifact["url"], timeout=60) as response:
        wheel_bytes = response.read(32 * 1024 * 1024 + 1)
    if len(wheel_bytes) > 32 * 1024 * 1024:
        raise GateError("Backend wheel exceeds the expected size limit")
    release["wheel"] = artifact
    verify_backend_wheel(wheel_bytes, release, version, expected)
    return "Published drdid " + version + " wheel matches every reviewed Python source byte and PyPI SHA256"


def output(message, ready=None, sha=None):
    print(message)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as stream:
            stream.write(message + "\n")
    target = os.environ.get("GITHUB_OUTPUT")
    if target:
        with open(target, "a") as stream:
            if ready is not None:
                stream.write("ready=" + str(ready).lower() + "\n")
            if sha:
                stream.write("sha=" + sha + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("ci", "backend", "lock"))
    parser.add_argument("--lock")
    parser.add_argument("--peer")
    args = parser.parse_args()
    try:
        if args.mode == "lock":
            lock = read_lock(args.lock)
            with open(os.environ["GITHUB_OUTPUT"], "a") as stream:
                for key in ("repository", "commit", "version"):
                    stream.write(key + "=" + lock[key] + "\n")
            output("Reviewed backend release: " + lock["version"] + " at " + lock["commit"])
            return 0
        if args.mode == "backend":
            output(backend_ready(args.lock, args.peer))
            return 0
        event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
        repo, repo_id = os.environ["GITHUB_REPOSITORY"], os.environ["GITHUB_REPOSITORY_ID"]
        sha = candidate(os.environ["GITHUB_EVENT_NAME"], event, repo, repo_id,
                        os.environ["GITHUB_REF"], os.environ["GITHUB_SHA"])
        checkout = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
        if checkout != sha:
            raise GateError("Release checkout does not match the tested commit")
        token = os.environ["GH_TOKEN"]
        ready, message = ci_ready(repo, repo_id, sha,
            lambda path: fetch_json("https://api.github.com/repos/" + repo + path, token))
        output(message, ready, sha)
        return 0
    except Exception as error:
        output("Release blocked: " + str(error), False)
        return 1


if __name__ == "__main__":
    sys.exit(main())
