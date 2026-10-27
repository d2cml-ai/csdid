#!/usr/bin/env python3
"""Validate exact CSDID/DRDID Git revisions as jointly installed wheels.

No package source checkout is on the test or bridge import path. A fresh venv,
plain test snapshots and module/file origin guards prevent an editable checkout
or an older same-version distribution from supplying the tested implementation.
"""
import argparse
import ast
import datetime
import email.parser
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
import zipfile


DEFAULT_CSDID = next((parent for parent in Path(__file__).resolve().parents
                      if (parent / "csdid" / "att_gt.py").is_file()), Path.cwd())
DEFAULT_LOCK = Path(__file__).with_name("python_drdid.json")
PACKAGES = ("csdid", "drdid")
BRIDGE_COUNT = 72
THREADS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")


class ChildFailure(RuntimeError):
    def __init__(self, stage, returncode):
        super().__init__("{} exited {}".format(stage, returncode))
        self.returncode = returncode
        self.stage = stage


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tree_hashes(folder):
    folder = Path(folder)
    return {path.relative_to(folder).as_posix(): sha256(path)
            for path in sorted(folder.rglob("*"))
            if path.is_file() and "__pycache__" not in path.parts
            and ".pytest_cache" not in path.parts and path.suffix != ".pyc"}


def load_lock(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate peer-lock field: " + key)
            result[key] = value
        return result
    lock = json.loads(Path(path).read_text(), object_pairs_hook=unique)
    if (not isinstance(lock, dict)
            or set(lock) != {"schema_version", "repository", "commit", "version"}
            or type(lock["schema_version"]) is not int or lock["schema_version"] != 1
            or lock["repository"] != "d2cml-ai/DRDIDpy"
            or not isinstance(lock["commit"], str)
            or re.fullmatch(r"[0-9a-f]{40}", lock["commit"]) is None
            or not isinstance(lock["version"], str) or not lock["version"].strip()):
        raise ValueError("Invalid reviewed Python DRDID lock")
    return lock


def verify_peer_ref(lock, requested):
    if requested is not None and requested.lower() != lock["commit"]:
        raise ValueError("Requested DRDID revision differs from reviewed peer lock")
    return lock["commit"]


def bind_archived_peer_lock(csdid_source, supplied, csdid_ref):
    path = Path(csdid_source) / "csdid/test_csdid/joint_sync/python_drdid.json"
    if path.is_file():
        committed = load_lock(path)
        if committed != supplied:
            raise ValueError("Supplied peer lock differs from the committed CSDID release configuration")
        return dict(mode="committed_peer_lock", archived_lock_sha256=sha256(path))
    if csdid_ref == "HEAD":
        raise ValueError("Default CSDID HEAD must contain its reviewed Python DRDID lock")
    return dict(mode="historical_explicit_ref_without_committed_peer_lock")


def _git(repository, *args):
    return subprocess.check_output(["git", "-C", str(repository), *args], text=True).strip()


def export_commit(repository, requested, destination, archive):
    if requested != "HEAD" and re.fullmatch(r"[0-9a-fA-F]{40}", requested) is None:
        raise ValueError("Explicit source revisions must be full 40-hex commit IDs")
    commit = _git(repository, "rev-parse", "--verify", requested + "^{commit}")
    if requested != "HEAD" and commit != requested.lower():
        raise ValueError("Resolved Git revision differs from requested commit")
    destination.mkdir(parents=True)
    with archive.open("wb") as stream:
        subprocess.run(["git", "-C", str(repository), "archive", "--format=tar", commit],
                       stdout=stream, check=True)
    with tarfile.open(str(archive)) as stream:
        for member in stream.getmembers():
            target = (destination / member.name).resolve()
            if (destination.resolve() not in target.parents and target != destination.resolve()
                    or not (member.isfile() or member.isdir())):
                raise ValueError("Source archive contains an unsafe path or link")
        # Members are regular files/directories and every target is contained.
        stream.extractall(str(destination))
    return dict(commit=commit, git_tree=_git(repository, "rev-parse", commit + "^{tree}"),
                archive_sha256=sha256(archive))


def runtime_hashes(source, package):
    return {path.relative_to(source).as_posix(): sha256(path)
            for path in sorted((source / package).rglob("*.py"))
            if "test_csdid" not in path.parts}


def wheel_identity(path, package, expected_files, expected_version=None):
    with zipfile.ZipFile(str(path)) as archive:
        names = archive.namelist()
        metadata_names = [name for name in names if name.endswith(".dist-info/METADATA")]
        if len(metadata_names) != 1:
            raise ValueError("Wheel must contain exactly one distribution metadata record")
        metadata = email.parser.BytesParser().parsebytes(archive.read(metadata_names[0]))
        if metadata["Name"].lower().replace("-", "_") != package:
            raise ValueError("Unexpected wheel package identity")
        if expected_version is not None and metadata["Version"] != expected_version:
            raise ValueError("Wheel version differs from reviewed source or peer lock")
        files = {name: hashlib.sha256(archive.read(name)).hexdigest() for name in names
                 if name.startswith(package + "/") and name.endswith(".py")
                 and "/test_csdid/" not in name}
    if files != expected_files:
        raise ValueError("Wheel runtime files differ from exact Git source archive: " + package)
    return dict(filename=path.name, version=metadata["Version"],
                requires_dist=metadata.get_all("Requires-Dist") or [], sha256=sha256(path),
                bytes=path.stat().st_size, runtime_files_sha256=files)


def collect_and_validate_origins(site_packages, expected, modules=None):
    """Reject source imports, same-version stale wheels and modified runtime files."""
    site_packages = Path(site_packages).resolve()
    modules = sys.modules if modules is None else modules
    origins = {}
    for name, module in list(modules.items()):
        if module is None or not any(name == pkg or name.startswith(pkg + ".") for pkg in PACKAGES):
            continue
        filename = getattr(module, "__file__", None)
        if filename is None:
            raise ValueError("Loaded runtime module has no verifiable origin: " + name)
        path = Path(filename).resolve()
        if site_packages not in path.parents:
            raise ValueError("Runtime module imported outside installed wheels: {} at {}".format(name, path))
        relative = path.relative_to(site_packages).as_posix()
        if relative not in expected or sha256(path) != expected[relative]:
            raise ValueError("Runtime module differs from reviewed wheel/source: " + name)
        origins[name] = str(path)
    if not all(package in origins for package in PACKAGES):
        raise ValueError("Both installed runtime packages must be loaded and verified")
    return origins


def verify_locked_dependency(requirements, version):
    backend = [value for value in requirements
               if re.match(r"drdid(?=\s|[<=>!~(\[;]|$)", value, re.IGNORECASE)]
    exact = r"drdid\s*(?:==\s*{0}|\(\s*==\s*{0}\s*\))".format(re.escape(version))
    if len(backend) != 1 or re.fullmatch(exact, backend[0], re.IGNORECASE) is None:
        raise ValueError("CSDID wheel must require exactly the reviewed DRDID version")


def verify_installed_files(site_packages, expected):
    for name, expected_hash in expected.items():
        path = Path(site_packages) / name
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError("Installed wheel runtime differs from source: " + name)


def source_version(path, variable):
    for statement in ast.parse(path.read_text()).body:
        if isinstance(statement, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == variable for target in statement.targets):
            value = ast.literal_eval(statement.value)
            if isinstance(value, str):
                return value
    raise ValueError("Cannot establish package version from reviewed source")


def _output_path(output, protected):
    output = Path(output).resolve()
    for source in protected:
        source = Path(source).resolve()
        if output == source or source in output.parents or output in source.parents:
            raise ValueError("Output must be outside source checkouts")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Output must be new or empty")
    return output


def _junit_counts(path):
    document = ET.parse(str(path))
    root = document.getroot()
    suites = [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
    result = {key: sum(int(suite.attrib.get(key, 0)) for suite in suites)
              for key in ("tests", "failures", "errors", "skipped")}
    if result["tests"] <= 0:
        raise ValueError("Empty JUnit test collection cannot establish joint validation")
    result["by_package"] = dict(csdid=0, drdid=0)
    for case in root.iter("testcase"):
        name = case.attrib.get("classname", "")
        if "plain_csdid" in name or "csdid_tests" in name:
            result["by_package"]["csdid"] += 1
        elif "plain_peer" in name:
            result["by_package"]["drdid"] += 1
    if not all(result["by_package"].values()):
        raise ValueError("JUnit must establish execution of both complete test snapshots")
    return result


def clean_environment(base=None):
    """Do not let inherited selectors, plugins or fixture overrides weaken CI."""
    environment = dict(os.environ if base is None else base)
    for key in ("PYTHONPATH", "PYTHONHOME", "PYTHONOPTIMIZE", "PYTEST_ADDOPTS", "PYTEST_PLUGINS",
                "CSDID_R_REFERENCE_DIR", "CSDID_R_REFERENCE_SOURCES",
                "DRDID_R_REFERENCE", "DRDID_R_REVISION"):
        environment.pop(key, None)
    environment.update({name: "1" for name in THREADS})
    environment["MPLBACKEND"] = "Agg"
    environment["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    return environment


def run(csdid, drdid, output, lock_path=DEFAULT_LOCK, csdid_ref="HEAD", drdid_ref=None,
        python=sys.executable, skip_slow=False):
    csdid, drdid, lock_path = Path(csdid).resolve(), Path(drdid).resolve(), Path(lock_path).resolve()
    output = _output_path(output, (csdid, drdid, DEFAULT_CSDID))
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    logs.mkdir()
    receipt = dict(schema_version=1, status="running", skip_slow=bool(skip_slow),
                   started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), commands=[])
    receipt_path = output / "receipt.json"
    def save():
        receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    environment = clean_environment()
    stage = "peer-lock"
    def command(arguments, name, env=None, checked=True):
        nonlocal stage
        stage = name
        destination = logs / (name + ".log")
        with destination.open("w") as stream:
            child = subprocess.run([str(arg) for arg in arguments], cwd=str(output),
                                   env=environment if env is None else env,
                                   stdout=stream, stderr=subprocess.STDOUT)
        receipt["commands"].append(dict(stage=name, arguments=[str(arg) for arg in arguments],
                                        returncode=child.returncode, log=str(destination.relative_to(output))))
        save()
        if checked and child.returncode:
            raise ChildFailure(name, child.returncode)
        return child.returncode
    save()
    try:
        lock = load_lock(lock_path)
        peer_ref = verify_peer_ref(lock, drdid_ref)
        receipt["peer_lock"] = lock
        receipt["peer_lock_sha256"] = sha256(lock_path)
        receipt["driver_sha256"] = sha256(__file__)
        stage = "export-exact-sources"
        sources, archives = output / "sources", output / "archives"
        sources.mkdir()
        archives.mkdir()
        receipt["sources"] = {}
        for package, repository, revision in (("csdid", csdid, csdid_ref), ("drdid", drdid, peer_ref)):
            target = sources / package
            receipt["sources"][package] = export_commit(repository, revision, target, archives / (package + ".tar"))
            receipt["sources"][package]["runtime_files_sha256"] = runtime_hashes(target, package)
            if package == "csdid":
                receipt["configuration_binding"] = bind_archived_peer_lock(target, lock, csdid_ref)
        if source_version(sources / "drdid/drdid/_version_.py", "version") != lock["version"]:
            raise ValueError("DRDID source version differs from reviewed peer lock")
        save()
        virtual = output / "venv"
        command([python, "-m", "venv", virtual], "create-fresh-venv")
        executable = virtual / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
        environment["VIRTUAL_ENV"] = str(virtual)
        environment["PATH"] = str(executable.parent) + os.pathsep + environment.get("PATH", "")
        command([executable, "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel", "pytest"],
                "install-build-and-test-tools")
        wheels = output / "wheels"
        wheels.mkdir()
        command([executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
                 "--wheel-dir", wheels, sources / "csdid", sources / "drdid"], "build-exact-wheels")
        receipt["wheels"] = {}
        selected_wheels = []
        for package in PACKAGES:
            matches = list(wheels.glob(package + "-*.whl"))
            if len(matches) != 1:
                raise ValueError("Missing or duplicate candidate wheel: " + package)
            selected_wheels.extend(matches)
            receipt["wheels"][package] = wheel_identity(matches[0], package,
                receipt["sources"][package]["runtime_files_sha256"],
                lock["version"] if package == "drdid" else
                source_version(sources / "csdid/csdid/_version.py", "__version"))
        if receipt["configuration_binding"]["mode"] == "committed_peer_lock":
            verify_locked_dependency(receipt["wheels"]["csdid"]["requires_dist"], lock["version"])
        command([executable, "-m", "pip", "install", *selected_wheels], "install-both-candidate-wheels")
        command([executable, "-m", "pip", "check"], "verify-dependency-consistency")
        command([executable, "-c", "import importlib.metadata as m,json,sys,sysconfig; "
                 "print(json.dumps({'python':sys.version,'site_packages':sysconfig.get_paths()['purelib'],"
                 "'packages':{d.metadata['Name']:d.version for d in m.distributions()}}))"], "installed-metadata")
        installed = json.loads((logs / "installed-metadata.log").read_text())
        receipt["installed"] = installed
        if installed["packages"].get("drdid") != lock["version"]:
            raise ValueError("Installed DRDID distribution differs from reviewed version")
        expected = {}
        for identity in receipt["wheels"].values():
            expected.update(identity["runtime_files_sha256"])
        verify_installed_files(installed["site_packages"], expected)
        stage = "copy-plain-test-snapshots"
        csdid_plain, peer_plain = output / "plain_csdid", output / "plain_peer"
        csdid_plain.mkdir()
        peer_plain.mkdir()
        tests = csdid_plain / "csdid_tests"
        peer_tests = peer_plain / "tests"
        ignore = shutil.ignore_patterns("__pycache__", "*.pyc", ".pytest_cache")
        shutil.copytree(sources / "csdid/csdid/test_csdid", tests, ignore=ignore)
        shutil.copytree(sources / "csdid/csdid/data", csdid_plain / "data")
        shutil.copytree(sources / "drdid/tests", peer_tests, ignore=ignore)
        extra_test_roots = []
        if (sources / "csdid/scripts").is_dir():
            shutil.copytree(sources / "csdid/scripts", csdid_plain / "scripts", ignore=ignore)
            if (csdid_plain / "scripts/tests").is_dir():
                extra_test_roots.append(csdid_plain / "scripts/tests")
        (peer_plain / "scripts").mkdir()
        bridge = peer_plain / "scripts/check_csdid_integration.py"
        shutil.copyfile(sources / "drdid/scripts/check_csdid_integration.py", bridge)
        bridge_cases = json.loads((peer_tests / "r_ref/cases.json").read_text())
        bridge_cases = [case for case in bridge_cases if case["trim_level"] == .995]
        if len(bridge_cases) != BRIDGE_COUNT:
            raise ValueError("Reviewed direct bridge must contain exactly 72 default-trim cases")
        receipt["bridge"] = dict(cases=len(bridge_cases), script_sha256=sha256(bridge),
                                  cases_sha256=sha256(peer_tests / "r_ref/cases.json"),
                                  checks="Unchanged raw ATT/full IF and declared original-unit comparisons")
        snapshots = {"csdid": tree_hashes(csdid_plain), "drdid": tree_hashes(peer_plain)}
        receipt["test_snapshot_sha256"] = {
            name: hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
            for name, hashes in snapshots.items()}
        config = dict(site_packages=installed["site_packages"], expected=expected)
        (output / "origin-config.json").write_text(json.dumps(config, indent=2) + "\n")
        shutil.copyfile(__file__, output / "_joint_gate.py")
        (output / "origin_guard.py").write_text(
            "import json,pathlib\nimport csdid,drdid\n"
            "from _joint_gate import collect_and_validate_origins\n"
            "config=json.loads(pathlib.Path('origin-config.json').read_text())\n"
            "def verify():\n"
            " return collect_and_validate_origins(config['site_packages'],config['expected'])\n"
            "def pytest_sessionstart(session):\n verify()\n"
            "def pytest_sessionfinish(session,exitstatus):\n"
            " pathlib.Path('pytest-module-origins.json').write_text(json.dumps(verify(),indent=2)+'\\n')\n")
        test_environment = dict(environment, PYTHONPATH=os.pathsep.join(
            (str(output), str(tests), str(peer_tests), str(csdid_plain / "scripts"))))
        junit = output / "joint-junit.xml"
        pytest_config = output / "pytest.ini"
        pytest_config.write_text("[pytest]\naddopts =\n")
        pytest_command = [executable, "-m", "pytest", "-c", pytest_config,
                          "--rootdir", output, "--import-mode=importlib", "-p", "origin_guard",
                          tests, peer_tests, *extra_test_roots, "-q", "--tb=short", "--junitxml=" + str(junit)]
        if skip_slow:
            pytest_command.extend(("-m", "not slow"))
        pytest_code = command(pytest_command, "joint-pytest", env=test_environment, checked=False)
        receipt["pytest_returncode"] = pytest_code
        save()
        bridge_wrapper = (
            "import json,pathlib,runpy; import csdid,drdid; "
            "from _joint_gate import collect_and_validate_origins; "
            "c=json.loads(pathlib.Path('origin-config.json').read_text()); "
            "collect_and_validate_origins(c['site_packages'],c['expected']); "
            "runpy.run_path({!r},run_name='__main__'); "
            "pathlib.Path('bridge-module-origins.json').write_text(json.dumps("
            "collect_and_validate_origins(c['site_packages'],c['expected']),indent=2)+'\\n')"
        ).format(str(bridge))
        bridge_code = command([executable, "-c", bridge_wrapper], "direct-72-case-bridge",
                              env=test_environment, checked=False)
        receipt["bridge_returncode"] = bridge_code
        save()
        stage = "verify-final-evidence"
        receipt["test_counts"] = _junit_counts(junit)
        receipt["module_origins"] = {}
        for name in ("pytest", "bridge"):
            path = output / (name + "-module-origins.json")
            receipt["module_origins"][name] = json.loads(path.read_text())
            if not all(package in receipt["module_origins"][name] for package in PACKAGES):
                raise ValueError("Missing installed runtime origin evidence")
        verify_installed_files(installed["site_packages"], expected)
        final_snapshots = {"csdid": tree_hashes(csdid_plain), "drdid": tree_hashes(peer_plain)}
        receipt["test_snapshots_unchanged"] = snapshots == final_snapshots
        if not receipt["test_snapshots_unchanged"]:
            raise ValueError("A plain test snapshot or frozen reference was modified")
        receipt["evidence_sha256"] = {str(path.relative_to(output)): sha256(path) for path in (
            junit, logs / "joint-pytest.log", logs / "direct-72-case-bridge.log")}
        receipt["pytest_returncode"], receipt["bridge_returncode"] = pytest_code, bridge_code
        if pytest_code or bridge_code:
            raise ChildFailure("joint-pytest" if pytest_code else "direct-72-case-bridge",
                               pytest_code or bridge_code)
        if receipt["test_counts"]["failures"] or receipt["test_counts"]["errors"]:
            raise ValueError("JUnit contains failed tests despite a successful child exit")
        receipt["status"] = "passed"
        returncode = 0
    except Exception as error:
        receipt["status"] = "failed"
        receipt["error"] = dict(stage=error.stage if isinstance(error, ChildFailure) else stage,
                                type=type(error).__name__, message=str(error))
        returncode = error.returncode if isinstance(error, ChildFailure) else 2
        if not 0 < returncode < 126:
            returncode = 1
    receipt["completed_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    receipt["returncode"] = returncode
    save()
    return returncode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csdid", type=Path, default=DEFAULT_CSDID)
    parser.add_argument("--csdid-ref", default="HEAD")
    parser.add_argument("--drdid", type=Path, required=True)
    parser.add_argument("--drdid-ref")
    parser.add_argument("--peer-lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-slow", action="store_true", help="CI only: deselect tests explicitly marked slow")
    args = parser.parse_args(argv)
    try:
        return run(args.csdid, args.drdid, args.output, args.peer_lock, args.csdid_ref,
                   args.drdid_ref, args.python, args.skip_slow)
    except Exception as error:
        parser.error(str(error))


if __name__ == "__main__":
    sys.exit(main())
