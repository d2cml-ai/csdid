# Development environments

`setup.py` is the runtime dependency manifest. It requires the reviewed
`drdid==1.1.7` release; ordinary package installs resolve that exact backend.
Joint development checks use the source revision in
[`python_drdid.json`](../csdid/test_csdid/joint_sync/python_drdid.json).

CI explicitly tests Python 3.8 and 3.12. For a modern local source installation:

```sh
python3.12 -m venv .venv/py312
.venv/py312/bin/python -m pip install pytest -e .
```

The Python 3.8 lane first installs the legacy dependency set, including
NumPy 1.24.3:

```sh
python3.8 -m venv .venv/py38
.venv/py38/bin/python -m pip install -r requirements.txt -e .
```

Use the intended interpreter when creating the environment; with pyenv, set
`PYENV_VERSION` to your installed development interpreter for these commands.
The root
`.python-version` selects the Python 3.11 line for dependency scanning of the
legacy requirements; `docs/.python-version` selects the Python 3.12 line for the
documentation dependencies. The exact selectors, 3.11.14 and 3.12.12, are
recognized by the scanner's pinned pyenv catalog; its runtime manager may use a
newer installed patch within that minor line. These tooling selections do not
define a package minimum or change the explicit Python 3.8/3.12 CI matrix.

README rendering has its own reviewed pins in
[`requirements-readme.txt`](requirements-readme.txt) and runs through
[`gen_readme.yml`](../.github/workflows/gen_readme.yml).

The retired Pipenv files described a June 2023 environment outside the current
documented development and CI paths.
Their historical version/hash snapshot remains available in Git:

```sh
git show eb9d320:Pipfile
git show eb9d320:Pipfile.lock
```

The current pip development paths use the package metadata and lane-specific
requirements; they do not provide a fully hashed transitive environment lock.
