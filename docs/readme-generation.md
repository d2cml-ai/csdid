# Maintaining the generated README

Edit `README.qmd`. It is the source for authored prose, badges, examples and
citation instructions. `README.md` is generated; a second `readme.md` must not
be tracked because the two paths alias on case-insensitive filesystems.

The renderer uses Quarto 1.8.25, Python 3.12, the packages pinned in
`docs/requirements-readme.txt`, the checked-out csdid source, and the exact
Python DRDID peer selected by `csdid/test_csdid/joint_sync/python_drdid.json`.
The workflow checks out both repositories into sibling directories and installs
the pair together with the pinned rendering dependencies. It therefore works
before the reviewed backend is published, while retaining the exact package
requirement. A future peer upgrade must update the rendering pin too. Examples use
the checked-in `data/mpdta.csv` and set NumPy's seed before bootstrap estimation.
Numerical outputs can change after a reviewed estimator correction; random
render-to-render variation should not create documentation updates.

From the repository root, with Quarto installed:

```sh
python3.12 -m venv /tmp/csdid-readme-venv
# Keep the Python DRDID checkout outside this repository. Fetch the exact SHA
# recorded in csdid/test_csdid/joint_sync/python_drdid.json and check it out.
/tmp/csdid-readme-venv/bin/python -m pip install \
  -r docs/requirements-readme.txt /path/to/reviewed/python-drdid .
/tmp/csdid-readme-venv/bin/python -m pip check
QUARTO_PYTHON=/tmp/csdid-readme-venv/bin/python \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONHASHSEED=0 quarto render README.qmd
python scripts/normalize_readme_png.py
python scripts/check_readme.py
python -m unittest discover -s scripts/tests -p 'test_*.py'
```

The README workflow renders twice and compares Markdown and PNG hashes before
opening a PR. It may commit only `README.md` and generated figure PNGs. The
PNG normalizer runs after those two-render checks. It restores tracked HEAD
bytes only when every non-IDAT chunk and the inflated IDAT stream match
exactly, preventing compression-only Linux/macOS differences from opening a
redundant PR. It validates PNG chunk framing/CRCs and complete zlib streams,
with 32 MiB file and 64 MiB inflated-data limits, before changing any figure.
Pixel/scanline or metadata changes, new valid figures and removals remain.
The helper uses only the standard library and does not decode pixels; different
scanline/filter encodings are conservatively retained. Its focused tests run
with the other pre-install guards in both Python 3.8 and 3.12 CI lanes. The
content guard also runs in both mandatory Python CI jobs; it rejects the
observed package-name regression, missing badges/examples/citation, duplicate
case-variant README paths, broken figure links and the corrupted Downloads
badge query. Maintainers should still review proposed wording and results.
The joint installed-wheel gate copies these root documentation fixtures into
its immutable plain test snapshot, so the README tests run there too while the
runtime packages continue to come only from the installed wheels.

The existing manual workflow trigger can validate a feature branch before
merging: `gh workflow run gen_readme.yml --ref BRANCH`. Such a run executes the
renders and guards but skips publishing a bot PR. Automatic main-branch runs
retain PR publication.

The bot uses `GITHUB_TOKEN`. GitHub does not start push-triggered workflows
from its automated pushes, so a generated PR still needs the required checks
to run before it can merge. This change does not add credentials or change the
repository's push-only upstream parity trigger.
