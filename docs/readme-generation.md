# Maintaining the generated README

Edit `README.qmd`. It is the source for authored prose, badges, examples and
citation instructions. `README.md` is generated; a second `readme.md` must not
be tracked because the two paths alias on case-insensitive filesystems.

The renderer uses Quarto 1.8.25, Python 3.12, the packages pinned in
`docs/requirements-readme.txt`, and the checked-out csdid source. Examples use
the checked-in `data/mpdta.csv` and set NumPy's seed before bootstrap estimation.
Numerical outputs can change after a reviewed estimator correction; random
render-to-render variation should not create documentation updates.

From the repository root, with Quarto installed:

```sh
python3.12 -m venv /tmp/csdid-readme-venv
/tmp/csdid-readme-venv/bin/python -m pip install -r docs/requirements-readme.txt
/tmp/csdid-readme-venv/bin/python -m pip install -c docs/requirements-readme.txt .
QUARTO_PYTHON=/tmp/csdid-readme-venv/bin/python \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONHASHSEED=0 quarto render README.qmd
python scripts/check_readme.py
python -m unittest discover -s scripts/tests -p 'test_*.py'
```

The README workflow renders twice and compares Markdown and PNG hashes before
opening a PR. It may commit only `README.md` and generated figure PNGs. The
content guard also runs in both mandatory Python CI jobs; it rejects the
observed package-name regression, missing badges/examples/citation, duplicate
case-variant README paths, broken figure links and the corrupted Downloads
badge query. Maintainers should still review proposed wording and results.

The bot uses `GITHUB_TOKEN`. GitHub does not start push-triggered workflows
from its automated pushes, so a generated PR still needs the required checks
to run before it can merge. This change does not add credentials or change the
repository's push-only upstream parity trigger.
