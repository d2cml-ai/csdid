# Installed Python package integration gate

This gate builds exact CSDID and reviewed Python DRDID Git commits, installs
both candidate wheels together in a fresh virtual environment, and runs their
complete test suites. It then runs the reviewed 72-case direct-kernel bridge,
which compares raw ATT/full influence output and the declared original outcome
units for the six shared estimators.

`python_drdid.json` identifies the reviewed peer repository, full commit and
package version. A requested peer SHA must match it. The driver verifies source,
wheel and installed runtime file hashes; every loaded `csdid`/`drdid` module in
the tests and bridge must come from those installed wheels. Plain test snapshots
contain fixtures/helpers but no runtime source packages. Optional CSDID
`scripts/tests` and their helpers are included when present.

Use clean or working Git checkouts containing the required commits. Working
tree changes are not packaged: the driver exports the selected committed Git
trees. No release is published, dependency requirement is weakened, R reference
is regenerated, or numerical tolerance is changed.

```sh
python csdid/test_csdid/joint_sync/run.py \
  --csdid /path/to/csdid --drdid /path/to/DRDIDpy \
  --output /tmp/new-empty-joint-package-run
```

The default CSDID revision is its checkout's `HEAD`; the peer revision comes
from the lock. `--csdid-ref` and `--drdid-ref` accept full commit IDs, and an
explicit peer revision must equal the lock. `--python` selects the interpreter
used to create the fresh environment. CI uses `--skip-slow` to deselect only
explicitly marked slow tests; local runs include them by default.

The output must be new or empty and outside the source repositories. It retains
Git archives, built wheels, logs, JUnit output, module origins, immutable test
snapshots and `receipt.json`. Build/install failures, test failures, a failed
bridge, source imports, mismatched versions/hashes or modified snapshots produce
a nonzero exit and a failed receipt. Normal pip caches remain available; the
fresh environment resolves the packages' actual declared dependencies.

The small failure probes reject an unreviewed peer before any child process,
a source-checkout import and a modified same-version wheel. The gate supplements
both packages' independent R parity suites: agreement between Python packages
alone does not establish agreement with R.
