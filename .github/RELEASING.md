# Automatic PyPI releases

The `Publish to PyPI` workflow publishes the version reviewed in
`csdid/_version.py` after the required CI runs on `main` succeed. A release PR
must increment that version when packaged Python source, build inputs or
dependencies change. PyPI versions are immutable; the workflow does not choose
or increment release versions.

The publisher accepts only the latest `main` commit and verifies the actual
GitHub Actions workflow, repository, event, commit and latest attempt behind
each of the three required checks. Pull request results, feature branch results,
other apps and superseded or failed attempts cannot authorize a release. The
checks are revalidated after any delay before the publishing job starts.

The reviewed Python DRDID release in
`csdid/test_csdid/joint_sync/python_drdid.json` must also be available on PyPI.
Before building CSDID, the publisher checks the backend wheel's PyPI SHA256,
version, purity and every Python source byte against the exact reviewed backend
commit. `setup.py` must require that exact backend version.

If DRDID has not reached PyPI yet, the workflow reports a waiting state. It
retries automatically at minute 17 of each hour, using the same CI and source
checks. No cross-repository token or manual dispatch is needed when DRDID is
published later. `workflow_dispatch` on `main` remains available for diagnosis.

All publication attempts share one concurrency group. Before building, the
workflow checks any existing CSDID wheel and source archive against reviewed
Python source, build inputs, dependency metadata and PyPI hashes. A complete
matching release skips the build and upload. Either artifact can be uploaded
first: an interrupted upload is completed only if its existing artifact matches.
A reused version with changed source, unexpected build configuration or yanked
files fails visibly and needs a new version.

After an upload, a separate clean checkout of the tested commit verifies the
public wheel and source archive again, including source, metadata and PyPI
hashes. Their filenames, hashes and download URLs must also appear in PyPI's
public Simple index. This catches a stale preflight 404 followed by an upload
that skipped an already existing artifact with different contents. It also
waits for the release to become visible to ordinary package installers.
The existing-release skip path performs the same bounded public confirmation in
the read-only gate, so a later scheduled retry can recover from an earlier
visibility timeout without a build, upload or publishing token.

Missing or partially visible publication data can retry six times, with waits
of 5, 10, 20, 40 and 60 seconds (135 seconds total backoff). Visible source,
metadata, hash or index identity mismatches fail immediately. The verification
step has an eight-minute hard timeout including network calls. A timeout or
mismatch fails the publisher even when PyPI accepted the upload; it never
deletes or replaces immutable files. Scheduled retries can finish a matching
partial release after visibility recovers. Client-side checks cannot make the
preflight, upload and public index update one atomic transaction, or prevent a
later owner action from changing publication availability.

## One-time PyPI configuration

The PyPI `csdid` project must register this GitHub Trusted Publisher:

| Field | Value |
| --- | --- |
| Owner | `d2cml-ai` |
| Repository | `csdid` |
| Workflow filename | `pypi.yml` |
| GitHub environment | `pypi` |

Publishing uses GitHub OIDC instead of a saved PyPI password or API token. The
GitHub `pypi` environment must permit automatic publication from `main` without
a required manual approval if unattended publication is desired. Required main
CI remains a separate merge and release condition.

## Platform limits and validation

GitHub may delay or drop an individual scheduled run; the following hourly run
tries again. In a public repository, GitHub disables schedules after 60 days of
inactivity. See [GitHub's schedule documentation](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#schedule).
The catch-up trigger reuses required current-main push CI results. It does not
rerun current R upstream parity hourly; upstream drift is checked by the parity
workflow on the next CSDID push.

The offline counterexamples run in both required Python CI lanes through the
jointly installed wheel gate. Run them directly with:

```sh
python -m unittest discover -s scripts/tests -v
```

They cover event provenance, spoofed or stale CI results, missing backend retry,
matching-release skips, both interrupted-upload orders, source/dependency drift,
archive build configuration and failed artifact hashes. Workflow changes should
also pass `actionlint .github/workflows/pypi.yml`.
