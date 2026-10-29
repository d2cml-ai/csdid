# Following upstream R changes

The reviewed numerical fixtures remain pinned. A separate GitHub Actions
[workflow](../../../.github/workflows/r-upstream-parity.yml) checks the latest
default branches of `bcallaway11/did` and `pedrohcgs/DRDID`, generates temporary
R references, and compares them with this Python checkout. Every unreviewed
tracked R Git-tree change requires review, even when the fixed numerical grid
still passes. It never commits reference updates or ports R assertions
automatically.

## When checks run

On pushes to Python `main`, pull requests targeting `main`, manual dispatch,
and the first of each month at 09:17 UTC. Main and pull-request checks include
documentation-only changes. Superseded runs are cancelled within the same
event type and ref; monthly or manual checks cannot cancel a main-push check.

The workflow is also reusable by the reviewed R repository callers. An R caller
supplies its exact pushed commit through `did-ref` or `r-drdid-ref`; omitted
inputs select the current R default branches. Each report records the commits
actually checked. A monthly check validates the current R upstreams even when
Python source has not changed. Validation does not publish Python packages;
scheduled PyPI catch-up runs have been removed.

The workflow needs only `contents: read`; results appear in the job summary and
an artifact retained for 14 days. Setup failures also retain available logs.

## What each run establishes

1. Run the ordinary Python suite, including the unchanged pinned numerical gate
   and automation failure tests.
2. Check out the two official R default branches once and retain their exact
   commits. Install their dependencies and selected package sources.
3. Compare each R Git tree with the reviewed baseline. Validate baseline source
   hashes and the complete did/DRDID test-file crosswalk before reporting drift.
   Every unreviewed tracked-tree difference requires review: runtime R/native code,
   tests, helpers, fixtures, snapshots, documentation and build metadata.
   Test-tree and other source changes are reported separately, so a source-only
   change does not imply that R tests changed or that new behavior is covered.
   The two exact, reviewed downstream-caller workflow additions and one bounded
   did version metadata edit are recorded separately as described below; their
   raw tree changes remain visible.
4. Write a strict `source-lock.json` containing the selected commits and package
   versions. Copy the five frozen input/multiplier CSVs to a new directory.
5. Generate all 25 R configurations in both R modes, then run all 50 Python
   numerical workflows and the comparator/provenance checks against those new
   references. Numerical tolerances, support and nonfinite-mask rules are
   unchanged. Recheck source identities and the reviewed fixtures after running.

The source-lock override is opt-in and requires a separate reference directory.
It changes expected source identity only. Ordinary Python CI still requires the
reviewed source pins; a lock cannot bypass artifact or numerical comparisons.
R is needed for the upstream workflow, not for ordinary pinned tests.

| Exit | Meaning | Required response |
|---:|---|---|
| 0 | Numerical comparisons passed; no unreviewed upstream tracked-tree changes | Retain the recorded source identities and comparison evidence |
| 1 | Python and the generated R output disagree | Diagnose the discrepancy; preserve the evidence |
| 2 | Evidence or execution failed, including source, schema, generation or environment errors | Fix the failure and rerun; this is not a parity pass |
| 3 | Numbers passed, but the upstream tracked tree has unreviewed changes | Review source and test changes, port applicable assertions, and update the reviewed baseline explicitly |

Numerical comparisons still run when the change audit returns 3, including for
source-only changes. If upstream drift and numerical failures both occur, exit 1
takes precedence and the audit remains in the report. Invalid evidence or failed
execution remains exit 2. The runner retains its existing
`coverage_review_required` status name for exit 3. A passing fixed grid does not
demonstrate coverage of new R behavior, tests or estimators outside csdid's
supported routes.

## Local replay

Use clean R Git checkouts with enough history to contain the reviewed commits.
The CLI uses supplied revisions without fetching, advancing or modifying them.
Install the Python package/test dependencies, R package dependencies and R
`pkgload`, `jsonlite`, `digest`, and `testthat`, then run from this repository:

```sh
python csdid/test_csdid/upstream_sync/run.py \
  --did /path/to/did --drdid /path/to/DRDID \
  --output /tmp/csdid-r-upstream-new-run
```

The output must be new or empty and separate from source checkouts and reviewed
fixtures. It contains the source lock, upstream change audit, generated references,
R/Python logs, JUnit output, `report.json` and `summary.md`. The workflow stores
these under `comparison/`, alongside checkout identities and Python self-test
logs. One numerical thread is used throughout.

For an offline upstream review without R or Python estimator dependencies:

```sh
python csdid/test_csdid/upstream_sync/audit.py \
  --did /path/to/did --drdid /path/to/DRDID \
  --output /tmp/r-test-changes.json --summary /tmp/r-test-changes.md
```

## Reviewed automation and version metadata

The [review record](reviewed_automation.json) accepts only the exact added
`.github/workflows/python-downstream.yml` file reviewed on 2026-10-29 in each R
repository. These callers run the Python compatibility workflow after an R
push; they do not change R package sources, tests, build inputs or numerical
fixtures. This limited review does not advance the pinned numerical baseline
or claim coverage of additional R assertions.

Each automation record binds the original numerical baseline commit, added path, SHA256,
regular-blob kind and `100644` mode. The path must be absent from that baseline.
Different bytes or file mode, any extra path, and every test-tree change still
require review. There is no directory ignore rule. Missing, malformed, duplicate
or baseline-mismatched review records fail as invalid evidence.

One separate did `version_metadata` record anticipates the existing R
`bump-version.yml` policy after the caller PR merges: `2.5.1.904` becomes
`2.5.1.905` in `DESCRIPTION`. It binds both complete file hashes and both version
strings to the original baseline. Both files must be regular `100644` blobs
with exactly one Version field, the first three version components must remain
identical, and the fourth must increase by one. Replacing only the baseline
Version value must reproduce every current file byte. Dependency edits, any
other metadata change, mode changes, runtime/test changes and later version
bumps remain unreviewed. The workflow-only checkout before the bump also
passes; this record does not require a bump or change the existing R policy.

The audit retains `source_changes` and its raw `source_change_count`, annotates
each source change with its review status, and reports reviewed additions,
reviewed version metadata and unreviewed paths with separate counts. Exact
reviewed additions produce `reviewed_automation_only`; an accepted version
edit is labeled `reviewed_version_metadata` and produces `reviewed_changes_only`.
These statuses have audit exit 0 only when there are no unreviewed changes.
The numerical runner must still
generate fresh R output and pass the unchanged comparisons before its overall
exit can be 0. The audit report also records the review file's SHA256.

## Advancing the reviewed baseline

Treat an update as a code review. Inspect the upstream source/test diff, port
applicable assertions and explain exclusions, then diagnose numerical changes.
If accepted, update the explicit source/version pins in `generate.R` and
`test_r_current_parity.py`, and revise both assertion crosswalks to describe the
new commits and all their test files. Regenerate the committed references using
the new clean sources and the existing inputs, following the
[fixture instructions](../r_ref/current/README.md), then run the full Python
suite and the upstream audit. Commit the evidence and any justified runtime
correction together. Do not replace input data, remove cases, widen tolerances
or mark new R tests covered merely to obtain a passing run.

## Validation record

The implementation adds independently reviewable automation to `main` after
#84 merged. It changes no estimator runtime code or statistical target. Acceptance
requires real R generation/comparison, preserved pinned numerical fixtures,
failure tests for numerical/test/schema/source drift, and a successful Linux
workflow. Unit tests use temporary Git histories and injected subprocess
failures; the live run uses actual R sources. These are distinct checks.

The initial 2026-10-24 validation below is historical. The current reviewed
baseline is did `48863b6` (2.5.1.904) and DRDID `eaac3b4` (1.4.0).

The first local live run used did `74b88eb` and DRDID `eaac3b4`, independently
confirmed as GitHub HEADs on 2026-10-24. It generated 25 configurations and passed
all 125 current-reference tests, including 50 numerical workflows. The test-tree
audit found no drift; reviewed artifacts remained unchanged. Generator replay
also confirmed that the source-lock interface changes no existing input, case
or numerical reference payload. The full local suite passed **2,021 tests**, with
25 warnings and no failures or skips. This adds 68 tests: 31 source-lock checks,
22 runner failure cases and 15 temporary-Git audit cases. The Linux workflow
results are recorded in the PR.

After #84 merged on 2026-10-26, #85 was synchronized with `main`. The generator
and test defaults now use the updated reviewed baseline while retaining the
explicit source-lock override and the corrected retained-data RC ID mapping.
The merge changes no estimator runtime files relative to `main`.

The complete integrated suite passes **2,310 tests**, with 25 warnings and no
failures, errors or skips. A normal pinned R replay preserves all five frozen
input/draw CSVs, case definitions, numerical references and R mode checks
byte-for-byte; the manifest changes only its generation timestamp and exporter
hash. A separate live source-lock replay against freshly verified official R
HEADs passes all **125 comparison/provenance checks**, audits all **44 did and
29 DRDID test files**, and reports zero source/test drift. It leaves the reviewed
fixture tree unchanged. Actionlint and an independent integration review pass.

A new end-to-end probe adds one R test in an isolated clean did commit while
leaving runtime R sources unchanged. All 125 numerical checks still pass, but
the runner correctly exits **3 / coverage_review_required**, identifies the
added test and leaves the reviewed fixtures unchanged. This verifies drift
detection from the expanded 44-file did baseline rather than relying only on
the earlier 34-file audit or mocked failure tests.

On 2026-10-27 the review gate was extended to the entire tracked upstream tree.
Runtime source changes now require review even when no R tests change and the
fixed numerical grid still passes. Source and test changes remain separate
evidence categories. The combined pinned-comparison, audit and runner suite
passes **169 tests**, including **44 audit/runner tests**. New temporary-Git
cases exercise runtime R code, native source, documentation and build metadata,
plus numerical-failure precedence when source-only review is required.

A separate live probe starts from did `48863b6` and commits only an executable,
equivalent positive-support comparison change in `R/utility_functions.R`
(`normalized > 0` becomes `0 < normalized`). No R tests change. Actual R
generation succeeds and all **125 numerical/provenance checks pass**, while
the audit and runner both return **3** and identify exactly one source change
and zero test changes. The reviewed fixtures remain unchanged. This live R
probe is distinct from the runner unit tests, which use real Git audits but
stub the R/numerical subprocesses to isolate control flow.
