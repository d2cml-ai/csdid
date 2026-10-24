# Following upstream R changes

The reviewed numerical fixtures remain pinned. A separate GitHub Actions
[workflow](../../../.github/workflows/r-upstream-parity.yml) checks the latest
default branches of `bcallaway11/did` and `pedrohcgs/DRDID`, generates temporary
R references, and compares them with this Python checkout. It also inventories
upstream test changes. It never commits reference updates or ports R assertions
automatically.

## When checks run

On every push to this Python repository, with no branch or changed-file filters.
This includes documentation-only changes. The pushed revision must contain the
workflow. Each run checks the latest R GitHub revisions available at checkout.

There is no daily schedule or separate pull-request trigger. A newer push to
the same branch or tag cancels its superseded run. A push to an R repository
alone does not start this workflow; its changes are checked on the next Python
repository push.

The workflow needs only `contents: read`; results appear in the job summary and
an artifact retained for 14 days. Setup failures also retain available logs.

## What each run establishes

1. Run the ordinary Python suite, including the unchanged pinned numerical gate
   and automation failure tests.
2. Check out the two official R default branches once and retain their exact
   commits. Install their dependencies and selected package sources.
3. Compare each R Git tree with the reviewed baseline. Validate baseline source
   hashes and the complete did/DRDID test-file crosswalk before reporting drift.
   Every changed file under `tests/` requires coverage review, including helper
   code, fixtures and snapshots. Other tracked changes are reported separately.
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
| 0 | Numerical comparisons passed; no upstream test-tree changes | Review any reported source changes before advancing pins |
| 1 | Python and the generated R output disagree | Diagnose the discrepancy; preserve the evidence |
| 2 | Evidence or execution failed, including source, schema, generation or environment errors | Fix the failure and rerun; this is not a parity pass |
| 3 | Numbers passed, but upstream tests changed | Review and port applicable assertions before updating coverage records |

Numerical comparisons still run when the test-change audit returns 3. If both
test changes and numerical failures occur, exit 1 takes precedence and the
coverage audit remains in the report. A passing fixed grid does not demonstrate
coverage of new R features, tests or estimators outside csdid's supported routes.

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
fixtures. It contains the source lock, test-change audit, generated references,
R/Python logs, JUnit output, `report.json` and `summary.md`. The workflow stores
these under `comparison/`, alongside checkout identities and Python self-test
logs. One numerical thread is used throughout.

For a coverage-only review without R or Python estimator dependencies:

```sh
python csdid/test_csdid/upstream_sync/audit.py \
  --did /path/to/did --drdid /path/to/DRDID \
  --output /tmp/r-test-changes.json --summary /tmp/r-test-changes.md
```

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

The implementation's scope is an independently reviewable automation PR on top
of #84. It changes no estimator runtime code or statistical target. Acceptance
requires real R generation/comparison, preserved pinned numerical fixtures,
failure tests for numerical/test/schema/source drift, and a successful Linux
workflow. Unit tests use temporary Git histories and injected subprocess
failures; the live run uses actual R sources. These are distinct checks.

The first local live run used did `74b88eb` and DRDID `eaac3b4`, independently
confirmed as GitHub HEADs on 2026-10-24. It generated 25 configurations and passed
all 125 current-reference tests, including 50 numerical workflows. The test-tree
audit found no drift; reviewed artifacts remained unchanged. Generator replay
also confirmed that the source-lock interface changes no existing input, case
or numerical reference payload. The full local suite passed **2,021 tests**, with
25 warnings and no failures or skips. This adds 68 tests: 31 source-lock checks,
22 runner failure cases and 15 temporary-Git audit cases. The Linux workflow
results are recorded in the PR.
