# R did synchronization, 2026-10-26

The supported Python estimators now incorporate the relevant changes merged
in R did PRs #276, #278 and #279. Validation uses clean did
[`48863b6`](https://github.com/bcallaway11/did/commit/48863b646e657f2aaa5cc9e9d2e7b9d764a6a44e)
(2.5.1.904) and unchanged DRDID
[`eaac3b4`](https://github.com/pedrohcgs/DRDID/commit/eaac3b49cd225d0d7706c633858555a5f204e21e)
(1.4.0). Both were verified against the official upstream repositories on
2026-10-26. The existing Python PR #84 remains the delivery branch.

The objective was to identify changes since the previous `74b88eb` audit,
correct affected Python behavior in both execution modes, and compare the
complete numerical results with the new R source. Passing older snapshots
alone, checking point estimates alone, dropping difficult cases, or widening
tolerances would not satisfy that objective. This audit establishes
implementation agreement within its tested scope; it does not establish
identification assumptions or repeated-sample statistical coverage.

| Change | Python behavior and evidence |
|---|---|
| Anticipation window | All affected cells use the last observed unaffected period under either base-period choice; Wald restrictions exclude declared anticipation. Independent weighted difference-in-means ATT/SE/full-IF calculations check the intercept-only reduction. |
| Weighted feasibility | Overlap uses positive locally normalized support; regression feasibility uses the weighted control Gram. The unweighted preliminary path remains separate, and cached overlap checks include the selected weight period. |
| Small RC groups | Every outcome regression that DR or regression actually fits must have at least as many positive normalized-weight observations as design columns. Tests cover backend non-entry, square designs, underflow, and panel routes using RC kernels. |
| Inference units | Finite positive SEs survive changes of units. Exact-zero/nonfinite scales remain unavailable; one multiplier draw preserves dimensions, cannot estimate a scale, and falls back to pointwise critical values. |
| Balanced event studies | `balance_e` uses common originally observed post-treatment support before missing-cell removal or display limits. Explicit ratio-estimator ATT/IF/SE oracles cover 200 randomized cases. |
| Column names | Direct indexing replaces query-string time filtering and counting. Names containing spaces/punctuation and argument-named columns preserve complete results across panel, RC, missing-row and balancing routes. |

R's legacy conditional CvM pre-test has no Python counterpart. Its new
IPW-only restriction and resolved-default fixes are recorded as not applicable;
the Python Wald test is a different statistic. PR #277 changes the R BMisc
version requirement and does not require a Python dependency change.

The complete local suite passed **2,242 tests**, with zero failures, errors
or skips and 25 warnings, on Python 3.12.12/macOS arm64. This adds 289
targeted regressions to the prior 1,953-test branch. The numerical libraries
were NumPy 2.5.4, pandas 3.1.0, SciPy 1.18.1 and statsmodels 0.15.0; the
installed Python `drdid` version was 1.1.6. Numerical threads were restricted
to one. Portable tests use independent calculations and invariances alongside
R comparisons; agreement between the two Python modes alone is insufficient.

The permanent **25 R configurations / 50 Python mode comparisons** retain
their original five frozen input/multiplier CSVs, case manifest and numerical
reference payloads unchanged. Source pins and provenance were refreshed;
all 94 reference/comparator/provenance tests passed. The new source produces
the same answers on these older, moderate-scale, zero-anticipation cases.

A fresh grid contains 144 configurations: three sampling designs, all three
estimators, both base-period choices, both control groups, anticipation
0/1/2, and small/large outcome-scale probes. All public outputs, unit/cell/event
support, nonfinite masks and every cell/aggregate IF component are compared.
The 252 ordinary/small-scale Python mode comparisons pass. A separate
24-configuration irregular-clock/missing-cell `balance_e` grid passes all
48 Python comparisons and the R mode checks. Another 48 weighted-support
comparisons pass, with maximum absolute differences of approximately
`8.0e-15` for ATT, `5.6e-16` for SE and `2.1e-13` for cell IF components.

The large-scale probes retain **36 failed comparisons** at outcome multiplier
`1e10`: near-zero aggregate IF components exceed the fixed componentwise
absolute tolerance after floating-point cancellation. All reported estimates,
SEs, intervals, critical values and Wald outputs pass; support and nonfinite
masks agree. R's own fast/standard comparisons also fail the same tolerance
in 42 IF fields on these probes. A separate 72/72 original-unit diagnostic
supports the cancellation diagnosis. It does not clear the raw-scale failures:
the stress harness preserves them and exits unsuccessfully. The policy remains
`atol=1e-7, rtol=1e-6`, and `atol=rtol=1e-5` for Wald outputs.

The expanded grid also exposed a reference-export defect: R's synthetic RC
observation numbers are assigned after preprocessing filters. Indexing the
original CSV by those numbers mislabels IF rows when early-treated cohorts
drop. The exporter now resolves IDs through retained `DIDparams$data`.
This removed 24 apparent RC support failures without changing estimator
samples, input data, tolerances or permanent numerical references. The first
full-suite attempt overlapped fixture regeneration and read an unfinished
reference file; it is superseded by the complete successful run above.

The updated [R assertion crosswalk](r_test_alignment_did.json) records all
44 current R test files and their hashes. It preserves older gaps and names
the remaining new partial coverage, including Monte Carlo rejection-rate
assertions, R-only pretest internals and some exact diagnostic/backend probes.
It does not claim one-to-one coverage of all R APIs or assertions.

Five serial alternating-source timing rounds compare the updated branch with
the prior PR head `2a36b48`, under the same interpreter/dependencies and one
numerical thread. Eight healthy weighted DR workloads cover both modes,
construction, fitting and all four aggregations. Across 80 timed executions,
complete outputs, full IFs, warnings and RNG states match the corresponding
baseline snapshots exactly. This does not describe result-changing edge cases.

| Workload | Mode | Median before → after | Increase |
|---|---|---:|---:|
| Balanced panel | fast | 117.0 → 120.4 ms | 2.9% |
| Balanced panel | standard | 383.7 → 387.2 ms | 0.9% |
| 12 covariates | fast | 197.8 → 202.7 ms | 2.5% |
| 12 covariates | standard | 684.8 → 687.8 ms | 0.4% |
| Repeated cross sections | fast | 246.2 → 248.6 ms | 1.0% |
| Repeated cross sections | standard | 241.6 → 245.3 ms | 1.5% |
| Unbalanced panel | fast | 189.9 → 193.6 ms | 1.9% |
| Unbalanced panel | standard | 190.4 → 192.9 ms | 1.3% |

These are small measured costs of the correctness guards, not evidence of
unchanged speed. For example, the five fast wide-design measurements range
from 195.5–199.4 ms before and 202.2–203.5 ms after, without overlap.
[Raw timings, ranges, exact-output checks and source hashes](upstream_20261026_performance.json)
are retained. These eight workloads do not bound performance on every design.

Reproduce the permanent gate with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -m pytest csdid/test_csdid/ -q
```

The [reference README](r_ref/current/README.md) gives exact R regeneration
and fresh stress commands. The tracked
[stress harness](r_ref/current/stress_compare.py),
[compressed stress evidence](r_ref/current/stress_comparison_20261026.json.gz),
[weighted-support harness](r_ref/current/support_compare.py), and
[weighted-support results](r_ref/current/support_comparison_20261026.json)
retain the case definitions, source/input hashes, comparisons and failures.
Large generated R arrays are reproducible outputs rather than package data.
The statistical target and API scope remain the supported Python estimators;
no new conditional-testing procedure was introduced.
