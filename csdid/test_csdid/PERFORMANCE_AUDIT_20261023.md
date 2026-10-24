# Performance audit: PyFixest, R did, and csdid

Date: 2026-10-23. Python baseline:
`972dbd0ff67f8a8b9d6099b7778774bd40c0f862`, which already includes the previous
0.4.3 performance improvements. This audit measures additional gains against
that baseline. The candidate is identified by per-file SHA256 hashes in the
benchmark evidence, independently of the eventual commit identifier.
The [compact evidence file](performance_audit_results_20261023.json) retains
raw timing repeats, environments, source hashes, and validation summaries.

## Objective and decision rules

Find transferable techniques in current PyFixest and R `did`/`DRDID`, locate
actual bottlenecks with profiles, implement justified improvements, and test
speed, memory, and numerical behavior on representative and adverse workloads.
Completion requires reproducible before/after evidence and explicit limits on
cross-package comparisons. A source inspection or an unrelated faster estimator
alone does not establish a performance improvement for this package.

The statistical target, sample selection, weighting, nuisance solvers, trimming,
influence functions, inference, warnings, and RNG behavior remain controlling
requirements. No scientific decision or estimator substitution is introduced.
Exact Python output comparison is the acceptance gate; R solver differences
use declared numerical tolerances. A universal "fastest DiD package" claim is
not supported by a finite workload grid with different estimator APIs.

## Implemented changes

1. **Reuse identical panel nuisance fits.** Faster-mode DR/IPW now reuse the
   unweighted overlap result and weighted propensity fit when their inputs
   repeat. Exact sample membership, covariate values/period, and weight
   values/period determine reuse. Outcome regression and influence arithmetic
   still run for every ATT cell. Each cohort retains at most its most recent
   propensity fit, and caches cannot survive a public fit call. Fits that emit
   warnings accepted by the caller's filters, fail to converge, or have
   nonfinite coefficients are not retained. Intercept-only models and
   overridden estimators, guards, or solver methods retain the uncached path.
2. **Group bootstrap clusters once.** Stable grouping replaces a full
   observation scan for every cluster. Original row order within each cluster,
   reduction dtype, empty NaN/NaT groups, bootstrap draw order, and scaling are
   preserved. At most 16 clusters retain the faster original mask loop;
   object labels also retain that loop to preserve object-NaN uniqueness
   semantics. This introduces no dependency or compiled extension.

The first analytical profile had 162 GLM fits for 81 ATT cells, consuming
0.151 of 0.317 profiled seconds. The candidate uses 18 GLM fits for those cells.
The original 24-period profile used 1,058 fits for 529 cells. These are call
counts and profiler diagnostics; the runtime tables use separate unprofiled
executions.

The approaches are motivated by
[R's guard cache](https://github.com/bcallaway11/did/blob/74b88eb07f1faa644df1271055a0f77848f6550c/R/compute.att_gt2.R#L689)
and PyFixest's
[reuse of compatible transformed data](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/pyfixest/estimation/internals/demean_.py#L189)
and [cluster grouping](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/src/crv1.rs#L5).
The Python changes were implemented independently; no upstream code was copied.

## Core before/after measurements

Median of five serial rounds, milliseconds, including preprocessing,
estimation, and all four public aggregations. Source order and mode order
alternate. Both sources use the same interpreter/dependencies and one numerical
thread. Data generation, imports, warmup, serialization, and profiling are
outside measured intervals.

| Workload | Standard before → after | Fast before → after | Fast speedup |
|---|---:|---:|---:|
| 8,000 units × 10 periods, analytical inference | 383.1 → 387.9 | 220.4 → 114.4 | 1.93× |
| Same panel, 400 clusters | 432.2 → 436.4 | 262.2 → 156.3 | 1.68× |
| 2,000 units × 8 periods, 399 bootstrap draws/bands | 268.3 → 269.6 | 204.6 → 164.0 | 1.25× |

Standard-mode differences are about 1% or less and do not establish an
improvement. The nuisance optimization applies to the faster mode, which was
already the package default before this audit.

Environment: macOS arm64, Python 3.12.12, NumPy 2.5.4, pandas 3.1.0, SciPy
1.18.1, statsmodels 0.15.0, Accelerate BLAS. Thread variables were set to one;
no other agent CPU jobs ran during timing. Exact versions, source hashes,
every repeat, ranges, and phase timings are retained in the evidence.

## Broader stress workloads

The extended benchmark runs each source/case in its own process, with five
serial alternating-source rounds and the faster mode explicitly enabled.
Times include construction, fit, and all four aggregations. Except where
specified, cases use 8,000 units, 10 periods, one covariate, unequal weights,
DR, never-treated controls, and analytical inference.

| Variation | Before | After | Speedup |
|---|---:|---:|---:|
| Baseline DR panel | 219.7 ms | 116.0 ms | 1.89× |
| IPW panel | 204.0 ms | 96.3 ms | 2.12× |
| Regression panel | 91.4 ms | 91.0 ms | 1.00× |
| Intercept only, unweighted | 125.8 ms | 124.4 ms | 1.01× |
| Time-varying covariates and weights | 223.2 ms | 170.2 ms | 1.31× |
| Not-yet-treated controls | 370.7 ms | 350.9 ms | 1.06× |
| Universal base period | 223.6 ms | 116.2 ms | 1.92× |
| 12 covariates | 481.8 ms | 187.6 ms | 2.57× |
| 4,000 units, 24 periods, 529 ATT cells | 772.8 ms | 284.2 ms | 2.72× |
| 50,000 units, 500,000 panel rows | 825.2 ms | 503.2 ms | 1.64× |
| RC, 16,000 rows | 243.7 ms | 242.9 ms | 1.00× |
| Unbalanced 2,000-unit, 8-period source panel | 190.4 ms | 188.0 ms | 1.01× |
| Bootstrap, 2,000 units, 8 periods, 399 draws | 202.8 ms | 163.7 ms | 1.24× |
| Bootstrap, 20,000 units, 8 periods, 999 draws | 2,779.0 ms | 2,708.8 ms | 1.03× |
| Bootstrap, 6,000 units, 8 periods, 199 draws, 5 clusters | 201.1 ms | 140.6 ms | 1.43× |
| Same workload, 3,000 clusters | 1,905.2 ms | 496.2 ms | 3.84× |

All 160 original broad executions passed exact complete-output parity. The
two clustered bootstrap rows above use the final implementation's separate
retiming, with another 40 executions covering both modes and exact parity.
In standard mode, five clusters took 287.5 → 286.7 ms (1.00×), and 3,000
clusters took 1,995.6 → 645.5 ms (3.09×). The small-cluster fallback therefore
removed the earlier regression while retaining the large-cluster gain.

Regression,
intercept-only, RC, and unbalanced workloads show essentially unchanged speed;
the panel cache deliberately does not optimize those routes. Changing
not-yet-treated comparison samples permits less reuse. The large unclustered
bootstrap is dominated by multiplier generation/products, limiting the gain
from nuisance caching.

Process peak RSS was also recorded before snapshot serialization. Median
baseline → candidate measurements were 345.8 → 250.2 MiB for 12 covariates,
650.7 → 525.1 MiB for 500,000 rows, and 623.7 → 561.2 MiB for the large
bootstrap. Not-yet-treated controls rose slightly, 318.9 → 324.2 MiB. These
are whole-process high-water marks including imports, data, warmup and retained
results, not isolated estimator allocations or a universal memory guarantee.

## R comparison

Current upstream source checkouts were loaded explicitly with `pkgload`, since
installed packages were older:

- [`did` 2.5.1.902, `74b88eb`](https://github.com/bcallaway11/did/commit/74b88eb07f1faa644df1271055a0f77848f6550c).
- [`DRDID` 1.4.0, `eaac3b4`](https://github.com/pedrohcgs/DRDID/commit/eaac3b49cd225d0d7706c633858555a5f204e21e).

Both languages use the same exported data and aligned estimator/inference
settings. R 4.6.0 uses Accelerate BLAS; numerical and `data.table` threads are
one. Five rounds per mode follow warmups, with separate Rprof runs. The R
medians are from a serial block preceding the Python timing block, not
interleaved language pairs. Python includes small result-extraction/formatting
costs within its public workflow timer; R snapshot extraction is outside it.

| Core workload | Python candidate fast | R standard | R fast |
|---|---:|---:|---:|
| Analytical | 114.4 ms | 152 ms | 181 ms |
| Analytical, 400 clusters | 156.3 ms | 190 ms | 197 ms |
| Bootstrap, 399 draws | 164.0 ms | 141 ms | 151 ms |

Python is faster on these analytical examples; R remains faster on this
unclustered bootstrap workflow. R's `faster_mode` is also not the fastest R
route on these particular datasets. These measurements cannot rank all
estimators, sample shapes, or inference requests.

All compared analytical public outputs and their complete influence functions
pass `atol=1e-7, rtol=1e-6`; Wald fields use `1e-5`. Maximum ATT and analytical
SE differences are `1.08e-10` and `4.46e-12`. The smaller bootstrap dataset's
unit influence matrix has a pre-existing solver discrepancy up to `5.80e-7`
that exceeds the componentwise threshold in both modes. It is reported as a
diagnostic failure, not hidden behind public-output agreement. Native R/Python
bootstrap draws differ, so cross-language bootstrap SE/band equality is not
asserted; same-language before/after bootstrap results are exact.

See [R source findings and replay](performance_audit_r/README.md) and the
[PyFixest source audit](PYFIXEST_SOURCE_AUDIT.md) for pinned implementation links.

## Restricted PyFixest comparison

PyFixest's general DiD menu does not implement the same CS doubly robust
group-time estimator. The fair numerical comparison here is a balanced
two-period panel, one treated cohort, no covariates or weights, and analytical
ATT/SE without small-sample corrections. PyFixest TWFE uses unit-clustered
covariance; its first-difference route uses HC0 and includes panel reshaping
inside the timer. P-values and confidence-interval critical values are not
compared.

PyFixest **0.60.0** was installed in an isolated environment. All four methods
ran in that same environment, which resolved the same numerical package
versions listed above. They each had a warmup, then five rounds with rotating
method order. Default public copying/storage and diagnostics remain included.

| Units | csdid standard DR | csdid fast DR | PyFixest TWFE | PyFixest first differences |
|---|---:|---:|---:|---:|
| 2,000 | 7.44 ms | 6.78 ms | 7.22 ms | 4.54 ms |
| 10,000 | 11.67 ms | 11.50 ms | 12.37 ms | 6.17 ms |
| 50,000 | 33.10 ms | 32.26 ms | 38.17 ms | 13.94 ms |

All 60 measured results passed the independent difference-in-means ATT/SE
oracle; maximum error among warmup results was `7.29e-14`. csdid's faster mode beats the
TWFE route on these examples, while PyFixest's first-difference route is
faster. This is evidence for a future specialized two-period implementation,
not a general ranking of staggered-treatment estimators. The inspected
PyFixest master revision and the benchmark wheel are recorded separately.

## Remaining priorities and rejected routes

| Priority | Measured evidence in candidate profiles | Proposed work and required check |
|---|---|---|
| Large unclustered bootstrap | Multiplier routine 2.422/2.704 s (90%) at 20,000 units and 999 draws | Investigate allocation and multiplication without changing draw streams, chunk shape semantics, or inference. Current 500-million-element cap can permit about 4 GB for an int64 multiplier buffer alone; larger jobs need a deliberate memory-budget design. |
| Many clusters | Cluster grouping/reduction 0.315/0.601 s (52%), about 105,000 sum calls | Reuse cluster metadata, then evaluate batched or compiled reductions. Preserve dtype, missing-label behavior, and NumPy's within-group summation order. |
| Large analytical panels | `wif` 0.169/0.585 s (29%) at 50,000 units | Reduce temporary influence-weight matrices and repeated time mapping. Check full aggregation influence vectors, rounding, and degeneracy guards. |
| RC and unbalanced panels | RC cell routine 75% and 87% of total time, approximately 1,000 DataFrame selections | Precompute period/sample arrays and unit mappings as in R. Preserve period-specific samples, fixed-weight missingness, and unit-level IF aggregation. |
| Simple 2×2 / intercept only | Equivalent first differences 13.94 vs 32.26 ms at 50,000 units | Test a specialized closed-form route against all weights, trimming, inference and diagnostic contracts before enabling it. |

These are hypotheses for subsequent improvements, not claimed speedups.
PyFixest's HDFE/GPU machinery does not directly solve the current low-dimensional
nuisance problem. Its old compressed-regression route is unavailable in the
audited revision; importing that advice would be stale.

A concrete additional bootstrap proposal was tested and **rejected**: mapping
the unchanged integer Bernoulli draws through a float lookup table instead of
the existing integer transform. Across 32 combinations of sizes, matrix
layouts, global/Generator RNGs and draw counts, results and RNG states matched
bitwise. Median speedup was only `0.9792×` (about 2.1% slower), range
`0.9662–1.0026×`, with essentially unchanged traced allocation peaks. No code
from that experiment was adopted. Script and raw results are retained under
`stress/speed_audit_20261023/bootstrap_experiment.*`.

## Correctness and review

The evidence supports unchanged results relative to `972dbd0` for the tested
inputs and recorded environments. It is not a formal proof or a risk-free
guarantee for every input, platform, dependency version, or user replacement
of numerical internals. The cache relies on deterministic nuisance fits with
identical inputs; grouped reductions preserve the original row order and NumPy
sum operation. Conservative fallbacks and permanent regressions protect the
edge cases discovered during this audit.

This scope is distinct from the earlier R-alignment fixes in the same release
PR, which intentionally change some inference and sample-handling results.
Those fixes are already present in this audit's baseline.

The PR's CI also exercises Linux x64 on Python 3.8 with the legacy requirements
and Python 3.12 with current compatible dependencies. Both lanes run all
non-slow tests, including the new regressions. Their status is reported in the
PR separately from the local macOS results below; two CI lanes do not cover
every dependency combination allowed by the package's unbounded requirements.

- **1,581 tests passed**, 16 warnings, no skips or failures.
- **84 source-to-source configurations and 336 aggregation requests** matched
  exactly for all stored public outputs, full unit and aggregate influence
  functions, warnings, and final NumPy RNG state. The grid crosses panel,
  unbalanced panel, RC, both controls/bases/modes, and DR/regression/IPW, with
  weighted factors, irregular dates, varying inputs, and clustered bootstrap.
- A separate independent probe checked **34 adverse configurations plus two
  pre-import override cases**, also exactly. These include zero weights with
  positive group totals, alternating inputs, signed zeros, shuffled IDs above
  `2**53`, separation, singular control designs, and clustered bootstrap.
- The new regressions include 98 cache tests and 55 cluster-grouping tests.
  They check solver call counts, warning filters and exception ordering,
  custom hooks, cache lifetime/memory bounds, float32/float64 cancellation,
  array layout, missing labels, and complete bootstrap/RNG outputs.
- All 60 core benchmark executions passed complete-output comparisons. The
  comparator requires matching keys, shapes and finite/NaN/infinity masks
  before comparing numerical values.

Review changed the implementation before final measurement: retained model
state was bounded to one entry; intercept-only reuse was excluded to preserve
warning suppression; original-function markers were moved beside their
definitions to protect hooks replaced before importing the fast module. The
earlier core run is retained separately and superseded by `core-final`.

The broad run then exposed a small-cluster cost: standard mode with five
clusters took 289.3 → 297.3 ms. A separate 20-scenario experiment found the
original loop faster in all 16 tested configurations with at most 16 clusters;
32-cluster cases were mixed. The final helper therefore keeps the original
loop through 16 clusters. This adds a second uniqueness scan for larger
cluster counts; its measured helper overhead was 2.5–6% in a separate probe,
so both affected public cluster workloads were re-timed in both modes. All
non-clustered function bodies remain unchanged (checked by comparing ASTs),
and their original timed source hashes are preserved rather than
rewritten to claim those runs used later source.
The evidence collector verified that the original core/broad source differs
from the final source only in this helper, and that final cluster measurements
and both final parity grids match the current package hashes. Earlier cluster
timings remain explicitly marked as superseded in the compact evidence.

A further counterexample with object arrays containing NaNs showed that
unique-label ordering can differ when requesting inverse codes. Object labels
now use the original loop, and three seeded regressions reproduce and reject
the faulty route. The original helper and all failed/superseded experiments
remain under `pre-small-cluster/` and `small_cluster_experiment/` in the local
evidence directory. The full suite and source-to-source grids were repeated
after these changes.

## Reproduce

Use the recorded dependency versions. Preserve a baseline checkout, then run
each CPU workload serially:

```sh
git worktree add --detach /tmp/csdid-speed-baseline 972dbd0ff67f8a8b9d6099b7778774bd40c0f862
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python csdid/test_csdid/benchmark_v043.py \
  --baseline /tmp/csdid-speed-baseline --output /tmp/csdid-core --reps 5 --profiles
python csdid/test_csdid/benchmark_speed_audit.py \
  --baseline /tmp/csdid-speed-baseline --output /tmp/csdid-broad --reps 5 --profiles
python csdid/test_csdid/benchmark_speed_audit.py \
  --baseline /tmp/csdid-speed-baseline --output /tmp/csdid-clusters --reps 5 \
  --cases few_clusters,many_clusters --modes true,false --profiles
python csdid/test_csdid/check_speed_parity.py \
  --baseline /tmp/csdid-speed-baseline --output /tmp/csdid-parity
python -m pytest csdid/test_csdid -q --tb=short
```

Use the linked R and PyFixest instructions for their separate comparisons.
Large snapshots, profiles and logs are preserved locally under the ignored
`stress/speed_audit_20261023/` directory and excluded from distributions.
