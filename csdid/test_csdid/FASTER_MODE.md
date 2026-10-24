# `faster_mode` and 0.4.3 performance validation

The subsequent [PyFixest/R performance audit](PERFORMANCE_AUDIT_20261023.md)
measures additional improvements against commit `972dbd0`, including nuisance-fit
reuse, bootstrap cluster grouping, larger stress workloads, and current-source
R comparisons. The measurements below describe the earlier optimization round.

`ATTgt` defaults to `faster_mode=True` in 0.4.3. Pass `faster_mode=False` to
select the standard computation path. Both modes remain tested across panel,
repeated cross sections, unbalanced panels, estimators, control groups, weights,
base periods and clustered inference. Results across the two paths agree to
numerical tolerance; they need not have identical floating-point representations.

```python
from csdid.att_gt import ATTgt

result = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=df,
               xformla="Y~X").fit(est_method="dr", bstrap=False)
# Set faster_mode=False explicitly to use the standard path.
```

The performance changes vectorize balanced-panel selection, reuse sorted
panel reshaping in standard mode, align period arrays once in fast mode, and
reuse time mappings, cohort-weight calculations and cluster membership inside
each aggregation call. Plotting imports are deferred until requested. The
estimators, overlap decisions, Wald filtering, bootstrap draw sequence and
bootstrap chunk sizes are preserved.

## Measured before/after results

These are actual implemented-source measurements on October 23, 2026, against
an aligned 0.4.3 source snapshot before the performance changes
(`fafac8dc8c3cabf38444296abbb332db58136abe`). Both mode flags were explicit in both sources, so
the table measures implementation changes independently of the new default.

Each number is the median of five serial, paired rounds in milliseconds for
construction/preprocessing, fitting, and all four public aggregations
(`simple`, `group`, `dynamic`, `calendar`).

| Workflow | Standard before | Standard after | Speedup | Fast before | Fast after | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| 80,000 panel rows, analytical inference | 623.4 | 402.3 | 1.55x | 405.7 | 227.9 | 1.78x |
| Same panel, 400 clusters | 752.6 | 436.8 | 1.72x | 540.7 | 268.0 | 2.02x |
| 16,000 panel rows, 399 bootstrap draws and simultaneous bands | 329.4 | 269.5 | 1.22x | 257.4 | 203.0 | 1.27x |

The analytical panel has 8,000 units, 10 periods and 81 ATT(g,t) cells; the
bootstrap panel has 2,000 units, eight periods and 49 cells. Both use seeded
synthetic data, the DR estimator, one covariate, unequal positive sampling
weights and never-treated controls. These gains describe these workloads;
repeated-cross-section or unbalanced-panel fits can benefit less because the
panel-specific changes do not apply to those branches.

The environment was macOS arm64, Python 3.12.12, NumPy 2.5.4, pandas 3.1.0,
SciPy 1.18.1 and statsmodels 0.15.0, with Accelerate BLAS and inherited thread
settings. Sources ran in separate processes with alternating source and mode
order. Dataset generation, imports, warmups, serialization and comparisons
were outside the timers. Profiling was a separate untimed call. No other agent
CPU jobs ran during the benchmark. Other machines and BLAS providers can yield
different results.

[benchmark_results_v043.json](benchmark_results_v043.json) records every raw
run, phase medians/ranges, environment, source hashes and comparison summary.

## Correctness checks

All 60 measured workflow/mode executions matched the aligned baseline exactly
for the same explicit mode: ATT, SE, bounds and critical values, Wald results,
complete unit and aggregation influence functions, warnings and final NumPy RNG
state. The comparator checks symmetric field sets, array shapes and finite,
NaN, positive-infinite and negative-infinite masks before numerical values.

A separate functional comparison was exact for 144 estimator configurations
and 2,304 aggregation requests, plus 48 bootstrap configurations using identical
multiplier matrices. The full package suite passed 1,428 tests. Tests include
shuffled and non-dense IDs, large integer and nullable IDs, time-varying
covariates/weights, factor covariates, missing cells, event windows, cluster
metadata, helper overrides and default-mode behavior.

These before/after Python checks are distinct from numerical agreement with
R's nuisance solvers. [UPSTREAM_SYNC_v0.4.3.md](UPSTREAM_SYNC_v0.4.3.md)
records the upstream revisions, independent oracles, adversarial checks and
remaining R solver-tolerance differences. No R/Python bit-identical claim is
made.

## Reproduce

Create a checkout of the aligned baseline commit shown above, then run the
current benchmark script with the same interpreter and BLAS settings for both
source trees:

```sh
git worktree add --detach /tmp/csdid-v043-before fafac8dc8c3cabf38444296abbb332db58136abe

python csdid/test_csdid/benchmark_v043.py \
  --baseline /tmp/csdid-v043-before \
  --output /tmp/csdid-v043-benchmark \
  --reps 5 --profiles

python -m pytest csdid/test_csdid -q
```

The benchmark always runs both mode flags and has no pass/fail speed threshold.
It fails on a complete-output parity difference. Its output directory preserves
raw timings, comparisons, source hashes, environment, first-pair full reference
arrays and optional cProfile files. Repeated identical reference arrays from
later rounds are removed to limit disk use.

Additional raw R and benchmark evidence is kept locally under the ignored
`csdid/test_csdid/stress/release_v043_20261023/` directory. Those large artifacts
are excluded from the package distribution and the PR.
