# csdid 0.4.3

This release aligns the supported estimators with recent fixes in R `did`
and `DRDID`, speeds up both computation paths, and makes `faster_mode=True`
the default. Pass `faster_mode=False` to select the standard path explicitly.

## Changes that can affect results

- Correct the signs of the pre-period outcome-regression terms in the
  repeated-cross-section DR influence function. ATT point estimates are
  unchanged, while analytical and bootstrap inference can change.
- Apply `min_e` and `max_e` to simple and group aggregations in original
  event-time units, including irregular observation dates. Handle empty
  windows, excluded groups, and samples with only pretreatment cells explicitly.
- Recheck cohorts and available controls after balancing a panel. Resolve
  panel balance before routing varying weights through RC estimators.
- Require unique observed IDs for genuine repeated cross sections. For units
  observed repeatedly, use `panel=True` and, if appropriate,
  `allow_unbalanced_panel=True`; for true repeated cross sections, supply
  distinct row IDs or omit `idname`.
- Preserve preprocessing's two-period `cband=False` decision and use a
  pointwise critical value when `cband=False`, including with bootstrap SEs.
- Validate sampling weights on the retained finite sample before balancing.
  Reject retained negative or all-zero weights while preserving missing-row
  exclusion.
- Give NaN estimator outputs missing influence functions throughout both
  computation paths. Regression-only RC estimation checks only the control
  outcome designs it actually fits; DR also checks treated outcome designs.

The reference R sources are `did` master `74b88eb` (2.5.1.902) and `DRDID`
master `eaac3b4` (1.4.0). These changes include the relevant behavior in
[DID #273](https://github.com/bcallaway11/did/pull/273),
[DID #275](https://github.com/bcallaway11/did/pull/275), and
[DRDID #20](https://github.com/pedrohcgs/DRDID/pull/20).
Unrelated R estimators and anticipation changes absent from the controlling
`did` master branch are not introduced.

## Performance

- Vectorize the per-unit panel-balancing filter shared by both modes.
- Avoid repeated sorting and grouping of validated standard-path panel cells.
- Align fast-path period arrays once and use masks for each comparison,
  preserving unit order and covariate matrix layout.
- Reuse identical balanced-panel DR/IPW overlap checks and propensity fits in
  faster mode. Cache keys include the exact sample, covariates, and weights;
  retained model state is bounded to one fit per cohort and public call.
- Group bootstrap influence rows by cluster once while preserving their
  original summation order. Small cluster counts and object labels retain
  the original loop.
- Use direct original-time lookups, reuse aggregation weight-influence
  calculations, and cache cluster metadata within each aggregation call.
- Generate the same bootstrap multiplier draws with fewer temporary arrays,
  preserving chunk shapes and random streams.
- Load pyplot when plotting is requested while retaining the plotting helpers.

The estimator formulas, overlap checks, Wald-test filtering, bootstrap draw
schedule, and explicit standard-mode option are preserved. Performance
measurements and reproduction instructions are in
[FASTER_MODE.md](csdid/test_csdid/FASTER_MODE.md).

## Validation

The complete suite passes 1,581 tests. The initial performance changes improved
the measured workflows by 1.22x to 2.02x. A further audit against that already
optimized baseline measured faster-mode gains of 1.93x for the core analytical
panel, 2.57x with 12 covariates, 2.72x with 24 periods, and 3.84x with 3,000
bootstrap clusters. These are medians of five serial rounds on the recorded
machine, not universal speed guarantees.

The additional optimizations preserve exact outputs, complete influence
functions, warnings, and random state in the 84-configuration comparison grid
and 36 independent adverse cases. Regression tests cover cache eligibility,
warning filters, custom hooks, missing cluster labels, and numerical reduction
order. This is empirical evidence within the tested environments, not a
proof of equivalence for all possible inputs or future dependencies. The
result-changing fixes listed above are separate from these optimizations.
See the [performance audit](csdid/test_csdid/PERFORMANCE_AUDIT_20261023.md) for
raw measurements, source hashes, R/PyFixest comparisons, and remaining limits.

Current-source R comparisons and their remaining solver precision differences
are documented in
[UPSTREAM_SYNC_v0.4.3.md](csdid/test_csdid/UPSTREAM_SYNC_v0.4.3.md).
