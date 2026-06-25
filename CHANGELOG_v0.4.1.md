# csdid 0.4.1

Differential-testing audit follow-up to 0.4.0, against R `did` 2.5.0 / DRDID 1.3.0.
Plus a new self-contained unit-test suite.

## ⚠️ Behavior changes — re-validate results

- **Parallel bootstrap is now seed-reproducible.** For large samples (`n > 2500`) with
  `pl=True, cores>1`, the multiplier bootstrap previously produced different draws run-to-run
  under a fixed seed (parallel workers did not inherit the seeded RNG). It now derives a
  deterministic per-chunk seed from the parent generator, so results are reproducible *and*
  the draws remain independent. **SEs / uniform bands on that path change** (they are now the
  correct, reproducible values) — re-run if you relied on the previous parallel output.
- **Graceful degradation instead of hard errors** (matching R): aggregating with an
  unavailable cluster variable now **warns and proceeds unclustered** (was a `ValueError`);
  a universal `base_period` with a group lacking a pre-period now **warns and drops** that
  group (was a hard error). Calls that previously raised now return results with a warning.

## Fixes

- **Confidence bands honor `alp`** in `plot_aggte` (previously hardcoded a ~95% critical
  value regardless of `alp`).
- **Diagnostics are emitted on the warnings channel**, not `print()`: the "dropped rows due
  to missing data" notice and the aggregation simultaneous-band messages now use
  `warnings.warn` consistently (so they can be caught/filtered).
- **`att_gt` accepts R-style keyword arguments** (`print_details`, `pl`, `cores`, and extra
  `**kwargs`) instead of raising `TypeError`, for compatibility with code ported from R call
  sites (no effect on numeric results).

## Testing

- **New self-contained unit-test suite** in `csdid/test_csdid/` (no R required at runtime):
  per-module tests for the estimator math (`drdid`), bootstrap (`mboot`), aggregation,
  conditioning/overlap guards, the `att_gt` driver, and preprocessing — plus tests ported
  from R `did`'s own `testthat` suite. These target the package's numeric core with
  hand-derived / analytic / property-based expectations.
