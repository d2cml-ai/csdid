# csdid Test Results — `csdid/test_csdid`

_Generated: 2026-06-20 · `pytest csdid/test_csdid`_

> **Update:** all features discussed below have since been implemented —
> `est_method`/treatment-reversal validation, custom callable `est_method`,
> `fix_weights`, and `faster_mode` (see `FASTER_MODE.md`). The suite now runs
> **564 passed, 0 skipped, 0 failures**. The skip analysis in this file is
> retained as the historical record that drove that work.

## Summary (current)

| Result | Count |
|--------|------:|
| **Passed** | **564** |
| Skipped | 0 |
| Failed | 0 |
| **Total collected** | **564** |

## Summary (at the time of this report)

| Result | Count |
|--------|------:|
| **Passed** | **460** |
| Skipped | 31 |
| Failed | 0 |
| xfailed | 0 |
| xpassed | 0 |
| **Total collected** | **491** |

Runtime: ~70 s. The suite is the canonical R-`did`-translated test set plus the
new tests added during the v2.5.x sync. **0 failures, 0 unexpected passes.**

## Live R parity

The package is validated against a live install of **R `did` v2.5.0** (R 4.4.2).
Reference values are produced by `r_ref/generate_reference.R` and checked by
`test_r_parity.py` (25 tests) across 5 scenarios (mpdta nevertreated /
notyettreated / covariates / ipw, and simulated data):

| Quantity | Agreement with R |
|----------|------------------|
| ATT(g,t) point estimates | **~1e-16 (machine precision)** |
| Analytical SEs (`bstrap=False`) | within **~0.1–1%** (max abs diff 3.7e-5) |
| Aggregations (simple/group/dynamic/calendar), overall + per-event | match |
| `nevertreated` coercion warning text | matches R verbatim |

## Per-file breakdown

| File | Tests | Notes |
|------|------:|-------|
| test_parametric_combinations.py | 120 | est_method × control × base_period × aggregation grid |
| test_att_gt.py | 66 | core `att_gt` (incl. 21 skipped — see below) |
| test_aggte_comprehensive.py | 46 | all four aggregation types |
| test_inference.py | 43 | bootstrap vs analytical SEs |
| test_edge_cases.py | 27 | unbalanced panels, no-never-treated, first-period |
| test_r_parity.py | 25 | **live R `did` 2.5.0 comparison** |
| test_validation.py | 22 | input validation guards |
| test_glance.py | 19 | tidy/glance summaries (2 skipped) |
| test_cluster_analytic.py | 18 | analytical clustered SEs (3 skipped) |
| test_ggdid.py | 16 | plotting (single-group bug fixed) |
| test_review_fixes.py | 16 | round-1/2 review regression tests |
| test_tidy.py | 16 | tidy output |
| test_error_handling.py | 11 | error/warn paths (4 skipped) |
| test_user_bug_fixes.py | 9 | GitHub-issue regressions |
| test_notyettreated.py | 6 | not-yet-treated control |
| test_clustered.py | 5 | cluster bootstrap |
| test_compute_inffunc.py | 5 | point-estimates-only mode |
| test_integration.py | 5 | end-to-end pipeline |
| test_jel_replication.py | 5 | JEL article replication |
| test_percell_failure.py | 4 | graceful per-cell NA |
| test_analytical_cluster_se.py | 3 | analytical clustered SE |
| test_mboot_cluster.py | 3 | clustered multiplier bootstrap |
| test_faster_mode_consistency.py | 1 | (skipped — see below) |

---

## The 31 skipped tests

All skips are **R features intentionally not (yet) ported to Python**, not
failures. They fall into six groups.

| # | Category | Tests | Why skipped |
|--:|----------|------:|-------------|
| 1 | `faster_mode` | **18** | R's `faster_mode=TRUE` uses `data.table` internals for a 2.5–3× speed-up. csdid has no equivalent, so tests that exercise it (or compare `faster_mode=TRUE` vs `FALSE`) cannot run. |
| 2 | `fix_weights` | **8** | R's `fix_weights` argument (`NULL`/`"varying"`/`"base_period"`/`"first_period"`) controls time-varying weight handling. Not implemented in csdid. |
| 3 | Bootstrap-heavy (`skip_on_cran`) | **2** | R marks these `skip_on_cran` because they run many bootstrap iterations; translated as `skip`. |
| 4 | Custom `est_method` as R function | **1** | The R test passes an R closure as `est_method`. csdid accepts a Python callable, but this specific test wraps an R-only estimator. |
| 5 | `est_method` string validation | **1** | csdid does not reject unknown `est_method` strings with a dedicated error. |
| 6 | Treatment-reversal detection | **1** | R warns when treatment switches off again; csdid has no such check. |

### Where they live
- **faster_mode (18):** `test_att_gt.py` `TestFasterMode` (8) + `TestIFConsistency` (5); `test_faster_mode_consistency.py` (1); `test_cluster_analytic.py` (1); `test_error_handling.py` (1); `test_glance.py` (2).
- **fix_weights (8):** `test_att_gt.py` `TestFixWeights` (7); `test_error_handling.py` (1).
- **Bootstrap-heavy (2):** `test_cluster_analytic.py` (2).
- **Custom est_method (1):** `test_att_gt.py:535`.
- **est_method validation (1):** `test_error_handling.py:141`.
- **Treatment reversal (1):** `test_error_handling.py:146`.

---

## Recommendations

| Category | Priority | Recommendation |
|----------|----------|----------------|
| **`fix_weights` (8)** | **High** | Real feature gap that changes results with time-varying weights. Implement the four modes (`NULL`/`varying`/`base_period`/`first_period`) in `preprocess_did.py` + `compute_att_gt.py`, then un-skip. Until then, add a guard that **errors** when weights vary within unit-period so users aren't silently misled. |
| **`est_method` string validation (1)** | **Low effort, do now** | Reject unknown `est_method` strings (anything not `dr`/`reg`/`ipw`/callable) with a clear `ValueError` in `att_gt`/`compute_att_gt`, then un-skip the test. Quick R-parity win. |
| **Custom `est_method` (1)** | **Low** | csdid already supports a Python callable `est_method`. Rewrite the skipped test with a Python estimator (signature `f(y, post, D, covariates, i_weights) -> (att, inf_func)`) and un-skip; no production code change needed. |
| **Treatment-reversal detection (1)** | **Low–Medium** | Add a preprocessing check that warns when a unit's treatment turns off after turning on (staggered-adoption assumption violated), matching R. Then un-skip. |
| **Bootstrap-heavy (2)** | **Low** | Convert from permanent `skip` to opt-in `@pytest.mark.slow` and run in nightly CI (not on every commit), so they still provide periodic coverage. |
| **`faster_mode` (18)** | **Out of scope / defer** | This is a performance optimization, **not** a correctness feature — results are identical to the standard path. csdid already vectorizes the multiplier bootstrap (~10× faster). Keep skipped and document as non-goal unless a pandas/polars fast path is later pursued; revisit only if profiling shows a bottleneck. |

### Suggested next step
Knock out the three **low-effort** items first (`est_method` validation,
custom-callable test rewrite, mark bootstrap-heavy as `slow`) to drop the skip
count from 31 → ~27 with no risk, then schedule **`fix_weights`** as the next
real feature PR. Treat `faster_mode` as an explicit non-goal.

---

## Can the skipped features be ported to Python?

**None are impossible** — every skipped feature is plain algorithm/validation
logic with no R-only magic. They differ only in effort and value. (Assessment is
grounded in the R `did` 2.5.0 source inspected via `getNamespace("did")`.)

| Feature | Tests | Portable? | Effort | Verdict |
|---------|------:|-----------|--------|---------|
| **Custom `est_method` (callable)** | 1 | **Already works** | none | csdid already supports `callable(est_method)` (`compute_att_gt.py:164,256`); verified a Python estimator reproduces built-in `reg` to 1e-10. Just rewrite the stub with a Python callable. |
| **`est_method` string validation** | 1 | **Yes** | tiny | Reject anything not in `{dr, reg, ipw}` or callable with a clear `ValueError`. ~3 lines. |
| **Treatment-reversal detection** | 1 | **Yes** | tiny | R simply requires `gname` to be time-invariant per unit ("treatment must be irreversible"). One `groupby(id)[gname].nunique()>1` check. |
| **Bootstrap-heavy (`skip_on_cran`)** | 2 | **Yes** | tiny | Not a feature — functionality already works. Un-skip behind `@pytest.mark.slow` for nightly CI. |
| **`fix_weights`** | 8 | **Yes** | moderate | Pure logic: pick each cell's weights from `NULL`/`base_period`/`first_period`, or `varying` (which routes panel data to RC estimators). No R internals. Needs time-varying-weight detection + per-cell weight selection + RC routing in `preprocess_did.py`/`compute_att_gt.py`. |
| **`faster_mode`** | 18 | **Technically yes, not worth it** | very high | R's `faster_mode=TRUE` is a ~540-line tensor reimplementation (`compute.att_gt2`, `pre_process_did2`, `get_did_tensors`, `run_att_gt_estimation`) **purely for speed** — it returns identical numbers to the standard path. A line-by-line port of R's `data.table` code would not be fast in Python; only a NumPy-native rewrite would help. The 5 "IF consistency" tests merely compare `faster_mode=TRUE` vs `FALSE`, so they are moot without it. |

### Bottom line
- **5 tests** are essentially free to recover today (1 already works, 4 are trivial): custom callable, `est_method` validation, treatment-reversal, and the 2 bootstrap-heavy tests.
- **8 tests** (`fix_weights`) are a genuine, self-contained feature — a good next PR; it changes results under time-varying weights, so it has real user value.
- **18 tests** (`faster_mode`) are portable in principle but **low value** (no new correctness) and **high cost**; recommend keeping as an explicit non-goal and instead profiling/vectorizing the existing path if speed ever matters.

Recovering the first two groups would take the suite from **460 passed / 31 skipped** to roughly **473 passed / 18 skipped**, with the remaining skips all attributable to `faster_mode`.

