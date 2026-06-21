# csdid Test Suite — Description & R-Reference Coverage

_Suite: `csdid/test_csdid` · 561 collected · **561 passed**, 0 skipped, 0 failures._

## How many tests use R reference values?

**65 of the 561 passing tests (~12%) compare csdid output directly against
numerical values produced by R `did`.** The rest are *translated from* R's
`testthat` suite (or authored for csdid) and validated via simulated ground
truth, internal consistency (incl. `faster_mode == standard`), or structure.

| R-reference group | Tests | Source of the R numbers |
|-------------------|------:|-------------------------|
| `test_r_parity.py` | 36 | **Live** R `did` 2.5.0 — `r_ref/generate_*.R` (ATT/SE, 4 aggregations, fix_weights, factor covariate, **RC/universal/anticipation/weighted/clustered-SE** on mpdta + sim_data) |
| `test_sim_parity.py` | 24 | **Live** R `did` 2.5.0 — `r_ref/generate_sim_reference.R` (6 simulated datasets × 2 controls × 2 methods) |
| `test_jel_replication.py` | 5 | Published Callaway–Sant'Anna **JEL** article values (computed in R) |
| **Total using R reference values** | **65** | |

### Are the "simulated ground truth" tests actually correct vs R?
**Yes — verified.** Much of the suite checks `att ≈ te` (the known simulated
effect) within a loose tolerance. To confirm those checks are backed by *exact*
R agreement, the underlying simulated datasets were run through R `did` 2.5.0 and
compared cell-by-cell:

| Simulated parity check | Result |
|------------------------|--------|
| Scenarios | 6 datasets (periods 2/4/5/8/10, constant & dynamic effects) × {nevertreated, notyettreated} × {dr, reg} = 24 |
| ATT(g,t) cells compared | 660 |
| Max ATT difference | **6.0e-10** (dr, optimizer tolerance) / **~5e-15** (reg, closed form) |
| Max analytical-SE difference | **1.1e-4** |

So the ground-truth `att ≈ te` assertions (≈23 tests assert this explicitly) are
far looser than the true ~1e-10 agreement with R. This is now locked in as the
24 permanent tests in `test_sim_parity.py`.

---

## Full breakdown by file (Source · Category · What it tests)

| File | Tests | Source | Category | What it tests |
|------|------:|--------|----------|---------------|
| test_parametric_combinations.py | 120 | csdid (authored) | Ground-truth + consistency | est_method × control × base_period × aggregation grid |
| test_att_gt.py | 66 | R testthat `test-att_gt.R` | Ground-truth | core `att_gt`; fix_weights; custom callable est_method |
| test_aggte_comprehensive.py | 46 | R testthat `test-aggte-comprehensive.R` | Ground-truth + structure | simple/group/dynamic/calendar aggregations |
| test_inference.py | 43 | csdid (authored) | Consistency | bootstrap SE ≈ analytical SE; CI coverage |
| test_r_parity.py | 36 | **Live R `did` 2.5.0** | **R reference** | ATT/SE/aggregations/fix_weights/factor/**RC/universal/anticipation/weighted/clustered-SE** |
| test_edge_cases.py | 27 | R testthat `test-edge-cases.R` | Behavior | unbalanced panels, no-never-treated, first-period |
| test_sim_parity.py | 24 | **Live R `did` 2.5.0** | **R reference** | ATT/SE on 6 simulated datasets (2 controls × 2 methods) |
| test_validation.py | 22 | csdid R-sync (Phase 1) | Error/validation | input guards (missing/reserved cols, dtypes, params) |
| test_glance.py | 19 | R testthat `test-glance.R` | Structure | `glance`/DIDparams keys & values |
| test_cluster_analytic.py | 18 | R testthat `test-cluster-analytic.R` | Consistency | analytical vs bootstrap clustered SE |
| test_ggdid.py | 16 | R testthat `test-ggdid.R` | Structure | plotting returns a `Figure` |
| test_review_fixes.py | 28 | csdid R-sync (review) | Regression | round-1/2 fixes, factor covariates, `allow_unbalanced_panel=False`, idname-numeric, universal-base NaN, `'post'` key alias |
| test_tidy.py | 16 | R testthat `test-tidy.R` | Structure | `tidy` output columns/keys |
| test_error_handling.py | 11 | R testthat `test-error-handling.R` | Error/validation | error/warn behavior (est_method, reversal, fix_weights) |
| test_user_bug_fixes.py | 9 | R testthat `test-user_bug_fixes.R` | Regression + ground-truth | GitHub-issue regressions |
| test_notyettreated.py | 6 | csdid R-sync (Phase 2.1) | Ground-truth + behavior | not-yet-treated control |
| test_jel_replication.py | 5 | **JEL article (R)** | **R reference** | published JEL values (2×2, 2×T, G×T) |
| test_clustered.py | 5 | csdid R-sync (Phase 3.2) | Behavior | clustered multiplier bootstrap |
| test_compute_inffunc.py | 5 | csdid R-sync (Phase 4.1) | Consistency | fast-mode == full-mode point estimates |
| test_integration.py | 5 | csdid R-sync (round-1) | Ground-truth | end-to-end pipeline |
| test_percell_failure.py | 4 | csdid R-sync (Phase 2.2) | Behavior | graceful per-cell NA |
| test_analytical_cluster_se.py | 3 | csdid R-sync (Phase 3.1) | Consistency | analytical vs bootstrap clustered SE + formula |
| test_mboot_cluster.py | 3 | R testthat `test-mboot-cluster.R` | Behavior | clustered bootstrap runs (balanced/unbalanced) |
| test_faster_mode_consistency.py | 24 | R testthat `test-faster-mode-consistency.R` | Consistency | `faster_mode` == standard (panel/RC × control × est × base period) |

### Source legend
- **Live R `did` 2.5.0** — reference values generated locally from the installed R package (regenerate via `r_ref/*.R`); tests skip if the CSVs are absent.
- **JEL article (R)** — published values from Callaway & Sant'Anna's JEL replication (computed with R `did`).
- **R testthat `test-*.R`** — translated assertion-for-assertion from the named file in the R `did` package's `tests/testthat/`.
- **csdid (authored)** — written for csdid; not derived from a specific R test file.
- **csdid R-sync (Phase N / review)** — added during the v2.5.x sync work (input validation, bug fixes, clustered inference, compute_inffunc, code-review regressions).

---

## Counts by Source and by Category

| Source | Tests | R-reference? |
|--------|------:|:------------:|
| R testthat (translated from `test-*.R`) | 255 | No (logic only) |
| csdid (authored: parametric grid, inference) | 163 | No |
| csdid R-sync (Phase / review additions) | 78 | No |
| Live R `did` 2.5.0 (parity CSVs) | 60 | ✅ Yes |
| JEL article (R, published) | 5 | ✅ Yes |
| **Total** | **561** | **65 reference** |

| Category | What it checks |
|----------|----------------|
| **R reference values** (58) | csdid == R `did` numbers (live or published) |
| Ground-truth | csdid recovers the known data-generating effect (`att ≈ te`) — now R-verified via `test_sim_parity` |
| Consistency | bootstrap ≈ analytical SE; fast-mode == full-mode; analytical-cluster ≈ bootstrap-cluster |
| Structure / behavior | shapes, dict keys, types, returns/finiteness, plotting |
| Error / validation | correct `ValueError`/`UserWarning` on bad input |

> Source ≠ Category: a file translated from R `testthat` (Source) often validates
> via ground-truth or structure (Category), **not** by comparing to R's numbers.
> Only the 58 "R reference" tests assert equality to actual R output.

### Why this split is appropriate
- The **58 R-reference tests** pin the estimator's *numbers* to R at machine
  precision (including on the simulated datasets the rest of the suite relies on).
- The **behavioral/consistency tests** exercise far more code paths, edge cases,
  and public APIs than would be practical to pin to R numbers, and they run
  **without requiring R**. Together they give depth (numerical fidelity) and
  breadth (coverage).

### To increase R-reference coverage further
Add scenarios to `r_ref/generate_reference.R` / `generate_sim_reference.R`
(e.g., universal base period, anticipation, clustered/weighted, repeated cross
sections) and extend `test_r_parity.py` / `test_sim_parity.py` to load them.
Each added scenario converts a class of behavioral tests into machine-precision
R-parity tests.
