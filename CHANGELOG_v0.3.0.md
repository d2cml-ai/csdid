# csdid v0.3.0 — Sync with R `did` v2.5.1

> **v0.3.1 (patch):** two edge-case fixes from a final GPT-5.5 review of v0.3.0
> (universal base cell identified by position not by an all-zero influence
> function; numeric dtype validation robust to pandas nullable integers). See
> §16. Validated by an independent Opus 4.8 review. 564 tests passing.

## Summary

This release brings the Python `csdid` package into full alignment with the R `did` package v2.5.x. Key additions: comprehensive input validation, critical bug fixes, proper clustered inference (bootstrap **and** analytical), a point-estimates-only mode, a ~10× faster bootstrap, the **`fix_weights`** parameter, and a vectorized **`faster_mode`** (~3× faster, bit-identical). Pre-existing Python defects fixed during this work (unrelated to the R sync) are tracked separately in **`PYTHON_BUGFIXES.md`**.

**564 tests** (`csdid/test_csdid`, 0 skipped) cover all changes; **65 compare directly against live R `did` 2.5.0** (ATT(g,t) to ~1e-16, analytical SEs within ~0.1–1%, across all aggregations, `fix_weights` modes, and factor covariates) or published JEL values. See §10–§16 and `csdid/test_csdid/test_suite_desc.md`.

---

## Changes Made

### 1. Input Validation (commit `0efb95f`)

**Files:** `csdid/attgt_fnc/preprocess_did.py`, `csdid/att_gt.py`  
**Tests:** `test/test_validation.py` (22 tests)

Added validation guards matching R `did` v2.5.1:
- Missing columns → clear error listing available columns
- Reserved internal column names (`w`, `rowid`, `G_m`, `C`, etc.) → rejected
- Non-numeric outcome/time/group variables → rejected
- Negative `gname` values → rejected with explanation
- Duplicate `(idname, tname)` rows → rejected
- `alp` must be in (0, 1); `biters` must be positive integer
- `anticipation` must be non-negative
- Strict `control_group` matching (only 'nevertreated' or 'notyettreated')
- Negative/zero-mean weights → rejected

### 2. Bug Fix: `notyettreated` Last-Cohort Retention (commit `5cf446a`)

**Files:** `csdid/attgt_fnc/preprocess_did.py`, `csdid/attgt_fnc/compute_att_gt.py`  
**Tests:** `test/test_notyettreated.py` (6 tests)

Fixed: when `control_group='notyettreated'` and no never-treated group exists, the last-treated cohort (which serves as the not-yet-treated comparison) was accidentally deleted from the data along with always-treated units. Now retained in data while excluded from estimation groups.

Also fixed: panel data path was not passing `post_treat` to the results, causing the post-treatment indicator to always be 0.

### 3. Bug Fix: Per-Cell Estimation Failure → NA (commit `9c5d4a5`)

**Files:** `csdid/attgt_fnc/compute_att_gt.py`  
**Tests:** `test/test_percell_failure.py` (4 tests)

When a 2×2 DiD estimation fails for a specific (g,t) cell (singular design matrix, too few observations), the package now warns and sets that cell's ATT to `NaN` instead of crashing the entire computation. Both panel and repeated cross-section paths are protected.

### 4. Fix: Cluster Bootstrap Aggregation (commit `8f0ea6f`)

**Files:** `csdid/utils/mboot.py`  
**Tests:** `test/test_clustered.py` (5 tests)

Rewrote the clustered multiplier bootstrap to match R `did` v2.5.1:
- Uses cluster **sums** of influence functions (not means) per Callaway & Sant'Anna (2021) Remark 10
- Fixed `clustervars` normalization (string/list/None handling)
- Fixed time-varying cluster detection (checks full data, not just first period)
- Fixed parallel bootstrap chunk calculation (prevents negative chunk sizes when `biters < cores`)
- Handles 1-d influence function arrays

### 5. Feature: Analytical Cluster-Robust SEs (commit `200054f`)

**Files:** `csdid/aggte_fnc/utils.py`  
**Tests:** `test/test_analytical_cluster_se.py` (3 tests)

When `clustervars` is set and `bstrap=FALSE`, cluster-robust standard errors are now computed analytically from cluster sums of influence functions. Previously only i.i.d. analytical SEs were available.

### 6. Feature: `compute_inffunc` Parameter (commit `4c3291d`)

**Files:** `csdid/att_gt.py`, `csdid/attgt_fnc/compute_att_gt.py`  
**Tests:** `test/test_compute_inffunc.py` (5 tests)

New `compute_inffunc=False` parameter for point-estimates-only runs:
- Returns identical ATT point estimates without influence functions or SEs
- Uses less memory (no n×k IF matrix)
- Runs faster (skips bootstrap entirely)
- `aggte()` is blocked with a clear error message
- Default remains `True` (backward compatible)

### 7. Performance: Vectorized Multiplier Bootstrap (commit `7bc09e0`)

**Files:** `csdid/utils/bmisc.py`

Replaced Python for-loop over `biters` iterations with a single matrix multiplication:
```python
# Before: O(biters) Python loop
for b in range(biters):
    Ub = np.random.choice([1, -1], size=(n, 1))
    outMat[b] = np.mean(inf_func * Ub, axis=0)

# After: single vectorized operation
Ub = np.random.choice([1, -1], size=(biters, n))
outMat = Ub @ inf_func / n
```
~10x faster for typical problem sizes.

---

### 8. Code Review Bug Fixes (commit `de8e0b5`)

**Files:** `csdid/utils/mboot.py`, `csdid/aggte_fnc/utils.py`, `csdid/attgt_fnc/preprocess_did.py`, `csdid/attgt_fnc/compute_att_gt.py`, `csdid/utils/bmisc.py`  
**Tests:** `test/test_review_fixes.py` (13 tests)

GPT-5.5 code review found 10 bugs:

| Bug | Severity | Fix |
|-----|----------|-----|
| Clustered bootstrap SE scaled by `1/√n_clusters` instead of R's `√n_clusters/n` | **Critical** | Corrected to `bSigma * sqrt(n_clusters) / n` |
| Analytical cluster SE used `√(mean(S²)/n_c)` instead of `√(sum(S²))/n` | **Critical** | Corrected formula |
| Panel path overwrote global `n` with cell-specific `len(disdat)` | **Critical** | Renamed to `n_cell` local variable |
| `asif_nev_treated` threshold missing `+ anticipation` | High | Added `+ anticipation` |
| `glist` not recomputed after recoding groups to 0 | High | Added recompute step |
| `treated_fp` threshold missing `+ anticipation` | High | Added `+ anticipation` |
| Panel/RCS count branches swapped | High | Swapped to match R |
| Unbalanced panels not switched to RC estimators | High | Detect unbalanced → set `panel=False` |
| Bad formula silently fell back to intercept | High | Raise `ValueError` on bad user formula |
| Boolean outcomes rejected (R allows logical) | Medium | Accept `np.bool_` dtype |
| Bootstrap chunking could produce wrong iteration count | Medium | Use `np.linspace` split |
| Vectorized bootstrap could OOM on large inputs | Medium | Added chunked allocation |

---

### 9. Second-Round Review Fixes (commit `4ebec68`)

**Files:** `csdid\attgt_fnc\preprocess_did.py`, `csdid\utils\mboot.py`  
**Tests:** `test/test_review_fixes.py` (+6), `test/test_notyettreated.py` (updated)

A second deep review against the R source (`R/pre_process_did.R`) found behavioral divergences:

| Issue | Severity | Fix |
|-------|----------|-----|
| `control_group='nevertreated'` with no never-treated group raised an error | **High** (behavior change) | R v2.5.1 warns + coerces the last cohort to never-treated and truncates late periods. Python now matches. |
| First-period-treated units dropped by cohort membership (`gname in glist_in`) | **High** | R drops by row identity (`~treated_fp`); the old approach could delete the not-yet-treated comparison cohort. Refactored to match R. |
| Unbalanced-panel detection used `nunique(counts) > 1` | Medium | R uses `nrow != n_units * n_periods`. The old heuristic missed panels where every unit has the same (but incomplete) period count. |
| `mboot` could crash when all bootstrap dimensions are degenerate | Low | Added `ncol(bres) == 0` guard returning NaN SE/crit_val (matches R). |

**Validation against canonical `mpdta` (Callaway & Sant'Anna):**
- ATT(g,t) point estimates match the R `did` vignette **exactly** (e.g., ATT(2004,2006) = -0.1373).
- Simple aggregation ATT = **-0.04** (R: -0.0400).
- Group aggregation overall ATT = **-0.031** (R: -0.0310).

**Note (behavior change):** Code relying on the old `ValueError("no available never-treated group")` should now expect a `UserWarning` and a successful run with the last cohort coerced to never-treated.

---

### 10. Test-Suite Consolidation & Live R Parity (commits `…`)

**Files:** `csdid/test_csdid/*`, `csdid/att_gt.py`, `csdid/attgt_fnc/preprocess_did.py`

The new tests were consolidated into the canonical R-translated suite at `csdid/test_csdid/`, and the package was validated against a **live install of R `did` v2.5.0** (R 4.4.2).

**Stale divergence markers removed** (my round-1/2 fixes made these match R — they were previously `xfail`/`xpass`):
- mboot clustervars-as-string bug → clustered bootstrap now runs
- unknown formula variables now raise (not silent intercept)
- anticipation-window coercion of groups treated after `max(t)`
- truly unbalanced panels route through RC estimators (no `vstack` failure)
- `nevertreated` with no never-treated group → warns + coerces (not error)
- first-period-treated drops emit a `UserWarning`

**New behavior/fixes:**
- `weights_name` may equal an internal name like `w` (R stores weights in `.w`); weights are captured before any internal columns are created.
- **Implemented att_gt-level analytical cluster-robust SE** (`bstrap=False` + `clustervars`), matching R; verified equal to the cluster bootstrap SE.
- Fixed `plot_attgt()` for a single group (`plt.subplots(..., squeeze=False)`).

**Live R parity (`csdid/test_csdid/test_r_parity.py`, 25 tests):**
Reference values generated by `r_ref/generate_reference.R` from R `did` 2.5.0 across 5 scenarios (mpdta nevertreated/notyettreated/covariates/ipw, simulated data):
- **ATT(g,t) point estimates match R to ~1e-16** (machine precision).
- Analytical (`bstrap=FALSE`) SEs match within ~0.1–1%.
- All four aggregation types (simple/group/dynamic/calendar) — overall and per-event — match R.
- R's `nevertreated`-coercion warning text matches csdid's wording.

**Full `csdid/test_csdid` suite: 460 passed, 31 skipped, 0 failures, 0 xfails.**

---

### 11. est_method / base_period validation, treatment-reversal, custom callable

**Files:** `csdid/att_gt.py`, `csdid/attgt_fnc/preprocess_did.py`

- `est_method` must be `dr`/`reg`/`ipw` or a callable; `base_period` must be `varying`/`universal` — both now validated with clear errors (matching R).
- Treatment must be irreversible: a time-varying `gname` within a unit now errors (matching R's message).
- Confirmed a **custom Python callable** `est_method` reproduces the built-in estimator (already supported; recovered the skipped test).

### 12. `fix_weights` parameter (`NULL`/`varying`/`base_period`/`first_period`)

**Files:** `csdid/att_gt.py`, `csdid/attgt_fnc/preprocess_did.py`, `csdid/attgt_fnc/compute_att_gt.py`  
**Tests:** `test_att_gt.py::TestFixWeights`, `test_r_parity.py::test_fix_weights_matches_r`

Controls which period's weights each (g,t) cell uses when weights are
time-varying. `base_period` uses the group's last pre-treatment period (R's
`pret_g`); `first_period` uses `tlist[0]`; `varying` routes panel data to the RC
estimators with per-period weights; the RC path drops units not observed in the
target period (with a warning), like R. A "time-varying weights detected"
warning is emitted for the `NULL` default. **All four modes match R `did` 2.5.0
to ~1e-16** on a shared time-varying-weight dataset.

### 13. `faster_mode` — vectorized compute (~3× faster, bit-identical)

**Files:** `csdid/attgt_fnc/compute_att_gt2.py` (new), `csdid/att_gt.py`  
**Tests:** `TestFasterMode`, `TestIFConsistency`, `test_faster_mode_consistency.py`  
**Doc:** `csdid/test_csdid/FASTER_MODE.md`

`ATTgt(..., faster_mode=True)` builds the covariate design **once per period**
(panel) / **once globally** (RC) and **pre-pivots** outcomes/weights, removing
the per-(g,t) `patsy` + `panel2cs2` cost. Results are **identical** to the
standard path (max |std−fast| = 1.4e-17 for ATT & SE across 72 configs).
Python(Fast) is ~3.3× faster than Python(Standard), roughly matching
R(Standard) in absolute time.

### 14. Factor / transformed covariates (R v2.5.1 "factor variables")

**Files:** `csdid/attgt_fnc/preprocess_did.py`, `compute_att_gt.py`, `compute_att_gt2.py`  
**Tests:** `test_review_fixes.py::TestFactorCovariates`, `test_r_parity.py::test_factor_covariate_matches_r`

Categorical/transformed covariates (`Y~C(cat)`, `Y~C(cat)+x`, `poly(...)`, etc.)
were **silently dropped**: `preprocess` expands the formula into design columns
(`C(cat)[T.b]`, …) and drops the raw `cat` column, then the (g,t) loop **re-ran
patsy with the original formula** on data that no longer had `cat` → it failed
and fell back to an intercept-only design (numeric `Y~x` only worked because the
design column is also named `x`). Now the design is **built once** in
`preprocess` (factor-consistent, global levels) and the loop **selects those
columns** (`dp['xcov_cols']`) — matching R `did` (which builds `model.matrix()`
once) and removing `patsy` from the inner loop entirely (a further speed-up for
the standard path too).

**Validated vs R `did` 2.5.0** on a shared factor dataset: ATT(g,t) match to
**4.4e-15** in both standard and `faster_mode`.

---

### 15. Polish: idname validation, universal-base SE, faster_mode dedup

**Files:** `csdid/attgt_fnc/preprocess_did.py`, `csdid/att_gt.py`,
`csdid/attgt_fnc/compute_att_gt.py`, `compute_att_gt2.py`,
`compute_att_gt_shared.py` (new)  
**Tests:** `test_review_fixes.py` (`TestIdnameNumericValidation`,
`TestUniversalBaseNaN`, `TestPostKeyAlias`), `test_r_parity.py` (strengthened)

- **`idname` numeric validation (R parity).** R `did` v2.5.1 requires the id
  variable to be numeric (`stop("The id variable '…' must be numeric…")`). csdid
  now validates this in `_validate_inputs` with the same message intent instead
  of failing obscurely later. Removed the stale `#todo: idname must be numeric`.
- **Universal base-period SE → `NaN` (R parity).** For `base_period="universal"`
  the base cell has `att = 0` by construction and an all-zero influence function,
  so its SE is undefined. R reports `NA`; csdid previously reported a misleading
  `0`. It now reports `NaN`, matching R's `ref_gaps` reference exactly (e.g.
  `universal, g=2004, t=2003 → att=0, se=NA`). `test_r_parity` was tightened to
  assert csdid is non-finite wherever R is `NA`.
- **`faster_mode` de-duplication (maintainability).** The estimator dispatch and
  the per-(g,t) pretreatment-period / control-flow logic were copy-pasted across
  the standard and `faster_mode` paths. They now live once in
  `compute_att_gt_shared.py` (`select_estimators`, `last_pretreatment_index`,
  `plan_cell`) and are consumed by both, so the cell-selection logic cannot
  drift. Behavior is unchanged — the bit-identical consistency suite (72 configs)
  and all R-parity tests still pass.

---

### 16. Final review (GPT-5.5) follow-ups

**Files:** `csdid/att_gt.py`, `csdid/attgt_fnc/preprocess_did.py`  
**Tests:** `test_review_fixes.py` (`TestUniversalBaseNaN`, `TestIdnameNumericValidation`, +3)

A final GPT-5.5 code review of §15 surfaced two real edge-case defects, both fixed
and verified:

- **Universal base cell mis-identified by all-zero IF (Medium).** §15 detected the
  universal base cell via `np.all(inffunc == 0, axis=1)`. A *degenerate but valid*
  estimated cell (e.g. a perfectly deterministic outcome with a true zero effect)
  also has an all-zero influence function and was wrongly assigned `SE=NaN` instead
  of `0`. The base cell is now identified by **position** — `year` equal to the
  cohort's last pre-treatment period (`last_pretreatment_index`) — exactly R's
  definition, independent of IF contents.
- **Numeric dtype check crashed on pandas nullable ints (Medium).** The
  `np.issubdtype(data[col].dtype, np.number)` checks raised
  `TypeError: Cannot interpret 'Int64Dtype()'` for pandas nullable integer columns
  (`Int64`), and pandas extension dtypes generally. Replaced with
  `pandas.api.types.is_numeric_dtype` (excluding `is_bool_dtype` for
  `tname`/`gname`/`idname`), so nullable numeric ids/keys are accepted and boolean
  ones are still rejected (matching R). String ids remain rejected for R parity.

---

## Pre-existing Python bug fixes

Defects in the Python port itself (crashes / wrong output, independent of the R
update) are documented separately in **`PYTHON_BUGFIXES.md`**: `plot_attgt`
single-group crash, panel post-treatment flag always 0, panel global-`n`
overwrite, the broken `allow_unbalanced_panel=False` branch, and `raise "<str>"`
masking error messages.

---

## Remaining Work (Future PRs)

> All headline R `did` v2.5.1 items are implemented: input validation,
> `notyettreated`, clustered inference (bootstrap + analytical), `compute_inffunc`,
> `fix_weights`, `faster_mode`, and **factor/transformed covariates**. Possible
> future work: 2-way clustering, additional `est_method` plug-ins, trimming
> options.

---

## Test suite

**`csdid/test_csdid`: 564 passed, 0 skipped, 0 failures.** 65 tests compare
directly to live/published R values (see `test_csdid/test_suite_desc.md`);
`faster_mode` details in `test_csdid/FASTER_MODE.md`.

## Release

- `csdid/_version.py` bumped `0.2.9` → **`0.3.0`**.
- `setup.py` now sets `long_description` (from `readme.md`) + `long_description_content_type`
  so PyPI renders the README; `MANIFEST.in` ships the changelog/license and prunes
  test fixtures from the sdist.

## Running Tests

```bash
# Full R-translated suite (no R needed; R-parity tests skip if references absent)
pytest csdid/test_csdid

# Regenerate R reference values (requires R with did >= 2.5.0):
Rscript csdid/test_csdid/r_ref/generate_reference.R csdid/test_csdid/r_ref csdid/test_csdid
Rscript csdid/test_csdid/r_ref/generate_sim_reference.R
Rscript csdid/test_csdid/r_ref/generate_fixweights_reference.R csdid/test_csdid/r_ref
```

