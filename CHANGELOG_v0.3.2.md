# csdid v0.3.2

Fixes from a differential-testing audit of `csdid` against R `did` 2.5.0 / DRDID 1.3.0. This is the
public successor to v0.3.1; it consolidates the internal `0.3.2`–`0.3.5` work milestones and the v3
coverage round into one release.

## Audit structure

The audit is a **stress test**: `csdid` — the Python port of R's `did` — and R `did` 2.5.0 / DRDID 1.3.0
are run on identical scenarios, and their outputs are compared at floating-point tolerance (≈1e-9). Two
axes — the scenarios fed in (**A, stimulus**) and what is checked in the result (**B, examination**).
Across more than 130 scenarios, each compared on every output channel, nine defects were identified and
corrected, and R's parallel-trends pre-test — previously not implemented in `csdid` — was added to
bit-for-bit parity.

### A — Stimulus

A scenario is a choice on each of several axes of variation. The stress test sweeps each axis across its
*complete* set of categories, then crosses them — so the coverage argument is enumeration and partition,
not a sample.

- **Estimator settings** — a closed, documented set; every value is exercised:
  - control group — never-treated *and* not-yet-treated;
  - base period — varying *and* universal;
  - estimation method — all three: doubly-robust, IPW, outcome regression;
  - aggregation — all four: simple, event-study, group, calendar;
  - anticipation — none, one period, more than one;
  - weighting — unweighted and time-varying (every `fix_weights` mode).
- **Data structure** — a partition of the structural possibilities; one case in each category:
  - panel type — balanced, unbalanced, repeated cross-section;
  - covariate form — none, single continuous, multiple, categorical, transformed;
  - cohort / boundary structure — standard, tiny cohort, last-period-only cohort, no never-treated group,
    missing values.
- **Numerical conditioning** — how hard the underlying outcome regressions and propensity-score fits are
  to compute. This is a property of the data *values*, not the function arguments, so it is invisible from
  the API; it is empirically where difference-in-differences implementations diverge (Cunningham, 2026),
  and it is layered on top of any setting × data-structure case:
  - well-conditioned (benign);
  - scale disparity — covariates on very different scales, giving ill-conditioned design matrices;
  - near-separation — a covariate that almost perfectly predicts treatment, driving propensity scores
    toward 0/1 and the implied weights toward instability;
  - both combined.

### B — Examination

Every deterministic output is compared, not only the point estimate:

- per-cell and aggregated ATT(g,t);
- **all** standard errors — analytical, clustered, and bootstrap (under shared draws);
- t-statistics, uniform confidence bands, and the parallel-trends Wald pre-test (`W` / `Wpval`);
- the influence-function matrices — the object every standard error derives from, so a standard-error
  defect is caught at its source;
- NA / error status.

---

## Fixes

### Standard-error correctness (well-conditioned data)

**1. Per-cell analytical SE — population, not sample variance** · `att_gt.py` · High
- **Symptom:** every i.i.d. analytical SE inflated by `sqrt(N/(N-1))` (≈0.05–0.13%; larger as N falls).
- **Cause:** computed as `np.std(inffunc, axis=1, ddof=1)/sqrt(n)` (sample variance); the CS / DRDID form
  is the population `sqrt(mean(IF^2))/sqrt(n)`. Only the i.i.d. path — the clustered and aggregation paths
  already used the population form.
- **Fix:** population form.

**2. By-group aggregation overall SE — incorrect group index** · `aggte_fnc/compute_aggte.py` · High
- **Symptom:** by-group overall SE understated ≈0.06–0.5% (anti-conservative); all estimation methods,
  including no covariates; other aggregations unaffected.
- **Cause:** the per-(g,t) `group` vector (e.g. `[2,2,2,2]`) passed to the weight correction `wif()`
  instead of the distinct cohort list (`[2,3,4,5]`), counting the first cohort repeatedly.
- **Fix:** pass `glist`.

**3. `fix_weights="varying"` — influence function understated 2×** · `compute_att_gt.py`,
`compute_att_gt2.py` · High
- **Symptom:** point estimate correct, but every SE / t-stat / band / aggregation SE understated by a
  factor of two — undetectable by a point-estimate comparison. `base_period`/`first_period` already correct.
- **Cause:** in the RC routing of varying weights, the `n/n1` rescaling used `n1 = #rows` (2× units)
  instead of `#units`.
- **Fix:** rescale by unique units (a no-op for genuine repeated cross-sections).

**10. `fix_weights="varying"` standard errors over-inflated on genuinely unbalanced panels** ·
`compute_att_gt.py`, `compute_att_gt2.py` · High
- **Symptom:** the #3 fix (scale by unit count) was right for *balanced* panels routed through the RC
  estimator, but was also applied to *genuinely unbalanced* panels — where R stays in repeated-cross-
  section form and scales by row count — over-inflating every SE / band / aggregation SE by
  `n_rows/n_units`. Point estimates unaffected.
- **Cause:** pre-processing forces `panel=False` for `fix_weights="varying"` *before* its
  balanced/unbalanced check, so the unit-count scaling was applied regardless of actual balance.
- **Fix:** gate the unit-count scaling on actual balance (`nrow == n_units·n_periods`); otherwise use the
  row count — matching R on both balanced and unbalanced panels.

### Robustness — NA where R returns NA

R returns NA (with a warning) for a 2×2 cell that is not reliably estimable; `csdid` returned a finite
value. Now matched exactly.

**4. Singular / ill-conditioned covariate designs** · `compute_att_gt*.py`, `compute_att_gt_shared.py` · Robustness
- **Fix:** added R's guard `rcond(crossprod(control_covs)) < eps` (LAPACK 1-norm) on the control design
  for `dr`/`reg`, both compute paths.

**5. Under-identified treated cohorts** · `compute_att_gt*.py`, `compute_att_gt_shared.py` · Robustness
- **Symptom:** a treated cohort too small for the design (e.g. 3 post-treatment obs vs a 4-parameter
  design) makes the treated regression singular; `csdid` checked only the control design.
- **Cause:** DRDID's locally efficient DR fits outcome regressions for both groups.
- **Fix:** mirror DRDID's per-group, path-dependent checks — control only (balanced panel); control and
  treated, pre and post (RC / unbalanced); plus the propensity-score design.

**6. Propensity-overlap violations** · `compute_att_gt*.py`, `compute_att_gt_shared.py` · Robustness
- **Symptom:** where fitted PS ≥ 0.999 R returns NA for the `dr`/`ipw` cell; `csdid` returned a value
  from diverging IPW weights.
- **Fix:** implement R's `overlap_check_fail` (per-cell logit; NA + warning at 0.999), both paths, before
  the singular guard. Well-conditioned cells unaffected.

### Point-estimate correctness

**8. `notyettreated` × anticipation, no never-treated group — incorrect ATT** · `attgt_fnc/preprocess_did.py` · High
- **Symptom:** ATT(g,t) incorrect by up to ≈0.1; exact at anticipation 0; divergent only when the
  not-yet-treated set is the binding comparison group.
- **Cause:** *not* the control cutoff (`gname > tlist[max(t_i,pret)+tfac]+anticipation`, identical to R) —
  upstream in pre-processing: with no never-treated group, R retains only cohorts in `c(0, glist)`, but
  `csdid` dropped only first-period-treated rows, leaving the last-treated cohort in as a comparison unit.
- **Fix:** restore R's cohort-membership filter in shared pre-processing — equivalent whenever a
  never-treated group is present (other cases unchanged), applied once so standard and `faster_mode`
  cannot diverge.

### Packaging

**7. Undeclared runtime dependencies** · `setup.py` · Packaging
- `drdid` and `statsmodels` are imported at runtime but absent from `install_requires`; a fresh
  `pip install csdid` could fail to import. **Fix:** both declared.

### Feature parity (newly implemented)

**9. Parallel-trends Wald pre-test** · `att_gt.py` · Feature parity
- **Was missing:** R `did` reports a Wald pre-test that the pre-treatment ATT(g,t) are jointly zero
  (`W`, `Wpval`); `csdid` computed neither, so the diagnostic was unavailable and could not be checked.
- **Added:** `csdid` now computes both on the `MP` object — `W = n · preatt' solve(preV) preatt`,
  `Wpval = 1 − chi2.cdf(W, q)` — from the same analytic influence-function covariance the SEs use, with
  all of R's guard semantics (no usable pre-cells / NA or singular pre-covariance / beyond-unit
  clustering with a clustered bootstrap → not computed).

---

## Differential finding — a difference in R, not a `csdid` defect

**Not-yet-treated controls with anticipation = 2 and no never-treated group.** In this corner almost
every ATT(g,t) is not estimable — only a couple of pre-treatment periods remain. `csdid`'s ATT(g,t)
**match R cell-for-cell** (confirming the #8 fix). When R then aggregates, it hits an *acknowledged bug
inside R itself* and errors out; `csdid` instead returns the partial pre-treatment event-study (the only
estimable cells), which is the more useful result. Matching R's crash would be wrong, so `csdid`'s
behavior is kept. Recorded for transparency — a difference in R, not a `csdid` defect.
