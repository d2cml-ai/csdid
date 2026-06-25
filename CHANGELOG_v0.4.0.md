# csdid v0.4.0

Public successor to v0.3.1, consolidating the internal `0.3.2`–`0.3.5` work milestones, the v3
coverage round, the v4 boundary/refusal-parity round, and the v4.5 full-channel re-read into one
release. Fixes from a differential-testing audit of `csdid` against R `did` 2.5.0 / DRDID 1.3.0.

## Summary of changes

This changeset corrects the numeric, robustness, and point-estimate defects below and adds
input-validation and diagnostic parity with R `did`. It changes outputs relative to the current port,
so it is proposed as a minor (`0.4.0`) rather than a patch release:

- **Standard errors** — the per-cell analytic SE no longer carries the spurious `sqrt(N/(N−1))`
  inflation (Fix 1).
- **Confidence bands** — the pointwise critical value follows `alp` (`qnorm(1−alp/2)`) instead of a
  hardcoded `1.96` (Fix V4-V1); significance stars can change at non-default `alp`.
- **Point estimates** — the `notyettreated × anticipation>0 × no-never-treated` corner is corrected
  (Fix 8).
- **Input validation** — inputs R rejects that the port had silently accepted or crashed on now raise
  the same class of error: time-varying cluster var, `fix_weights="varying"` + callable `est_method`,
  non-logical `compute_inffunc`, all-NA aggregation, all-rows-dropped, multiple cluster vars, missing
  formula column, missing `idname` (Fixes V4-S1/S3/S4, V4-U1–U5).
- **Diagnostics** — the parallel-trends Wald pre-test (`W`/`Wpval`, Fix 9), plus a warning when no
  pre-treatment cells are available (Fix V4-V2). The pre-test is now computed over the valid
  pre-treatment cells even when some are not estimable (Fix V4.5-1), and the time-varying-weights
  warning fires regardless of `fix_weights` (Fix V4.5-2).

---

## Audit structure

The audit is a **stress test**: `csdid` — the Python port of R's `did` — and R `did` 2.5.0 / DRDID 1.3.0
are run on identical scenarios, and their outputs are compared at floating-point tolerance (≈1e-9). Two
axes — the scenarios fed in (**A, stimulus**) and what is checked in the result (**B, examination**).
Across more than 130 scenarios, each compared on every output channel, the numeric, robustness, and
point-estimate defects below were identified and corrected, and R's parallel-trends pre-test —
previously not implemented in `csdid` — was added to bit-for-bit parity.

A final round then added a **new stimulus dimension and a new examination channel** — *boundary /
refusal parity* — orthogonal to the three valid-scenario axes below. Where those axes vary a
well-formed scenario that produces output, this dimension drives the API to its edges: every argument
taken to its domain extremes and invalid values, and inputs constructed to trip each of R's documented
rejection guards. The matching examination is not "did `csdid` also error" but whether it refuses **the
same inputs** R refuses, **for the same reason** — the same *class* of message. This round surfaced ten
further divergences: inputs R rejects with a clear error that `csdid` had silently accepted or crashed
on, a confidence-band defect, and a missing diagnostic warning (one accepted as a deliberate deviation;
see below).

A further round then re-read the **valid-scenario corpus through the full examination channel suite** —
the per-cell critical value, pointwise band, significance stars, the Wald pre-test (`W`/`Wpval`), the
aggregation influence functions, and the warning/message text — channels the earlier valid-scenario
rounds had captured but never compared. With the confidence-band defect already fixed, the critical
value matched R on every case; this round instead surfaced two further divergences on the
diagnostic/communication surface: the Wald pre-test was silently skipped for an entire fit whenever a
single pre-treatment cell was not estimable (where R computes it over the rest), and R's
time-varying-weights warning was suppressed whenever `fix_weights` was set.

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
- **Argument domain & input validity** *(boundary / refusal dimension)* — orthogonal to the three axes
  above, which all assume a well-formed scenario that returns a result. This axis instead targets the
  edges of the API and the inputs R is built to *reject*:
  - argument domain extremes & invalid values — `alp` across `(0,1)` and out of range; non-logical
    `compute_inffunc`; a callable `est_method` crossed with `fix_weights="varying"`; a *list* of cluster
    variables; `panel=True` with no `idname`;
  - inputs that should trip R's documented rejection guards — time-varying cluster variable, more than
    one cluster variable, a formula naming a column not in the data, data that drops to empty, an
    aggregation in which every ATT(g,t) is NA.

  Coverage here is by **enumeration of the guards** (every documented rejection guard given a case —
  33/33 covered) and **partition of each argument's domain**, again not a sample.

### B — Examination

Every deterministic output is compared, not only the point estimate:

- per-cell and aggregated ATT(g,t);
- **all** standard errors — analytical, clustered, and bootstrap (under shared draws);
- t-statistics, uniform confidence bands, and the parallel-trends Wald pre-test (`W` / `Wpval`);
- the influence-function matrices — the object every standard error derives from, so a standard-error
  defect is caught at its source;
- NA / error status;
- **refusal parity** *(boundary dimension)* — whether `csdid` refuses the inputs R refuses, and does so
  with the same **class of message** ("same reason", not merely "both errored"), so a port that crashes
  with an opaque exception where R returns a clean, classified error is flagged;
- **channel behaviour at argument edges** *(boundary dimension)* — output channels checked as the
  arguments move, not only at defaults: e.g. the pointwise critical value as `alp` varies, and the
  presence of the pre-test warning when no pre-treatment cells are available.

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

### Confidence-band correctness

**V4-V1. Pointwise critical value ignored `alp`** · `att_gt.py` · High
- **Symptom:** the per-cell critical value `MP$c` — and the pointwise band and significance stars
  derived from it — was a hardcoded `1.96` on the analytic-SE path, ignoring the `alp` argument. Wrong
  at any `alp ≠ 0.05` (e.g. `alp=0.10` → R's `qnorm(0.95)=1.6449` vs `csdid`'s `1.96`; stars can flip),
  and off by ≈3.6e-5 even at the default `0.05`.
- **Scope:** confined to the analytic `att_gt` band — the bootstrap critical value and the aggregation
  critical value (`crit_val_egt`) already followed `alp` and matched R.
- **Fix:** `crit_val = qnorm(1 − alp/2)` (`scipy.stats.norm.ppf`).

### Input validation & error parity

R rejects a number of ill-formed inputs up front with a clear, classified error; `csdid` either
accepted them silently (a wrong or misleading result) or crashed with an opaque low-level exception.
Each now raises the **same class of error as R, at the same point**.

| id | Input | Before (`csdid`) | Now |
|----|-------|------------------|-----|
| V4-S1 | time-varying cluster variable | silently used the first period's labels | refuses ("time-varying cluster variables are not supported"), in pre-processing |
| V4-S3 | `fix_weights="varying"` + callable `est_method` (panel) | invoked the callable with the wrong (RC) signature | refuses before estimation |
| V4-S4 | non-logical `compute_inffunc` (e.g. `"yes"`) | `bool()`-coerced (any truthy value → True) | refuses ("must be a single logical") |
| V4-U1 | aggregation when every ATT(g,t) is NA (`na_rm=True`) | crashed: `max() arg is an empty sequence` | refuses ("all att_gt() estimates are NA") |
| V4-U2 | a model column entirely missing → all rows dropped | crashed: zero-size-array reduction | refuses ("all observations dropped due to missing data") |
| V4-U3 | two cluster variables (a list) | crashed: `unhashable type: 'list'` | refuses ("at most one cluster variable"); a single-element list is unwrapped |
| V4-U4 | a formula naming a column not in the data | patsy `NameError` ("error processing formula …") | reports the missing variable ("… not in data: …") |
| V4-U5 | `panel=True` with no `idname` | "column 'None' not found in data" | "must provide idname when panel=True", before the column check |

Fixes live in `attgt_fnc/preprocess_did.py` (S1, U2–U5), `att_gt.py` (S3, S4), and
`aggte_fnc/compute_aggte.py` (U1).

**V4-V2. Missing pre-test warning** · `att_gt.py` · Diagnostic parity
- **Symptom:** when no usable pre-treatment cells exist (e.g. a two-period panel), the Wald pre-test
  cannot be computed; R emits a warning, `csdid` returned `W=None` silently.
- **Fix:** warn (matching R) before returning, so the absent diagnostic is signalled, not silent.

### Diagnostic parity — full-channel re-read of the valid corpus

**V4.5-1. Wald pre-test silently skipped when some pre-treatment cells are not estimable** ·
`att_gt.py` · Medium
- **Symptom:** when *some* (not all) pre-treatment ATT(g,t) cells are NA — a single overlap/singular/
  near-separated cell, universal-base reference cells, or a not-yet-treated × anticipation corner — the
  port returned `W=None`/`Wpval=None` for the **whole** fit, silently dropping the parallel-trends
  diagnostic. R instead drops the non-estimable pre cells and computes the Wald over the valid subset.
- **Cause:** the pre-cell filter `zero_na = se <= eps` did not catch NaN-variance cells (`NaN <= x` is
  `False`), so an NA cell survived the filter, left `NaN` in the pre-covariance, and forced the pre-test
  to bail — even when other pre cells were perfectly estimable. The same path also bypassed the V4-V2
  warning in the all-NA case (`pre.size` was non-zero because the NaN cells were retained).
- **Fix:** also drop non-finite-variance pre cells (`| ~np.isfinite(se)`). The port now reproduces R's
  `W`/`Wpval` over the valid pre cells (verified to ~1e-9), and the V4-V2 warning fires correctly when
  *every* pre cell is non-estimable.

**V4.5-2. Time-varying-weights warning suppressed when `fix_weights` is set** ·
`attgt_fnc/preprocess_did.py` · Low
- **Symptom:** with time-varying sampling weights on a panel, R emits "Time-varying weights detected … use
  the 'fix_weights' argument to control this behavior" *regardless* of `fix_weights`; the port emitted it
  only when `fix_weights` was unset, so a user who explicitly chose a `fix_weights` mode saw no warning.
- **Cause:** the warning was gated on `fix_weights is None`.
- **Fix:** warn whenever time-varying weights are detected on a panel, matching R.

---

## Differential finding — a difference in R, not a `csdid` defect

**Not-yet-treated controls with anticipation = 2 and no never-treated group.** In this corner almost
every ATT(g,t) is not estimable — only a couple of pre-treatment periods remain. `csdid`'s ATT(g,t)
**match R cell-for-cell** (confirming the #8 fix). When R then aggregates, it hits an *acknowledged bug
inside R itself* and errors out; `csdid` instead returns the partial pre-treatment event-study (the only
estimable cells), which is the more useful result. Matching R's crash would be wrong, so `csdid`'s
behavior is kept. Recorded for transparency — a difference in R, not a `csdid` defect.

**Reserved column name `.w` — deliberate deviation (audit id V4-S2).** R reserves dotted internal
names (`.w`, `.y0`, …) and refuses a user weights column named `.w`. `csdid`'s internal columns are
*non-dotted* (`w`, not `.w`), so a user `.w` column does not actually collide with anything internal;
`csdid` keeps accepting it. Accepted as a deliberate deviation rather than mirroring R's self-imposed
naming restriction — recorded for transparency.

---

Builds on **v0.3.1**; proposed version **0.4.0** — the changeset alters numeric outputs (standard
errors, confidence bands, some point estimates) and adds input-validation errors, so it is a minor
bump rather than a patch. Merged in #NN.
