# csdid 0.4.2

Themed differential audit against R `did` **2.5.1** / DRDID 1.3.0 (10 mechanism-based themes
run to convergence on the pristine 0.4.1 baseline). This release folds in the 2.5.0→2.5.1 sync
work (N1 + `fix_weights`) and everything the themed audit surfaced. Channel-complete parity was
re-verified against R 2.5.1 via the differential comparator; all 703 packaged unit tests pass.

## ⚠️ Behavior changes — re-validate results

- **`allow_unbalanced_panel` now defaults to `False`** (was `True`), matching R `att_gt`. With
  the default and an **unbalanced** panel, the data is now coerced to a balanced panel (dropping
  units not observed in every period, with a warning) instead of silently routing to the
  repeated-cross-section estimators. **Point estimates, `n`, and SEs change** for any
  default-argument call on unbalanced data — pass `allow_unbalanced_panel=True` explicitly to keep
  the previous behavior. *(theme O2)*
- **`fix_weights="varying"` on a balanced panel:** the analytic SE was exactly **2× too large**
  (the RC influence function was not folded ÷2 per unit). SEs on this path are now correct —
  inference was previously twice too conservative. With **time-varying covariates** on this path
  the **point estimate** also changes: stacked rows are now conditioned on the unit's
  earlier-period covariate (R: "fix_weights changes weights, not the covariate conditioning set").
  *(themes O1, C1)*
- **No-never-treated + `control_group="notyettreated"` + a first-period-treated cohort (N1):**
  the first-period drop now removes only the already-treated **units** (by row identity), not every
  member of their **cohort**. The previous cohort-membership filter also deleted the latest cohort
  that was deliberately retained as a not-yet-treated control, corrupting `ATT(g,t)` for the other
  groups. Estimates on this path change (now correct). *(theme O3; reverses the earlier 2.5.0-matched
  port behavior, per R 2.5.1.)*
- **`aggte` corrected aggregates** *(theme O6)*:
  - `type="calendar", na.rm=True` with a calendar period whose post-treatment cells are all NA: the
    period is now dropped (R `has_post` guard) instead of being kept as a phantom `att=0.0` that
    silently biased the overall calendar ATT.
  - `type="group", na.rm=True` with a finite `max_e` and a group whose only non-NA post cell is past
    the window: that group is now dropped (R's `t <= group + max_e` filter), instead of yielding
    `overall_att = NaN`.
- **`alp=NaN` is now rejected** (`att_gt` and `aggte`). It previously slipped past the `0<alp<1`
  guard and produced all-NaN critical values / confidence bands. *(theme C1)*

## Graceful failures instead of crashes (now matching R's clean refusal)

- `biters=1` with `bstrap=True` no longer raises `IndexError` (degenerate single-draw bootstrap now
  degrades gracefully). *(themes O7, C1)*
- An empty event-window (`min_e>max_e`, out-of-range window, etc.) now raises a clear
  "No event times fall within the requested window…" instead of a numpy
  `need at least one array to concatenate`. *(themes O6, C3)*
- Aggregating when **every** `att(g,t)` is NA, or when cohort/anticipation trimming empties the
  data, now raises the clean R-style message instead of a `TypeError`/`IndexError`. *(themes O5, O3, C3)*
- Container-valued or column-aliasing name arguments (e.g. `gname=['x']`, `clustervars==idname`,
  `yname==gname`) are unwrapped/de-duplicated or refused cleanly instead of crashing with
  `unhashable type` / `AttributeError`. *(theme C2)*

## New input validation (matching R `did` 2.5.1's `validate_*` hardening)

- **Logical-scalar** validation for `panel, allow_unbalanced_panel, bstrap, cband, faster_mode,
  print_details, pl` (and `na_rm`): a mistyped flag (`panel=0`, `bstrap="TRUE"`, `(False,)`) is now
  refused instead of being silently truthy-coerced to the wrong meaning. *(theme C1)*
- **`cores`** must be a positive whole number (was silently `int()`-truncated / accepted 0/neg).
- **`anticipation`** must be a finite, non-negative, whole number (was accepting `0.5`/`NaN`/`Inf`/bool).
- **`alp` / `anticipation`** now also accept idiomatic numpy scalars and length-1 numeric vectors.
- **`aggte`-stage args** (`alp`, `min_e`, `max_e`, `balance_e`, `na_rm`) are validated like R.
- **`idname=None` is allowed when `panel=False`** (repeated cross-sections); a bool outcome column is
  cast to numeric; `clustervars=[idname, x]` drops `idname` before the one-cluster check; ±Inf rows in
  required columns are dropped (R `complete_finite_cases`) instead of kept. *(themes C1, C2)*

## Warnings / messages and argument aliases (matching R)

- `att_gt` / `aggte` now accept R's argument spellings **`clustervars`** (→ `clustervar`) and
  **`weightsname`** (→ `weights_name`); previously these were silently dropped to `**kwargs`,
  **disabling clustering / weighting without notice**. *(themes C2, C3)*
- Added the missing notices: the `anticipation > 0` note; the "`min_e`/`max_e`/`balance_e` ignored for
  `type="calendar"`" warning; the aggregate-cluster-unavailable fallback warning; and a warning when
  the Wald pre-test is suppressed due to a singular / NA pre-treatment covariance. *(themes C1, O3, O6, C3)*

## Known cross-implementation differences (documented, not fixed)

These are floating-point / linear-algebra differences between R's compiled LAPACK/DRDID and the
port's numpy/scipy, not logic bugs; they vanish or stay within numerical tolerance:

- dr/ipw analytic SE on covariate models differs at the ~1e-9–1e-7 level and **scales away with `n`**
  (DRDID-vs-numpy nuisance-term noise). *(themes O1, O2, O4, O5)*
- The near-singular `rcond` keep/drop boundary (port exact `1/cond` vs R LAPACK `dgecon`) can differ
  only within ~3× machine-eps; near-collinear/affine designs and intercept-suppressed (`~0+x`)
  outcome regressions differ structurally on degenerate input. *(themes O4, O5)*

## Testing

- 703 packaged unit tests pass; channel-complete parity re-verified against R `did` 2.5.1 on the
  themed audit's case corpora (the result-changing fixes above flip from divergent to matching;
  faithful paths — RCS, explicit unbalanced routing, factors, base-period — remain identical).
</content>
