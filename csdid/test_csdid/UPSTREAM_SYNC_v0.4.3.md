# R alignment and adversarial validation for 0.4.3

Follow-up: the [2026-10-24 test audit](R_TEST_ALIGNMENT_20261024.md) adds
assertion crosswalks, permanent pinned current-R fixtures and four further
correctness fixes, including final-score propensity information. Results and
solver residuals below describe the earlier source stage, not a rerun after
that correction. See the follow-up for current validation and remaining gaps.

The release targets the supported Python `dr`, `reg`, and stabilized `ipw`
estimators, their preprocessing, and aggregation. It also changes the Python
default to `faster_mode=True`; `False` remains supported. The objective is to
port relevant fixes with small runtime changes, preserving estimator and solver
choices while checking both computation paths independently.

## Controlling sources

The audit used current source checkouts rather than the older installed R
packages. GitHub and local source commits were checked on 2026-10-23:

- [`did` 2.5.1.902, `74b88eb`](https://github.com/bcallaway11/did/commit/74b88eb07f1faa644df1271055a0f77848f6550c).
- [`DRDID` 1.4.0, `eaac3b4`](https://github.com/pedrohcgs/DRDID/commit/eaac3b49cd225d0d7706c633858555a5f204e21e).

[`did` #273](https://github.com/bcallaway11/did/pull/273) and
[#275](https://github.com/bcallaway11/did/pull/275) supply the relevant
aggregation-window, repeated-cross-section ID, inference, and panel-balancing
behavior. [`DRDID` #20](https://github.com/pedrohcgs/DRDID/pull/20) corrects
pre-period terms in the traditional RC DR influence function. Retained-sample
weight validation was checked against current R preprocessing and
[`DRDID` #34](https://github.com/pedrohcgs/DRDID/pull/34). The RC regression
feasibility guard now checks only models actually fitted by the chosen estimator.
See the [release notes](../../CHANGELOG_v0.4.3.md) for user-visible consequences.

At the audited commits, `did` #276's anticipation changes belong to a separate
feature branch and are absent from master; #278 is open and concerns entry to
R's native `fastglm`. R dependency-only changes and DRDID's TWFE, nonstabilized
IPW, IPT/improved DR, wrapper, and weighted-refitting bootstrap changes do not
describe estimator routes dispatched by this Python package. They were not
ported merely to match unrelated source edits.

## Checks and results

- `pytest csdid/test_csdid/ -q --tb=short`: **1,428 passed**, no failures or
  skips, with 16 plotting/small-group warnings, on Python 3.12/macOS arm64.
- The original Python suite passed 703 tests. Alignment added 160 regressions;
  performance and default-mode checks bring the final total to 1,428. Existing
  genuine RC fixtures now omit IDs or give rows distinct IDs, retaining their
  mathematical assertions.
- Independent empirical derivatives and period reversal check the corrected
  RC DR influence function, including misspecified propensity models, unequal
  weights, zero weights, and active trimming. The expanded stress check compared
  1,920 observation derivatives (maximum absolute error `2.31e-6`).
- Independent ratio/IF calculations check aggregation weights, full influence
  vectors, and SEs. Original-code rejection checks and seven surgical mutations
  show that the new alignment tests detect the defects they target.
- A current-source R grid covers **144 fits, 1,944 ATT cells, and 2,304
  aggregation requests**, including 720 symmetric expected errors. It crosses
  both explicit `faster_mode` values, panel/unbalanced/RC data, regular/irregular
  dates, three estimators, two bases, and two controls. Weights, clustering,
  and anticipation rotate through these cases rather than forming another full
  factorial. Public analytic outputs pass the declared tolerances; maximum
  ATT/SE differences are `2.53e-11`/`1.02e-8`.
- A separate shared-multiplier grid covers **48 DR fits, 408 cells, and 288
  aggregation requests**, using 399 identical signed draws. It includes both
  modes, two/five periods, unequal clusters, and both confidence-band requests.
  Unit IDs align influence rows to draws before comparison.
- Complete before/after comparisons for all 144 analytic and 48 bootstrap fits
  are **exact** for public outputs, full influence functions, and warnings.
  Performance regression tests additionally exercise IDs beyond `2**53`, missing
  and nullable IDs, mixed integer/float dates, reordered frames, helper overrides,
  unusual direct helper inputs, bootstrap chunk boundaries, and random states.
- Wheel and source distribution build successfully as 0.4.3 and exclude tests
  and stress artifacts. An installed-wheel smoke check confirms the default,
  both computation paths, and aggregation. Performance measurement and replay
  instructions are in [FASTER_MODE.md](FASTER_MODE.md).

## Numerical limits and preserved failed routes

The R comparison uses `atol=1e-7, rtol=1e-6` (Wald statistics use `1e-5` for
both). With production solver defaults, full influence functions exceed the
stricter componentwise threshold in 96 analytic and 32 bootstrap case checks,
with maximum absolute differences `3.19e-5` and `2.95e-5`. Four bootstrap
critical-value fields in two equivalent clustered-RC cases also exceed it
(maximum `1.15e-5`); supports and finite masks match.

A separately labeled sensitivity run tightening GLM convergence to `1e-12`
makes all public analytic and bootstrap outputs pass, including those four
critical values. Full influence functions still have 14 analytic and 16
bootstrap case exceedances, with maxima `6.88e-6` and `2.85e-6`. These residuals
are present before the speedups. Production solver tolerances were preserved;
the evidence supports numerical implementation fidelity, not bitwise R equality
or complete R API/default equivalence.

A proposed Wald covariance shortcut was rejected after adversarial boundaries
changed guard decisions and p-values. Bootstrap chunk shapes, overlap checks,
and weighted versus unweighted propensity fits were retained. An initial R
bootstrap harness with incorrectly paired influence rows was corrected and its
failed results retained separately. The original lower-coverage R grids remain
historical evidence; the 144-fit comparison above supersedes them.

Detailed scripts, source hashes, arrays, logs, and rejected-route evidence are
preserved locally under the ignored `stress/upstream_sync_20261023/`,
`stress/performance_20261023/`, and `stress/release_v043_20261023/` directories.
Large arrays and archives are excluded from the PR and distribution. R replay
scripts retain recorded local checkout paths; restore those paths or adjust the
recorded root constants. Self-contained regression tests require no R runtime.
