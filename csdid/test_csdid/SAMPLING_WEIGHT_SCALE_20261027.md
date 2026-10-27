# Cross-package sampling-weight audit

The goal is equivalent estimates and inference for equivalent specifications in
R DRDID, R did, Python DRDID and Python CSDID. Matching point estimates alone is
insufficient: the complete influence function and inference must agree too.

## Where changes belong

The standalone Python distribution is `drdid`, maintained in
[DRDIDpy](https://github.com/d2cml-ai/DRDIDpy). CSDID currently imports its
regression kernels, while its traditional DR and normalized IPW kernels live in
`attgt_fnc/drdid_trim.py`, with additional cached-fit support. Therefore a fix in
only one repository cannot certify all consumers.

[DRDIDpy PR #21](https://github.com/d2cml-ai/DRDIDpy/pull/21) corrects the six
existing standalone kernels and adds numerical R and downstream CSDID CI. This
CSDID correction works with the currently published `drdid==1.1.6` and with the
proposed upstream implementation. A published upstream release is still needed
before raising CSDID's minimum dependency or consolidating its local kernels;
that consolidation must retain cached-fit behavior and the group-time tests.

Compare the same sample, outcome, design, weights, estimator and inference
settings. R DRDID's high-level default is improved/IPT; R did and CSDID's built-in
DR route use the traditional estimator. Improved/IPT, alternative RC estimators,
TWFE, nonnormalized IPW and R's weighted-refit bootstrap are outside these six
shared Python kernels.

## Failures found and corrected

An ordinary mean reduction can overflow when the true mean of finite positive
sampling weights is finite. All 12 public method/design/mode combinations could
then return NA estimates under a harmless common weight scale of `1e307`.
Scale-safe means now retain the ordinary arithmetic when it succeeds and use a
bounded reduction only when needed. The selected two-period weights are repaired
before invoking the existing dependency. Raw user columns, unrelated cohorts and
custom-estimator inputs are preserved.

Aggregation also formed a squared cohort-mass denominator. Large scales remove
the denominator derivative through overflow; small scales yield nonfinite or
inaccurate derivatives through underflow and subnormal rounding. The corrected
fallback computes the same ratio derivative with dimensionless selected-cohort
shares. Selected mass totals and averages also use bounded arithmetic when
necessary. Ordinary-scale calculations retain their existing path.

The same aggregation defect exists in R did's slow path; its faster path
normalizes weights earlier. [R did PR #280](https://github.com/bcallaway11/did/pull/280)
corrects that source. For a single cohort and treatment cell, all four
aggregations must reproduce the cell's IF. In a shared example the old aggregate
SE was `0.1909288`, versus the correct cell/aggregate SE `0.1630024`. This is a
behavioral correction supported by an analytic nesting identity and independent
empirical derivatives.

## Checks and limits

The full Python suite passes 2,469 tests, including 159 new sampling-weight tests
and the two slow tests. The same new tests fail 120 times on the baseline, with
39 passing. The tests cover all public methods, panels/RC, both modes and all
four aggregations, complete IF and SE, scales from `1e-200` to `1e307`, float32,
duplicate cohort keys, zero mass, and local support in the presence of an
unrelated large-weight cohort. They also check input immutability.

The standalone proposed wheel passes 379 tests, including 80 full R numerical
configurations and independent formula/derivative/transformation contracts. Its
four upstream CI jobs pass, including downstream CSDID integration.

The live common-input audit compares 120 Python fits with 120 corrected R fits
across single- and multiple-cohort panel/RC designs, all three methods, both
modes and five weight scales. All 552 complete-payload cross-language, mode and
unit checks pass with unchanged tolerances. It checks row identity explicitly.
An independent before/after audit finds 240 healthy workflows and 960 aggregates
bit-identical, including all numeric fields, warnings and RNG state. Replay
scripts and compact records are in [audit20261027](audit20261027/).

Reviewed ordinary R fixtures, source pins, tolerances, workflow triggers and
required checks remain unchanged. The newly discovered R slow-path discrepancy
is reported explicitly; it is not concealed by replacing those fixtures. After
the R companion merges, the existing upstream drift gate requires review of its
new tests and source pin before accepting that new default-branch baseline.

The R upstream gate runs on Python pushes. An R-only push does not itself
trigger that workflow, and passing a numerical grid never ports an upstream change
automatically. Cross-package equivalence applies to supported, matched
specifications and the numerical checks recorded here, not to every R API or
arbitrary unrepresentable floating-point input.
