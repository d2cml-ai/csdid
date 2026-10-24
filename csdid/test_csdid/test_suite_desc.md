# Python test coverage and R alignment

The canonical suite is `csdid/test_csdid`. Run it with `python -m pytest
csdid/test_csdid -q`; R is not required. GitHub CI runs the non-slow tests on
Python 3.8 with legacy dependencies and Python 3.12 with current dependencies.
The two slow tests are also included in the full local run.

The 2026-10-24 full run passed **1,953 tests**, with no failures or skips.

Python and R do **not** have one-to-one unit-test coverage. A test translated
from R may check behavior or simulated ground truth without comparing any R
numbers. The assertion-level crosswalks and current validation results are in
[the October 2026 audit](R_TEST_ALIGNMENT_20261024.md).

## Numerical reference tests

| Test module | Numerical workflows | Reference and scope |
|---|---:|---|
| `test_r_current_parity.py` | 50 | 25 pinned current-source R configurations in both Python modes; full ATT/SE/bands/Wald output, unit influences, and all four aggregations with full influences |
| `test_r_parity.py` | 36 | Historical R did 2.5.0 fixtures, including weights, factors, controls, bases, anticipation and clustering |
| `test_sim_parity.py` | 24 | Historical R did 2.5.0 fixtures: six simulated datasets, two controls and DR/regression |
| `test_jel_replication.py` | 5 | Published JEL replication values |

The current-source references pin did `74b88eb` (2.5.1.902) and DRDID `eaac3b4`
(1.4.0). They cover weighted balanced panels, true repeated cross sections and
unbalanced panels; DR/regression/IPW; two control/base combinations; varying
weights; analytical clustering; and two shared-draw bootstrap configurations.
The bootstrap tests use identical frozen multipliers, not matching RNG seeds.

Every numerical component uses `atol=1e-7, rtol=1e-6`, except Wald outputs
(`atol=rtol=1e-5`). Support, dimensions and nonfinite masks must agree exactly.
Separate corruption tests check that the comparison rejects malformed results;
a provenance test verifies source commits, artifact hashes and R mode checks.
These are additional tests, not additional statistical configurations.
See [fixture documentation](r_ref/current/README.md) for regeneration and limits.

The historical fixtures retain their original, broader tolerances, including
3–5% relative SE tolerances. They are useful regressions but do not establish
machine-precision agreement. Missing historical or current fixture files fail
the suite rather than silently skipping comparisons.

## Behavioral and independent numerical tests

The remaining tests cover public validation, sample selection, error and warning
behavior, estimator identities, aggregation formulas, influence derivatives,
bootstrap behavior, plotting, execution-mode parity and performance safeguards.
The October audit adds 114 did-alignment cases, 152 DRDID/aggregation cases and
12 final-score Hessian cases. Their independent oracles include perturbed
sampling-weight derivatives, closed-form aggregation calculations and explicitly
retained/releveled samples.

The [did crosswalk](r_test_alignment_did.json) reviews all 34 upstream test files:
9 covered, 22 partial, 2 not applicable and 1 gap. The
[DRDID crosswalk](r_test_alignment_drdid.json) reviews all 29 files: 12 covered,
3 partial and 14 not applicable. These are qualitative assertion-intent
classifications, not coverage percentages or claims that every R assertion was
ported. Each partial or excluded contract is described in the crosswalk.

R-only conditional CvM tests, improved/IPT estimators, TWFE, non-normalized IPW
and weighted-refit bootstrap are outside csdid's current supported estimator
routes. Formula-language and some diagnostic differences remain. Finite test
coverage cannot guarantee identical answers for every input or dependency.
