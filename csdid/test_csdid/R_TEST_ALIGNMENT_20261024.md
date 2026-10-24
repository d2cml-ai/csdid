# R test alignment and permanent numerical references

Audit completed 2026-10-24, following performance commit `dab3876`. The question
was whether Python's tests reproduce the applicable R assertions and whether
current R numerical output is checked permanently. They did not provide full
coverage. This change adds explicit coverage records, strict committed R
references and regression tests, and fixes four defects revealed by those tests.

## Scope and acceptance rules

The target is csdid's existing DR, normalized IPW and regression estimators,
preprocessing, inference and four aggregation types, in both execution modes.
The audit reads assertion intent rather than treating similarly named files or
test counts as evidence of equivalence. A reference comparison must retain all
cells, unit identities and influence components. Missing results, mismatched
nonfinite masks and absent fixtures fail; numerical disagreement is investigated
without changing data, removing cases or widening tolerances.

The controlling sources, verified as upstream heads during this audit, are
[did `74b88eb`](https://github.com/bcallaway11/did/commit/74b88eb07f1faa644df1271055a0f77848f6550c)
(2.5.1.902) and
[DRDID `eaac3b4`](https://github.com/pedrohcgs/DRDID/commit/eaac3b49cd225d0d7706c633858555a5f204e21e)
(1.4.0). “Current” means these pinned revisions, not an automatically moving
upstream target. The fixture manifest records clean source trees, source and
artifact SHA256 hashes, loaded package versions and the R environment.

## Assertion crosswalk

| Upstream suite | Files reviewed | Covered | Partial | Not applicable | Gap |
|---|---:|---:|---:|---:|---:|
| did | 34 | 9 | 22 | 2 | 1 |
| DRDID | 29 | 12 | 3 | 14 | 0 |

These are qualitative file classifications, not percentages of R assertions
ported. [did](r_test_alignment_did.json) and
[DRDID](r_test_alignment_drdid.json) records identify source assertions, actual
Python tests and remaining distinctions. Additions include 114 did behavioral
cases, 152 DRDID/aggregation cases and 12 final-score inference cases.
Independent oracles differentiate refitted estimators after perturbing sampling
mass, calculate weighted aggregation ratios and influences, and compare against
explicit retained/releveled samples.

The earlier numerical references use R did 2.5.0, with some relative SE
tolerances of 3–5%. They remain useful historical tests, but do not establish
machine-precision agreement with current R. Their 60 tests now require the
committed fixtures and exact cell support and compare missing/nonfinite values
instead of silently skipping them. Earlier coverage prose claiming full
assertion translation or machine precision has been replaced.

## Current R numerical gate

There are 25 configurations, each evaluated in both Python modes: 23 analytical
cases and two cases with frozen shared bootstrap multipliers. They cover
weighted balanced panels, repeated cross sections and unbalanced panels;
DR/regression/IPW; never-treated/varying-base and not-yet-treated/universal-base
settings; weight conventions; and analytical clustering. All four aggregations
are checked for every configuration.

Comparisons include ATT, SE, critical values, interval endpoints, Wald outputs,
full unit influence matrices and every aggregate influence component. Public
outputs and influence functions use componentwise `atol=1e-7, rtol=1e-6`; Wald
outputs use `atol=rtol=1e-5`. Dimensions, support and NaN/infinity masks agree
exactly. Fields are validated before permutations so extra rows or columns
cannot be discarded by alignment. Corrupted-output tests independently check
the comparison gate, and provenance hashes prevent silent fixture replacement.

R's standard and fast paths agree on support/masks for all 25 cases; maximum
absolute difference is `2.416e-13`. Normal regeneration reuses the committed
input and multiplier CSVs and rejects missing inputs. Explicit
`--regenerate-inputs` is required to regenerate data. Replay preserved all
25 reference payloads and frozen inputs exactly. See
[fixture documentation](r_ref/current/README.md) for commands and schema.

The two bootstrap cases test identical unit/cluster draws and downstream SE/band
calculations. They do not equate R and NumPy random seeds, and 99 draws are a
test configuration rather than an empirical-analysis recommendation.

## Defects found and corrected

1. **Unused categorical levels.** Globally unused levels, including levels lost
   through missing-data exclusion or balancing, could leave singular designs.
   Both modes returned zero finite ATT cells out of ten in a fixture whose
   explicitly releveled data yielded ten. Rebuilding categorical designs on
   the retained sample now reproduces the cleaned sample's ATT, SE and full IF,
   without modifying the caller's data.
2. **Missing custom influences.** Per-unit folding replaced an observed NaN
   influence with zero. The regression injected one missing contribution per
   cell and previously obtained finite SEs throughout. Missing observed
   contributions now propagate; units absent from a cell retain structural
   zero, and healthy contributions remain unchanged.
3. **Balanced panels with varying weights.** Routing estimation through RC
   kernels inadvertently switched aggregation to average-period cohort weights.
   R retains the balanced-panel sampling design and uses first-period cohort
   weights. Eighteen independent ATT/full-IF/SE oracles cover the correction and
   true-RC/unbalanced boundaries. Four affected current-R comparisons now pass.
4. **Final propensity information.** Statsmodels can retain preceding-iteration
   covariance weights when IRLS returns its final fitted probabilities. R
   DRDID explicitly rebuilds the propensity information from final capped
   probabilities and normalized weights. All four owned DR/IPW kernels now do
   this, with the existing condition cutoff and a Cholesky inverse. GLM calls,
   default convergence tolerance, fitted coefficients, ATT arithmetic, trimming
   and propensity-cache behavior are unchanged. IFs and SEs intentionally change.

The last discrepancy initially looked like ordinary cross-solver precision.
Diagnosis of six failing IPW cases found 287 componentwise violations across
24 fields, with maximum IF difference `3.3513e-6`. Recomputing only final-score
information removed every violation and reduced the maximum to `3.1144e-8`.
Tightening GLM tolerance was examined but not adopted: fitted coefficients
already differed by at most `2.59e-12`. The compact
[diagnostic evidence](r_test_alignment_results_20261024.json) retains these
outcomes. Frozen data, expected outputs and numerical tolerances were unchanged.

## Validation and limits

The final complete suite passes **1,953 tests** in 45.84 seconds, with 25 warnings
and no failures or skips, on Python 3.12/macOS arm64. This adds 372 tests to the
performance-stage suite. The 94 current-reference tests comprise 50 numerical
workflows, 43 corruption/support rejection cases and one provenance check.
Independent review found and closed a comparison weakness: shape validation
must precede key-based reordering, which otherwise can discard excess rows or
columns. Thirty-four of the rejection cases exercise that full comparison.

Fresh wheel and source archives exclude tests, stress outputs and R references;
all 18 installed source modules match the checkout. An isolated installed-wheel
smoke outside the source tree ran six DR/IPW/regression fits and 24 aggregations.
All 108 compared numeric fields, including full IFs, match exactly across modes,
with no runtime warnings. The pre-existing declared `twine` dependency is absent
from the local smoke environment; all eight runtime dependencies are present.
The companion [evidence file](r_test_alignment_results_20261024.json) records the
environment, runtime hashes, packaging results and numerical diagnosis. GitHub
CI on the PR independently runs the non-slow suite on Python 3.8 and 3.12/Linux.

This is not a claim of complete R API parity or risk-free universal equality.
R-only conditional CvM tests, improved/IPT estimators, TWFE, non-normalized IPW
and exponential-weight nuisance-refit bootstrap are outside Python's dispatched
routes. Remaining partial coverage includes diagnostics, formula semantics,
some invariance variants and native multiway-bootstrap oracles. The current R
fixtures use numeric `~X` and zero anticipation; separate historical/behavioral
tests cover other formulas and anticipation. Captured warning text is not part
of the numerical fixture comparison. The standalone
unbalanced-panel diagnostic-message file remains a named gap. A categorical
predictor with only one observed level can still be accepted by Patsy where R
rejects its contrasts; this is an explicit formula-language difference.
Condition decisions near machine precision can differ across numerical
libraries. These cases do not justify relaxing the committed reference gate.

The earlier [performance audit](PERFORMANCE_AUDIT_20261023.md) measured
optimizations through `dab3876`; its source hashes, timings and exact-output
claims refer to that stage. The four fixes here intentionally correct results
and must not be described as preserving every old Python number. Broader earlier
R grids are historical evidence, not reruns of this final change.

A serial follow-up measures the inference correction's cost against `dab3876`
on the same machine/dependencies, with five alternating rounds and one numerical
thread. Timers include construction, fit and all four aggregations; inputs,
imports and warmup are excluded. Fast-mode medians are 114.37 → 117.76 ms for the
core panel (+3.0%), 190.29 → 197.07 ms with 12 covariates (+3.6%), and 280.27 →
297.87 ms with 24 periods (+6.3%). The extra final-score matrix calculation
adds modest fitting cost. The correction is retained for valid inference;
this follow-up does not establish a new speedup. All 30 executions preserve
point estimates, support, warnings and RNG state exactly. Inference is excluded
from that old-source equality check because its correction is intentional;
the permanent current-R gate validates it instead. Raw timing repeats, source
hashes, environment and checks are in the companion evidence file.
