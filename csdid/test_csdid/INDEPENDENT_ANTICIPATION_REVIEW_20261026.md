# Independent anticipation and regression audit, 2026-10-26

The candidate `81fa180` was checked against the clean prior PR head `2a36b48`
in separate Python processes under the same dependencies. No runtime source
was changed during this audit. The evidence supports the specific improvements
and unchanged workflows below; it cannot establish absence of every regression.

## Evidence independent of R

A noiseless two-cohort example sets anticipation effects to exactly **0.8**
in periods 3 and 4, treatment effects to **1.5** from period 5, and the last
unaffected reference to period 2. Under the old varying base, period 4 compared
two already anticipating periods and returned approximately zero. The corrected
code returns 0.8. DR, regression and IPW, both execution modes, and both base
choices give **12/12 correct candidate fits**, with maximum ATT error
`6.7e-16`. The old code fails all six varying-base fits by 0.8; its six universal
fits were already correct. Positive varying sampling weights are included.
This example checks effects; its deterministic residuals are not evidence for
inferential coverage.

An additional **96 noisy balanced-panel fits** are checked against an oracle
built directly from weighted means and centered residual contributions. The
oracle independently selects treatment, control and reference samples, derives
the full unit influence vector and analytical SE, and computes the Wald
quadratic form using only `t < g - anticipation`. It calls neither R nor the
package's cell-selection or estimator helpers. It crosses all methods, both
execution modes and bases, both control choices, weighted/unweighted fits and
anticipation 1/2. Every candidate comparison passes at `atol=1e-10, rtol=1e-9`;
the maximum absolute discrepancy across the checked quantities is `7.2e-14`.
The old source fails the independent oracle in all 96 fits because its Wald
restrictions include declared anticipation; 24 fits additionally change window
ATTs under the varying base with anticipation 2.

The oracle initially compared unrounded chi-square tail probabilities with
the existing public `Wpval`, which is rounded to five decimals. This caused
48 presentation-only mismatches. The oracle now applies that documented-in-code
public presentation convention. The inputs, estimator and numeric tolerance
were not changed; ATT, SE, full IF and Wald-statistic comparisons already passed.
This correction remains recorded in the machine-readable evidence.

## Before/after regression bounds

The **576 anticipation-zero workflows** cross three sampling designs, three
methods, both execution modes, both base choices, both control groups,
invariant/varying covariates, and four preprocessing patterns: clean data,
missing outcomes, missing rows, and absence of never-treated units. Healthy
varying weights and clustering are rotated through the grid, rather than
claimed as a complete cross of weight and cluster choices. Weight-anchor and
extreme-support grids are covered by a separate audit.

All public estimates, SEs, intervals, critical values, Wald outputs and every
cell/aggregate IF component agree **exactly** with the prior source, including
unit/cell/event support and nonfinite masks. These comparisons cover
13,266,912 numeric elements. NumPy's RNG state is unchanged by every analytic
fit and aggregation, and the caller's input frame is preserved exactly.

For the 96 nonzero-anticipation cases, true pre-period and post-treatment cells
retain their old ATT/SE/interval/full-IF outputs. The simple, group and calendar
aggregations also remain unchanged. The anticipation window, Wald test and
dynamic aggregation are intentionally excluded from before/after equality;
the first two are checked against the independent oracle instead.

Six additional seeded bootstrap workflows cover all three sampling designs and
both execution modes, with weights, clustering, varying covariates and missing
outcomes. Their complete outputs, masks, full IFs and final RNG states match
the old source exactly over another 156,756 numeric elements. Each uses 63
draws. These six cases do not exhaust bootstrap degeneracies or RNG backends.

The complete audit has **678/678 passing checks**, covering 14,201,076 numeric
comparisons, plus the separately reported 12-fit noiseless counterexample.

## Independent source review

- The shared cell planner changes the reference only when
  `t >= g - anticipation`; it selects the last *observed* period strictly before
  that threshold. The post-treatment flag remains `t >= g`, so anticipation
  cells do not enter post-treatment-only aggregations.
- The overlap cache now includes the selected weight period. The existing
  cohort-local cache and exact unit mask prevent reuse across changing controls;
  covariate and weight keys distinguish time-varying inputs. Healthy cached and
  uncached behavior is covered by both modes and covariate profiles above.
- Direct column indexing preserves the original sample predicates while fixing
  names that query strings interpret as expressions. The audit confirms normal
  names, missing-data filtering, no-never fallback and caller-data preservation;
  dedicated portable tests separately cover unusual names.
- The reference exporter correctly resolves synthetic RC row numbers through
  retained data. Synthetic row numbers cannot index the original CSV after
  preprocessing deletes observations. The earlier R grid explicitly covers
  that regression at anticipation 2.
- The earlier 24-case R `balance_e` grid exercises irregular clocks and display
  limits. Deliberately missing ATT cells are checked by the portable balance
  tests; they should not be described as a separate cross-language missing-ATT
  grid. Raw large-scale IF failures remain disclosed in the earlier evidence.

No additional runtime defect was found in this review of anticipation,
preprocessing, cache interaction or evidence claims. Weighted numerical guards,
bootstrap degeneracies and balance selection also have separate specialized
reviews. Agreement with older code establishes preservation on the tested
grid; the direct oracle supplies the independent correctness check for the
intercept-only balanced-panel reduction. Neither establishes identification
assumptions, repeated-sample coverage, or correctness at arbitrary magnitudes.

Reproduce with the same Python environment and clean checkouts:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python csdid/test_csdid/r_ref/current/anticipation_regression_compare.py \
  --baseline /clean/2a36b48 --candidate /clean/81fa180 \
  --output /tmp/csdid-anticipation-proof
```

[The compressed evidence](r_ref/current/anticipation_regression_comparison_20261026.json.gz)
retains all case definitions and outcomes, earlier-source oracle failures,
counterexample outputs, exact runtime-source hashes, clean-source checks,
dependency versions, numeric tolerances and the oracle correction.
