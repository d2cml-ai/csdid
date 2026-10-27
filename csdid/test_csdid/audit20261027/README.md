# Sampling-weight units: cross-package validation

The corrected Python implementation and corrected R `did` sources agree on
complete public output in **240 cross-language comparisons**, using 120 fits in
each language. Every comparison passes the existing strict comparator without
changing tolerances, support rules, nonfinite masks or frozen reference files.

The grid contains two-period single-cohort data and four-period data with a
never-treated group and three treated cohorts. Both panel and repeated-cross-
section routes use the same authoritative CSV rows in both languages. Each
uses DR, normalized IPW and regression, both execution modes, and weight scales
`1`, `1e-200`, `1e-160`, `1e155` and `1e307`.

The check compares all ATT cells, confidence limits, analytical SEs, Wald
statistics and full influence functions. All four aggregations are compared,
including their component estimates, SEs, critical values and complete overall
and component influence functions. Retained original observation IDs determine
row alignment, including R's synthetic RC IDs. It also checks 60 mode pairs in
each language and 96 weight-unit comparisons in each language: **552 complete
output checks**, with zero failures.

Maximum observed cross-language absolute errors are `3.77e-12` for ATT,
`7.58e-10` for cell/aggregate IFs and `5.83e-13` for SEs. Numerical tolerances
remain `atol=1e-7, rtol=1e-6`, with the existing `1e-5` Wald tolerances. Shapes,
unit/cell/event support and nonfinite masks must match exactly. This audit uses
analytical inference; it makes no cross-language RNG or bootstrap claim.

The reviewed R `did` correction is commit `75168454d6aa7faa2f9167b05ae7ad4f38f9e27b`;
R DRDID is official commit `eaac3b49cd225d0d7706c633858555a5f204e21e`. Python
uses the existing PyPI `drdid` 1.1.6 dependency, demonstrating that this CSDID
repair works before the standalone upstream correction is released. The compact
[validation record](weight_scale_validation.json) retains exact runtime-file,
input, exporter, comparator and payload hashes. Its Python Git head is the
worktree's parent revision; the recorded runtime-file hashes identify the
reviewed changes before their commit. Sources and frozen fixtures were unchanged
throughout the run.

The new `test_sampling_weight_scale.py` supplies separate regression evidence:
159 tests pass with the correction; the identical tests produce 120 failures
and 39 passes on the original runtime. Its ratio-of-cohort-masses derivative
oracle and single-cell aggregation nesting check do not use R output as an
expected answer. It also protects raw input columns and locally positive
weights when an unrelated cohort has vastly larger weights.

## Replay

Use clean R checkouts at the commits above, with their R dependencies,
`pkgload`, `jsonlite` and `data.table` available. Install this Python checkout and
its ordinary test dependencies. From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 python \
  csdid/test_csdid/audit20261027/compare_weight_scale.py \
  --did /path/to/reviewed/did --drdid /path/to/official/DRDID \
  --output /tmp/new-empty-weight-scale-audit
```

The output must be new or empty and outside the repository. It contains shared
CSV inputs, case specifications, complete compressed R/Python payloads,
`sessionInfo.txt`, the R generation log and `result.json`. These large payloads
are intentionally not committed. Any execution/comparison failure prevents a
passing result. The script hashes the unchanged committed fixture tree before
and after execution.

By default the script regenerates both seeded input designs. A separate replay
of those generators reproduced all four input hashes recorded in the completed
run. The original audit used existing shared CSV authorities via `--input-dir`;
that directory must contain `panel.csv`, `rc.csv`, `multi_panel.csv` and
`multi_rc.csv`. No input or case was selected in response to a numerical result.

The R exporter adapts the numeric schema and retained-ID mapping in
`r_ref/current/generate.R`, but never sources or executes that generator. The
Python runner reuses `_run_case` and `_assert_case_output` from the permanent
parity test, so this audit does not introduce a weaker numerical gate.
