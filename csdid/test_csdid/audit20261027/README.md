# Joint installed-package validation

The final tested pair is **CSDID 0.4.4** at
`e8354bc3454e5ed5bc78eb5f2df3c8b8422b24aa` and **Python DRDID 1.1.7** at
`eb896640c7fb74e2baac06ec0e208ee19f290d7c`. CSDID declares `drdid==1.1.7`;
its committed peer lock selects that exact reviewed source. These are tested
release candidates. Validation did not publish either distribution.

The [joint installed receipt](joint_installed_validation.json) records
**2,935 tests plus 19 subtests**, with zero failures, errors or skips. Slow tests
and publisher-helper tests were included. The 72-case direct-kernel bridge also
passes its unchanged raw-output and original-unit checks. Both wheels were built
from exact Git archives and installed together in a fresh virtual environment.
Source/wheel/installed-file fingerprints, module origins and unchanged plain
test snapshots were verified. JUnit reports 2,954 cases because it includes the
19 subtests; those are not additional test functions.

The final installed pair also passes **648 complete R comparison checks**:
552 checks on two-period and multiple-cohort panel/RC inputs, plus 96 checks on
sparse float32 weights. Every grid uses DR, normalized IPW and regression in
both execution modes. The main scales are `1`, `1e-200`, `1e-160`, `1e155` and
`1e307`. The sparse grid compares actual float32 unit masks and masks scaled by
`1.401298464324817e-45`; their float32 mean rounds to zero, whereas the exact
represented values have double mean `2.8025969286496343e-46`.

All ATT cells, confidence limits, analytical SEs, Wald statistics and full IFs
are compared, along with all four aggregations' component estimates, SEs,
critical values and full IFs. Original retained IDs determine alignment,
including R's synthetic RC IDs. Tolerances remain `atol=1e-7, rtol=1e-6`, with
the existing `1e-5` Wald tolerances; shapes, support and nonfinite masks match
exactly. Maximum IF error is `7.58e-10` on the main grid and `6.07e-11` on the
sparse grid. No fixture or source pin was advanced to obtain a pass.

The numerical receipts are [main grid](weight_scale_validation.json) and
[sparse float32](sparse_float32_validation.json). Their respective
[installed guard](weight_scale_installed_guard.json) and
[installed guard](sparse_float32_installed_guard.json) bind them to the joint
wheels: both runtime packages were preloaded under Python `-I`, and every loaded
runtime origin/file hash was verified before and after the source-based audit
helpers ran. The R sources are corrected did commit
`75168454d6aa7faa2f9167b05ae7ad4f38f9e27b` and official DRDID commit
`eaac3b49cd225d0d7706c633858555a5f204e21e`.

The [final installed healthy-identity proof](joint_installed_healthy_identity.json)
retains 240 workflows, 960 aggregations and 1,920 field hashes that are
bit-identical to the baseline. The older [healthy receipt](healthy_identity_summary.json),
[main-grid receipt](historical_weight_scale_validation_afa674b.json),
[sparse receipt](historical_sparse_float32_validation_afa674b.json) and
[initial summary](historical_validation_summary_afa674b.json) remain historical.
The [float32 before/after record](sparse_float32_before_after.json) preserves
its original failing cases and intermediate correction; no failure was removed
or tolerance widened.

## Preserved numerical boundary

The final installed pair replays the original 144-case outcome-unit grid:
**252 of 288 raw comparisons pass; 36 large-outcome comparisons still fail**
the declared full aggregate-IF tolerance at outcome scale `1e10`. All non-IF
public fields pass. All 72 original-unit diagnostic comparisons pass, but those
diagnostics do not override the raw failures or their exit status **1**. The
preserved R fast/slow grid also contains 42 failed IF fields.

The [installed raw-outcome receipt](raw_outcome_installed_replay.json) records
these failures, exact affected cases and artifact hashes. It reuses the original
R references (did `48863b6`, DRDID `eaac3b4`), rather than claiming fresh R
export. This evidence supports the reviewed pair and tested grids; it does not
establish universal full-IF equivalence in arbitrary outcome units. The scope
is the six shared estimators under equivalent settings, as requested.

## Replay

First run the [joint wheel gate](../joint_sync/README.md) using the reviewed
commits. It retains candidate wheels, full JUnit/logs, hashes and `receipt.json`
in a new external directory. Then use that gate's Python for the R audits:

```sh
/tmp/joint/venv/bin/python -I \
  csdid/test_csdid/audit20261027/installed_run.py \
  --joint /tmp/joint --audit compare_weight_scale \
  --guard-receipt /tmp/new-weight-audit/installed_guard.json -- \
  --did /path/to/reviewed/did --drdid /path/to/official/DRDID \
  --output /tmp/new-weight-audit
```

Use `--audit compare_sparse_float32` with another new output/guard directory for
the sparse proof. R needs the checkout dependencies, `pkgload`, `jsonlite` and
`data.table`. Numerical work uses one thread. The main audit regenerates seeded
inputs by default; a separate replay reproduced all four original CSV hashes.
`--input-dir` can instead supply the original `panel.csv`, `rc.csv`,
`multi_panel.csv` and `multi_rc.csv` authorities.

Full payloads stay outside the repository. The R exporter adapts the numeric
schema and retained-ID mapping from `r_ref/current/generate.R` without executing
that generator. The Python audit reuses the permanent `_run_case` and
`_assert_case_output` gate. Comparisons use analytical inference and make no
cross-language RNG/bootstrap claim. The compact [summary](validation_summary.json)
links the current receipts and retained history.
