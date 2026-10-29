# Current-source R performance comparison

These scripts compare the same weighted synthetic panels as
`../benchmark_v043.py`: 8,000 units × 10 periods with analytical inference,
the same panel with 400 clusters, and 2,000 units × 8 periods with 399 bootstrap
draws and simultaneous bands. They measure the complete public workflow:
construction/preprocessing, fitting, and simple/group/dynamic/calendar
aggregation. R's `att_gt` includes preprocessing, so compare its `att_gt_s` with
Python's `construct_s + fit_s`, and compare `total_s` directly.

Run all CPU benchmarks serially, with the same thread environment. Inputs,
imports, warmups, garbage collection, snapshot extraction and serialization are
outside R timers. Each mode gets one warmup; mode order alternates across five
measured repetitions; medians and every raw repetition are retained. R profiles
run separately after timing. The profiler's memory columns describe allocations,
not peak process RSS. Python's existing benchmark includes small result-extraction
costs inside its timing window; R snapshot extraction is outside the timer.

```sh
python csdid/test_csdid/performance_audit_r/prepare_inputs.py \
  --output /tmp/csdid-r-perf-input

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
Rscript csdid/test_csdid/performance_audit_r/benchmark.R \
  --did /path/to/did --drdid /path/to/DRDID \
  --input /tmp/csdid-r-perf-input --output /tmp/csdid-r-perf \
  --reps 5 --profiles

python csdid/test_csdid/performance_audit_r/compare.py \
  --python-snapshot /path/to/python-benchmark/00-candidate/snapshot \
  --r-output /tmp/csdid-r-perf --output /tmp/csdid-r-perf/comparison.json
```

`--cases analytic,clustered` can limit R execution; the Python snapshot must have
the same case support. Source package paths are explicit and loaded with
`pkgload::load_all`, avoiding stale installed packages. Output records R/package
versions, thread settings, source commits/status, per-file source MD5 hashes,
input MD5/SHA256 hashes, raw timings, medians, full snapshots and optional Rprof
call/allocation profiles. Repeated R results must match exactly within each mode.
Influence rows are aligned by authoritative unit IDs before exporting.

The comparison requires analytical public outputs at `atol=1e-7, rtol=1e-6`
(Wald fields `1e-5`), and reports full influence-function discrepancies at the
same threshold separately. Add `--require-influence` to make those discrepancies
a failing gate too. R and Python retain different GLM solvers; residual IF
differences must be reported, not interpreted as an optimization regression.
Bootstrap comparisons deliberately check only point estimates and full influence
functions. Equal seeds in R and Python do **not** produce equal multiplier draws,
so this benchmark cannot assert bootstrap SE/band equality. The earlier shared
multiplier audit is the appropriate test for that requirement.

## Source inspection, 2026-10-23

The connected GitHub API confirmed local checkouts matched current upstream:

- `did` 2.5.1.902, commit
  [`74b88eb07f1faa644df1271055a0f77848f6550c`](https://github.com/bcallaway11/did/commit/74b88eb07f1faa644df1271055a0f77848f6550c).
- `DRDID` 1.4.0, commit
  [`eaac3b49cd225d0d7706c633858555a5f204e21e`](https://github.com/pedrohcgs/DRDID/commit/eaac3b49cd225d0d7706c633858555a5f204e21e).

R mechanisms that may transfer to Python, subject to profiling and parity checks:

- [Guard memoization](https://github.com/bcallaway11/did/blob/74b88eb07f1faa644df1271055a0f77848f6550c/R/compute.att_gt2.R#L689):
  cache overlap/regression-feasibility results only when all guard inputs repeat;
  preserve per-cell warning behavior. R restricts this to balanced panels,
  never-treated controls and compatible weight rules.
- [RC masks and unit aggregation](https://github.com/bcallaway11/did/blob/74b88eb07f1faa644df1271055a0f77848f6550c/R/compute.att_gt2.R#L651):
  precompute period/cohort masks, fixed-weight joins, and a sparse row-to-unit
  aggregation operator. These avoid repeated full-data scans and hashing.
- [Sparse influence assembly](https://github.com/bcallaway11/did/blob/74b88eb07f1faa644df1271055a0f77848f6550c/R/compute.att_gt2.R#L789):
  build sparse triplets directly; optional point-estimate-only calls omit
  influence computation entirely. Sparse storage is useful only if downstream
  covariance/aggregation/bootstrap code avoids immediately densifying it.
- [Low-level GLM entry point](https://github.com/pedrohcgs/DRDID/blob/eaac3b49cd225d0d7706c633858555a5f204e21e/R/fastglm_fit.R):
  bypass repeated high-level model wrappers while preserving solver defaults.
  The source comment's 25–30% per-fit claim is historical R evidence, not a
  measured Python speedup.
- [Bootstrap work suppression](https://github.com/bcallaway11/did/blob/74b88eb07f1faa644df1271055a0f77848f6550c/R/mboot.R):
  `return_V=FALSE` avoids a covariance matrix when callers only need standard
  errors and critical values; group influence functions before multiplier draws.

The older R profiling vignette compares `did` 2.0.0 to 2.1.0 and shows large
data-management savings, but its reported runtimes are not evidence about
current R/Python performance. The separate `did_speed_checks` project benchmarks
different estimators (CS, SA, BJS, Gardner, dCDH), so that ranking cannot establish
that one implementation is fastest for the same estimand and inference contract.
