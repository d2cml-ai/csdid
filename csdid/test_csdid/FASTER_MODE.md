# csdid `faster_mode` — Implementation, Correctness & Benchmark

`faster_mode=True` is the Python equivalent of R `did`'s `faster_mode=TRUE`. It
produces results **numerically identical** to the standard path but eliminates
the two dominant per-(g,t) costs (profiled in the standard loop):

- per-cell `patsy` covariate construction → the design matrix is built **once
  per period** (panel) or **once globally** (repeated cross sections);
- per-cell `panel2cs2` reshaping → outcomes/weights are **pre-pivoted** to a
  per-period lookup (panel).

Implementation: `csdid/attgt_fnc/compute_att_gt2.py`, routed from
`ATTgt(..., faster_mode=True)`.

## Correctness — identical to the standard path

Validated across **72 configurations** (3 datasets × {panel, RC} ×
{nevertreated, notyettreated} × {dr, reg, ipw} × {varying, universal}):

| Quantity | Max |standard − fast| |
|----------|----------------------|
| ATT(g,t) | **1.4e-17** (bit-identical) |
| Standard error | **1.4e-17** |
| Influence functions | < 1e-9 (see `TestIFConsistency`) |

This equivalence is locked in by the suite: `TestFasterMode` (8),
`TestIFConsistency` (5), and `test_faster_mode_consistency.py` (24 parametrized)
all assert `faster_mode == standard`.

## Speed — 4-way comparison

`dr` estimator, `bstrap=False`, best of 5 runs, on this machine
(R `did` 2.5.0 / R 4.4.2; Python 3.13). Reproduce with
`r_ref/benchmark.py` and `r_ref/benchmark.R`.

| Dataset | rows | R (Standard) | R (Fast) | Python (Standard) | **Python (Fast)** |
|---------|-----:|-------------:|---------:|------------------:|------------------:|
| tp8_dyn (2k units × 8 periods) | 16,000 | 115 ms | 133 ms | 505 ms | **138 ms** |
| tp10_const (1k units × 10 periods) | 10,000 | 118 ms | 114 ms | 726 ms | **220 ms** |
| big_tp10 (8k units × 10 periods) | 80,000 | 396 ms | 453 ms | 1387 ms | **423 ms** |

### Takeaways
- **Python `faster_mode` ≈ 3.3× faster** than Python standard (3.65× / 3.31× /
  3.28×) — matching R's documented "2.5–3×" headline.
- **Python (Fast) ≈ R (Standard)** in absolute time (e.g. 138 ms vs 115 ms),
  largely closing the Python↔R gap that comes from R's compiled `DRDID`.
- **R `faster_mode` is roughly break-even** at these sizes (its tensor-setup
  overhead ≈ its savings); R's gains appear on much larger panels (more
  units × periods × cohorts). The Python port pays off earlier because the
  cost it removes — per-cell `patsy`/`pandas` reshaping — is proportionally
  larger in Python.

## Usage

```python
from csdid.att_gt import ATTgt
res = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=df,
            faster_mode=True).fit(est_method="dr")   # identical results, ~3x faster
```

`faster_mode` composes with everything else (panel/RC, control groups, all
est_methods, `base_period`, `fix_weights`, clustering, `bstrap`).
