# PyFixest performance audit

Audit date: 2026-10-23. Source inspected at
[`4245747fad805ec32c013ea4c0170a020e9ceded`](https://github.com/py-econometrics/pyfixest/tree/4245747fad805ec32c013ea4c0170a020e9ceded).
The separately installed benchmark wheel is **PyFixest 0.60.0**; it is not assumed
to be identical to this source revision.

## Mechanisms worth transferring

| Priority | Upstream mechanism | Application to csdid | Status |
|---|---|---|---|
| High | [Count/group cluster IDs once, then reduce score rows](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/src/crv1.rs#L5) | Replace bootstrap's full observation scan for every cluster with stable grouped indices. Preserve original row order and NumPy reduction dtype. | Implemented independently in `utils/mboot.py`; exact floating-point and RNG parity tests in `test_speed_cluster_groups.py`. |
| High | [Cache transformed columns and compute only missing columns](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/pyfixest/estimation/internals/demean_.py#L189) | Cache identical nuisance fits across related ATT cells. Keys must include sample, treatment/control membership, covariate period/design, and weights. Same-cohort never-treated post-period fits can sometimes share propensity estimates. | Implemented independently for supported balanced-panel faster-mode DR/IPW fits; see `test_performance_nuisance_cache.py`. |
| Medium | [Limit caches to compatible blocks and release unneeded model state](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/pyfixest/estimation/runner.py#L86) | Ensure point-estimate-only calls skip influence-function algebra inside the estimator, not just storage in the outer loop. Limit future caches to prevent memory growth. | Recommendation; no inference contract change made here. |
| Conditional | [Precompute group weights and reuse per-thread scratch arrays](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/pyfixest/estimation/numba/demean_nb.py#L65) | Consider compiled grouped reductions or fused influence arithmetic only after profiling shows NumPy/Python overhead remains important. | No compiled dependency added. |

Cache reuse is invalid when control sets, sample membership, covariates, weights,
or estimator settings change. In particular, cohort identity alone is not a
sufficient cache key. Read-only cache arrays, as used upstream, help prevent
accidental mutation across cells.

The implemented propensity cache holds only the most recent fit within each
cohort, remains local to a public fit call, and retains converged finite fits
without captured warnings. It stores the fitted propensity model, so each cell
still applies its own trimming and influence-function algebra. Custom estimators,
intercept-only designs, replaced fitting functions, and unsupported designs keep
the original route. Exact equality checks identify invariant covariates/weights;
changing control membership is included in the key.

The bootstrap grouping change is independently implemented with NumPy. **No
upstream implementation code was copied.** Its tests cover unsorted/noncontiguous
labels, strings, NaN/NaT, empty/singleton groups, C/F array layouts, float32/float64,
and unchanged bootstrap outputs and subsequent RNG draws. A measured
small-cluster fallback retains the original loop for at most 16 clusters.
Object labels also retain that loop because object NaNs can change the unique
label sequence when inverse codes are requested. Low and high cluster counts
are benchmarked separately; raw measurements retain the extra uniqueness-scan
cost of this dispatch for large cluster counts.

## Mechanisms that do not transfer directly

- PyFixest [event-study choices](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/pyfixest/did/estimation.py#L28)
  are TWFE, Gardner DID2S, and saturated cohort-specific regression. Its
  [saturated model](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/pyfixest/did/saturated_twfe.py#L100)
  follows Sun and Abraham. These are not substitutes for Callaway–Sant'Anna
  doubly robust group-time ATT in general staggered designs.
- High-dimensional fixed-effect demeaning, sparse LSMR/preconditioning, and GPU
  backends solve a different computational problem from csdid's current
  low-dimensional two-period nuisance regressions.
- Old compressed/Mundlak regression advice is obsolete for this revision:
  [`use_compression=True` raises `NotImplementedError`](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/pyfixest/estimation/api/feols.py#L547).
  Collapsing observations also need not retain the covariate, sampling-weight,
  cluster, or influence-function information required by CS inference.

## Restricted numerical comparison and timing

`benchmark_pyfixest_2x2.py` compares four public-API routes on a balanced panel
with two periods, one treated cohort, never-treated controls, no covariates, and
no sampling weights:

1. csdid standard DR;
2. csdid faster-mode DR;
3. PyFixest TWFE with unit-clustered covariance and both small-sample corrections
   disabled (`k_adj=False`, `G_adj=False`);
4. PyFixest first-difference regression with HC0 covariance, including reshaping
   inside the timed call.

All four must match a direct difference-in-means ATT and its analytical
influence-function standard error before timing is reported. Confidence-interval
critical values and p-values are not compared. This restricted equivalence does
not establish equivalence or a speed ranking for staggered treatment, covariate
adjustment, bootstrap bands, repeated cross sections, or arbitrary weights.

The script uses a common environment, fixes numerical/Rayon thread counts,
excludes data generation and imports, warms every method, rotates the timing
order, preserves default PyFixest copy/storage settings, and saves every repeat
plus versions and csdid source hashes. Upstream's
[modular benchmark](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/benchmarks/modular/feols_benchmarkers.py#L245)
turns off input copying and data storage, so its published results should not be
read as default-API comparisons against csdid.

From the csdid checkout, create a separate environment without modifying project
dependencies:

```sh
uv venv /tmp/csdid-pyfixest-compare-venv --python 3.12
uv pip install --python /tmp/csdid-pyfixest-compare-venv/bin/python \
  --only-binary :all: pyfixest==0.60.0 statsmodels drdid patsy threadpoolctl
PYTHONPATH=. /tmp/csdid-pyfixest-compare-venv/bin/python \
  csdid/test_csdid/benchmark_pyfixest_2x2.py --units 400 --parity-only
PYTHONPATH=. /tmp/csdid-pyfixest-compare-venv/bin/python \
  csdid/test_csdid/benchmark_pyfixest_2x2.py --units 2000,10000,50000 \
  --repeats 5 --output /tmp/csdid-pyfixest-2x2.json
```

Use an otherwise idle machine. The output records resolved dependency versions;
pin those versions as well when reproducing a particular recorded run. The
400-unit initial smoke check matched ATT/SE to at most `4.44e-16`; it was not a
timing run.

## Dependencies and license

Audited master uses [maturin/Rust, `within-py`, and optional Numba/Torch extras](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/pyproject.toml#L1).
Adding these only for grouped sums would create build/runtime costs without an
established benefit over NumPy. PyFixest is
[MIT licensed](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/LICENSE.md),
with [additional MIT notices for adapted Formulaic code](https://github.com/py-econometrics/pyfixest/blob/4245747fad805ec32c013ea4c0170a020e9ceded/THIRD_PARTY_NOTICES.md).
Future source reuse must retain applicable notices. This change adds no package
dependency and copies no upstream code.
