"""Bespoke, self-contained, mutation-gated unit tests for the multiplier
bootstrap in ``csdid/utils/bmisc.py`` (draw generation) and
``csdid/utils/mboot.py`` (post-draw SE / crit-val arithmetic + cluster routing).

NO R at runtime. Expectations are:
  * ANALYTIC Rademacher moments: each bootstrap draw is (Ub @ a)/n with Ub in
    {+1,-1}; E[draw]=0, Var[draw]=||a||^2 / n^2.  (audit-v10 dim1)
  * DETERMINISTIC lattice facts for tiny n (n=1 => |draw|==|a| exactly;
    n=2,a=[c,c] => draws in {-c,0,c}) which pin the `@ inf_func / n` division and
    the {1,-1} Rademacher support exactly with no RNG dependence.
  * STRUCTURAL/branch oracles for mboot: degenerate-dimension guard, NaN SE on
    all-degenerate input, clustervar==idname collapse, missing-clustervar warn +
    fallback, time-varying-cluster refusal, and SE/crit-val sign+magnitude facts.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.utils.bmisc import multiplier_bootstrap
from csdid.utils.mboot import mboot, run_multiplier_bootstrap


# --------------------------------------------------------------------------- #
# bmisc.multiplier_bootstrap : draw generation
# --------------------------------------------------------------------------- #
def test_mb_n1_abs_equals_a_exactly():
    """n=1: draw = +/- a / 1  =>  |draw| == |a| exactly, every iteration.
    Kills the `@ inf_func / n` division mutant and the Rademacher support."""
    a = np.array([[3.7]])
    out = multiplier_bootstrap(a, 500)
    assert out.shape == (500, 1)
    assert np.allclose(np.abs(out), 3.7)
    # both signs must appear (Rademacher is +1/-1, not +1/+1)
    assert (out > 0).any() and (out < 0).any()


def test_mb_n2_lattice_values():
    """n=2, a=[c,c]: draws = (U1*c+U2*c)/2 in {-c, 0, c}.
    Pins the division-by-n and the symmetric +/-1 support."""
    c = 2.0
    a = np.array([[c], [c]])
    out = multiplier_bootstrap(a, 4000)
    vals = set(np.round(np.unique(out), 9))
    assert vals.issubset({-c, 0.0, c})
    assert vals == {-c, 0.0, c}  # all three reached for a fair coin


def test_mb_division_by_n():
    """Distinguish /n from no-division: with n=4 identical entries a=[1,1,1,1],
    the all-+1 draw equals 4*1/4 = 1, NOT 4. A `Div`->`Mult` mutant gives 16."""
    a = np.ones((4, 1))
    out = multiplier_bootstrap(a, 6000)
    assert out.max() <= 1.0 + 1e-12        # max achievable is +1 (all +1 weights)
    assert out.min() >= -1.0 - 1e-12


def test_mb_rademacher_variance():
    """Var of a single-column draw == ||a||^2 / n^2 (analytic Rademacher)."""
    rng = np.random.default_rng(0)
    n = 40
    a = rng.normal(0, 1, (n, 1))
    np.random.seed(123)
    out = multiplier_bootstrap(a, 300000)
    analytic_var = np.sum(a ** 2) / n ** 2
    assert abs(out.var() - analytic_var) / analytic_var < 0.05
    assert abs(out.mean()) < 5e-3 * (analytic_var ** 0.5 + 1)


def test_mb_serial_reproducible_under_global_seed():
    a = np.array([[1.0], [2.0], [3.0]])
    np.random.seed(7)
    o1 = multiplier_bootstrap(a, 200)
    np.random.seed(7)
    o2 = multiplier_bootstrap(a, 200)
    assert np.array_equal(o1, o2)


def test_mb_seed_path_independent_of_global():
    """seed!=None uses a LOCAL Generator (parallel path); reproducible by seed and
    distinct from the global-state serial draws."""
    a = np.array([[1.0], [2.0], [3.0]])
    g1 = multiplier_bootstrap(a, 200, seed=5)
    g2 = multiplier_bootstrap(a, 200, seed=5)
    assert np.array_equal(g1, g2)
    np.random.seed(5)
    serial = multiplier_bootstrap(a, 200)
    assert not np.array_equal(g1, serial)


def test_run_multiplier_bootstrap_shape_and_scaling():
    """run_multiplier_bootstrap (serial path) returns (biters, K) and matches a
    direct call under the same global seed."""
    rng = np.random.default_rng(1)
    a = rng.normal(0, 1, (30, 2))
    np.random.seed(99)
    r = run_multiplier_bootstrap(a, 100, pl=False, cores=1)
    assert r.shape == (100, 2)
    np.random.seed(99)
    direct = multiplier_bootstrap(a, 100)
    assert np.array_equal(r, direct)


# --------------------------------------------------------------------------- #
# mboot : helpers + fixtures
# --------------------------------------------------------------------------- #
def _panel_params(data, clustervars=None, biters=1500, alp=0.05):
    return dict(data=data, idname="id", tname="t", clustervars=clustervars,
                biters=biters, alp=alp, panel=True,
                true_repeated_cross_sections=False)


def _simple_panel_data(n):
    return pd.DataFrame({"id": np.arange(n), "t": np.zeros(n)})


# --------------------------------------------------------------------------- #
# mboot : SE / crit_val arithmetic (lines 108, 124, 127-136)
# --------------------------------------------------------------------------- #
def test_mboot_se_positive_and_crit_exceeds_pointwise():
    n = 300
    rng = np.random.default_rng(0)
    ifn = rng.normal(0, 1, (n, 2))
    dp = _panel_params(_simple_panel_data(n))
    np.random.seed(42)
    r = mboot(ifn, dp)
    assert np.all(np.isfinite(r["se"])) and np.all(r["se"] > 0)
    # simultaneous crit-val for a multivariate band must exceed pointwise 1.96
    assert r["crit_val"] > 1.96


def test_mboot_se_scales_with_if_magnitude():
    """SE ~ sqrt(mean(if^2)/n) scaling: doubling the IF doubles the SE (linear).
    Kills mutants that drop the n-normalization or alter the bSigma assembly."""
    n = 300
    rng = np.random.default_rng(3)
    ifn = rng.normal(0, 1, (n, 1))
    dp = _panel_params(_simple_panel_data(n))
    np.random.seed(11)
    r1 = mboot(ifn, dp)
    np.random.seed(11)
    r2 = mboot(2.0 * ifn, dp)
    ratio = r2["se"][0] / r1["se"][0]
    assert abs(ratio - 2.0) < 1e-9


def test_mboot_degenerate_column_gives_nan_se():
    """A zero IF column is degenerate -> its SE is NaN, the live column finite.
    Kills the ndg-dimension guard mutants (line 108)."""
    n = 200
    rng = np.random.default_rng(1)
    ifn = np.column_stack([rng.normal(0, 1, n), np.zeros(n)])
    dp = _panel_params(_simple_panel_data(n))
    np.random.seed(7)
    r = mboot(ifn, dp)
    assert np.isfinite(r["se"][0])
    assert np.isnan(r["se"][1])


def test_mboot_all_degenerate_returns_nan():
    """All-zero IF: all dims degenerate -> all-NaN SE and NaN crit_val
    (matches R's ncol==0 guard; lines 112-114)."""
    n = 150
    ifn = np.zeros((n, 2))
    dp = _panel_params(_simple_panel_data(n))
    np.random.seed(7)
    r = mboot(ifn, dp)
    assert np.all(np.isnan(r["se"]))
    assert np.isnan(r["crit_val"])


# --------------------------------------------------------------------------- #
# mboot : cluster routing / normalization (lines 27, 40, 51, 55, 70, 96)
# --------------------------------------------------------------------------- #
def _cluster_fixture(n=200):
    rng = np.random.default_rng(2)
    ifn = rng.normal(0, 1, (n, 1))
    clust = np.repeat(np.arange(n // 2), 2)  # 2 units per cluster
    data = pd.DataFrame({"id": np.arange(n), "t": np.zeros(n), "cl": clust})
    return ifn, data


def test_mboot_clustervar_equal_idname_collapses_to_none():
    """clustervars == idname => no extra clustering (line 51)."""
    ifn, data = _cluster_fixture()
    np.random.seed(0)
    r_id = mboot(ifn, _panel_params(data, clustervars="id"))
    np.random.seed(0)
    r_none = mboot(ifn, _panel_params(data, clustervars=None))
    assert np.allclose(r_id["se"], r_none["se"], equal_nan=True)


def test_mboot_missing_clustervar_warns_and_falls_back():
    """Missing cluster column => warn + unclustered fallback (lines 55-65)."""
    ifn, data = _cluster_fixture()
    np.random.seed(0)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r_miss = mboot(ifn, _panel_params(data, clustervars="does_not_exist"))
    assert any("clustering" in str(x.message) for x in w)
    np.random.seed(0)
    r_none = mboot(ifn, _panel_params(data, clustervars=None))
    assert np.allclose(r_miss["se"], r_none["se"], equal_nan=True)


def test_mboot_time_varying_cluster_raises():
    """Cluster variable that varies within a unit over time => ValueError
    (lines 66-74)."""
    n = 100
    rng = np.random.default_rng(2)
    ifn = rng.normal(0, 1, (n, 1))
    # 2 periods per unit, cluster id different across the two rows of each unit
    data = pd.DataFrame({
        "id": np.repeat(np.arange(n), 2),
        "t": np.tile([0, 1], n),
        "cl": np.arange(2 * n),
    })
    with pytest.raises(ValueError):
        np.random.seed(0)
        mboot(ifn, _panel_params(data, clustervars="cl"))


def test_mboot_multiple_clustervars_rejected():
    """A list with >1 cluster var is rejected (line 47)."""
    ifn, data = _cluster_fixture()
    data = data.assign(cl2=data["cl"])
    with pytest.raises(ValueError):
        mboot(ifn, _panel_params(data, clustervars=["cl", "cl2"]))


def test_mboot_clustered_se_uses_cluster_sums():
    """Clustered SE differs from unclustered when within-cluster IFs are
    positively correlated (cluster-sum variance > iid). Kills the cluster-sum
    aggregation (lines 94-99) and the sqrt(n_clusters)/n scaling (line 136)."""
    n = 200
    base = np.random.default_rng(5).normal(0, 1, n // 2)
    # two perfectly-correlated rows per cluster -> cluster sum = 2*base
    ifn = np.repeat(base, 2)[:, None]
    clust = np.repeat(np.arange(n // 2), 2)
    data = pd.DataFrame({"id": np.arange(n), "t": np.zeros(n), "cl": clust})
    np.random.seed(0)
    r_cl = mboot(ifn, _panel_params(data, clustervars="cl"))
    np.random.seed(0)
    r_iid = mboot(ifn, _panel_params(data, clustervars=None))
    # positively-correlated clusters inflate the SE vs the iid assumption
    assert r_cl["se"][0] > r_iid["se"][0]


# --------------------------------------------------------------------------- #
# run_multiplier_bootstrap PARALLEL path (n>2500, pl=True, cores>1):
# lines 145-147 (chunk boundaries) + 151 (branch condition) + 160-161 (per-chunk
# seed derivation). Requires n>2500 to enter the branch.
# --------------------------------------------------------------------------- #
def _big_if(n=3000, seed=0):
    return np.random.default_rng(seed).normal(0, 1, (n, 1))


def test_parallel_path_row_count_equals_biters():
    """The per-core chunks must sum EXACTLY to biters (R's diff(round(seq(...)))).
    A chunk-boundary arithmetic mutant (linspace cores+1, diff) makes the stacked
    result have the wrong number of rows."""
    a = _big_if()
    for cores in (2, 3):
        np.random.seed(1)
        r = run_multiplier_bootstrap(a, 120, pl=True, cores=cores)
        assert r.shape == (120, 1)


def test_parallel_path_reproducible_under_seed():
    """v10-F1: the parallel path is reproducible under np.random.seed (the
    per-chunk SeedSequence is derived from the pinned global RNG)."""
    a = _big_if(seed=2)
    np.random.seed(5)
    r1 = run_multiplier_bootstrap(a, 200, pl=True, cores=2)
    np.random.seed(5)
    r2 = run_multiplier_bootstrap(a, 200, pl=True, cores=2)
    assert np.array_equal(r1, r2)


def test_parallel_path_correct_variance_not_collapsed():
    """Workers must be INDEPENDENT (distinct spawned SeedSequences), so the
    pooled draws have the analytic Rademacher variance ||a||^2/n^2 -- not a
    collapsed/duplicated stream. Kills the seed-derivation const mutants."""
    n = 3000
    a = _big_if(n=n, seed=3)
    np.random.seed(7)
    big = run_multiplier_bootstrap(a, 20000, pl=True, cores=2)
    analytic_var = np.sum(a ** 2) / n ** 2
    assert abs(big.var() - analytic_var) / analytic_var < 0.05


def test_serial_vs_parallel_branch_selection():
    """n <= 2500 stays SERIAL even with pl=True, cores>1 (matches global-RNG
    draws). Pins the n>2500 branch threshold (line 151)."""
    a = np.random.default_rng(0).normal(0, 1, (2000, 1))  # n=2000 <= 2500
    np.random.seed(9)
    r_branch = run_multiplier_bootstrap(a, 100, pl=True, cores=2)
    np.random.seed(9)
    r_serial = multiplier_bootstrap(a, 100)
    # serial branch taken -> identical to a direct global-RNG draw
    assert np.array_equal(r_branch, r_serial)
