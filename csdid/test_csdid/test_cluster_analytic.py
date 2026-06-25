"""
Tests for analytical (no-bootstrap) cluster-robust standard errors.
Translated from R package 'did' tests/testthat/test-cluster-analytic.R (4 test_that blocks).

Tests that require faster_mode are skipped — not in Python csdid.
Tests that require internal cluster_vector or raw inffunc matrix exposure
are adapted to test observable behavior (clustered SE != unclustered SE).
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt


def _make_clustered_shocks(seed, G=50):
    """Panel where units within a cluster share a common shock each period."""
    rng = np.random.default_rng(seed)
    sz = np.tile([1, 2, 4, 10], G // 4 + 1)[:G]
    N = int(sz.sum())
    cl = np.repeat(np.arange(1, G + 1), sz)
    alpha = rng.standard_normal(G)[cl - 1]
    nu = rng.standard_normal(N)
    per = np.arange(1, 5)
    eta = rng.standard_normal((G, 4)) * 1.5
    g = rng.choice([2, 3, 0], size=G, p=[0.34, 0.33, 0.33])[cl - 1]

    rows = []
    for i in range(N):
        for t in per:
            te = 1.0 if (g[i] != 0 and t >= g[i]) else 0.0
            y = alpha[i] + nu[i] + 0.3 * t + eta[cl[i] - 1, t - 1] + te + rng.standard_normal()
            rows.append((i + 1, t, cl[i], g[i], y))

    df = pd.DataFrame(rows, columns=["id", "t", "cl", "g", "y"])
    return df.sort_values(["id", "t"]).reset_index(drop=True)


# ─────────────────────────────────────────────────────────────
# Analytical clustered SE differs from iid SE
# ─────────────────────────────────────────────────────────────

def test_clustered_se_differs_from_iid():
    """Analytical clustered SE differs from unclustered SE when units within a cluster are dependent."""
    d = _make_clustered_shocks(404)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res_cl = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated", clustervar="cl",
        ).fit(est_method="reg", bstrap=False, base_period="varying")

        res_iid = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated",
        ).fit(est_method="reg", bstrap=False, base_period="varying")

    se_cl = np.asarray(res_cl.results["se"], dtype=float)
    se_iid = np.asarray(res_iid.results["se"], dtype=float)
    finite = np.isfinite(se_cl) & np.isfinite(se_iid) & (se_iid > 0)
    # Clustered and iid SEs should not be identical
    assert not np.allclose(se_cl[finite], se_iid[finite], rtol=0.01), \
        "Clustered SE should differ from iid SE"


def test_clustered_se_runs_without_errors():
    """ATTgt with clustervar runs successfully and produces finite SEs."""
    d = _make_clustered_shocks(505)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated", clustervar="cl",
        ).fit(est_method="reg", bstrap=False, base_period="varying")

    se = np.asarray(res.results["se"], dtype=float)
    assert np.any(np.isfinite(se) & (se > 0))


# ─────────────────────────────────────────────────────────────
# Clustered SE flows through aggte
# ─────────────────────────────────────────────────────────────

def test_clustered_se_through_aggte():
    """Analytical clustered SE propagates through aggte for simple/group/dynamic."""
    d = _make_clustered_shocks(505)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res_cl = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated", clustervar="cl",
        ).fit(est_method="reg", bstrap=False, base_period="varying")

        res_iid = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated",
        ).fit(est_method="reg", bstrap=False, base_period="varying")

    for typec in ["simple", "group", "dynamic"]:
        import copy
        obj_cl = copy.copy(res_cl)
        obj_cl.MP = dict(res_cl.MP)
        obj_cl.did_object = dict(res_cl.did_object)
        obj_iid = copy.copy(res_iid)
        obj_iid.MP = dict(res_iid.MP)
        obj_iid.did_object = dict(res_iid.did_object)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obj_cl.aggte(typec=typec, bstrap=False)
            obj_iid.aggte(typec=typec, bstrap=False)

        se_cl = obj_cl.atte["overall_se"]
        se_iid = obj_iid.atte["overall_se"]
        assert np.isfinite(se_cl) and se_cl > 0, f"typec={typec}"


# ─────────────────────────────────────────────────────────────
# Skipped: analytical cluster SE vs bootstrap (skip_on_cran in R)
# ─────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────
# Analytical cluster SE vs cluster bootstrap (slow; skip_on_cran in R)
# ─────────────────────────────────────────────────────────────

@pytest.mark.slow
def test_analytical_cluster_se_vs_bootstrap():
    """Analytical clustered SE agrees with the cluster bootstrap SE (panel)."""
    d = _make_clustered_shocks(505)
    common = dict(yname="y", tname="t", idname="id", gname="g", data=d,
                  control_group="nevertreated", clustervar="cl")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        se_an = np.asarray(ATTgt(**common).fit(
            est_method="reg", bstrap=False, base_period="varying").results["se"], dtype=float)
        np.random.seed(0)
        se_bs = np.asarray(ATTgt(**common, biters=3000).fit(
            est_method="reg", bstrap=True, base_period="varying").results["se"], dtype=float)
    valid = np.isfinite(se_an) & np.isfinite(se_bs) & (se_bs > 0)
    assert valid.any()
    np.testing.assert_allclose(se_an[valid], se_bs[valid], rtol=0.25)


# ─────────────────────────────────────────────────────────────
# Skipped: repeated cross-sections clustered SE (requires faster_mode)
# ─────────────────────────────────────────────────────────────

def test_clustered_se_repeated_cross_sections():
    """Clustered analytical SE on repeated cross sections; faster_mode agrees with standard."""
    d = _make_clustered_shocks(808)
    common = dict(yname="y", tname="t", idname="id", gname="g", data=d,
                  control_group="nevertreated", clustervar="cl", panel=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        std = ATTgt(**common, faster_mode=False).fit("reg", bstrap=False, base_period="varying")
        fast = ATTgt(**common, faster_mode=True).fit("reg", bstrap=False, base_period="varying")
    se_s = np.asarray(std.results["se"], dtype=float)
    se_f = np.asarray(fast.results["se"], dtype=float)
    valid = np.isfinite(se_s) & np.isfinite(se_f)
    assert valid.any()
    np.testing.assert_allclose(se_s[valid], se_f[valid], atol=1e-9)


@pytest.mark.slow
def test_clustered_bootstrap_vs_analytical_rcs():
    """Analytical vs bootstrap clustered SE agree for repeated cross-sections."""
    d = _make_clustered_shocks(606)
    common = dict(yname="y", tname="t", idname="id", gname="g", data=d,
                  control_group="nevertreated", clustervar="cl", panel=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        se_an = np.asarray(ATTgt(**common).fit(
            est_method="reg", bstrap=False, base_period="varying").results["se"], dtype=float)
        np.random.seed(0)
        se_bs = np.asarray(ATTgt(**common, biters=3000).fit(
            est_method="reg", bstrap=True, base_period="varying").results["se"], dtype=float)
    valid = np.isfinite(se_an) & np.isfinite(se_bs) & (se_bs > 0)
    assert valid.any()
    np.testing.assert_allclose(se_an[valid], se_bs[valid], rtol=0.3)


# ─────────────────────────────────────────────────────────────
# Parametrized: clustered SE runs with all est_methods
# ─────────────────────────────────────────────────────────────

EST_METHODS = ["dr", "reg", "ipw"]


@pytest.mark.parametrize("em", EST_METHODS)
def test_clustered_se_runs_all_methods(em):
    """ATTgt with clustervar runs with all est_methods."""
    d = _make_clustered_shocks(606)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated", clustervar="cl",
        ).fit(est_method=em, bstrap=False, base_period="varying")
    se = np.asarray(res.results["se"], dtype=float)
    assert np.any(np.isfinite(se) & (se > 0))


@pytest.mark.parametrize("em", EST_METHODS)
@pytest.mark.parametrize("typec", ["simple", "group", "dynamic"])
def test_clustered_se_through_aggte_all_methods(em, typec):
    """Clustered SE propagates through aggte for all methods and types."""
    d = _make_clustered_shocks(707)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated", clustervar="cl",
        ).fit(est_method=em, bstrap=False, base_period="varying")

    import copy
    obj = copy.copy(res)
    obj.MP = dict(res.MP)
    obj.did_object = dict(res.did_object)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.aggte(typec=typec, bstrap=False)
    assert np.isfinite(obj.atte["overall_se"])
    assert obj.atte["overall_se"] > 0
