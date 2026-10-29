"""
Tests for cluster-robust multiplier bootstrap.
Translated from R package 'did' tests/testthat/test-mboot-cluster.R (3 test_that blocks).

The bootstrap-heavy tests are marked slow (skip_on_cran in R).
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt


def _make_clustered(seed, bal, G=40):
    """Small staggered panel with units grouped into clusters."""
    rng = np.random.default_rng(seed)
    sz = np.full(G, 4) if bal else np.tile([1, 2, 3, 8], G // 4 + 1)[:G]
    N = int(sz.sum())
    cl = np.repeat(np.arange(1, G + 1), sz)
    a = rng.standard_normal(G)[cl - 1]
    nu = rng.standard_normal(N)
    g = np.where(rng.random(N) < 0.5, 2, 0)

    rows = []
    for i in range(N):
        for t in [1, 2, 3]:
            te = 1.0 if (g[i] == 2 and t >= 2) else 0.0
            y = a[i] + nu[i] + 0.5 * t + te + rng.standard_normal()
            rows.append((i + 1, t, cl[i], g[i], y))

    df = pd.DataFrame(rows, columns=["id", "t", "cl", "g", "y"])
    return df.sort_values(["id", "t"]).reset_index(drop=True)


# ─────────────────────────────────────────────────────────────
# Clustered bootstrap with unbalanced clusters
# ─────────────────────────────────────────────────────────────

def test_clustered_mboot_unbalanced_runs():
    """Clustered multiplier bootstrap runs for unbalanced clusters."""
    d = _make_clustered(101, bal=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated", clustervar="cl",
        ).fit(est_method="reg", bstrap=True, base_period="varying")

    se = np.asarray(res.results["se"], dtype=float)
    assert np.any(np.isfinite(se) & (se > 0)), "Bootstrap with clustering should produce finite SEs"


# ─────────────────────────────────────────────────────────────
# Clustered bootstrap with balanced clusters
# ─────────────────────────────────────────────────────────────

def test_clustered_mboot_balanced_runs():
    """Clustered multiplier bootstrap runs for balanced clusters."""
    d = _make_clustered(202, bal=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="y", tname="t", idname="id", gname="g", data=d,
            control_group="nevertreated", clustervar="cl",
        ).fit(est_method="reg", bstrap=True, base_period="varying")

    se = np.asarray(res.results["se"], dtype=float)
    assert np.any(np.isfinite(se) & (se > 0)), "Bootstrap with balanced clustering should produce finite SEs"


# ─────────────────────────────────────────────────────────────
# Clustering validation
# ─────────────────────────────────────────────────────────────

def test_clustering_validation_rejects_multiple_clustervars():
    """At most one cluster variable beyond idname is allowed."""
    d = _make_clustered(303, bal=False)
    d["cl2"] = d["cl"]

    # Python csdid takes clustervar as a single string, not a list.
    # Passing an invalid clustervar should error or be rejected.
    # The validation may happen in mboot when bstrap=True.
    with pytest.raises(Exception):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ATTgt(
                yname="y", tname="t", idname="id", gname="g", data=d,
                control_group="nevertreated", clustervar=["cl", "cl2"],
            ).fit(est_method="reg", bstrap=True, base_period="varying")
