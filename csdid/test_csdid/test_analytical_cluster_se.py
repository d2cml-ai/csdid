"""Tests for analytical clustered SEs (Phase 3.1 of R sync).

R did v2.5.1 added analytical cluster-robust SEs without requiring
the bootstrap (bstrap=FALSE with clustervars set).
"""
import pytest
import numpy as np
import pandas as pd
from csdid.att_gt import ATTgt


@pytest.fixture
def clustered_data():
    """Panel data with cluster structure."""
    np.random.seed(42)
    n_clusters = 10
    units_per_cluster = 5
    n_units = n_clusters * units_per_cluster
    n_periods = 4
    
    cluster_id = np.repeat(np.arange(1, n_clusters + 1), units_per_cluster)
    groups_cluster = np.concatenate([
        np.full(3, 2005), np.full(3, 2006), np.full(4, 0),
    ])
    groups_unit = np.repeat(groups_cluster, units_per_cluster)
    
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(2003, 2007), n_units)
    groups = np.repeat(groups_unit, n_periods)
    clusters = np.repeat(cluster_id, n_periods)
    
    cluster_effect = np.random.randn(n_clusters)
    y = np.random.randn(n_units * n_periods)
    y += np.repeat(np.repeat(cluster_effect, units_per_cluster), n_periods)
    treated = (times >= groups) & (groups > 0)
    y[treated] += 0.5
    
    return pd.DataFrame({
        'id': ids, 'year': times, 'group': groups,
        'y': y, 'cluster': clusters,
    })


class TestAnalyticalClusteredSE:
    def test_analytical_clustered_se_runs(self, clustered_data):
        """Analytical clustered SEs should work without bootstrap."""
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data, clustervar="cluster", biters=100,
        )
        result = obj.fit(est_method='dr', bstrap=False)
        assert result is not None
        se_values = np.array(result.results['se'])
        # SEs should be non-zero
        assert np.all(se_values[~np.isnan(se_values)] > 0)

    def test_analytical_clustered_se_matches_bootstrap(self, clustered_data):
        """Analytical clustered SEs should match the cluster bootstrap SEs.

        Both are R's cluster-robust SE: sqrt(sum_c S_c^2)/n. (Comparing clustered
        vs i.i.d. is not meaningful here because DiD differences away the
        time-invariant cluster effect, so clustered can be smaller than i.i.d.)
        """
        obj_clust = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data, clustervar="cluster", biters=3000,
        )
        np.random.seed(0)
        result_an = obj_clust.fit(est_method='dr', bstrap=False)

        obj_bs = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data, clustervar="cluster", biters=3000,
        )
        np.random.seed(0)
        result_bs = obj_bs.fit(est_method='dr', bstrap=True)

        se_an = np.array(result_an.results['se'], dtype=float)
        se_bs = np.array(result_bs.results['se'], dtype=float)
        valid = ~np.isnan(se_an) & ~np.isnan(se_bs) & (se_bs > 0)
        # Analytical and bootstrap cluster-robust SEs should agree within ~20%.
        np.testing.assert_allclose(se_an[valid], se_bs[valid], rtol=0.2)

    def test_aggte_with_analytical_clustered_se(self, clustered_data):
        """aggte should work with analytical clustered SEs."""
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data, clustervar="cluster", biters=100,
        )
        result = obj.fit(est_method='dr', bstrap=False)
        agg = result.aggte(typec='simple')
        assert agg.atte is not None
        assert agg.atte['overall_se'] is not None
        assert agg.atte['overall_se'] > 0
