"""Tests for cluster bootstrap fix (Phase 3.2 of R sync).

Verifies that the clustered multiplier bootstrap correctly aggregates
influence functions by cluster before drawing bootstrap weights.
"""
import pytest
import numpy as np
import pandas as pd
from csdid.att_gt import ATTgt


@pytest.fixture
def clustered_data():
    """Panel data with a cluster variable (e.g., state-level clustering)."""
    np.random.seed(42)
    n_clusters = 10
    units_per_cluster = 5
    n_units = n_clusters * units_per_cluster
    n_periods = 4
    
    # Assign units to clusters
    cluster_id = np.repeat(np.arange(1, n_clusters + 1), units_per_cluster)
    
    # Groups: first 3 clusters treated at 2005, next 3 at 2006, rest never-treated
    groups_cluster = np.concatenate([
        np.full(3, 2005),
        np.full(3, 2006),
        np.full(4, 0),
    ])
    groups_unit = np.repeat(groups_cluster, units_per_cluster)
    
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(2003, 2003 + n_periods), n_units)
    groups = np.repeat(groups_unit, n_periods)
    clusters = np.repeat(cluster_id, n_periods)
    
    # Add cluster-level random effect
    cluster_effect = np.random.randn(n_clusters)
    y = np.random.randn(n_units * n_periods)
    y += np.repeat(np.repeat(cluster_effect, units_per_cluster), n_periods)
    # Treatment effect
    treated = (times >= groups) & (groups > 0)
    y[treated] += 0.5
    
    return pd.DataFrame({
        'id': ids,
        'year': times,
        'group': groups,
        'y': y,
        'cluster': clusters,
    })


class TestClusteredBootstrap:
    def test_clustered_bootstrap_runs(self, clustered_data):
        """Basic clustered bootstrap should run without errors."""
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data,
            clustervar="cluster",
            biters=200,
        )
        result = obj.fit(est_method='dr')
        assert result is not None
        assert len(result.results['att']) > 0
        se_values = np.array(result.results['se'])
        # SEs should be non-zero for valid estimates
        assert np.all(se_values[~np.isnan(se_values)] > 0)

    def test_clustered_se_differs_from_unclustered(self, clustered_data):
        """Clustered SEs should generally differ from non-clustered SEs."""
        # Without clustering
        obj_noclust = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data,
            biters=500,
        )
        result_noclust = obj_noclust.fit(est_method='dr')
        
        # With clustering
        obj_clust = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data,
            clustervar="cluster",
            biters=500,
        )
        result_clust = obj_clust.fit(est_method='dr')
        
        se_noclust = np.array(result_noclust.results['se'])
        se_clust = np.array(result_clust.results['se'])
        
        # They should not be identical (within-cluster correlation inflates SEs)
        assert not np.allclose(se_noclust, se_clust, rtol=0.1), \
            "Clustered and non-clustered SEs should differ when there's within-cluster correlation"

    def test_clustered_se_differs_in_magnitude(self, clustered_data):
        """Clustered SEs should differ from non-clustered SEs (properly scaled)."""
        obj_noclust = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data,
            biters=500,
        )
        result_noclust = obj_noclust.fit(est_method='dr')
        
        obj_clust = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data,
            clustervar="cluster",
            biters=500,
        )
        result_clust = obj_clust.fit(est_method='dr')
        
        se_noclust = np.array(result_noclust.results['se'])
        se_clust = np.array(result_clust.results['se'])
        
        # Clustered SEs should be finite and positive
        valid = ~np.isnan(se_clust)
        assert np.all(se_clust[valid] > 0), "Clustered SEs should be positive"
        # SEs should be different (clustering changes variance estimation)
        assert not np.allclose(se_clust[valid], se_noclust[~np.isnan(se_noclust)]), \
            "Clustered SEs should differ from non-clustered"

    def test_unit_level_cluster_similar_to_no_cluster(self, clustered_data):
        """Clustering at unit level (each unit = own cluster) should give similar SEs."""
        # Create a unique-per-unit cluster variable
        data = clustered_data.copy()
        data['unit_cluster'] = data['id']
        
        obj_noclust = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=clustered_data,
            biters=500,
        )
        result_noclust = obj_noclust.fit(est_method='dr')
        
        # Cluster on a per-unit variable (= no intra-cluster correlation)
        obj_unit_clust = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=data,
            clustervar="unit_cluster",
            biters=500,
        )
        result_unit_clust = obj_unit_clust.fit(est_method='dr')
        
        se_noclust = np.array(result_noclust.results['se'])
        se_unit = np.array(result_unit_clust.results['se'])
        
        # Should be approximately equal (both are unit-level)
        valid = ~np.isnan(se_noclust) & ~np.isnan(se_unit)
        assert np.allclose(se_noclust[valid], se_unit[valid], rtol=0.3), \
            "Unit-level clustering should give similar SEs to no clustering"

    def test_time_varying_cluster_rejected(self):
        """Time-varying cluster variables should be rejected."""
        np.random.seed(42)
        n_units, n_periods = 20, 4
        ids = np.repeat(np.arange(1, n_units + 1), n_periods)
        times = np.tile(np.arange(2003, 2007), n_units)
        groups_unit = np.concatenate([np.full(10, 2005), np.full(10, 0)])
        groups = np.repeat(groups_unit, n_periods)
        y = np.random.randn(n_units * n_periods)
        # Time-varying cluster: changes every period for each unit
        cluster = (ids * 10 + times) % 5
        
        data = pd.DataFrame({
            'id': ids, 'year': times, 'group': groups,
            'y': y, 'tv_cluster': cluster,
        })
        
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=data, clustervar="tv_cluster", biters=50,
        )
        with pytest.raises(ValueError, match="varies over time"):
            obj.fit(est_method='dr')
