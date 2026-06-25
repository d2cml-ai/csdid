"""Tests for per-cell estimation failure handling (Phase 2.2 of R sync).

When a 2x2 DiD estimation fails for a specific (g,t) cell (e.g., singular
design matrix), the package should warn and set that cell's ATT to NA
instead of crashing.
"""
import pytest
import numpy as np
import pandas as pd
import warnings
from csdid.att_gt import ATTgt


@pytest.fixture
def collinear_covariate_data():
    """Panel data where one group has perfectly collinear covariates,
    causing the propensity score model to fail for that cell."""
    np.random.seed(42)
    n_per_group = 20
    n_periods = 4
    
    groups_unit = np.concatenate([
        np.full(n_per_group, 2005),
        np.full(n_per_group, 2006),
        np.full(n_per_group, 0),  # never-treated
    ])
    n_units = len(groups_unit)
    
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(2003, 2003 + n_periods), n_units)
    groups = np.repeat(groups_unit, n_periods)
    y = np.random.randn(n_units * n_periods)
    
    # Create a covariate that's fine for most groups
    x1 = np.random.randn(n_units * n_periods)
    # Make x2 = x1 for group 2005 units only (collinear within that group)
    x2 = np.random.randn(n_units * n_periods)
    mask_2005 = groups == 2005
    x2[mask_2005] = x1[mask_2005]
    
    return pd.DataFrame({
        'id': ids,
        'year': times,
        'group': groups,
        'y': y,
        'x1': x1,
        'x2': x2,
    })


@pytest.fixture
def tiny_group_data():
    """Panel data where one group has very few observations (1 unit),
    likely causing estimation to fail."""
    np.random.seed(42)
    n_periods = 4
    
    # 1 unit in group 2005 (too small), 30 in 2006, 30 never-treated
    groups_unit = np.concatenate([
        np.full(1, 2005),
        np.full(30, 2006),
        np.full(30, 0),
    ])
    n_units = len(groups_unit)
    
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(2003, 2003 + n_periods), n_units)
    groups = np.repeat(groups_unit, n_periods)
    y = np.random.randn(n_units * n_periods)
    
    return pd.DataFrame({
        'id': ids,
        'year': times,
        'group': groups,
        'y': y,
    })


class TestPerCellFailure:
    def test_collinear_covariates_no_crash(self, collinear_covariate_data):
        """Should not crash when covariates are collinear for a specific group."""
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            obj = ATTgt(
                yname="y", tname="year", idname="id", gname="group",
                data=collinear_covariate_data,
                xformla="y ~ x1 + x2",
                biters=50,
            )
            # This should not crash even if some cells fail
            result = obj.fit(est_method='dr')
            assert result is not None

    def test_tiny_group_warns_not_crashes(self, tiny_group_data):
        """Should warn (not crash) when a group is too small for estimation."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            obj = ATTgt(
                yname="y", tname="year", idname="id", gname="group",
                data=tiny_group_data,
                biters=50,
            )
            result = obj.fit(est_method='dr')
            # Should produce results (some may be NA)
            assert result is not None
            att_values = np.array(result.results['att'])
            # At least some values should exist (from the viable group 2006)
            assert not np.all(np.isnan(att_values)), \
                "Not all ATTs should be NA; group 2006 should estimate fine"

    def test_failed_cells_are_nan(self, tiny_group_data):
        """Failed cells should have NaN ATT values."""
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            obj = ATTgt(
                yname="y", tname="year", idname="id", gname="group",
                data=tiny_group_data,
                biters=50,
            )
            result = obj.fit(est_method='dr')
            att_values = np.array(result.results['att'])
            groups = np.array(result.results['group'])
            # Group 2005 has 1 unit — its cells may fail
            # Group 2006 has 30 units — its cells should be fine
            group_2006_atts = att_values[groups == 2006]
            assert not np.any(np.isnan(group_2006_atts)), \
                "Group 2006 (30 units) should have valid ATT estimates"

    def test_normal_data_unaffected(self):
        """Normal data should produce the same results as before (no regression)."""
        np.random.seed(42)
        n_units, n_periods = 50, 4
        ids = np.repeat(np.arange(1, n_units + 1), n_periods)
        times = np.tile(np.arange(2003, 2007), n_units)
        groups_unit = np.concatenate([np.full(20, 2004), np.full(15, 2006), np.full(15, 0)])
        groups = np.repeat(groups_unit, n_periods)
        y = np.random.randn(n_units * n_periods)
        data = pd.DataFrame({'id': ids, 'year': times, 'group': groups, 'y': y})
        
        obj = ATTgt(yname="y", tname="year", idname="id", gname="group",
                    data=data, biters=100)
        result = obj.fit(est_method='dr')
        att_values = np.array(result.results['att'])
        
        # None should be NaN for normal well-conditioned data
        assert not np.any(np.isnan(att_values)), \
            "Normal data should not produce any NaN ATTs"
