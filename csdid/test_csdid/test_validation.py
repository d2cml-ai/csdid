"""Tests for input validation in csdid package (Phase 1 of R sync)."""
import pytest
import numpy as np
import pandas as pd
from csdid.att_gt import ATTgt


@pytest.fixture
def sample_data():
    """Create minimal valid panel data for testing."""
    np.random.seed(42)
    n_units = 50
    n_periods = 4
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(2003, 2003 + n_periods), n_units)
    # Groups: 20 units treated at 2004, 15 at 2006, 15 never-treated (0)
    groups_unit = np.concatenate([
        np.full(20, 2004),
        np.full(15, 2006),
        np.full(15, 0)
    ])
    groups = np.repeat(groups_unit, n_periods)
    y = np.random.randn(n_units * n_periods) + 0.5 * (times >= groups) * (groups > 0)
    
    return pd.DataFrame({
        'id': ids,
        'year': times,
        'group': groups,
        'y': y,
    })


class TestColumnValidation:
    def test_missing_yname(self, sample_data):
        with pytest.raises(ValueError, match="Column.*not found"):
            ATTgt(yname="nonexistent", tname="year", idname="id",
                  gname="group", data=sample_data)

    def test_missing_tname(self, sample_data):
        with pytest.raises(ValueError, match="Column.*not found"):
            ATTgt(yname="y", tname="nonexistent", idname="id",
                  gname="group", data=sample_data)

    def test_missing_idname(self, sample_data):
        with pytest.raises(ValueError, match="Column.*not found"):
            ATTgt(yname="y", tname="year", idname="nonexistent",
                  gname="group", data=sample_data)

    def test_missing_gname(self, sample_data):
        with pytest.raises(ValueError, match="Column.*not found"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="nonexistent", data=sample_data)

    def test_missing_weights_name(self, sample_data):
        with pytest.raises(ValueError, match="Column.*not found"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=sample_data, weights_name="nonexistent")

    def test_missing_clustervar(self, sample_data):
        with pytest.raises(ValueError, match="Column.*not found"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=sample_data, clustervar="nonexistent")


class TestReservedNames:
    def test_reserved_column_w(self, sample_data):
        data = sample_data.copy()
        data['w'] = 1.0
        with pytest.raises(ValueError, match="conflict with names used internally"):
            ATTgt(yname="w", tname="year", idname="id",
                  gname="group", data=data)

    def test_reserved_column_rowid(self, sample_data):
        data = sample_data.copy()
        data['rowid'] = range(len(data))
        with pytest.raises(ValueError, match="conflict with names used internally"):
            ATTgt(yname="y", tname="year", idname="rowid",
                  gname="group", data=data)


class TestNumericValidation:
    def test_non_numeric_outcome(self, sample_data):
        data = sample_data.copy()
        data['y_str'] = 'a'
        with pytest.raises(ValueError, match="must be numeric"):
            ATTgt(yname="y_str", tname="year", idname="id",
                  gname="group", data=data)

    def test_non_numeric_tname(self, sample_data):
        data = sample_data.copy()
        data['year_str'] = data['year'].astype(str)
        with pytest.raises(ValueError, match="must be numeric"):
            ATTgt(yname="y", tname="year_str", idname="id",
                  gname="group", data=data)

    def test_non_numeric_gname(self, sample_data):
        data = sample_data.copy()
        data['group_str'] = data['group'].astype(str)
        with pytest.raises(ValueError, match="must be numeric"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group_str", data=data)


class TestNegativeGname:
    def test_negative_gname_rejected(self, sample_data):
        data = sample_data.copy()
        data.loc[data['group'] == 0, 'group'] = -1
        with pytest.raises(ValueError, match="negative values"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=data)


class TestDuplicateRows:
    def test_duplicate_id_time(self, sample_data):
        # Add a duplicate row
        dup_row = sample_data.iloc[[0]].copy()
        data = pd.concat([sample_data, dup_row], ignore_index=True)
        with pytest.raises(ValueError, match="duplicate.*idname.*tname"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=data)


class TestParameterValidation:
    def test_alp_too_low(self, sample_data):
        with pytest.raises(ValueError, match="alp.*must be.*between 0 and 1"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=sample_data, alp=0)

    def test_alp_too_high(self, sample_data):
        with pytest.raises(ValueError, match="alp.*must be.*between 0 and 1"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=sample_data, alp=1.5)

    def test_biters_negative(self, sample_data):
        with pytest.raises(ValueError, match="biters.*must be.*positive"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=sample_data, biters=-10)

    def test_biters_zero(self, sample_data):
        with pytest.raises(ValueError, match="biters.*must be.*positive"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=sample_data, biters=0)

    def test_negative_anticipation(self, sample_data):
        with pytest.raises(ValueError, match="anticipation.*non-negative"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=sample_data, anticipation=-1)

    def test_invalid_control_group(self, sample_data):
        with pytest.raises(ValueError, match="control_group.*must be one of"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=sample_data, control_group="invalid")


class TestWeightsValidation:
    def test_negative_weights(self, sample_data):
        data = sample_data.copy()
        data['wt'] = 1.0
        data.iloc[0, data.columns.get_loc('wt')] = -1.0
        with pytest.raises(ValueError, match="negative values"):
            ATTgt(yname="y", tname="year", idname="id",
                  gname="group", data=data, weights_name="wt")


class TestValidDataPasses:
    def test_valid_data_no_error(self, sample_data):
        """Ensure valid data doesn't raise during initialization."""
        obj = ATTgt(yname="y", tname="year", idname="id",
                    gname="group", data=sample_data)
        assert obj.dp is not None

    def test_valid_data_with_fit(self, sample_data):
        """Ensure valid data can complete fit."""
        obj = ATTgt(yname="y", tname="year", idname="id",
                    gname="group", data=sample_data, biters=100)
        result = obj.fit(est_method='dr')
        assert result.results is not None
        assert len(result.results['att']) > 0
