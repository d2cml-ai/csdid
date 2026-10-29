"""Integration test: full end-to-end workflow."""
import numpy as np
import pandas as pd
import pytest
from csdid.att_gt import ATTgt


@pytest.fixture
def panel_data():
    np.random.seed(42)
    n_units, n_periods = 100, 5
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(2003, 2008), n_units)
    groups_unit = np.concatenate([np.full(30, 2005), np.full(30, 2006), np.full(40, 0)])
    groups = np.repeat(groups_unit, n_periods)
    y = np.random.randn(n_units * n_periods) + 0.8 * (times >= groups) * (groups > 0)
    return pd.DataFrame({'id': ids, 'year': times, 'group': groups, 'y': y})


class TestEndToEnd:
    def test_full_pipeline(self, panel_data):
        """Complete pipeline: ATTgt → fit → aggte for all types."""
        obj = ATTgt(yname="y", tname="year", idname="id", gname="group",
                    data=panel_data, biters=200)
        result = obj.fit(est_method='dr')
        
        # Check ATTs are computed
        att = np.array(result.results['att'])
        assert len(att) == 8  # 2 groups × 4 time periods
        assert np.all(np.isfinite(att))
        
        # Check SEs are computed
        se = np.array(result.results['se'])
        assert np.all(se > 0)
        
        # Run all aggregation types
        for typec in ['simple', 'group', 'dynamic', 'calendar']:
            agg = result.aggte(typec=typec)
            assert agg.atte['overall_att'] is not None
            assert agg.atte['overall_se'] is not None

    def test_treatment_effect_detected(self, panel_data):
        """With true effect=0.8, should detect positive ATT."""
        obj = ATTgt(yname="y", tname="year", idname="id", gname="group",
                    data=panel_data, biters=200)
        result = obj.fit(est_method='dr')
        agg = result.aggte(typec='simple')
        assert agg.atte['overall_att'] > 0.3  # True effect is 0.8

    def test_fast_mode_consistency(self, panel_data):
        """Fast mode produces same ATTs as full mode."""
        obj_full = ATTgt(yname="y", tname="year", idname="id", gname="group",
                         data=panel_data, biters=200, compute_inffunc=True)
        result_full = obj_full.fit(est_method='dr')
        
        obj_fast = ATTgt(yname="y", tname="year", idname="id", gname="group",
                         data=panel_data, compute_inffunc=False)
        result_fast = obj_fast.fit(est_method='dr')
        
        assert np.allclose(result_full.results['att'], result_fast.results['att'])

    def test_all_estimation_methods(self, panel_data):
        """All estimation methods should work."""
        for method in ['dr', 'ipw', 'reg']:
            obj = ATTgt(yname="y", tname="year", idname="id", gname="group",
                        data=panel_data, biters=100)
            result = obj.fit(est_method=method)
            att = np.array(result.results['att'])
            assert np.all(np.isfinite(att)), f"Method '{method}' produced non-finite ATTs"

    def test_notyettreated_control(self, panel_data):
        """notyettreated control group should work."""
        obj = ATTgt(yname="y", tname="year", idname="id", gname="group",
                    data=panel_data, control_group='notyettreated', biters=100)
        result = obj.fit(est_method='dr')
        att = np.array(result.results['att'])
        assert len(att) > 0
        assert np.all(np.isfinite(att))
