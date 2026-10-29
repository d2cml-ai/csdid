"""Tests for compute_inffunc parameter (Phase 4.1 of R sync).

R did v2.5.1 added compute_inffunc=FALSE for point-estimates-only runs
that skip influence function computation for speed/memory savings.
"""
import pytest
import numpy as np
import pandas as pd
from csdid.att_gt import ATTgt


@pytest.fixture
def sample_data():
    np.random.seed(42)
    n_units, n_periods = 50, 4
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(2003, 2007), n_units)
    groups_unit = np.concatenate([np.full(20, 2004), np.full(15, 2006), np.full(15, 0)])
    groups = np.repeat(groups_unit, n_periods)
    y = np.random.randn(n_units * n_periods) + 0.5 * (times >= groups) * (groups > 0)
    return pd.DataFrame({'id': ids, 'year': times, 'group': groups, 'y': y})


class TestComputeInffunc:
    def test_point_estimates_match(self, sample_data):
        """Point estimates should be identical with and without IF computation."""
        obj_full = ATTgt(yname="y", tname="year", idname="id", gname="group",
                         data=sample_data, biters=100, compute_inffunc=True)
        result_full = obj_full.fit(est_method='dr')
        
        obj_fast = ATTgt(yname="y", tname="year", idname="id", gname="group",
                         data=sample_data, biters=100, compute_inffunc=False)
        result_fast = obj_fast.fit(est_method='dr')
        
        assert np.allclose(result_full.results['att'], result_fast.results['att'])

    def test_se_is_nan_without_inffunc(self, sample_data):
        """SEs should be NaN when compute_inffunc=False."""
        obj = ATTgt(yname="y", tname="year", idname="id", gname="group",
                    data=sample_data, biters=100, compute_inffunc=False)
        result = obj.fit(est_method='dr')
        assert np.all(np.isnan(result.results['se']))

    def test_aggte_blocked_without_inffunc(self, sample_data):
        """aggte() should raise when compute_inffunc=False."""
        obj = ATTgt(yname="y", tname="year", idname="id", gname="group",
                    data=sample_data, biters=100, compute_inffunc=False)
        obj.fit(est_method='dr')
        
        with pytest.raises(ValueError, match="compute_inffunc=False"):
            obj.aggte(typec='simple')

    def test_faster_without_inffunc(self, sample_data):
        """compute_inffunc=False must skip influence-function/bootstrap work.

        Rather than rely on noisy wall-clock timing, assert the deterministic
        contract that makes fast mode cheaper: the influence-function matrix is
        not materialized and no bootstrap-based SEs/critical values are produced.
        """
        obj_fast = ATTgt(yname="y", tname="year", idname="id", gname="group",
                         data=sample_data, biters=500, compute_inffunc=False)
        res_fast = obj_fast.fit(est_method='dr')

        obj_full = ATTgt(yname="y", tname="year", idname="id", gname="group",
                         data=sample_data, biters=500, compute_inffunc=True)
        res_full = obj_full.fit(est_method='dr')

        # Fast mode produces no influence functions (the expensive artifact).
        assert obj_fast.compute_inffunc is False
        assert res_fast.MP['inffunc'] is None

        # Full mode does compute them.
        assert obj_full.compute_inffunc is True
        assert res_full.MP['inffunc'] is not None

        # Fast mode yields no finite bootstrap SEs (it skips the bootstrap entirely).
        se_fast = np.asarray(res_fast.results['se'], dtype=float)
        assert np.all(np.isnan(se_fast)), "Fast mode should not compute standard errors"

    def test_default_is_true(self, sample_data):
        """Default should compute influence functions."""
        obj = ATTgt(yname="y", tname="year", idname="id", gname="group",
                    data=sample_data, biters=100)
        assert obj.compute_inffunc is True
