"""Tests for notyettreated bug fix (Phase 2.1 of R sync).

The R `did` v2.5.1 fixed a bug where control_group='notyettreated' with no
never-treated group would accidentally delete the last-treated cohort
(which serves as the comparison group), biasing remaining estimates or
producing NAs.
"""
import pytest
import numpy as np
import pandas as pd
from csdid.att_gt import ATTgt


@pytest.fixture
def notyettreated_data():
    """Panel data with NO never-treated group.
    
    Groups: treated at 2005, treated at 2006, treated at 2007.
    The 2007 cohort should serve as the not-yet-treated comparison
    for the 2005 and 2006 cohorts.
    """
    np.random.seed(123)
    n_per_group = 30
    n_periods = 5  # 2003, 2004, 2005, 2006, 2007
    
    groups_unit = np.concatenate([
        np.full(n_per_group, 2005),  # treated at 2005
        np.full(n_per_group, 2006),  # treated at 2006
        np.full(n_per_group, 2007),  # treated at 2007 (last cohort = comparison)
    ])
    n_units = len(groups_unit)
    
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(2003, 2003 + n_periods), n_units)
    groups = np.repeat(groups_unit, n_periods)
    
    # Generate outcomes with treatment effect
    y = np.random.randn(n_units * n_periods)
    treated = (times >= groups) & (groups > 0)
    y[treated] += 1.0  # positive treatment effect
    
    return pd.DataFrame({
        'id': ids,
        'year': times,
        'group': groups,
        'y': y,
    })


class TestNotyettreatedLastCohort:
    def test_no_crash_notyettreated_no_nevertreated(self, notyettreated_data):
        """Should not crash when there's no never-treated group."""
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=notyettreated_data,
            control_group="notyettreated",
            biters=100,
        )
        result = obj.fit(est_method='dr')
        assert result is not None
        assert len(result.results['att']) > 0

    def test_last_cohort_retained_in_data(self, notyettreated_data):
        """The last-treated cohort (2007) should remain in the data as control."""
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=notyettreated_data,
            control_group="notyettreated",
            biters=100,
        )
        # Check that 2007 group units still exist in preprocessed data
        processed_data = obj.dp['data']
        assert (processed_data['group'] == 2007).any(), \
            "Last-treated cohort (2007) should be retained in data as comparison group"

    def test_last_cohort_not_in_glist(self, notyettreated_data):
        """The last-treated cohort should NOT be in glist (not an estimated group)."""
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=notyettreated_data,
            control_group="notyettreated",
            biters=100,
        )
        glist = obj.dp['glist']
        assert 2007 not in glist, \
            "Last-treated cohort (2007) should not be in estimation groups"

    def test_estimates_not_nan(self, notyettreated_data):
        """ATT estimates should not be all NaN."""
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=notyettreated_data,
            control_group="notyettreated",
            biters=100,
        )
        result = obj.fit(est_method='dr')
        att_values = np.array(result.results['att'])
        assert not np.all(np.isnan(att_values)), \
            "ATT estimates should not be all NaN when last cohort serves as control"

    def test_nevertreated_coerces_without_control(self, notyettreated_data):
        """R did v2.5.1: control_group='nevertreated' with no never-treated group
        warns and coerces the last cohort to never-treated (it no longer errors)."""
        with pytest.warns(UserWarning, match="No never-treated group"):
            obj = ATTgt(
                yname="y", tname="year", idname="id", gname="group",
                data=notyettreated_data,
                control_group="nevertreated",
                biters=100,
            )
        result = obj.fit(est_method='dr')
        # Coercion should yield at least some finite ATT estimates.
        att_values = np.array(result.results['att'])
        assert np.isfinite(att_values).any(), \
            "Coerced never-treated run should produce finite ATT estimates"

    def test_positive_treatment_effect_detected(self, notyettreated_data):
        """With a true effect of 1.0, post-treatment ATTs should be positive on average."""
        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=notyettreated_data,
            control_group="notyettreated",
            biters=200,
        )
        result = obj.fit(est_method='dr')
        att_values = np.array(result.results['att'])
        post_mask = np.array(result.results['post']) == 1
        post_atts = att_values[post_mask]
        # With n=30/group and effect=1.0, mean should be clearly positive
        assert np.nanmean(post_atts) > 0.3, \
            f"Post-treatment ATTs should detect positive effect, got mean={np.nanmean(post_atts):.3f}"
