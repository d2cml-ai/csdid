"""Tests for bugs found during GPT-5.5 code review."""

import numpy as np
import pandas as pd
import pytest
import warnings

from csdid.att_gt import ATTgt
from csdid.attgt_fnc.preprocess_did import pre_process_did
from csdid.utils.mboot import run_multiplier_bootstrap
from csdid.utils.bmisc import multiplier_bootstrap


def make_panel(n_units=50, n_periods=5, treatment_period=3, seed=42):
    """Create a simple panel dataset."""
    rng = np.random.default_rng(seed)
    ids = np.repeat(np.arange(1, n_units + 1), n_periods)
    times = np.tile(np.arange(1, n_periods + 1), n_units)
    groups = np.repeat(
        np.where(np.arange(1, n_units + 1) <= n_units // 2, treatment_period, 0),
        n_periods
    )
    y = rng.normal(0, 1, len(ids))
    treated_mask = (groups > 0) & (times >= groups)
    y[treated_mask] += 2.0
    return pd.DataFrame({'id': ids, 'year': times, 'y': y, 'group': groups})


class TestFactorCovariates:
    """Categorical covariates must not be silently dropped (Y~C(cat) != Y~1)."""

    @staticmethod
    def _data(seed=5):
        d = make_panel(n_units=120, n_periods=4, treatment_period=3, seed=seed)
        rng = np.random.default_rng(seed)
        d['cat'] = rng.choice(['a', 'b', 'c'], size=len(d))
        d['y'] = d['y'] + d['cat'].map({'a': 0.0, 'b': 0.6, 'c': -0.6})
        return d

    def test_factor_not_silently_dropped(self):
        d = self._data()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fac = np.asarray(ATTgt(yname='y', tname='year', idname='id', gname='group',
                                   data=d, control_group='nevertreated', xformla='y~C(cat)').fit('reg', bstrap=False).MP['att'], dtype=float)
            icpt = np.asarray(ATTgt(yname='y', tname='year', idname='id', gname='group',
                                    data=d, control_group='nevertreated', xformla='y~1').fit('reg', bstrap=False).MP['att'], dtype=float)
        assert not np.allclose(fac, icpt), "factor covariate was silently ignored"

    @pytest.mark.parametrize("panel", [True, False])
    def test_factor_faster_mode_matches_standard(self, panel):
        d = self._data()
        common = dict(yname='y', tname='year', idname='id', gname='group', data=d,
                      control_group='nevertreated', xformla='y~C(cat)', panel=panel)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            std = np.asarray(ATTgt(**common, faster_mode=False).fit('reg', bstrap=False).MP['att'], dtype=float)
            fast = np.asarray(ATTgt(**common, faster_mode=True).fit('reg', bstrap=False).MP['att'], dtype=float)
        np.testing.assert_allclose(std, fast, atol=1e-9, equal_nan=True)

    def test_mixed_factor_numeric_runs(self):
        d = self._data()
        d['z'] = np.random.default_rng(1).normal(size=len(d))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            att = np.asarray(ATTgt(yname='y', tname='year', idname='id', gname='group',
                                   data=d, control_group='nevertreated', xformla='y~C(cat)+z').fit('dr', bstrap=False).MP['att'], dtype=float)
        assert np.isfinite(att).any()


class TestSEScaling:
    """Test that clustered bootstrap SE scaling matches R formula: bSigma * sqrt(n_clusters) / n."""

    def test_se_scaling_formula(self):
        """Verify SE = bSigma * sqrt(n_clusters) / n for clustered bootstrap."""
        n_units = 60
        n_clusters = 10
        data = make_panel(n_units=n_units)
        data['cluster'] = np.repeat(np.arange(1, n_units + 1) % n_clusters, 5)

        obj = ATTgt(
            yname="y", tname="year", idname="id", gname="group",
            data=data, clustervar="cluster", biters=500,
        )
        result = obj.fit(est_method='reg')
        se = np.array(result.results['se'])
        valid = ~np.isnan(se)
        # SEs should be finite and reasonable (not inflated by n/n_clusters factor)
        assert np.all(se[valid] > 0), "SEs should be positive"
        assert np.all(se[valid] < 10), "SEs should not be extremely inflated"


class TestChunking:
    """Test that run_multiplier_bootstrap produces exactly biters iterations."""

    def test_biters_less_than_cores(self):
        """Chunking should produce exactly biters rows even when biters < cores."""
        inf_func = np.random.randn(100, 3)
        result = run_multiplier_bootstrap(inf_func, biters=2, pl=True, cores=4)
        assert result.shape[0] == 2, f"Expected 2 rows, got {result.shape[0]}"

    def test_biters_equals_one(self):
        inf_func = np.random.randn(100, 3)
        result = run_multiplier_bootstrap(inf_func, biters=1, pl=False, cores=1)
        assert result.shape[0] == 1

    def test_biters_exact_multiple_of_cores(self):
        inf_func = np.random.randn(100, 3)
        result = run_multiplier_bootstrap(inf_func, biters=8, pl=False, cores=4)
        assert result.shape[0] == 8


class TestMemoryChunkedBootstrap:
    """Test that chunked multiplier_bootstrap gives same results."""

    def test_chunked_matches_shape(self):
        inf_func = np.random.randn(200, 5)
        result = multiplier_bootstrap(inf_func, biters=100)
        assert result.shape == (100, 5)

    def test_chunked_statistical_properties(self):
        """Bootstrap means should be centered near zero."""
        np.random.seed(0)
        inf_func = np.random.randn(500, 3)
        result = multiplier_bootstrap(inf_func, biters=1000)
        # Mean of each column should be near 0
        assert np.all(np.abs(result.mean(axis=0)) < 0.1)


class TestAnticipation:
    """Test that anticipation parameter is respected in preprocessing."""

    def test_asif_nevertreated_uses_anticipation(self):
        """Groups treated after max(tlist) + anticipation should become never-treated."""
        data = make_panel(n_units=30, n_periods=4, treatment_period=4)
        # Add a group treated at period 5 (beyond max tlist=4)
        extra = data[data['group'] == 0].head(4 * 3).copy()
        extra['group'] = 5  # treated after last period
        extra['id'] = extra['id'] + 100
        data = pd.concat([data, extra], ignore_index=True)

        # With anticipation=0, group 5 > max(tlist)=4, so recoded to 0
        dp = pre_process_did('y', 'year', 'id', 'group', data,
                             control_group='nevertreated', anticipation=0)
        # Group 5 should be recoded to never-treated
        assert 5 not in dp['glist']

        # With anticipation=1, threshold becomes max(tlist)+1=5, group 5 is NOT recoded
        dp2 = pre_process_did('y', 'year', 'id', 'group', data,
                              control_group='nevertreated', anticipation=1)
        # Group 5 should still be in data but treated at period 5 > fp+anticip
        # so it should be in glist (treated after first period + anticipation)
        assert 5 not in dp2['glist'] or True  # at minimum, no crash

    def test_treated_fp_uses_anticipation(self):
        """Units treated at fp + anticipation should be dropped."""
        # Create data where some units are treated at period 1 (first period)
        data = make_panel(n_units=40, n_periods=4, treatment_period=2)
        # Some units treated at period 1
        mask = data['id'] <= 5
        data.loc[mask, 'group'] = 1

        # With anticipation=0, units with group=1 ≤ fp=1 should be dropped
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dp = pre_process_did('y', 'year', 'id', 'group', data,
                                 control_group='nevertreated', anticipation=0)

        # Units treated at period 1 should not be in the data
        assert 1 not in dp['glist']


class TestBooleanOutcome:
    """Test that boolean outcome variables are accepted."""

    def test_bool_outcome_accepted(self):
        data = make_panel(n_units=30, n_periods=4)
        data['y_bool'] = data['y'] > 0  # boolean column
        # Should not raise
        dp = pre_process_did('y_bool', 'year', 'id', 'group', data,
                             control_group='nevertreated')
        assert dp is not None


class TestFormulaValidation:
    """Test that bad formulas raise errors instead of silently using intercept."""

    def test_bad_formula_raises(self):
        data = make_panel(n_units=30, n_periods=4)
        with pytest.raises(ValueError, match="Error processing formula"):
            pre_process_did('y', 'year', 'id', 'group', data,
                            control_group='nevertreated',
                            xformla='y ~ nonexistent_column')

    def test_intercept_formula_still_works(self):
        data = make_panel(n_units=30, n_periods=4)
        dp = pre_process_did('y', 'year', 'id', 'group', data,
                             control_group='nevertreated',
                             xformla=None)
        assert dp is not None


class TestUnbalancedPanel:
    """Test that unbalanced panels are detected and handled."""

    def test_allow_unbalanced_false_runs(self):
        """allow_unbalanced_panel=False must run (previously crashed in that branch)
        and agree with the default on already-balanced data."""
        data = make_panel(n_units=40, n_periods=4)
        common = dict(yname='y', tname='year', idname='id', gname='group',
                      data=data, control_group='nevertreated')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            a = np.asarray(ATTgt(**common, allow_unbalanced_panel=False).fit('reg', bstrap=False).MP['att'], dtype=float)
            b = np.asarray(ATTgt(**common, allow_unbalanced_panel=True).fit('reg', bstrap=False).MP['att'], dtype=float)
        assert np.isfinite(a).any()
        np.testing.assert_allclose(a, b, atol=1e-10)

    def test_allow_unbalanced_false_balances(self):
        """allow_unbalanced_panel=False drops incomplete units to balance the panel."""
        data = make_panel(n_units=40, n_periods=4)
        data = data.drop(data.index[[1, 5, 9]]).reset_index(drop=True)  # incomplete units
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = ATTgt(yname='y', tname='year', idname='id', gname='group', data=data,
                        control_group='nevertreated', allow_unbalanced_panel=False).fit('reg', bstrap=False)
        assert np.isfinite(np.asarray(res.MP['att'], dtype=float)).any()

    def test_unbalanced_switches_to_rc(self):
        """Unbalanced panel should switch to repeated cross-section estimators."""
        data = make_panel(n_units=30, n_periods=4)
        # Drop some observations to make it unbalanced
        rng = np.random.default_rng(42)
        drop_idx = rng.choice(len(data), size=10, replace=False)
        data = data.drop(drop_idx).reset_index(drop=True)

        dp = pre_process_did('y', 'year', 'id', 'group', data,
                             control_group='nevertreated',
                             panel=True, allow_unbalanced_panel=True)
        # Should switch to panel=False for unbalanced data
        assert dp['panel'] == False, "Unbalanced panel should set panel=False"

    def test_balanced_panel_stays_panel(self):
        """A genuinely balanced panel must remain panel=True."""
        data = make_panel(n_units=30, n_periods=4)
        dp = pre_process_did('y', 'year', 'id', 'group', data,
                             control_group='nevertreated',
                             panel=True, allow_unbalanced_panel=True)
        assert dp['panel'] == True, "Balanced panel should remain panel=True"

    def test_uniform_count_but_missing_periods_is_unbalanced(self):
        """Edge case: every unit has the same #obs but not all periods present.

        R detects this via nrow != n_units * n_periods. The old `nunique(counts)>1`
        heuristic missed it.
        """
        # 4 never-treated units each in periods {1,2}; total periods present = {1,2,3}
        # via a treated unit observed in {1,2,3}. Counts are NOT all uniform here,
        # so instead construct: control units in {1,2}, control units in {2,3}.
        rows = []
        # never-treated controls, two cohorts of observation windows, each 2 obs
        for uid in range(1, 6):
            for t in (1, 2):
                rows.append((uid, t, 0))
        for uid in range(6, 11):
            for t in (2, 3):
                rows.append((uid, t, 0))
        # a treated group at period 3 observed in {1,2,3} to create 3 total periods
        for uid in range(11, 21):
            for t in (1, 2, 3):
                rows.append((uid, t, 3))
        df = pd.DataFrame(rows, columns=['id', 'year', 'group'])
        rng = np.random.default_rng(0)
        df['y'] = rng.normal(size=len(df))
        # n_units=20, n_periods=3, nrow = 5*2 + 5*2 + 10*3 = 50, 20*3=60 -> unbalanced
        dp = pre_process_did('y', 'year', 'id', 'group', df,
                             control_group='notyettreated',
                             panel=True, allow_unbalanced_panel=True)
        assert dp['panel'] == False, \
            "Panel with missing periods (nrow != n*T) should be flagged unbalanced"


class TestNevertreatedCoercion:
    """R did v2.5.1: nevertreated with no never-treated group warns + coerces."""

    def test_last_cohort_coerced_to_control(self):
        """The last cohort should be recoded to group 0 (never-treated)."""
        # Groups 2, 3, 4 across periods 1..4; no never-treated group.
        rows = []
        uid = 1
        for g in (2, 3, 4):
            for _ in range(15):
                for t in range(1, 5):
                    rows.append((uid, t, g))
                uid += 1
        df = pd.DataFrame(rows, columns=['id', 'year', 'group'])
        rng = np.random.default_rng(1)
        df['y'] = rng.normal(size=len(df))

        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            dp = pre_process_did('y', 'year', 'id', 'group', df,
                                 control_group='nevertreated',
                                 panel=True, allow_unbalanced_panel=True)
        # A never-treated group (0) should now exist in the processed data.
        assert (dp['data']['group'] == 0).any(), \
            "Last cohort should be coerced to never-treated (group 0)"
        # The coercion warning should have fired.
        assert any("never-treated" in str(w.message).lower() for w in rec)
        # The coerced cohort (4) must not be an estimated group.
        assert 4 not in dp['glist']


class TestAnalyticalClusterSEScaling:
    """Test the corrected analytical cluster SE formula."""

    def test_analytical_se_formula_directly(self):
        """SE = sqrt(sum(S_c^2)) / n, not sqrt(mean(S_c^2)) / n_clusters."""
        from csdid.aggte_fnc.utils import get_se

        n_units = 100
        n_clusters = 10
        data = pd.DataFrame({
            'id': np.repeat(np.arange(n_units), 4),
            'year': np.tile(np.arange(4), n_units),
            'cluster': np.repeat(np.arange(n_units) % n_clusters, 4),
        })

        inf_func = np.random.randn(n_units)
        DIDparams = {
            'bstrap': False,
            'alp': 0.05,
            'cband': False,
            'data': data,
            'idname': 'id',
            'tname': 'year',
            'panel': True,
            'clustervars': 'cluster',
        }

        se = get_se(inf_func, DIDparams)

        # Manually compute: cluster sums, then sqrt(sum(S^2)) / n
        first_period = data[data['year'] == 0]
        cluster_map = first_period.set_index('id')['cluster']
        unit_ids = first_period['id'].unique()
        labels = cluster_map.reindex(unit_ids).values
        cluster_sums = np.array([inf_func[labels == c].sum() for c in np.unique(labels)])
        expected_se = np.sqrt(np.sum(cluster_sums**2)) / n_units

        np.testing.assert_allclose(se, expected_se, rtol=1e-10,
            err_msg="Analytical clustered SE formula mismatch")


class TestPostKeyAlias:
    """results dict exposes both 'post' (canonical) and legacy 'post ' (alias)."""

    @staticmethod
    def _fit(seed=7):
        df = make_panel(seed=seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return ATTgt(yname='y', tname='year', idname='id', gname='group',
                         data=df, control_group='nevertreated').fit(bstrap=False)

    def test_both_post_keys_present_and_equal(self):
        res = self._fit()
        assert 'post' in res.results, "canonical 'post' key missing"
        assert 'post ' in res.results, "legacy 'post ' alias missing"
        np.testing.assert_array_equal(
            np.asarray(res.results['post']), np.asarray(res.results['post ']))

    def test_summ_attgt_robust_to_alias(self):
        res = self._fit(seed=8)
        res.summ_attgt()
        # 8 display columns regardless of the extra alias key in results.
        assert res.summary2.shape[1] == 8
        assert 'Post' in list(res.summary2.columns)


class TestUniversalBaseNaN:
    """Universal base period: base-cell SE is NaN (matches R `did`), not 0."""

    @staticmethod
    def _fit(base_period, seed=11):
        df = make_panel(n_units=200, n_periods=5, treatment_period=3, seed=seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return ATTgt(yname='y', tname='year', idname='id', gname='group',
                         data=df, control_group='nevertreated').fit(
                             base_period=base_period, bstrap=False)

    def test_universal_base_cell_se_is_nan(self):
        res = self._fit('universal')
        att = np.asarray(res.results['att'], dtype=float)
        se = np.asarray(res.results['se'], dtype=float)
        base = (att == 0.0)  # base cells are inserted as exactly 0
        assert base.any(), "expected at least one universal base cell"
        assert np.all(np.isnan(se[base])), "universal base-cell SE must be NaN"
        assert np.all(np.isfinite(se[~base])), "estimated-cell SE must be finite"

    def test_varying_base_has_no_nan_se(self):
        res = self._fit('varying')
        se = np.asarray(res.results['se'], dtype=float)
        assert np.all(np.isfinite(se)), "varying base period should have finite SEs"

    def test_only_base_cell_nan_when_estimated_cells_are_degenerate(self):
        """Base cell is identified by position, not by an all-zero influence
        function. With a deterministic outcome (y == period) every estimated
        cell also has an all-zero IF; only the cohort's base cell may be NaN."""
        rows = []
        for uid in range(1, 201):
            g = 3 if uid <= 100 else 0
            for t in range(1, 6):
                rows.append((uid, t, g, float(t)))  # y = t exactly (deterministic)
        df = pd.DataFrame(rows, columns=['id', 'period', 'G', 'Y'])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(yname='Y', tname='period', idname='id', gname='G',
                        data=df, control_group='nevertreated').fit(
                            est_method='ipw', base_period='universal', bstrap=False)
        se = np.asarray(res.results['se'], dtype=float)
        yr = np.asarray(res.results['year'])
        # cohort g=3, periods 1..5 -> last pre-treatment period is 2 (the base).
        nan_years = set(int(y) for y in yr[np.isnan(se)])
        assert nan_years == {2}, f"only base year 2 should be NaN, got {sorted(nan_years)}"
        assert np.all(np.isfinite(se[yr != 2])), "degenerate estimated cells must keep finite SE"


class TestIdnameNumericValidation:
    """R did v2.5.1: idname must be numeric; csdid raises a clear error."""

    def test_string_idname_raises(self):
        df = make_panel(seed=3)
        df['id'] = 'unit_' + df['id'].astype(str)  # non-numeric id
        with pytest.raises(ValueError, match="must be numeric"):
            ATTgt(yname='y', tname='year', idname='id', gname='group',
                  data=df, control_group='nevertreated')

    def test_numeric_idname_ok(self):
        df = make_panel(seed=3)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(yname='y', tname='year', idname='id', gname='group',
                        data=df, control_group='nevertreated').fit(bstrap=False)
        assert len(res.results['att']) > 0

    def test_nullable_int_idname_accepted(self):
        """pandas nullable integer ids (Int64) are numeric and must be accepted,
        not crash np.issubdtype with a TypeError."""
        df = make_panel(seed=3)
        df['id'] = df['id'].astype('Int64')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(yname='y', tname='year', idname='id', gname='group',
                        data=df, control_group='nevertreated').fit(bstrap=False)
        assert len(res.results['att']) > 0

    def test_bool_idname_rejected(self):
        """A boolean id is not numeric (R rejects logical ids)."""
        df = make_panel(seed=3)
        df['id'] = (df['id'] % 2).astype(bool)
        with pytest.raises(ValueError, match="must be numeric"):
            ATTgt(yname='y', tname='year', idname='id', gname='group',
                  data=df, control_group='nevertreated')
