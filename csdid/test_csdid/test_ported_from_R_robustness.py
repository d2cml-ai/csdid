"""Ported from R `did` tests/testthat/ (bcallaway11/did @ 9aba07d, v2.5.1).

Self-contained ports (no R at runtime) of three R regression intents that the
Python suite did not already cover:

  * test-mutation-safety.R      -> att_gt must not mutate the caller's DataFrame.
  * test-aggte-edge-coverage.R  -> single-group dynamic aggregation respects the
                                   min_e/max_e window and len(att_egt)==len(egt).
  * test-robustness-guards.R    -> gname=+Inf is a valid never-treated sentinel
                                   (preserved, results identical to coding it 0).

Data is regenerated in numpy (R's RNG stream is irrelevant to these structural
invariants); expected values are analytic / structural, not R-RNG-baked.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt


def _panel(seed=1, cohorts=(2, 3, 0), n_per=100, periods=4, te=2.0, drift=0.5):
    rng = np.random.default_rng(seed)
    rows = []
    uid = 0
    for g in cohorts:
        for _ in range(n_per):
            uid += 1
            fe = rng.normal(0, 1)
            for t in range(1, periods + 1):
                post = 1 if (g > 0 and t >= g) else 0
                y = fe + drift * t + te * post + rng.normal(0, 0.3)
                rows.append((uid, t, g, y))
    return pd.DataFrame(rows, columns=["id", "period", "G", "Y"])


# --- test-mutation-safety.R --------------------------------------------------
class TestMutationSafety:
    """R: att_gt must not modify the caller's data (columns or values)."""

    def test_caller_dataframe_unchanged_panel(self):
        df = _panel(seed=3, cohorts=(2, 3, 0))
        snap = df.copy(deep=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ATTgt(yname="Y", tname="period", idname="id", gname="G", data=df,
                  xformla="Y~1", control_group="nevertreated").fit(est_method="reg")
        assert list(df.columns) == list(snap.columns)
        pd.testing.assert_frame_equal(df, snap)

    def test_caller_dataframe_unchanged_rcs(self):
        df = _panel(seed=4, cohorts=(2, 3, 0))
        snap = df.copy(deep=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ATTgt(yname="Y", tname="period", idname="id", gname="G", data=df,
                  xformla="Y~1", control_group="nevertreated",
                  panel=False).fit(est_method="reg")
        assert list(df.columns) == list(snap.columns)
        pd.testing.assert_frame_equal(df, snap)


# --- test-aggte-edge-coverage.R (block 2) ------------------------------------
class TestSingleGroupDynamicWindow:
    """R: single treated cohort, dynamic aggregation with min_e/max_e must return
    an event-time set fully within [min_e, max_e] with len(att_egt)==len(egt)."""

    def test_dynamic_window_bounded_and_aligned(self):
        df = _panel(seed=7, cohorts=(3, 0), n_per=200, periods=5, te=2.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(yname="Y", tname="period", idname="id", gname="G",
                        data=df.copy(), xformla="Y~1",
                        control_group="nevertreated").fit(est_method="reg")
            agg = res.aggte(typec="dynamic", min_e=-1, max_e=1)
        at = agg.atte
        egt = np.asarray(at["egt"])
        att_egt = np.asarray(at["att_egt"])
        assert len(egt) == len(att_egt)
        assert egt.min() >= -1
        assert egt.max() <= 1
        # the on-impact effect (e=0) should recover the simulated te (~2)
        e0 = att_egt[np.where(egt == 0)[0][0]]
        assert abs(e0 - 2.0) < 0.3


# --- test-robustness-guards.R (Inf gname blocks) -----------------------------
class TestInfGnameAsNeverTreated:
    """R: gname=+Inf is a valid never-treated sentinel: such units are preserved
    (not dropped) and yield results identical to coding them as 0."""

    def _att_keyed(self, df):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(yname="Y", tname="period", idname="id", gname="G",
                        data=df.copy(), xformla="Y~1",
                        control_group="nevertreated").fit(est_method="reg").results
        return {
            f"{int(g)}_{int(t)}": float(a)
            for g, t, a in zip(res["group"], res["year"], res["att"])
        }

    def test_inf_gname_equals_zero_coding(self):
        df0 = _panel(seed=5, cohorts=(2, 3, 0), n_per=150, periods=4)
        df_inf = df0.copy()
        df_inf["G"] = df_inf["G"].astype(float)
        df_inf.loc[df_inf.G == 0, "G"] = np.inf
        a0 = self._att_keyed(df0)
        ainf = self._att_keyed(df_inf)
        # never-treated units must NOT become an estimated group under either coding
        assert "0_2" not in a0 and "inf_2" not in {k for k in ainf}
        # identical estimated cells and identical ATT values
        assert set(a0) == set(ainf)
        for k in a0:
            if not np.isnan(a0[k]):
                assert abs(a0[k] - ainf[k]) <= 1e-10, f"cell {k}: {a0[k]} vs {ainf[k]}"
