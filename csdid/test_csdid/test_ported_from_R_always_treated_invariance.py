"""Ported from R `did` tests/testthat/test-always-treated-invariance.R
(bcallaway11/did @ 9aba07d, v2.5.1).

R INTENT (regression for the latest-control-cohort deletion bug): att_gt output
must be INVARIANT to the presence of effectively-always-treated units (cohorts
with g <= the first period). Dropping those units must not change ATT(g,t) on the
cells common to both runs (including the NA pattern).

PORT NOTES (self-contained; no R at runtime):
- R's data is RNG-specific; we regenerate an equivalent balanced panel in numpy.
  The invariance assertions are DATA-INDEPENDENT (dropping units that should never
  be read must not move the other estimates), so the exact RNG stream is irrelevant.
- The csdid result dict keys point estimates under "att" with group key "group" and
  time key "year".
- Covers the nevertreated invariance + scaling oracle + structural group-selection
  (always-treated cohort dropped, latest cohort kept as control).
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt


# --- helper: build a balanced panel from a vector of cohorts -----------------
# Mirrors R's .mk_design: y = b + 50*(t-1) + 200*(t>=g & g>0) + noise, with a
# large heterogeneous level `b` (so leakage of always-treated units is visible).
def _mk_design(cohorts, n_periods, n=25, seed=1, het=1.2):
    rng = np.random.default_rng(seed)
    recs = []
    uid = 0
    for g in cohorts:
        b = np.exp(rng.normal(np.log(1e5), het, size=n))
        for i in range(n):
            uid += 1
            recs.append((uid, g, b[i]))
    rows = []
    for uid, g, b in recs:
        for t in range(1, n_periods + 1):
            post = 200.0 * (1 if (t >= g and g > 0) else 0)
            y = b + 50.0 * (t - 1) + post + rng.normal(0, 5)
            rows.append((uid, t, g, y))
    return pd.DataFrame(rows, columns=["uid", "t", "g", "y"])


def _att_keyed(data, control_group, drop_uid=(), est_method="reg", **kw):
    """Return {"<g>_<t>": att} like R's .att_keyed."""
    if len(drop_uid):
        data = data[~data["uid"].isin(drop_uid)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="y", tname="t", idname="uid", gname="g", data=data.copy(),
            control_group=control_group, bstrap=False,
        ).fit(est_method=est_method, **kw).results
    return {
        f"{int(g)}_{int(t)}": float(a)
        for g, t, a in zip(res["group"], res["year"], res["att"])
    }


def _assert_invariant(full, drop, tol=1e-8):
    common = set(full) & set(drop)
    assert len(common) > 0
    # NA pattern must match on common cells
    for k in common:
        assert np.isnan(full[k]) == np.isnan(drop[k]), f"NA mismatch at {k}"
    for k in common:
        if not np.isnan(full[k]):
            assert abs(full[k] - drop[k]) <= tol, f"cell {k}: {full[k]} vs {drop[k]}"


# --- nevertreated branch: R says it was always immune -----------------------
def test_nevertreated_invariant_to_always_treated_presence():
    """R: 'control_group=nevertreated is invariant to always-treated presence'."""
    d = _mk_design([0, 1, 2, 3, 4], n_periods=5, seed=2024)
    eat = sorted(d.loc[d.g == 1, "uid"].unique())
    full = _att_keyed(d, "nevertreated")
    drop = _att_keyed(d, "nevertreated", drop_uid=eat)
    _assert_invariant(full, drop)


def test_scaling_oracle_always_treated_outcomes_never_read_nevertreated():
    """R 'scaling oracle': scaling the always-treated cohort's outcomes by 1e6 must
    not change ATT(g,t) for the other groups (those outcomes are never read).
    Ported on the nevertreated branch (where the port honors the oracle)."""
    d = _mk_design([0, 1, 2, 3, 4], n_periods=5, seed=2024)
    d_scaled = d.copy()
    d_scaled.loc[d_scaled.g == 1, "y"] *= 1e6
    a1 = _att_keyed(d, "nevertreated")
    a2 = _att_keyed(d_scaled, "nevertreated")
    other = [k for k in (set(a1) & set(a2)) if not k.startswith("1_")]
    assert other
    for k in other:
        if not np.isnan(a1[k]):
            assert abs(a1[k] - a2[k]) <= 1e-8, f"cell {k} moved: {a1[k]} vs {a2[k]}"


# --- structural: latest cohort kept as control, always-treated dropped ------
def test_latest_cohort_kept_as_control_always_treated_dropped():
    """R: latest cohort retained as control (no ATT of its own); always-treated
    cohort dropped. This STRUCTURAL behavior is correct in the port."""
    d = _mk_design([1, 2, 3, 4, 5], n_periods=5, seed=5)
    keyed = _att_keyed(d, "notyettreated")
    est_groups = sorted({int(k.split("_")[0]) for k in keyed})
    assert 1 not in est_groups       # always-treated cohort dropped
    assert 5 not in est_groups       # latest cohort kept as control, not estimated
    assert est_groups == [2, 3, 4]
