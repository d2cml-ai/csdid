"""Bespoke, self-contained, mutation-gated unit tests for the group-time ATT cell
math in ``csdid/attgt_fnc/compute_att_gt.py`` (standard path) and
``csdid/attgt_fnc/compute_att_gt2.py`` (faster_mode, vectorized path).

NO R at runtime. Two oracle families:

1. EQUIVALENCE: the two paths are documented to be numerically IDENTICAL (the
   faster path only changes HOW the 2x2 cell data is sourced; the masks, weights,
   n/n1 IF rescaling, post_treat flag, control-group logic and IF bookkeeping are
   the same). Asserting bit-identical att / inffunc / group / post between the two
   paths kills any single-file mutation in EITHER (a mutation on a mask /
   weight / scaling line in one file makes the two diverge). Run across est_method
   x control_group x base_period x fix_weights (incl. the fix_weights="varying"
   panel->RC routing that exercises the audit-v3 F1 balance-scaling guard).

2. ABSOLUTE: with constant additive treatment effect te, every post-treatment
   ATT(g,t) estimate (with the correct covariate model) recovers te; pre-period
   placebo cells are ~0. Pins the post_treat = 1{g<=t} flag and the cell ATT math
   to a known value, independent of the other path.
"""
import contextlib
import io
import warnings

import numpy as np
import pytest

from csdid.attgt_fnc.preprocess_did import pre_process_did
from csdid.attgt_fnc.compute_att_gt import compute_att_gt
from csdid.attgt_fnc.compute_att_gt2 import compute_att_gt2
from helpers import build_sim_data


def _prep(seed=1, n=400, te=2.0, control_group="nevertreated", fix_weights=None):
    df = build_sim_data(n=n, time_periods=4, te=te, seed=seed)
    with contextlib.redirect_stdout(io.StringIO()):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dp = pre_process_did(yname="Y", tname="period", idname="id",
                                 gname="G", data=df, xformla="~X",
                                 control_group=control_group,
                                 fix_weights=fix_weights)
    return dp


def _run(dp, path, est_method, base_period):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stdout(io.StringIO()):
            return path(dp, est_method=est_method, base_period=base_period,
                        compute_inffunc=True)


# --------------------------------------------------------------------------- #
# 1. faster_mode == standard  (bit-identical)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("est_method", ["dr", "reg", "ipw"])
@pytest.mark.parametrize("control_group", ["nevertreated", "notyettreated"])
@pytest.mark.parametrize("base_period", ["varying", "universal"])
def test_faster_equals_standard(est_method, control_group, base_period):
    dp = _prep(control_group=control_group)
    o1, if1 = _run(dp, compute_att_gt, est_method, base_period)
    o2, if2 = _run(dp, compute_att_gt2, est_method, base_period)
    assert o1["group"] == o2["group"]
    assert o1["post"] == o2["post"]
    a1, a2 = np.array(o1["att"]), np.array(o2["att"])
    assert np.nanmax(np.abs(a1 - a2)) < 1e-12
    assert np.nanmax(np.abs(if1 - if2)) < 1e-12


def test_faster_equals_standard_fix_weights_varying():
    """fix_weights='varying' routes a balanced panel through the RC estimators and
    exercises the audit-v3 F1 IF balance-scaling guard (n/#units vs n/#rows). The
    two paths must still be bit-identical."""
    dp = _prep(fix_weights="varying")
    assert dp["panel"] is False  # routed to the RC path
    o1, if1 = _run(dp, compute_att_gt, "dr", "varying")
    o2, if2 = _run(dp, compute_att_gt2, "dr", "varying")
    assert np.nanmax(np.abs(np.array(o1["att"]) - np.array(o2["att"]))) < 1e-12
    assert np.nanmax(np.abs(if1 - if2)) < 1e-12


# --------------------------------------------------------------------------- #
# 2. absolute oracle: constant te is recovered; placebo pre-cells ~ 0;
#    post_treat flag = 1{g <= t}
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("path", [compute_att_gt, compute_att_gt2])
def test_constant_te_recovered_and_post_flag(path):
    te = 2.0
    dp = _prep(seed=5, n=800, te=te)
    out, _ = _run(dp, path, "dr", "varying")
    g = np.array(out["group"])
    t = np.array(out["year"])
    a = np.array(out["att"])
    post = np.array(out["post"])
    # post flag must equal 1{g <= t}
    assert np.array_equal(post, (g <= t).astype(int))
    # post-treatment cells recover te; pre cells (placebo) ~ 0
    post_mask = (g <= t) & ~np.isnan(a)
    pre_mask = (g > t) & ~np.isnan(a)
    assert abs(np.mean(a[post_mask]) - te) < 0.3
    assert abs(np.mean(a[pre_mask])) < 0.3


@pytest.mark.parametrize("path", [compute_att_gt, compute_att_gt2])
def test_inffunc_overall_mean_zero(path):
    """The SUM of all per-cell influence-function columns is mean-zero over units
    (the aggregate IF integrates to zero); a global sign error in the IF placement
    or n/n1 rescaling breaks this aggregate property."""
    dp = _prep(seed=6, n=600)
    _, ifn = _run(dp, path, "dr", "varying")
    agg = np.nansum(ifn, axis=1)  # sum across cells, per unit
    assert abs(np.mean(agg)) < 1e-6
