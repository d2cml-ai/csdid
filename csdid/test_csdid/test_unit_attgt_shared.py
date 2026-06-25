"""Bespoke, self-contained, mutation-gated unit tests for the shared cell-selection
and feasibility-guard helpers in ``csdid/attgt_fnc/compute_att_gt_shared.py``.

NO R at runtime. These are pure functions with STRUCTURAL / closed-form oracles:
  * last_pretreatment_index = tail(which(tlist+anticipation < g))  (R formula)
  * plan_cell control flow: ESTIMATE / BASE_ZERO / BREAK with the documented
    pretreatment-index resolution.
  * rcond_check_fail / overlap_check_fail / drdid_design_singular: singular vs
    full-rank designs, the 0.999 overlap threshold, and the dr/reg vs ipw
    design-set distinction (control-only on balanced panel; 4 group*period
    designs on RC; PS full-sample for dr/ipw).
"""
import warnings

import numpy as np
import pytest

from csdid.attgt_fnc import drdid_trim
from drdid import reg_did
from csdid.attgt_fnc.compute_att_gt_shared import (
    last_pretreatment_index,
    plan_cell,
    rcond_check_fail,
    overlap_check_fail,
    drdid_design_singular,
    select_estimators,
    BREAK,
    BASE_ZERO,
    ESTIMATE,
)

TLIST = [1, 2, 3, 4]


# --------------------------------------------------------------------------- #
# last_pretreatment_index  (R: tail(which((tlist+anticipation) < g), 1))
# --------------------------------------------------------------------------- #
def test_lpi_basic():
    # tlist < 3 -> indices {0,1} -> last = 1
    assert last_pretreatment_index(3, TLIST, 0) == 1


def test_lpi_no_pretreatment_returns_none():
    # nothing in tlist < 1
    assert last_pretreatment_index(1, TLIST, 0) is None


def test_lpi_anticipation_shifts():
    # tlist + 1 < 3  ->  tlist < 2  ->  {0}  -> last = 0
    assert last_pretreatment_index(3, TLIST, 1) == 0
    # anticipation makes the cutoff strictly smaller; with a=2, tlist<1 -> none
    assert last_pretreatment_index(3, TLIST, 2) is None


# --------------------------------------------------------------------------- #
# plan_cell : control flow
# --------------------------------------------------------------------------- #
def test_plan_cell_pre_period_estimate():
    # g=3, t_i=0 (tn=2): pre-period (g>tn) -> ESTIMATE, pret stays t_i=0
    pret, action = plan_cell(3, 0, TLIST, 0, "varying", 1)
    assert action == ESTIMATE
    assert pret == 0


def test_plan_cell_post_period_uses_last_pretreatment():
    # g=3, t_i=2 (tn=4): post-period (g<=tn) -> pret = last_pretreatment = 1
    pret, action = plan_cell(3, 2, TLIST, 0, "varying", 1)
    assert action == ESTIMATE
    assert pret == 1


def test_plan_cell_universal_base_zero():
    # universal base: g=3, t_i=1 (tfac=0, tn=tlist[1]=2); pret=lpi(3)=1,
    # pret_year = tlist[1] = 2 == tn -> BASE_ZERO
    pret, action = plan_cell(3, 1, TLIST, 0, "universal", 0)
    assert action == BASE_ZERO


def test_plan_cell_universal_no_pretreatment_breaks():
    # g=1 has no pre-treatment period under universal base -> warn + BREAK
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pret, action = plan_cell(1, 0, TLIST, 0, "universal", 0)
    assert action == BREAK


def test_plan_cell_post_no_pretreatment_breaks():
    # g=1, post period, no pre period -> BREAK
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, action = plan_cell(1, 1, TLIST, 0, "varying", 1)
    assert action == BREAK


# --------------------------------------------------------------------------- #
# select_estimators
# --------------------------------------------------------------------------- #
def test_select_estimators_dr():
    p, r = select_estimators("dr")
    assert p is drdid_trim.drdid_panel
    assert r is drdid_trim.drdid_rc


def test_select_estimators_ipw():
    p, r = select_estimators("ipw")
    assert p is drdid_trim.std_ipw_did_panel
    assert r is drdid_trim.std_ipw_did_rc


def test_select_estimators_reg():
    p, r = select_estimators("reg")
    assert p is reg_did.reg_did_panel
    assert r is reg_did.reg_did_rc


def test_select_estimators_callable():
    f = lambda *a, **k: 0
    p, r = select_estimators(f)
    assert p is f and r is f


# --------------------------------------------------------------------------- #
# rcond_check_fail
# --------------------------------------------------------------------------- #
def test_rcond_full_rank_passes():
    cov = np.column_stack([np.ones(5), np.arange(5.0)])
    assert not rcond_check_fail(cov)


def test_rcond_collinear_fails():
    cov = np.column_stack([np.ones(5), np.ones(5)])
    assert rcond_check_fail(cov)


def test_rcond_no_control_fails():
    assert rcond_check_fail(np.zeros((0, 2)))


# --------------------------------------------------------------------------- #
# overlap_check_fail (intercept-only: pbar >= 0.999)  -- lines 97, 100, 101
# --------------------------------------------------------------------------- #
def test_overlap_balanced_ok():
    # pbar = 0.5 < 0.999 -> no violation
    assert overlap_check_fail(np.ones((10, 1)), np.array([0, 1] * 5)) is False


def test_overlap_all_treated_violates():
    # pbar = 1.0 >= 0.999 -> violation
    assert overlap_check_fail(np.ones((10, 1)), np.ones(10)) is True


def test_overlap_threshold_is_0999():
    """Just below 0.999 must pass; at/above must fail. Pins the 0.999 constant
    and the >= comparison (lines 100-101)."""
    n = 1000
    # 998 treated -> pbar = 0.998 < 0.999 -> ok
    y_below = np.array([1] * 998 + [0] * 2)
    assert overlap_check_fail(np.ones((n, 1)), y_below) is False
    # 999 treated -> pbar = 0.999 >= 0.999 -> fail
    y_at = np.array([1] * 999 + [0] * 1)
    assert overlap_check_fail(np.ones((n, 1)), y_at) is True


def test_overlap_empty_no_violation():
    # pbar defaults to 0.0 for empty -> no violation
    assert overlap_check_fail(np.zeros((0, 1)), np.array([])) is False


# --------------------------------------------------------------------------- #
# drdid_design_singular : design-set selection by method/path
# --------------------------------------------------------------------------- #
def _panel_designs():
    cov = np.column_stack([np.ones(6), np.array([1, 2, 3, 4, 5, 6.0])])
    G = np.array([0, 0, 0, 1, 1, 1.0])
    w = np.ones(6)
    return cov, w, G


def test_design_singular_panel_full_rank():
    cov, w, G = _panel_designs()
    assert drdid_design_singular(cov, w, G, "reg", post=None) is False


def test_design_singular_panel_control_singular():
    cov = np.column_stack([np.ones(6), np.ones(6)])
    G = np.array([0, 0, 0, 1, 1, 1.0])
    w = np.ones(6)
    assert drdid_design_singular(cov, w, G, "reg", post=None) is True


def test_design_singular_ipw_checks_ps_fullsample():
    cov = np.column_stack([np.ones(6), np.ones(6)])  # singular full-sample design
    G = np.array([0, 0, 0, 1, 1, 1.0])
    w = np.ones(6)
    # ipw doesn't fit outcome regressions but DOES fit the PS logit -> singular
    assert drdid_design_singular(cov, w, G, "ipw", post=None) is True


def test_design_singular_rc_full_rank():
    cov = np.column_stack([np.ones(8), np.arange(8.0)])
    G = np.array([0, 0, 0, 0, 1, 1, 1, 1.0])
    post = np.array([0, 1, 0, 1, 0, 1, 0, 1.0])
    w = np.ones(8)
    assert drdid_design_singular(cov, w, G, "reg", post=post) is False


def test_design_singular_rc_tiny_treated_post_fails():
    """RC dr/reg checks all 4 group*period designs; a treated-post cell with a
    single obs cannot support a 2-column design -> singular. Pins the RC mask set
    (lines 174-176)."""
    cov = np.column_stack([np.ones(8), np.arange(8.0)])
    G = np.array([0, 0, 0, 0, 1, 1, 1, 1.0])
    # treated rows 4-7 have post = [0,0,0,1] -> only ONE treated-post obs
    post = np.array([0, 1, 0, 1, 0, 0, 0, 1.0])
    w = np.ones(8)
    assert drdid_design_singular(cov, w, G, "dr", post=post) is True
    # ipw on RC only checks the full-sample PS design (full rank here) -> OK
    assert drdid_design_singular(cov, w, G, "ipw", post=post) is False
