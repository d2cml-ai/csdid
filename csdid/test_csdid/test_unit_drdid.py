"""Bespoke, self-contained, mutation-gated unit tests for the DR/IPW estimator
math in ``csdid/attgt_fnc/drdid_trim.py``.

NO R at runtime. Expected values are ANALYTIC / closed-form for the intercept-only
(no-covariate) case, where the doubly-robust and IPW estimators collapse to a
simple 2x2 (or 2-group) difference of weighted means, and the influence function
has a known closed form. Any sign/operator flip on the ATT or IF arithmetic, or a
change to the trim threshold / weight normalization / NA guards, breaks one of
these identities -> kills the mutant.

Key closed forms (intercept-only, equal/given weights):
  PANEL (drdid_panel, std_ipw_did_panel):
     att = mean(dy | treated) - mean(dy | control)
     IF_i = D_i/p * (dy_i - eta_t) - (1-D_i)/(1-p) * (dy_i - eta_c)
            with p = mean(D), eta_t = mean(dy|treat), eta_c = mean(dy|control)
  RC (std_ipw_did_rc, drdid_rc):
     att = (mu11 - mu10) - (mu01 - mu00)         (mu_dt = mean(y|D=d,post=t))
     IF_i = I(1,1)/lam11*(y-mu11) - I(1,0)/lam10*(y-mu10)
            - ( I(0,1)/lam01*(y-mu01) - I(0,0)/lam00*(y-mu00) )
The DR and IPW estimators are IDENTICAL in the intercept-only case (the OR term
drops out), giving an extra cross-check oracle.
"""
import numpy as np
import pytest

from csdid.attgt_fnc.drdid_trim import (
    _trim_pscore,
    drdid_panel,
    std_ipw_did_panel,
    std_ipw_did_rc,
    drdid_rc,
)

TOL = 1e-9


# --------------------------------------------------------------------------- #
# helpers: analytic closed forms
# --------------------------------------------------------------------------- #
def _panel_closed_form(dy, D, w=None):
    D = np.asarray(D, float)
    dy = np.asarray(dy, float)
    if w is None:
        w = np.ones_like(D)
    w = w / w.mean()
    wt = w[D == 1]
    wc = w[D == 0]
    eta_t = np.sum(wt * dy[D == 1]) / np.sum(wt)
    eta_c = np.sum(wc * dy[D == 0]) / np.sum(wc)
    att = eta_t - eta_c
    return att, eta_t, eta_c


def _panel_if_closed_form(dy, D):
    """IF for the equal-weight intercept-only case."""
    D = np.asarray(D, float)
    dy = np.asarray(dy, float)
    p = D.mean()
    eta_t = dy[D == 1].mean()
    eta_c = dy[D == 0].mean()
    return D / p * (dy - eta_t) - (1 - D) / (1 - p) * (dy - eta_c)


def _rc_means(y, D, post):
    def mu(dd, tt):
        m = (D == dd) & (post == tt)
        return y[m].mean()
    return mu(1, 1), mu(1, 0), mu(0, 1), mu(0, 0)


def _rc_if_closed_form(y, D, post):
    y = np.asarray(y, float)
    D = np.asarray(D)
    post = np.asarray(post)

    def lam(dd, tt):
        return np.mean((D == dd) & (post == tt))

    def ind(dd, tt):
        return ((D == dd) & (post == tt)).astype(float)

    mu11, mu10, mu01, mu00 = _rc_means(y, D, post)
    return (ind(1, 1) / lam(1, 1) * (y - mu11)
            - ind(1, 0) / lam(1, 0) * (y - mu10)
            - (ind(0, 1) / lam(0, 1) * (y - mu01)
               - ind(0, 0) / lam(0, 0) * (y - mu00)))


def _panel_data(seed=11, n=400, att=3.0):
    rng = np.random.default_rng(seed)
    D = np.array([0, 1] * (n // 2), float)
    dy = rng.normal(0.0, 1.0, n) + att * D
    return dy, D


def _rc_data(seed=12, n=1200, att=2.0):
    rng = np.random.default_rng(seed)
    D = np.array([0, 1] * (n // 2), float)
    post = rng.integers(0, 2, n).astype(float)
    y = rng.normal(0.0, 1.0, n) + att * (D * post)
    return y, D, post


# --------------------------------------------------------------------------- #
# PANEL: ATT closed form  (kills att-arithmetic mutants: lines 55-57, 114-116)
# --------------------------------------------------------------------------- #
def test_drdid_panel_att_closed_form():
    dy, D = _panel_data()
    att, ifn = drdid_panel(dy, np.zeros_like(dy), D, covariates=None)
    exp, _, _ = _panel_closed_form(dy, D)
    assert abs(att - exp) < TOL
    assert len(ifn) == len(D)


def test_ipw_panel_att_closed_form():
    dy, D = _panel_data(seed=7)
    att, _ = std_ipw_did_panel(dy, np.zeros_like(dy), D, covariates=None)
    exp, _, _ = _panel_closed_form(dy, D)
    assert abs(att - exp) < TOL


def test_panel_dr_equals_ipw_intercept_only():
    """DR == IPW when there are no covariates (OR term vanishes)."""
    dy, D = _panel_data(seed=99)
    a_dr, if_dr = drdid_panel(dy, np.zeros_like(dy), D, covariates=None)
    a_ipw, if_ipw = std_ipw_did_panel(dy, np.zeros_like(dy), D, covariates=None)
    assert abs(a_dr - a_ipw) < TOL
    assert np.max(np.abs(if_dr - if_ipw)) < 1e-8


# --------------------------------------------------------------------------- #
# PANEL: full IF closed form (kills lines 59-70, 109-141 arithmetic/signs)
# --------------------------------------------------------------------------- #
def test_drdid_panel_if_closed_form():
    dy, D = _panel_data(seed=1, n=600)
    _, ifn = drdid_panel(dy, np.zeros_like(dy), D, covariates=None)
    exp = _panel_if_closed_form(dy, D)
    assert np.max(np.abs(ifn - exp)) < 1e-8


def test_ipw_panel_if_closed_form():
    dy, D = _panel_data(seed=2, n=600)
    _, ifn = std_ipw_did_panel(dy, np.zeros_like(dy), D, covariates=None)
    exp = _panel_if_closed_form(dy, D)
    assert np.max(np.abs(ifn - exp)) < 1e-8


def test_panel_if_mean_zero():
    """IF must be mean-zero (a defining property of an influence function).
    Any single sign flip in the IF assembly breaks the mean-zero property."""
    dy, D = _panel_data(seed=3)
    for fn in (drdid_panel, std_ipw_did_panel):
        _, ifn = fn(dy, np.zeros_like(dy), D, covariates=None)
        assert abs(ifn.mean()) < 1e-8


# --------------------------------------------------------------------------- #
# PANEL: weight normalization + non-negativity guard (lines 31-37, 80-86)
# --------------------------------------------------------------------------- #
def test_panel_weight_scale_invariance():
    """i_weights/mean(i_weights) -> ATT invariant to a positive rescale.
    Kills the `Div`->`Mult` mutant on the normalization line."""
    dy, D = _panel_data(seed=5, n=300)
    rng = np.random.default_rng(0)
    w = rng.uniform(0.5, 2.0, len(D))
    a1, _ = drdid_panel(dy, np.zeros_like(dy), D, i_weights=w)
    a2, _ = drdid_panel(dy, np.zeros_like(dy), D, i_weights=10.0 * w)
    assert abs(a1 - a2) < 1e-9


def test_panel_weighted_att_matches_closed_form():
    dy, D = _panel_data(seed=6, n=300)
    rng = np.random.default_rng(1)
    w = rng.uniform(0.5, 2.0, len(D))
    a, _ = drdid_panel(dy, np.zeros_like(dy), D, i_weights=w)
    exp, _, _ = _panel_closed_form(dy, D, w)
    assert abs(a - exp) < 1e-9


def test_panel_negative_weight_raises():
    dy, D = _panel_data(seed=8, n=100)
    bad = np.ones(len(D))
    bad[0] = -1.0
    with pytest.raises(ValueError):
        drdid_panel(dy, np.zeros_like(dy), D, i_weights=bad)
    with pytest.raises(ValueError):
        std_ipw_did_panel(dy, np.zeros_like(dy), D, i_weights=bad)


# --------------------------------------------------------------------------- #
# RC: ATT closed form  (kills att-arithmetic: lines 175-203 / 269-289)
# --------------------------------------------------------------------------- #
def test_ipw_rc_att_closed_form():
    y, D, post = _rc_data(seed=2)
    att, _ = std_ipw_did_rc(y, post, D, covariates=None)
    mu11, mu10, mu01, mu00 = _rc_means(y, D, post)
    exp = (mu11 - mu10) - (mu01 - mu00)
    assert abs(att - exp) < TOL


def test_drdid_rc_att_closed_form():
    y, D, post = _rc_data(seed=3)
    att, _ = drdid_rc(y, post, D, covariates=None)
    mu11, mu10, mu01, mu00 = _rc_means(y, D, post)
    exp = (mu11 - mu10) - (mu01 - mu00)
    assert abs(att - exp) < TOL


def test_rc_dr_equals_ipw_intercept_only():
    y, D, post = _rc_data(seed=4)
    a_dr, if_dr = drdid_rc(y, post, D, covariates=None)
    a_ipw, if_ipw = std_ipw_did_rc(y, post, D, covariates=None)
    assert abs(a_dr - a_ipw) < 1e-8
    assert np.max(np.abs(if_dr - if_ipw)) < 1e-7


# --------------------------------------------------------------------------- #
# RC: full IF closed form (kills lines 186-203 / 315-358 arithmetic/signs)
# --------------------------------------------------------------------------- #
def test_ipw_rc_if_closed_form():
    y, D, post = _rc_data(seed=12, n=1600)
    _, ifn = std_ipw_did_rc(y, post, D, covariates=None)
    exp = _rc_if_closed_form(y, D, post)
    assert np.max(np.abs(ifn - exp)) < 1e-7


def test_drdid_rc_if_closed_form():
    y, D, post = _rc_data(seed=13, n=1600)
    _, ifn = drdid_rc(y, post, D, covariates=None)
    exp = _rc_if_closed_form(y, D, post)
    assert np.max(np.abs(ifn - exp)) < 1e-7


def test_rc_if_mean_zero():
    y, D, post = _rc_data(seed=14)
    for fn in (std_ipw_did_rc, drdid_rc):
        _, ifn = fn(y, post, D, covariates=None)
        assert abs(ifn.mean()) < 1e-7


def test_rc_negative_weight_raises():
    y, D, post = _rc_data(seed=15, n=200)
    bad = np.ones(len(D))
    bad[0] = -1.0
    with pytest.raises(ValueError):
        std_ipw_did_rc(y, post, D, i_weights=bad)


# --------------------------------------------------------------------------- #
# _trim_pscore: threshold + clamp + control selection (lines 17-22, 151-157)
# --------------------------------------------------------------------------- #
def test_trim_pscore_clamp_upper():
    """ps is clamped to 1 - 1e-6. Kills the `1 - 1e-6` const mutants."""
    ps = np.array([0.5, 1.0, 2.0, 0.999999999])
    out, _ = _trim_pscore(ps.copy(), np.zeros(4), 0.995)
    assert np.all(out <= 1 - 1e-6 + 1e-15)
    assert abs(out.max() - (1 - 1e-6)) < 1e-12


def test_trim_pscore_controls_threshold():
    """Control units with ps >= trim_level are trimmed (False);
    those below are kept (True). Kills the `< trim_level` compare + `D==0`."""
    ps = np.array([0.10, 0.50, 0.99, 0.996, 0.9999])
    D = np.zeros(5)  # all controls
    _, trim = _trim_pscore(ps.copy(), D, 0.995)
    # < 0.995 -> True (kept); >= 0.995 -> False (trimmed)
    assert list(trim) == [True, True, True, False, False]


def test_trim_pscore_treated_never_trimmed():
    """Treated units use the ps < 1.01 guard (always True after clamp)."""
    ps = np.array([0.10, 0.50, 0.99, 0.9999])
    D = np.ones(4)  # all treated
    _, trim = _trim_pscore(ps.copy(), D, 0.995)
    assert np.all(trim)


def test_trim_pscore_mixed_only_controls_affected():
    ps = np.array([0.9999, 0.9999])
    D = np.array([1.0, 0.0])  # treated kept, control trimmed
    _, trim = _trim_pscore(ps.copy(), D, 0.995)
    assert trim[0] == True   # treated kept
    assert trim[1] == False  # control trimmed


# --------------------------------------------------------------------------- #
# WITH-COVARIATES: finite-difference Gateaux-derivative IF oracle.
#
# The influence function is the Gateaux derivative of the estimator wrt the
# empirical distribution: IF_i = d ATT / d (eps * delta_i), recoverable by
# upweighting observation i by eps and central-differencing the ATT. This is a
# fully self-contained oracle (NO R) that exercises the FULL propensity-score
# correction machinery (score_ps, Hessian_ps, asy_lin_rep_*, M1/M2/M3 moment
# terms) which VANISHES in the intercept-only case. Any sign/operator flip on the
# PS-correction or IF-assembly lines breaks the FD<->IF match.
# --------------------------------------------------------------------------- #
def _fd_if_max_error(fn, y, D, cov, post=None, eps=1e-5, idxs=(0, 5, 17, 99)):
    n = len(D)

    def att_with_weights(w):
        if post is None:
            a, _ = fn(y, np.zeros(n), D, covariates=cov, i_weights=w)
        else:
            a, _ = fn(y, post, D, covariates=cov, i_weights=w)
        return a

    if post is None:
        _, ifn = fn(y, np.zeros(n), D, covariates=cov)
    else:
        _, ifn = fn(y, post, D, covariates=cov)

    worst = 0.0
    for i in idxs:
        wp = np.ones(n); wp[i] += eps * n
        wm = np.ones(n); wm[i] -= eps * n
        fd = (att_with_weights(wp) - att_with_weights(wm)) / (2 * eps)
        worst = max(worst, abs(fd - ifn[i]))
    return worst


def _cov_panel_data(seed=7, n=400):
    rng = np.random.default_rng(seed)
    X = rng.normal(0, 1, n)
    ps = 1.0 / (1.0 + np.exp(-0.5 * X))
    D = (rng.uniform(size=n) < ps).astype(float)
    dy = 2.0 * D + 0.8 * X + rng.normal(0, 0.5, n)
    return dy, D, X[:, None]


def test_drdid_panel_if_fd_with_covariates():
    """The DR panel estimator is Neyman-orthogonal, so its influence function
    equals the (central) finite-difference Gateaux derivative EXACTLY -- including
    the propensity-score correction (score_ps, Hessian_ps, asy_lin_rep_wols/ps,
    M1/M2/M3) which vanishes in the intercept-only case. Kills the sign/operator
    mutants on the WITH-covariate DR PS-correction and IF-assembly lines.

    (Only DR is asserted: the IPW estimator is NOT orthogonal, so its finite
    finite-sample IF is not the simple Gateaux derivative of the plugin estimator;
    those IPW PS-correction mutants are documented residuals in the Phase-2
    report -- they hinge on statsmodels GLM covariance internals that have no
    self-contained analytic oracle without R.)
    """
    dy, D, cov = _cov_panel_data(seed=7)
    assert _fd_if_max_error(drdid_panel, dy, D, cov) < 1e-3


# --------------------------------------------------------------------------- #
# trim_level default (0.995): plant a few control units at very high covariate
# values so their fitted propensity score exceeds 0.995. Then the DEFAULT
# trim_level trims them, so the ATT DIFFERS from the no-trim (trim_level=2.0)
# ATT. Pins the 0.995 default-arg constant on each estimator.
# --------------------------------------------------------------------------- #
def _planted_panel(seed=0, n=400):
    rng = np.random.default_rng(seed)
    X = rng.normal(0, 1, n)
    D = (X > 0).astype(float)
    ctrl = np.where(D == 0)[0][:5]
    X[ctrl] = 4.0  # high-X controls -> fitted ps >= 0.995 -> trimmed by default
    y = 2.0 * D + X + rng.normal(0, 0.3, n)
    return y, D, X[:, None]


def test_panel_trim_default_is_active():
    y, D, cov = _planted_panel(seed=0)
    a_default, _ = std_ipw_did_panel(y, np.zeros_like(y), D, covariates=cov)
    a_notrim, _ = std_ipw_did_panel(y, np.zeros_like(y), D, covariates=cov,
                                    trim_level=2.0)
    assert abs(a_default - a_notrim) > 1e-4
    a_default_dr, _ = drdid_panel(y, np.zeros_like(y), D, covariates=cov)
    a_notrim_dr, _ = drdid_panel(y, np.zeros_like(y), D, covariates=cov,
                                 trim_level=2.0)
    assert abs(a_default_dr - a_notrim_dr) > 1e-4


def test_rc_trim_default_is_active():
    rng = np.random.default_rng(0)
    n = 400
    X = rng.normal(0, 1, n)
    D = (X > 0).astype(float)
    ctrl = np.where(D == 0)[0][:5]
    X[ctrl] = 4.0
    post = rng.integers(0, 2, n).astype(float)
    y = 2.0 * (D * post) + X + rng.normal(0, 0.3, n)
    cov = X[:, None]
    a_default, _ = std_ipw_did_rc(y, post, D, covariates=cov)
    a_notrim, _ = std_ipw_did_rc(y, post, D, covariates=cov, trim_level=2.0)
    assert abs(a_default - a_notrim) > 1e-4
    a_default_dr, _ = drdid_rc(y, post, D, covariates=cov)
    a_notrim_dr, _ = drdid_rc(y, post, D, covariates=cov, trim_level=2.0)
    assert abs(a_default_dr - a_notrim_dr) > 1e-4


# --------------------------------------------------------------------------- #
# default i_weights=None path must SUCCEED on every estimator. Kills the
# `i_weights is None` -> `is not None` mutant (which would call None.flatten()).
# --------------------------------------------------------------------------- #
def test_default_weights_path_succeeds_all_estimators():
    dy, D = _panel_data(seed=21, n=200)
    for fn in (drdid_panel, std_ipw_did_panel):
        att, ifn = fn(dy, np.zeros_like(dy), D, covariates=None)
        assert np.isfinite(att) and np.all(np.isfinite(ifn))
    y, D2, post = _rc_data(seed=22, n=400)
    for fn in (std_ipw_did_rc, drdid_rc):
        att, ifn = fn(y, post, D2, covariates=None)
        assert np.isfinite(att) and np.all(np.isfinite(ifn))
