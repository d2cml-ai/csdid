"""Shared helpers for the standard (``compute_att_gt``) and faster_mode
(``compute_att_gt2``) group-time ATT paths.

Centralizing the estimator dispatch and the per-(g, t) pretreatment-period /
control-flow logic keeps the two paths in lock-step. The only thing that
legitimately differs between them is *how* the 2x2 cell data is sourced
(per-cell subsetting + ``panel2cs2`` in the standard path vs. precomputed
per-period lookups in the faster path); the cell *selection* logic below must
stay identical, so it lives here once.
"""
import warnings

import numpy as np

from drdid import reg_did
from csdid.attgt_fnc import drdid_trim


def select_estimators(est_method):
    """Return ``(panel_estimator, rc_estimator)`` for the requested method.

    A callable ``est_method`` is used for both panel and repeated-cross-section
    cells (matches R `did`, which accepts a custom estimator).
    """
    if callable(est_method):
        return est_method, est_method
    panel = {
        "reg": reg_did.reg_did_panel,
        "ipw": drdid_trim.std_ipw_did_panel,
        "dr": drdid_trim.drdid_panel,
    }[est_method]
    rc = {
        "reg": reg_did.reg_did_rc,
        "ipw": drdid_trim.std_ipw_did_rc,
        "dr": drdid_trim.drdid_rc,
    }[est_method]
    return panel, rc


def rcond_check_fail(covariates_control):
    """Per-(g, t) regression-feasibility guard, mirroring R `did`'s
    ``rcond_check_fail`` (DRDID >= 1.3.0 era): return ``True`` when the
    control-unit design's Gram matrix ``Xc'Xc`` is singular or numerically
    ill-conditioned, i.e. ``rcond(crossprod(control_covs)) < .Machine$double.eps``.

    ``covariates_control`` is the design matrix (intercept column included, as the
    port already builds it) restricted to the comparison's control units (treatment
    indicator ``== 0``). R uses the LAPACK 1-norm reciprocal condition number, so we
    match with ``1 / np.linalg.cond(M, 1)`` and the ``np.finfo(float).eps`` cutoff
    (= R's ``.Machine$double.eps``). A singular Gram (``inv`` fails / non-finite
    cond) gives reciprocal condition number 0 -> fail, as in R; with no control
    units the Gram is empty and the check fails, mirroring R's intercept-only
    ``!any(y == 0)`` closed form. Applied only for outcome-regression methods
    (``dr``/``reg``); ``ipw`` and custom estimators are exempt, matching R.
    """
    xc = np.asarray(covariates_control, dtype=float)
    if xc.ndim == 1:
        xc = xc[:, None]
    if xc.shape[0] == 0:
        return True  # no control units
    gram = xc.T @ xc
    try:
        rcond = 1.0 / np.linalg.cond(gram, 1)
    except np.linalg.LinAlgError:
        return True  # singular Gram -> reciprocal condition number 0
    if not np.isfinite(rcond):
        return True
    return rcond < np.finfo(float).eps


def overlap_check_fail(int_cov, treated, intercept_only=None):
    """Per-(g, t) propensity-overlap guard, mirroring R `did`'s
    ``overlap_check_fail`` (DRDID >= 1.3.0 era): NA a ``dr``/``ipw`` cell whose
    propensity-score model is (near-)separated, rather than letting exploding IPW
    weights return a garbage estimate.

    R fits the treatment logit ``treated ~ int_cov`` on the 2x2 cell sample
    (``fastglm::fastglmPure``, IRLS) and fails when ``max(fitted) >= 0.999``. We
    use ``statsmodels`` GLM (Binomial / IRLS) -- the same fit the audit harness
    uses to MEASURE the overlap mechanism -- so the port's NA decision tracks R's.

    ``int_cov`` is the design matrix (intercept column included, as the port
    builds it); ``treated`` is the treatment indicator (``G``/``D``). For an
    intercept-only model R short-circuits to the treated fraction
    (``mean(treated)``); we detect intercept-only as a single-column design (the
    port always carries an intercept), matching R's ``xformla == ~1``. Applied
    only to ``dr``/``ipw`` (the PS-fitting methods); ``reg`` and custom estimators
    are exempt, matching R.
    """
    y = np.asarray(treated, dtype=float).ravel()
    x = np.asarray(int_cov, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    if intercept_only is None:
        intercept_only = (x.shape[1] == 1)
    if intercept_only:
        pbar = float(np.mean(y)) if y.size else 0.0
        # R: only short-circuits when not right at the 0.999 boundary; for an
        # intercept-only logit every fitted value equals pbar anyway.
        if abs(pbar - 0.999) > 1e-6:
            return pbar >= 0.999
        return pbar >= 0.999
    import statsmodels.api as sm
    import warnings as _w
    try:
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            res = sm.GLM(y, x, family=sm.families.Binomial()).fit()
            fitted = np.asarray(res.fittedvalues, dtype=float)
    except sm.tools.sm_exceptions.PerfectSeparationError:
        # Perfect separation -> fitted PS reaches 1 -> overlap violated.
        return True
    except Exception:
        # A non-separation fit failure is not an overlap violation; let the
        # downstream estimator (or other guards) handle it.
        return False
    if fitted.size == 0 or not np.isfinite(fitted).any():
        return False
    return float(np.nanmax(fitted)) >= 0.999


def _gram_singular(int_cov, i_weights, mask):
    """``rcond(crossprod(w*mask * int_cov, int_cov)) < eps`` -- the weighted Gram of
    the design restricted to ``mask`` (a units-subset boolean). Mirrors one of
    DRDID's internal ``rcond`` design checks. Returns ``True`` if singular / empty.
    The 1/n scaling R applies is omitted (scale-invariant for ``rcond``)."""
    m = np.asarray(mask, dtype=bool)
    if not m.any():
        return True
    x = np.asarray(int_cov, dtype=float)[m]
    w = np.asarray(i_weights, dtype=float)[m]
    gram = (x * w[:, None]).T @ x
    try:
        rcond = 1.0 / np.linalg.cond(gram, 1)
    except np.linalg.LinAlgError:
        return True
    if not np.isfinite(rcond):
        return True
    return rcond < np.finfo(float).eps


def drdid_design_singular(int_cov, i_weights, treated, est_method, post=None):
    """Replicate DRDID's *internal* per-cell design-feasibility checks
    (DRDID >= 1.3.0), which NA a cell whose 2x2 estimating designs are singular:

      * ``dr``/``reg`` fit outcome regressions for BOTH the control and treated
        groups, so each group's weighted design must be non-singular. On the
        repeated-cross-section / unbalanced path (``post`` given) the control and
        treated outcome regressions are fit separately pre and post, so all four
        period-by-group designs are checked (this is what catches a tiny treated
        cohort with too few post observations to support the covariate design --
        e.g. 3 treated obs against a 4-column design). On the balanced panel path
        (``post=None``) the differenced design is checked once per group.
      * ``dr``/``ipw`` additionally fit a propensity-score logit on the full
        sample, so the full-sample design is checked too.

    ``treated`` is the treatment indicator (``D``/``G_m``). Returns ``True`` if any
    applicable design is singular, matching R's NA-the-cell behavior. Custom
    (callable) estimators are exempt and handled by the caller.
    """
    treated = np.asarray(treated).astype(bool)
    control = ~treated
    masks = []
    if est_method in ("dr", "reg"):
        if post is None:
            # Balanced panel (drdid_panel/reg_did_panel): ONLY the control-group
            # outcome design is checked -- there is NO treated-design check.
            masks += [control]
        else:
            # Repeated cross-section / unbalanced (drdid_rc/reg_did_rc): control
            # AND treated outcome regressions are fit separately pre and post, so
            # all four period-by-group designs are checked. The treated-design
            # checks are what NA a tiny treated cohort with too few observations.
            post = np.asarray(post).astype(bool)
            pre = ~post
            masks += [control & post, control & pre, treated & post, treated & pre]
    for m in masks:
        if _gram_singular(int_cov, i_weights, m):
            return True
    if est_method in ("dr", "ipw"):
        # Propensity-score design (full sample), checked by dr/ipw on both paths.
        if _gram_singular(int_cov, i_weights, np.ones(len(treated), dtype=bool)):
            return True
    return False


def last_pretreatment_index(g, tlist, anticipation):
    """Index of the last period that is pre-treatment for cohort ``g``.

    Mirrors R's ``tail(which((tlist + anticipation) < g), 1)``. Returns ``None``
    when the cohort has no pre-treatment period.
    """
    idx = np.where((np.asarray(tlist) + anticipation) < g)[0]
    return int(idx[-1]) if len(idx) else None


# Actions returned by ``plan_cell`` describing what the caller should do.
ESTIMATE = "estimate"   # estimate ATT(g, t) for this cell
BASE_ZERO = "base_zero"  # universal base cell: append att=0 with a zero IF
BREAK = "break"          # no pre-treatment period: stop scanning periods for g


def plan_cell(g, t_i, tlist, anticipation, base_period, tfac):
    """Resolve the pretreatment period and control flow for one (g, t) cell.

    Shared verbatim by the standard and faster_mode paths so they cannot drift.
    Returns ``(pret, action)`` where ``pret`` is the (final) pretreatment period
    index and ``action`` is one of :data:`ESTIMATE`, :data:`BASE_ZERO`,
    :data:`BREAK`. Raises ``ValueError`` for a universal base period with no
    pre-treatment period (matches the standard path).

    The caller is expected to derive ``tn``/``pret_year``/``post_treat`` from its
    own ``tlist`` so that stored values keep their original dtypes.
    """
    tl = np.asarray(tlist)
    pret = t_i
    tn = tl[t_i + tfac]

    # Universal base period: fixed reference = last pre-treatment period.
    if base_period == "universal":
        idx = last_pretreatment_index(g, tl, anticipation)
        if idx is None:
            raise ValueError(
                f"There are no pre-treatment periods for the group first treated at {g}. "
                f"Units from this group are dropped."
            )
        pret = idx

    # Post-treatment period: reference is the last pre-treatment period.
    if g <= tn:
        idx = last_pretreatment_index(g, tl, anticipation)
        if idx is None:
            warnings.warn(
                f"There are no pre-treatment periods for the group first treated at {g}\n"
                f"Units from this group are dropped"
            )
            return pret, BREAK
        pret = idx

    pret_year = tl[pret]
    if base_period == "universal" and pret_year == tn:
        return pret, BASE_ZERO
    return pret, ESTIMATE
