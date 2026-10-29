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
from csdid.utils.weights import _weight_mean, _normalize_weights


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


def rcond_check_fail(covariates_control, treated=None, i_weights=None):
    """Check the actual normalized-weight control Gram, as in R ``did``.

    The one-argument form accepts an already restricted control design. Passing
    ``treated`` allows normalization on the complete cell before restricting its
    controls, which matters when a positive weight rounds to zero locally.
    """
    x = np.asarray(covariates_control, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    control = (np.ones(len(x), dtype=bool) if treated is None
               else np.asarray(treated).ravel() == 0)
    if i_weights is None:
        xc = x[control]
        if len(xc) == 0:
            return True
        if x.shape[1] == 1 and np.all(x[:, 0] == 1):
            return False
        gram = xc.T @ xc
    else:
        w = np.asarray(i_weights, dtype=float).ravel()
        with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
            w = _normalize_weights(w)
        control &= np.isfinite(w) & (w > 0)
        if not control.any():
            return True
        if x.shape[1] == 1 and np.all(x[:, 0] == 1):
            return False
        xc = x[control]
        gram = (xc * w[control, None]).T @ xc / len(x)
    try:
        rcond = 1.0 / np.linalg.cond(gram, 1)
    except np.linalg.LinAlgError:
        return True
    return bool(not np.isfinite(rcond) or rcond < np.finfo(float).eps)


def overlap_check_fail(int_cov, treated, intercept_only=None, i_weights=None):
    """Per-(g, t) propensity-overlap guard, mirroring R `did`'s
    ``overlap_check_fail`` (DRDID >= 1.3.0 era): NA a ``dr``/``ipw`` cell whose
    propensity-score model is (near-)separated, rather than letting exploding IPW
    weights return a garbage estimate.

    R fits the treatment logit ``treated ~ int_cov`` on the 2x2 cell sample
    (``fastglm::fastglmPure``, IRLS) on positive normalized-weight support and
    fails when ``max(fitted) >= 0.999``. The preliminary fit stays unweighted. We
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
    if i_weights is not None:
        w = np.asarray(i_weights, dtype=float).ravel()
        with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
            normalized = _normalize_weights(w)
        supported = np.isfinite(normalized) & (normalized > 0)
        if not supported.all():
            x, y = x[supported], y[supported]
    if y.size == 0:
        return True
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
    positive = np.isfinite(w) & (w > 0)
    if not positive.any():
        return True
    x, w = x[positive], w[positive]
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

      * Panel ``dr``/``reg`` fit an outcome regression only for controls.
        Repeated-cross-section ``reg`` fits controls separately pre and post;
        ``dr`` additionally fits treated outcome regressions in both periods.
        Check only designs actually fitted by the selected estimator, so a tiny
        treated cell does not incorrectly invalidate a regression-only estimate.
      * ``dr``/``ipw`` additionally fit a propensity-score logit on the full
        sample, so the full-sample design is checked too.

    ``treated`` is the treatment indicator (``D``/``G_m``). Returns ``True`` if any
    applicable design is singular, matching R's NA-the-cell behavior. Custom
    (callable) estimators are exempt and handled by the caller.
    """
    treated = np.asarray(treated).astype(bool)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
        i_weights = _normalize_weights(np.asarray(i_weights, dtype=float))
    control = ~treated
    masks = []
    if est_method in ("dr", "reg"):
        if post is None:
            # Balanced panel (drdid_panel/reg_did_panel): ONLY the control-group
            # outcome design is checked -- there is NO treated-design check.
            masks += [control]
        else:
            # Only DR fits treated regressions in addition to control regressions.
            post = np.asarray(post).astype(bool)
            pre = ~post
            masks += [control & post, control & pre]
            if est_method == "dr":
                masks += [treated & post, treated & pre]
    for m in masks:
        if _gram_singular(int_cov, i_weights, m):
            return True
    if est_method in ("dr", "ipw"):
        # Propensity-score design (full sample), checked by dr/ipw on both paths.
        if _gram_singular(int_cov, i_weights, np.ones(len(treated), dtype=bool)):
            return True
    return False


def rc_outcome_has_too_few_obs(treated, post, i_weights, n_covariates,
                               include_treated=True):
    """Necessary physical-size check before each fitted RC outcome design.

    Count finite positive weights after the estimator's local normalization;
    even a positive subnormal counts, but a weight rounded to zero does not.
    Passing this check does not certify full rank.
    """
    w = np.asarray(i_weights, dtype=float).ravel()
    with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
        scale = _weight_mean(w) if w.size else np.nan
        if not np.isfinite(scale) or scale <= 0:
            return True
        normalized = w / scale
    positive = np.isfinite(normalized) & (normalized > 0)
    labels = np.asarray(post, dtype=int) + 2 * np.asarray(treated, dtype=int)
    counts = np.bincount(labels[positive], minlength=4)
    required = counts if include_treated else counts[:2]
    return bool(np.any(required < n_covariates))


_ORIGINAL_PANEL_GUARDS = (overlap_check_fail, rcond_check_fail, drdid_design_singular)


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
    :data:`BREAK`. For a universal base period with no pre-treatment period it
    warns and returns :data:`BREAK` (warn-and-drop), mirroring R `did` and the
    varying-base path (v7-DIV2).

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
            # v7-DIV2: match R `did` -- warn and drop the cohort (warn-and-drop),
            # mirroring the varying-base path below, rather than raising.
            warnings.warn(
                f"There are no pre-treatment periods for the group first treated at {g}\n"
                f"Units from this group are dropped"
            )
            return pret, BREAK
        pret = idx

    # Treatment and anticipation periods share the last unaffected reference.
    # Within the anticipation window, adjacent-period differences would remove
    # an anticipation effect that persists over more than one period.
    if g - anticipation <= tn:
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
