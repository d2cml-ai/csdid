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
