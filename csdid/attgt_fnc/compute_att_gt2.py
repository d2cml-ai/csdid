"""Vectorized (faster_mode) group-time ATT computation.

This is the Python equivalent of R `did`'s ``faster_mode=TRUE`` path. It produces
results numerically identical to ``compute_att_gt`` (the standard path) but avoids
the two dominant per-cell costs profiled in the standard loop:

* per-(g,t) ``patsy`` covariate construction — the design matrix is built **once
  per period** (panel) or **once globally** (repeated cross sections); and
* per-(g,t) ``panel2cs2`` reshaping — outcomes/weights are pre-pivoted to a
  per-period lookup (panel).

The 2x2 estimators, the pre-period/control-group logic, the fix_weights handling
and the influence-function bookkeeping are identical to ``compute_att_gt``.
"""
import numpy as np, pandas as pd
from drdid import reg_did
from csdid.attgt_fnc import drdid_trim
from csdid.attgt_fnc.compute_att_gt_shared import (
    select_estimators, last_pretreatment_index, plan_cell, rcond_check_fail,
    drdid_design_singular, overlap_check_fail, rc_outcome_has_too_few_obs, BREAK, BASE_ZERO,
    _ORIGINAL_PANEL_GUARDS)
import warnings


_CACHEABLE_PANEL_ESTIMATORS = drdid_trim._ORIGINAL_PANEL_ESTIMATORS


def _period_cache_keys(period_arrays):
    """Share the first period's key only for exactly unchanged arrays.

    Period identities already permit post-treatment reuse with time-varying
    data. Comparing each period to the first also covers time-invariant inputs
    in linear time, without hashes or tolerance-based equality.
    """
    first = period_arrays[0]
    return {p: (0 if values.dtype == first.dtype and np.array_equal(values, first)
                else p) for p, values in period_arrays.items()}


def compute_att_gt2(dp, est_method="dr", base_period="varying", compute_inffunc=True):
    yname = dp['yname']
    tname = dp['tname']
    idname = dp['idname']
    xformla = dp['xformla']
    data = dp['data']
    weights_name = dp['weights_name']
    panel = dp['panel']
    control_group = dp['control_group']
    anticipation = dp['anticipation']
    gname = dp['gname']
    n = dp['n']
    tlist = np.asarray(dp['tlist'])
    glist = dp['glist']
    fix_weights = dp.get('fix_weights')
    xcov_cols = dp.get('xcov_cols')

    tlist_len = len(tlist) - 1 if base_period != "universal" else len(tlist)
    tfac = 1 if base_period != "universal" else 0

    never_treated = control_group == 'nevertreated'

    inf_func = []
    att_est, group, year, post_array = [], [], [], []

    # current (g, tn) captured by add_att_data via closure variables
    state = {'g': None, 'tn': None}

    def add_att_data(att=0, pst=0, inf_f=None):
        if compute_inffunc:
            inf_func.append(inf_f)
        att_est.append(att)
        group.append(state['g'])
        year.append(state['tn'])
        post_array.append(pst)

    est_panel, est_rc = select_estimators(est_method)

    # Regression-feasibility guards (match R): `did`'s control pre-check (dr/reg)
    # plus DRDID's per-group design checks (dr/reg outcome designs, dr/ipw PS).
    apply_rcond = (not callable(est_method)) and est_method in ("dr", "reg")
    apply_guard = not callable(est_method)
    # Propensity-overlap guard (matches R): dr/ipw only.
    apply_overlap = (not callable(est_method)) and est_method in ("dr", "ipw")
    cache_nuisance = (panel and apply_overlap and len(xcov_cols) > 1
                      and any(est_panel is original for original in _CACHEABLE_PANEL_ESTIMATORS)
                      and all(current is original for current, original in zip(
                          (overlap_check_fail, rcond_check_fail, drdid_design_singular),
                          _ORIGINAL_PANEL_GUARDS))
                      and drdid_trim._panel_pscore_cache_available())

    # ─── Precompute period-indexed structures ───────────────────────────────
    if panel:
        # per-unit cohort (from the first period) and a fixed global unit order
        first_frame = data[data[tname] == tlist[0]]
        unit_ids = first_frame[idname].to_numpy()
        unit_g = first_frame[gname].to_numpy()
        global_pos = pd.Series(np.arange(len(unit_ids)), index=unit_ids)

        period_y, period_w, period_cov = {}, {}, {}
        for p, tp in enumerate(tlist):
            fr = data[data[tname] == tp]
            fid = fr[idname].to_numpy()
            # Align once to the global unit order, including when a caller has
            # reordered the preprocessed frame. Cells then select NumPy rows.
            period_y[p] = pd.Series(fr[yname].to_numpy(), index=fid).reindex(unit_ids).to_numpy()
            period_w[p] = pd.Series(fr['w'].to_numpy(), index=fid).reindex(unit_ids).to_numpy()
            period_cov[p] = pd.DataFrame(fr[xcov_cols].to_numpy(), index=fid).reindex(unit_ids).to_numpy()
        if cache_nuisance:
            cov_keys = _period_cache_keys(period_cov)
            weight_keys = _period_cache_keys(period_w)
    else:
        # repeated cross sections / unbalanced panel: one global design matrix
        full_cov = pd.DataFrame(data[xcov_cols].to_numpy(), index=data.index)
        rowid_all = data['rowid'].to_numpy()
        rowid_unique = data['rowid'].unique()

    # ─── Main loop over groups × periods ────────────────────────────────────
    for g_index, g in enumerate(glist):
        state['g'] = g
        _pretg = last_pretreatment_index(g, tlist, anticipation)
        pret_g = _pretg if _pretg is not None else 0
        # Keep only the most recent design, so retained fit memory is bounded
        # independently of the number of periods. Entries cannot leak across
        # treatment cohorts or public fit calls. Intercept-only designs retain
        # the uncached path because their overlap guard has no warning context.
        nuisance_cache = ({'overlap_key': None, 'pscore_key': None, 'cov_keys': cov_keys,
                           'weight_keys': weight_keys} if cache_nuisance else None)

        if not panel:
            data = data.assign(G_m=1 * (data[gname] == g))

        for t_i in range(tlist_len):
            tn = tlist[t_i + tfac]
            state['tn'] = tn
            pret, action = plan_cell(g, t_i, tlist, anticipation, base_period, tfac)
            if action == BREAK:
                break

            tn_idx = t_i + tfac
            post_treat = 1 * (g <= tn)
            if action == BASE_ZERO:
                add_att_data(att=0, pst=post_treat, inf_f=np.zeros(n))
                continue

            pret_year = tlist[pret]

            if panel:
                _att_gt_panel_cell(
                    g, tn, tn_idx, pret, pret_g, pret_year, post_treat,
                    unit_ids, unit_g, never_treated, anticipation, tlist, tfac,
                    period_y, period_w, period_cov, global_pos, n,
                    fix_weights, est_panel, add_att_data, t_i,
                    apply_rcond, apply_guard, est_method, apply_overlap, nuisance_cache,
                    weighted=weights_name is not None)
            else:
                _att_gt_rc_cell(
                    g, tn, tn_idx, pret, pret_g, pret_year, post_treat,
                    data, tname, idname, yname, gname, never_treated, anticipation,
                    tlist, tfac, full_cov, rowid_unique, n, fix_weights,
                    est_rc, add_att_data, t_i, apply_rcond, apply_guard, est_method,
                    apply_overlap, xcov_cols, weighted=weights_name is not None)

    output = {'group': group, 'year': year, "att": att_est, 'post': post_array}
    if compute_inffunc:
        return (output, np.vstack(inf_func))
    return (output, None)


def _att_gt_panel_cell(g, tn, tn_idx, pret, pret_g, pret_year, post_treat,
                       unit_ids, unit_g, never_treated, anticipation, tlist, tfac,
                       period_y, period_w, period_cov, global_pos, n,
                       fix_weights, est_att_f, add_att_data, t_i, apply_rcond=False,
                       apply_guard=False, est_method="dr", apply_overlap=False,
                       nuisance_cache=None, weighted=False):
    Gm = (unit_g == g)
    if never_treated:
        Cm = (unit_g == 0)
    else:
        thresh = tlist[max(t_i, pret) + tfac] + anticipation
        Cm = (unit_g == 0) | ((unit_g > thresh) & (unit_g != g))
    mask = Gm | Cm
    uids = unit_ids[mask]
    if len(uids) == 0:
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return
    G = Gm[mask].astype(float)

    earlier_idx = min(tn_idx, pret)
    ypre = period_y[pret][mask]
    ypost = period_y[tn_idx][mask]
    # Pandas' former row reindexing yielded column-contiguous covariates.
    # Retain that layout so estimator reductions keep their rounding behavior.
    covariates = np.asfortranarray(period_cov[earlier_idx][mask])

    if fix_weights == "base_period":
        weight_idx = pret_g
    elif fix_weights == "first_period":
        weight_idx = 0
    else:
        weight_idx = earlier_idx
    w = period_w[weight_idx][mask]

    pscore_cache = None
    if nuisance_cache is not None:
        # Treatment labels are fixed within the cohort; the exact mask captures
        # changing not-yet-treated controls, including gaps between cohorts.
        overlap_key = (nuisance_cache['cov_keys'][earlier_idx], mask.tobytes(),
                       nuisance_cache['weight_keys'][weight_idx])
        pscore_key = (overlap_key, nuisance_cache['weight_keys'][weight_idx])
        if pscore_key != nuisance_cache['pscore_key']:
            nuisance_cache['pscore_key'] = pscore_key
            nuisance_cache['pscore'] = {}
        pscore_cache = nuisance_cache['pscore']

    # Propensity-overlap guard (matches R), checked before the singular guard.
    overlap_failed = False
    if apply_overlap:
        if nuisance_cache is None or overlap_key != nuisance_cache['overlap_key']:
            overlap_failed = overlap_check_fail(covariates, G, i_weights=w if weighted else None)
            if nuisance_cache is not None:
                nuisance_cache['overlap_key'] = overlap_key
                nuisance_cache['overlap_failed'] = overlap_failed
        else:
            overlap_failed = nuisance_cache['overlap_failed']
            # The original guard enters this context even for a clean fit.
            # Keep its warning-registry behavior when skipping that fit.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
    if overlap_failed:
        warnings.warn(
            f"overlap condition violated for group {g} in time period {tn}")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    # Regression-feasibility guard (matches R): control pre-check + DRDID's
    # per-group designs (incl. the treated outcome regression).
    if (apply_rcond and rcond_check_fail(covariates, G, i_weights=w if weighted else None)) or (
            apply_guard and drdid_design_singular(
                covariates, w, G, est_method, post=None)):
        warnings.warn(
            f"Covariate matrix for control units is singular or numerically "
            f"ill-conditioned for group {g} in time period {tn}; consider "
            f"centering/rescaling covariates or removing collinear terms")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    try:
        if pscore_cache is None:
            att_gt, att_inf_func = est_att_f(ypost, ypre, G, i_weights=w, covariates=covariates)
        else:
            att_gt, att_inf_func = est_att_f(
                ypost, ypre, G, i_weights=w, covariates=covariates,
                _pscore_cache=pscore_cache)
    except Exception as e:
        warnings.warn(
            f"Estimation failed for group {g} in time period {tn}: {e}. "
            f"Setting ATT to NA for this cell.")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    if np.isnan(att_gt):
        warnings.warn(f"ATT for (g,t) = ({g},{tn}) is NaN; setting it to NA")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    n1 = len(uids)
    inf_zeros = np.zeros(n)
    positions = np.flatnonzero(mask)
    inf_zeros[positions] = (n / n1) * att_inf_func
    add_att_data(att_gt, pst=post_treat, inf_f=inf_zeros)


def _att_gt_rc_cell(g, tn, tn_idx, pret, pret_g, pret_year, post_treat,
                    data, tname, idname, yname, gname, never_treated, anticipation,
                    tlist, tfac, full_cov, rowid_unique, n, fix_weights,
                    est_att_f, add_att_data, t_i, apply_rcond=False,
                    apply_guard=False, est_method="dr", apply_overlap=False,
                    xcov_cols=None, weighted=False):
    if never_treated:
        C_main = (data[gname] == 0)
    else:
        thresh = tlist[max(t_i, pret) + tfac] + anticipation
        C_main = (data[gname] == 0) | ((data[gname] > thresh) & (data[gname] != g))
    data = data.assign(C=1 * C_main)

    two_period = data[(data[tname] == tlist[tn_idx]) | (data[tname] == tlist[pret])]
    right_ids = two_period.loc[two_period.G_m.eq(1) | two_period.C.eq(1), 'rowid'].to_numpy()
    dis_idx = (data['rowid'].isin(right_ids)) & (data[tname].isin([tlist[tn_idx], tlist[pret]]))
    disdat = data.loc[dis_idx]

    G = disdat.G_m.to_numpy()
    C = disdat.C.to_numpy()
    Y = disdat[yname].to_numpy()
    post = 1 * (disdat[tname] == tlist[tn_idx]).to_numpy()
    w = disdat.w.to_numpy()

    if fix_weights in ("base_period", "first_period"):
        target_period = tlist[pret_g] if fix_weights == "base_period" else tlist[0]
        sub = data[data[tname] == target_period]
        wmap = sub.drop_duplicates(subset=[idname]).set_index(idname)['w']
        w = disdat[idname].map(wmap).to_numpy()
        valid = ~np.isnan(w)
        if not valid.all():
            n_dropped = disdat.loc[~valid, idname].nunique()
            warnings.warn(
                f"Dropped {n_dropped} units not observed in {fix_weights} "
                f"(period {target_period}) for group {g} in time period {tn}.")
            disdat = disdat[valid]
            G, C, Y, post, w = G[valid], C[valid], Y[valid], post[valid], w[valid]
            right_ids = disdat['rowid'].to_numpy()

    n1 = sum(G + C)

    skip = (np.sum(G * post) == 0 or np.sum(G * (1 - post)) == 0
            or np.sum(C * post) == 0 or np.sum(C * (1 - post)) == 0)
    if skip:
        warnings.warn(f"No (treated or control) units for group {g} in time period {tn}; setting ATT to NA.")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    # Balanced panel force-routed through RC by fix_weights='varying' (balanced iff
    # every unit observed in every period: nrow == n_units * n_periods). Computed
    # here (not just below) because it also gates the covariate remap.
    varying_forced_rc_balanced = (
        fix_weights == "varying" and len(data) == n * data[tname].nunique())

    # F-O1-2 / compute.att_gt.R:526-540: fix_weights only changes weights, not the
    # covariate conditioning set. In the forced-RC "varying" balanced branch, R
    # remaps every stacked row to its unit's earlier-period -- min(pret, t+tfac) --
    # covariate (earlier_idx_v; cov_early[id_map,]). The default RC path keeps each
    # row's own-period covariate, which diverges from R when the covariate is
    # time-varying. With a time-invariant covariate (or none) earlier-period ==
    # own-period, so this is a no-op there.
    if varying_forced_rc_balanced:
        earlier_idx = min(pret, tn_idx)
        earlier_period = tlist[earlier_idx]
        cov_early = (
            data[data[tname] == earlier_period]
            .drop_duplicates(subset=[idname])
            .set_index(idname)[xcov_cols])
        covariates = cov_early.reindex(disdat[idname]).to_numpy()
    else:
        covariates = full_cov.loc[disdat.index].to_numpy()

    if apply_rcond and rc_outcome_has_too_few_obs(
            G, post, w, covariates.shape[1], include_treated=est_method == "dr"):
        warnings.warn(
            f"Outcome regression design is singular (fewer positive-weight observations "
            f"than covariates) for group {g} in time period {tn}; "
            f"the ATT for this cell will be set to NA.")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    # Propensity-overlap guard (matches R), checked before the singular guard.
    if apply_overlap and overlap_check_fail(covariates, G, i_weights=w if weighted else None):
        warnings.warn(
            f"overlap condition violated for group {g} in time period {tn}")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    # Regression-feasibility guard (matches R): control pre-check + DRDID's
    # per-period, per-group designs (incl. the treated outcome regression).
    if (apply_rcond and rcond_check_fail(covariates, G, i_weights=w if weighted else None)) or (
            apply_guard and drdid_design_singular(
                covariates, w, G, est_method, post=post)):
        warnings.warn(
            f"Covariate matrix for control units is singular or numerically "
            f"ill-conditioned for group {g} in time period {tn}; consider "
            f"centering/rescaling covariates or removing collinear terms")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    try:
        att_gt, att_inf_func = est_att_f(y=Y, post=post, D=G, i_weights=w, covariates=covariates)
    except Exception as e:
        warnings.warn(
            f"Estimation failed for group {g} in time period {tn}: {e}. "
            f"Setting ATT to NA for this cell.")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    if np.isnan(att_gt):
        warnings.warn(f"ATT for (g,t) = ({g},{tn}) is NaN; setting it to NA")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    # n/n1 rescaling (matches R `did`; mirrors compute_att_gt). Scale the per-cell
    # IF by #units only for a BALANCED panel force-routed through the RC estimator
    # by fix_weights='varying' (R reshapes wide; two rows per unit summed back per
    # unit below, so the row count would halve every SE/band). A genuinely
    # unbalanced panel must scale by #rows (= n1) like R's long/RC path; using
    # #units there over-inflates the IF by n_rows/n_units (audit-v3 F1).
    # (varying_forced_rc_balanced computed above, where it also gates the cov remap.)
    n1_scale = len(np.unique(right_ids)) if varying_forced_rc_balanced else n1
    # R `did` 2.5.1 (commit 4e9de53): fold pre+post per unit must divide by 2 (the
    # RC IF is normalized over 2*n_units stacked rows) so the balanced-panel
    # fix_weights='varying' SE matches the panel normalization; without it the SE is
    # exactly 2x too large (R 2.5.0's bug). Mirrors compute_att_gt.
    fold_factor = 0.5 if varying_forced_rc_balanced else 1.0
    att_inf_func = fold_factor * (n / n1_scale) * att_inf_func
    inf_func_df = pd.DataFrame({"inf_func": att_inf_func, "right_ids": right_ids})
    inf_zeros = np.zeros(n)
    aggte_if = inf_func_df.groupby('right_ids').inf_func.sum()
    # Preserve missing observed contributions; absent units remain zero below.
    missing_if = pd.isna(att_inf_func)
    if np.any(missing_if):
        aggte_if.loc[np.unique(right_ids[missing_if])] = np.nan
    try:
        dis_idx1 = np.isin(data['rowid'].unique(), aggte_if.index.to_numpy())
    except Exception:
        dis_idx1 = np.isin(data['rowid'].unique().to_numpy(), aggte_if.index.to_numpy())
    inf_zeros[dis_idx1] = np.array(aggte_if)
    add_att_data(att_gt, pst=post_treat, inf_f=inf_zeros)
