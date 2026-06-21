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
    select_estimators, last_pretreatment_index, plan_cell, BREAK, BASE_ZERO)
import warnings


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
            period_y[p] = pd.Series(fr[yname].to_numpy(), index=fid)
            period_w[p] = pd.Series(fr['w'].to_numpy(), index=fid)
            period_cov[p] = pd.DataFrame(fr[xcov_cols].to_numpy(), index=fid)
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
                    fix_weights, est_panel, add_att_data, t_i)
            else:
                _att_gt_rc_cell(
                    g, tn, tn_idx, pret, pret_g, pret_year, post_treat,
                    data, tname, idname, yname, gname, never_treated, anticipation,
                    tlist, tfac, full_cov, rowid_unique, n, fix_weights,
                    est_rc, add_att_data, t_i)

    output = {'group': group, 'year': year, "att": att_est, 'post': post_array}
    if compute_inffunc:
        return (output, np.vstack(inf_func))
    return (output, None)


def _att_gt_panel_cell(g, tn, tn_idx, pret, pret_g, pret_year, post_treat,
                       unit_ids, unit_g, never_treated, anticipation, tlist, tfac,
                       period_y, period_w, period_cov, global_pos, n,
                       fix_weights, est_att_f, add_att_data, t_i):
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
    ypre = period_y[pret].reindex(uids).to_numpy()
    ypost = period_y[tn_idx].reindex(uids).to_numpy()
    covariates = period_cov[earlier_idx].reindex(uids).to_numpy()

    if fix_weights == "base_period":
        w = period_w[pret_g].reindex(uids).to_numpy()
    elif fix_weights == "first_period":
        w = period_w[0].reindex(uids).to_numpy()
    else:
        w = period_w[earlier_idx].reindex(uids).to_numpy()

    try:
        att_gt, att_inf_func = est_att_f(ypost, ypre, G, i_weights=w, covariates=covariates)
    except Exception as e:
        warnings.warn(
            f"Estimation failed for group {g} in time period {tn}: {e}. "
            f"Setting ATT to NA for this cell.")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    n1 = len(uids)
    inf_zeros = np.zeros(n)
    positions = global_pos.reindex(uids).to_numpy()
    inf_zeros[positions] = (n / n1) * att_inf_func
    add_att_data(att_gt, pst=post_treat, inf_f=inf_zeros)


def _att_gt_rc_cell(g, tn, tn_idx, pret, pret_g, pret_year, post_treat,
                    data, tname, idname, yname, gname, never_treated, anticipation,
                    tlist, tfac, full_cov, rowid_unique, n, fix_weights,
                    est_att_f, add_att_data, t_i):
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

    covariates = full_cov.loc[disdat.index].to_numpy()

    try:
        att_gt, att_inf_func = est_att_f(y=Y, post=post, D=G, i_weights=w, covariates=covariates)
    except Exception as e:
        warnings.warn(
            f"Estimation failed for group {g} in time period {tn}: {e}. "
            f"Setting ATT to NA for this cell.")
        add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
        return

    att_inf_func = (n / n1) * att_inf_func
    inf_func_df = pd.DataFrame({"inf_func": att_inf_func, "right_ids": right_ids}).fillna(0)
    inf_zeros = np.zeros(n)
    aggte_if = inf_func_df.groupby('right_ids').inf_func.sum()
    try:
        dis_idx1 = np.isin(data['rowid'].unique(), aggte_if.index.to_numpy())
    except Exception:
        dis_idx1 = np.isin(data['rowid'].unique().to_numpy(), aggte_if.index.to_numpy())
    inf_zeros[dis_idx1] = np.array(aggte_if)
    add_att_data(att_gt, pst=post_treat, inf_f=inf_zeros)
