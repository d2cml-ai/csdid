import numpy as np, pandas as pd
import patsy 
from drdid import reg_did
from csdid.attgt_fnc import drdid_trim
from csdid.attgt_fnc.compute_att_gt_shared import (
    select_estimators, last_pretreatment_index, plan_cell, rcond_check_fail,
    drdid_design_singular, overlap_check_fail, BREAK, BASE_ZERO)

from csdid.utils.bmisc import panel2cs2
import warnings


fml = patsy.dmatrices
# Initialize a list to store data for each iteration
results_list = []

def compute_att_gt(dp, est_method = "dr", base_period = 'varying', compute_inffunc = True):
    yname = dp['yname']
    tname = dp['tname']
    idname = dp['idname']
    xformla = dp['xformla']
    data = dp['data']
    weights_name = dp['weights_name']
    # base_period = dp['base_period']
    panel = dp['panel']
    true_rep_cross_section = dp['true_rep_cross_section']
    control_group = dp['control_group']
    anticipation = dp['anticipation']
    gname = dp['gname']
    n = dp['n']
    nT = dp['nT']
    nG = dp['nG']
    tlist = dp['tlist']
    glist = dp['glist']
    fix_weights = dp.get('fix_weights')
    # Pre-built (globally factor-consistent) design column names; selected per
    # cell instead of re-running patsy (which broke transformed/factor covariates).
    xcov_cols = dp.get('xcov_cols')

    # Per-period weight lookup (idname -> weight in a given period), used by
    # fix_weights='base_period'/'first_period' to fix each unit's weight to a
    # specific period rather than the per-cell earlier period (matches R `did`).
    def _period_weight_map(target_period):
        sub = data[data[tname] == target_period]
        return sub.drop_duplicates(subset=[idname]).set_index(idname)['w']

    # Calculate time periods and adjustment factor
    tlist_len = len(tlist) - 1 if base_period != "universal" else len(tlist)
    tfac = 1 if base_period != "universal" else 0

    inf_func = []

    att_est, group, year, post_array = [], [], [], []

    def add_att_data(att = 0, pst = 0, inf_f = []):
        if compute_inffunc:
            inf_func.append(inf_f)
        att_est.append(att)
        group.append(g)
        year.append(tn)
        post_array.append(pst)

    # Handle never treated case
    never_treated = control_group == 'nevertreated'
    if never_treated:
        data['C'] = (data[gname] == 0).astype(int)
    data['y_main'] = data[yname]

    # 2x2 estimators (panel + repeated cross sections), selected once.
    est_panel, est_rc = select_estimators(est_method)

    # Distinguish a BALANCED panel force-routed through the RC estimators by
    # fix_weights="varying" (preprocess sets panel=False *before* its imbalance
    # check, so the balanced/unbalanced distinction is otherwise lost) from a
    # GENUINELY unbalanced panel. R `did` scales the per-cell influence function
    # by #units in the former (it reshapes wide) and by #rows in the latter (it
    # stays long / RC). A balanced panel has every unit observed in every period,
    # i.e. nrow == n_units * n_periods (R's own imbalance test, preprocess L552).
    # Without this, fix_weights="varying" on a real unbalanced panel over-scaled
    # the IF by n_rows/n_units, inflating every SE/band (audit-v3 F1).
    varying_forced_rc_balanced = (
        fix_weights == "varying" and not true_rep_cross_section
        and len(data) == n * data[tname].nunique())

    # Regression-feasibility guards (match R). `did`'s pooled control pre-check
    # applies to outcome-regression methods (dr/reg); DRDID's internal per-group
    # design checks apply to any built-in method (dr/reg outcome designs, dr/ipw
    # PS design). Custom (callable) estimators are exempt.
    apply_rcond = (not callable(est_method)) and est_method in ("dr", "reg")
    apply_guard = not callable(est_method)
    # Propensity-overlap guard (matches R): dr/ipw fit a PS logit, so a
    # (near-)separated cell must be NA'd rather than estimated from exploding IPW
    # weights. reg / custom estimators are exempt.
    apply_overlap = (not callable(est_method)) and est_method in ("dr", "ipw")

    # Loop over groups
    for g_index, g in enumerate(glist):  
        # Set up .G once
        # Create a binary column 'G_m' to indicate if a row belongs to the current group 'g'
        G_main = (data[gname] == glist[g_index])
        data = data.assign(G_m=1 * G_main)

        # Group's last pre-treatment period index (R's pret_g), used as the
        # fixed reference period for fix_weights='base_period'.
        _pretg = last_pretreatment_index(g, tlist, anticipation)
        pret_g = _pretg if _pretg is not None else 0

        # Loop over time periods
        for t_i in range(tlist_len):

            # Resolve the pretreatment period and control flow. This
            # selection logic is shared with the faster_mode path
            # (compute_att_gt2) via plan_cell so the two cannot drift.
            tn = tlist[t_i + tfac]
            pret, _action = plan_cell(g, t_i, tlist, anticipation, base_period, tfac)
            if _action == BREAK:
                break

            pret_year = tlist[pret]
            post_treat = 1 * (g <= tn)
            if _action == BASE_ZERO:
                add_att_data(att=0, pst=post_treat, inf_f=np.zeros(n))
                continue

            # Control-group indicator for not-yet-treated comparisons. Using
            # the final pret matches the legacy ordering: max(t_i, pret)
            # collapses to t_i for post-treatment cells, and pret is unchanged
            # for pre-treatment cells.
            if not never_treated:
                n1 = (data[gname] == 0)
                n2 = (data[gname] > (tlist[max(t_i, pret) + tfac] + anticipation))
                n3 = (data[gname] != glist[g_index])
                row_eval = n1 | (n2 & n3)
                data = data.assign(C=1 * row_eval)

            # Subset the data for the current and pretreatment periods
            disdat = data[(data[tname] == tn) | (data[tname] == tlist[pret])]
            # print("Shape of disdat:", disdat.shape)

        # results for the case with panel data
        #-----------------------------------------------------------------------------

            if panel:
                disdat = panel2cs2(disdat, yname, idname, tname)
                disdat = disdat.dropna()
                n_cell = len(disdat)
                dis_idx = np.array(disdat.G_m == 1) | np.array(disdat.C == 1)
                disdat = disdat.loc[dis_idx, :]
                n1 = len(disdat)
                G = disdat.G_m
                C = disdat.C
                w = disdat.w

                # fix_weights: replace the per-cell (earlier-period) weights with
                # each unit's weight from a fixed period (matches R `did`). The
                # panel path only runs on balanced panels, so every unit is present
                # in the target period.
                if fix_weights in ("base_period", "first_period"):
                    target_period = tlist[pret_g] if fix_weights == "base_period" else tlist[0]
                    wmap = _period_weight_map(target_period)
                    w = disdat[idname].map(wmap)

                ypre = disdat.y0 if tn > pret_year else disdat.y1
                ypost = disdat.y0 if tn < pret_year else disdat.y1
                covariates = disdat[xcov_cols].to_numpy()

                G, C, w, ypre = map(np.array, [G, C, w, ypre])
                ypost, covariates = map(np.array, [ypost, covariates])

                est_att_f = est_panel

                # Propensity-overlap guard (matches R, checked BEFORE the singular
                # guard so its precedence/message match): NA a dr/ipw cell whose PS
                # logit is (near-)separated (max fitted PS >= 0.999).
                if apply_overlap and overlap_check_fail(covariates, G):
                    warnings.warn(
                        f"overlap condition violated for group {g} in time period {tn}"
                    )
                    add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
                    continue

                # Regression-feasibility guard (matches R): NA the cell when a 2x2
                # estimating design is singular/ill-conditioned, rather than letting
                # the solver return an untrustworthy estimate. Covers the control
                # design (`did` pre-check) and DRDID's per-group designs (incl. the
                # treated outcome regression -- a tiny treated cohort can be singular
                # even when the control design is fine).
                if (apply_rcond and rcond_check_fail(covariates[G == 0])) or (
                        apply_guard and drdid_design_singular(
                            covariates, w, G, est_method, post=None)):
                    warnings.warn(
                        f"Covariate matrix for control units is singular or numerically "
                        f"ill-conditioned for group {g} in time period {tn}; consider "
                        f"centering/rescaling covariates or removing collinear terms"
                    )
                    add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
                    continue

                try:
                    att_gt, att_inf_func = est_att_f(ypost, ypre, G, i_weights=w, covariates=covariates)
                except Exception as e:
                    warnings.warn(
                        f"Estimation failed for group {g} in time period {tn}: {e}. "
                        f"Setting ATT to NA for this cell."
                    )
                    add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
                    continue

                inf_zeros = np.zeros(n)
                att_inf = n_cell / n1 * att_inf_func
                inf_zeros[dis_idx] = att_inf

                add_att_data(att_gt, pst=post_treat, inf_f=inf_zeros)

        #-----------------------------------------------------------------------------
        # results for the case with no panel data
        #-----------------------------------------------------------------------------
            
            if not panel:
                # Fixed right_ids selection
                right_ids = disdat.loc[disdat.G_m.eq(1) | disdat.C.eq(1), 'rowid'].to_numpy()
                # print("Shape of ids:",right_ids.shape)  # Show dimensions
                # Consistent time period selection
                dis_idx = (data['rowid'].isin(right_ids)) & \
                            (data[tname].isin([tlist[t_i + tfac], tlist[pret]]))

                disdat = data.loc[dis_idx]

                G = disdat.G_m.to_numpy()
                C = disdat.C.to_numpy()
                Y = disdat[yname].to_numpy()
                post = 1 * (disdat[tname] == tlist[t_i + tfac]).to_numpy()
                w = disdat.w.to_numpy()

                # fix_weights for the RC path: 'varying' keeps per-period weights
                # (already in disdat.w); 'base_period'/'first_period' fix each row to
                # its unit's weight from the target period (matches R `did`).
                if fix_weights in ("base_period", "first_period"):
                    target_period = tlist[pret_g] if fix_weights == "base_period" else tlist[0]
                    wmap = _period_weight_map(target_period)
                    w = disdat[idname].map(wmap).to_numpy()
                    # Drop units not observed in the target period (matches R `did`).
                    valid = ~np.isnan(w)
                    if not valid.all():
                        n_dropped = disdat.loc[~valid, idname].nunique()
                        warnings.warn(
                            f"Dropped {n_dropped} units not observed in {fix_weights} "
                            f"(period {target_period}) for group {g} in time period {tn}."
                        )
                        disdat = disdat[valid]
                        G, C, Y, post, w = G[valid], C[valid], Y[valid], post[valid], w[valid]
                        right_ids = disdat['rowid'].to_numpy()

                n1 = sum(G + C)

                # Store the current iteration's data
                current_data = {
                    'Y': Y,
                    'post': post,
                    'G': G,
                    'group': g,
                    'time_period': tn
                }
                
                # print(f"Lengths: Y={len(Y)}, post={len(post)}, G={len(G)}, group={type(g)}, time_period={type(tn)}")
                # results_list.append(pd.DataFrame(current_data))

                #-----------------------------------------------------------------------------
                # checks to make sure that we have enough observations
                
                skip_this_att_gt = False

                if np.sum(G * post) == 0:
                    warnings.warn(f"No units in group {g} in time period {tn}; setting ATT to NA.")
                    skip_this_att_gt = True 

                if np.sum(G * (1 - post)) == 0:
                    warnings.warn(f"No units in group {g} in the base period for time {tn}; setting ATT to NA.")
                    skip_this_att_gt = True 

                if np.sum(C * post) == 0:
                    warnings.warn(f"No available control units for group {g} in time period {tn}; setting ATT to NA.")
                    skip_this_att_gt = True 

                if np.sum(C * (1 - post)) == 0:
                    warnings.warn(f"No available control units for group {g} in the base period for time {tn}; setting ATT to NA.")
                    skip_this_att_gt = True 

                if skip_this_att_gt:
                    # Append results with missing ATT and NA influence function
                    add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
                    #add_att_data()
                    continue

                # return (inf_func)
                # F-O1-2 / compute.att_gt.R:526-540: fix_weights only changes
                # weights, not the covariate conditioning set. In the forced-RC
                # "varying" balanced branch, R remaps every stacked row to its
                # unit's earlier-period -- min(pret, t+tfac) -- covariate
                # (earlier_idx_v; cov_early[id_map,]). The default RC path keeps
                # each row's own-period covariate, which diverges from R when the
                # covariate is time-varying. With a time-invariant covariate (or
                # none) earlier-period == own-period, so this is a no-op there.
                if varying_forced_rc_balanced:
                    earlier_idx = min(pret, t_i + tfac)
                    earlier_period = tlist[earlier_idx]
                    cov_early = (
                        data[data[tname] == earlier_period]
                        .drop_duplicates(subset=[idname])
                        .set_index(idname)[xcov_cols])
                    covariates = cov_early.reindex(disdat[idname]).to_numpy()
                else:
                    covariates = disdat[xcov_cols].to_numpy()

                #-----------------------------------------------------------------------------
                # code for actually computing att(g,t)
                #-----------------------------------------------------------------------------
                
                est_att_f = est_rc

                # Propensity-overlap guard (matches R), checked before the singular
                # guard. On the RC / varying path the PS logit is fit on the stacked
                # long sample (covariates, G as built here).
                if apply_overlap and overlap_check_fail(covariates, G):
                    warnings.warn(
                        f"overlap condition violated for group {g} in time period {tn}"
                    )
                    add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
                    continue

                # Regression-feasibility guard (matches R): on the RC / unbalanced
                # path DRDID fits control AND treated outcome regressions separately
                # pre and post, so each period-by-group design is checked (this is
                # what NAs a tiny treated cohort with too few post observations).
                if (apply_rcond and rcond_check_fail(covariates[G == 0])) or (
                        apply_guard and drdid_design_singular(
                            covariates, w, G, est_method, post=post)):
                    warnings.warn(
                        f"Covariate matrix for control units is singular or numerically "
                        f"ill-conditioned for group {g} in time period {tn}; consider "
                        f"centering/rescaling covariates or removing collinear terms"
                    )
                    add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
                    continue

                try:
                    att_gt, att_inf_func = est_att_f(y=Y, post=post, D = G, i_weights=w, covariates=covariates)
                except Exception as e:
                    warnings.warn(
                        f"Estimation failed for group {g} in time period {tn}: {e}. "
                        f"Setting ATT to NA for this cell."
                    )
                    add_att_data(att=np.nan, pst=post_treat, inf_f=np.full(n, np.nan))
                    continue

                # n/n1 rescaling (matches R `did`). Two regimes, keyed on whether
                # this is a BALANCED panel force-routed through the RC estimator by
                # fix_weights='varying' (preprocess sets panel=False), or a genuine
                # RC / unbalanced panel:
                #   * balanced-through-RC: R reshapes wide and scales by #units. The
                #     stacked long 2x2 has TWO rows per unit; the per-row IF is
                #     summed back per unit (groupby right_ids below), so scaling by
                #     the row count would halve every SE/band. Use #units.
                #   * genuine RC / unbalanced: R stays long and scales by #rows; a
                #     unit may contribute one OR two rows, so #units != #rows and
                #     scaling by #units over-inflates the IF by n_rows/n_units
                #     (audit-v3 F1). Use #rows (= n1).
                n1_scale = len(np.unique(right_ids)) if varying_forced_rc_balanced else n1
                # R `did` 2.5.1 (commit 4e9de53): the RC influence function is
                # normalized over the 2*n_units stacked rows; folding pre+post per
                # unit (the groupby-sum below) must divide by 2 so the balanced-panel
                # fix_weights='varying' SE matches the panel normalization. Without
                # the 0.5 the SE is exactly 2x too large -- R 2.5.0's bug, which the
                # port had inherited.
                fold_factor = 0.5 if varying_forced_rc_balanced else 1.0
                att_inf_func = fold_factor * (n/n1_scale) * att_inf_func

                inf_func_df = pd.DataFrame(
                {
                    "inf_func": att_inf_func,
                    "right_ids": right_ids
                }
                ).fillna(0)

                inf_zeros = np.zeros(n)
                aggte_infffuc = inf_func_df.groupby('right_ids').inf_func.sum()
                try:
                    dis_idx1 = np.isin(data['rowid'].unique(), aggte_infffuc.index.to_numpy())
                except Exception:
                    dis_idx1 = np.isin(data['rowid'].unique().to_numpy(), aggte_infffuc.index.to_numpy())
                
                inf_zeros[dis_idx1] = np.array(aggte_infffuc)

                add_att_data(att_gt, pst = post_treat, inf_f=inf_zeros)
                # print(att_est)

    output = {
    'group': group ,
    'year': year,
    "att" : att_est,
    'post': post_array
    }

    if compute_inffunc:
        return (output, np.vstack(inf_func))
    else:
        return (output, None)


# import numpy as np, pandas as pd
# import patsy 
# from drdid import drdid, reg_did, ipwd_did

# from csdid.utils.bmisc import panel2cs2

# fml = patsy.dmatrices


# # data Struct
# # output = {
# #   "att" : []
# #   'group': []
# #   'year': []
# #   'post ': []
# # }


# def compute_att_gt(dp, est_method = "dr", base_period = 'varying'):
#   yname = dp['yname']
#   tname = dp['tname']
#   idname = dp['idname']
#   xformla = dp['xformla']
#   data = dp['data']
#   weights_name = dp['weights_name']
#   # base_period = dp['base_period']
#   panel = dp['panel']
#   true_rep_cross_section = dp['true_rep_cross_section']
#   control_group = dp['control_group']
#   anticipation = dp['anticipation']
#   gname = dp['gname']
#   n = dp['n']
#   nT = dp['nT']
#   nG = dp['nG']
#   tlist = dp['tlist']
#   glist = dp['glist']
#   tlist_len, tfac = len(tlist), 0
#   if base_period != 'universal':
#     tlist_len = tlist_len - 1
#     tfac = 1

#   inf_func = []

#   att_est, group, year, post_array = [], [], [], []

#   def add_att_data(att = 0, pst = 0, inf_f = []):
#     inf_func.append(inf_f)
#     att_est.append(att)
#     group.append(g)
#     year.append(tn)
#     post_array.append(pst)

#   never_treated = control_group == 'nevertreated'
#   if never_treated:
#     data = data.assign(C = 1 * (data[gname] == 0))
#   data = data.assign(y_main = data[yname])

#   # g, t = glist[0], tlist[0]

#   # for _, g, in enumerate(glist):
#   for g_index, g in enumerate(glist):

#     # g = glist[1]
#     G_main = (data[gname] == g)
#     data = data.assign(G_m = 1 * G_main)

#     for t_i in range(tlist_len):
#       pret = t_i
#       tn = tlist[t_i + tfac]
#       if base_period == 'universal' or g < tn:
#         try:
#           pret = np.where(tlist + anticipation < g)[0][-1]
#         except:
#           raise f"There are no pre-treatment periods for the group first treated at {g}\nUnits from this group are dropped"
#           # break


#       if base_period == 'universal':
#         if pret == tn:
#           add_att_data()

#       if not never_treated:
#         n1 = data[gname] == 0
#         n2 = (data[gname] > (tlist[np.max([t_i, pret]) + tfac]) + anticipation)
#         #n3 = np.where(data[gname] != glist[g], True, False)
#         n3 = np.where(data[gname] != glist[g_index], True, False)

#         row_eval = n1 | n2 & n3
#         data = data.assign(C = 1 * row_eval)

#       post_treat = 1 * (g <= tn)
#       disdat = data[(data[tname] == tn) | (data[tname] == tlist[pret])]

#       if panel: 
#         disdat = panel2cs2(disdat, yname, idname, tname)
#         disdat = disdat.dropna()
#         n = len(disdat)
#         dis_idx = np.array(disdat.G_m == 1) | np.array(disdat.C == 1)
#         disdat = disdat.loc[dis_idx, :]
#         n1 = len(disdat)
#         G = disdat.G_m
#         C = disdat.C
#         w = disdat.w

#         ypre = disdat.y0 if tn > pret else disdat.y1
#         ypost = disdat.y0 if tn < pret else disdat.y1
#         _, covariates = fml(xformla, data = disdat, return_type = 'dataframe')

#         G, C, w, ypre = map(np.array, [G, C, w, ypre])
#         ypost, covariates = map(np.array, [ypost, covariates])

#         if callable(est_method):
#           est_att_f = est_method
#         elif est_method == "reg":
#           est_att_f = reg_did.reg_did_panel
#         elif est_method == "ipw":
#           est_att_f = ipwd_did.std_ipw_did_panel
#         elif est_method == "dr":
#           est_att_f = drdid.drdid_panel

#         att_gt, att_inf_func = est_att_f(ypost, ypre, G, i_weights=w, covariates=covariates)

#         inf_zeros = np.zeros(n)
#         att_inf = n / n1 * att_inf_func
#         inf_zeros[dis_idx] = att_inf

#         add_att_data(att_gt, inf_f=inf_zeros)

#       if not panel:
#         right_ids = np.array(disdat.query('(G_m == 1) or (C == 1)').rowid.to_numpy())
#         dis_idx = (data['rowid'].isin(right_ids)) &\
#           ((data[tname] == tlist[t_i + tfac]) |\
#             (data[tname] == tlist[pret]))

#         disdat = data.loc[dis_idx]

#         G = disdat.G_m.to_numpy()
#         C = disdat.C.to_numpy()
#         Y = disdat[yname].to_numpy()
#         post = 1 * (disdat[tname] == tlist[t_i + tfac]).to_numpy()
#         w = disdat.w.to_numpy()

#         # G, C, Y, post, w = map(np.array, [G, C, Y, post, w])


#         n1 = sum(G + C)

#         skip_this_att_gt = False
#         if np.sum(G * post) == 0:
#           print(f"No units in group {g} in time period {tn}")
#           skip_this_att_gt = True 

#         if np.sum(G * (1 - post)) == 0:
#           print(f"No units in group {g} in time period {tn}")
#           skip_this_att_gt = True 

#         if np.sum(C * post) == 0:
#           print(f"No available control units for group {g} in time period {tn}")
#           skip_this_att_gt = True 

#         if np.sum(C * (1 - post)) == 0:
#           print(f"No available control units for group {g} in time period {tn}")
#           skip_this_att_gt = True 

#         if skip_this_att_gt:
#           add_att_data()

#         try:
#           _, covariates = fml(xformla, data = disdat, return_type = 'dataframe')
#           covariates = np.array(covariates)
#         except:
#           y_str, x_str = xformla.split("~")
#           xs1 = x_str.split('+')
#           xs1_col_names = [x.strip() for x in xs1 if x.strip() != '1']
#           n_dis = len(disdat)
#           ones = np.ones((n_dis, 1))
#           try:
#             covariates = disdat[xs1_col_names].to_numpy()
#             covariates = np.append(covariates, ones, axis=1)
#           except:
#             covariates = ones

        

#         if callable(est_method):
#           est_att_f = est_method
#         elif est_method == "reg":
#           est_att_f = reg_did.reg_did_rc
#         elif est_method == "ipw":
#           est_att_f = ipwd_did.std_ipw_did_rc
#         elif est_method == "dr":
#           est_att_f = drdid.drdid_rc
#         att_gt, att_inf_func = est_att_f(y=Y, post=post, D = G, i_weights=w, covariates=covariates)

#         inf_func_df = pd.DataFrame(
#           {
#             "inf_func": att_inf_func,
#             "right_ids": right_ids
#           }
#         )
#         inf_zeros = np.zeros(n)
#         aggte_infffuc = inf_func_df.groupby('right_ids').inf_func.sum()
#         try:
#           dis_idx1 = np.isin(data['rowid'].unique(), aggte_infffuc.index.to_numpy())
#         except:
#           dis_idx1 = np.isin(data['rowid'].unique().to_numpy(), aggte_infffuc.index.to_numpy())
#         inf_zeros[dis_idx1] = np.array(aggte_infffuc)

#         add_att_data(att_gt, pst = post_treat, inf_f=inf_zeros)
#         # print(att_est)

#   output = {
#     'group': group ,
#     'year': year,
#     "att" : att_est,
#     'post ': post_array
#   }
#   return (output, np.array(inf_func))
