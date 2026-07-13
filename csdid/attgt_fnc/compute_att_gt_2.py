import numpy as np
import pandas as pd
from scipy import sparse
import patsy
from drdid import drdid, reg_did, ipwd_did
from csdid.utils.bmisc import panel2cs2
import warnings
from typing import Dict, List, Union, Optional

def compute_att_gt(dp: Dict) -> pd.DataFrame:
    """
    Compute Group-Time Average Treatment Effects.
    
    Args:
        dp: Dictionary containing DID parameters
            Required keys: data, yname, tname, idname, xformla, weights_name,
            est_method, base_period, panel, control_group, anticipation,
            gname, n, nT, nG, tlist, glist
    
    Returns:
        pd.DataFrame: Results containing treatment effects and related data
    """
    # Extract parameters from dp dictionary
    data = dp['data'].copy()
    yname = dp['yname']
    tname = dp['tname']
    idname = dp['idname']
    xformla = dp['xformla']
    weights_name = dp.get('weights_name')
    est_method = dp.get('est_method', 'dr')
    base_period = dp.get('base_period', 'varying')
    panel = dp['panel']
    control_group = dp['control_group']
    anticipation = dp['anticipation']
    gname = dp['gname']
    n = dp['n']
    nT = dp['nT']
    nG = dp['nG']
    tlist = dp['tlist']
    glist = dp['glist']

    # Initialize results DataFrame with empty structure
    results_df = pd.DataFrame(columns=['Y', 'post', 'G', 'group', 'time_period'])
    
    # Calculate time periods and adjustment factor
    tlist_len = len(tlist) - 1 if base_period != "universal" else len(tlist)
    tfac = 1 if base_period != "universal" else 0
    
    # Initialize influence function matrix
    inffunc = sparse.csr_matrix((n, nG*(nT-tfac)), dtype=float)
    inffunc_updates = []
    update_counter = 1
    counter = 1

    # Handle never treated case
    never_treated = control_group == 'nevertreated'
    if never_treated:
        data['C'] = (data[gname] == 0).astype(int)
    data['y_main'] = data[yname]
    
    # # Ensure rowid exists for tracking
    # if 'rowid' not in data.columns:
    #     data['rowid'] = range(1, len(data) + 1)  # 1-based indexing like R

    # Main loop over groups
    for g_idx, g in enumerate(glist, 1):
        # Set group indicator
        data['G'] = (data[gname] == g).astype(int)
        
        # Time period loop
        for t in range(1, tlist_len + 1):
            # Set base pre-treatment period
            pret = t
            current_period = tlist[t + tfac - 1]  # Adjust for 0-based indexing
            
            # Handle universal base period
            if base_period == "universal":
                pret_mask = (np.array(tlist) + anticipation) < g
                if not any(pret_mask):
                    warnings.warn(
                        f"There are no pre-treatment periods for the group first treated at {g}\n"
                        "Units from this group are dropped"
                    )
                    break
                pret = np.where(pret_mask)[0][-1]
            
            # Update control group for not-never-treated case
            if not never_treated:
                max_t = max(t, pret)
                threshold = tlist[max_t + tfac - 1] + anticipation
                future_treat = data[gname] > threshold
                not_current = data[gname] != g
                data['C'] = ((data[gname] == 0) | (future_treat & not_current)).astype(int)

            # Check if in post-treatment period
            if g <= current_period:
                pret_mask = (np.array(tlist) + anticipation) < g
                if not any(pret_mask):
                    warnings.warn(
                        f"There are no pre-treatment periods for the group first treated at {g}\n"
                        "Units from this group are dropped"
                    )
                    break
                pret = np.where(pret_mask)[0][-1]
                print("This is alex", pret)

            # Handle universal base period normalization
            if base_period == "universal" and tlist[pret] == current_period:
                # Add zero treatment effect and continue
                results_df = pd.concat([
                    results_df,
                    pd.DataFrame({
                        'Y': [0],
                        'post': [0],
                        'G': [1],
                        'group': [g],
                        'time_period': [current_period]
                    })
                ], ignore_index=True)
                continue

            # Get relevant time period data
            period_mask = data[tname].isin([current_period, tlist[pret]])
            disdat = data[period_mask].copy()

            if panel:
                # Handle panel data case
                disdat = panel2cs2(disdat, yname, idname, tname)
                n = len(disdat)
                treat_control_mask = (disdat['G'] == 1) | (disdat['C'] == 1)
                disdat = disdat[treat_control_mask].copy()
                n1 = len(disdat)

                # Extract variables
                G = disdat['G'].values
                C = disdat['C'].values
                Ypre = disdat['y0'] if current_period > tlist[pret] else disdat['y1']
                Ypost = disdat['y1'] if current_period > tlist[pret] else disdat['y0']
                w = disdat[weights_name].values if weights_name else np.ones(len(disdat))

                # Get covariates
                try:
                    _, covariates = patsy.dmatrices(xformla, data=disdat)
                    covariates = np.asarray(covariates)
                    
                    # Check estimation method conditions
                    if not callable(est_method):
                        if est_method in ['dr', 'ipw']:
                            from sklearn.linear_model import LogisticRegression
                            logit = LogisticRegression().fit(covariates, G)
                            pscores = logit.predict_proba(covariates)[:, 1]
                            if np.max(pscores) >= 0.999:
                                warnings.warn(f"Overlap condition violated for {g} in time period {current_period}")
                                raise ValueError("Overlap condition violated")

                except Exception:
                    covariates = np.ones((len(disdat), 1))

                # Select estimation method
                est_func = {
                    'reg': reg_did.reg_did_panel,
                    'ipw': ipwd_did.std_ipw_did_panel,
                    'dr': drdid.drdid_panel
                }.get(est_method, est_method if callable(est_method) else drdid.drdid_panel)

                # Compute treatment effects
                att_gt, att_inf_func = est_func(
                    Ypost, Ypre, G, i_weights=w, covariates=covariates
                )
                
                # Adjust influence function
                att_inf_func = (n/n1) * att_inf_func

            else:
                # Handle non-panel data case
                right_ids = disdat.loc[
                    (disdat['G'] == 1) | (disdat['C'] == 1), 'rowid'
                ].values
                
                dis_idx = data['rowid'].isin(right_ids) & period_mask
                disdat = data[dis_idx].copy()

                # Extract variables
                G = disdat['G'].values
                C = disdat['C'].values
                Y = disdat[yname].values
                post = (disdat[tname] == current_period).astype(int).values
                w = disdat[weights_name].values if weights_name else np.ones(len(disdat))
                n1 = np.sum(G + C)

                # Store current iteration data
                current_results = pd.DataFrame({
                    'Y': Y,
                    'post': post,
                    'G': G,
                    'group': g,
                    'time_period': current_period
                })
                results_df = pd.concat([results_df, current_results], ignore_index=True)

                # Check for sufficient observations
                skip_conditions = [
                    (np.sum(G * post) == 0, f"No units in group {g} in time period {current_period}"),
                    (np.sum(G * (1-post)) == 0, f"No units in group {g} in time period {tlist[pret]}"),
                    (np.sum(C * post) == 0, f"No control units for group {g} in time period {current_period}"),
                    (np.sum(C * (1-post)) == 0, f"No control units for group {g} in time period {tlist[pret]}")
                ]

                if any(cond[0] for cond in skip_conditions):
                    for condition, message in skip_conditions:
                        if condition:
                            warnings.warn(message)
                    continue

                # Get covariates
                try:
                    _, covariates = patsy.dmatrices(xformla, data=disdat)
                    covariates = np.asarray(covariates)
                except Exception as e:
                    try:
                        y_str, x_str = xformla.split("~")
                        xs1 = [x.strip() for x in x_str.split('+') if x.strip() != '1']
                        covariates = disdat[xs1].values
                        covariates = np.column_stack([covariates, np.ones(len(disdat))])
                    except:
                        covariates = np.ones((len(disdat), 1))

                # Select estimation method
                est_func = {
                    'reg': reg_did.reg_did_rc,
                    'ipw': ipwd_did.std_ipw_did_rc,
                    'dr': drdid.drdid_rc
                }.get(est_method, est_method if callable(est_method) else drdid.drdid_rc)

                # Compute treatment effects
                att_gt, att_inf_func = est_func(
                    y=Y, post=post, D=G, i_weights=w, covariates=covariates
                )

                # Handle NaN values
                if np.isnan(att_gt):
                    att_gt = np.nan
                    att_inf_func = np.zeros_like(att_inf_func)

                # Aggregate influence functions
                agg_inf = (pd.DataFrame({'inf': att_inf_func, 'id': right_ids})
                          .groupby('id')['inf'].sum())
                
                # Update influence function
                inf_zeros = np.zeros(n)
                unique_ids = data['rowid'].unique()
                dis_idx_inf = np.isin(unique_ids, agg_inf.index)
                inf_zeros[dis_idx_inf] = agg_inf.values

            # Update influence function matrix
            if panel:
                inffunc_updates.append({
                    'indices': treat_control_mask,
                    'values': att_inf_func
                })
            else:
                inffunc_updates.append({
                    'indices': dis_idx_inf,
                    'values': inf_zeros
                })

            counter += 1
            update_counter += 1

    # Update final influence function matrix
    for i, update in enumerate(inffunc_updates):
        idx = update['indices']
        inffunc[idx, i] = update['values'][idx]

    return results_df