import pickle
import pandas as pd
import numpy as np
from scipy.stats import norm


from csdid.aggte_fnc.utils import get_agg_inf_func, get_se, wif, AGGTEobj
from csdid.utils.mboot import mboot
import warnings
import numbers


# =============================================================================
# Scalar-argument validators (AGG-validate: O6-F5 / C1-6 / C1-7).
# Mirror R `did`'s utility_functions.R validators so the aggte stage refuses
# invalid alp / balance_e / min_e / max_e / na_rm cleanly instead of accepting
# nonsense (alp in {0,1,1.5,-0.1,NaN}, fractional/negative balance_e, truthy
# na_rm) or crashing with a cryptic error on string/vector args.
# =============================================================================
def _is_numeric_scalar(x):
    # Accept Python/NumPy real scalars (including whole-valued floats); reject
    # bool (R is.numeric(TRUE) is FALSE), strings, vectors/arrays, None.
    if isinstance(x, bool):
        return False
    if isinstance(x, numbers.Number):
        return True
    if isinstance(x, np.generic):
        return isinstance(x.item(), numbers.Number) and not isinstance(x.item(), bool)
    return False


def _validate_logical_scalar(x, name):
    # R: validate_logical_scalar -> single non-NA logical. Python bools (and
    # numpy bool_) qualify; truthy ints/strings do NOT.
    if isinstance(x, np.bool_):
        return
    if isinstance(x, bool):
        return
    raise ValueError(f"{name} must be a single logical (TRUE or FALSE).")


def _validate_numeric_scalar(x, name):
    # R: validate_numeric_scalar -> single non-missing number (Inf allowed).
    if (not _is_numeric_scalar(x)) or (isinstance(x, float) and np.isnan(x)) or \
       (isinstance(x, np.generic) and np.isnan(x)):
        raise ValueError(f"{name} must be a single non-missing number.")


def _validate_alp(alp, name="alp"):
    # R: validate_alp -> single number strictly between 0 and 1, non-NA.
    if (not _is_numeric_scalar(alp)) or \
       (isinstance(alp, (float, np.generic)) and np.isnan(alp)) or \
       alp <= 0 or alp >= 1:
        raise ValueError(f"{name} must be a single number strictly between 0 and 1.")


def _validate_nonnegative_whole_number(x, name):
    # R: validate_nonnegative_whole_number -> single non-NA finite non-negative
    # whole number (whole-valued floats accepted; Inf rejected).
    if (not _is_numeric_scalar(x)) or \
       (isinstance(x, (float, np.generic)) and np.isnan(x)) or \
       (not np.isfinite(x)) or x < 0 or x != round(x):
        raise ValueError(f"{name} must be a single non-negative whole number.")

def compute_aggte(MP,
                      typec         = "group",
                      balance_e     = None,
                      min_e         = float('-inf'),
                      max_e         = float('inf'),
                      na_rm         = False,
                      bstrap        = None,
                      biters        = None,
                      cband         = None,
                      alp           = None,
                      clustervars   = None,
                      call          = None):
  
# =============================================================================
# unpack MP object
# =============================================================================
    group       = np.array( MP['group'] )
    t           = np.array( MP['t']     )
    att         = np.array( MP['att'] )
    dp          = MP['DIDparams']
    tlist       = np.array( dp['tlist'] )
    glist       = np.array( dp['glist'] )
    data        = dp['data']
    inffunc     = MP['inffunc']['inffunc']
    n           = MP['n']
    gname       = dp['gname']
    tname       = dp['tname']
    idname      = dp['idname']
    # typec       = MP['type']  
    panel       = dp['panel']    


    
    if clustervars is None:
        clustervars = dp['clustervars']
    # AGG-clust-list (C3-F4): accept a list (faithful translation of R `c(...)`)
    # in addition to a bare string. att_gt's analytic SE path keys on a single
    # cluster variable; unwrap a length-1 list, take the first of a longer list
    # (R semantics: at most idname + one extra cluster var). Avoids the
    # `TypeError: unhashable type: 'list'` crash in get_se's column lookup.
    if isinstance(clustervars, (list, tuple, np.ndarray)):
        cl = list(clustervars)
        clustervars = cl[0] if len(cl) > 0 else None

    # AGG-clust-warn (O6-F2): if an extra cluster variable (beyond idname) is
    # requested but is not available in the retained att_gt data, warn and fall
    # back to non-clustered SE rather than silently doing so. Mirrors R
    # compute.aggte.R:87-107. (get_se's `clustervars in dta.columns` check is the
    # port's availability test; replicate it here so the warning surfaces.)
    if clustervars is not None and clustervars not in ("", idname):
        if clustervars not in data.columns:
            warnings.warn(
                "Clustered standard errors were requested in aggte() "
                f"(clustervars = '{clustervars}'), but the cluster information "
                "needed is not available in this att_gt object. Reporting "
                "standard errors that do NOT account for clustering; re-run "
                f"att_gt() with clustervars = '{clustervars}' for clustered "
                "inference."
            )
            clustervars = None

    if bstrap is None:
        bstrap = dp['bstrap']
    if biters is None:
        biters = dp['biters']
    if alp is None:
        alp = dp['alp']
    if cband is None:
        cband = dp['cband']



    # Overwrite MP objects (to compute bootstrap)
    MP['DIDparams']['clustervars'] = clustervars
    MP['DIDparams']['bstrap'] = bstrap
    MP['DIDparams']['biters'] = biters
    MP['DIDparams']['alp'] = alp
    MP['DIDparams']['cband'] = cband
    
# =============================================================================
#  Treat data
# =============================================================================

    if typec not in ["simple", "dynamic", "group", "calendar"]:
        raise ValueError("`typec` must be one of ['simple', 'dynamic', 'group', 'calendar']")

    # AGG-validate (O6-F5 / C1-6 / C1-7): mirror R compute.aggte.R:42-51,122
    # scalar-arg validation. alp is validated after defaulting (above).
    _validate_logical_scalar(na_rm, "na_rm")
    _validate_numeric_scalar(min_e, "min_e")
    _validate_numeric_scalar(max_e, "max_e")
    if balance_e is not None:
        _validate_nonnegative_whole_number(balance_e, "balance_e")
    _validate_alp(alp, "alp")

    # Removing missing values
    if na_rm:
        notna = ~np.isnan(att)
        group = group[notna]
        t = t[notna]
        att = att[notna]
        inffunc = inffunc[:, notna]
        glist = np.sort(np.unique(group))
    
        if typec == "group":
            gnotna = []
            for g in glist:
                # AGG-grp-maxe (O6-F4): restrict to the SAME max_e window the
                # group-specific estimate uses below; mirrors R
                # compute.aggte.R:169-181 `(t <= (group + max_e))`. Without it, a
                # group whose only non-NA post cell is PAST max_e survives this
                # filter but then has an empty in-window selection -> att.egt=NaN
                # and overall_att=NaN. No-op when max_e=Inf.
                indices = np.where((group == g) & (g <= t) & (t <= (g + max_e)))
                is_not_na = np.any(~np.isnan(att[indices]))
                gnotna.append(is_not_na)
            
            gnotna = np.array(gnotna)
            glist = glist[gnotna]
            not_all_na = np.isin(group, glist)
            group = group[not_all_na]
            t = t[not_all_na]
            att = att[not_all_na]
            inffunc = inffunc[:, not_all_na]
            glist = np.sort(np.unique(group))

        # All att_gt() estimates were NA -> na_rm stripped every cell, leaving
        # nothing to aggregate. Downstream `max(t)` would crash on the empty
        # sequence (V4-U1); R returns a clean refusal instead.
        if len(att) == 0:
            raise ValueError(
                "All att_gt() estimates are NA. Cannot compute aggregated "
                "treatment effects."
            )


    if (not na_rm) and np.any(np.isnan(att)):
        raise ValueError("Missing values at att_gt found. If you want to remove these, set `na_rm = True`.")

    if panel:
        dta = data[data[tname] == tlist[0]]
    else:
        dta_cols = [idname, gname, "w1"]
        dta = data.loc[:, dta_cols].groupby(idname, as_index=False).mean()
        dta = dta.drop(columns=[idname], errors="ignore")

# =============================================================================
#  Treat data 2
# =============================================================================

    
    originalt = t
    originalgroup = group
    originalglist = glist
    originaltlist = tlist
    # In case g's are not part of tlist
    originalgtlist = np.sort(np.unique(np.concatenate((originaltlist, originalglist))))
    uniquet = list(range(1, len(originalgtlist) + 1))
    
    # Function to switch from "new" t values to original t values
    def t2orig(t):
        return originalgtlist[uniquet.index(t) if t in uniquet else -1]
    
    # Function to switch between "original" t values and new t values
    def orig2t(orig):
        new_t = [uniquet[i] for i in range(len(originalgtlist)) if originalgtlist[i] == orig]
        out = new_t[0] if new_t else None
        return out
    
    t     = [orig2t(orig) for orig in originalt]
    group = [orig2t(orig) for orig in originalgroup]
    glist = [orig2t(orig) for orig in originalglist]
    tlist = np.asarray(list(set(t)))
    maxT  = max(t)
        
    # Set the weights
    # return data.columns
    weights_ind = dta['w1'].to_numpy()
    
    # We can work in overall probabilities because conditioning will cancel out
    # since it shows up in numerator and denominator
    pg = np.array([np.mean(weights_ind * (dta[gname].to_numpy() == g)) for g in originalglist])
    
    # Length of this is equal to the number of groups
    pgg = pg

    # Same but length is equal to the number of ATT(g,t)
    pg = [pg[glist.index(g)] for g in group]  
    
    # Which group time average treatment effects are post-treatment
    keepers = [i for i in range(len(group)) if group[i] <= t[i] <= (group[i] + max_e)] ### added second condition to allow for limit on longest period included in att
    
    # n x 1 vector of group variable
    G = [orig2t(g) for g in dta[gname].to_numpy()]

# =============================================================================
#  simple
# =============================================================================


    if typec == "simple":
        # AGG-empty-window (O6-F3 / C3-F6): an out-of-range max_e (or window)
        # can leave the post-treatment selection EMPTY. R refuses via
        # get_agg_inf_func()'s "No valid att_gt() estimates" stop; mirror that
        # cleanly here rather than returning a phantom att=None/se=None object.
        if len(keepers) == 0:
            raise ValueError(
                "No valid att_gt() estimates found for this aggregation. "
                "This may happen if all estimates for a particular group or "
                "time period are NA, or the requested max_e window excludes "
                "every post-treatment period."
            )
        # Simple ATT
        # Averages all post-treatment ATT(g,t) with weights given by group size
        pg = np.array(pg)
        simple_att = np.sum(att[keepers] * pg[keepers]) / np.sum(pg[keepers])
        if np.isnan(simple_att):
            simple_att = None
    
        # Get the part of the influence function coming from estimated weights
        simple_wif = wif(keepers, pg, weights_ind, G, group)
    
        # Get the overall influence function
        simple_if = get_agg_inf_func(att = att , 
                                     inffunc = inffunc , 
                                     whichones = keepers ,
                                     weights_agg = np.array(pg)[keepers]/np.sum(np.array(pg)[keepers]) , 
                                     wif = simple_wif )[:, None]
    
        # Get standard errors from the overall influence function
        simple_se = get_se(simple_if, dp)
        
        # AGG-sentinel (O5-F-A / C3-F7): use np.nan, not None, for the
        # degenerate-SE sentinel so downstream `overall_att + cval*overall_se`
        # in AGGTEobj degrades to NaN gracefully (matching dynamic/calendar)
        # instead of crashing with `float * None`. Mirrors R `simple.se <- NA`.
        if simple_se is not None:
            if simple_se <= np.sqrt(np.finfo(float).eps) * 10:
                simple_se = np.nan

        AGGTEobj_print = AGGTEobj(overall_att=simple_att, 
                                  overall_se=simple_se, 
                                  typec=typec,
                                  inf_function={'simple_att': simple_if}, 
                                  call=call, DIDparams=dp)
        
        return AGGTEobj_print


# =============================================================================
#  GRoup
# =============================================================================

    if typec == "group":
        group = np.array(group)
        t = np.array(t) 
        pg = np.array(pg) 
        selective_att_g = [np.mean(att[( group== g) & (t >= g) & (t <= (group + max_e))]) for g in glist]
        selective_att_g = np.asarray(selective_att_g)
        selective_att_g[np.isnan(selective_att_g)] = None
    
        selective_se_inner = [None] * len(glist)
        for i, g in enumerate(glist):
            whichg = np.where(np.logical_and.reduce((group == g, g <= t, t <= (group + max_e))))[0]
            weightsg =  pg[whichg] / np.sum(pg[whichg])
            inf_func_g = get_agg_inf_func(att = att , 
                                            inffunc = inffunc , 
                                            whichones = whichg ,
                                            weights_agg = weightsg , 
                                            wif = None)[:, None]
            se_g = get_se(inf_func_g, dp)
            selective_se_inner[i] = {'inf_func': inf_func_g, 'se': se_g}
            
        # recover standard errors separately by group   
        selective_se_g = np.asarray([item['se'] for item in selective_se_inner]).T
        
        selective_se_g[selective_se_g <= np.sqrt(np.finfo(float).eps) * 10] = None
        
        selective_inf_func_g = np.column_stack([elem["inf_func"] for elem in selective_se_inner])
  
        # use multiplier bootstrap (across groups) to get critical value
        # for constructing uniform confidence bands   
        selective_crit_val = norm.ppf(1 - alp/2)
        
        if dp['cband']:
            if not dp['bstrap']:
                # v7-CF2: warn (matches calendar branch + R `did`), not print().
                warnings.warn("Used bootstrap procedure to compute simultaneous confidence band")

            selective_crit_val = mboot(selective_inf_func_g, dp)['crit_val']

            if np.isnan(selective_crit_val) or np.isinf(selective_crit_val):
                warnings.warn("Simultaneous critical value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.")
                selective_crit_val = norm.ppf(1 - alp/2)
                dp['cband'] = False

            if selective_crit_val < norm.ppf(1 - alp/2):
                warnings.warn("Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals")
                selective_crit_val = norm.ppf(1 - alp/2)
                dp['cband'] = False

            if selective_crit_val >= 7:
                warnings.warn("Simultaneous critical value is arguably 'too large' to be reliable. This usually happens when the number of observations per group is small and/or there is not much variation in outcomes.")
  
        # get overall att under selective treatment timing
        # (here use pgg instead of pg because we can just look at each group)            
        selective_att = np.sum(selective_att_g * pgg) / np.sum(pgg)
        
        # account for having to estimate pgg in the influence function.
        # `keepers` indexes the DISTINCT groups (pgg/selective_att_g are per-group),
        # so `wif` must compare each unit's group against the distinct group values
        # (`glist`), NOT the per-(g,t) `group` vector whose first len(glist) entries
        # repeat the first group. Matches R `did`'s compute.aggte (passes the group
        # list). Fixes a systematic, anti-conservative error in the group-overall SE.
        selective_wif = wif(keepers = np.arange(1, len(glist)+1)-1,
                            pg  = pgg,
                            weights_ind = weights_ind,
                            G = G,
                            group = glist)
        
        # get overall influence function   
        selective_inf_func = get_agg_inf_func(att = selective_att_g, 
                                              inffunc = selective_inf_func_g,
                                              whichones = np.arange(1, len(glist)+1)-1, 
                                              weights_agg = pgg/np.sum(pgg),
                                              wif = selective_wif)[:, None]    
        
        # get overall standard error        
        selective_se = get_se(selective_inf_func, dp)
        # AGG-sentinel (O5-F-A / C3-F7): np.nan, not None, so AGGTEobj's
        # `overall_att + cval*overall_se` degrades to NaN instead of crashing
        # with `float * None`. Mirrors R `selective.se <- NA`.
        if not np.isnan(selective_se):
            if selective_se <= np.sqrt(np.finfo(float).eps) * 10:
                selective_se = np.nan
    
        AGGTEobj_print = AGGTEobj(overall_att = selective_att, 
                            overall_se = selective_se, 
                            typec = typec,
                            egt = originalglist,
                            att_egt = selective_att_g,
                            se_egt = selective_se_g,
                            crit_val_egt = selective_crit_val,
                            inf_function = {'selective_inf_func_g': selective_inf_func_g, 
                                            'selective_inf_func': selective_inf_func},
                            call = call, 
                            DIDparams = dp)

        return AGGTEobj_print


# =============================================================================
#  Dynamic
# =============================================================================

    if typec == "dynamic":
        # event times
        # this looks at all available event times
        # note: event times can be negative here.
        # note: event time = 0 corresponds to "on impact"
        eseq = np.unique(np.array(originalt) - np.array(originalgroup) ) # Subtract corresponding elements and convert to NumPy array
        eseq = np.sort(eseq)  # Sort the unique values in ascending order
    
        # if the user specifies balance_e, then we are going to
        # drop some event times and some groups; if not, we just
        # keep everything (that is what this variable is for)
        originalt = np.array(originalt)
        originalgroup = np.array(originalgroup)
        pg = np.array(pg)
        include_balanced_gt = np.repeat(True, len(originalgroup))

        if balance_e is not None:
            include_balanced_gt = (t2orig(maxT) - originalgroup >= balance_e)        
            eseq = np.unique(originalt[include_balanced_gt] - originalgroup[include_balanced_gt])
            eseq = np.sort(eseq)      
            eseq = eseq[(eseq <= balance_e) & (eseq >= balance_e - t2orig(maxT) + t2orig(1))]
        eseq = eseq[(eseq >= min_e) & (eseq <= max_e)]

        # AGG-empty-window (O6-F3 / C3-F6): a window (min_e>max_e, out-of-range
        # min_e/max_e, over-large balance_e) can leave NO event times. Without
        # this guard `np.column_stack([])` later raises a cryptic
        # "need at least one array to concatenate". Mirror R
        # compute.aggte.R:472-476 with a clean refusal.
        if len(eseq) == 0:
            raise ValueError(
                "No event times fall within the requested window. "
                "Adjust 'min_e'/'max_e' (and 'balance_e') so at least one event "
                "time is included."
            )

        dynamic_att_e = []
        for e in eseq:
            whiche = np.where((originalt - originalgroup == e) & include_balanced_gt)
            atte = att[whiche]
            pge = pg[whiche] / np.sum(pg[whiche])
            dynamic_att_e.append(np.sum(atte * pge))
        
        dynamic_se_inner = []
        for e in eseq:
            whiche = np.where((originalt - originalgroup == e) & (include_balanced_gt) )[0]
            pge = pg[whiche] / sum(pg[whiche])
            wif_e = wif(whiche, 
                        pg, 
                        weights_ind, 
                        G, 
                        group)
            inf_func_e = get_agg_inf_func(att         = att, 
                                          inffunc     = inffunc, 
                                          whichones   = whiche, 
                                          weights_agg = pge, 
                                          wif         = wif_e)[:, None]
            se_e = get_se(inf_func_e, dp)
            dynamic_se_inner.append({'inf_func': inf_func_e, 'se': se_e})

        dynamic_se_e = np.array([item['se'] for item in dynamic_se_inner]).T
        
        dynamic_se_e[dynamic_se_e <= np.sqrt(np.finfo(float).eps) * 10] = np.nan
        
        dynamic_inf_func_e = np.column_stack([item['inf_func'] for item in dynamic_se_inner])
                
        dynamic_crit_val = norm.ppf(1 - alp/2)
        if dp['cband']:
            if not dp['bstrap']:
                # v7-CF2: warn (matches calendar branch + R `did`), not print().
                warnings.warn('Used bootstrap procedure to compute simultaneous confidence band')
            dynamic_crit_val = mboot(dynamic_inf_func_e, dp)['crit_val']

            if np.isnan(dynamic_crit_val) or np.isinf(dynamic_crit_val):
                warnings.warn('Simultaneous critical value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.')
                dynamic_crit_val = norm.ppf(1 - alp/2)
                dp['cband'] = False

            if dynamic_crit_val < norm.ppf(1 - alp/2):
                warnings.warn('Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals')
                dynamic_crit_val = norm.ppf(1 - alp/2)
                dp['cband'] = False

            if dynamic_crit_val >= 7:
                warnings.warn("Simultaneous critical value is arguably 'too large' to be reliable. This usually happens when the number of observations per group is small and/or there is not much variation in outcomes.")

        epos = eseq >= 0
        dynamic_att = np.mean(np.array(dynamic_att_e)[epos])
        dynamic_inf_func = get_agg_inf_func(att         = np.array(dynamic_att_e)[epos],
                                            inffunc     = np.array(dynamic_inf_func_e[:, epos]),
                                            whichones   = np.arange(1, np.sum(epos)+1)-1,
                                            weights_agg = np.repeat(1 / np.sum(epos), np.sum(epos)),
                                            wif=None)[:, None]
        
        dynamic_se = get_se(dynamic_inf_func, dp)
        if not np.isnan(dynamic_se):
            if dynamic_se <= np.sqrt(np.finfo(float).eps) * 10:
                dynamic_se = np.nan

        AGGTEobj_print = AGGTEobj(overall_att=dynamic_att,
                                overall_se=dynamic_se,
                                typec=typec,
                                egt=eseq,
                                att_egt=dynamic_att_e,
                                se_egt=dynamic_se_e,
                                crit_val_egt=dynamic_crit_val,
                                inf_function={'dynamic_inf_func_e': dynamic_inf_func_e,
                                              'dynamic_inf_func': dynamic_inf_func},
                                call=call,
                                min_e=min_e,
                                max_e=max_e,
                                balance_e=balance_e,
                                DIDparams=dp)
            
        return AGGTEobj_print



# =============================================================================
#  Calendar
# =============================================================================

 # np.array(group)
    if typec == "calendar":
        # AGG-cal-warn (O6-F1 / C3-F5): min_e/max_e/balance_e have no effect on
        # calendar aggregation; warn so the (correct) unrestricted result is not
        # mistaken for a windowed one. Mirrors R compute.aggte.R:579-582.
        if np.isfinite(max_e) or np.isfinite(min_e) or (balance_e is not None):
            warnings.warn(
                "`min_e`, `max_e`, and `balance_e` are ignored for type = "
                "\"calendar\"; returning the unrestricted calendar-time effects."
            )
        minG = min(group)
        calendar_tlist = tlist[tlist >= minG]
        pg = np.array(pg)
        group = np.array(group)
        t = np.array(t)

        # AGG-cal-haspost (O6-F6): drop calendar periods with no non-missing
        # post-treatment ATT(g,t) cell (e.g. after na_rm stripped all of a
        # period's cells). Without this guard the port keeps a phantom att=0.0
        # period (empty `whicht` -> np.sum(pgt*attt)=0.0), corrupting the overall
        # calendar ATT. Mirrors R compute.aggte.R:593-600 `has_post`.
        has_post = np.array([np.any((t == t1) & (group <= t)) for t1 in calendar_tlist])
        calendar_tlist = calendar_tlist[has_post]
        if len(calendar_tlist) == 0:
            raise ValueError(
                "No calendar periods have non-missing post-treatment att_gt() "
                "estimates. Cannot compute calendar aggregation. Check your "
                "att_gt() results."
            )

        calendar_att_t = []
        for t1 in calendar_tlist:
            whicht = np.where((t == t1) & (group <= t))[0]
            attt = att[whicht]
            pgt = pg[whicht] / np.sum(pg[whicht])
            calendar_att_t.append(np.sum(pgt * attt))
            
        # get standard errors and influence functions
        # for each time specific att
        calendar_se_inner = []
        for t1 in calendar_tlist:
            which_t = np.where((t == t1) & (group <= t))[0]
            pgt = pg[which_t] / np.sum(pg[which_t])
            wif_t = wif(keepers=which_t, 
                        pg=pg, 
                        weights_ind=weights_ind, 
                        G=G, 
                        group=group)
            inf_func_t = get_agg_inf_func(att=att, 
                                            inffunc=inffunc, 
                                            whichones=which_t, 
                                            weights_agg=pgt, 
                                            wif=wif_t)[:, None]
            se_t = get_se(inf_func_t, dp)
            calendar_se_inner.append({"inf_func": inf_func_t, "se": se_t})
    
    
    
        # recover standard errors separately by time
        calendar_se_t = np.array([se["se"] for se in calendar_se_inner]).T
        calendar_se_t[calendar_se_t <= np.sqrt(np.finfo(float).eps) * 10] = np.nan
        
        # recover influence function separately by time
        calendar_inf_func_t = np.column_stack([se["inf_func"] for se in calendar_se_inner])
    
        # use multiplier boostrap (across groups) to get critical value
        # for constructing uniform confidence bands
        calendar_crit_val = norm.ppf(1 - alp/2)
        
        if dp['cband']:
            if not dp['bstrap']:
                warnings.warn('Used bootstrap procedure to compute simultaneous confidence band')
        
            # mboot function is not provided, please define it separately
            calendar_crit_val = mboot(calendar_inf_func_t, dp)['crit_val']
        
            if np.isnan(calendar_crit_val) or np.isinf(calendar_crit_val):
                warnings.warn('Simultaneous critical value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.')
                calendar_crit_val = norm.ppf(1 - alp/2)
                dp['cband'] = False
        
            if calendar_crit_val < norm.ppf(1 - alp/2):
                warnings.warn('Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals.')
                calendar_crit_val = norm.ppf(1 - alp/2)
                dp['cband'] = False
        
            if calendar_crit_val >= 7:
                warnings.warn("Simultaneous critical value is arguably 'too large' to be reliable. This usually happens when the number of observations per group is small and/or there is not much variation in outcomes.")
    
    
        # get overall att under calendar time effects
        # this is just average over all time periods
        calendar_att = np.mean(calendar_att_t)
        
        # get overall influence function
        calendar_inf_func = get_agg_inf_func(att=calendar_att_t,
                                             inffunc=calendar_inf_func_t,
                                             whichones=range(len(calendar_tlist)),
                                             weights_agg=np.repeat(1/len(calendar_tlist), len(calendar_tlist)),
                                             wif=None)[:, None]
        calendar_inf_func = np.array(calendar_inf_func)
        
        # get overall standard error
        calendar_se = get_se(calendar_inf_func, dp)
        if not np.isnan(calendar_se):
            if calendar_se <= np.sqrt(np.finfo(float).eps) * 10:
                calendar_se = np.nan
        
        AGGTEobj_print = AGGTEobj(overall_att=calendar_att,
                                overall_se=calendar_se,
                                typec=typec,
                                egt=list(map(t2orig, calendar_tlist)),
                                att_egt=calendar_att_t,
                                se_egt=calendar_se_t,
                                crit_val_egt=calendar_crit_val,
                                inf_function={"calendar_inf_func_t": calendar_inf_func_t,
                                              "calendar_inf_func": calendar_inf_func},
                                call=call,
                                DIDparams=dp)
    
        return AGGTEobj_print
