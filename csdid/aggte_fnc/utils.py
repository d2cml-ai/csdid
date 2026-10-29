from csdid.utils.bmisc import TorF
from csdid.utils.mboot import mboot
from csdid.utils.weights import _weight_mean

import numpy as np
import scipy.stats as stats
import pandas as pd

def wif(keepers, pg, weights_ind, G, group):
    """Influence of estimated cohort weights, reusing the numerator matrix."""
    pg_array = np.asarray(pg)
    weights = np.asarray(weights_ind)
    unit_groups = np.asarray(G)
    groups = np.asarray(group)
    selected = np.asarray(keepers)
    # Preserve the scalar casting, errors and warnings of unusual direct calls.
    # Normal aggregation inputs are nonempty one-dimensional float64 arrays.
    if (type(weights_ind) is not np.ndarray or type(G) not in (list, tuple, np.ndarray) or
            selected.ndim != 1 or selected.size == 0 or
            selected.dtype.kind not in "iu" or pg_array.ndim != 1 or
            groups.ndim != 1 or groups.dtype.kind not in "iuf" or
            weights.ndim != 1 or unit_groups.ndim != 1 or
            weights.shape != unit_groups.shape or
            pg_array.dtype != np.float64 or weights.dtype != np.float64):
        return _wif_legacy(keepers, pg, weights_ind, G, group)
    selected_pg = pg_array[selected]
    if (not np.all(np.isfinite(selected_pg)) or
            not np.all(np.isfinite(weights))):
        return _wif_legacy(keepers, pg, weights_ind, G, group)
    try:
        # If arithmetic would warn, retain the original per-column warnings.
        with np.errstate(all="raise"):
            denominator = np.sum(selected_pg)
            if denominator == 0:
                raise FloatingPointError
            numerator = (weights[:, None] *
                         (unit_groups[:, None] == groups[selected]) - selected_pg)
            row_sum = np.sum(numerator, axis=1)
            # Keep builtin sum: the original reduction is sequential.
            squared_denominator = sum(selected_pg) ** 2
            if (not np.isfinite(squared_denominator) or
                    squared_denominator < np.finfo(float).tiny):
                raise FloatingPointError
            multiplier = selected_pg / squared_denominator
            numerator /= denominator
            numerator -= np.outer(row_sum, multiplier)
    except FloatingPointError:
        stable = _wif_scaled(selected, selected_pg, weights, unit_groups, groups)
        if stable is not None:
            return stable
        return _wif_legacy(keepers, pg, weights_ind, G, group)
    return numerator


def _valid_cohort_labels(values):
    if values.dtype.kind in "iuf":
        return np.all(np.isfinite(values))
    if values.dtype.kind == "O":
        # orig2t maps never-treated units to None. Those units legitimately
        # match none of the selected treated cohorts.
        return all(value is None or
                   (isinstance(value, (int, float, np.integer, np.floating))
                    and np.isfinite(value)) for value in values)
    return False


def _wif_scaled(selected, selected_pg, weights, unit_groups, groups):
    """Equivalent selected-cohort algebra without squaring its total mass."""
    # Only use the repair for actual cohort means. Unusual direct calls retain
    # the historical casting, errors and warnings of the legacy helper.
    if (not _valid_cohort_labels(unit_groups) or
            not np.all(np.isfinite(groups[selected])) or
            np.any(weights < 0) or np.any(selected_pg < 0)):
        return None
    precision = max([1e-12] + [8 * np.finfo(values.dtype).eps
                               for values in (selected_pg, weights)
                               if values.dtype.kind == "f"])
    # Promotion belongs only to the repair. Ordinary float32/pandas calls keep
    # their historical casting and reductions in the legacy implementation.
    original_pg = selected_pg
    selected_pg = np.asarray(selected_pg, dtype=float)
    weights = np.asarray(weights, dtype=float)
    matches = unit_groups[:, None] == groups[selected]
    empirical_pg = np.array([_weight_mean(weights * matches[:, i])
                             for i in range(len(selected))])
    if (original_pg.dtype.kind == "f" and original_pg.dtype.itemsize < 8
            and np.any((empirical_pg > 0)
                       & (empirical_pg < np.finfo(original_pg.dtype).tiny))):
        # A caller's native float32 cohort mean may already have rounded to
        # zero. Repair only masses consistent with rounding the actual means
        # back to that dtype, leaving malformed direct-call diagnostics intact.
        with np.errstate(under="ignore"):
            rounded_pg = empirical_pg.astype(original_pg.dtype).astype(float)
        if np.allclose(selected_pg, rounded_pg, rtol=precision, atol=0):
            selected_pg = empirical_pg
    if not np.max(selected_pg) > 0:
        return None
    if not np.allclose(selected_pg, empirical_pg, rtol=precision, atol=0):
        return None
    with np.errstate(all="ignore"):
        scale = np.max(selected_pg)
        scaled_pg = selected_pg / scale
        # Mask before dividing: weights in excluded cohorts can have a very
        # different scale and have no contribution to this conditional average.
        numerator = (weights[:, None] * matches) / scale - scaled_pg
        numerator /= np.sum(scaled_pg)
        result = numerator - np.outer(np.sum(numerator, axis=1),
                                      scaled_pg / sum(scaled_pg))
    return result if np.all(np.isfinite(result)) else None


def _wif_legacy(keepers, pg, weights_ind, G, group):
    # Valid float32 and pandas inputs also need the numerical repair, but retain
    # the original implementation whenever its denominator arithmetic is safe.
    pg_array, weights = np.asarray(pg), np.asarray(weights_ind)
    unit_groups, groups = np.asarray(G), np.asarray(group)
    selected = np.asarray(keepers)
    if (isinstance(weights_ind, (np.ndarray, pd.Series)) and
            pg_array.ndim == weights.ndim == unit_groups.ndim == groups.ndim == 1 and
            selected.ndim == 1 and selected.size and selected.dtype.kind in "iu" and
            pg_array.dtype.kind == "f" and weights.dtype.kind in "fiu" and
            _valid_cohort_labels(unit_groups) and groups.dtype.kind in "iuf" and
            weights.shape == unit_groups.shape and
            np.all((selected >= 0) & (selected < len(pg_array))) and
            np.all(selected < len(groups))):
        selected_pg = pg_array[selected]
        with np.errstate(all="ignore"):
            total = sum(selected_pg)
            squared = total ** 2
            multiplier = selected_pg / squared
        tiny = np.finfo(np.asarray(squared).dtype).tiny
        low_precision_mass = (
            weights.dtype.kind == "f" and weights.dtype.itemsize < 8
            and np.any((selected_pg > 0)
                       & (selected_pg < np.finfo(weights.dtype).tiny)))
        if (not np.isfinite(squared) or squared < tiny or low_precision_mass or
                not np.all(np.isfinite(multiplier))):
            stable = _wif_scaled(selected, selected_pg, weights, unit_groups, groups)
            if stable is not None:
                return stable
    return _wif_legacy_unscaled(keepers, pg, weights_ind, G, group)


def _wif_legacy_unscaled(keepers, pg, weights_ind, G, group):
    # note: weights are all of the form P(G=g|cond)/sum_cond(P(G=g|cond))
    # this is equal to P(G=g)/sum_cond(P(G=g)) which simplifies things here
    pg = np.array(pg)
    group = np.array(group)
    
    # effect of estimating weights in the numerator
    if1 = np.empty((len(weights_ind), len(keepers)))
    for i, k  in enumerate(keepers):
        numerator = (weights_ind * 1 * TorF(G == group[k])) - pg[k]
        # denominator = sum(np.array(pg)[keepers]) )[:, None]  
        denominator = np.sum(pg[keepers])

        result = numerator[:, None]  / denominator
        if1[:, i] = result.squeeze()
    
    # effect of estimating weights in the denominator
    if2 = np.empty((len(weights_ind), len(keepers)))
    for i, k  in enumerate(keepers):
        numerator = ( weights_ind * 1 * TorF(G == group[k]) ) - pg[k]
        # result = numerator.to_numpy()[:, None]  @ multipler[:, None].T
        if2[:, i] = numerator.squeeze()
    if2 = np.sum(if2, axis=1)    
    multiplier = ( pg[keepers] / sum( pg[keepers] ) ** 2 )   
    if2 = np.outer( if2 , multiplier)

    # if1 = [((weights_ind * 1*TorF(G==group[k])) - pg[k]) / sum(pg[keepers]) for k in keepers]
    # if2 = np.dot(np.array([weights_ind*1*TorF(G==group[k]) - pg[k] for k in keepers]).T, pg[keepers]/(sum(pg[keepers])**2))
    wif_factor = if1 - if2
    return wif_factor

def get_agg_inf_func(att, inffunc, whichones, weights_agg, wif=None):
    # enforce weights are in matrix form
    weights_agg = np.asarray(weights_agg)

    # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
    thisinffunc = np.dot(inffunc[:, whichones], weights_agg)

    # Incorporate influence function of the weights
    if wif is not None:
        thisinffunc = thisinffunc + np.dot(wif, np.array(att[whichones]))
        
    # return influence function
    return thisinffunc


def _positive_finite_se(se):
    """Only positive finite scales support inference, in any outcome units."""
    values = np.asarray(se)
    result = np.where(np.isfinite(values) & (values > 0), values, np.nan)
    return result[()] if result.ndim == 0 else result


def get_se(thisinffunc, DIDparams=None):
    alpha = 0.05
    bstrap = False
    if DIDparams is not None:
        bstrap = DIDparams['bstrap']
        alpha = DIDparams['alp']
        cband = DIDparams['cband']
        n = len(thisinffunc)

    if bstrap:
        bout = mboot(thisinffunc, DIDparams)
        return _positive_finite_se(bout['se'])
    else:
        # Analytical standard errors
        thisinffunc = np.asarray(thisinffunc).flatten()
        n = len(thisinffunc)
        clustervars = DIDparams.get('clustervars') if DIDparams else None
        
        if clustervars is not None and DIDparams is not None:
            # Cluster-robust analytical SEs:
            # SE = sqrt(Var(cluster_sums) / n_clusters^2)
            data = DIDparams['data']
            idname = DIDparams['idname']
            tname = DIDparams['tname']
            panel = DIDparams['panel']
            
            # Get unit-to-cluster mapping
            if panel:
                tlist = np.sort(data[tname].unique())
                dta = data[data[tname] == tlist[0]]
            else:
                dta = data.drop_duplicates(subset=[idname])
            
            if clustervars in dta.columns:
                cluster_map = dta[[idname, clustervars]].drop_duplicates().set_index(idname)[clustervars]
                unit_ids = dta[idname].unique()
                cluster_labels = cluster_map.reindex(unit_ids).values
                
                # Sum IFs by cluster
                unique_clusters = np.unique(cluster_labels)
                n_clusters = len(unique_clusters)
                cluster_sums = np.zeros(n_clusters)
                for i, c in enumerate(unique_clusters):
                    mask = cluster_labels == c
                    cluster_sums[i] = np.sum(thisinffunc[mask])
                
                # Cluster-robust variance: SE = sqrt(sum(S_c^2)) / n
                # where n is the total number of units (not clusters),
                # matching R did v2.5.1
                n_units = len(thisinffunc)
                return _positive_finite_se(np.sqrt(np.sum(cluster_sums**2)) / n_units)
            
        # Standard (iid) analytical SE
        return _positive_finite_se(np.sqrt(np.mean(thisinffunc**2) / n))


_ORIGINAL_GET_SE = get_se


def _make_cached_get_se(get_se_func):
    """Cache analytical cluster membership for one compute_aggte invocation.

    Bootstrap calls retain their original dispatch and random draw sequence.
    An overridden public helper is also called unchanged.
    """
    if get_se_func is not _ORIGINAL_GET_SE:
        return get_se_func
    cache = {}

    def cached_get_se(thisinffunc, DIDparams=None):
        if DIDparams is None or DIDparams['bstrap']:
            return get_se_func(thisinffunc, DIDparams)
        clustervars = DIDparams.get('clustervars')
        if clustervars is None:
            return get_se_func(thisinffunc, DIDparams)
        data = DIDparams['data']
        idname, tname = DIDparams['idname'], DIDparams['tname']
        key = (id(data), idname, tname, DIDparams['panel'], clustervars)
        if key not in cache:
            if DIDparams['panel']:
                tlist = np.sort(data[tname].unique())
                dta = data[data[tname] == tlist[0]]
            else:
                dta = data.drop_duplicates(subset=[idname])
            if clustervars not in dta.columns:
                cache[key] = None
            else:
                cluster_map = (dta[[idname, clustervars]].drop_duplicates()
                               .set_index(idname)[clustervars])
                labels = cluster_map.reindex(dta[idname].unique()).values
                # Ascending indices preserve the original boolean-mask sum order.
                rows = [np.flatnonzero(labels == c) for c in np.unique(labels)]
                cache[key] = (len(labels), rows)
        metadata = cache[key]
        influence = np.asarray(thisinffunc).flatten()
        if metadata is None or metadata[0] != len(influence):
            return get_se_func(thisinffunc, DIDparams)
        rows = metadata[1]
        # Match the original float64 accumulator even for float32 direct calls.
        cluster_sums = np.zeros(len(rows))
        for i, index in enumerate(rows):
            cluster_sums[i] = np.sum(influence[index])
        return _positive_finite_se(np.sqrt(np.sum(cluster_sums**2)) / len(influence))

    return cached_get_se

def AGGTEobj(overall_att=None,
             overall_se=None,
             typec="simple",
             egt=None,
             att_egt=None,
             se_egt=None,
             crit_val_egt=None,
             inf_function=None,
             min_e=None,
             max_e=None,
             balance_e=None,
             call=None,
             DIDparams=None):

    out = {
        "overall_att": overall_att,
        "overall_se": overall_se,
        "type": typec,
        "egt": egt,
        "att_egt": att_egt,
        "se_egt": se_egt,
        "crit_val_egt": crit_val_egt,
        "inf_function": inf_function,
        "min_e": min_e,
        "max_e": max_e,
        "balance_e": balance_e,
        "call": call,
        "DIDparams": DIDparams
    }
    
    
    # Overall estimates
    alp = out["DIDparams"]["alp"]
    pointwise_cval = stats.norm.ppf(1 - alp / 2)
    # AGG-sentinel (O5-F-A / C3-F7): a degenerate-variance or all-NA aggregation
    # can leave overall_att / overall_se as None. `float * None` would crash with
    # a TypeError; coerce a None/non-finite SE (and a None att) to np.nan so the
    # confidence band degrades to NaN gracefully, like R returns NA, for ALL four
    # aggregation types.
    _oa = out["overall_att"]
    _ose = out["overall_se"]
    _oa = np.nan if _oa is None else _oa
    _ose = np.nan if _ose is None else _ose
    overall_cband_upper = _oa + pointwise_cval * _ose
    overall_cband_lower = _oa - pointwise_cval * _ose
    out1 = np.column_stack((_oa, _ose, overall_cband_lower, overall_cband_upper))
    out1 = np.round(out1, 4)
    overall_sig = (overall_cband_upper < 0) | (overall_cband_lower > 0)
    overall_sig = np.asarray(overall_sig)
    overall_sig = np.where(np.isnan(overall_sig), False, overall_sig)
    overall_sig_text = np.atleast_1d(np.where(overall_sig, "*", ""))
    out1 = np.column_stack((out1, overall_sig_text))
    
    print("\n")
    if out["type"] == "dynamic":
        print("Overall summary of ATT's based on event-study/dynamic aggregation:")
    elif out["type"] == "group":
        print("Overall summary of ATT's based on group/cohort aggregation:")
    elif out["type"] == "calendar":
        print("Overall summary of ATT's based on calendar time aggregation:")
    colnames = ["ATT", "Std. Error", f"[{100 * (1 - out['DIDparams']['alp'])}%"," Conf. Int.]", ""]
    print(pd.DataFrame(out1, columns=colnames).to_string(index=False))
    print("\n")
    
    # Handle cases depending on type
    if out["type"] in ["group", "dynamic", "calendar"]:
        if out["type"] == "dynamic":
            c1name = "Event time"
            print("Dynamic Effects:")
        elif out["type"] == "group":
            c1name = "Group"
            print("Group Effects:")
        elif out["type"] == "calendar":
            c1name = "Time"
            print("Time Effects (calendar):")
    
        cband_text1a = f"{100 * (1 - out['DIDparams']['alp'])}% "
        cband_text1b = "Simult. " if out["DIDparams"]["bstrap"] else "Pointwise "
        cband_text1 = f"[{cband_text1a}{cband_text1b}"
    
        cband_lower = out["att_egt"] - out["crit_val_egt"] * out["se_egt"]
        cband_upper = out["att_egt"] + out["crit_val_egt"] * out["se_egt"]
    
        sig = (cband_upper < 0) | (cband_lower > 0)
        sig[np.isnan(sig)] = False
        sig_text = np.where(sig, "*", "").T

        out2 = pd.DataFrame([out["egt"], 
                             out["att_egt"], 
                             out["se_egt"].flatten(), 
                             np.hstack(cband_lower),
                             np.hstack(cband_upper)]).T
        
        out2 = out2.round(4)
        out2[0] = out2[0].astype(int)
        out2 = pd.concat([out2, pd.DataFrame(sig_text, columns=['sig_text']) ], axis=1)    
        
        out2.columns = [c1name, "Estimate", "Std. Error", cband_text1, "Conf. Band", ""]
        print(out2)
    
    
    

    
    print("---")
    print("Signif. codes: `*' confidence band does not cover 0")
    
    # Set control group text
    control_group = out["DIDparams"]["control_group"]
    control_group_text = None
    if control_group == "nevertreated":
        control_group_text = "Never Treated"
    elif control_group == "notyettreated":
        control_group_text = "Not Yet Treated"
    
    if control_group:
        print("Control Group: ", control_group_text, ", ")
    
    # Anticipation periods
    print("Anticipation Periods: ", out["DIDparams"]["anticipation"])
    
    # Estimation method text
    est_method = out["DIDparams"]["est_method"]
    if isinstance(est_method, str):
        est_method_text = est_method
        if est_method == "dr":
            est_method_text = "Doubly Robust"
        elif est_method == "ipw":
            est_method_text = "Inverse Probability Weighting"
        elif est_method == "reg":
            est_method_text = "Outcome Regression"
    
        print("Estimation Method: ", est_method_text)
        print("\n")
        
    return out

