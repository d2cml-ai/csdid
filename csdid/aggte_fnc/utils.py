from csdid.utils.bmisc import TorF
from csdid.utils.mboot import mboot

import numpy as np
import scipy.stats as stats
import pandas as pd

def wif(keepers, pg, weights_ind, G, group):
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
        return bout['se']
    else:
        return np.sqrt(np.mean((thisinffunc)**2) / n)

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
    overall_cband_upper = out["overall_att"] + pointwise_cval * out["overall_se"]
    overall_cband_lower = out["overall_att"] - pointwise_cval * out["overall_se"]
    out1 = np.column_stack((out["overall_att"], out["overall_se"], overall_cband_lower, overall_cband_upper))
    out1 = np.round(out1, 4)
    overall_sig = (overall_cband_upper < 0) | (overall_cband_lower > 0)
    overall_sig[np.isnan(overall_sig)] = False
    overall_sig_text = np.where(overall_sig, "*", "")
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


