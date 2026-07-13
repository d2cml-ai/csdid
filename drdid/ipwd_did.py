from .utils import *

def std_ipw_did_panel(y1, y0, D, covariates, i_weights = None):
    # print("ipw: panel")
    D = np.asarray(D).flatten()
    n = len(D)
    delta_y = np.asarray(y1 - y0).flatten()
    int_cov = np.ones((n, 1))
    i_weights = np.asarray(i_weights).flatten()
    if covariates is not None:
        covariates = np.asarray(covariates)
        if np.all(covariates[:, 0] == 1):
            int_cov = covariates
        else:
            int_cov = np.column_stack((np.ones(n), covariates))
    
    if i_weights is None:
        i_weights = np.ones(n)
    elif np.min(i_weights) < 0:
        raise ValueError("i_weights must be non-negative")
    
    i_weights = i_weights / np.mean(i_weights)
    pscore_model = sm.GLM(D, int_cov, family=sm.families.Binomial(), freq_weights=i_weights)
    pscore_results = pscore_model.fit()
    # print(D.mean())
    # print(pscore_results.summary2())
    if not pscore_results.converged:
        print("Warning: glm algorithm did not converge")
    if np.any(np.isnan(pscore_results.params)):
        raise ValueError("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
    ps_fit = pscore_results.predict()
    ps_fit = np.minimum(ps_fit, 1 - 1e-16)

    w_treat = i_weights * D
    w_cont = i_weights * ps_fit * (1 - D) / (1 - ps_fit)
    
    att_treat = w_treat * delta_y
    att_cont = w_cont * delta_y

    eta_treat = mean(att_treat) / mean(w_treat)
    eta_cont = mean(att_cont) / mean(w_cont)

    ipw_att = eta_treat - eta_cont

    score_ps = i_weights[:, np.newaxis] * (D - ps_fit)[:, np.newaxis] * int_cov
    Hessian_ps = pscore_results.cov_params() * n
    asy_lin_rep_ps = np.dot(score_ps, Hessian_ps)

    inf_treat = (att_treat - w_treat * eta_treat) / mean(w_treat)
    inf_cont_1 = att_cont - w_cont * eta_cont
    pre_m2 = w_cont * (delta_y - eta_cont)
    M2 = np.mean(pre_m2[:, np.newaxis] * int_cov, axis = 0)
    # print(M2)
    inf_cont_2 = np.dot(asy_lin_rep_ps, M2)

    inf_control = (inf_cont_1 + inf_cont_2) / np.mean(w_cont)
    att_inf_func = inf_treat - inf_control
    se = np.std(att_inf_func, ddof=1) / np.sqrt(n)
    # print(f"att: {ipw_att} \t se: {se}")    
    # print(np.std(att_inf_func) / np.sqrt(n))
    return ipw_att, att_inf_func

def std_ipw_did_rc(y, post, D, covariates, i_weights = None):
    # print("ipw: rc")
    D = np.asarray(D).flatten()
    y = np.asarray(y).flatten()
    post = np.asarray(post).flatten()
    n = len(D)
    if covariates is None:
        int_cov = np.ones((n, 1))
    else:
        covariates = np.asarray(covariates)
        if np.all(covariates[:, 0] == 1):
            int_cov = covariates
        else:
            int_cov = np.column_stack((np.ones(n), covariates))
    
    # Pesos
    if i_weights is None:
        i_weights = np.ones(n)
    else:
        i_weights = np.asarray(i_weights)
        if np.min(i_weights) < 0:
            raise ValueError("i_weights must be non-negative")
    
    # Normalizar pesos
    i_weights = np.asarray(i_weights).flatten()
    i_weights = i_weights / np.mean(i_weights)

    pscore_model = sm.GLM(D, int_cov, family=sm.families.Binomial(), freq_weights=i_weights)
    pscore_results = pscore_model.fit()
    if not pscore_results.converged:
        print("Warning: glm algorithm did not converge")
    if np.any(np.isnan(pscore_results.params)):
        raise ValueError("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
    ps_fit = pscore_results.predict()
    ps_fit = np.minimum(ps_fit, 1 - 1e-16)


    w_treat_pre = i_weights * D * (1 - post)
    w_treat_post = i_weights * D * post
    # print(np.mean(w_treat_pre))

    w_cont_pre = i_weights * ps_fit * (1 - D) * (1 - post)/(1 - ps_fit)
    w_cont_post = i_weights * ps_fit * (1 - D) * post/(1 - ps_fit)

    # Elements of the influence function (summands)
    eta_treat_pre = w_treat_pre * y / np.mean(w_treat_pre)
    eta_treat_post = w_treat_post * y / np.mean(w_treat_post)
    # print(eta_treat_pre)

    eta_cont_pre = w_cont_pre * y / np.mean(w_cont_pre)
    eta_cont_post = w_cont_post * y / np.mean(w_cont_post)

    # Estimator of each component
    att_treat_pre = np.mean(eta_treat_pre)
    att_treat_post = np.mean(eta_treat_post)
    att_cont_pre = np.mean(eta_cont_pre)
    att_cont_post = np.mean(eta_cont_post)
    ipw_att = (att_treat_post - att_treat_pre) - (att_cont_post - att_cont_pre)

    score_ps = (i_weights * (D - ps_fit))[:, np.newaxis] * int_cov
    Hessian_ps = pscore_results.cov_params() * n
    asy_lyn_rep_ps = np.dot(score_ps, Hessian_ps)

    inf_treat_pre = eta_treat_pre - w_treat_pre * att_treat_pre/np.mean(w_treat_pre)
    inf_treat_post = eta_treat_post - w_treat_post * att_treat_post/np.mean(w_treat_post)
    inf_treat = inf_treat_post - inf_treat_pre
    # Now, get the influence function of control component
    # Leading term of the influence function: no estimation effect
    inf_cont_pre = eta_cont_pre - w_cont_pre * att_cont_pre/np.mean(w_cont_pre)
    inf_cont_post = eta_cont_post - w_cont_post * att_cont_post/np.mean(w_cont_post)
    inf_cont = inf_cont_post - inf_cont_pre

    # Estimation effect from gamma hat (pscore)
    # Derivative matrix (k x 1 vector)
  
    M2_pre = np.mean((w_cont_pre *(y - att_cont_pre))[:, np.newaxis] * int_cov, axis = 0)/np.mean(w_cont_pre)
    M2_post = np.mean((w_cont_post *(y - att_cont_post))[:, np.newaxis] * int_cov, axis = 0)/np.mean(w_cont_post)

    # Now the influence function related to estimation effect of pscores
    M2 = M2_post - M2_pre
    # print()

    inf_cont_ps = np.dot(asy_lyn_rep_ps, M2)

    # Influence function for the control component
    inf_cont = inf_cont + inf_cont_ps

    #get the influence function of the DR estimator (put all pieces together)
    att_inf_func = inf_treat - inf_cont
    # print(np.std(att_inf_func) / np.sqrt(n))
    se = np.std(att_inf_func, ddof=1) / np.sqrt(n)
    # print(f"att: {ipw_att} \t se: {se}")    
    return ipw_att, att_inf_func





