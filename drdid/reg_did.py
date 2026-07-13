import numpy as np
import statsmodels.api as sm
from numpy import ndarray

from .utils import * 
import numpy as np
import statsmodels.api as sm
from numpy import ndarray

from .utils import * 

def asy_lin_wols(d, post, y, out_y, int_cov, i_w, n):
  weigths_ols = i_w * d * post
  # weigths_ols_pre
  wols_x = weigths_ols[:, n_x] * int_cov
  wols_ex = (weigths_ols * (y - out_y))[:, n_x] * int_cov
  cr = np.dot(wols_x.T, int_cov) / n
  xpx_inv = qr_solver(cr)
  asy_lin_rep_ols = np.dot(wols_ex, xpx_inv)
  return asy_lin_rep_ols

import numpy as np
import statsmodels.api as sm

def reg_did_panel(y1, y0, D, covariates=None, i_weights=None):
    # Convert inputs to numpy arrays
    # print("reg: panel")
    y1 = np.array(y1)
    y0 = np.array(y0)
    D = np.array(D)
    
    # Sample size
    n = len(D)
    
    # Generate deltaY
    deltaY = y1 - y0
    
    # Add constant to covariate vector
    if covariates is None:
        int_cov = np.ones((n, 1))
    else:
        covariates = np.asarray(covariates)
        if np.all(covariates[:, 0] == 1):
            int_cov = covariates
        else:
            int_cov = np.column_stack((np.ones(n), covariates))
    
    
    # Weights
    if i_weights is None:
        i_weights = np.ones(n)
    else:
        i_weights = np.array(i_weights)
        if np.min(i_weights) < 0:
            raise ValueError("i_weights must be non-negative")
    
    # Normalize weights
    i_weights = i_weights / np.mean(i_weights)
    
    # Compute the Outcome regression for the control group using OLS
    model = sm.WLS(deltaY[D == 0], int_cov[D == 0], weights=i_weights[D == 0])
    results = model.fit()
    reg_coeff = results.params
    # print(reg_coeff)
    if np.any(np.isnan(reg_coeff)):
        raise ValueError("Outcome regression model coefficients have NA components. Multicollinearity (or lack of variation) of covariates is probably the reason for it.")
    
    out_delta = np.dot(int_cov, reg_coeff)
    # print(out_delta)
    
    # Compute the OR-DiD estimator
    w_treat = i_weights * D
    w_cont = i_weights * D
    reg_att_treat = w_treat * deltaY
    reg_att_cont = w_cont * out_delta
    eta_treat = np.mean(reg_att_treat) / np.mean(w_treat)
    eta_cont = np.mean(reg_att_cont) / np.mean(w_cont)
    # print(np.mean(reg_att_cont), np.mean(w_cont))
    reg_att = eta_treat - eta_cont
    
    # Compute influence function
    weights_ols = i_weights * (1 - D)
    wols_x = weights_ols[:, np.newaxis] * int_cov
    wols_eX = weights_ols[:, np.newaxis] * (deltaY - out_delta)[:, np.newaxis] * int_cov
    XpX_inv = np.linalg.inv(np.dot(wols_x.T, int_cov) / n)
    asy_lin_rep_ols = np.dot(wols_eX, XpX_inv)
    
    inf_treat = (reg_att_treat - w_treat * eta_treat) / np.mean(w_treat)
    
    inf_cont_1 = reg_att_cont - w_cont * eta_cont
    M1 = np.mean(w_cont[:, np.newaxis] * int_cov, axis=0)
    inf_cont_2 = np.dot(asy_lin_rep_ols, M1)
    inf_control = (inf_cont_1 + inf_cont_2) / np.mean(w_cont)
    
    reg_att_inf_func = inf_treat - inf_control
    se = np.std(reg_att_inf_func, ddof=1) / np.sqrt(n)
    # print(f"att: {reg_att} \t se: {se}")
    
    return reg_att, reg_att_inf_func




def reg_did_rc(y, post, D, covariates, i_weights=None):
#   print("reg_rc")
  D = np.asarray(D).flatten()
  post = np.asarray(post).flatten()
  n = len(D)
  y = np.asarray(y).flatten()
  i_weights = np.asarray(i_weights).flatten()
  int_cov = np.ones((n, 1))
  
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
  
  # Pre-treatment regression
  mask_pre = (D == 0) & (post == 0)
  X_pre = int_cov[mask_pre]
  y_pre = y[mask_pre]
  w_pre = i_weights[mask_pre]
  model_pre = sm.WLS(y_pre, X_pre, weights=w_pre)
  results_pre = model_pre.fit()
  reg_coeff_pre = results_pre.params
  
  if np.any(np.isnan(reg_coeff_pre)):
      raise ValueError("Outcome regression model coefficients have NA components. \n Multicollinearity of covariates is probably the reason for it.")
  
  out_y_pre = np.dot(int_cov, reg_coeff_pre)
  
  # Post-treatment regression
  mask_post = (D == 0) & (post == 1)
  X_post = int_cov[mask_post]
  y_post = y[mask_post]
  w_post = i_weights[mask_post]
  model_post = sm.WLS(y_post, X_post, weights=w_post)
  results_post = model_post.fit()
  reg_coeff_post = results_post.params
  
  if np.any(np.isnan(reg_coeff_post)):
      raise ValueError("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is probably the reason for it.")
  
  out_y_post = np.dot(int_cov, reg_coeff_post)
  
  w_treat_pre = i_weights * D * (1 - post)
  w_treat_post = i_weights * D * post
  w_cont = i_weights * D
  reg_att_treat_pre = w_treat_pre * y
  reg_att_treat_post = w_treat_post * y
  reg_att_cont = w_cont * (out_y_post - out_y_pre)
  eta_treat_pre = np.mean(reg_att_treat_pre) / np.mean(w_treat_pre)
  eta_treat_post = np.mean(reg_att_treat_post) / np.mean(w_treat_post)
  eta_cont = np.mean(reg_att_cont) / np.mean(w_cont)
  reg_att = (eta_treat_post - eta_treat_pre) - eta_cont
  
  weights_ols_pre = i_weights * (1 - D) * (1 - post)
  # print(weights_ols_pre.reshape((n, 1)))
  wols_x_pre = weights_ols_pre[:, np.newaxis] * int_cov
  wols_eX_pre = weights_ols_pre[:, np.newaxis] * (y - out_y_pre)[:, np.newaxis] * int_cov
  XpX_inv_pre = np.linalg.inv(np.dot(wols_x_pre.T, int_cov) / n)
  asy_lin_rep_ols_pre = np.dot(wols_eX_pre, XpX_inv_pre)
  
  weights_ols_post = i_weights * (1 - D) * post
  wols_x_post = weights_ols_post[:, np.newaxis] * int_cov
  wols_eX_post = weights_ols_post[:, np.newaxis] * (y - out_y_post)[:, np.newaxis] * int_cov
  XpX_inv_post = np.linalg.inv(np.dot(wols_x_post.T, int_cov) / n)
  asy_lin_rep_ols_post = np.dot(wols_eX_post, XpX_inv_post)
  
  inf_treat_pre = (reg_att_treat_pre - w_treat_pre * eta_treat_pre) / np.mean(w_treat_pre)
  inf_treat_post = (reg_att_treat_post - w_treat_post * eta_treat_post) / np.mean(w_treat_post)
  inf_treat = inf_treat_post - inf_treat_pre
  
  inf_cont_1 = (reg_att_cont - w_cont * eta_cont)
  M1 = np.mean(w_cont[:, np.newaxis] * int_cov, axis=0)
  inf_cont_2_post = np.dot(asy_lin_rep_ols_post, M1)
  inf_cont_2_pre = np.dot(asy_lin_rep_ols_pre, M1)
  inf_control = (inf_cont_1 + inf_cont_2_post - inf_cont_2_pre) / np.mean(w_cont)
  
  reg_att_inf_func = (inf_treat - inf_control)
#   se_reg_att = np.std(reg_att_inf_func) / np.sqrt(n)
  se = np.std(reg_att_inf_func, ddof=1) / np.sqrt(n)
#   print(f"att: {reg_att} \t se: {se}")

  return reg_att, reg_att_inf_func


