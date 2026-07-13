import numpy as np
import statsmodels.api as sm

lm = sm.WLS
glm = sm.GLM
n_x = np.newaxis
qr_solver = np.linalg.pinv
binomial = sm.families.Binomial()
mean = np.mean

def has_intercept(covariates, n):
  int_cov = np.ones(n)
  if covariates is not None:
    if np.all(covariates[:, 0] == int_cov):
      int_cov = covariates
    else:
      int_cov = np.concatenate((np.ones((n, 1)), covariates), axis=1)
  
  return int_cov

def has_weights(weights, n):
  if weights is None:
    weights = np.ones(n)
  elif np.min(weights) < 0:
    raise ValueError('weights must be non-negative')

def asy_lin_wols(w, d, x, y, out_y, pst = 1):
  n = len(y)
  w_ols = w * d * pst
  wols_x = w_ols[:, n_x] * x
  w_diff = w_ols * (y - out_y)
  wols_ex = w_diff[:, n_x] * x
  cr = np.dot(wols_x.T, x) / n
  xpx_inv = qr_solver(cr)
  asy_lin = np.dot(wols_ex, xpx_inv)
  return asy_lin

def out_wols(y, d, x, rows, wgts, pst=1):
  ols_coef = lm(
    y[rows], x[rows], weights=wgts[rows]
  ).fit().params
  out_y = np.dor(ols_coef, x.T)

  ### asy_lin_rep
  asy_lin_rep = asy_lin_wols(
    wgts, d, x, y, out_y, pst
  )


  return out_y, asy_lin_rep

def w_tc_val(w, d, pst=1):
  return w * d * pst

def eta_val(att, w_tc=1, y=None):
  if y is None:
    eta_r = mean(att) / mean(w_tc)
  else:
    eta_r = att * y / mean(att)


def inf_treat_f(att, w, eta):
  inf_f_u = att - w * eta
  inf_f = inf_f_u / mean(w)
  return inf_f

def bstrap_se():
  # replace with mboot
  pass

def fit_ps(d, x, w, post = 1):
  n = len(d)
  glm_fit = glm(
    d, x, family=binomial, freq_weights=w
  ).fit()
  ps_fit = glm_fit.fittedvalues
  ps_fit = np.minimum(ps_fit, 1 - 1e-16)

  rest_cont = w * ps_fit * (1 - d)
  post1 = (1 - post)
  ps_fit1 = (1 - ps_fit)
  # Rc
  w_cont_pre = rest_cont * post1 / ps_fit1 
  w_cont_post = rest_cont * post / ps_fit
  # panel
  w_cont = rest_cont / ps_fit1

  w_ps = w * (d - ps_fit)
  score_ps = w_ps[:, n_x] * x
  hessian = glm_fit.cov_params() * n

  asy_lin_ps = np.dot(score_ps, hessian)

  return (ps_fit, w_cont_pre, w_cont, w_cont_post, 0, asy_lin_ps)
