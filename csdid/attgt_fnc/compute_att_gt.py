import numpy as np, pandas as pd
import patsy 
from drdid import drdid, reg_did, ipwd_did

from csdid.utils.bmisc import panel2cs2

fml = patsy.dmatrices


# data Struct
# output = {
#   "att" : []
#   'group': []
#   'year': []
#   'post ': []
# }


def compute_att_gt(dp, est_method = "dr", base_period = 'varying'):
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

  tlist_len, tfac = len(tlist), 0
  if base_period != 'universal':
    tlist_len = tlist_len - 1
    tfac = 1

  inf_func = []

  att_est, group, year, post_array = [], [], [], []

  def add_att_data(att = 0, pst = 0, inf_f = []):
    inf_func.append(inf_f)
    att_est.append(att)
    group.append(g)
    year.append(tn)
    post_array.append(pst)

  never_treated = control_group == 'nevertreated'
  if never_treated:
    data = data.assign(C = 1 * (data[gname] == 0))
  data = data.assign(y_main = data[yname])

  return (dp,inf_func)