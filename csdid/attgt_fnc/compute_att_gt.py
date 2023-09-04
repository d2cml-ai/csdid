import numpy as np, pandas as pd
import patsy 
from drdid import drdid, reg_did

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

  # g, t = glist[0], tlist[0]

  for _, g, in enumerate(glist):
    # g = glist[1]
    G_main = (data[gname] == g)
    data = data.assign(G_m = 1 * G_main)

    for t_i in range(tlist_len):
      pret = t_i
      tn = tlist[t_i + tfac]
      if base_period == 'universal' or g < tn:
        try:
          pret = np.where(tlist + anticipation < g)[0][-1]
        except:
          raise f"There are no pre-treatment periods for the group first treated at {g}\nUnits from this group are dropped"
          # break


      if base_period == 'universal':
        if pret == tn:
          add_att_data()

      if not never_treated:
        n1 = data[gname] == 0
        n2 = (data[gname] > (tlist[np.max([t_i, pret]) + tfac]) + anticipation)
        n3 = np.where(data[gname] != glist[g], True, False)
        row_eval = n1 | n2 & n3
        data = data.assign(C = 1 * row_eval)

      post_treat = 1 * (g <= tn)
      disdat = data[(data[tname] == tn) | (data[tname] == tlist[pret])]

      if panel: 
        disdat = panel2cs2(disdat, yname, idname, tname)
        disdat = disdat.dropna()
        n = len(disdat)
        dis_idx = np.array(disdat.G_m == 1) | np.array(disdat.C == 1)
        disdat = disdat.loc[dis_idx, :]
        n1 = len(disdat)
        G = disdat.G_m
        C = disdat.C
        w = disdat.w

        ypre = disdat.y0 if tn > pret else disdat.y1
        ypost = disdat.y0 if tn < pret else disdat.y1
        _, covariates = fml(xformla, data = disdat, return_type = 'dataframe')

        G, C, w, ypre = map(np.array, [G, C, w, ypre])
        ypost, covariates = map(np.array, [ypost, covariates])

        if callable(est_method):
          est_att_f = est_method
        elif est_method == "reg":
          est_att_f = reg_did.reg_did_panel
        elif est_method == "dr":
          est_att_f = drdid.drdid_panel

        att_gt, att_inf_func = est_att_f(ypost, ypre, G, i_weights=w, covariates=covariates)

        inf_zeros = np.zeros(n)
        att_inf = n / n1 * att_inf_func
        inf_zeros[dis_idx] = att_inf

        add_att_data(att_gt, inf_f=inf_zeros)

      if not panel:
        right_ids = np.array(disdat.query('(G_m == 1) or (C == 1)').rowid.to_numpy())
        dis_idx = (data['rowid'].isin(right_ids)) &\
          ((data[tname] == tlist[t_i + tfac]) |\
            (data[tname] == tlist[pret]))

        disdat = data.loc[dis_idx]

        G = disdat.G_m.to_numpy()
        C = disdat.C.to_numpy()
        Y = disdat[yname].to_numpy()
        post = 1 * (disdat[tname] == tlist[t_i + tfac]).to_numpy()
        w = disdat.w.to_numpy()

        # G, C, Y, post, w = map(np.array, [G, C, Y, post, w])


        n1 = sum(G + C)

        skip_this_att_gt = False
        if np.sum(G * post) == 0:
          print(f"No units in group {g} in time period {tn}")
          skip_this_att_gt = True 

        if np.sum(G * (1 - post)) == 0:
          print(f"No units in group {g} in time period {t}")
          skip_this_att_gt = True 

        if np.sum(C * post) == 0:
          print(f"No available control units for group {g} in time period {tn}")
          skip_this_att_gt = True 

        if np.sum(C * (1 - post)) == 0:
          print(f"No available control units for group {g} in time period {t}")
          skip_this_att_gt = True 

        if skip_this_att_gt:
          add_att_data()

        try:
          _, covariates = fml(xformla, data = disdat, return_type = 'dataframe')
          covariates = np.array(covariates)
        except:
          y_str, x_str = xformla.split("~")
          xs1 = x_str.split('+')
          xs1_col_names = [x.strip() for x in xs1 if x.strip() != '1']
          n_dis = len(disdat)
          ones = np.ones((n_dis, 1))
          try:
            covariates = disdat[xs1_col_names].to_numpy()
            covariates = np.append(covariates, ones, axis=1)
          except:
            covariates = ones

        

        if callable(est_method):
          est_att_f = est_method
        elif est_method == "reg":
          est_att_f = reg_did.reg_did_rc
        elif est_method == "dr":
          est_att_f = drdid.drdid_rc
        att_gt, att_inf_func = est_att_f(y=Y, post=post, D = G, i_weights=w, covariates=covariates)

        inf_func_df = pd.DataFrame(
          {
            "inf_func": att_inf_func,
            "right_ids": right_ids
          }
        )
        inf_zeros = np.zeros(n)
        aggte_infffuc = inf_func_df.groupby('right_ids').inf_func.sum()
        try:
          dis_idx1 = np.isin(data['rowid'].unique(), aggte_infffuc.index.to_numpy())
        except:
          dis_idx1 = np.isin(data['rowid'].unique().to_numpy(), aggte_infffuc.index.to_numpy())
        inf_zeros[dis_idx1] = np.array(aggte_infffuc)

        add_att_data(att_gt, pst = post_treat, inf_f=inf_zeros)
        # print(att_est)

  output = {
    'group': group ,
    'year': year,
    "att" : att_est,
    'post ': post_array
  }
  return (output, np.array(inf_func))