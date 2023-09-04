
import pandas as pd, numpy as np
import patsy 
from csdid.utils.bmisc import makeBalancedPanel

fml = patsy.dmatrices

def pre_process_did(yname, tname, idname, gname, data: pd.DataFrame, 
  control_group = ['nevertreated', 'notyettreated'], 
  anticipation = 0, xformla : str = None,
  panel = True, allow_unbalanced_panel = True, cband = False,
  clustervar = None,  weights_name = None
  ) -> dict:

  n, t = data.shape
  control_group = control_group[0]
  columns = [idname, tname, yname, gname]
  # print(columns)
  # Columns
  if clustervar is not None:
    columns += [clustervar]
  if weights_name is not None:
    columns += [weights_name]
    w = data[weights_name]
  else:
    w = np.ones(n)


  if xformla is None:
    xformla = f'{yname} ~ 1'

  # if xformla is None:
  try:
    xformla = f'{yname} ~ 1'
    _, x_cov = fml(xformla, data = data, return_type='dataframe')
    _, n_cov = x_cov.shape
    data = pd.concat([data[columns], x_cov], axis=1)
    data = data.assign(w = w)
  except:
    data = data.assign(intercept = 1)
    clms = columns + ['intercept']
    n_cov = len(data.columns)
    # patsy dont work with pyspark
    data = data[clms]
    if weights_name is None:
      data = data.assign(w = 1)
    else:
      data = data.assign(w = lambda x: x[weights_name] * 1)


  data = data.dropna()
  ndiff = n - len(data) 
  if ndiff != 0: 
    print(f'dropped, {ndiff}, rows from original data due to missing data')
  try:

    tlist = np.sort(data[tname].unique())
    glist = np.sort(data[gname].unique())
  except:
    tlist = np.sort(data[tname].unique().to_numpy())
    glist = np.sort(data[gname].unique().to_numpy())

  asif_nev_treated = data[gname] > np.max(tlist)
  asif_nev_treated.fillna(False, inplace=True)
  data.loc[asif_nev_treated, gname] = 0

  if len(glist[glist == 0]) == 0:
    if control_group == "nevertreated":
      raise ValueError("There is no available never-treated group")
    else:
      value = np.max(glist) - anticipation
      data = data.query(f'{tname} < @value')
      tlist = np.sort(data[tname].unique())
      glist = np.sort(data[gname].unique())
      glist = glist[glist < np.max(glist)]

  glist = glist[glist > 0]
  # first prerios 
  fp = tlist[0]
  glist = glist[glist > fp + anticipation]

  treated_fp = (data[gname] <= fp) & ~(data[gname] == 0)
  treated_fp.fillna(False, inplace=True)

  try:

    nfirst_period = np.sum(treated_fp) if panel \
      else len(data.loc[treated_fp, idname].unique())
  except:
    nfirst_period = treated_fp.sum() if panel \
      else len(data.loc[treated_fp, idname].unique())

  if nfirst_period > 0:
    warning_message = f"Dropped {nfirst_period} units that were already treated in the first period."
    print(warning_message)
    glist_in = np.append(glist, [0])
    data = data.query(f'{gname} in @glist_in')
    tlist = np.sort(data[tname].unique())
    glist = np.sort(data[gname].unique())
    glist = glist[glist > 0]
    fp = tlist[0]
    glist = glist[glist > fp + anticipation]

  #todo: idname must be numeric
  true_rep_cross_section = False
  if not panel:
    true_rep_cross_section = True

  if panel:
    if allow_unbalanced_panel: 
      panel = False
      true_rep_cross_section = False
    else:
      keepers = data.dropna().index
      n = len(data[idname].unique)
      print(n)
      n_keep = len(data.iloc[keepers, idname].unique())

      if len(data.loc[keepers] < len(data)):
        print(f"Dropped {n-n_keep} observations that had missing data.")
        data = data.loc[keepers]
      # make balanced data set
      n_old = len(data[idname].unique())
      data = makeBalancedPanel(data, idname=idname, tname=tname)
      n = len(data[idname].unique())
      if len(data) == 0:
        raise ValueError("All observations dropped to convert data to balanced panel. Consider setting `panel=False` and/or revisit 'idname'.")
      if n < n_old:
        warnings.warn(f"Dropped {n_old-n} observations while converting to balanced panel.")
      tn = tlist[0]
      n = len(data.query(f'{tname} == @tn'))
  # add rowid
  if not panel:

    keepers = data.dropna().index.to_numpy()
    ndiff = len(data.loc[keepers]) - len(data)
    if len(keepers) == 0:
      raise "All observations dropped due to missing data problems."
    if ndiff < 0:
      mssg = f"Dropped {ndiff} observations that had missing data."
      data = data.loc[keepers]
    if true_rep_cross_section: 
      # fix: posible error
      data = data.assign(rowid = range(len(data)))
      idname = 'rowid'
    else:
      # r_id = np.array(data[idname])
      data = data.assign(rowid = lambda x: x[idname] * 1)
    
    n = len(data[idname].unique())

  data = data.sort_values([idname, tname])
  data = data.assign(w1 = lambda x: x['w'] * 1)
  # data.loc[:, ".w"] = data['w']
  if len(glist) == 0:
    raise f"No valid groups. The variable in '{gname}' should be expressed as the time a unit is first treated (0 if never-treated)."
  if len(tlist) == 2:
    cband = False
  gsize = data.groupby(data[gname]).size().reset_index(name="count")
  gsize["count"] /= len(tlist)

  reqsize = n_cov + 5
  gsize = gsize[gsize["count"] < reqsize]

  if len(gsize) > 0:
    gpaste = ",".join(map(str, gsize[gname]))
    warnings.warn(f"Be aware that there are some small groups in your dataset.\n  Check groups: {gpaste}.")

    if 0 in gsize[gname].to_numpy() and control_group == "nevertreated":
      raise "Never-treated group is too small, try setting control_group='notyettreated'."
  nT, nG = map(len, [tlist, glist])
  did_params = {
    'yname' : yname, 'tname': tname,
    'idname' : idname, 'gname': gname,
    'xformla': xformla, 'data': data,
    'tlist': tlist, 'glist': glist,
    'n': n, 'nG': nG, 'nT': nT,
    'control_group': control_group, 'anticipation': anticipation,
    'weights_name': weights_name, 'panel': panel,
    'true_rep_cross_section': true_rep_cross_section,
    'clustervars': clustervar
  }
  return did_params

