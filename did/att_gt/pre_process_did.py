import pandas as pd, numpy as np
from DID_params import DIDparams

def pre_process_did(yname, tname, idname, gname, data: pd.DataFrame, 
  control_group = ['nevertreated'], est_method = 'dr', base_period = 'varying',
  alp = 0.05, anticipation = 0, biters = 1000, cores = 1,
  panel = True, allow_unbalanced_panel = True, bstrap = False, cband = False,
  print_details = True, pl = False, 
  clustervars = None, call = None, xformla : str = None, weights_name = None
  ) -> dict:

  control_group = control_group[0]
  if control_group not in ["nevertreated", "notyettreated"]:
    raise ValueError("control_group must be either 'nevertreated' or 'notyettreated'")

  # make sure dataset is a DataFrame
  # this gets around Python's default of reading data as a different data structure
  # if not isinstance(data, pd.DataFrame):
  #   data = pd.DataFrame(data)

  # make sure time periods are numeric
  if not np.issubdtype(data[tname].dtype, np.number):
    raise ValueError("data[tname] must be numeric")

  x_len = 0
  if xformla is not None:
    y_str, x_str = xformla.split("~")
    x_str_cols = x_str.strip().split("+")
    x_str_cols = [x.strip() for x in x_str_cols]
    x_len = len(x_str_cols)

  # make sure gname is numeric
  # if not data[gname].dtype == "float64":
  if not np.issubdtype(data[gname].dtype, np.number):
    raise ValueError("data[gname] must be numeric")

  # put in blank xformla if no covariates
  if xformla is None:
    xformla = "~1"


  columns = [idname, tname, yname, gname]

  if clustervars is not None:
    columns = columns + [clustervars] 


  if weights_name is None: 
    w = np.ones(data.shape[0])
  else:
    w = data[weights_name]
    columns = columns + [weights_name]

  if "w" in data.columns:
    raise ValueError("`did` tried to use column named 'w' internally, but there was already a column with this name")
 
  #  drop irrelevant columns from data
  # data <- cbind.data.frame(data[,c(idname, tname, yname, gname, weightsname, clustervars)], model.frame(xformla, data=data, na.action=na.pass))

  data = data[columns]
  n_orig, c_orig = data.shape

  data = data.dropna()
  n_diff = n_orig - data.shape[0]


  data = data.assign(w = w)

  tlist = np.unique(data[tname])
  tlist = np.sort(tlist)

  asif_never_treated = data[gname] > max(tlist)
  asif_never_treated.fillna(False, inplace=True)
  data.loc[asif_never_treated, gname] = 0

  glist = np.unique(data[gname])
  glist = np.sort(glist)
  if len(glist[glist == 0]) == 0:
    if control_group == "nevertreated":
        raise ValueError("There is no available never-treated group")
    else:
        # Drop all time periods with time periods >= latest treated
        data = data[data[tname] < (max(glist) - anticipation)]
        
        # Replace last treated time with zero
        # lines_gmax = data[gname] == max(glist, na_rm=True)
        # data.loc[lines_gmax, gname] = 0

        tlist = sorted(data[tname].unique())
        glist = sorted(data[gname].unique())

        # Don't compute ATT(g,t) for groups that are only treated at end
        # and only play a role as a comparison group
        glist = [val for val in glist if val < max(glist)]

  glist = glist[glist > 0]

  # drop groups treated in the first period or before
  first_period = tlist[0]
  glist = glist[glist > first_period + anticipation]

  # check for groups treated in the first period and drop these
  treated_first_period = (data[gname] <= first_period) & ~(data[gname] == 0)
  treated_first_period.fillna(False, inplace=True)
  nfirstperiod = len(data.loc[treated_first_period, idname].unique()) if panel else len(data.loc[treated_first_period])
  if nfirstperiod > 0:
      warning_message = f"Dropped {nfirstperiod} units that were already treated in the first period."
      print(warning_message)
      data = data[data[gname].isin([0] + glist)]
      # update tlist and glist
      tlist = sorted(data[tname].unique())
      glist = sorted(data[gname].unique())
      glist = [val for val in glist if val > 0]

      # drop groups treated in the first period or before
      first_period = tlist[0]
      glist = [val for val in glist if val > first_period + anticipation]

  if idname is not None:
    # make sure id is numeric
    if not data[idname].dtype.kind in ["i", "f"]:
      raise ValueError("data[idname] must be numeric")

  true_repeated_cross_sections = False
  if not panel:
      true_repeated_cross_sections = True

  #-----------------------------------------------------------------------------
  # setup data in panel case
  #-----------------------------------------------------------------------------
  if panel:

    # check for unbalanced panel
    if allow_unbalanced_panel:

      # code will run through repeated cross sections, so set panel to be False
      panel = False
      true_repeated_cross_sections = False

      if not data[idname].dtype.kind in ["i", "f"]:
          raise ValueError("Must provide a numeric id")

    else:

      # this is the case where we coerce balanced panel

      # check for complete cases
      keepers = data.dropna().index
      n = len(data[idname].unique())
      n_keep = len(data.loc[keepers, idname].unique())
      if len(data.loc[keepers, :]) < len(data):
        warnings.warn(f"Dropped {n-n_keep} observations that had missing data.")
        data = data.loc[keepers, :]

      # make it a balanced data set
      n_old = len(data[idname].unique())
      data = data.groupby([idname, tname]).first().reset_index()
      n = len(data[idname].unique())
      if n < n_old:
        warnings.warn(f"Dropped {n_old-n} observations while converting to balanced panel.")

      # If drop all data, you do not have a panel.
      if len(data) == 0:
        raise ValueError("All observations dropped to convert data to balanced panel. Consider setting `panel=False` and/or revisit 'idname'.")

      n = len(data[data[tname] == tlist[0]])
  if not panel:

    # check for complete cases
    keepers = data.dropna().index
    if len(data.loc[keepers, :]) < len(data):
      warnings.warn(f"Dropped {len(data) - len(data.loc[keepers, :])} observations that had missing data.")
      data = data.loc[keepers, :]

    # If drop all data, you do not have a panel.
    if len(data) == 0:
      raise ValueError("All observations dropped due to missing data problems.")

    # n-row data.frame to hold the influence function
    if true_repeated_cross_sections:
      data[".rowid"] = range(1, len(data) + 1)
      idname = ".rowid"
    else:
      # set rowid to idname for repeated cross section/unbalanced
      data[".rowid"] = data[idname]

    # n is unique number of cross section observations
    # this is different for repeated cross sections and unbalanced panel
    n = len(data[idname].unique())
  if len(glist) == 0:
    raise ValueError("No valid groups. The variable in 'gname' should be expressed as the time a unit is first treated (0 if never-treated).")

  # if there are only two time periods, then uniform confidence
  # bands are the same as pointwise confidence intervals
  if len(tlist) == 2:
      cband = False

  #-----------------------------------------------------------------------------
  # more error handling after we have balanced the panel

  # check against very small groups
  gsize = data.groupby(data[gname]).size().reset_index(name="count")
  gsize["count"] /= len(tlist)

  # how many in each group before give warning
  # 5 is just a buffer, could pick something else, but seems to work fine
  # print(gsize)
  reqsize = x_len + 5

  # which groups to warn about
  gsize = gsize[gsize["count"] < reqsize]

  # warn if some groups are small
  if len(gsize) > 0:
    gpaste = ",".join(map(str, gsize[gname]))
    warnings.warn(f"Be aware that there are some small groups in your dataset.\n  Check groups: {gpaste}.")

    if 0 in gsize[gname].tolist() and control_group == "nevertreated":
      raise ValueError("Never-treated group is too small, try setting control_group='notyettreated'.")
  
  nT = len(tlist)
  nG = len(glist)

  data = data.sort_values([idname, tname])

  did_params = {
    'yname' : yname,
    'tname' : tname,
    'idname' : idname,
    'gname' : gname,
    'xformla' : xformla,
    'data' : data,
    'control_group' : control_group,
    'anticipation' : anticipation,
    'weights_name' : weights_name,
    'alp' : alp,
    'bstrap' : bstrap,
    'biters' : biters,
    'clustervars' : clustervars,
    'cband' : cband,
    'print_details' : print_details,
    'pl' : pl,
    'cores' : cores,
    'est_method' : est_method,
    'base_period' : base_period,
    'panel' : panel,
    'true_repeated_cross_sections' : true_repeated_cross_sections,
    'n' : n,
    'nG' : nG,
    'nT' : nT,
    'tlist' : tlist,
    'glist' : glist,
    'call' : call
  }
  return did_params

  return 


data = pd.read_csv("../../data/mpdta.csv")

yname = "lemp"
gname = "first.treat"
idname = "countyreal"
tname = "year"
est_method = "reg"

# print(pre_process_did(yname, tname, idname, gname, data))