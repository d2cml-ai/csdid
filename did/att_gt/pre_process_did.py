import pandas as pd, numpy as np
from DID_params import DIDparams

def pre_process_did(
  yname, tname, idname, gname, data: pd.DataFrame, 
  xformla = None, weights_names = None, 
  control_group = ['nevertreated']) -> list:

  control_group = control_group[0]
  if control_group not in ["nevertreated", "notyettreated"]:
    raise ValueError("control_group must be either 'nevertreated' or 'notyettreated'")

  # make sure dataset is a DataFrame
  # this gets around Python's default of reading data as a different data structure
  # if not isinstance(data, pd.DataFrame):
  #   data = pd.DataFrame(data)

  # make sure time periods are numeric
  if not data[tname].dtype == "float64":
    raise ValueError("data[tname] must be numeric")

  # make sure gname is numeric
  if not data[gname].dtype == "float64":
    raise ValueError("data[gname] must be numeric")

  # put in blank xformla if no covariates
  if xformla is None:
    xformla = "~1"
  
  columns = [idname, tname, yname, gname, weightsname, clustervars]
  #  drop irrelevant columns from data
  # data <- cbind.data.frame(data[,c(idname, tname, yname, gname, weightsname, clustervars)], model.frame(xformla, data=data, na.action=na.pass))

  data = data[columns]
  
  n_orig, c_orig = data.shape

  data = data.dropna()
  n_diff = n_orig - data.shape[0]
  if weights_names is None: 
    w = np.ones(data.shape[0])
  else:
    w = data[weights_names]


  if "w" in data.columns:
    raise ValueError("`did` tried to use column named 'w' internally, but there was already a column with this name")

  data = data.assign(w = w)
  tlist = data[tname].unique()
  tlist = sorted(tlist)

  asif_never_treated = data[gname] > max(tlist)
  asif_never_treated.fillna(False, inplace=True)
  data.loc[asif_never_treated, gname] = 0

  glist = data[gname].unique()
  glist = sorted(glist)
  