
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

  return dp

