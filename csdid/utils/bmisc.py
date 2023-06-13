import pandas as pd
def makeBalancedPanel(data, idname, tname):
  data = data.sort_values([idname, tname]).reset_index(drop = True)
  nt = len(data[tname].unique())
  data = data.groupby(idname)\
    .filter(lambda x: len(x) == nt)
  return data


def panel2cs2(data: pd.DataFrame, yname, idname, tname):
  if len(data[tname].unique()) != 2:
    raise ValueError('panel2cs2 only for 2 perios of apnel data')

  data = data.sort_values([idname, tname])
  y1 = data.groupby(idname)[yname].shift(-1)
  y0 = data[yname]
  dy = y1 - y0
  data = data.assign(
    y1 = y1, y0 = y0, dy = dy
  )
  return data.dropna()

# -*- coding: utf-8 -*-
"""
Created on Wed May 31 18:58:35 2023

@author: Carlos Guevara
"""
import numpy as np

def TorF(cond, use_isTRUE=False):

    if not isinstance(cond, np.ndarray) or cond.dtype != bool:
        raise ValueError("cond should be a logical vector")
    if use_isTRUE:
        cond = np.array([x is True for x in cond])
    else:
        cond[np.isnan(cond)] = False
    return cond

def multiplier_bootstrap(inf_func, biters): # This function comes from c++
    n, K = inf_func.shape
    biters = int(biters)
    innerMat = np.zeros((n, K))
    Ub = np.zeros(n)
    outMat = np.zeros((biters,K))

    for b in range(biters):
        # draw Rademechar weights
        # Ub = ( np.ones(n) - 2 * np.round(np.random.rand(n)) )[:, np.newaxis]
        Ub = np.random.choice([1, -1], size=(n, 1))
        innerMat = inf_func * Ub
        outMat[b] = np.mean(innerMat, axis=0)

    return outMat