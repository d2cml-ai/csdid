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

def multiplier_bootstrap(inf_func, biters):
    """Multiplier bootstrap using Rademacher weights (vectorized, chunked for large inputs)."""
    n, K = inf_func.shape
    biters = int(biters)
    # Chunk to limit peak memory: each chunk allocates (chunk_size x n) array
    max_chunk = max(1, min(biters, int(5e8 / max(n, 1))))  # ~500M elements max
    outMat = np.empty((biters, K))
    pos = 0
    while pos < biters:
        chunk = min(max_chunk, biters - pos)
        Ub = np.random.choice([1, -1], size=(chunk, n))
        outMat[pos:pos+chunk] = Ub @ inf_func / n
        pos += chunk
    return outMat