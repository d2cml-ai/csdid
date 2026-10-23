import pandas as pd
def makeBalancedPanel(data, idname, tname):
  data = data.sort_values([idname, tname]).reset_index(drop = True)
  nt = len(data[tname].unique())
  counts = data.groupby(idname)[tname].transform('size')
  return data.loc[counts.eq(nt)]


def panel2cs2(data: pd.DataFrame, yname, idname, tname, assume_sorted=False):
  """Reshape a two-period panel, optionally reusing its unit/time ordering."""
  if len(data[tname].unique()) != 2:
    raise ValueError('panel2cs2 only for 2 perios of apnel data')

  if not assume_sorted:
    data = data.sort_values([idname, tname])
  ids = data[idname]
  # Compare adjacent values without a shift that casts large integer IDs to
  # float64, or index alignment that can pair rows with duplicate labels.
  adjacent = pd.array(ids.array[:-1] == ids.array[1:], dtype='boolean')
  next_same_id = np.r_[adjacent.fillna(False).to_numpy(dtype=bool), False]
  next_same_id &= ids.notna().to_numpy()
  y1 = data[yname].shift(-1).where(next_same_id)
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

def multiplier_bootstrap(inf_func, biters, seed=None):
    """Multiplier bootstrap using Rademacher weights (vectorized, chunked for large inputs).

    ``seed`` (v10-F1): when ``None`` the draws come from the GLOBAL ``np.random``
    state (so the serial path stays reproducible under ``np.random.seed`` and
    bit-identical to the legacy behavior). When an integer/``SeedSequence`` is
    given, a LOCAL ``np.random.Generator`` is used instead -- this lets the
    parallel path seed each worker deterministically (derived from the parent
    global RNG) so the parallel bootstrap is reproducible AND its workers stay
    independent.
    """
    n, K = inf_func.shape
    biters = int(biters)
    rng = np.random.default_rng(seed) if seed is not None else None
    # Chunk to limit peak memory: each chunk allocates (chunk_size x n) array
    max_chunk = max(1, min(biters, int(5e8 / max(n, 1))))  # ~500M elements max
    outMat = np.empty((biters, K))
    pos = 0
    while pos < biters:
        chunk = min(max_chunk, biters - pos)
        if rng is not None:
            Ub = rng.integers(0, 2, size=(chunk, n))
        else:
            Ub = np.random.randint(0, 2, size=(chunk, n))
        Ub *= -2
        Ub += 1
        outMat[pos:pos+chunk] = Ub @ inf_func / n
        pos += chunk
    return outMat
