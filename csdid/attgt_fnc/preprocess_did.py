# from typing import Union, List, Optional, Dict
# import pandas as pd, numpy as np
# import patsy 
# from csdid.utils.bmisc import makeBalancedPanel

# fml = patsy.dmatrices

# def pre_process_did(yname, tname, idname, gname, data: pd.DataFrame, 
#   control_group: Union[str, List[str]] = "nevertreated", 
#   anticipation = 0, xformla : str = None,
#   panel = True, allow_unbalanced_panel = True, cband = False,
#   clustervar = None,  weights_name = None
#   ) -> dict:

#   n, t = data.shape
#   # control_group = "nevertreated"
#   columns = [idname, tname, yname, gname]
#   # print(columns)
#   # Columns
#   if clustervar is not None:
#     columns += [clustervar]
#   if weights_name is not None:
#     columns += [weights_name]
#     w = data[weights_name]
#   else:
#     w = np.ones(n)


#   if xformla is None:
#     xformla = f'{yname} ~ 1'

#   # if xformla is None:
#   try:
#     _, x_cov = fml(xformla, data = data, return_type='dataframe')
#     _, n_cov = x_cov.shape
#     data = pd.concat([data[columns], x_cov], axis=1)
#     data = data.assign(w = w)
#   except:
#     data = data.assign(intercept = 1)
#     clms = columns + ['intercept']
#     n_cov = len(data.columns)
#     # patsy dont work with pyspark
#     data = data[clms]
#     if weights_name is None:
#       data = data.assign(w = 1)
#     else:
#       data = data.assign(w = lambda x: x[weights_name] * 1)


#   data = data.dropna()
#   ndiff = n - len(data) 
#   if ndiff != 0: 
#     print(f'dropped, {ndiff}, rows from original data due to missing data')
#   try:

#     tlist = np.sort(data[tname].unique())
#     glist = np.sort(data[gname].unique())
#   except:
#     tlist = np.sort(data[tname].unique().to_numpy())
#     glist = np.sort(data[gname].unique().to_numpy())

#   asif_nev_treated = data[gname] > np.max(tlist)
#   asif_nev_treated.fillna(False, inplace=True)
#   data.loc[asif_nev_treated, gname] = 0

#   if len(glist[glist == 0]) == 0:
#     if control_group == "nevertreated":
#       raise ValueError("There is no available never-treated group")
#     else:
#       value = np.max(glist) - anticipation
#       data = data.query(f'{tname} < @value')
#       tlist = np.sort(data[tname].unique())
#       glist = np.sort(data[gname].unique())
#       glist = glist[glist < np.max(glist)]

#   glist = glist[glist > 0]
#   # first prerios 
#   fp = tlist[0]
#   glist = glist[glist > fp + anticipation]

#   treated_fp = (data[gname] <= fp) & ~(data[gname] == 0)
#   treated_fp.fillna(False, inplace=True)

#   try:

#     nfirst_period = np.sum(treated_fp) if panel \
#       else len(data.loc[treated_fp, idname].unique())
#   except:
#     nfirst_period = treated_fp.sum() if panel \
#       else len(data.loc[treated_fp, idname].unique())

#   if nfirst_period > 0:
#     warning_message = f"Dropped {nfirst_period} units that were already treated in the first period."
#     print(warning_message)
#     glist_in = np.append(glist, [0])
#     data = data.query(f'{gname} in @glist_in')
#     tlist = np.sort(data[tname].unique())
#     glist = np.sort(data[gname].unique())
#     glist = glist[glist > 0]
#     fp = tlist[0]
#     glist = glist[glist > fp + anticipation]

#   # (idname numeric check now performed in _validate_inputs)
#   true_rep_cross_section = False
#   if not panel:
#     true_rep_cross_section = True



#   #-----------------------------------------------------------------------------
#   # setup data in panel case
#   #-----------------------------------------------------------------------------

#   # # Check if data is a balanced panel if panel = True and allow_unbalanced_panel = True
#   # if panel and allow_unbalanced_panel:
#   #     # First, focus on complete cases
#   #     keepers = data.dropna()
#   #     data_comp = keepers.copy()

#   #     # Make it a balanced dataset
#   #     n_all = len(data_comp[idname].unique())
#   #     data_bal = makeBalancedPanel(data_comp, idname, tname)
#   #     n_bal = len(data_bal[idname].unique())

#   #     if n_bal < n_all:
#   #         # If fewer unique IDs in the balanced panel, it means the panel is unbalanced
#   #         print("You have an unbalanced panel. Proceeding as such.")
#   #         allow_unbalanced_panel = True
#   #     else:
#   #         # If the number of unique IDs remains the same, it means the panel is balanced
#   #         print("You have a balanced panel. Setting allow_unbalanced_panel = False.")
#   #         allow_unbalanced_panel = False


#   if panel:
#     if allow_unbalanced_panel: 
#       panel = False
#       true_rep_cross_section = False
#     else:
#       keepers = data.dropna().index
#       n = len(data[idname].unique)
#       print(n)
#       n_keep = len(data.iloc[keepers, idname].unique())

#       if len(data.loc[keepers] < len(data)):
#         print(f"Dropped {n-n_keep} observations that had missing data.")
#         data = data.loc[keepers]
#       # make balanced data set
#       n_old = len(data[idname].unique())
#       data = makeBalancedPanel(data, idname=idname, tname=tname)
#       n = len(data[idname].unique())
#       if len(data) == 0:
#         raise ValueError("All observations dropped to convert data to balanced panel. Consider setting `panel=False` and/or revisit 'idname'.")
#       if n < n_old:
#         warnings.warn(f"Dropped {n_old-n} observations while converting to balanced panel.")
#       tn = tlist[0]
#       n = len(data.query(f'{tname} == @tn'))

#   # # Add rowid
#   # if not panel:
#   #     # Handle missing data
#   #     keepers = data.dropna().index
#   #     dropped_count = len(data) - len(keepers)
#   #     if dropped_count > 0:
#   #         print(f"Dropped {dropped_count} observations that had missing data.")
      
#   #     # Drop incomplete rows
#   #     data = data.loc[keepers]

#   #     # Check if all rows are dropped
#   #     if len(data) == 0:
#   #         raise ValueError("All observations dropped due to missing data problems.")

#   #     # Add rowid based on true_repeated_cross_section
#   #     if true_rep_cross_section:
#   #         data = data.assign(rowid=range(len(data)))
#   #         idname = 'rowid'
#   #     else:
#   #         if idname not in data.columns:
#   #             raise ValueError(f"Column {idname} is not in the dataset.")
#   #         data = data.assign(rowid=lambda x: x[idname] * 1)

#   #     # Number of unique observations
#   #     n = data[idname].nunique()

#   if not panel:
#     # Check for complete cases (equivalent to complete.cases in R)
#     keepers = data.notna().all(axis=1)
#     dropped_count = len(data) - keepers.sum()
    
#     if dropped_count > 0:
#       import warnings
#       warnings.warn(f"Dropped {dropped_count} observations that had missing data.")
#       data = data[keepers].copy()
    
#     # If drop all data, raise error
#     if len(data) == 0:
#       raise ValueError("All observations dropped due to missing data problems.")
    
#     # Add rowid column
#     if true_rep_cross_section:
#       data.loc[:, 'rowid'] = range(1, len(data) + 1)  # 1-based indexing like R
#       idname = 'rowid'
#     else:
#       # Set rowid to idname for repeated cross section/unbalanced
#       if idname not in data.columns:
#         raise ValueError(f"Column {idname} not found in dataset")
#       data.loc[:, 'rowid'] = data[idname]
    
#     # Calculate n as unique number of cross section observations
#     n = data[idname].nunique()

#   data = data.sort_values([idname, tname])
#   data = data.assign(w1 = lambda x: x['w'] * 1)
#   # data.loc[:, ".w"] = data['w']
#   if len(glist) == 0:
#     raise f"No valid groups. The variable in '{gname}' should be expressed as the time a unit is first treated (0 if never-treated)."
#   if len(tlist) == 2:
#     cband = False
#   gsize = data.groupby(data[gname]).size().reset_index(name="count")
#   gsize["count"] /= len(tlist)

#   reqsize = n_cov + 5
#   gsize = gsize[gsize["count"] < reqsize]

#   if len(gsize) > 0:
#     gpaste = ",".join(map(str, gsize[gname]))
#     warnings.warn(f"Be aware that there are some small groups in your dataset.\n  Check groups: {gpaste}.")

#     if 0 in gsize[gname].to_numpy() and control_group == "nevertreated":
#       raise "Never-treated group is too small, try setting control_group='notyettreated'."
#   nT, nG = map(len, [tlist, glist])
#   did_params = {
#     'yname' : yname, 'tname': tname,
#     'idname' : idname, 'gname': gname,
#     'xformla': xformla, 'data': data,
#     'tlist': tlist, 'glist': glist,
#     'n': n, 'nG': nG, 'nT': nT,
#     'control_group': control_group, 'anticipation': anticipation,
#     'weights_name': weights_name, 'panel': panel,
#     'true_rep_cross_section': true_rep_cross_section,
#     'clustervars': clustervar
#   }
#   return did_params


import re
import pandas as pd, numpy as np
import patsy
from pandas.api.types import is_numeric_dtype, is_bool_dtype
from csdid.utils.bmisc import makeBalancedPanel
import warnings

fml = patsy.dmatrices

# Internal column names that must not collide with user data.
# P-reserved (C2-F2): union the port's genuinely-used internal columns
# (undotted w/w1/rowid/... -- the port really creates these, so colliding on them
# is unsafe and being stricter than R there is fine) with R's reserved set
# (check_reserved_did_names: .w/.rowid/.G/.C/post/asif_never_treated/
# treated_first_period). Adding R's names closes the PORT-SILENT-GAP where the
# port accepted R-reserved names that R refuses. Residual: the port stays stricter
# than R on w/w1/G_m/C/y_main/rowid/y0/y1/dy (R reserves only the dotted forms).
_RESERVED_NAMES = {
  'w', 'w1', 'G_m', 'C', 'y_main', 'rowid', 'y0', 'y1', 'dy',
  '.w', '.rowid', '.G', '.C', 'post', 'asif_never_treated', 'treated_first_period',
}


def _normalize_name_arg(value, argname):
  """Mirror R `validate_character_scalar`: unwrap a length-1 list/tuple/set/ndarray
  to its scalar element and accept it; refuse length>=2 (or other multi-element /
  unhashable containers, e.g. dict/set) with R's wording BEFORE the later
  `value not in data.columns` membership test, which would otherwise crash with
  `unhashable type` (dict/set) or an ambiguous-truth error (ndarray). ``None``
  passes through unchanged so the caller's allow-null logic still applies
  (P-name-scalar / C2-F1)."""
  if value is None or isinstance(value, str):
    return value
  # Any non-string container: a length-1 list/tuple/set/ndarray unwraps to its
  # element; everything else (length!=1, dict, mapping, ...) is refused cleanly.
  if isinstance(value, (list, tuple, set, frozenset, np.ndarray)):
    seq = list(value)
    if len(seq) == 1:
      return seq[0]
    raise ValueError(f"{argname} must be a single non-missing character string.")
  if isinstance(value, dict):
    raise ValueError(f"{argname} must be a single non-missing character string.")
  return value


def _strictly_numeric(series):
  """True for numeric dtypes (including pandas nullable Int64/Float64), but not
  boolean. Uses pandas dtype checks so extension dtypes are handled and do not
  crash ``np.issubdtype`` (which raises on e.g. ``Int64Dtype``)."""
  return is_numeric_dtype(series) and not is_bool_dtype(series)


def _validate_inputs(yname, tname, idname, gname, data, control_group,
                     anticipation, panel, clustervar, weights_name):
  """Validate inputs to pre_process_did, raising clear errors on bad data.

  Returns a tuple of the normalized name args
  ``(yname, tname, gname, idname, weights_name, clustervar)``. Each name arg has a
  length-1 container unwrapped to its scalar (P-name-scalar / C2-F1); ``clustervar``
  additionally has a single-element list/tuple unwrapped to a scalar and an
  idname-only / empty one collapsed to ``None`` (P-clust-id), so callers store
  scalar identifiers consistently.
  """

  # --- Normalize name args (P-name-scalar / C2-F1): a length-1 list/tuple/ndarray
  # is unwrapped to its scalar (R `validate_character_scalar` accepts a length-1
  # character vector); a length>=2 / unhashable container is refused cleanly with
  # R's wording BEFORE the `value not in data.columns` membership test, which would
  # otherwise crash with `unhashable type` (list/set/dict) or an ambiguous-truth
  # error (ndarray). ---
  yname = _normalize_name_arg(yname, "yname")
  tname = _normalize_name_arg(tname, "tname")
  gname = _normalize_name_arg(gname, "gname")
  idname = _normalize_name_arg(idname, "idname")
  weights_name = _normalize_name_arg(weights_name, "weightsname")

  # --- panel requires an idname (R `did`: missing_idname). Must precede the
  # column-existence check so a None idname is not reported as a missing column
  # (V4-U5). ---
  if panel and idname is None:
    raise ValueError(
      "Must provide idname when panel=True. Set panel=False for repeated cross sections."
    )

  # --- Normalize clustervar: at most one cluster variable is supported beyond
  # the unit (R `did`). A raw list reaches the column-existence check below and
  # would crash with `unhashable type: 'list'` (V4-U3); refuse >1 cleanly and
  # unwrap a single-element list/tuple to a scalar. ---
  # P-clust-id (C2-F7): R removes idname from the cluster set BEFORE the
  # "at most one" count check (clustering at the unit level via idname is
  # implicit), so `clustervars = [idname, one_other]` reduces to 1 effective
  # cluster var and is accepted. Strip idname from a list/tuple clustervar first.
  if isinstance(clustervar, (list, tuple, set, frozenset, np.ndarray)):
    # Normalize any cluster-var container (list/tuple/set/ndarray) the same way; an
    # ndarray otherwise reaches a `.any()`-style truth test and raises the
    # ambiguous-truth error (P-name-scalar sibling).
    cv_list = list(clustervar)
    cv_check = [c for c in cv_list if c != idname] if idname is not None else cv_list
    if len(cv_check) > 1:
      raise ValueError(
        "At most one cluster variable (beyond 'idname') is supported. "
        "Please reduce to one."
      )
    clustervar = cv_check[0] if len(cv_check) == 1 else None
  elif clustervar is not None and idname is not None and clustervar == idname:
    # A scalar clustervar equal to idname is no effective cluster var (R drops it).
    clustervar = None

  # --- Check columns exist ---
  required = {
    'yname': yname, 'tname': tname, 'gname': gname
  }
  # P-idnone (C1-11 / C2-F5): R validates idname with allow_null = !panel, so a
  # None idname is allowed for repeated cross sections (panel=False) -- the
  # true-RCS path synthesizes `.rowid` below. Only require the idname column to
  # exist when idname is non-None (panel=True + None already refused above).
  if idname is not None:
    required['idname'] = idname
  if weights_name is not None:
    required['weights_name'] = weights_name
  if clustervar is not None:
    required['clustervar'] = clustervar

  missing = [f"{k}='{v}'" for k, v in required.items() if v not in data.columns]
  if missing:
    raise ValueError(
      f"Column(s) not found in data: {', '.join(missing)}. "
      f"Available columns: {list(data.columns)}"
    )

  # --- Time-varying cluster variable: only the first period's labels are used
  # internally, so a within-unit time-varying cluster var is silently mis-applied
  # (V4-S1). R refuses it up front; the fix belongs here in _validate_inputs, not
  # in the SE code (the N1 lesson). ---
  # P-panelgate (C2-F4): R runs this check whenever idname is non-null (n>=2),
  # NOT gated on panel (pre_process_did.R lines 329-376). Ungate it from `panel`
  # so it also fires on the RCS path; require idname for the per-unit grouping.
  if clustervar is not None and idname is not None and len(data) >= 2 and \
      data.groupby(idname)[clustervar].nunique().gt(1).any():
    raise ValueError(
      "Time-varying cluster variables are not supported. "
      "Please provide a time-invariant cluster variable."
    )

  # --- Reserved column name checks ---
  # The weights column is read into the internal 'w' column before any internal
  # columns are created, so a weights column named like an internal column (e.g. 'w')
  # is safe and is allowed (matches R 'did', which stores weights in '.w').
  user_cols = {yname, tname, idname, gname}
  if clustervar:
    user_cols.add(clustervar)
  conflicts = user_cols & _RESERVED_NAMES
  if conflicts:
    raise ValueError(
      f"Column name(s) {conflicts} conflict with names used internally by csdid. "
      f"Please rename these columns before calling ATTgt."
    )

  # --- Non-numeric outcome (allow bool, which R treats as logical) ---
  if not (is_numeric_dtype(data[yname]) or is_bool_dtype(data[yname])):
    raise ValueError(
      f"Outcome variable '{yname}' must be numeric, got dtype '{data[yname].dtype}'."
    )

  # --- tname / gname / idname must be numeric (not boolean). pandas dtype
  # checks accept nullable integers (Int64) and avoid np.issubdtype crashing on
  # extension dtypes. ---
  if not _strictly_numeric(data[tname]):
    raise ValueError(f"Time variable '{tname}' must be numeric, got dtype '{data[tname].dtype}'.")
  if not _strictly_numeric(data[gname]):
    raise ValueError(f"Group variable '{gname}' must be numeric, got dtype '{data[gname].dtype}'.")

  # --- Non-numeric idname (R did v2.5.1: idname must be numeric so the
  # influence function can be indexed by unit). Guarded on idname being provided
  # (P-idnone: panel=False allows idname=None, where .rowid is synthesized). ---
  if idname is not None and not _strictly_numeric(data[idname]):
    raise ValueError(
      f"The id variable '{idname}' must be numeric, got dtype '{data[idname].dtype}'. "
      f"Please convert it to numeric."
    )

  # --- Negative gname ---
  if (data[gname] < 0).any():
    raise ValueError(
      f"Group variable '{gname}' contains negative values. "
      f"Values must be 0 (never-treated) or a positive treatment time."
    )

  # --- Treatment irreversibility and (id, time) uniqueness ---
  # P-panelgate (C2-F4): R runs BOTH checks whenever idname is non-null and n>=2,
  # NOT gated on `panel` (pre_process_did.R lines 329-376). Ungating from `panel`
  # makes the RCS path (panel=False with a real idname) refuse duplicate (id,t)
  # rows and treatment reversal, matching R -- and is verified safe for the
  # faithful shared-id RCS cases (O2_E_rcs_*), whose data have no (id,t)
  # duplicates and time-invariant gname, so they still run. R checks
  # irreversibility before uniqueness, so order them the same way.
  if idname is not None and len(data) >= 2:
    g_tv = data.groupby(idname)[gname].nunique()
    if (g_tv > 1).any():
      raise ValueError(
        "The value of gname (treatment variable) must be the same across all "
        "periods for each particular unit. The treatment must be irreversible."
      )
    dups = data.duplicated(subset=[idname, tname], keep=False)
    if dups.any():
      n_dups = dups.sum()
      raise ValueError(
        f"Found {n_dups} duplicate (idname, tname) rows. "
        f"Each unit must appear at most once per time period. "
        f"Check for duplicate entries in your data."
      )

  # --- anticipation (P-anticipation / C1-2) ---
  # Mirror R `validate_anticipation`: numeric, length-1, finite, non-NA,
  # non-negative WHOLE number; reject bool (R: is.numeric(TRUE) is FALSE). Accept
  # numpy numeric scalars and length-1 numeric vectors / 0-d arrays (R accepts
  # those) by unwrapping first. The old guard accepted 0.5/1.5/NaN/Inf/True.
  _anti = anticipation
  if isinstance(_anti, (list, tuple, np.ndarray)):
    _seq = list(np.asarray(_anti).ravel())
    if len(_seq) != 1:
      raise ValueError(
        "anticipation must be a single finite non-missing number. "
        "Please check your arguments."
      )
    _anti = _seq[0]
  if isinstance(_anti, (bool, np.bool_)):
    raise ValueError("anticipation must be numeric. Please convert it.")
  if not isinstance(_anti, (int, float, np.integer, np.floating)):
    raise ValueError("anticipation must be numeric. Please convert it.")
  if not np.isfinite(_anti):
    raise ValueError(
      "anticipation must be a single finite non-missing number. "
      "Please check your arguments."
    )
  if _anti < 0 or float(_anti) != round(float(_anti)):
    raise ValueError(
      "anticipation must be a non-negative whole number. "
      "Please check your arguments."
    )
  # Use the unwrapped scalar downstream (a len-1 list/array would otherwise break
  # the broadcasted `data[gname] > max(tlist) + anticipation` comparison).
  anticipation = _anti

  # --- control_group ---
  valid_cg = ['nevertreated', 'notyettreated']
  if control_group not in valid_cg:
    raise ValueError(
      f"'control_group' must be one of {valid_cg}, got '{control_group}'."
    )

  # --- weights validation ---
  if weights_name is not None:
    w = data[weights_name]
    if (w < 0).any():
      raise ValueError(f"Weights column '{weights_name}' contains negative values.")
    if w.mean() <= 0:
      raise ValueError(f"Weights column '{weights_name}' has non-positive mean.")

  # Return the normalized name args + unwrapped anticipation alongside clustervar
  # so the caller uses the unwrapped scalars (P-name-scalar / C2-F1; P-anticipation).
  return yname, tname, gname, idname, weights_name, clustervar, anticipation


def pre_process_did(yname, tname, idname, gname, data: pd.DataFrame, 
  control_group = ['nevertreated', 'notyettreated'], 
  anticipation = 0, xformla : str = None,
  panel = True, allow_unbalanced_panel = True, cband = False,
  clustervar = None,  weights_name = None, fix_weights = None
  ) -> dict:

  n, t = data.shape
  # P-control-vec (C1-5): the att_gt default is the 2-list
  # ['nevertreated','notyettreated'] (a sentinel meaning "default"); element-0
  # unwrap of that default is correct. A USER-supplied multi-element control_group
  # that differs from the default is already refused upstream in att_gt.__init__
  # (the att_gt builder owns that file), so here we only need the element-0 unwrap
  # of the default / any length-1 container. (See DEFERRED note in report.)
  if isinstance(control_group, (list, tuple)):
    control_group = control_group[0]

  # fix_weights = 'base_period'/'first_period' are not supported for true
  # repeated cross sections (matches R `did`).
  input_panel = panel
  if (not input_panel) and fix_weights in ("base_period", "first_period"):
    raise ValueError(
      f"fix_weights = '{fix_weights}' is not supported for repeated cross sections. "
      f"Use fix_weights = 'varying' or None instead."
    )

  # Validate inputs early with clear error messages. Returns the normalized name
  # args (length-1 containers unwrapped; P-name-scalar / C2-F1) and clustervar
  # (single-element list unwrapped to a scalar; idname-only collapsed to None;
  # see V4-U3 / P-clust-id).
  yname, tname, gname, idname, weights_name, clustervar, anticipation = _validate_inputs(
    yname, tname, idname, gname, data, control_group,
    anticipation, panel, clustervar, weights_name)

  # P-msg-antic (C1-17 / O3-F4 / C3-F1): when anticipation > 0, R emits a message
  # (the K channel) explaining never-treated units are assumed not to anticipate.
  # The port emitted nothing; print to stdout (the port's K channel, like the
  # aggte summary). Uses the actual anticipation value.
  if anticipation and anticipation > 0:
    print(
      f"Note: anticipation = {anticipation}. Never-treated units (with group "
      f"status 0 or Inf) are assumed to never anticipate treatment. Anticipation "
      f"only applies to eventually-treated units."
    )

  # P-idnone (C1-11 / C2-F5): idname may be None for a true repeated cross section
  # (panel=False); drop it from the selected columns (R's c(idname,...) drops NULL)
  # -- the .rowid synthesis below supplies the unit index.
  columns = [c for c in (idname, tname, yname, gname) if c is not None]
  # print(columns)
  # Columns
  if clustervar is not None:
    columns += [clustervar]
  if weights_name is not None:
    columns += [weights_name]
    w = data[weights_name]
  else:
    w = np.ones(n)

  # P-dedup (C2-F6): de-duplicate the column list (preserving order). When an arg
  # aliases another column (e.g. clustervar==idname, yname==gname,
  # weightsname==gname), `data[columns]` would otherwise select a repeated column
  # and downstream `data[col]` would return a DataFrame -> AttributeError on
  # `.unique()/.to_numpy()`. R selects unique(keep_cols), so de-dup here too.
  columns = list(dict.fromkeys(columns))


  if xformla is None:
    xformla = f'{yname} ~ 1'

  # Pre-check formula variables against the data columns. A bare name in the
  # formula that is not a column otherwise surfaces downstream as a patsy
  # NameError wrapped in "Error processing formula ..." (unclassified); R reports
  # a clean missing-column error naming the variable(s) (V4-U4). Only flag bare
  # identifiers (e.g. `Ghost`) -- transformed/call terms like `np.log(X)` or
  # `C(x)` are left for patsy to resolve.
  # P-reqsize (O3-F3 / C2-F10): count the number of RHS *variables* in xformla
  # (0 for intercept-only ~1), mirroring R `length(BMisc::rhs_vars(xformla))`.
  # This excludes the patsy intercept, which the old `n_cov = x_cov.shape[1]`
  # wrongly counted (reqsize was 1 too high). Default to 0 so the fallback path
  # (intercept-only) is also correct.
  n_rhs = 0
  try:
    _desc = patsy.ModelDesc.from_formula(xformla)
    _fnames = {f.name() for side in (_desc.lhs_termlist, _desc.rhs_termlist)
               for term in side for f in term.factors}
    _missing_vars = sorted(
      nm for nm in _fnames
      if re.fullmatch(r'[A-Za-z_]\w*', nm) and nm not in data.columns
    )
    # Distinct RHS variable names, excluding the Intercept term (R's rhs_vars
    # likewise returns only the named variables, not the intercept).
    _rhs_vars = {f.name() for term in _desc.rhs_termlist for f in term.factors}
    n_rhs = len(_rhs_vars)
  except Exception:
    _missing_vars = []
  if _missing_vars:
    raise ValueError(
      "The following variables are not in data: " + ", ".join(_missing_vars) + "."
    )

  # if xformla is None:
  try:
    try:
      _, x_cov = fml(xformla, data=data, return_type='dataframe')
    except Exception:
      x_cov = patsy.dmatrix(xformla, data=data, return_type='dataframe')
    _, n_cov = x_cov.shape
    # Names of the (globally built, factor-consistent) design columns. These are
    # selected directly in the (g,t) loop instead of re-running patsy per cell,
    # which both fixes transformed/factor covariates and removes patsy from the
    # inner loop (matches R `did`, which builds model.matrix() once).
    xcov_cols = list(x_cov.columns)
    data = pd.concat([data[columns], x_cov], axis=1)
    data = data.assign(w=w)
  except Exception as e:
    if xformla is not None and xformla != f'{yname} ~ 1':
      raise ValueError(
        f"Error processing formula '{xformla}': {e}. "
        f"Check that all formula variables exist in the data."
      ) from e
    # Only fall back to intercept for the default formula case
    data = data.assign(intercept = 1)
    xcov_cols = ['intercept']
    clms = columns + ['intercept']
    n_cov = len(data.columns)
    # patsy dont work with pyspark
    data = data[clms]
    if weights_name is None:
      data = data.assign(w = 1)
    else:
      data = data.assign(w = lambda x: x[weights_name] * 1)


  # P-bool-y (C2-F8): R allows a logical/0-1 outcome (is.logical) and treats it as
  # numeric; the port's estimator calls np.isnan/arithmetic on a bool array, which
  # raises "ufunc 'isnan' not supported" -> all-NaN ATT. Cast a bool outcome to
  # float (0/1) here, after validation and before estimation. Only yname is cast;
  # gname/tname/id are already numeric-validated.
  if yname in data.columns and is_bool_dtype(data[yname]):
    data[yname] = data[yname].astype(float)

  # P-finite (C2-F3): R `complete_finite_cases(data, finite_exclude = gname)` drops
  # every row with a non-finite (+/-Inf or NaN) value in ANY numeric column EXCEPT
  # gname (gname=Inf is a valid never-treated code, handled by the gname>max recode
  # below). The old `data.dropna()` removed only NaN and kept +/-Inf, which then
  # poisoned the estimator. Replicate the finite filter: keep rows finite in every
  # numeric column other than gname, and also drop NaN in gname.
  _finite_mask = pd.Series(True, index=data.index)
  for _col in data.columns:
    if not is_numeric_dtype(data[_col]):
      continue
    _vals = data[_col]
    if _col == gname:
      # gname: NA/NaN dropped, but Inf preserved (valid never-treated code).
      _finite_mask &= _vals.notna()
    else:
      _finite_mask &= np.isfinite(_vals.to_numpy(dtype='float64', na_value=np.nan))
  # Non-numeric columns (e.g. a cluster label): drop NaN only (mirrors
  # complete.cases, which is not a finiteness test for non-numeric columns).
  for _col in data.columns:
    if not is_numeric_dtype(data[_col]):
      _finite_mask &= data[_col].notna()
  data = data[_finite_mask]
  # All rows dropped (e.g. a fully-missing model column) -> downstream reductions
  # like np.max(tlist) crash on the empty frame (V4-U2). The clean all-dropped
  # guard previously existed only on the balanced-panel branch.
  if len(data) == 0:
    raise ValueError(
      "All observations dropped due to missing data. Check your outcome, group, "
      "time, and covariate variables for missing values."
    )
  ndiff = n - len(data)
  if ndiff != 0:
    # v7-CF1: surface on the warnings channel (matches R `did`), not stdout.
    warnings.warn(f"dropped {ndiff} rows from original data due to missing or non-finite data")
  try:

    tlist = np.sort(data[tname].unique())
    glist = np.sort(data[gname].unique())
  except Exception:
    tlist = np.sort(data[tname].unique().to_numpy())
    glist = np.sort(data[gname].unique().to_numpy())

  asif_nev_treated = data[gname] > np.max(tlist) + anticipation
  asif_nev_treated.fillna(False, inplace=True)
  data.loc[asif_nev_treated, gname] = 0

  # Recompute glist after modifying gname (some cohorts may have been recoded to 0)
  try:
    glist = np.sort(data[gname].unique())
  except Exception:
    glist = np.sort(data[gname].unique().to_numpy())

  # Handle the case with no never-treated group.
  # R did v2.5.1 warns and coerces (does not error) when control_group='nevertreated'
  # and no never-treated units exist: it coerces the last cohort to never-treated and
  # truncates periods at/after that cohort.
  latest_g = None
  if not (glist == 0).any():
    latest_g = np.max(glist)
    cutoff_t = latest_g - anticipation
    if control_group == "nevertreated":
      warnings.warn(
        "No never-treated group is available. The last treated cohort is being "
        "coerced as 'never-treated' units, and data from periods at/after that cohort "
        "is being filtered out (no available comparison groups)."
      )
      data = data.query(f'{tname} < @cutoff_t')
      data.loc[data[gname] == latest_g, gname] = 0
    else:
      # notyettreated: drop late periods; the last cohort stays in the data as the
      # not-yet-treated comparison group.
      data = data.query(f'{tname} < @cutoff_t')
    tlist = np.sort(data[tname].unique())
    glist = np.sort(data[gname].unique())
    # The last cohort only serves as a comparison; exclude it from estimated groups.
    if control_group != "nevertreated":
      glist = glist[glist < latest_g]

  glist = glist[glist > 0]
  # P-empty (O3-F2): the cohort/anticipation/control trimming above can empty the
  # data (and hence tlist/glist). `tlist[0]` / np.max(...) on a zero-length array
  # raises IndexError; R degrades to a clean stop instead. Guard here and raise the
  # existing "No valid groups" error (the same terminal condition R reaches via its
  # length(glist)==0 stop).
  if len(tlist) == 0 or len(data) == 0:
    raise ValueError(
      f"No valid groups. The variable in '{gname}' should be expressed as the "
      f"time a unit is first treated (0 if never-treated)."
    )
  # first period
  fp = tlist[0]
  glist = glist[glist > fp + anticipation]

  # Identify units already treated in the first period (accounting for anticipation).
  treated_fp = (data[gname] <= fp + anticipation) & ~(data[gname] == 0)
  treated_fp = treated_fp.fillna(False)

  try:
    nfirst_period = len(data.loc[treated_fp, idname].unique()) if panel \
      else int(np.sum(treated_fp))
  except Exception:
    nfirst_period = len(data.loc[treated_fp, idname].unique()) if panel \
      else int(treated_fp.sum())

  if nfirst_period > 0:
    warnings.warn(
      "Dropped {n} units that were already treated in the first period{anti}.".format(
        n=nfirst_period,
        anti=(f" (accounting for anticipation = {anticipation})" if anticipation > 0 else ""),
      )
    )
    # Drop ONLY the first-period-treated units, by row identity. Matches R `did`
    # 2.5.1 (commit 3246281). The previous `data[gname].isin(glist ∪ {0})` dropped
    # by cohort membership in glist, which -- when there is no never-treated group --
    # also deleted the latest cohort that was deliberately removed from glist above
    # (the `glist < latest_g` trim) so it could serve as a not-yet-treated control.
    # That silently deleted a valid comparison cohort and corrupted ATT(g,t) for the
    # other groups. The `treated_fp` mask removes exactly the already-treated units,
    # nothing else. (R 2.5.0 shipped the buggy membership filter; 2.5.1 reversed it
    # to this row-identity drop -- so this REVERSES the earlier port "N1 fix" that
    # had matched 2.5.0 bug-for-bug.)
    data = data[~treated_fp]
    tlist = np.sort(data[tname].unique())
    glist = np.sort(data[gname].unique())
    glist = glist[glist > 0]
    # P-empty (O3-F2): dropping the first-period-treated units can empty the data
    # (e.g. every unit is treated in period 1). Guard before `tlist[0]` so the
    # clean "No valid groups" error is raised instead of an IndexError.
    if len(tlist) == 0 or len(data) == 0:
      raise ValueError(
        f"No valid groups. The variable in '{gname}' should be expressed as the "
        f"time a unit is first treated (0 if never-treated)."
      )
    fp = tlist[0]
    glist = glist[glist > fp + anticipation]
    # The latest cohort stays in the data as a not-yet-treated control but, when
    # there is still no never-treated group, must remain excluded from glist (it
    # gets no ATT of its own) -- mirroring the nfirst_period == 0 case above.
    if control_group != "nevertreated" and not (data[gname] == 0).any():
      glist = glist[glist < latest_g]

  # idname is validated to be numeric in _validate_inputs (R did v2.5.1).
  true_rep_cross_section = False
  if not panel:
    true_rep_cross_section = True

  # fix_weights="varying" with panel data: use the repeated-cross-section
  # estimators with per-period weights (matches R `did`).
  if fix_weights == "varying" and panel:
    panel = False
    true_rep_cross_section = False

  if panel:
    if allow_unbalanced_panel:
      try:
        n = data[idname].nunique()
      except Exception:
        n = len(pd.unique(data[idname]))
      # Detect imbalance exactly as R did v2.5.1: a panel is balanced iff
      # nrow(data) == n_units * n_periods. If unbalanced, switch to the
      # repeated-cross-section estimators.
      n_periods = data[tname].nunique()
      if len(data) != n * n_periods:
        panel = False
        true_rep_cross_section = False
    else:
      keepers = data.dropna().index
      n = len(data[idname].unique())
      n_keep = len(data.loc[keepers, idname].unique())

      if len(data.loc[keepers]) < len(data):
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
      raise ValueError("All observations dropped due to missing data problems.")
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

  # Warn about the handling of time-varying weights (matches R `did`). R emits this
  # whenever time-varying weights are detected on a panel, regardless of fix_weights
  # (the message points the user at fix_weights to control the behavior) -- so the
  # warning must NOT be gated on `fix_weights is None` (F45-2).
  if weights_name is not None and input_panel:
    try:
      w_tv = data.groupby(idname)['w'].nunique()
      if (w_tv > 1).any():
        warnings.warn(
          "Time-varying weights detected. For balanced panel data, the default "
          "behavior uses the weight from the earlier of the two time periods in each "
          "2x2 comparison (the base period for post-treatment cells). Use the "
          "'fix_weights' argument to control this behavior."
        )
    except Exception:
      pass

  # data.loc[:, ".w"] = data['w']
  if len(glist) == 0:
    raise ValueError(f"No valid groups. The variable in '{gname}' should be expressed as the time a unit is first treated (0 if never-treated).")
  if len(tlist) == 2:
    cband = False
  gsize = data.groupby(data[gname]).size().reset_index(name="count")
  gsize["count"] /= len(tlist)

  # P-reqsize (O3-F3 / C2-F10): R uses length(BMisc::rhs_vars(xformla)) + 5, which
  # counts RHS variables only (0 for ~1) -- NOT the patsy intercept column. The old
  # `n_cov + 5` was 1 too high (e.g. ~1: port reqsize 6 vs R 5), over-warning and
  # over-refusing groups of size 5. Use n_rhs (computed from the formula above).
  reqsize = n_rhs + 5
  gsize = gsize[gsize["count"] < reqsize]

  if len(gsize) > 0:
    gpaste = ",".join(map(str, gsize[gname]))
    warnings.warn(f"Be aware that there are some small groups in your dataset.\n  Check groups: {gpaste}.")

    if 0 in gsize[gname].to_numpy() and control_group == "nevertreated":
      raise ValueError("Never-treated group is too small, try setting control_group='notyettreated'.")
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
    'clustervars': clustervar, 'fix_weights': fix_weights,
    'input_panel': input_panel,
    'xcov_cols': xcov_cols
  }
  return did_params
