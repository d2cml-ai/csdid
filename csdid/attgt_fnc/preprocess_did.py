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

# Internal column names that must not collide with user data
_RESERVED_NAMES = {'w', 'w1', 'G_m', 'C', 'y_main', 'rowid', 'y0', 'y1', 'dy'}


def _strictly_numeric(series):
  """True for numeric dtypes (including pandas nullable Int64/Float64), but not
  boolean. Uses pandas dtype checks so extension dtypes are handled and do not
  crash ``np.issubdtype`` (which raises on e.g. ``Int64Dtype``)."""
  return is_numeric_dtype(series) and not is_bool_dtype(series)


def _validate_inputs(yname, tname, idname, gname, data, control_group,
                     anticipation, panel, clustervar, weights_name):
  """Validate inputs to pre_process_did, raising clear errors on bad data.

  Returns the normalized ``clustervar`` (a single-element list/tuple is unwrapped
  to a scalar; an empty one becomes ``None``), so callers store a scalar cluster
  variable consistently.
  """

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
  if isinstance(clustervar, (list, tuple)):
    if len(clustervar) > 1:
      raise ValueError(
        "At most one cluster variable (beyond 'idname') is supported. "
        "Please reduce to one."
      )
    clustervar = clustervar[0] if len(clustervar) == 1 else None

  # --- Check columns exist ---
  required = {
    'yname': yname, 'tname': tname, 'idname': idname, 'gname': gname
  }
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
  if panel and clustervar is not None and \
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
  # influence function can be indexed by unit). ---
  if not _strictly_numeric(data[idname]):
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

  # --- Duplicate (id, time) ---
  if panel:
    dups = data.duplicated(subset=[idname, tname], keep=False)
    if dups.any():
      n_dups = dups.sum()
      raise ValueError(
        f"Found {n_dups} duplicate (idname, tname) rows. "
        f"Each unit must appear at most once per time period. "
        f"Check for duplicate entries in your data."
      )

  # --- Treatment must be irreversible (gname time-invariant per unit) ---
  if panel:
    g_tv = data.groupby(idname)[gname].nunique()
    if (g_tv > 1).any():
      raise ValueError(
        "The value of gname (treatment variable) must be the same across all "
        "periods for each particular unit. The treatment must be irreversible."
      )

  # --- anticipation ---
  if not isinstance(anticipation, (int, float)) or anticipation < 0:
    raise ValueError(f"'anticipation' must be a non-negative number, got {anticipation}.")

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

  return clustervar


def pre_process_did(yname, tname, idname, gname, data: pd.DataFrame, 
  control_group = ['nevertreated', 'notyettreated'], 
  anticipation = 0, xformla : str = None,
  panel = True, allow_unbalanced_panel = True, cband = False,
  clustervar = None,  weights_name = None, fix_weights = None
  ) -> dict:

  n, t = data.shape
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

  # Validate inputs early with clear error messages. Returns the normalized
  # clustervar (single-element list unwrapped to a scalar; see V4-U3).
  clustervar = _validate_inputs(yname, tname, idname, gname, data, control_group,
                                anticipation, panel, clustervar, weights_name)

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

  # Pre-check formula variables against the data columns. A bare name in the
  # formula that is not a column otherwise surfaces downstream as a patsy
  # NameError wrapped in "Error processing formula ..." (unclassified); R reports
  # a clean missing-column error naming the variable(s) (V4-U4). Only flag bare
  # identifiers (e.g. `Ghost`) -- transformed/call terms like `np.log(X)` or
  # `C(x)` are left for patsy to resolve.
  try:
    _desc = patsy.ModelDesc.from_formula(xformla)
    _fnames = {f.name() for side in (_desc.lhs_termlist, _desc.rhs_termlist)
               for term in side for f in term.factors}
    _missing_vars = sorted(
      nm for nm in _fnames
      if re.fullmatch(r'[A-Za-z_]\w*', nm) and nm not in data.columns
    )
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


  data = data.dropna()
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
    warnings.warn(f"dropped {ndiff} rows from original data due to missing data")
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
    # Match R `did` (pre_process_did) EXACTLY: keep only rows whose cohort is in
    # {0} ∪ glist. At this point glist already excludes the first-period-treated
    # cohort(s) (via `glist > fp + anticipation`) AND -- in the no-never-treated
    # case -- the last cohort (via `glist < latest_g` above). R's
    #   data <- data[data[, gname] %in% c(0, glist), ]
    # therefore drops the already-treated units *and* the last-treated cohort's
    # units. The latter is essential and was the N1 bug: with
    # control_group="notyettreated", anticipation>0 and no never-treated group,
    # keeping the last cohort left its units in the data, where the (correct)
    # not-yet-treated control cutoff then admitted them as valid controls,
    # contaminating the comparison group and biasing the ATT. The earlier behavior
    # (`data = data[~treated_fp]`) dropped only the first-period rows and so kept
    # the last cohort. When a never-treated group IS present this filter is
    # equivalent to dropping `treated_fp` (glist = all cohorts > fp+anticipation,
    # plus 0), so non-N1 cases are unchanged.
    glist_in = np.append(glist, 0)
    data = data[data[gname].isin(glist_in)]
    tlist = np.sort(data[tname].unique())
    glist = np.sort(data[gname].unique())
    glist = glist[glist > 0]
    fp = tlist[0]
    glist = glist[glist > fp + anticipation]

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

  reqsize = n_cov + 5
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
