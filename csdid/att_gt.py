# from aggte import AGGte
from csdid.aggte_fnc.aggte import aggte as agg_te

from csdid.attgt_fnc.preprocess_did import pre_process_did
from csdid.attgt_fnc.compute_att_gt import compute_att_gt
from csdid.attgt_fnc.compute_att_gt2 import compute_att_gt2
from csdid.attgt_fnc.compute_att_gt_shared import last_pretreatment_index

from csdid.utils.mboot import mboot

from csdid.plots.gplot import gplot, splot


import matplotlib.pyplot as plt

import warnings

import numpy as np, pandas as pd

from scipy.stats import norm

import math
from decimal import Decimal


def _unwrap_len1(x):
  """Unwrap a length-1 container / 0-d numpy array to its scalar element, matching
  R's treatment of length-1 vectors as scalars. Returns ``(scalar, True)`` when an
  unwrap happened, else ``(x, False)``. Multi-element containers are returned as-is
  with the second flag False so the caller can reject them by length. A python str /
  bytes is NOT treated as a container (it is itself a scalar value)."""
  if isinstance(x, np.ndarray):
    if x.ndim == 0:
      return x[()], True
    if x.size == 1:
      return x.reshape(-1)[0], True
    return x, False
  if isinstance(x, (list, tuple)):
    if len(x) == 1:
      return x[0], True
    return x, False
  return x, False


def _validate_logical_scalar(x, name):
  """Mirror R `validate_logical_scalar`: accept only a single logical (python bool /
  np.bool_), unwrapping a length-1 logical vector first (F-C1-1; so ``(False,)`` ->
  False rather than being read as truthy). Reject 'yes'/0/1/2/None/1.0/[True,False]/
  ''/2.5/etc. (C1-1)."""
  x, _ = _unwrap_len1(x)
  if not isinstance(x, (bool, np.bool_)):
    raise ValueError(f"{name} must be a single logical (True or False).")
  return bool(x)


def _validate_positive_whole_number(x, name):
  """Mirror R `validate_positive_whole_number`: a single finite positive whole number
  (C1-3). Accepts python/np ints and whole floats (4.0), unwrapping a length-1 vector
  first. Rejects bool, str, non-finite, non-whole (1.5), and < 1 (0/neg)."""
  x, _ = _unwrap_len1(x)
  # bool is an int subclass in python; R logicals are not whole numbers here.
  if isinstance(x, (bool, np.bool_)):
    raise ValueError(f"{name} must be a single positive whole number.")
  if not isinstance(x, (int, float, np.integer, np.floating, Decimal)):
    raise ValueError(f"{name} must be a single positive whole number.")
  xf = float(x)
  if not math.isfinite(xf) or xf < 1 or xf != round(xf):
    raise ValueError(f"{name} must be a single positive whole number.")
  return int(round(xf))


def _coerce_alp(alp):
  """Generalize the `alp` scalar type check to mirror R `validate_alp` (C1-4, C1-8):
  accept length-1 numeric vectors / 0-d arrays / np.number / Decimal / python
  int|float (unwrapping a length-1 container to its scalar), and reject NaN/non-finite
  (NaN<=0 and NaN>=1 are both False, so the old guard let NaN through -> all-NaN
  inference) as well as length!=1 containers, bool, str, and out-of-range values.
  Returns a python float in (0, 1)."""
  alp, _ = _unwrap_len1(alp)
  # bool would pass the numeric isinstance below but is not a valid significance level.
  if isinstance(alp, (bool, np.bool_)):
    raise ValueError(f"'alp' must be a number strictly between 0 and 1, got {alp}.")
  if not isinstance(alp, (int, float, np.integer, np.floating, Decimal)):
    raise ValueError(f"'alp' must be a number strictly between 0 and 1, got {alp}.")
  alpf = float(alp)
  if not math.isfinite(alpf) or alpf <= 0 or alpf >= 1:
    raise ValueError(f"'alp' must be a number strictly between 0 and 1, got {alp}.")
  return alpf


def _cluster_sums(inffunc, dp):
    """Sum the influence function over the units in each cluster.

    Returns ``(cluster_sums, n)`` where ``cluster_sums`` is ``(n_clusters x n_cells)``
    and ``n`` is the number of units, or ``(None, n)`` when there is no effective
    (beyond-unit) cluster variable. Shared by the cluster-robust analytical SE and
    the clustered analytical Wald pre-test (both match R `did`, bstrap=FALSE).
    """
    clustervars = dp.get('clustervars')
    idname = dp['idname']
    cv = clustervars
    if isinstance(cv, (list, tuple)):
        cv = cv[0] if len(cv) else None
    inffunc = np.asarray(inffunc)
    n = inffunc.shape[1] if inffunc.ndim == 2 else len(inffunc)
    if cv is None or cv == idname:
        return None, n

    data = dp['data']
    tname = dp['tname']
    panel = dp['panel']
    if panel:
        tlist = np.sort(data[tname].unique())
        dta = data[data[tname] == tlist[0]]
    else:
        dta = data.drop_duplicates(subset=[idname])
    if cv not in dta.columns:
        return None, n

    cluster_map = dta[[idname, cv]].drop_duplicates().set_index(idname)[cv]
    unit_ids = dta[idname].unique()
    labels = np.asarray(cluster_map.reindex(unit_ids).values)

    uniq = np.unique(labels)
    codes = np.searchsorted(uniq, labels)

    # cluster_sums[cluster, cell] = sum of IF over units in that cluster
    cluster_sums = np.zeros((len(uniq), inffunc.shape[0]))
    np.add.at(cluster_sums, codes, inffunc.T)
    return cluster_sums, n


def _analytical_cluster_se(inffunc, dp):
    """Cluster-robust analytical SE per (g,t) cell, matching R `did` (bstrap=FALSE).

    For each cell the SE is ``sqrt(sum_c S_c^2) / n`` where ``S_c`` is the sum of the
    influence function over the units in cluster ``c`` and ``n`` is the number of
    units. Returns ``None`` when there is no effective cluster variable (so the
    caller falls back to the i.i.d. analytical SE).
    """
    cluster_sums, n = _cluster_sums(inffunc, dp)
    if cluster_sums is None:
        return None
    return np.sqrt((cluster_sums ** 2).sum(axis=0)) / n


def _wald_pretest(inffunc, att, group, year, dp):
    """Parallel-trends Wald pre-test ``(W, Wpval)``, matching R `did` 2.5.0 `att_gt`.

    Wald test that the pre-treatment ``ATT(g,t)`` (the cells with cohort ``g`` first
    treated strictly after period ``t``) are jointly zero, using the ANALYTIC
    influence-function covariance ``V`` -- the same ``V`` the analytic SEs derive from
    (the pre-test always uses the analytic variance, never the bootstrap one, even
    when ``bstrap=True``). Returns ``(W, Wpval)`` or ``(None, None)`` when not
    computable, with R's guard semantics:

      * clustering BEYOND the unit level together with a clustered bootstrap -> the
        analytic ``V`` ignores between-cluster correlation, so the pre-test is
        suppressed (use the bootstrap confidence intervals instead);
      * no usable pre-treatment cells (every pre cell has zero/NA variance -- e.g. a
        universal-base base cell -- or every cohort is treated from the first period);
      * any NA in the pre covariance, or a singular pre covariance
        (``rcond(preV) <= eps``).

    ``W = n * preatt' solve(preV) preatt`` and
    ``Wpval = round(1 - chi2.cdf(W, q), 5)`` with ``q`` the number of usable pre
    cells -- byte-for-byte R's formula.
    """
    from scipy.stats import chi2

    inffunc = np.asarray(inffunc, dtype=float)
    if inffunc.ndim != 2 or inffunc.size == 0:
        return None, None
    att = np.asarray(att, dtype=float)
    group = np.asarray(group, dtype=float)
    year = np.asarray(year, dtype=float)

    clustervars = dp.get('clustervars')
    cv = (clustervars[0] if isinstance(clustervars, (list, tuple)) and len(clustervars)
          else clustervars)
    extra_cluster = cv is not None and cv != dp.get('idname')
    bstrap = bool(dp.get('bstrap', False))

    # R `wald_invalid`: beyond-unit clustering with a clustered bootstrap -> the
    # analytic variance does not account for between-cluster correlation -> suppress.
    if extra_cluster and bstrap:
        return None, None

    # Analytic V: clustered crossprod when a beyond-unit clustervar is set and the
    # bootstrap is off (R `cluster_analytic`); i.i.d. otherwise.
    if extra_cluster:
        cluster_sums, n = _cluster_sums(inffunc, dp)
        if cluster_sums is None:
            n = inffunc.shape[1]
            V = inffunc @ inffunc.T / n
        else:
            V = cluster_sums.T @ cluster_sums / n
    else:
        n = inffunc.shape[1]
        V = inffunc @ inffunc.T / n

    se = np.sqrt(np.clip(np.diag(V), 0.0, None) / n)
    # Drop pre cells with zero OR NA/inf variance. NaN must be tested explicitly:
    # `NaN <= x` is False, so a NaN-variance cell would otherwise survive this filter,
    # leave NaN in preV, and force the pre-test to bail with None even when other pre
    # cells are valid -- R instead computes the Wald over the valid subset (F45-1).
    zero_na = (se <= np.sqrt(np.finfo(float).eps) * 10) | ~np.isfinite(se)

    pre = np.where(group > year)[0]
    pre = pre[~zero_na[pre]]
    if pre.size == 0:
        # No usable pre-treatment cells -> the Wald pre-test of parallel trends
        # cannot be computed. R `did` warns here rather than silently dropping it
        # (V4-V2); e.g. a two-period (T=2) panel has no pre-treatment periods.
        warnings.warn(
            "No pre-treatment periods available for the Wald pre-test of parallel "
            "trends. The pre-test will not be reported."
        )
        return None, None
    preatt = att[pre]
    preV = V[np.ix_(pre, pre)]
    if np.isnan(preV).any():
        # A-wald-warn (C3-F9): R warns before suppressing (att_gt.R:667) rather than
        # returning None silently.
        warnings.warn(
            "Not returning pre-test Wald statistic due to NA pre-treatment values"
        )
        return None, None
    # R: rcond(preV) <= .Machine$double.eps  (singular pre-covariance)
    try:
        rcond = 1.0 / np.linalg.cond(preV, p=1)
    except Exception:
        rcond = 0.0
    if not np.isfinite(rcond) or rcond <= np.finfo(float).eps:
        # A-wald-warn (C3-F9): R warns before suppressing (att_gt.R:672).
        warnings.warn(
            "Not returning pre-test Wald statistic due to singular covariance matrix"
        )
        return None, None
    try:
        W = float(n * preatt @ np.linalg.solve(preV, preatt))
    except np.linalg.LinAlgError:
        return None, None
    q = int(pre.size)
    Wpval = round(float(1.0 - chi2.cdf(W, q)), 5)
    return W, Wpval


# class ATTgt(AGGte):
class ATTgt:
  def __init__(self, yname, tname, idname, gname, data, control_group = ['nevertreated', 'notyettreated'],
  xformla: str = None, panel = True, allow_unbalanced_panel = False,
  clustervar = None, weights_name = None, anticipation = 0,
  cband = False, biters = 1000, alp = 0.05, compute_inffunc = True,
  fix_weights = None, faster_mode = False,
  print_details = False, pl = False, cores = 1, **kwargs
  ):
    # v7-O1..O4: accept R `att_gt`'s extra formals for call-site compatibility.
    #   print_details : verbosity toggle -- accepted as a no-op (the port is quiet).
    #   pl, cores     : parallel multiplier-bootstrap controls, threaded through to
    #                   mboot in fit() (R defaults pl=FALSE/cores=1 reproduce the
    #                   serial behavior exactly).
    #   **kwargs      : extra args forwarded to a custom est_method in R; with a
    #                   built-in est_method R warns "Extra arguments ... are ignored"
    #                   (att_gt__02). Mirror that warning here instead of raising.
    # A-argalias (C2 / C3-F3/F4): R's formals are `clustervars` (plural) and
    # `weightsname` (no underscore). When an R-faithful call passes those, the port
    # would otherwise drop them via **kwargs and only emit the generic "Extra
    # arguments ... ignored" warning -- silently disabling clustering / weights.
    # Consume the aliases here BEFORE the extra-args warning, preferring the
    # explicitly-set canonical arg and erroring on a genuine conflict (both given
    # with different values).
    if 'clustervars' in kwargs:
      _alias = kwargs.pop('clustervars')
      if clustervar is not None and clustervar != _alias:
        raise ValueError(
          "Both 'clustervar' and its alias 'clustervars' were supplied with "
          "different values; pass only one."
        )
      if clustervar is None:
        clustervar = _alias
    if 'weightsname' in kwargs:
      _alias = kwargs.pop('weightsname')
      if weights_name is not None and weights_name != _alias:
        raise ValueError(
          "Both 'weights_name' and its alias 'weightsname' were supplied with "
          "different values; pass only one."
        )
      if weights_name is None:
        weights_name = _alias

    # A-logical (C1-1): validate logical scalars rather than truthy/bool()-coercing.
    self.print_details = _validate_logical_scalar(print_details, "print_details")
    self.pl = _validate_logical_scalar(pl, "pl")
    panel = _validate_logical_scalar(panel, "panel")
    allow_unbalanced_panel = _validate_logical_scalar(
      allow_unbalanced_panel, "allow_unbalanced_panel")
    cband = _validate_logical_scalar(cband, "cband")
    faster_mode = _validate_logical_scalar(faster_mode, "faster_mode")

    # A-cores (C1-3): validate `cores` as a positive whole number (was int(cores),
    # which silently truncated 1.5, accepted 0/neg, coerced '4'/True).
    self.cores = _validate_positive_whole_number(cores, "cores")

    # A-control-vec (C1-5): R `validate_choice_scalar` requires a length-1
    # control_group. The port's default is the LIST ['nevertreated','notyettreated']
    # (a sentinel meaning "use the first / default"), which preprocess_did reduces to
    # element[0]. Refuse a USER-supplied multi-element control_group that is NOT that
    # exact default sentinel, instead of silently reducing it to element[0]. The
    # default list (and any length-1 container) is left untouched for preprocess_did.
    _default_cg = ['nevertreated', 'notyettreated']
    if isinstance(control_group, (list, tuple)):
      if len(control_group) > 1 and list(control_group) != _default_cg:
        raise ValueError(
          "control_group must be either 'nevertreated' or 'notyettreated'."
        )

    # A-argalias: only warn about genuinely-unknown kwargs (aliases already consumed).
    if kwargs:
      warnings.warn(
        "Extra arguments " + ", ".join(repr(k) for k in kwargs)
        + " are ignored when using a built-in est_method."
      )
    # A-alp-nan (C1-4) / A-np-scalar (C1-8): accept length-1 numeric vectors / 0-d
    # arrays / np.number / Decimal for alp, and reject NaN/non-finite (the old
    # `alp<=0 or alp>=1` guard let NaN through -> norm.ppf(NaN)=NaN -> all-NaN
    # inference).
    alp = _coerce_alp(alp)
    # biters: validate at construction (a single positive whole number). R defers
    # this to the bootstrap branch (C1-9 flags the port as marginally stricter), but
    # validating early is benign and strictly safer -- an invalid biters is a user
    # error regardless of bstrap -- so the port keeps its construction-time guard.
    # biters=1 is valid (the degenerate-bootstrap crash is fixed in mboot, O7-F1).
    biters = _validate_positive_whole_number(biters, "biters")

    # Validate fix_weights (matches R `did`).
    if fix_weights is not None and fix_weights not in ("varying", "base_period", "first_period"):
      raise ValueError(
        "fix_weights must be None or one of 'varying', 'base_period', or 'first_period'."
      )

    # Point-estimates-only mode: skip IF computation. compute_inffunc must be a
    # single logical (matches R `did`); a truthy non-bool (e.g. "yes") was silently
    # bool()-coerced to True (V4-S4), masking a user error.
    # A-inffunc-list (C1-10): accept a length-1 logical vector ([True]/(True,)) by
    # unwrapping before the bool check (R accepts; the port's bare isinstance rejected
    # a list). _validate_logical_scalar does the unwrap-then-validate.
    self.compute_inffunc = _validate_logical_scalar(
      compute_inffunc, "compute_inffunc")
    if not self.compute_inffunc:
      cband = False

    dp = pre_process_did(
      yname=yname, tname=tname, idname=idname, gname=gname,
      data=data, control_group=control_group, anticipation=anticipation,
      xformla=xformla, panel=panel, allow_unbalanced_panel=allow_unbalanced_panel,
      cband=cband, clustervar=clustervar, weights_name=weights_name,
      fix_weights=fix_weights
    )

    dp['biters'] = biters
    dp['alp'] = alp
    dp['true_repeated_cross_sections'] = dp['true_rep_cross_section']
    dp['cband'] = cband
    self.faster_mode = bool(faster_mode)
    self.dp = dp

  def fit(self, est_method = 'dr', base_period = 'varying', bstrap = True):
    dp = self.dp

    # A-logical (C1-1): bstrap is a logical scalar (was bool()-coerced).
    bstrap = _validate_logical_scalar(bstrap, "bstrap")

    # Validate est_method and base_period (matches R `did`).
    if not (callable(est_method) or est_method in ('dr', 'reg', 'ipw')):
      raise ValueError(
        f"'est_method' must be one of 'dr', 'reg', 'ipw', or a callable, got {est_method!r}."
      )
    if base_period not in ('varying', 'universal'):
      raise ValueError(
        f"'base_period' must be 'varying' or 'universal', got {base_period!r}."
      )

    # fix_weights='varying' on panel data switches the engine to the
    # repeated-cross-section estimators (preprocess flips panel->False), whose
    # callable est_method signature differs from the panel one. R refuses this
    # combination up front rather than invoking the callable with the wrong
    # signature (V4-S3). dp['input_panel'] is the user's original panel flag,
    # before the varying-weights flip.
    if (callable(est_method) and dp.get('fix_weights') == 'varying'
        and dp.get('input_panel')):
      raise ValueError(
        'fix_weights = "varying" is not currently supported with a custom '
        '(callable) est_method.'
      )

    # Point-estimates-only mode
    if not self.compute_inffunc:
      bstrap = False

    _compute = compute_att_gt2 if self.faster_mode else compute_att_gt
    result, inffunc = _compute(dp, est_method = est_method, base_period = base_period,
                                     compute_inffunc = self.compute_inffunc)
    att = result['att']
    
    if self.compute_inffunc:
      n_len = list(map(len, inffunc))
      # Analytical i.i.d. SE = sqrt(mean(IF^2)) / sqrt(n): the population (1/n),
      # non-centered influence-function variance, matching R `did`/DRDID and the
      # port's own clustered (`_analytical_cluster_se`) and aggregation
      # (`aggte_fnc/utils.py`) paths. (Previously np.std(..., ddof=1), the sample
      # variance, which inflated every i.i.d. SE by exactly sqrt(n/(n-1)).)
      inffunc_arr = np.asarray(inffunc)
      # Pointwise critical value MP$c = qnorm(1 - alp/2) (matches R `did`). Was a
      # hardcoded 1.96 literal (V4-V1), which ignored alp entirely on the analytic
      # path -- wrong at any alp != 0.05 (and off by ~3.6e-5 even at the default).
      # The bootstrap path below overwrites crit_val from mboot.
      crit_val, se, V = (
              norm.ppf(1 - dp['alp'] / 2),
              np.sqrt(np.mean(inffunc_arr ** 2, axis = 1)) / np.sqrt(n_len),
              np.zeros(len(att)),
          )
      if bstrap:
        # v7-O2/O3: thread pl/cores through to the multiplier bootstrap.
        ref_se = mboot(inffunc.T, dp, pl=self.pl, cores=self.cores)
        crit_val, se = ref_se['crit_val'], ref_se['se']
        V = ref_se['V']
      else:
        # Analytical cluster-robust SEs when a cluster variable is set (matches R).
        cl_se = _analytical_cluster_se(inffunc, dp)
        if cl_se is not None:
          se = cl_se
      # Universal base period: each cohort's base cell has att=0 by construction
      # and an undefined SE; R `did` reports NA there. Identify that cell by its
      # position -- year == the cohort's last pre-treatment period -- NOT by an
      # all-zero influence function, since a degenerate-but-valid estimated cell
      # can also have an all-zero IF and must keep its (zero) SE.
      if base_period == 'universal':
        tl = np.asarray(dp['tlist'])
        ant = dp['anticipation']
        groups = np.asarray(result['group'])
        years = np.asarray(result['year'])
        se = np.array(se, dtype=float)
        for _i in range(len(se)):
          _b = last_pretreatment_index(groups[_i], tl, ant)
          if _b is not None and years[_i] == tl[_b]:
            se[_i] = np.nan
      inf_fnc = {'inffunc': inffunc.T}
    else:
      # No influence functions: no SEs, no bootstrap
      att = np.array(att)
      se = np.full(len(att), np.nan)
      crit_val = np.nan
      V = np.zeros(len(att))
      inf_fnc = None
    
    ############# aggte input
    group = result['group']
    att = result['att']
    tt = result['year']

    dp['bstrap'] = bstrap
    dp['est_method'] = est_method
    dp['base_period'] = base_period
    self.dp = dp
    n = dp['n']

    # Parallel-trends Wald pre-test (R `did` 2.5.0 att_gt $W / $Wpval). Needs the
    # influence functions; in point-estimates-only mode there is nothing to test.
    if self.compute_inffunc:
      W, Wpval = _wald_pretest(inffunc_arr, result['att'],
                               result['group'], result['year'], dp)
    else:
      W, Wpval = None, None

    mp = {
      'group': group, 'att': att, 't': tt,
      'DIDparams': dp, 'inffunc': inf_fnc,
      'n': n, 'W': W, 'Wpval': Wpval
    }
    self.MP = mp

    att = np.array(att)
    se = np.array(se)
    cband_lower = att - crit_val * se
    cband_upper = att + crit_val * se
    sig = (cband_upper < 0) | (cband_lower > 0)
    sig[np.isnan(sig)] = False
    sig_text = np.where(sig, "*", "")

    result.update(
      {
        'se': se, 'l_se': cband_lower,
        'c': crit_val,
        'u_se': cband_upper, 'sig': sig_text
       })
    
    self.results = result
    # Backward-compatible alias for the legacy trailing-space key. The canonical
    # key is now 'post'; 'post ' is kept so older code does not hit a KeyError.
    result['post '] = result['post']

    rst = result
    did_object = {
      'group': mp['group'],
      't': mp['t'],
      'att': rst['att'],
      'se': rst['se'],
      'c': rst['c'],
    }
    self.did_object = did_object
    return self
  def summ_attgt(self, n = 4):
    result = self.results
    cols = ['group', 'year', 'att', 'post', 'se', 'l_se', 'u_se', 'sig']
    att_gt = pd.DataFrame({c: result[c] for c in cols})
    name_attgt_df = ['Group', 'Time', 'ATT(g, t)', 'Post', "Std. Error", "[95% Pointwise", 'Conf. Band]', '']
    att_gt.columns = name_attgt_df
    att_gt = att_gt.round(n)
    self.summary2 = att_gt
    return self

  def aggte(
    self, 
    typec         = "group",
    balance_e     = None,
    min_e         = float('-inf'),
    max_e         = float('inf'),
    na_rm         = False,
    bstrap        = None,
    biters        = None,
    cband         = None,
    alp           = None,
    clustervars   = None,
    ):
    if not self.compute_inffunc:
      raise ValueError(
        "Cannot run aggte() when compute_inffunc=False. "
        "Re-run ATTgt with compute_inffunc=True to compute influence functions."
      )
    mp = self.MP
    did_object = self.did_object
    
    did_object.update({
      'type': typec
      }
    )


    atte = agg_te(
      mp, typec=typec, balance_e=balance_e, 
      min_e=min_e, max_e=max_e, na_rm=na_rm, bstrap=bstrap, 
      biters=biters, cband=cband, alp=alp, clustervars=clustervars
    )

    self.atte = atte
    return self
  def plot_attgt(self, ylim=None, 
                xlab=None, 
                ylab=None, 
                title="Group",
                xgap=1, 
                ncol=1, 
                legend=True, 
                group=None, 
                ref_line=0,
                theming=True, 
                grtitle="Group"
                ):

    did_object = self.did_object

    grp = did_object['group']
    t_i = did_object['t']

    G = len(np.unique(grp))
    Y = len(np.unique(t_i))
    g = np.unique(grp)[np.argsort(np.unique(grp))].astype(int)
    y = np.unique(t_i)

    results = pd.DataFrame({'year': np.tile(y, G)})
    results['group'] = np.repeat(g, Y)
    results['grtitle'] = grtitle + ' ' + results['group'].astype(str)
    results['att'] = did_object['att']
    results['att_se'] = did_object['se']
    results['post'] = np.where(results['year'] >= grp, 1, 0)
    results['year'] = results['year']
    results['c'] = did_object['c']

    self.results_plot_df_attgt = results

    if group is None:
      group = g
      if any(group not in g for group in group):
        raise ValueError("Some of the specified groups do not exist in the data. Reporting all available groups.")


    legend_1 = False    # for multiple subplots, legend outside 
    # squeeze=False keeps axes as a 2-D array even for a single group, so axes[i, 0]
    # works uniformly (plt.subplots(nrows=1) would otherwise return a bare Axes).
    fig, axes = plt.subplots(nrows=len(group), ncols=1, figsize=(10, 5), squeeze=False)  # Adjust the figsize as needed
    handles = []
    labels = []
    for i, group_cat in enumerate(group):
        group_data = results.loc[results['group'] == group_cat]
        title = group_data['grtitle'].unique()[0]
        ax = axes[i, 0]
        ax = gplot(group_data, ax, ylim, xlab, ylab, title, xgap, legend_1, ref_line, theming)
    plt.tight_layout()
    if legend is True:
        handles_ax, labels_ax = ax.get_legend_handles_labels()
        handles.extend(handles_ax)
        labels.extend(labels_ax)
        fig.legend(handles, labels, loc='lower center', fontsize='small', bbox_to_anchor=(0.545, -0.075), ncol=2)
    
    plt.show()
    return fig 

  def plot_aggte(self, ylim=None, 
                   xlab=None, 
                   ylab=None, 
                   title="", 
                   xgap=1, 
                   legend=True, 
                   ref_line=0, 
                   theming=True,
                   **kwargs):

    did_object = self.atte

    post_treat = 1 * (np.asarray(did_object["egt"]).astype(int) >= 0)
    
    results = {
        "year": list(map(int, did_object["egt"])),
        "att": did_object["att_egt"],
        "att_se": did_object["se_egt"][0],
        "post": post_treat
    }
    
    results = pd.DataFrame(results)
    self.results_plot_df_aggte = results
    
    if did_object['crit_val_egt'] is None:
        # v6-F-PLOTC: respect alp rather than hardcoding qnorm(0.975) ~ 1.96.
        # R `did` uses qnorm(1 - alp/2) for the pointwise band crit value.
        results['c'] = norm.ppf(1 - did_object['DIDparams']['alp'] / 2)
    else:
        results['c'] = did_object['crit_val_egt']

    if title == "":
        title = "Average Effect by Group" if\
          did_object["type"] == "group" else\
            "Average Effect by Length of Exposure"


    if did_object["type"] == "group":
        fig, ax = plt.subplots(figsize=(10, 5))
        p = splot(results, ax, ylim, xlab, ylab, title, legend, ref_line, theming)
        plt.tight_layout()
        plt.show()

    else:
        fig, ax = plt.subplots(figsize=(10, 5))
        p = gplot(results, ax, ylim, xlab, ylab, title, xgap, legend, ref_line, theming)
        plt.tight_layout()
        plt.show() 
        
    return p
