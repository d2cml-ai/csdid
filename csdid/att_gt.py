# from aggte import AGGte
from csdid.aggte_fnc.aggte import aggte as agg_te

from csdid.attgt_fnc.preprocess_did import pre_process_did
from csdid.attgt_fnc.compute_att_gt import compute_att_gt
from csdid.attgt_fnc.compute_att_gt2 import compute_att_gt2

from csdid.utils.mboot import mboot

from csdid.plots.gplot import gplot, splot


import matplotlib.pyplot as plt

import warnings

import numpy as np, pandas as pd


def _analytical_cluster_se(inffunc, dp):
    """Cluster-robust analytical SE per (g,t) cell, matching R `did` (bstrap=FALSE).

    For each cell the SE is ``sqrt(sum_c S_c^2) / n`` where ``S_c`` is the sum of the
    influence function over the units in cluster ``c`` and ``n`` is the number of
    units. Returns ``None`` when there is no effective cluster variable (so the
    caller falls back to the i.i.d. analytical SE).
    """
    clustervars = dp.get('clustervars')
    idname = dp['idname']
    cv = clustervars
    if isinstance(cv, (list, tuple)):
        cv = cv[0] if len(cv) else None
    if cv is None or cv == idname:
        return None

    data = dp['data']
    tname = dp['tname']
    panel = dp['panel']
    if panel:
        tlist = np.sort(data[tname].unique())
        dta = data[data[tname] == tlist[0]]
    else:
        dta = data.drop_duplicates(subset=[idname])
    if cv not in dta.columns:
        return None

    cluster_map = dta[[idname, cv]].drop_duplicates().set_index(idname)[cv]
    unit_ids = dta[idname].unique()
    labels = np.asarray(cluster_map.reindex(unit_ids).values)

    inffunc = np.asarray(inffunc)
    n = inffunc.shape[1]
    uniq = np.unique(labels)
    codes = np.searchsorted(uniq, labels)

    # cluster_sums[cell, cluster] = sum of IF over units in that cluster
    cluster_sums = np.zeros((len(uniq), inffunc.shape[0]))
    np.add.at(cluster_sums, codes, inffunc.T)
    return np.sqrt((cluster_sums ** 2).sum(axis=0)) / n


# class ATTgt(AGGte):
class ATTgt:
  def __init__(self, yname, tname, idname, gname, data, control_group = ['nevertreated', 'notyettreated'], 
  xformla: str = None, panel = True, allow_unbalanced_panel = True, 
  clustervar = None, weights_name = None, anticipation = 0, 
  cband = False, biters = 1000, alp = 0.05, compute_inffunc = True,
  fix_weights = None, faster_mode = False
  ):
    # Validate alp and biters
    if not isinstance(alp, (int, float)) or alp <= 0 or alp >= 1:
      raise ValueError(f"'alp' must be a number strictly between 0 and 1, got {alp}.")
    if not isinstance(biters, (int, float)) or biters < 1 or biters != int(biters):
      raise ValueError(f"'biters' must be a positive integer, got {biters}.")
    biters = int(biters)

    # Validate fix_weights (matches R `did`).
    if fix_weights is not None and fix_weights not in ("varying", "base_period", "first_period"):
      raise ValueError(
        "fix_weights must be None or one of 'varying', 'base_period', or 'first_period'."
      )

    # Point-estimates-only mode: skip IF computation
    self.compute_inffunc = bool(compute_inffunc)
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

    # Validate est_method and base_period (matches R `did`).
    if not (callable(est_method) or est_method in ('dr', 'reg', 'ipw')):
      raise ValueError(
        f"'est_method' must be one of 'dr', 'reg', 'ipw', or a callable, got {est_method!r}."
      )
    if base_period not in ('varying', 'universal'):
      raise ValueError(
        f"'base_period' must be 'varying' or 'universal', got {base_period!r}."
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
      crit_val, se, V = (
              1.96,
              np.std(inffunc, axis=1, ddof = 1) / np.sqrt(n_len),
              np.zeros(len(att)),
          )
      if bstrap:
        ref_se = mboot(inffunc.T, dp)
        crit_val, se = ref_se['crit_val'], ref_se['se']
        V = ref_se['V']
      else:
        # Analytical cluster-robust SEs when a cluster variable is set (matches R).
        cl_se = _analytical_cluster_se(inffunc, dp)
        if cl_se is not None:
          se = cl_se
      # Universal base period: the base cell has att=0 by construction and an
      # all-zero influence function, so its SE is undefined. R `did` reports NA
      # there, so set those SEs to NaN (instead of a misleading 0) to match.
      if base_period == 'universal':
        base_mask = np.all(np.asarray(inffunc) == 0, axis=1)
        se = np.array(se, dtype=float)
        se[base_mask] = np.nan
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

    mp = {
      'group': group, 'att': att, 't': tt,
      'DIDparams': dp, 'inffunc': inf_fnc, 
      'n': n
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
        results['c'] = abs(norm.ppf(0.025))
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
