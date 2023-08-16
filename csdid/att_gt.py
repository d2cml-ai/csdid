# from aggte import AGGte
from csdid.aggte_fnc.aggte import aggte as agg_te

from csdid.attgt_fnc.preprocess_did import pre_process_did
from csdid.attgt_fnc.compute_att_gt import compute_att_gt

from csdid.utils.mboot import mboot

from csdid.plots.gplot import gplot, splot


import matplotlib.pyplot as plt


import numpy as np, pandas as pd

# class ATTgt(AGGte):
class ATTgt:
  def __init__(self, yname, tname, idname, gname, data, control_group = ['nevertreated', 'notyettreated'], 
  xformla: str = None, panel = True, allow_unbalanced_panel = True, 
  clustervar = None, weights_name = None, anticipation = 0, 
  cband = False, biters = 1000, alp = 0.05
  ):
    dp = pre_process_did(
      yname=yname, tname = tname, idname=idname, gname = gname,
      data = data, control_group=control_group, anticipation=anticipation,
      xformla=xformla, panel=panel, allow_unbalanced_panel=allow_unbalanced_panel, cband=cband, clustervar=None, weights_name=None
    )

    dp['biters'] = biters
    dp['alp'] = alp
    dp['true_repeated_cross_sections'] = dp['true_rep_cross_section']
    dp['cband'] = cband
    self.dp = dp

  def fit(self, est_method = 'dr', base_period = 'varying', bstrap = True):
    # print(self.dp)
    dp = self.dp
    result, inffunc = compute_att_gt(dp)
    att = result['att']
    crit_val, se, V = np.zeros(len(att)), np.zeros(len(att)), np.zeros(len(att))
    if bstrap:
      ref_se = mboot(inffunc.T, dp)
      crit_val, se = ref_se['crit_val'], ref_se['se']
      V = ref_se['V']
    
    ############# aggte input
    group = result['group']
    att = result['att']
    tt = result['year']
    inf_fnc = {'inffunc': inffunc.T}

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
    att_gt = pd.DataFrame(result)
    att_gt = att_gt.drop('c', axis=1)
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
    fig, axes = plt.subplots(nrows=len(group), ncols=1, figsize=(10, 5))  # Adjust the figsize as needed
    handles = []
    labels = []
    for i, group_cat in enumerate(group):
        group_data = results.loc[results['group'] == group_cat]
        title = group_data['grtitle'].unique()[0]
        ax = axes[i]
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
