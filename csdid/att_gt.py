# from aggte import AGGte
from csdid.aggte_fnc.aggte import aggte as agg_te

from csdid.attgt_fnc.preprocess_did import pre_process_did
from csdid.attgt_fnc.compute_att_gt import compute_att_gt

from csdid.utils.mboot import mboot

import numpy as np

# class ATTgt(AGGte):
class ATTgt:
  def __init__(self, yname, tname, idname, gname, data, control_group = ['nevertreated', 'notyertreated'], 
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

    MP = {
      'group': group, 'att': att, 't': tt,
      'DIDparams': dp, 'inffunc': inf_fnc, 
      'n': n
    }
    self.MP = MP


    cband_lower = att - crit_val * se
    cband_upper = att + crit_val * se
    sig = (cband_upper < 0) | (cband_lower > 0)
    sig[np.isnan(sig)] = False
    sig_text = np.where(sig, "*", "")

    result.update(
      {
        'std': se, 'l_se': cband_lower,
        'u_se': cband_upper, 'sig': sig_text
       })

    self.results = result
    return self
  def summ_attgt(self, n = 4):
    result = self.results
    att_gt = pd.DataFrame(result)
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
    agg_te(
      mp, typec=typec, balance_e=balance_e, 
      min_e=min_e, max_e=max_e, na_rm=na_rm, bstrap=bstrap, 
      biters=biters, cband=cband, alp=alp, clustervars=clustervars
    )
    return self
    
  

# # print(aggte(b.MP, typec='dynamic'))
# # data = pd.read_csv(dt['simdata'])

# # yname = "Y"
# # tname = "period"
# # idname = "id"
# # gname = "G"
# # data = data
# # xformla = "Y~1"


# a = ATTgt(yname, idname, gname, data, xformla=xformla).fit().summ_attgt()
# mp = a.MP


# b = aggte(mp, typec='simple')
# # print(b)

# print(b)

# print(a)