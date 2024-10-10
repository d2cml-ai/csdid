# from aggte import AGGte
from csdid.aggte_fnc.aggte import aggte as agg_te

from csdid.attgt_fnc.preprocess_did import pre_process_did
from csdid.attgt_fnc.compute_att_gt import compute_att_gt

from csdid.utils.mboot import mboot

from csdid.plots.gplot import gplot, splot


import matplotlib.pyplot as plt


import numpy as np, pandas as pd

# aggte: Función que calcula efectos promedio por grupo o por el tiempo de exposición
# Pre_process_did: Preprocesa los datos para preparar el análisis
# Compute_att_gt: Calcula el ATT para cada grupo y periodo de tiempo
# nboot: Realiza una re-muestreo bootstrap para calcular intervalos de confianza
# gplot y splot: Generan graficos

# class ATTgt(AGGte): 
class ATTgt:
  def __init__(self, yname, tname, idname, gname, data, control_group = ['nevertreated', 'notyettreated'], 
  xformla: str = None, panel = True, allow_unbalanced_panel = True, 
  clustervar = None, weights_name = None, anticipation = 0, 
  cband = False, biters = 1000, alp = 0.05
  ):
    # yname: Nombre de la variable dependiente
    # tname: Nombre de la variable temporal
    # idname: Identificador para cada unidad de análisis
    # gname: Variable que indica los grupos tratados
    # control_group: Grupo que nunca recibió el tratamiento
    # xformula: Formula que describe las covariables
    # panel: Indica si los datos estan formato panel
    # cband: Intervalos de confianza
    # biters: Iteraciones bootstrap
    # alp: Nivel de significancia
    # Otros: Grupo de control, anticipación de tratamiento, otros
    dp = pre_process_did(
      yname=yname, tname = tname, idname=idname, gname = gname,
      data = data, control_group=control_group, anticipation=anticipation,
      xformla=xformla, panel=panel, allow_unbalanced_panel=allow_unbalanced_panel, cband=cband, clustervar=None, weights_name=None
    )
    #dp = pre_process_did(
    #  yname="lemp", tname = "year", idname="countyreal", gname = "first.treat",
    #  data = data, control_group="notyettreated", anticipation='anticipation',
    #  xformla= "lemp~1", panel=True, allow_unbalanced_panel=True, cband=False, clustervar=None, weights_name=None
    #)

    dp['biters'] = biters
    dp['alp'] = alp
    dp['true_repeated_cross_sections'] = dp['true_rep_cross_section']
    dp['cband'] = cband
    dp['panel'] = panel
    self.dp = dp
  
  #Calcula el ATT para cada grupo
  def fit(self, est_method = 'dr', base_period = 'varying', bstrap = True): 
    dp = self.dp
    # Calcula el ATT y la función de influencia
    result, inffunc = compute_att_gt(dp, est_method = est_method, base_period = base_period)
    att = result['att']
    # Calcula los errores estándar y el valor crítico para los intervalos de confianza
    n_len = list(map(len, inffunc))
    crit_val, se, V = (
            1.96,
            np.std(inffunc, axis=1, ddof = 1) / np.sqrt(n_len),
            np.zeros(len(att)),
        )
    # Realiza un remuestreo bootstrap para ajustar los errores estándar
    if bstrap:
      ref_se = mboot(inffunc.T, dp)
      crit_val, se = ref_se['crit_val'], ref_se['se']
      V = ref_se['V']

    return dp