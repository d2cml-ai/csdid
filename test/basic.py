

import yaml, pandas as pd

with open('configs/data.yml') as f:
  dt = yaml.safe_load(f)

data = pd.read_csv(dt['mpdata'])


yname = "lemp"
gname = "first.treat"
idname = "countyreal"
tname = "year"
xformla = f"lemp~1"

from csdid.attgt_fnc import preprocess_did
from csdid.attgt_fnc import compute_att_gt

dp = preprocess_did(yname, tname, idname, gname, data = data, xformla=xformla)


# data = mpdta
# print(data)
# print(tname)

# from csdid.att_gt import ATTgt

# b = ATTgt(yname, tname, idname, gname, data = data, xformla=xformla).fit()
# c = b.summ_attgt(n = 12).summary2

# # print(dir(b))
# # print(b.MP)
# # print(b.dp)
# # print(b.results)

# # print(c)

# c = b.aggte(balance_e=1)
# print(c)
# b.sdplot()
# b.dplto()
# algo()