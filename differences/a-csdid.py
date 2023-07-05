# csdid.py
# pip install csdid
import pandas as pd
from memory_profiler import profile
from csdid.att_gt import ATTgt
import time
data = pd.read_csv("https://raw.githubusercontent.com/d2cml-ai/csdid/function-aggte/data/mpdta.csv")

start_time = time.time()
@profile
def csdid_p():
  out = ATTgt(yname = "lemp",
              gname = "first.treat",
              idname = "countyreal",
              tname = "year",
              xformla = f"lemp~1",
              data = data,
              )
  out.fit(est_method = 'dr')

if __name__ == '__main__':
  csdid_p()
  end_time = time.time()
  count = end_time - start_time
  print(f'Execution time (csdid): {count}')