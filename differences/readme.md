## Coments:

- `csdid` has the same syntax as `R`'s `did` package, making it easy to get started and achieve the same structure in the output.
- `Difference` has poor documentation and is not very intuitive for new databases.

## Memory Profiler:

- In order to compare memory usage and execution speed between the "csdid" package and "difference," we will use `Memory Profiler`, which displays line-by-line memory usage and also calculates the execution time of Python code.

- To view this information, the corresponding codes will be executed in two different Python files.

### `csdid`

Data: The original R package's "mpdata" database is being used, which contains 2500 rows.

```py
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
```

```sh
# python csdid.py
```

```
Filename: .\a-csdid.py

Line #    Mem usage    Increment  Occurrences   Line Contents
=============================================================
    10    161.9 MiB    161.9 MiB           1   @profile
    11                                         def csdid_p():
    12    162.5 MiB      0.6 MiB           2     out = ATTgt(yname = "lemp",
    13    161.9 MiB      0.0 MiB           1                 gname = "first.treat",
    14    161.9 MiB      0.0 MiB           1                 idname = "countyreal",
    15    161.9 MiB      0.0 MiB           1                 tname = "year",
    16    161.9 MiB      0.0 MiB           1                 xformla = f"lemp~1",
    17    161.9 MiB      0.0 MiB           1                 data = data,
    18                                                       )
    19    165.1 MiB      2.6 MiB           1     out.fit(est_method = 'dr')


Execution time (csdid): 0.899705171585083
```

### `Difference`

- Data:
  - `simulate_data(nentity=313)`: 2504 rows
  - We were unable to replicate the structure of `difference` using the R examples. However, we attempted to simulate a similar context.

```py
# difference.py
# pip install difference
from differences import ATTgt, simulate_data
from memory_profiler import profile
import time

# almost 2500 rows
df = simulate_data(nentity=313)


start_time = time.time()

@profile
def diffe_rence():
	att_gt = ATTgt(data=df, cohort_name='cohort')
	att_gt.fit(formula='y')

if __name__ == '__main__':
	diffe_rence()
	end_time = time.time()
	count = end_time - start_time
	print(f'Execution time (difference): {count}')
```

```sh
python b_difference.py
```

```
Computing ATTgt [workers=1]   100%|████████████████████| 21/21 [00:00<00:00, 41.48it/s]
Filename: .\b_difference.py

Line #    Mem usage    Increment  Occurrences   Line Contents
=============================================================
    12    163.0 MiB    163.0 MiB           1   @profile
    13                                         def diffe_rence():
    14    163.6 MiB      0.6 MiB           1    att_gt = ATTgt(data=df, cohort_name='cohort')
    15    165.1 MiB      1.5 MiB           1    att_gt.fit(formula='y')

Execution time (difference): 0.8982996940612793
```

### Result

Output interpretation (Memory):

> The first column represents the line number of the code that has been profiled, the second column (Mem usage) the memory usage of the Python interpreter after that line has been executed. The third column (Increment) represents the difference in memory of the current line with respect to the last one. The last column (Line Contents) prints the code that has been profiled.

In terms of **memory usage**, `csdid` uses 161.9 MiB before the execution of the main code, 0.6 MiB to create the class, and 2.6 MiB for estimation. At the end of the execution, the total memory usage is 165.1 MiB.
`differences` uses 163 MiB before the execution of the main code, 0.6 MiB to create the class, and 1.5 MiB for estimation. The total memory usage at the end is also 165.1 MiB.

Time:

- `differences` takes 0.89830 seconds to execute the code, while `csdid` takes 0.8997 seconds.

> The differences in memory usage and execution speed between both packages are similar. In terms of memory, there is no difference, while in terms of execution speed, "difference" is 0.0014 seconds faster.

## PySpark

Currently, the package does not support a PySpark object. Part of the error is due to the `patsy` package, which does not support PySpark and is included in the data preparation function to be estimated.

[See error](https://colab.research.google.com/drive/1EWbvYxXbRaofSeyVVrOO80XWP5dgvgMC?usp=sharing)

```
---------------------------------------------------------------------------
PandasNotImplementedError                 Traceback (most recent call last)
<ipython-input-10-b05bdaddb778> in <cell line: 1>()
----> 1 out = ATTgt(yname = "lemp",
      2               gname = "first.treat",
      3               idname = "countyreal",
      4               tname = "year",
      5               xformla = f"lemp~1",

9 frames
/usr/local/lib/python3.10/dist-packages/pyspark/pandas/missing/__init__.py in unsupported_function(*args, **kwargs)
     21 def unsupported_function(class_name, method_name, deprecated=False, reason=""):
     22     def unsupported_function(*args, **kwargs):
---> 23         raise PandasNotImplementedError(
     24             class_name=class_name, method_name=method_name, reason=reason
     25         )

PandasNotImplementedError: The method `pd.Series.__iter__()` is not implemented. If you want to collect your data as an NumPy array, use
'to_numpy()' instead.
```
