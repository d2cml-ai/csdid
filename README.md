# Difference in Difference in Python


The **csdid** package contains tools for computing average treatment
effect parameters in a Difference-in-Differences setup allowing for

- More than two time periods

- Variation in treatment timing (i.e., units can become treated at
  different points in time)

- Treatment effect heterogeneity (i.e, the effect of participating in
  the treatment can vary across units and exhibit potentially complex
  dynamics, selection into treatment, or time effects)

- The parallel trends assumption holds only after conditioning on
  covariates

The main parameters are **group-time average treatment effects**. These
are the average treatment effect for a particular group (group is
defined by treatment timing) in a particular time period. These
parameters are a natural generalization of the average treatment effect
on the treated (ATT) which is identified in the textbook case with two
periods and two groups to the case with multiple periods.

Group-time average treatment effects are also natural building blocks
for more aggregated treatment effect parameters such as overall
treatment effects or event-study-type estimands.

## Getting Started

There has been some recent work on DiD with multiple time periods. The
**did** package implements the framework put forward in

- [Callaway, Brantly and Pedro H.C. Sant’Anna.
  “Difference-in-Differences with Multiple Time Periods.” Journal of
  Econometrics, Vol. 225, No. 2, pp. 200-230,
  2021.](https://doi.org/10.1016/j.jeconom.2020.12.001) or
  \[arXiv\](https://arxiv.org/abs/1803.09015

This project is based on the original [did R
package](https://github.com/bcallaway11/did).

## Instalation

You can install **csdid** from `pypi` with:

    pip install csdid

or via github:

    pip install git+https://github.com/d2cml-ai/csdid/

### Dependencies

Additionally, I have created an additional library called `drdid`, which
can be installed via GitHub.

    pip install git+https://github.com/d2cml-ai/DRDID

## Basic Example

The following is a simplified example of the effect of states increasing
their minimum wages on county-level teen employment rates which comes
from [Callaway and Sant’Anna
(2021)](https://authors.elsevier.com/a/1cFzc15Dji4pnC).

- [More detailed examples are also
  available](https://bcallaway11.github.io/did/articles/did-basics.html)

A subset of the data is available in the package and can be loaded by

``` python
from csdid.att_gt import ATTgt
import pandas as pd
data = pd.read_csv("https://raw.githubusercontent.com/d2cml-ai/csdid/function-aggte/data/mpdta.csv")
```

The dataset contains 500 observations of county-level teen employment
rates from 2003-2007. Some states are first treated in 2004, some in
2006, and some in 2007 (see the paper for more details). The important
variables in the dataset are

- **lemp** This is the log of county-level teen employment. It is the
  outcome variable

- **first.treat** This is the period when a state first increases its
  minimum wage. It can be 2004, 2006, or 2007. It is the variable that
  defines *group* in this application

- **year** This is the year and is the *time* variable

- **countyreal** This is an id number for each county and provides the
  individual identifier in this panel data context

To estimate group-time average treatment effects, use the
**ATTgt().fit()** method

``` python
out = ATTgt(yname = "lemp",
              gname = "first.treat",
              idname = "countyreal",
              tname = "year",
              xformla = f"lemp~1",
              data = data,
              ).fit(est_method = 'dr')
```

Summary table

``` python
out.summ_attgt().summary2
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
&#10;    .dataframe tbody tr th {
        vertical-align: top;
    }
&#10;    .dataframe thead th {
        text-align: right;
    }
</style>

|     | Group | Time | ATT(g, t) | Post | Std. Error | \[95% Pointwise | Conf. Band\] |     |
|-----|-------|------|-----------|------|------------|-----------------|--------------|-----|
| 0   | 2004  | 2004 | -0.0105   | 1    | 0.0257     | -0.0809         | 0.0599       |     |
| 1   | 2004  | 2005 | -0.0704   | 1    | 0.0323     | -0.1589         | 0.0181       |     |
| 2   | 2004  | 2006 | -0.1373   | 1    | 0.0393     | -0.2449         | -0.0296      | \*  |
| 3   | 2004  | 2007 | -0.1008   | 1    | 0.0349     | -0.1963         | -0.0053      | \*  |
| 4   | 2006  | 2004 | 0.0065    | 0    | 0.0225     | -0.0552         | 0.0682       |     |
| 5   | 2006  | 2005 | -0.0028   | 0    | 0.0185     | -0.0533         | 0.0478       |     |
| 6   | 2006  | 2006 | -0.0046   | 1    | 0.0188     | -0.0561         | 0.0469       |     |
| 7   | 2006  | 2007 | -0.0412   | 1    | 0.0200     | -0.0961         | 0.0136       |     |
| 8   | 2007  | 2004 | 0.0305    | 0    | 0.0155     | -0.0119         | 0.0729       |     |
| 9   | 2007  | 2005 | -0.0027   | 0    | 0.0166     | -0.0481         | 0.0427       |     |
| 10  | 2007  | 2006 | -0.0311   | 0    | 0.0181     | -0.0805         | 0.0184       |     |
| 11  | 2007  | 2007 | -0.0261   | 1    | 0.0175     | -0.0739         | 0.0218       |     |

</div>

In the graphs, a semicolon `;` should be added to prevent printing the
class and the graph information.

``` python
out.plot_attgt();
```

![](README_files/figure-commonmark/cell-5-output-1.png)

``` python
out.aggte(typec='calendar');
```



    Overall summary of ATT's based on calendar time aggregation:
        ATT Std. Error [95.0%  Conf. Int.]  
    -0.0417     0.0165 -0.074      -0.0094 *


    Time Effects (calendar):
       Time  Estimate  Std. Error  [95.0% Simult.   Conf. Band   
    0  2004   -0.0105      0.0242          -0.0580      0.0370   
    1  2005   -0.0704      0.0313          -0.1318     -0.0090  *
    2  2006   -0.0488      0.0210          -0.0900     -0.0077  *
    3  2007   -0.0371      0.0143          -0.0650     -0.0091  *
    ---
    Signif. codes: `*' confidence band does not cover 0
    Control Group:  Never Treated , 
    Anticipation Periods:  0
    Estimation Method:  Doubly Robust

``` python
out.plot_aggte();
```

![](README_files/figure-commonmark/cell-7-output-1.png)

## Additional options

`ATTgt(...)` and `.fit(...)` support several options for inference,
performance, and weighting (all matching the R `did` package):

- **Control group**: `control_group="nevertreated"` (default) or
  `"notyettreated"`.
- **Estimation method**: `.fit(est_method="dr")` (doubly robust,
  default), `"reg"` (outcome regression), `"ipw"` (inverse-probability
  weighting), or a Python callable.
- **Base period**: `.fit(base_period="varying")` (default) or
  `"universal"`.
- **Clustered inference**: pass `clustervar="..."` for cluster-robust
  standard errors (analytical when `bstrap=False`, or via the multiplier
  bootstrap).
- **Point estimates only**: `compute_inffunc=False` skips
  influence-function and standard-error computation for faster, lighter
  exploratory runs.
- **Time-varying weights** (`weights_name`): `fix_weights` controls
  which period’s weights each 2×2 comparison uses — `None` (earlier
  period, default), `"base_period"`, `"first_period"`, or `"varying"`.
- **Faster computation**: `faster_mode=True` returns identical results
  to the default path but is ~3× faster by precomputing the covariate
  design once per period and pre-pivoting outcomes/weights.

``` python
# ~3x faster, identical results; cluster-robust SEs; point estimates only:
out = ATTgt(yname="lemp", gname="first.treat", idname="countyreal", tname="year",
            xformla="lemp~1", data=data, clustervar="countyreal",
            faster_mode=True).fit(est_method="dr")
```
