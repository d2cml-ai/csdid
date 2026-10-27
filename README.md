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
| 0   | 2004  | 2004 | -0.0105   | 1    | 0.0235     | -0.0566         | 0.0356       |     |
| 1   | 2004  | 2005 | -0.0704   | 1    | 0.0335     | -0.1360         | -0.0048      | \*  |
| 2   | 2004  | 2006 | -0.1373   | 1    | 0.0417     | -0.2189         | -0.0556      | \*  |
| 3   | 2004  | 2007 | -0.1008   | 1    | 0.0375     | -0.1743         | -0.0273      | \*  |
| 4   | 2006  | 2004 | 0.0065    | 0    | 0.0242     | -0.0409         | 0.0540       |     |
| 5   | 2006  | 2005 | -0.0028   | 0    | 0.0197     | -0.0414         | 0.0359       |     |
| 6   | 2006  | 2006 | -0.0046   | 1    | 0.0190     | -0.0417         | 0.0326       |     |
| 7   | 2006  | 2007 | -0.0412   | 1    | 0.0208     | -0.0819         | -0.0005      | \*  |
| 8   | 2007  | 2004 | 0.0305    | 0    | 0.0155     | 0.0001          | 0.0609       | \*  |
| 9   | 2007  | 2005 | -0.0027   | 0    | 0.0177     | -0.0374         | 0.0320       |     |
| 10  | 2007  | 2006 | -0.0311   | 0    | 0.0187     | -0.0677         | 0.0055       |     |
| 11  | 2007  | 2007 | -0.0261   | 1    | 0.0169     | -0.0591         | 0.0070       |     |

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
        ATT Std. Error  [95.0%  Conf. Int.]  
    -0.0417     0.0177 -0.0764       -0.007 *


    Time Effects (calendar):
       Time  Estimate  Std. Error  [95.0% Simult.   Conf. Band   
    0  2004   -0.0105      0.0235          -0.0566      0.0356   
    1  2005   -0.0704      0.0326          -0.1344     -0.0065  *
    2  2006   -0.0488      0.0204          -0.0889     -0.0088  *
    3  2007   -0.0371      0.0145          -0.0654     -0.0087  *
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
performance, and weighting:

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
- **Faster computation**: `faster_mode=True` is the default starting in
  0.4.3. It precomputes panel outcomes, weights, and covariates by
  period. Pass `faster_mode=False` to use the standard computation path.
  Both modes use the same estimators and inference; speed gains depend
  on the data and options.

``` python
# Cluster-robust SEs using the default faster computation path:
out = ATTgt(yname="lemp", gname="first.treat", idname="countyreal", tname="year",
            xformla="lemp~1", data=data,
            clustervar="countyreal").fit(est_method="dr")
```
