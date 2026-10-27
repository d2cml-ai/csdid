# Difference in Difference in Python


[![PyPI
version](https://img.shields.io/pypi/v/csdid.svg?color=blue)](https://pypi.org/project/csdid/)
<a href="https://pepy.tech/project/csdid"><img alt="Downloads" src="https://static.pepy.tech/personalized-badge/csdid?period=total&amp;units=international_system&amp;left_color=blue&amp;right_color=grey&amp;left_text=Downloads" /></a>
[![Last
commit](https://img.shields.io/github/last-commit/d2cml-ai/csdid.svg)](https://github.com/d2cml-ai/csdid/commits/main)
[![GitHub
stars](https://img.shields.io/github/stars/d2cml-ai/csdid.svg?style=social)](https://github.com/d2cml-ai/csdid/stargazers)
[![GitHub
issues](https://img.shields.io/github/issues/d2cml-ai/csdid.svg)](https://github.com/d2cml-ai/csdid/issues)
[![License](https://img.shields.io/github/license/d2cml-ai/csdid.svg)](https://github.com/d2cml-ai/csdid/blob/main/LICENSE)

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
**csdid** package implements the framework put forward in

- [Callaway, Brantly and Pedro H.C. Sant’Anna.
  “Difference-in-Differences with Multiple Time Periods.” Journal of
  Econometrics, Vol. 225, No. 2, pp. 200-230,
  2021.](https://doi.org/10.1016/j.jeconom.2020.12.001) or
  [arXiv](https://arxiv.org/abs/1803.09015).

This project is based on the original [did R
package](https://github.com/bcallaway11/did).

## Installation

You can install **csdid** from `pypi` with:

    pip install csdid

or via github:

    pip install git+https://github.com/d2cml-ai/csdid/

### Dependencies

The `drdid` dependency is installed automatically with `csdid`. CSDID
0.4.4 requires the reviewed DRDID 1.1.7 release. The Python backend is
maintained in [DRDIDpy](https://github.com/d2cml-ai/DRDIDpy);
development checkouts use the exact peer revision recorded in
[`python_drdid.json`](csdid/test_csdid/joint_sync/python_drdid.json).
Update both packages together using their matching reviewed versions.

## Basic Example

The following is a simplified example of the effect of states increasing
their minimum wages on county-level teen employment rates which comes
from [Callaway and Sant’Anna
(2021)](https://authors.elsevier.com/a/1cFzc15Dji4pnC).

- [More detailed examples are also
  available](https://bcallaway11.github.io/did/articles/did-basics.html)

A subset of the data is available at [data/mpdta.csv](data/mpdta.csv).
Run the following examples from the repository root:

``` python
from csdid.att_gt import ATTgt
import numpy as np
import pandas as pd
data = pd.read_csv("data/mpdta.csv")
np.random.seed(12345)
```

The dataset contains 2,500 observations for 500 counties, with teen
employment rates from 2003-2007. Some states are first treated in 2004,
some in 2006, and some in 2007 (see the paper for more details). The
important variables in the dataset are

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
| 0   | 2004  | 2004 | -0.0105   | 1    | 0.0243     | -0.0581         | 0.0371       |     |
| 1   | 2004  | 2005 | -0.0704   | 1    | 0.0348     | -0.1387         | -0.0022      | \*  |
| 2   | 2004  | 2006 | -0.1373   | 1    | 0.0412     | -0.2179         | -0.0566      | \*  |
| 3   | 2004  | 2007 | -0.1008   | 1    | 0.0356     | -0.1706         | -0.0311      | \*  |
| 4   | 2006  | 2004 | 0.0065    | 0    | 0.0249     | -0.0424         | 0.0554       |     |
| 5   | 2006  | 2005 | -0.0028   | 0    | 0.0197     | -0.0414         | 0.0359       |     |
| 6   | 2006  | 2006 | -0.0046   | 1    | 0.0173     | -0.0386         | 0.0294       |     |
| 7   | 2006  | 2007 | -0.0412   | 1    | 0.0197     | -0.0799         | -0.0026      | \*  |
| 8   | 2007  | 2004 | 0.0305    | 0    | 0.0146     | 0.0019          | 0.0591       | \*  |
| 9   | 2007  | 2005 | -0.0027   | 0    | 0.0170     | -0.0361         | 0.0306       |     |
| 10  | 2007  | 2006 | -0.0311   | 0    | 0.0177     | -0.0657         | 0.0035       |     |
| 11  | 2007  | 2007 | -0.0261   | 1    | 0.0181     | -0.0615         | 0.0094       |     |

</div>

In the graphs, a semicolon `;` should be added to prevent printing the
class and the graph information.

``` python
out.plot_attgt();
```

<div id="fig-attgt">

<img src="README_files/figure-commonmark/fig-attgt-output-1.png"
id="fig-attgt" />

Figure 1

</div>

``` python
out.aggte(typec='calendar');
```



    Overall summary of ATT's based on calendar time aggregation:
        ATT Std. Error  [95.0%  Conf. Int.]  
    -0.0417     0.0164 -0.0738      -0.0096 *


    Time Effects (calendar):
       Time  Estimate  Std. Error  [95.0% Simult.   Conf. Band   
    0  2004   -0.0105      0.0242          -0.0580      0.0370   
    1  2005   -0.0704      0.0325          -0.1342     -0.0067  *
    2  2006   -0.0488      0.0195          -0.0871     -0.0105  *
    3  2007   -0.0371      0.0146          -0.0657     -0.0084  *
    ---
    Signif. codes: `*' confidence band does not cover 0
    Control Group:  Never Treated , 
    Anticipation Periods:  0
    Estimation Method:  Doubly Robust

``` python
out.plot_aggte();
```

<div id="fig-calendar">

<img src="README_files/figure-commonmark/fig-calendar-output-1.png"
id="fig-calendar" />

Figure 2

</div>

## Event studies

Although in the current example it is pretty easy to directly interpret
the group-time average treatment effects, there are many cases where it
is convenient to aggregate the group-time average treatment effects into
a small number of parameters. A main type of aggregation is into an
*event study* plot.

To make an event study plot in the **csdid** package, one can use the
**aggte** function with the **dynamic** option:

``` python
out.aggte(typec='dynamic');
```



    Overall summary of ATT's based on event-study/dynamic aggregation:
        ATT Std. Error  [95.0%  Conf. Int.]  
    -0.0772     0.0212 -0.1187      -0.0357 *


    Dynamic Effects:
       Event time  Estimate  Std. Error  [95.0% Simult.   Conf. Band   
    0          -3    0.0305      0.0156          -0.0001      0.0612   
    1          -2   -0.0006      0.0128          -0.0257      0.0246   
    2          -1   -0.0245      0.0141          -0.0522      0.0033   
    3           0   -0.0199      0.0118          -0.0430      0.0032   
    4           1   -0.0510      0.0184          -0.0869     -0.0150  *
    5           2   -0.1373      0.0408          -0.2172     -0.0574  *
    6           3   -0.1008      0.0367          -0.1727     -0.0289  *
    ---
    Signif. codes: `*' confidence band does not cover 0
    Control Group:  Never Treated , 
    Anticipation Periods:  0
    Estimation Method:  Doubly Robust

The column `event time` is for each group relative to when they first
participate in the treatment. For example, `event time=0` corresponds to
the *on impact* effect, and `event time=-1` is the *effect* in the
period before a unit becomes treated. Checking that this is equal to
zero is potentially useful as a pre-test when no anticipation is
specified.

To plot the event study, use **plot_aggte**:

``` python
out.plot_aggte();
```

<div id="fig-dynamic">

<img src="README_files/figure-commonmark/fig-dynamic-output-1.png"
id="fig-dynamic" />

Figure 3

</div>

The figure here is very similar to the group-time average treatment
effects. Red dots are pre-treatment periods, blue dots are
post-treatment periods. The difference is that the x-axis is in event
time.

## Overall effect of participating in treatment

The event study above reported an overall effect of participating in the
treatment. This was computed by averaging the average effects computed
at each length of exposure.

In many cases, a more general purpose overall treatment effect parameter
is given by computing the average treatment effect for each group, and
then averaging across groups. This procedure provides an average
treatment effect parameter with a similar interpretation to the Average
Treatment Effect on the Treated (ATT) in the two-period, two-group case.

To compute this overall average treatment effect parameter, use:

``` python
out.aggte(typec='group');
```



    Overall summary of ATT's based on group/cohort aggregation:
       ATT Std. Error  [95.0%  Conf. Int.]  
    -0.031      0.013 -0.0565      -0.0055 *


    Group Effects:
       Group  Estimate  Std. Error  [95.0% Simult.   Conf. Band   
    0   2004   -0.0797      0.0293          -0.1371     -0.0224  *
    1   2006   -0.0229      0.0166          -0.0554      0.0096   
    2   2007   -0.0261      0.0163          -0.0581      0.0059   
    ---
    Signif. codes: `*' confidence band does not cover 0
    Control Group:  Never Treated , 
    Anticipation Periods:  0
    Estimation Method:  Doubly Robust

Of particular interest is the `Overall ATT` in the results. Here, we
estimate that increasing the minimum wage decreased teen employment by
approximately 3.1% in this example.

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

## How to cite

If you want to cite CSDID, you can use the following BibTeX entry:

``` bibtex
@software{csdid,
  author  = {Callaway, Brantly and Sant'Anna, Pedro HC and Quispe, Alexander and Guevara, Carlos},
  title   = {{csdid: Difference-in-Differences with Multiple Time Periods in Python}},
  year    = {2024},
  url     = {https://github.com/d2cml-ai/csdid}
}
```
