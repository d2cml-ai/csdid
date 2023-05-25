# Library Creation Process

## R DID PACKAGE

[DID](https://github.com/bcallaway11/did/blob/master/R/DIDparams.R)

## `agg_te`

- "`did/aggte/utils.py`" file
  - `wif, get_agg_inf_func, getSE -> compute_aggte`
    - [R ref](https://github.com/bcallaway11/did/blob/master/R/compute.aggte.R)
    - [Python ref](https://github.com/bernardodionisi/differences/blob/2b6acc94cd37105893ba9dc1f9a48e37fb916c4f/differences/attgt/aggregate.py#L258)
      - R `wif` -> Python `get_wif`
      - R `get_agg_inf_func` -> Python `get_agg_influence_func`
- "`did/aggte/mboot.py`" file
  - `sum_multiplier_bootstrap`, `mboot`
    - [R ref](https://github.com/bcallaway11/did/blob/master/R/mboot.R)
    - [python ref](https://github.com/bernardodionisi/differences/blob/2b6acc94cd37105893ba9dc1f9a48e37fb916c4f/differences/attgt/mboot.py#L15)
- "`did/aggte/aggte.py`" file
  - `aggte`
    - [R ref](https://github.com/bcallaway11/did/blob/master/R/aggte.R)
    - [Python ref]

![](did/aggte/aggte.png)

## `agg_gt`

- `pre_process_did`:
  - [R ref](https://github.com/bcallaway11/did/blob/master/R/pre_process_did.R)
- `process_att_gt`:
  - [R ref](https://github.com/bcallaway11/did/blob/master/R/process_attgt.R)
- `compute_att_gt`"
  - [R ref](https://github.com/bcallaway11/did/blob/master/R/compute.att_gt.R)
- `att_gt`:
  - [R ref](https://github.com/bcallaway11/did/blob/master/R/att_gt.R)

### `DRDID` package

[DRDID](https://github.com/pedrohcgs/DRDID/tree/master)

To Build o adapt:

- `std_ipw_did_panel`:
  - [R ref](https://github.com/pedrohcgs/DRDID/blob/master/R/std_ipw_did_rc.R)
- `reg_did_panel`:
  - [R ref](https://github.com/pedrohcgs/DRDID/blob/master/R/reg_did_panel.R)
  - [Python `Difference` package](https://github.com/bernardodionisi/differences/blob/main/differences/did/did_cal.py)
- `drdid_panel`:
  - [R ref](https://github.com/pedrohcgs/DRDID/blob/master/R/drdid_panel.R)
  - [Python `Difference` package](https://github.com/bernardodionisi/differences/blob/main/differences/did/did_cal.py)
- `std_ipw_did_rc`:
  - [R ref](https://github.com/pedrohcgs/DRDID/blob/master/R/std_ipw_did_rc.R)
  - [Python `Difference` package](https://github.com/bernardodionisi/differences/blob/main/differences/did/did_cal.py)
- `reg_did_rc`:
  - [R ref](https://github.com/pedrohcgs/DRDID/blob/master/R/reg_did_rc.R)
  - [Python `Difference` package](https://github.com/bernardodionisi/differences/blob/main/differences/did/did_cal.py)
- `drdid_rc`:
  - [R ref](https://github.com/pedrohcgs/DRDID/blob/master/R/drdid.R)
  - [Python `Difference` package](https://github.com/bernardodionisi/differences/blob/main/differences/did/did_cal.py)

![](did/att_gt/att_gt.png)

## `simulate data`

![](did/sim_data/sim.png)
