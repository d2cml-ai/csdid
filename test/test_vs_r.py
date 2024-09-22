import pytest

import numpy as np
import pandas as pd
import pytest
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# rpy2 imports
from rpy2.robjects.packages import importr
from csdid.att_gt import ATTgt


pandas2ri.activate()

did = importr("did")

@pytest.fixture
def data():

  return pd.read_csv("https://raw.githubusercontent.com/d2cml-ai/csdid/function-aggte/data/mpdta.csv")

def check_absolute_diff(x1, x2, tol, msg=None):
    msg = "" if msg is None else msg
    assert np.all(np.abs(x1 - x2) < tol), msg

def test_ate(data):

  "Test simple ATE via Py vs R."

  py_did = ATTgt(
    yname = "lemp",
    gname = "first.treat",
    idname = "countyreal",
    tname = "year",
    data = data,
    biters = 10_000,
  ).fit(est_method = 'dr')

  py_res = py_did.aggte("simple")
  py_coef = py_res.atte.get("overall_att")
  py_ses = py_res.atte.get("overall_se")

  r_did = did.att_gt(
    yname = "lemp",
    gname = "first.treat",
    idname = "countyreal",
    tname = "year",
    data = data,
    biters = 10_000
  )

  r_coef = did.aggte(r_did, type = "simple").rx2('overall.att')
  r_se = did.aggte(r_did, type = "simple").rx2('overall.se')

  check_absolute_diff(py_coef, r_coef, 1e-5, "ATEs are not equal.")
  check_absolute_diff(py_se, r_se, 1e-3, "SEs are not equal.")
