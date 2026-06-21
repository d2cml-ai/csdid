"""Tests for faster_mode consistency.

Translated from R package 'did' tests/testthat/test-faster-mode-consistency.R,
whose 87 tests all compare faster_mode=TRUE vs faster_mode=FALSE.

csdid now implements faster_mode (vectorized precompute of covariates and
per-period outcomes) which is numerically identical to the standard path. These
tests verify that equivalence across panel/RC, control groups, est_methods and
base periods.
"""
import warnings

import numpy as np
import pytest

from csdid.att_gt import ATTgt
from helpers import build_sim_data


def _att_se(data, faster_mode, **kw):
    common = dict(yname="Y", tname="period", idname="id", gname="G", data=data, xformla="Y~X")
    common.update({k: kw[k] for k in ("control_group", "panel", "allow_unbalanced_panel") if k in kw})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(**common, faster_mode=faster_mode).fit(
            est_method=kw.get("est_method", "dr"),
            base_period=kw.get("base_period", "varying"), bstrap=False)
    return (np.asarray(res.MP["att"], dtype=float),
            np.asarray(res.results["se"], dtype=float))


@pytest.mark.parametrize("panel", [True, False])
@pytest.mark.parametrize("control_group", ["nevertreated", "notyettreated"])
@pytest.mark.parametrize("est_method", ["dr", "reg", "ipw"])
@pytest.mark.parametrize("base_period", ["varying", "universal"])
def test_faster_mode_equals_standard(panel, control_group, est_method, base_period):
    """ATT and SE are identical between faster_mode and the standard path."""
    data = build_sim_data(n=1000, time_periods=4, te=0, te_e=[1, 2, 3, 4], seed=42)
    kw = dict(control_group=control_group, panel=panel,
              est_method=est_method, base_period=base_period)
    a_std, se_std = _att_se(data, False, **kw)
    a_fast, se_fast = _att_se(data, True, **kw)
    np.testing.assert_allclose(a_std, a_fast, atol=1e-9, equal_nan=True)
    fin = np.isfinite(se_std) & np.isfinite(se_fast)
    np.testing.assert_allclose(se_std[fin], se_fast[fin], atol=1e-9)
