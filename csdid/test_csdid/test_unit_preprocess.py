"""Bespoke, self-contained, mutation-gated unit tests for the input validation and
preprocessing logic in ``csdid/attgt_fnc/preprocess_did.py``.

NO R at runtime. ``_validate_inputs`` is a pure guard function: each test feeds a
specific bad input and asserts the corresponding ValueError, pinning every
threshold/comparison in the guard (>1 cluster var, time-varying cluster, negative
gname, negative / non-positive-mean weights, invalid control_group, negative
anticipation, duplicate (id,t), irreversible-treatment, panel-without-idname,
strict-numeric dtype). The structural pre_process_did paths (anticipation-text in
the dropped-first-period warning, time-varying-weights warning) are exercised
through public construction with crafted data.
"""
import contextlib
import io
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.attgt_fnc.preprocess_did import (
    _validate_inputs, _strictly_numeric, pre_process_did,
)
from helpers import build_sim_data


# --------------------------------------------------------------------------- #
# _validate_inputs : pure guard
# --------------------------------------------------------------------------- #
def _mk(**over):
    n = 40
    d = pd.DataFrame({
        "Y": np.random.default_rng(0).standard_normal(n * 2),
        "period": np.tile([1, 2], n),
        "id": np.repeat(np.arange(n), 2),
        "G": np.repeat([0, 2], n),
        "cl": np.repeat(np.arange(n), 2),
        "w": np.ones(n * 2),
    })
    for k, v in over.items():
        d[k] = v
    return d


_BASE = dict(yname="Y", tname="period", idname="id", gname="G",
             control_group="nevertreated", anticipation=0, panel=True,
             clustervar=None, weights_name=None)


def _validate(data, **over):
    # _validate_inputs now returns the normalized
    # (yname, tname, gname, idname, weights_name, clustervar, anticipation) tuple
    # (it also normalizes/unwraps the name + anticipation args, mirroring R
    # `validate_*`); the cluster-normalization contract these tests pin is slot 5.
    return _validate_inputs(data=data, **{**_BASE, **over})[5]


def test_validate_normal_returns_none_clustervar():
    assert _validate(_mk()) is None


def test_validate_unwraps_single_clustervar():
    assert _validate(_mk(), clustervar=["cl"]) == "cl"


def test_validate_drops_idname_from_clustervars():
    # C2-F7: R removes idname from the cluster set BEFORE the "at most one" count
    # check, so [other, idname] is a single effective cluster var -> accepted.
    d = _mk()
    d["cl2"] = np.repeat(np.arange(40), 2)
    assert _validate(d, clustervar=["cl", "id"]) == "cl"


def test_validate_rejects_multiple_clustervars():
    # Two genuine (non-idname) cluster vars still exceed the one-beyond-unit limit.
    d = _mk()
    d["cl2"] = np.repeat(np.arange(40), 2)
    with pytest.raises(ValueError):
        _validate(d, clustervar=["cl", "cl2"])


def test_validate_rejects_time_varying_cluster():
    d = _mk()
    d["cl"] = np.arange(80)  # unique per row -> varies within unit over time
    with pytest.raises(ValueError):
        _validate(d, clustervar="cl")


def test_validate_rejects_negative_gname():
    with pytest.raises(ValueError):
        _validate(_mk(G=np.repeat([0, -2], 40)))


def test_validate_rejects_negative_weights():
    w = np.ones(80)
    w[0] = -1.0
    with pytest.raises(ValueError):
        _validate(_mk(w=w), weights_name="w")


def test_validate_rejects_nonpositive_mean_weights():
    with pytest.raises(ValueError):
        _validate(_mk(w=np.zeros(80)), weights_name="w")


def test_validate_accepts_positive_weights():
    # a strictly-positive weight column must pass (kills weight-guard mutants that
    # would over-reject)
    assert _validate(_mk(w=np.full(80, 2.0)), weights_name="w") is None


def test_validate_rejects_invalid_control_group():
    with pytest.raises(ValueError):
        _validate(_mk(), control_group="bogus")


def test_validate_rejects_negative_anticipation():
    with pytest.raises(ValueError):
        _validate(_mk(), anticipation=-1)


def test_validate_rejects_duplicate_id_time():
    d = _mk()
    d.loc[1, "period"] = 1  # id 0 now appears twice in period 1
    with pytest.raises(ValueError):
        _validate(d)


def test_validate_rejects_reversible_treatment():
    d = _mk()
    d.loc[1, "G"] = 2  # id 0 has G in {0, 2} across periods
    with pytest.raises(ValueError):
        _validate(d)


def test_validate_panel_requires_idname():
    with pytest.raises(ValueError):
        _validate(_mk(), idname=None)


def test_validate_rejects_nonnumeric_time():
    d = _mk()
    d["period"] = d["period"].astype(str)
    with pytest.raises(ValueError):
        _validate(d)


# --------------------------------------------------------------------------- #
# _strictly_numeric
# --------------------------------------------------------------------------- #
def test_strictly_numeric_int_true():
    assert _strictly_numeric(pd.Series([1, 2, 3]))


def test_strictly_numeric_bool_false():
    assert not _strictly_numeric(pd.Series([True, False, True]))


def test_strictly_numeric_string_false():
    assert not _strictly_numeric(pd.Series(["a", "b"]))


# --------------------------------------------------------------------------- #
# pre_process_did : structural warning paths
# --------------------------------------------------------------------------- #
def test_time_varying_weights_warning():
    """Panel + within-unit time-varying weights -> the time-varying-weights
    warning fires regardless of fix_weights (F45-2). Pins the `w_tv > 1` guard
    (line 677) and the `input_panel` gate."""
    df = build_sim_data(n=200, time_periods=3, seed=2)
    df["wt"] = 1.0
    df.loc[df["period"] == 2, "wt"] = 2.0
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(io.StringIO()):
            pre_process_did(yname="Y", tname="period", idname="id", gname="G",
                            data=df, weights_name="wt", xformla="~X")
    assert any("Time-varying weights" in str(x.message) for x in w)


def test_constant_weights_no_tv_warning():
    """Constant weights must NOT trigger the time-varying-weights warning
    (kills a `nunique() > 1` -> `>= 1` / `> 0` mutant that would always warn)."""
    df = build_sim_data(n=200, time_periods=3, seed=2)
    df["wt"] = 1.0
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(io.StringIO()):
            pre_process_did(yname="Y", tname="period", idname="id", gname="G",
                            data=df, weights_name="wt", xformla="~X")
    assert not any("Time-varying weights" in str(x.message) for x in w)


def test_anticipation_text_in_dropped_first_period_warning():
    """When anticipation > 0 and units are dropped for being treated in the first
    period, the warning includes the anticipation note (line 579 `anticipation>0`)."""
    df = build_sim_data(n=200, time_periods=4, seed=3)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(io.StringIO()):
            pre_process_did(yname="Y", tname="period", idname="id", gname="G",
                            data=df, anticipation=1, xformla="~X")
    msgs = [str(x.message) for x in w]
    assert any("anticipation = 1" in m for m in msgs)
