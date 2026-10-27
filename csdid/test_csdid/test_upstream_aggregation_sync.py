"""DID #273: event windows in original time units, with analytic IF oracles."""
import copy

import numpy as np
import pandas as pd
import pytest

from csdid.aggte_fnc.compute_aggte import compute_aggte


def _mp(off_grid=False):
    # Unequal cohort sizes and nonconstant weights expose weight-IF mistakes.
    cohorts = np.array([0, 0, 0, 2004, 2004, 2007, 2007, 2007, 2007])
    if off_grid:
        cohorts[cohorts == 2007] = 2006
    years = np.array([2000, 2004, 2007, 2012])
    weights = np.arange(1, 10, dtype=float) / 5
    data = pd.DataFrame({
        "id": np.repeat(np.arange(9), 4), "year": np.tile(years, 9),
        "cohort": np.repeat(cohorts, 4), "w1": np.repeat(weights, 4),
    })
    groups = np.repeat(np.unique(cohorts[cohorts > 0]), 3)
    times = np.tile(years[1:], 2)
    inf = np.random.default_rng(873).normal(size=(9, 6))
    inf -= inf.mean(axis=0)
    return {
        "group": groups, "t": times, "att": np.array([2., 8., 20., -1., 6., 15.]),
        "inffunc": {"inffunc": inf}, "n": 9,
        "DIDparams": {
            "data": data, "idname": "id", "tname": "year", "gname": "cohort",
            "panel": True, "tlist": years, "glist": np.unique(groups),
            "bstrap": False, "cband": False, "biters": 99,
            "clustervars": None, "alp": .05, "control_group": "nevertreated",
            "est_method": "dr", "anticipation": 0,
        },
    }


def _ratio_oracle(values, influence, unit_mass):
    # Differentiate A = sum_k p_k a_k / sum_k p_k directly.
    mass = unit_mass.mean(axis=0)
    estimate = mass @ values / mass.sum()
    inf = influence @ (mass / mass.sum())
    inf += unit_mass @ (values - estimate) / mass.sum()
    return estimate, inf


def _expected(mp, kind, lo, hi, na_rm=False):
    g, t, a = (np.asarray(mp[k]) for k in ("group", "t", "att"))
    inf = mp["inffunc"]["inffunc"]
    first = mp["DIDparams"]["data"].iloc[::4]
    unit_g, w = first.cohort.to_numpy(), first.w1.to_numpy()
    selected = (t >= g) & (t - g >= lo) & (t - g <= hi)
    if na_rm:
        selected &= np.isfinite(a)
    if kind == "simple":
        unit_mass = w[:, None] * (unit_g[:, None] == g[selected])
        estimate, influence = _ratio_oracle(a[selected], inf[:, selected], unit_mass)
        return estimate, influence, None, None
    groups = np.unique(g[selected])
    means = np.array([a[selected & (g == gg)].mean() for gg in groups])
    group_inf = np.column_stack([inf[:, selected & (g == gg)].mean(axis=1)
                                 for gg in groups])
    unit_mass = w[:, None] * (unit_g[:, None] == groups)
    estimate, influence = _ratio_oracle(means, group_inf, unit_mass)
    return estimate, influence, means, group_inf


@pytest.mark.parametrize("kind", ["simple", "group"])
@pytest.mark.parametrize("off_grid", [False, True])
@pytest.mark.parametrize("lo,hi", [(-np.inf, np.inf), (1, np.inf), (0, 3), (3, 8)])
def test_original_time_windows_att_se_and_influence(kind, off_grid, lo, hi):
    mp = _mp(off_grid)
    estimate, inf, means, group_inf = _expected(mp, kind, lo, hi)
    result = compute_aggte(copy.deepcopy(mp), typec=kind, min_e=lo, max_e=hi)
    np.testing.assert_allclose(result["overall_att"], estimate, atol=1e-13)
    np.testing.assert_allclose(result["overall_se"], np.linalg.norm(inf) / 9, atol=1e-13)
    key = "simple_att" if kind == "simple" else "selective_inf_func"
    np.testing.assert_allclose(result["inf_function"][key].ravel(), inf, atol=1e-13)
    if kind == "group":
        np.testing.assert_allclose(result["att_egt"], means, atol=1e-13)
        np.testing.assert_allclose(result["se_egt"], np.linalg.norm(group_inf, axis=0) / 9)


@pytest.mark.parametrize("lo,hi", [(1, 3), (4, 3), (100, 101)])
def test_group_empty_window_error_and_na_rm(lo, hi):
    mp = _mp()
    with pytest.raises(ValueError, match="window.*group"):
        compute_aggte(copy.deepcopy(mp), typec="group", min_e=lo, max_e=hi)
    if lo == 1:
        expected, inf, means, _ = _expected(mp, "group", lo, hi, na_rm=True)
        result = compute_aggte(mp, typec="group", min_e=lo, max_e=hi, na_rm=True)
        np.testing.assert_array_equal(result["egt"], [2004])
        np.testing.assert_allclose(result["overall_att"], expected)
        np.testing.assert_allclose(result["overall_se"], np.linalg.norm(inf) / 9)
        np.testing.assert_allclose(result["att_egt"], means)
    else:
        with pytest.raises(ValueError):
            compute_aggte(mp, typec="group", min_e=lo, max_e=hi, na_rm=True)


@pytest.mark.parametrize("lo,hi", [(100, 101), (4, 3), (-10, -1)])
def test_simple_empty_window_error(lo, hi):
    with pytest.raises(ValueError, match="min_e.*max_e"):
        compute_aggte(_mp(), typec="simple", min_e=lo, max_e=hi)


def test_na_rm_drops_group_with_only_missing_in_window_cell():
    mp = _mp()
    mp["att"][5] = np.nan  # late cohort has only its e=0 cell left, outside window
    expected, inf, means, _ = _expected(mp, "group", 1, np.inf, na_rm=True)
    result = compute_aggte(mp, typec="group", min_e=1, na_rm=True)
    np.testing.assert_array_equal(result["egt"], [2004])
    np.testing.assert_allclose(result["att_egt"], means)
    np.testing.assert_allclose(result["overall_att"], expected)
    np.testing.assert_allclose(result["overall_se"], np.linalg.norm(inf) / 9)


@pytest.mark.parametrize("kind", ["simple", "group"])
def test_balance_e_warning_preserves_result(kind):
    baseline = compute_aggte(_mp(), typec=kind)
    with pytest.warns(UserWarning, match="balance_e.*ignored"):
        result = compute_aggte(_mp(), typec=kind, balance_e=1)
    np.testing.assert_equal(result["overall_att"], baseline["overall_att"])
    np.testing.assert_equal(result["overall_se"], baseline["overall_se"])


@pytest.mark.parametrize("kind", ["simple", "group", "dynamic", "calendar"])
@pytest.mark.parametrize("na_rm", [False, True])
def test_all_pretreatment_cells_refused(kind, na_rm):
    mp = _mp()
    mp["t"] = np.full(6, 2000)
    with pytest.raises(ValueError, match="No post-treatment ATT"):
        compute_aggte(mp, typec=kind, na_rm=na_rm)
