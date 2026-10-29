"""Bespoke, self-contained, mutation-gated unit tests for the aggregation math in
``csdid/aggte_fnc/compute_aggte.py`` and the pure helpers in
``csdid/aggte_fnc/utils.py`` (wif / get_agg_inf_func / get_se).

NO R at runtime. The aggregated ATTs are asserted against HAND-COMPUTED weighted
averages of the per-(g,t) ATTs, using the SAME pg = mean(w1 * 1{G==g}) group
weights the code uses. These oracles pin the keepers/event-time/calendar
selection logic and every aggregation sum/mean:

  simple   = sum(att[keepers]*pg) / sum(pg)             (keepers: g<=t)
  group_g  = mean(att | group==g, t>=g);  overall = pgg-weighted avg
  dynamic_e= pg-weighted avg within event time e=t-g;  overall = mean over e>=0
  calendar_t= pg-weighted avg within t (g<=t);          overall = mean over t

Pure helpers use closed forms: get_se = sqrt(mean(IF^2)/n);
get_agg_inf_func = inffunc[:,which] @ weights.
"""
import contextlib
import io
import warnings

import numpy as np
import pytest

from csdid.aggte_fnc.compute_aggte import compute_aggte
from csdid.aggte_fnc.utils import get_se, get_agg_inf_func, wif
from helpers import build_sim_data
from csdid.att_gt import ATTgt


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
def _fit_mp(seed=1, n=600, te=2.0):
    df = build_sim_data(n=n, time_periods=4, te=te, seed=seed)
    att = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=df)
    with contextlib.redirect_stdout(io.StringIO()):
        att.fit(est_method="dr", bstrap=False)
    return att.MP


def _fit_mp_cband(seed=1, n=600, te=2.0, biters=300):
    df = build_sim_data(n=n, time_periods=4, te=te, seed=seed)
    att = ATTgt(yname="Y", tname="period", idname="id", gname="G", data=df,
                cband=True, biters=biters)
    with contextlib.redirect_stdout(io.StringIO()):
        att.fit(est_method="dr", bstrap=True)
    return att.MP


def _aggte(MP, typec, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return compute_aggte(MP, typec=typec, **kw)


def _pg_by_group(MP):
    dp = MP["DIDparams"]
    data = dp["data"]
    tlist = np.array(dp["tlist"])
    dta = data[data[dp["tname"]] == tlist[0]]
    w1 = dta["w1"].to_numpy()
    gcol = dta[dp["gname"]].to_numpy()
    glist = np.array(dp["glist"])
    return {gg: np.mean(w1 * (gcol == gg)) for gg in glist}


def _cells(MP):
    return (np.array(MP["group"]), np.array(MP["t"]), np.array(MP["att"]))


# --------------------------------------------------------------------------- #
# SIMPLE  (lines 159 keepers, 173 weighted sum)
# --------------------------------------------------------------------------- #
def test_simple_overall_matches_weighted_average():
    MP = _fit_mp()
    g, t, a = _cells(MP)
    pgmap = _pg_by_group(MP)
    keepers = [i for i in range(len(g)) if g[i] <= t[i]]
    pg = np.array([pgmap[g[i]] for i in range(len(g))])
    exp = np.sum(a[keepers] * pg[keepers]) / np.sum(pg[keepers])
    got = _aggte(MP, "simple")["overall_att"]
    assert abs(got - exp) < 1e-9


# --------------------------------------------------------------------------- #
# GROUP  (lines 211 per-group, 260 overall)
# --------------------------------------------------------------------------- #
def test_group_per_group_and_overall():
    MP = _fit_mp()
    g, t, a = _cells(MP)
    pgmap = _pg_by_group(MP)
    glist = np.array(MP["DIDparams"]["glist"])
    att_g = []
    for gg in glist:
        cells = [i for i in range(len(g)) if g[i] == gg and t[i] >= gg]
        att_g.append(np.mean(a[cells]))
    att_g = np.array(att_g)
    pgg = np.array([pgmap[gg] for gg in glist])
    overall = np.sum(att_g * pgg) / np.sum(pgg)

    r = _aggte(MP, "group")
    assert np.allclose(np.asarray(r["att_egt"], float), att_g, atol=1e-9)
    assert abs(r["overall_att"] - overall) < 1e-9


# --------------------------------------------------------------------------- #
# DYNAMIC  (lines 327 window, 333 within-e weights, 380 overall mean)
# --------------------------------------------------------------------------- #
def test_dynamic_event_times_and_overall():
    MP = _fit_mp()
    g, t, a = _cells(MP)
    pgmap = _pg_by_group(MP)
    e = t - g
    dyn = {}
    for ev in sorted(set(e)):
        cells = np.where(e == ev)[0]
        pge = np.array([pgmap[g[i]] for i in cells])
        pge = pge / pge.sum()
        dyn[ev] = np.sum(a[cells] * pge)

    r = _aggte(MP, "dynamic")
    egt = list(r["egt"])
    assert egt == sorted(set(e))
    assert np.allclose(np.asarray(r["att_egt"], float),
                       [dyn[ev] for ev in egt], atol=1e-9)
    overall = np.mean([dyn[ev] for ev in egt if ev >= 0])
    assert abs(r["overall_att"] - overall) < 1e-9


def test_dynamic_max_e_window():
    MP = _fit_mp()
    r = _aggte(MP, "dynamic", max_e=1)
    assert max(r["egt"]) <= 1


def test_dynamic_min_e_window():
    MP = _fit_mp()
    r = _aggte(MP, "dynamic", min_e=0)
    assert min(r["egt"]) >= 0


# --------------------------------------------------------------------------- #
# CALENDAR  (lines 424 selection, 426 weights, 484 overall mean)
# --------------------------------------------------------------------------- #
def test_calendar_per_time_and_overall():
    MP = _fit_mp()
    g, t, a = _cells(MP)
    pgmap = _pg_by_group(MP)
    minG = g.min()
    cal_t = sorted([tt for tt in np.unique(t) if tt >= minG])
    cal = {}
    for tt in cal_t:
        cells = np.where((t == tt) & (g <= t))[0]
        pgt = np.array([pgmap[g[i]] for i in cells])
        pgt = pgt / pgt.sum()
        cal[tt] = np.sum(a[cells] * pgt)

    r = _aggte(MP, "calendar")
    assert np.allclose(np.asarray(r["att_egt"], float),
                       [cal[tt] for tt in cal_t], atol=1e-9)
    assert abs(r["overall_att"] - np.mean(list(cal.values()))) < 1e-9


# --------------------------------------------------------------------------- #
# validation + na_rm + SE positivity
# --------------------------------------------------------------------------- #
def test_invalid_type_raises():
    MP = _fit_mp()
    with pytest.raises(ValueError):
        _aggte(MP, "bogus")


def test_all_types_have_positive_overall_se():
    MP = _fit_mp()
    for tp in ("simple", "group", "dynamic", "calendar"):
        r = _aggte(MP, tp)
        se = float(r["overall_se"])
        assert np.isfinite(se) and se > 0


# --------------------------------------------------------------------------- #
# pure helpers (utils.py)
# --------------------------------------------------------------------------- #
def test_get_se_analytic_formula():
    """Analytical (non-bootstrap) SE = sqrt(mean(IF^2)/n)."""
    inf = np.array([1.0, -2.0, 3.0, -1.0, 0.5])
    se = get_se(inf, None)
    exp = np.sqrt(np.mean(inf ** 2) / len(inf))
    assert abs(se - exp) < 1e-12


def test_get_agg_inf_func_weighted_combination():
    """get_agg_inf_func = inffunc[:, which] @ weights_agg."""
    inffunc = np.arange(12.0).reshape(4, 3)
    w = np.array([0.5, 0.3, 0.2])
    out = get_agg_inf_func(att=None, inffunc=inffunc, whichones=[0, 1, 2],
                           weights_agg=w)
    assert np.allclose(out, inffunc @ w)


def test_get_agg_inf_func_adds_weight_if():
    """With wif given, the weight-IF term att[which] @ wif is added."""
    inffunc = np.zeros((4, 2))
    att = np.array([5.0, 7.0])
    wifmat = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]])
    out = get_agg_inf_func(att=att, inffunc=inffunc, whichones=[0, 1],
                           weights_agg=np.array([0.5, 0.5]), wif=wifmat)
    # inffunc contribution is 0; only the wif @ att[which] term remains
    assert np.allclose(out, wifmat @ att)


# --------------------------------------------------------------------------- #
# bootstrap simultaneous bands (cband=True): the simultaneous critical value
# exceeds the pointwise 1.96, and is finite (< 7). Exercises the cband mboot
# branches (lines 243/364/466) and the crit-val guard comparisons.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("typec", ["group", "dynamic", "calendar"])
def test_cband_crit_val_exceeds_pointwise(typec):
    MP = _fit_mp_cband()
    np.random.seed(0)
    r = _aggte(MP, typec)
    cv = float(r["crit_val_egt"])
    pointwise = 1.959963984540054  # norm.ppf(0.975)
    assert cv >= pointwise          # simultaneous >= pointwise
    assert cv < 7.0                 # finite / reliable (the >=7 warn-guard line)


# --------------------------------------------------------------------------- #
# balance_e windowing (dynamic): restricting to a balanced panel of event times
# bounds the event-time sequence at balance_e and drops longer horizons.
# Pins the balance_e selection arithmetic (lines 323/326).
# --------------------------------------------------------------------------- #
def test_dynamic_balance_e_bounds_event_times():
    MP = _fit_mp()
    full = _aggte(MP, "dynamic")
    bal = _aggte(MP, "dynamic", balance_e=1)
    assert max(bal["egt"]) <= 1
    # balance_e drops the longest event times present without balancing
    assert max(bal["egt"]) < max(full["egt"])


def test_wif_mean_zero_columns():
    """Each column of the weight influence function is mean-zero (it is a centered
    indicator minus its probability), a defining property. Any sign flip in the
    numerator/denominator breaks it."""
    n = 50
    rng = np.random.default_rng(0)
    G = rng.integers(0, 3, n).astype(float)
    group = np.array([0.0, 1.0, 2.0])
    pg = np.array([np.mean(G == k) for k in group])
    weights_ind = np.ones(n)
    keepers = [0, 1, 2]
    wf = wif(keepers, pg, weights_ind, G, group)
    # mean over units of each column ~ 0
    assert np.max(np.abs(wf.mean(axis=0))) < 1e-12
