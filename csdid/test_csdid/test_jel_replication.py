"""
Replication tests for the Callaway & Sant'Anna JEL article.
Translated from R package 'did' tests/testthat/test-jel_replication.R (39 assertions).

These tests verify that csdid Python results match the published article.
Tests are skipped if the JEL data file is not available.
The last R test (faster_mode consistency) is skipped — not in Python csdid.
"""
import os
import tempfile
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt


# ─────────────────────────────────────────────────────────────
# Data helpers
# ─────────────────────────────────────────────────────────────

def _get_jel_data_path():
    """Resolve JEL data path, downloading from GitHub if needed."""
    env_path = os.environ.get(
        "JEL_DID_DATA_PATH",
        os.path.join(os.path.expanduser("~"), "JEL-DiD", "data", "county_mortality_data.csv"),
    )
    if os.path.isfile(env_path):
        return env_path

    tmp_path = os.path.join(tempfile.gettempdir(), "county_mortality_data.csv")
    if os.path.isfile(tmp_path) and os.path.getsize(tmp_path) > 0:
        return tmp_path

    try:
        import urllib.request
        urllib.request.urlretrieve(
            "https://raw.githubusercontent.com/pedrohcgs/JEL-DiD/main/data/county_mortality_data.csv",
            tmp_path,
        )
        if os.path.isfile(tmp_path) and os.path.getsize(tmp_path) > 0:
            return tmp_path
    except Exception:
        pass

    return tmp_path  # skip_if_not will handle missing file


def _load_jel_data(path, filter_2xt=True):
    """Load and clean JEL data."""
    df = pd.read_csv(path, dtype={"county": str})
    df["state"] = df["county"].str[-2:]
    df = df[~df["state"].isin(["DC", "DE", "MA", "NY", "VT"])].copy()

    if filter_2xt:
        df = df[(df["yaca"] == 2014) | df["yaca"].isna() | (df["yaca"] > 2019)].copy()

    df["perc_white"] = df["population_20_64_white"] / df["population_20_64"] * 100
    df["perc_hispanic"] = df["population_20_64_hispanic"] / df["population_20_64"] * 100
    df["perc_female"] = df["population_20_64_female"] / df["population_20_64"] * 100
    df["unemp_rate"] = df["unemp_rate"] * 100
    df["median_income"] = df["median_income"] / 1000

    keep = ["county_code", "year", "population_20_64", "yaca", "crude_rate_20_64",
            "perc_female", "perc_white", "perc_hispanic", "unemp_rate", "poverty_rate",
            "median_income"]
    df = df[keep].copy()

    non_yaca = [c for c in keep if c != "yaca"]
    df = df.dropna(subset=non_yaca)

    both_years = df[df["year"].isin([2013, 2014])].groupby("county_code").size()
    valid = both_years[both_years == 2].index
    df = df[df["county_code"].isin(valid)].copy()

    full_panel = df.groupby("county_code").size()
    valid = full_panel[full_panel == 11].index
    df = df[df["county_code"].isin(valid)].copy()

    return df


@pytest.fixture(scope="module")
def jel_data_path():
    path = _get_jel_data_path()
    if not os.path.isfile(path):
        pytest.skip("JEL-DiD data not available")
    return path


COVS = "crude_rate_20_64 ~ perc_female + perc_white + perc_hispanic + unemp_rate + poverty_rate + median_income"


# ─────────────────────────────────────────────────────────────
# Table 7: 2x2 CS-DiD with Covariates
# ─────────────────────────────────────────────────────────────

def test_jel_table7_2x2(jel_data_path):
    """JEL Table 7: 2x2 CS-DiD point estimates match published values."""
    mydata = _load_jel_data(jel_data_path, filter_2xt=True)
    short = mydata[mydata["year"].isin([2013, 2014])].copy()
    short["treat_year"] = short["yaca"].apply(
        lambda x: 2014 if pd.notna(x) and x == 2014 else 0
    )
    short["county_code"] = short["county_code"].astype(float)

    wt = short.loc[short["year"] == 2013, ["county_code", "population_20_64"]].copy()
    wt.columns = ["county_code", "set_wt"]
    short = short.merge(wt, on="county_code")

    expected = {
        ("reg", None): -1.6154372119,
        ("ipw", None): -0.8585625501,
        ("dr", None): -1.2256473242,
        ("reg", "set_wt"): -3.4592200594,
        ("ipw", "set_wt"): -3.8416966846,
        ("dr", "set_wt"): -3.7561045985,
    }

    for (method, wt_name), exp_att in expected.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(
                yname="crude_rate_20_64", tname="year", idname="county_code",
                gname="treat_year", xformla=COVS, data=short,
                control_group="nevertreated", weights_name=wt_name,
            ).fit(est_method=method, base_period="universal", bstrap=False)
            res.aggte(typec="simple", na_rm=True, bstrap=False)

        assert abs(res.atte["overall_att"] - exp_att) < 0.01, (
            f"Table 7 {method} {'weighted' if wt_name else 'unweighted'}"
        )


# ─────────────────────────────────────────────────────────────
# 2xT Event Study
# ─────────────────────────────────────────────────────────────

def test_jel_2xt_event_study(jel_data_path):
    """JEL 2xT: event study ATT(g,t) point estimates match."""
    mydata = _load_jel_data(jel_data_path, filter_2xt=True)
    mydata["treat_year"] = mydata["yaca"].apply(
        lambda x: 2014 if pd.notna(x) and x == 2014 else 0
    )
    mydata["county_code"] = mydata["county_code"].astype(float)

    wt = mydata.loc[mydata["year"] == 2013, ["county_code", "population_20_64"]].copy()
    wt.columns = ["county_code", "set_wt"]
    mydata = mydata.merge(wt, on="county_code")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="crude_rate_20_64", tname="year", idname="county_code",
            gname="treat_year", data=mydata, weights_name="set_wt",
            control_group="nevertreated",
        ).fit(est_method="reg", base_period="universal", bstrap=False)

    expected_att = [
        4.1292044043, -0.5016807242, 2.7531791360, 2.7804626426,
        0.0, -2.5628745138, -1.6973291127, 0.2189009815,
        -0.8133358354, -1.1532954495, 1.7866564429,
    ]

    att = np.asarray(res.MP["att"], dtype=float)
    np.testing.assert_allclose(att, expected_att, atol=1e-4)

    # Dynamic aggregation
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic", bstrap=False)
    att_egt = np.asarray(res.atte["att_egt"], dtype=float)
    np.testing.assert_allclose(att_egt, expected_att, atol=1e-4)

    # Overall ATT for e in {0,5}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic", min_e=0, max_e=5, bstrap=False)
    assert abs(res.atte["overall_att"] - (-0.7035462478)) < 0.01


# ─────────────────────────────────────────────────────────────
# 2xT with covariates
# ─────────────────────────────────────────────────────────────

def test_jel_2xt_with_covariates(jel_data_path):
    """JEL 2xT: event study with covariates — base period e=-1 should be 0."""
    mydata = _load_jel_data(jel_data_path, filter_2xt=True)
    mydata["treat_year"] = mydata["yaca"].apply(
        lambda x: 2014 if pd.notna(x) and x == 2014 else 0
    )
    mydata["county_code"] = mydata["county_code"].astype(float)

    wt = mydata.loc[mydata["year"] == 2013, ["county_code", "population_20_64"]].copy()
    wt.columns = ["county_code", "set_wt"]
    mydata = mydata.merge(wt, on="county_code")

    for method in ["reg", "ipw", "dr"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ATTgt(
                yname="crude_rate_20_64", tname="year", idname="county_code",
                gname="treat_year", xformla=COVS, data=mydata,
                weights_name="set_wt", control_group="nevertreated",
            ).fit(est_method=method, base_period="universal", bstrap=False)
            res.aggte(typec="dynamic", na_rm=True, bstrap=False)

        egt = np.asarray(res.atte["egt"], dtype=float)
        att_egt = np.asarray(res.atte["att_egt"], dtype=float)

        # e=-1 (base period) should be 0
        base_idx = np.where(egt == -1)[0]
        if len(base_idx) > 0:
            assert abs(att_egt[base_idx[0]]) < 1e-6, f"method={method} e=-1"

        # All estimates should be finite
        non_na = ~np.isnan(att_egt)
        assert np.all(np.isfinite(att_egt[non_na])), f"method={method} finite"


# ─────────────────────────────────────────────────────────────
# GxT: Staggered event study (no covariates)
# ─────────────────────────────────────────────────────────────

def test_jel_gxt_no_covs(jel_data_path):
    """JEL GxT: staggered event study without covariates matches."""
    mydata = _load_jel_data(jel_data_path, filter_2xt=False)
    mydata["treat_year"] = mydata["yaca"].apply(
        lambda x: int(x) if pd.notna(x) and x <= 2019 else 0
    )
    mydata["county_code"] = mydata["county_code"].astype(float)

    wt = mydata.loc[mydata["year"] == 2013, ["county_code", "population_20_64"]].copy()
    wt.columns = ["county_code", "set_wt"]
    mydata = mydata.merge(wt, on="county_code")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="crude_rate_20_64", tname="year", idname="county_code",
            gname="treat_year", data=mydata, weights_name="set_wt",
            control_group="notyettreated",
        ).fit(est_method="reg", base_period="universal", bstrap=False)
        res.aggte(typec="dynamic", bstrap=False)

    expected = {
        -5: 2.2186783578, -4: 0.8579225109, -3: 1.9161499399,
        -2: 2.5644742860, -1: 0.0,
        0: -1.6545648988, 1: -0.2616435324, 2: 1.7055625922,
        3: -0.5405232028, 4: -0.5148819184, 5: 1.7866564429,
    }

    egt = np.asarray(res.atte["egt"], dtype=float)
    att_egt = np.asarray(res.atte["att_egt"], dtype=float)

    for e_val, exp in expected.items():
        idx = np.where(egt == e_val)[0]
        if len(idx) > 0:
            assert abs(att_egt[idx[0]] - exp) < 0.01, f"e={e_val}"

    # Overall ATT e in {0,5}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic", min_e=0, max_e=5, bstrap=False)
    assert abs(res.atte["overall_att"] - 0.0867675805) < 0.01


# ─────────────────────────────────────────────────────────────
# GxT: DR with covariates
# ─────────────────────────────────────────────────────────────

def test_jel_gxt_dr_covs(jel_data_path):
    """JEL GxT: staggered event study with DR covariates matches."""
    mydata = _load_jel_data(jel_data_path, filter_2xt=False)
    mydata["treat_year"] = mydata["yaca"].apply(
        lambda x: int(x) if pd.notna(x) and x <= 2019 else 0
    )
    mydata["county_code"] = mydata["county_code"].astype(float)

    wt = mydata.loc[mydata["year"] == 2013, ["county_code", "population_20_64"]].copy()
    wt.columns = ["county_code", "set_wt"]
    mydata = mydata.merge(wt, on="county_code")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ATTgt(
            yname="crude_rate_20_64", tname="year", idname="county_code",
            gname="treat_year", xformla=COVS, data=mydata,
            weights_name="set_wt", control_group="notyettreated",
        ).fit(est_method="dr", base_period="universal", bstrap=False)

    # Overall ATT for e in {0,5}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic", min_e=0, max_e=5, bstrap=False)
    assert abs(res.atte["overall_att"] - (-2.2469982988)) < 0.05

    # Dynamic at key event times
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res.aggte(typec="dynamic", bstrap=False)

    expected = {
        -5: 2.6684811691, -4: 2.1333836537, -3: 2.8574389179,
        -2: 2.9129673584, -1: 0.0,
        0: -1.4929553032, 1: -1.9693922262, 2: -2.7250659699,
        3: -5.0556625460, 4: -4.7161630370, 5: 2.4772492894,
    }

    egt = np.asarray(res.atte["egt"], dtype=float)
    att_egt = np.asarray(res.atte["att_egt"], dtype=float)

    for e_val, exp in expected.items():
        idx = np.where(egt == e_val)[0]
        if len(idx) > 0:
            assert abs(att_egt[idx[0]] - exp) < 0.05, f"e={e_val}"
