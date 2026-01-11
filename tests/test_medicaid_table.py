import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from csdid.att_gt import ATTgt

EXPECTED = {
    "unweighted": {
        "reg": {"att": -1.6154372, "se": 4.678972},
        "ipw": {"att": -0.8585626, "se": 4.611698},
        "dr": {"att": -1.2256473, "se": 4.942809},
    },
    "weighted": {
        "reg": {"att": -3.4592201, "se": 2.388022},
        "ipw": {"att": -3.8416967, "se": 3.394186},
        "dr": {"att": -3.7561046, "se": 3.240180},
    },
}


def load_data(path):
    data = pd.read_csv(path)
    data["state"] = data["county"].str[-2:]
    data = data[~data["state"].isin(["DC", "DE", "MA", "NY", "VT"])]
    data = data[(data["yaca"] == 2014) | (data["yaca"].isna()) | (data["yaca"] > 2019)]

    covs = [
        "perc_female",
        "perc_white",
        "perc_hispanic",
        "unemp_rate",
        "poverty_rate",
        "median_income",
    ]

    data = data.assign(
        perc_white=data["population_20_64_white"] / data["population_20_64"] * 100,
        perc_hispanic=data["population_20_64_hispanic"] / data["population_20_64"] * 100,
        perc_female=data["population_20_64_female"] / data["population_20_64"] * 100,
        unemp_rate=data["unemp_rate"] * 100,
        median_income=data["median_income"] / 1000,
    )

    keep_cols = [
        "state",
        "county",
        "county_code",
        "year",
        "population_20_64",
        "yaca",
    ] + [c for c in data.columns if c.startswith("perc_")] + [
        "crude_rate_20_64",
    ] + covs
    keep_cols = list(dict.fromkeys(keep_cols))
    data = data[keep_cols]

    cols_except_yaca = [c for c in data.columns if c != "yaca"]
    data = data.dropna(subset=cols_except_yaca)

    data = data.groupby("county_code").filter(
        lambda x: ((x["year"] == 2013) | (x["year"] == 2014)).sum() == 2
    )
    data = data.groupby("county_code").filter(
        lambda x: x["crude_rate_20_64"].notna().sum() == 11
    )

    short_data = data.copy()
    short_data["Treat"] = (
        (short_data["yaca"] == 2014) & (~short_data["yaca"].isna())
    ).astype(int)
    short_data["Post"] = (short_data["year"] == 2014).astype(int)
    short_data = short_data[short_data["year"].isin([2013, 2014])].copy()
    short_data["set_wt"] = short_data.groupby("county_code")["population_20_64"].transform(
        lambda s: s.loc[short_data.loc[s.index, "year"] == 2013].iloc[0]
    )

    short_data["treat_year"] = np.where(
        (short_data["yaca"] == 2014) & (~short_data["yaca"].isna()), 2014, 0
    )
    short_data["county_code"] = pd.to_numeric(short_data["county_code"], errors="coerce")

    return short_data, covs


def run_att_gt(data, covs, method, weights_name, biters=25000):
    np.random.seed(20240924)
    xformla = "~ " + " + ".join(covs)
    did = ATTgt(
        yname="crude_rate_20_64",
        tname="year",
        idname="county_code",
        gname="treat_year",
        data=data,
        xformla=xformla,
        panel=True,
        control_group=["nevertreated"],
        biters=biters,
        weights_name=weights_name,
    ).fit(est_method=method, base_period="universal", bstrap=True)

    np.random.seed(20240924)
    ag = did.aggte(typec="group", na_rm=True, bstrap=True, biters=biters)
    egt = list(ag.atte.get("egt"))
    idx = egt.index(2014)
    att = ag.atte.get("att_egt")[idx]
    se = ag.atte.get("se_egt")[idx]
    if isinstance(se, (list, tuple, np.ndarray)):
        se = se[0]
    return float(att), float(se)


def test_medicaid_expansion_table():
    data, covs = load_data("tests/fixtures/county_mortality_data.csv")

    for weight_key, weights_name in [("unweighted", None), ("weighted", "set_wt")]:
        for method in ["reg", "ipw", "dr"]:
            att, se = run_att_gt(data, covs, method, weights_name)
            exp_att = EXPECTED[weight_key][method]["att"]
            exp_se = EXPECTED[weight_key][method]["se"]

            assert abs(att - exp_att) < 1e-6, f"ATT mismatch for {weight_key}/{method}"
            assert abs(se - exp_se) < 0.05, f"SE mismatch for {weight_key}/{method}"
