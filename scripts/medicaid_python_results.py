# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Medicaid Expansion Results (Python)
# Compute ATT and SE for unweighted and weighted configurations, then write CSV output.

# %% [markdown]
# ## Imports and Paths

# %%
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    ROOT = Path(__file__).resolve().parents[1]
    SCRIPT_DIR = Path(__file__).resolve().parent
except NameError:
    ROOT = Path.cwd()
    SCRIPT_DIR = ROOT / "scripts"

sys.path.insert(0, str(ROOT))

from csdid.att_gt import ATTgt

# %% [markdown]
# ## Configuration

# %%
COVS = [
    "perc_female",
    "perc_white",
    "perc_hispanic",
    "unemp_rate",
    "poverty_rate",
    "median_income",
]
METHODS = ["reg", "ipw", "dr"]
BITERS = 25000
SEED = 20240924

# %% [markdown]
# ## Data Preparation

# %%
def load_data(path):
    data = pd.read_csv(path)
    data["state"] = data["county"].str[-2:]
    data = data[~data["state"].isin(["DC", "DE", "MA", "NY", "VT"])]
    data = data[(data["yaca"] == 2014) | (data["yaca"].isna()) | (data["yaca"] > 2019)]

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
    ] + COVS
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

    return short_data

# %% [markdown]
# ## Estimation Helper

# %%
def run_att_gt(data, method, weights_name):
    np.random.seed(SEED)
    xformla = "~ " + " + ".join(COVS)
    did = ATTgt(
        yname="crude_rate_20_64",
        tname="year",
        idname="county_code",
        gname="treat_year",
        data=data,
        xformla=xformla,
        panel=True,
        control_group=["nevertreated"],
        biters=BITERS,
        weights_name=weights_name,
    ).fit(est_method=method, base_period="universal", bstrap=True)

    np.random.seed(SEED)
    ag = did.aggte(typec="group", na_rm=True, bstrap=True, biters=BITERS)
    egt = list(ag.atte.get("egt"))
    idx = egt.index(2014)
    att = ag.atte.get("att_egt")[idx]
    se = ag.atte.get("se_egt")[idx]
    if isinstance(se, (list, tuple, np.ndarray)):
        se = se[0]
    return float(att), float(se)

# %% [markdown]
# ## CSV Writer

# %%
def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

# %% [markdown]
# ## Run and Save

# %%
def main():
    data_path = ROOT / "tests" / "fixtures" / "county_mortality_data.csv"
    out_path = SCRIPT_DIR / "medicaid_python_results.csv"

    data = load_data(data_path)
    rows = []
    for weight_key, weights_name in [("unweighted", None), ("weighted", "set_wt")]:
        for method in METHODS:
            att, se = run_att_gt(data, method, weights_name)
            rows.append(
                {
                    "weighting": weight_key,
                    "method": method,
                    "att": att,
                    "se": se,
                    "source": "python",
                }
            )

    write_csv(out_path, rows)
    print(f"Python results: {out_path}")


if __name__ == "__main__":
    main()
