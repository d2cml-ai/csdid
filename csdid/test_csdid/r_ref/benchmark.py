"""Benchmark csdid faster_mode vs standard. Prints a table and writes bench_py.csv.

Run: python csdid/test_csdid/r_ref/benchmark.py
Pairs with benchmark.R for the full R(Fast/Standard) vs Python(Fast/Standard) table.
"""
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir, os.pardir))
sys.path.insert(0, os.path.join(HERE, os.pardir))  # csdid/test_csdid -> helpers
sys.path.insert(0, REPO_ROOT)  # repo root -> csdid package
from helpers import build_sim_data  # noqa: E402
from csdid.att_gt import ATTgt  # noqa: E402

SIM = os.path.join(HERE, "sim")

DATASETS = [
    ("tp8_dyn", os.path.join(SIM, "tp8_dyn.csv")),
    ("tp10_const", os.path.join(SIM, "tp10_const.csv")),
    ("big_tp10", os.path.join(SIM, "big_tp10.csv")),
]


def _ensure_big():
    p = os.path.join(SIM, "big_tp10.csv")
    if not os.path.isfile(p):
        d = build_sim_data(n=8000, time_periods=10, te=0, te_e=list(range(1, 11)), seed=99)
        d.to_csv(p, index=False)


def _time(d, faster_mode, reps=5):
    best = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        ATTgt(yname="Y", tname="period", idname="id", gname="G", data=d,
              control_group="nevertreated", xformla="Y~X",
              faster_mode=faster_mode).fit("dr", bstrap=False)
        best = min(best, time.perf_counter() - t0)
    return best * 1000.0  # ms


def main():
    _ensure_big()
    rows = []
    print(f"{'dataset':14s} {'rows':>7s} {'Py(Standard)':>13s} {'Py(Fast)':>10s} {'speedup':>8s}")
    for name, path in DATASETS:
        d = pd.read_csv(path)
        std = _time(d, False)
        fast = _time(d, True)
        rows.append(dict(dataset=name, rows=len(d), py_standard_ms=round(std, 1),
                         py_fast_ms=round(fast, 1), py_speedup=round(std / fast, 2)))
        print(f"{name:14s} {len(d):7d} {std:13.1f} {fast:10.1f} {std/fast:7.2f}x")
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "bench_py.csv"), index=False)
    print("\nwrote bench_py.csv")


if __name__ == "__main__":
    main()
