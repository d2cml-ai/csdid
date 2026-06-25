"""Shared helper functions for csdid tests."""
import numpy as np
import pandas as pd


def build_sim_data(n=1000, time_periods=4, te=1.0, te_e=None, seed=9142024):
    """Build simulated DiD panel data with known treatment effects.

    Groups: 0 (never-treated), 2, 3, ..., time_periods.
    Treatment effect at exposure *e* is ``te_e[e]`` when *te_e* is given,
    otherwise the constant *te*.
    """
    rng = np.random.default_rng(seed)
    groups = [0] + list(range(2, time_periods + 1))
    n_per = n // len(groups)
    unit_g = np.concatenate([np.full(n_per, g, dtype=int) for g in groups])
    rem = n - len(unit_g)
    if rem > 0:
        unit_g = np.append(unit_g, np.zeros(rem, dtype=int))
    n_act = len(unit_g)
    ids = np.arange(1, n_act + 1)
    X = rng.standard_normal(n_act)
    alpha = rng.standard_normal(n_act) * 0.3
    cluster = (ids % 10) + 1

    frames = []
    for t in range(1, time_periods + 1):
        eps = rng.standard_normal(n_act)
        Y = alpha + 0.3 * X + 0.1 * t + eps
        for g in groups:
            if g > 0 and t >= g:
                mask = unit_g == g
                e = t - g
                eff = te_e[min(e, len(te_e) - 1)] if te_e is not None else te
                Y = Y.copy()
                Y[mask] += eff
        frames.append(pd.DataFrame({
            "id": ids, "period": t, "G": unit_g, "Y": Y,
            "X": X, "cluster": cluster,
        }))
    return pd.concat(frames, ignore_index=True)
