"""Known-DGP and independent Wald checks for the anticipation correction.

Run separately against old/new clean source roots, with one numerical thread::

    python csdid/test_csdid/audit_anticipation_calibration.py \
      --repo /path/to/source --output /tmp/calibration.json

The seeded panel has parallel untreated trends, an anticipation effect of one,
and a treatment effect of two. Ignoring anticipation intentionally violates
the pre-trend null. This one DGP is a diagnostic, not a general coverage proof.
"""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replications", type=int, default=200)
    args = parser.parse_args()
    repo = args.repo.resolve()
    sys.path.insert(0, str(repo))
    import csdid
    from csdid.att_gt import ATTgt
    assert Path(csdid.__file__).resolve().parent.parent == repo
    n, periods, cohort = 400, 8, 6
    treatment = np.arange(n) >= n // 2
    seeds = np.arange(95001, 95001 + args.replications)
    records = []
    for seed in seeds:
        rng = np.random.default_rng(int(seed))
        time = np.arange(1, periods + 1)
        y = (rng.normal(size=(n, 1)) + .2 * time + rng.normal(size=(n, periods))
             + treatment[:, None] * ((time >= 4) & (time < cohort))
             + 2 * treatment[:, None] * (time >= cohort))
        data = pd.DataFrame(dict(id=np.repeat(np.arange(n), periods),
                                 period=np.tile(time, n),
                                 G=np.repeat(cohort * treatment, periods), Y=y.ravel()))
        # Independent weighted-free DID and full empirical covariance, built
        # directly from the DGP arrays rather than the estimator's influence.
        delta = np.column_stack([y[:, 1] - y[:, 0], y[:, 2] - y[:, 1]])
        estimate = delta[treatment].mean(axis=0) - delta[~treatment].mean(axis=0)
        influence = np.empty_like(delta)
        influence[treatment] = 2 * (delta[treatment] - delta[treatment].mean(axis=0))
        influence[~treatment] = -2 * (delta[~treatment] - delta[~treatment].mean(axis=0))
        covariance = influence.T @ influence / n ** 2
        expected_wald = float(estimate @ np.linalg.solve(covariance, estimate))
        window = y[:, 3:5] - y[:, 2, None]
        expected_window = window[treatment].mean(axis=0) - window[~treatment].mean(axis=0)
        for anticipation in (0, 2):
            with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                warnings.simplefilter("ignore")
                fit = ATTgt("Y", "period", "id", "G", data=data,
                            anticipation=anticipation, faster_mode=True).fit(
                                est_method="reg", bstrap=False)
            t = np.asarray(fit.MP["t"])
            att = np.asarray(fit.MP["att"])
            record = dict(seed=int(seed), anticipation=anticipation,
                          p_value=fit.MP["Wpval"], wald=fit.MP["W"],
                          window_att=att[np.isin(t, [4, 5])].tolist(),
                          post_att_mean=float(att[t >= cohort].mean()))
            if anticipation == 2:
                record["independent_wald"] = expected_wald
                record["wald_abs_error"] = abs(fit.MP["W"] - expected_wald)
                record["window_oracle_max_abs_error"] = float(np.max(
                    np.abs(att[np.isin(t, [4, 5])] - expected_window)))
            records.append(record)
    declared = [r for r in records if r["anticipation"] == 2]
    undeclared = [r for r in records if r["anticipation"] == 0]
    summary = dict(
        declared_rejection_rate=float(np.mean([r["p_value"] < .05 for r in declared])),
        undeclared_rejection_rate=float(np.mean([r["p_value"] < .05 for r in undeclared])),
        declared_window_mean=np.mean([r["window_att"] for r in declared], axis=0).tolist(),
        declared_post_mean=float(np.mean([r["post_att_mean"] for r in declared])),
        max_wald_oracle_error=max(r["wald_abs_error"] for r in declared),
        max_window_oracle_error=max(r["window_oracle_max_abs_error"] for r in declared),
    )
    checks = dict(
        declared_size_below_predeclared_15_percent=summary["declared_rejection_rate"] <= .15,
        undeclared_power_above_predeclared_80_percent=summary["undeclared_rejection_rate"] >= .8,
        both_window_effects_within_predeclared_01=np.all(
            np.abs(np.asarray(summary["declared_window_mean"]) - 1) < .1).item(),
        treatment_effect_within_predeclared_01=abs(summary["declared_post_mean"] - 2) < .1,
        wald_matches_independent_formula=summary["max_wald_oracle_error"] < 1e-8,
        window_matches_independent_formula=summary["max_window_oracle_error"] < 1e-10,
    )
    payload = dict(
        schema_version=1, replications=args.replications, units=n, periods=periods,
        true_anticipation_effect=1, true_treatment_effect=2,
        seed_range=[int(seeds[0]), int(seeds[-1])],
        commit=subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        source_sha256={str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (repo / "csdid").rglob("*.py") if "test_csdid" not in p.parts},
        summary=summary, checks=checks, records=records,
        scope="One fixed DGP and analytical covariance oracle; no general finite-sample coverage claim.",
    )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(dict(summary=summary, checks=checks), indent=2))


if __name__ == "__main__":
    main()
