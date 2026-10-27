"""Permanent numerical regressions against pinned, current-source R output.

R is required only to regenerate the committed fixtures. Missing references,
different cell/unit support, nonfinite-mask changes and numerical discrepancies
fail these tests; no comparison is skipped because a result is unavailable.
"""
import contextlib
import copy
import gzip
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.att_gt import ATTgt
from r_reference_checks import assert_numeric, decode_numeric, support_order


REFERENCE_DIR = Path(os.environ.get(
    "CSDID_R_REFERENCE_DIR", str(Path(__file__).parent / "r_ref" / "current")))
CASES = json.loads((REFERENCE_DIR / "cases.json").read_text())
FIT_FIELDS = ("att", "se", "c", "l_se", "u_se", "W", "Wpval", "cband")
AGG_FIELDS = ("overall_att", "overall_se", "egt", "att_egt", "se_egt", "crit_val_egt")
KINDS = ("simple", "group", "dynamic", "calendar")
COMPONENT_IF = {"group": "selective_inf_func_g", "dynamic": "dynamic_inf_func_e",
                "calendar": "calendar_inf_func_t"}
OVERALL_IF = {"simple": "simple_att", "group": "selective_inf_func",
              "dynamic": "dynamic_inf_func", "calendar": "calendar_inf_func"}


def _vector(value):
    """Public vector fields can use a one-row/one-column NumPy representation."""
    if value is None:
        return None
    value = np.asarray(value)
    assert value.ndim <= 1 or (value.ndim == 2 and 1 in value.shape)
    return value.reshape(-1)


@pytest.fixture(scope="module")
def references():
    with gzip.open(REFERENCE_DIR / "references.json.gz", "rt") as stream:
        result = json.load(stream)
    names = [case["name"] for case in CASES]
    assert len(names) == len(set(names)) == 25
    assert set(result) == set(names), "R fixture cases are missing, duplicated or unexpected"
    return result


def _unit_ids(obj, case):
    dp = obj.MP["DIDparams"]
    data, internal_id = dp["data"], dp["idname"]
    ids = data[internal_id].unique()
    original_id = case["constructor"]["idname"]
    # True RC replaces the public ID with a synthetic rowid. Recover identity
    # from the retained frame, never from coincidentally matching row positions.
    if internal_id != original_id:
        mapping = data[[internal_id, original_id]].drop_duplicates()
        assert mapping[internal_id].is_unique
        ids = mapping.set_index(internal_id)[original_id].reindex(ids).to_numpy()
    return ids


def _run_case(case, faster_mode):
    data = pd.read_csv(REFERENCE_DIR / case["data"],
                       converters={key: float for key in ("X", "Y", "weight")})
    calls = []
    with contextlib.ExitStack() as stack:
        patch = stack.enter_context(pytest.MonkeyPatch.context())
        stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
        caught = stack.enter_context(warnings.catch_warnings(record=True))
        warnings.simplefilter("always")
        obj = ATTgt(data=data, faster_mode=faster_mode, **case["constructor"])
        if case["fit"].get("bstrap", False):
            multipliers = pd.read_csv(REFERENCE_DIR / case["bootstrap_draws"]["file"],
                                      dtype=str).astype(float)
            draws = multipliers.to_numpy(dtype=float)
            dp = obj.dp
            cluster = case["constructor"].get("clustervar")
            expected_labels = (np.sort(dp["data"][cluster].unique()) if cluster
                               else dp["data"][dp["idname"]].unique())
            np.testing.assert_array_equal(np.asarray(multipliers.columns, dtype=float),
                                          expected_labels)

            def shared_draws(inf_func, biters, pl=False, cores=1):
                influence = np.asarray(inf_func, dtype=float)
                assert draws.shape == (biters, influence.shape[0])
                assert not pl and cores == 1
                calls.append(influence.shape)
                return draws @ influence / influence.shape[0]

            patch.setattr(importlib.import_module("csdid.utils.mboot"),
                          "run_multiplier_bootstrap", shared_draws)
        obj.fit(**case["fit"])
        fit = {key: obj.results[key] for key in ("att", "se", "c", "l_se", "u_se")}
        fit.update({key: obj.MP[key] for key in ("group", "t", "W", "Wpval")})
        fit["inffunc"] = obj.MP["inffunc"]["inffunc"].copy()
        fit["cband"] = obj.MP["DIDparams"]["cband"]
        result = {"unit_ids": _unit_ids(obj, case), "fit": fit, "aggs": {}}
        for kind in KINDS:
            obj.aggte(typec=kind, **case["aggregation"])
            result["aggs"][kind] = copy.deepcopy(
                {key: obj.atte[key] for key in AGG_FIELDS + ("inf_function",)})
        if case["fit"].get("bstrap", False):
            assert calls, "Shared multiplier hook was never called"
        result["warnings"] = [(w.category.__name__, str(w.message)) for w in caught]
        return result


def _require_shape(value, shape, label):
    value = np.asarray(value)
    assert value.shape == shape, label + ": shape {} != {}".format(value.shape, shape)
    return value


def _fit_layout(fit, unit_ids, label):
    """Validate unsliced public arrays before key-based alignment can hide extras."""
    assert set(fit) == set(FIT_FIELDS) | {"group", "t", "inffunc"}
    ids = np.asarray(unit_ids)
    assert ids.ndim == 1, label + ": unit-ID shape must be one-dimensional"
    vectors = {key: _vector(value) for key, value in fit.items() if key != "inffunc"}
    assert vectors["group"] is not None, label + ": absent group labels"
    ncell = len(vectors["group"])
    for key in ("group", "t", "att", "se", "l_se", "u_se"):
        _require_shape(vectors[key], (ncell,), label + "." + key)
    for key in ("c", "W", "Wpval", "cband"):
        if vectors[key] is None:
            assert key in ("W", "Wpval"), label + ": absent " + key
        else:
            _require_shape(vectors[key], (1,), label + "." + key)
    vectors["inffunc"] = _require_shape(fit["inffunc"], (len(ids), ncell),
                                         label + ".inffunc")
    return vectors


def _aggregation_layout(aggregation, kind, nunit, label):
    """Only overall IF vectors may normalize to an n-unit by one matrix."""
    assert set(aggregation) == set(AGG_FIELDS) | {"inf_function"}
    normalized = {key: _vector(aggregation[key]) for key in AGG_FIELDS}
    for key in ("overall_att", "overall_se"):
        _require_shape(normalized[key], (1,), label + "." + key)
    if kind == "simple":
        for key in ("egt", "att_egt", "se_egt", "crit_val_egt"):
            assert normalized[key] is None, label + ": unexpected " + key
        n_event = None
    else:
        assert normalized["egt"] is not None, label + ": absent event labels"
        n_event = len(normalized["egt"])
        for key in ("egt", "att_egt", "se_egt"):
            _require_shape(normalized[key], (n_event,), label + "." + key)
        _require_shape(normalized["crit_val_egt"], (1,), label + ".crit_val_egt")
    expected_keys = {OVERALL_IF[kind]}
    if kind in COMPONENT_IF:
        expected_keys.add(COMPONENT_IF[kind])
    assert set(aggregation["inf_function"]) == expected_keys, label + ": IF fields"
    normalized["inf_function"] = {}
    for key, value in aggregation["inf_function"].items():
        value = np.asarray(value)
        overall = key == OVERALL_IF[kind]
        if overall and value.ndim == 1:
            value = value[:, None]
        normalized["inf_function"][key] = _require_shape(
            value, (nunit, 1 if overall else n_event), label + "." + key)
    return normalized


def _assert_case_output(actual, reference):
    expected_fit = {key: value if key == "cband" else decode_numeric(value)
                    for key, value in reference["fit"].items()}
    actual_fit = _fit_layout(actual["fit"], actual["unit_ids"], "actual fit")
    expected_fit = _fit_layout(expected_fit, reference["unit_ids"], "reference fit")
    cells = support_order(zip(actual_fit["group"], actual_fit["t"]),
                          zip(expected_fit["group"], expected_fit["t"]), "ATT cells")
    units = support_order(actual["unit_ids"], reference["unit_ids"], "IF unit IDs")
    issues = []

    def check(left, right, field, matrix=False):
        try:
            tolerance = {"atol": 1e-5, "rtol": 1e-5} if field in ("W", "Wpval") else {}
            assert_numeric(left, right, field, matrix=matrix, **tolerance)
        except AssertionError as exc:
            issues.append(str(exc))

    for field in FIT_FIELDS:
        value = actual_fit[field]
        if field in ("att", "se", "l_se", "u_se"):
            value = np.asarray(value)[cells]
        check(value, expected_fit[field], field)
    check(actual_fit["inffunc"][np.ix_(units, cells)],
          expected_fit["inffunc"], "inffunc", matrix=True)
    assert set(actual["aggs"]) == set(reference["aggs"]) == set(KINDS)
    for kind in KINDS:
        rr = reference["aggs"][kind]
        assert set(rr) == set(AGG_FIELDS) | {"inf_function"}
        decoded = {field: decode_numeric(rr[field]) for field in AGG_FIELDS}
        decoded["inf_function"] = {field: decode_numeric(value)
                                   for field, value in rr["inf_function"].items()}
        aa = _aggregation_layout(actual["aggs"][kind], kind, len(actual["unit_ids"]),
                                 "actual " + kind)
        decoded = _aggregation_layout(decoded, kind, len(reference["unit_ids"]),
                                      "reference " + kind)
        events = None
        if decoded["egt"] is not None:
            events = support_order(aa["egt"], decoded["egt"], kind + ": events")
        for field in AGG_FIELDS:
            value = aa[field]
            if events is not None and field in ("egt", "att_egt", "se_egt"):
                value = np.asarray(value)[events]
            check(value, decoded[field], kind + "." + field)
        for field, value in aa["inf_function"].items():
            value = np.asarray(value)[units]
            if field == COMPONENT_IF.get(kind):
                assert events is not None
                value = value[:, events]
            check(value, decoded["inf_function"][field],
                  kind + "." + field, matrix=True)
    assert not issues, "\n\n".join(issues)


@pytest.mark.parametrize("faster_mode", [False, True], ids=["standard", "fast"])
@pytest.mark.parametrize("case", CASES, ids=lambda case: case["name"])
def test_current_r_complete_numerical_output(case, faster_mode, references):
    _assert_case_output(_run_case(case, faster_mode), references[case["name"]])


@pytest.fixture(scope="module")
def numerical_gate_example(references):
    case = CASES[0]
    actual, reference = _run_case(case, True), references[case["name"]]
    _assert_case_output(actual, reference)
    return actual, reference


def _numeric_payload_for_corruption(value):
    """Re-encode a deliberately corrupted field without altering its NA masks."""
    value = np.asarray(value, dtype=float)
    flat = value.reshape(-1).copy()
    masks = {key: np.flatnonzero(mask(flat)).tolist()
             for key, mask in (("nan", np.isnan), ("posinf", np.isposinf),
                               ("neginf", np.isneginf))}
    flat[~np.isfinite(flat)] = 0
    return dict(shape=list(value.shape), values=flat.tolist(), **masks)


@pytest.mark.parametrize("side", ["actual", "reference"])
@pytest.mark.parametrize("path,axis", [
    (("fit", "group"), 0), (("fit", "t"), 0),
    (("fit", "att"), 0), (("fit", "se"), 0),
    (("fit", "l_se"), 0), (("fit", "u_se"), 0),
    (("fit", "inffunc"), 0), (("fit", "inffunc"), 1),
    (("aggs", "dynamic", "egt"), 0),
    (("aggs", "dynamic", "att_egt"), 0),
    (("aggs", "dynamic", "se_egt"), 0),
    (("aggs", "dynamic", "crit_val_egt"), 0),
    (("aggs", "dynamic", "overall_att"), 0),
    (("aggs", "dynamic", "inf_function", "dynamic_inf_func_e"), 0),
    (("aggs", "dynamic", "inf_function", "dynamic_inf_func_e"), 1),
    (("aggs", "dynamic", "inf_function", "dynamic_inf_func"), 0),
    (("aggs", "dynamic", "inf_function", "dynamic_inf_func"), 1),
])
def test_full_gate_rejects_extra_elements_before_alignment(numerical_gate_example, side, path, axis):
    """Exercise the full comparator: alignment must not slice away an extra value."""
    actual, reference = copy.deepcopy(numerical_gate_example)
    target = actual if side == "actual" else reference
    for key in path[:-1]:
        target = target[key]
    value = decode_numeric(target[path[-1]]) if side == "reference" else target[path[-1]]
    value = np.asarray(value)
    # Overall IF accepts n-vectors and n-by-one matrices; test extras in either axis.
    if value.ndim == 1 and axis == 1:
        value = value[:, None]
    value = np.atleast_1d(value)
    padding_shape = list(value.shape)
    padding_shape[axis] = 1
    corrupted = np.concatenate((value, np.zeros(padding_shape)), axis=axis)
    target[path[-1]] = (_numeric_payload_for_corruption(corrupted)
                       if side == "reference" else corrupted)
    with pytest.raises(AssertionError, match="shape"):
        _assert_case_output(actual, reference)


@pytest.mark.parametrize("defect", ["shape", "nan", "positive_infinity", "negative_infinity", "value"])
def test_numerical_gate_rejects_corrupted_results(defect):
    reference = np.array([0., 1., np.nan, np.inf, -np.inf])
    corrupted = reference.copy()
    if defect == "shape":
        corrupted = corrupted[:-1]
    elif defect == "nan":
        corrupted[2] = 0.
    elif defect == "positive_infinity":
        corrupted[3] = -np.inf
    elif defect == "negative_infinity":
        corrupted[4] = np.inf
    else:
        corrupted[0] = 1e-3
    with pytest.raises(AssertionError):
        assert_numeric(corrupted, reference, defect)


@pytest.mark.parametrize("actual", [[1], [1, 2, 3], [1, 1], [1, 3]])
def test_support_gate_rejects_missing_extra_duplicate_or_changed_keys(actual):
    with pytest.raises(AssertionError):
        support_order(actual, [1, 2], "cells")


def test_current_r_reference_provenance(references):
    manifest = json.loads((REFERENCE_DIR / "manifest.json").read_text())
    assert manifest["schema_version"] == 1
    expected_commits = {
        "did": "48863b646e657f2aaa5cc9e9d2e7b9d764a6a44e",
        "DRDID": "eaac3b49cd225d0d7706c633858555a5f204e21e",
    }
    for package, commit in expected_commits.items():
        assert manifest["sources"][package]["commit"] == commit
        assert manifest["sources"][package]["tracked_source_clean"] is True
    assert manifest["packages"]["did"] == "2.5.1.904"
    assert manifest["packages"]["DRDID"] == "1.4.0"
    assert manifest["case_count"] == len(CASES) == 25
    expected_artifacts = {
        "panel.csv", "rc.csv", "unbalanced.csv", "multipliers_unit.csv",
        "multipliers_cluster.csv", "cases.json", "references.json.gz",
        "r_mode_checks.json", "sessionInfo.txt", "generate.R",
    }
    assert set(manifest["artifacts_sha256"]) == expected_artifacts
    for name, expected in manifest["artifacts_sha256"].items():
        actual = hashlib.sha256((REFERENCE_DIR / name).read_bytes()).hexdigest()
        assert actual == expected, "R reference artifact changed: " + name
    policy = manifest["tolerance_policy"]
    assert policy["public"] == policy["influence_function"] == {"atol": 1e-7, "rtol": 1e-6}
    assert policy["wald"] == {"atol": 1e-5, "rtol": 1e-5}
    assert policy["missing_masks_and_support"] == "exact"
    checks = json.loads((REFERENCE_DIR / "r_mode_checks.json").read_text())
    assert set(checks) == set(references)
    for name, fields in checks.items():
        assert fields, name + ": missing R fast/standard check"
        assert all(field["same_support"] for field in fields), name
        assert all(field["max_abs"] <= 1e-7 for field in fields), name
