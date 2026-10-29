"""Regression checks for estimation speedups and unchanged bootstrap streams."""
import importlib
import warnings

import numpy as np
import pandas as pd
import pytest

from csdid.attgt_fnc.compute_att_gt2 import compute_att_gt2
from csdid.attgt_fnc.preprocess_did import pre_process_did
from csdid.att_gt import ATTgt
from csdid.utils.bmisc import makeBalancedPanel, multiplier_bootstrap, panel2cs2
from helpers import build_sim_data


standard_module = importlib.import_module("csdid.attgt_fnc.compute_att_gt")


def _legacy_balance(data, idname, tname):
    data = data.sort_values([idname, tname]).reset_index(drop=True)
    nt = len(data[tname].unique())
    return data.groupby(idname).filter(lambda x: len(x) == nt)


def _legacy_panel2cs2(data, yname, idname, tname, assume_sorted=False):
    data = data.sort_values([idname, tname])
    y1 = data.groupby(idname)[yname].shift(-1)
    y0 = data[yname]
    return data.assign(y1=y1, y0=y0, dy=y1-y0).dropna()


@pytest.mark.parametrize("n", [0, 1, 4, 37])
@pytest.mark.parametrize("periods", [1, 2, 5])
@pytest.mark.parametrize("keep", [0.1, 0.6, 1.0])
def test_balance_preserves_exact_frames(n, periods, keep):
    rng = np.random.default_rng(300)
    data = pd.DataFrame({
        "id": np.repeat(np.arange(n, dtype=float)*31-70, periods),
        "t": np.tile(np.arange(periods), n),
        "y": rng.normal(size=n*periods),
    })
    data = data.loc[rng.random(len(data)) < keep].sample(frac=1, random_state=0)
    before = data.copy(deep=True)
    pd.testing.assert_frame_equal(makeBalancedPanel(data, "id", "t"),
                                  _legacy_balance(data, "id", "t"))
    pd.testing.assert_frame_equal(data, before)


@pytest.mark.parametrize("data", [
    pd.DataFrame({"id": [1, 1, 1, 2, 2, np.nan], "t": [1, 1, 2, 1, 2, 1]}),
    pd.DataFrame({"id": pd.Categorical(["a", "b", "a", "b", "c"],
                                       categories=["a", "b", "c", "unused"]),
                  "t": [1, 1, 2, 2, 1]}),
    pd.DataFrame({"id": pd.Series([1, 1, None, None, 2, 2], dtype="Int64"),
                  "t": [1, 2, 1, 2, 1, 2]}),
])
def test_balance_preserves_duplicate_missing_and_categorical_ids(data):
    pd.testing.assert_frame_equal(makeBalancedPanel(data, "id", "t"),
                                  _legacy_balance(data, "id", "t"))


@pytest.mark.parametrize("data", [
    pd.DataFrame({"id": [1, 1, 2, 2, 3, 3], "t": [1, 2]*3,
                  "Y": [0, 1, 2, 3, 4, 5]}),
    pd.DataFrame({"id": [1, 1, 2, 2, 3, 3], "t": [1, 2]*3,
                  "Y": pd.Series([0, 1, 2, None, 4, 5], dtype="Int64")}),
    pd.DataFrame({"id": [1, 1, 2, 2, 3, 3], "t": [1, 2]*3,
                  "Y": pd.Series([0, 1, 2, 3, 4, None], dtype="Float32"),
                  "other": [0, 1, 2, np.nan, 4, 5]}),
    pd.DataFrame({"id": [1, 1, 1, 2, 2, np.nan, np.nan],
                  "t": [1, 1, 2, 1, 2, 1, 2], "Y": np.arange(7.0)}),
    pd.DataFrame({"id": pd.Categorical(["a", "a", "b", "b", "c"]),
                  "t": [1, 2, 1, 2, 2], "Y": np.arange(5.0)}),
    pd.DataFrame({"id": pd.Series([1, 1, 2, 2, None, None], dtype="Int64"),
                  "t": [1, 2]*3, "Y": np.arange(6.0)}),
])
def test_panel_shift_preserves_frames_and_never_crosses_unit_boundaries(data):
    data = data.sample(frac=1, random_state=9)
    before = data.copy(deep=True)
    expected = _legacy_panel2cs2(data, "Y", "id", "t")
    pd.testing.assert_frame_equal(panel2cs2(data, "Y", "id", "t"), expected)
    sorted_data = data.sort_values(["id", "t"])
    pd.testing.assert_frame_equal(
        panel2cs2(sorted_data, "Y", "id", "t", assume_sorted=True), expected)
    pd.testing.assert_frame_equal(data, before)


@pytest.mark.parametrize("periods", [0, 1, 3])
@pytest.mark.parametrize("assume_sorted", [False, True])
def test_panel_shift_rejects_other_than_two_periods(periods, assume_sorted):
    data = pd.DataFrame({"id": [1]*periods, "t": range(periods), "Y": range(periods)})
    with pytest.raises(ValueError, match="panel2cs2 only"):
        panel2cs2(data, "Y", "id", "t", assume_sorted=assume_sorted)


@pytest.mark.parametrize("dtype,offset", [
    ("int64", 2**53), ("uint64", 2**63), ("Int64", 2**53), (object, 2**64),
])
@pytest.mark.parametrize("duplicate_index", [False, True])
def test_panel_shift_keeps_adjacent_large_integer_ids_separate(dtype, offset, duplicate_index):
    ids = pd.Series(np.repeat([offset, offset+1, offset+2], 2), dtype=dtype)
    data = pd.DataFrame({"id": ids, "t": [1, 2]*3, "Y": np.arange(6.)})
    if duplicate_index:
        data.index = [0, 0, 1, 1, 0, 0]
    data = data.sample(frac=1, random_state=19)
    expected = _legacy_panel2cs2(data, "Y", "id", "t")
    assert len(expected) == 3
    pd.testing.assert_frame_equal(panel2cs2(data, "Y", "id", "t"), expected)
    pd.testing.assert_frame_equal(
        panel2cs2(data.sort_values(["id", "t"]), "Y", "id", "t", assume_sorted=True),
        expected)


@pytest.mark.parametrize("dtype", ["Int64", object])
def test_large_integer_shift_excludes_missing_ids(dtype):
    offset = 2**53
    data = pd.DataFrame({
        "id": pd.Series([offset, offset, pd.NA, None, offset+1, offset+1], dtype=dtype),
        "t": [1, 2]*3, "Y": np.arange(6.),
    }).sample(frac=1, random_state=49)
    expected = _legacy_panel2cs2(data, "Y", "id", "t")
    assert len(expected) == 2
    pd.testing.assert_frame_equal(panel2cs2(data, "Y", "id", "t"), expected)


@pytest.mark.parametrize("dtype,offset", [
    ("int64", 2**53), ("uint64", 2**63), ("Int64", 2**53),
])
@pytest.mark.parametrize("faster_mode", [False, True])
@pytest.mark.parametrize("base", ["varying", "universal"])
def test_weighted_public_fit_preserves_large_integer_unit_ids(dtype, offset, faster_mode, base):
    data = build_sim_data(n=100, time_periods=4, seed=24)
    data["weight"] = 0.2+(data["id"] % 11)/7
    data = data.sample(frac=1, random_state=39)
    huge_ids = data.copy(deep=True)
    huge_ids["id"] = pd.Series([int(value)+offset for value in data["id"]],
                                index=data.index, dtype=dtype)

    def fit(frame):
        return ATTgt(yname="Y", tname="period", idname="id", gname="G", data=frame,
                     xformla="~X", weights_name="weight", faster_mode=faster_mode).fit(
                         est_method="dr", base_period=base, bstrap=False)

    expected, actual = fit(data), fit(huge_ids)
    assert actual.dp["n"] == expected.dp["n"] == 100
    for field in ("att", "se", "c", "l_se", "u_se"):
        np.testing.assert_array_equal(actual.results[field], expected.results[field])
    np.testing.assert_array_equal(actual.MP["inffunc"]["inffunc"],
                                  expected.MP["inffunc"]["inffunc"])
    np.testing.assert_equal(actual.MP["W"], expected.MP["W"])
    np.testing.assert_equal(actual.MP["Wpval"], expected.MP["Wpval"])


def _preprocessed_data(fix_weights=None, control_group="notyettreated"):
    data = build_sim_data(n=100, time_periods=4, seed=17)
    data["id"] = 7*data["id"]-901
    data["X"] += 0.07*data["period"]
    data["category"] = np.where(data["id"] % 3 == 0, "a", "b")
    data["weight"] = 0.2+((data["id"] % 11)+data["period"])/7
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pre_process_did(
            yname="Y", tname="period", idname="id", gname="G",
            data=data.sample(frac=1, random_state=7), xformla="~X+C(category)",
            control_group=control_group, weights_name="weight", fix_weights=fix_weights)


def _run_standard(dp, method, base):
    with warnings.catch_warnings(record=True) as found:
        warnings.simplefilter("always")
        output, inf = standard_module.compute_att_gt(
            dict(dp, data=dp["data"].copy(deep=True)), est_method=method, base_period=base)
    return output, inf, [(type(w.message), str(w.message)) for w in found]


@pytest.mark.parametrize("method", ["dr", "ipw", "reg"])
@pytest.mark.parametrize("base", ["varying", "universal"])
@pytest.mark.parametrize("fix_weights", [None, "base_period", "first_period"])
@pytest.mark.parametrize("control_group", ["nevertreated", "notyettreated"])
def test_standard_reshape_is_bit_identical(monkeypatch, method, base, fix_weights, control_group):
    dp = _preprocessed_data(fix_weights, control_group)
    actual = _run_standard(dp, method, base)
    monkeypatch.setattr(standard_module, "panel2cs2", _legacy_panel2cs2)
    expected = _run_standard(dp, method, base)
    for key in ("group", "year", "post"):
        assert actual[0][key] == expected[0][key]
    np.testing.assert_array_equal(actual[0]["att"], expected[0]["att"])
    np.testing.assert_array_equal(actual[1], expected[1])
    assert actual[2] == expected[2]


@pytest.mark.parametrize("shuffle", [False, True])
@pytest.mark.parametrize("base", ["varying", "universal"])
@pytest.mark.parametrize("fix_weights", [None, "base_period", "first_period"])
def test_fast_arrays_preserve_alignment_layout_and_influence_positions(shuffle, base, fix_weights):
    dp = _preprocessed_data(fix_weights)
    data = dp["data"].copy(deep=True)
    data["Y"] = data["id"]*100+data["period"]
    if shuffle:
        data = data.sample(frac=1, random_state=2)
    dp["data"] = data
    first_ids = data.loc[data.period.eq(dp["tlist"][0]), "id"].to_numpy()
    positions = pd.Series(np.arange(len(first_ids)), index=first_ids)
    records = []

    def estimator(ypost, ypre, D, i_weights=None, covariates=None):
        ids = np.floor_divide(ypre, 100).astype(int)
        pret = int(ypre[0] % 100)
        tn = int(ypost[0] % 100)
        early = min(pret, tn)
        selected = data.loc[data.period.eq(early)].set_index("id").reindex(ids)
        g = selected.loc[D == 1, "G"].unique().item()
        cohorts = data.loc[data.period.eq(dp["tlist"][0])].set_index("id").loc[first_ids, "G"].to_numpy()
        mask = (cohorts == g) | (cohorts == 0) | ((cohorts > max(pret, tn)) & (cohorts != g))
        np.testing.assert_array_equal(ids, first_ids[mask])
        np.testing.assert_array_equal(ypost, ids*100+tn)
        np.testing.assert_array_equal(D, (selected.G == g).to_numpy(dtype=float))
        expected_cov = selected[dp["xcov_cols"]].to_numpy()
        np.testing.assert_array_equal(covariates, expected_cov)
        assert covariates.flags.f_contiguous
        weight_period = early
        if fix_weights == "first_period":
            weight_period = dp["tlist"][0]
        elif fix_weights == "base_period":
            weight_period = max(t for t in dp["tlist"] if t < g)
        expected_w = data.loc[data.period.eq(weight_period)].set_index("id").reindex(ids).w.to_numpy()
        np.testing.assert_array_equal(i_weights, expected_w)
        influence = np.arange(len(ids), dtype=float)+0.125
        expected_inf = np.zeros(dp["n"])
        expected_inf[positions.reindex(ids).to_numpy()] = (dp["n"]/len(ids))*influence
        records.append(expected_inf)
        return float(len(records)), influence

    output, influence = compute_att_gt2(dp, est_method=estimator, base_period=base)
    assert len(records) > 0
    for estimate, actual_inf in zip(output["att"], influence):
        expected_inf = np.zeros(dp["n"]) if estimate == 0 else records[int(estimate)-1]
        np.testing.assert_array_equal(actual_inf, expected_inf)


@pytest.mark.parametrize("seed_kind", ["global", "integer", "sequence", "generator"])
@pytest.mark.parametrize("shape", [(1, 1), (2, 4), (7, 5), (8, 0)])
@pytest.mark.parametrize("biters", [0, 1, 73])
def test_bootstrap_matches_choice_outputs_and_rng_state(seed_kind, shape, biters):
    influence = np.random.default_rng(9).normal(size=shape)
    np.random.seed(123)
    if seed_kind == "global":
        expected = np.random.choice([1, -1], size=(biters, shape[0])) @ influence / shape[0]
        expected_next = np.random.randint(0, 100000, size=31)
        np.random.seed(123)
        actual = multiplier_bootstrap(influence, biters)
        np.testing.assert_array_equal(np.random.randint(0, 100000, size=31), expected_next)
    else:
        expected_rng = np.random.default_rng(17)
        expected = expected_rng.choice([1, -1], size=(biters, shape[0])) @ influence / shape[0]
        seed = {"integer": 17, "sequence": np.random.SeedSequence(17),
                "generator": np.random.default_rng(17)}[seed_kind]
        before = np.random.get_state()
        actual = multiplier_bootstrap(influence, biters, seed=seed)
        after = np.random.get_state()
        for old, new in zip(before, after):
            np.testing.assert_array_equal(new, old)
        if seed_kind == "generator":
            np.testing.assert_array_equal(seed.integers(0, 100000, size=31),
                                          expected_rng.integers(0, 100000, size=31))
    np.testing.assert_array_equal(actual, expected)


def test_bootstrap_retains_500_million_element_chunk_boundary(monkeypatch):
    """Exercise the real chunk arithmetic without allocating enormous arrays."""
    calls = []

    class Influence:
        shape = (250_000_000, 1)

    class Draws:
        def __init__(self, size):
            self.size = size
            self.value = 0

        def __imul__(self, value):
            self.value *= value
            return self

        def __iadd__(self, value):
            self.value += value
            return self

        def __matmul__(self, influence):
            assert influence is source
            assert self.value == 1
            return np.full((self.size[0], 1), self.size[1], dtype=float)

    def draw(low, high, size):
        assert (low, high) == (0, 2)
        calls.append(size)
        return Draws(size)

    source = Influence()
    monkeypatch.setattr(np.random, "randint", draw)
    np.testing.assert_array_equal(multiplier_bootstrap(source, 5), np.ones((5, 1)))
    assert calls == [(2, 250_000_000), (2, 250_000_000), (1, 250_000_000)]
