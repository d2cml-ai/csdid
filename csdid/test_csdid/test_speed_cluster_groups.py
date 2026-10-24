"""Numerical and RNG parity for bootstrap cluster grouping."""
import importlib

import numpy as np
import pandas as pd
import pytest

from csdid.utils.mboot import _cluster_sums_by_label


def _reference_cluster_sums(inf_func, labels):
    """The former per-cluster mask, including its reduction order and dtype."""
    unique = np.unique(labels)
    result = np.zeros((len(unique), inf_func.shape[1]))
    for i, label in enumerate(unique):
        result[i] = inf_func[labels == label].sum(axis=0)
    return result


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("layout", ["C", "F"])
@pytest.mark.parametrize("label_kind", ["integer", "string", "missing", "datetime"])
@pytest.mark.parametrize("cluster_count", [5, 16, 17])
def test_cluster_sums_preserve_order_dtype_and_missing_labels(
        dtype, layout, label_kind, cluster_count):
    rng = np.random.default_rng(861)
    codes = rng.choice(cluster_count, 503)
    labels = codes * 1_000_003 - 17
    if label_kind == "string":
        labels = np.array([f"cluster-{code}" for code in codes])
    elif label_kind == "missing":
        labels = labels.astype(float)
        labels[::11] = np.nan
    elif label_kind == "datetime":
        labels = np.datetime64("2000-01-01") + codes.astype("timedelta64[D]")
        labels[::11] = np.datetime64("NaT", "D")

    # Cancellation makes a different accumulation order observable, especially
    # when a float32 input would otherwise be promoted before summation.
    values = rng.normal(size=(len(labels), 7))
    values[::3] *= 1e10
    values[1::3] *= -1e10
    values = np.array(values, dtype=dtype, order=layout)
    expected = _reference_cluster_sums(values, labels)
    actual = _cluster_sums_by_label(values, labels)
    np.testing.assert_array_equal(actual, expected)
    assert actual.dtype == expected.dtype


@pytest.mark.parametrize("n_columns", [1, 5])
def test_cluster_sums_empty_input_and_singleton_clusters(n_columns):
    empty = np.empty((0, n_columns), dtype=np.float32)
    np.testing.assert_array_equal(
        _cluster_sums_by_label(empty, np.array([], dtype=int)),
        _reference_cluster_sums(empty, np.array([], dtype=int)),
    )
    labels = np.array([99, -2, 4, 17])
    values = np.arange(len(labels) * n_columns).reshape(len(labels), n_columns)
    np.testing.assert_array_equal(
        _cluster_sums_by_label(values, labels),
        _reference_cluster_sums(values, labels),
    )


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_object_labels_with_nan_preserve_legacy_unique_sequence(seed):
    rng = np.random.default_rng(seed)
    labels = rng.integers(0, 30, 500).astype(object)
    labels[::11] = np.nan
    values = rng.normal(size=(500, 4))
    np.testing.assert_array_equal(
        _cluster_sums_by_label(values, labels),
        _reference_cluster_sums(values, labels),
    )


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_clustered_bootstrap_outputs_and_rng_match_previous_grouping(monkeypatch, dtype):
    module = importlib.import_module("csdid.utils.mboot")
    rng = np.random.default_rng(910)
    labels = rng.integers(0, 23, size=173)
    values = np.array(rng.normal(size=(len(labels), 4)), dtype=dtype, order="F")
    data = pd.DataFrame({"id": np.arange(len(labels)), "t": 0, "cluster": labels})
    params = dict(
        data=data, idname="id", tname="t", clustervars="cluster",
        biters=87, alp=0.05, panel=True, true_repeated_cross_sections=False,
    )
    np.random.seed(632)
    actual = module.mboot(values, params)
    actual_next_draw = np.random.random(8)

    monkeypatch.setattr(module, "_cluster_sums_by_label", _reference_cluster_sums)
    np.random.seed(632)
    expected = module.mboot(values, params)
    expected_next_draw = np.random.random(8)
    for key in ("bres", "V", "se", "crit_val"):
        np.testing.assert_array_equal(actual[key], expected[key])
    np.testing.assert_array_equal(actual_next_draw, expected_next_draw)
