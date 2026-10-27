"""Strict support, shape and numerical checks for committed R references."""
import numpy as np


def decode_numeric(value):
    """Decode explicit R nonfinite masks without conflating NA and infinity."""
    if value is None:
        return None
    assert set(value) == {"shape", "values", "nan", "posinf", "neginf"}
    shape = tuple(value["shape"])
    assert all(isinstance(n, int) and n >= 0 for n in shape)
    flat = np.asarray(value["values"], dtype=float)
    assert flat.ndim == 1 and len(flat) == int(np.prod(shape))
    assert np.isfinite(flat).all(), "Nonfinite reference values must use explicit masks"
    used = set()
    for key, replacement in (("nan", np.nan), ("posinf", np.inf), ("neginf", -np.inf)):
        indices = value[key]
        assert all(isinstance(i, int) and 0 <= i < len(flat) for i in indices)
        assert len(indices) == len(set(indices))
        assert not used.intersection(indices), "Overlapping nonfinite masks"
        used.update(indices)
        flat[indices] = replacement
    return flat.reshape(shape)


def support_order(actual, expected, label):
    """Return actual positions in expected order, rejecting missing/extra keys."""
    actual, expected = list(actual), list(expected)
    assert len(actual) == len(set(actual)), label + ": duplicate actual keys"
    assert len(expected) == len(set(expected)), label + ": duplicate reference keys"
    assert set(actual) == set(expected), label + ": different support"
    positions = {key: i for i, key in enumerate(actual)}
    return np.asarray([positions[key] for key in expected], dtype=int)


def assert_numeric(actual, expected, label, atol=1e-7, rtol=1e-6, matrix=False):
    """Compare every element, requiring identical shape and all nonfinite masks."""
    if actual is None or expected is None:
        assert actual is None and expected is None, label + ": absent-field mismatch"
        return
    actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    if matrix:
        if actual.ndim == 1:
            actual = actual.reshape(-1, 1)
        if expected.ndim == 1:
            expected = expected.reshape(-1, 1)
    else:
        # Public scalar fields may be scalars or length-one vectors. No other
        # broadcasting or dimensional coercion is allowed.
        actual, expected = np.atleast_1d(actual), np.atleast_1d(expected)
    assert actual.shape == expected.shape, (
        "{}: shape {} != {}".format(label, actual.shape, expected.shape))
    for name, mask in (("finite", np.isfinite), ("nan", np.isnan),
                       ("posinf", np.isposinf), ("neginf", np.isneginf)):
        np.testing.assert_array_equal(mask(actual), mask(expected),
                                      err_msg=label + ": " + name + " mask")
    finite = np.isfinite(expected)
    np.testing.assert_allclose(actual[finite], expected[finite], atol=atol,
                               rtol=rtol, err_msg=label)
