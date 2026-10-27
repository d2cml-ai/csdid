"""Scale-safe sampling-weight reductions with the ordinary path preserved."""

import numpy as np


def _low_precision_mean_is_unsafe(values, mean):
    """Float16/32 can lose a mean that is representable in R's doubles."""
    return (values.dtype.kind == "f" and values.dtype.itemsize < 8
            and np.isfinite(mean) and 0 <= mean < np.finfo(values.dtype).tiny
            and values.size and np.all(np.isfinite(values))
            and np.all(values >= 0) and np.any(values > 0))


def _weight_mean(weights):
    """Keep valid weight means representable in the double-precision contract."""
    values = np.asarray(weights)
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        result = np.mean(values)
    if _low_precision_mean_is_unsafe(values, result):
        return np.mean(values.astype(np.float64))
    if (not np.isfinite(result) and values.size
            and np.all(np.isfinite(values)) and np.all(values >= 0)):
        largest = np.max(values)
        if largest > 0:
            return largest * np.mean(values / largest)
    return result


def _normalize_weights(weights):
    values = np.asarray(weights)
    mean = _weight_mean(values)
    if _low_precision_mean_is_unsafe(values, mean):
        # NumPy 1.x can cast a scalar denominator back to the array's dtype.
        # Promote the array too so a positive double mean cannot become zero.
        values = values.astype(np.float64)
    return values / mean


def _rescale_weights_if_mean_overflows(weights):
    """Protect legacy external kernels on the selected 2x2 sample only.

    Repair sum overflow and a float16/32 mean in the subnormal range.
    Retaining raw weights elsewhere avoids rounding a small cohort's weights
    to zero because an unrelated cohort has a much larger weight scale.
    """
    values = np.asarray(weights)
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        ordinary_mean = np.mean(values)
    if (_low_precision_mean_is_unsafe(values, ordinary_mean) or
            (np.isinf(ordinary_mean) and values.size
             and np.all(np.isfinite(values)) and np.all(values >= 0))):
        return _normalize_weights(values)
    return weights


def _cohort_proportions(masses, sequential=False):
    """Normalize a selected set of cohort masses, preserving ordinary rounding."""
    values = np.asarray(masses)
    reduction = sum if sequential else np.sum
    try:
        with np.errstate(all="raise"):
            return values / reduction(values)
    except FloatingPointError:
        if (values.size and np.all(np.isfinite(values))
                and np.all(values >= 0) and np.max(values) > 0):
            with np.errstate(under="ignore"):
                scaled = values / np.max(values)
                return scaled / reduction(scaled)
        return values / reduction(values)


def _cohort_average(values, masses):
    """Avoid overflowing or underflowing a weighted-average numerator."""
    values, masses = np.asarray(values), np.asarray(masses)
    try:
        with np.errstate(all="raise"):
            return np.sum(values * masses) / np.sum(masses)
    except FloatingPointError:
        if (masses.size and np.all(np.isfinite(values))
                and np.all(np.isfinite(masses)) and np.all(masses >= 0)
                and np.max(masses) > 0):
            return np.sum(values * _cohort_proportions(masses))
        return np.sum(values * masses) / np.sum(masses)
