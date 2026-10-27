"""Scale-safe sampling-weight reductions with the ordinary path preserved."""

import numpy as np


def _weight_mean(weights):
    """Match a finite non-negative weight mean without overflowing its sum."""
    values = np.asarray(weights)
    with np.errstate(over="ignore", invalid="ignore"):
        result = np.mean(values)
    if (not np.isfinite(result) and values.size
            and np.all(np.isfinite(values)) and np.all(values >= 0)):
        largest = np.max(values)
        if largest > 0:
            return largest * np.mean(values / largest)
    return result


def _normalize_weights(weights):
    values = np.asarray(weights)
    return values / _weight_mean(values)


def _rescale_weights_if_mean_overflows(weights):
    """Protect legacy external kernels on the selected 2x2 sample only.

    Retaining raw weights elsewhere avoids rounding a small cohort's weights
    to zero because an unrelated cohort has a much larger weight scale.
    """
    values = np.asarray(weights)
    with np.errstate(over="ignore", invalid="ignore"):
        ordinary_mean = np.mean(values)
    if (np.isinf(ordinary_mean) and values.size
            and np.all(np.isfinite(values)) and np.all(values >= 0)):
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
