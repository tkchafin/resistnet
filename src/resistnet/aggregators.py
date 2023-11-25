import numpy as np
import scipy.stats
import pandas as pd


def aggregateDist(stuff, method):
    """
    Aggregate distances using a specified method.

    Args:
        stuff: A numpy array or list of distances.
        method: A string specifying the aggregation method.

    Returns:
        The aggregated distance value based on the specified method.

    Raises:
        ValueError: For invalid input types, zero values in 'HARM',
                    unsupported methods, or other issues with the input data.
    """
    # Validate input types
    if not isinstance(stuff, (list, np.ndarray, pd.Series)):
        raise ValueError("Input 'stuff' must be a list or numpy array.")
    if not isinstance(method, str):
        raise ValueError("Input 'method' must be a string.")

    try:
        if method == "HARM":
            if 0 in stuff:
                raise ValueError(
                    "Zero values found. Harmonic mean cannot be calculated. \
                        Try 'ADJHARM'."
                )
            return scipy.stats.hmean(stuff)
        elif method == "ARITH":
            return np.mean(stuff)
        elif method == "GEOM":
            if any(x < 0 for x in stuff):
                raise ValueError(
                    "Negative values cannot be used for geometric mean"
                )
            return scipy.stats.mstats.gmean(stuff)
        elif method == "MEDIAN":
            return np.median(stuff)
        elif method == "MAX":
            return np.max(stuff)
        elif method == "MIN":
            return np.min(stuff)
        elif method == "ADJHARM":
            return adjustedHarmonicMean(stuff)
        elif method == "SD":
            return np.std(stuff)
        elif method == "VAR":
            return np.var(stuff)
        elif method == "FIRST":
            if len(stuff) == 0:
                raise ValueError("Input 'stuff' is empty.")
            return stuff.flat[0]
        elif method == "CV":
            mean = np.mean(stuff)
            return np.std(stuff) / mean if mean != 0 else 0.0
        elif method == "SUM":
            return np.sum(stuff)
        else:
            raise ValueError(f"Unsupported aggregation method: {method}")

    except ValueError as e:
        raise ValueError(f"Error in aggregateDist: {e}")

    except Exception as e:
        raise ValueError(f"Unexpected error: {e}")


def adjustedHarmonicMean(stuff):
    """
    Compute an adjusted harmonic mean for non-positive values.

    This function calculates the harmonic mean for a dataset, correcting
    for non-positive values.

    Args:
        stuff: A numpy array or list of values.

    Returns:
        float: The adjusted harmonic mean of the values.
    """
    s = np.array(stuff)
    positive_vals = s[s > 0.0]
    num_positive = len(positive_vals)
    total_count = len(s)
    if num_positive == 0:
        return 0.0

    inverse_sum = np.sum([1.0 / x for x in positive_vals])
    adjustment_factor = (num_positive - (total_count - num_positive))
    mean_inverse = inverse_sum / adjustment_factor
    adjusted_mean = (1.0 / mean_inverse) * (num_positive / total_count)

    return adjusted_mean
