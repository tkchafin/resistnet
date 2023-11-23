import sys
import numpy as np
import scipy.stats


def aggregateDist(stuff, method):
    """
    Aggregate distances using a specified method.

    This function aggregates an array of distances using various statistical
    methods like harmonic mean, arithmetic mean, geometric mean, median, etc.

    Args:
        stuff: A numpy array or list of distances.
        method: A string specifying the aggregation method.

    Returns:
        The aggregated distance value based on the specified method.

    Raises:
        SystemExit: If 'HARM' method encounters a zero distance.
    """
    try:
        if method == "HARM":
            return scipy.stats.hmean(stuff)
        elif method == "ARITH":
            return np.mean(stuff)
        elif method == "GEOM":
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
            return stuff.flat[0]
        elif method == "CV":
            mean = np.mean(stuff)
            return np.std(stuff) / mean if mean != 0 else 0.0
        elif method == "SUM":
            return np.sum(stuff)
    except ValueError as e:
        print(e)
        print("ERROR (DivideByZero): Harmonic mean cannot be calculated using "
              "a zero distance. Try recomputing using the \"ADJHARM\" option.")
        sys.exit(1)


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
