import numpy as np
import resistnet.MLPE as mlpe
import resistnet.CFPT as cfpt_samc


def parsePairwise(points, inc_matrix, multi, gendist):
    """
    Parse pairwise resistance calculations.

    Args:
        points (dict): A dictionary of points.
        inc_matrix (numpy.ndarray): The incidence matrix.
        multi (pd.Series): List of edge resistances.
        gendist (float): The general distance.

    Returns:
        tuple: A tuple containing the effective resistance matrix (r) and the
               result (res).
    """
    r = effectiveResistanceMatrix(points, inc_matrix, multi)
    # res = mlpe_rga.MLPE_R(gendist, r, scale=True)
    res2 = mlpe.MLPE(gendist, r, scale=True)

    return r, res2


def effectiveResistanceMatrix(points, inc_matrix, edge_resistance):
    """
    Calculate the effective resistance matrix.

    Args:
        points (dict): A dictionary of points.
        inc_matrix (numpy.ndarray): The incidence matrix.
        edge_resistance (pd.Series): List of edge resistances.

    Returns:
        numpy.ndarray: The effective resistance matrix.
    """
    num_points = len(points)
    if inc_matrix.shape[0] != (num_points * (num_points - 1)) // 2:
        raise ValueError(
            "Incorrect dimensions of inc_matrix for the given number of points"
        )

    # Rows of inc_matrix follow itertools.combinations order, which is the
    # same as the upper-triangle (row-major) order
    d = np.asarray(inc_matrix) @ np.asarray(edge_resistance, dtype='float64')
    r = np.zeros((num_points, num_points))
    ia, ib = np.triu_indices(num_points, 1)
    r[ia, ib] = d
    r[ib, ia] = d

    return r


def conditionalFirstPassTime(Q, sites_i, gendist, absorption=None,
                             combine="mean"):
    """
    Conditional first passage times (CFPT) among sites, fitted to genetic
    distances with MLPE.

    Genetic distances are symmetric, so the model combines the two
    directions between each pair of sites (their mean or their maximum).

    Args:
        Q (scipy.sparse matrix): Transient transition matrix.
        sites_i (list): State index of each site.
        gendist (numpy.ndarray): Genetic distance matrix.
        absorption (array-like, optional): Total absorption of each state
            (see CFPT.cond_passage).
        combine (str): How the two directions between a pair of sites are
            combined: "mean" or "max".

    Returns:
        tuple: The CFPT matrix (row = origin, column = destination), and the
               MLPE results (None if some passage times are not finite).
    """
    cfpt = cfpt_samc.CFPT(Q, sites_i, absorption)
    if not np.all(np.isfinite(cfpt)):
        return cfpt, None
    res = mlpe.MLPE(gendist, combine_directions(cfpt, combine), scale=True)
    return cfpt, res


def combine_directions(T, how="mean"):
    """
    Combines the passage times of the two directions between each pair of
    sites into a symmetric distance: their mean, or their maximum.

    Args:
        T (numpy.ndarray): Passage times (row = origin, column =
                           destination).
        how (str): "mean" or "max".

    Returns:
        numpy.ndarray: The symmetric distances.
    """
    T = np.asarray(T, dtype=float)
    if how == "mean":
        return (T + T.T) / 2.0
    if how == "max":
        return np.maximum(T, T.T)
    raise ValueError(f"combine must be 'mean' or 'max', not {how!r}")
