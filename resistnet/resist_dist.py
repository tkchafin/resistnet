import itertools
import pandas as pd
import numpy as np

import resistnet.MLPE as mlpe_rga


def parsePairwise(points, inc_matrix, multi, gendist):
    """
    Parse pairwise resistance calculations.

    Args:
        points (dict): A dictionary of points.
        inc_matrix (numpy.ndarray): The incidence matrix.
        multi (list): List of edge resistances.
        gendist (float): The general distance.

    Returns:
        tuple: A tuple containing the effective resistance matrix (r) and the
               result (res).
    """
    r = effectiveResistanceMatrix(points, inc_matrix, multi)
    res = mlpe_rga.MLPE_R(gendist, r, scale=True)
    return r, res


def effectiveResistanceMatrix(points, inc_matrix, edge_resistance):
    """
    Calculate the effective resistance matrix.

    Args:
        points (dict): A dictionary of points.
        inc_matrix (numpy.ndarray): The incidence matrix.
        edge_resistance (list): List of edge resistances.

    Returns:
        numpy.ndarray: The effective resistance matrix.
    """
    r = pd.DataFrame(
        columns=list(points.values()), index=list(points.values())
    )
    inc_row = 0
    edge_resistance = np.array(edge_resistance)

    for ia, ib in itertools.combinations(range(0, len(points)), 2):
        inc = np.array(inc_matrix[inc_row, :])
        d = np.sum(np.multiply(edge_resistance, inc))
        inc_row += 1
        r.loc[list(points.values())[ia], list(points.values())[ib]] = d
        r.loc[list(points.values())[ib], list(points.values())[ia]] = d

    r = r.astype('float64').to_numpy()
    np.fill_diagonal(r, 0.0)

    return r
