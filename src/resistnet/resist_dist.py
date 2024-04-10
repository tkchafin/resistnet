import itertools
import pandas as pd
import sys
import numpy as np

import resistnet.MLPE as mlpe_rga
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
    res = mlpe_rga.MLPE_R(gendist, r, scale=True)
    return r, res


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

def conditionalFirstPassTime(Q, R, sites_i, gendist):

    # get cfpt matrix
    cfpt = cfpt_samc.CFPT(Q, R, sites_i)
    cfpt = np.array(cfpt)

    # fit MLPE
    res = mlpe_rga.MLPE_R(gendist, cfpt, scale=True)
    return cfpt, res


