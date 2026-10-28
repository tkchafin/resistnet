import numpy as np
import resistnet.MLPE as mlpe


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


# import resistnet.CFPT as cfpt_samc
# def conditionalFirstPassTime(Q, R, sites_i, gendist,
#                              rtol=0.00001, max_iter=1000, max_fail=1,
#                              solver="iterative"):

#     # get cfpt matrix
#     cfpt = cfpt_samc.CFPT(Q, R, sites_i, rtol, max_iter, max_fail, solver)

#     # fit MLPE
#     if cfpt is not None:
#         cfpt_a = np.array(cfpt)
#         if np.all(cfpt_a == cfpt_a[0]) or np.all(np.isnan(cfpt_a)):
#             return None, float('-inf')
#         try:
#             res = mlpe.MLPE(gendist, cfpt_a, scale=True)
#             return cfpt, res
#         except Exception as e:
#             # NOTE: This is often caused by most of the CFPT values
#             # being 0
#             print("Unexpected error fitting MLPE:", e)
#             return None, float('-inf')
#     else:
#         return None, float('-inf')
