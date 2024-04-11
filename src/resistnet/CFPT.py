import sys
import warnings
import numpy as np
from scipy.sparse import diags, eye, csr_matrix, SparseEfficiencyWarning
from scipy.sparse.linalg import spsolve, lgmres, cg
from scipy.sparse import SparseEfficiencyWarning

warnings.simplefilter('ignore', SparseEfficiencyWarning)

def CFPT(Q, R, edge_site_indices):
    N = len(edge_site_indices)
    cfpt_matrix = np.zeros((N, N))

    for i, dest in enumerate(edge_site_indices):
        Q_temp = Q.copy().tolil()
        absorption_factor = R[dest]

        # Get indices and data for the dest row
        for j in range(Q.shape[1]):
            if j != dest:  # Avoid altering the diagonal element at this step
                Q_temp[dest, j] *= (1 - absorption_factor)
        Q_temp = Q_temp.tocsr()
        # Ensure diagonals make row sums to 1
        Q_temp = _set_diags(Q_temp)

        for j, orig in enumerate(edge_site_indices):
            if orig != dest:
                # Extract Qj and qj for current destination
                mask = np.ones(Q.shape[0], dtype=bool)
                mask[dest] = False
                Qj = Q_temp[mask, :][:, mask]
                qj = Q_temp[mask, dest]
                
                # Apply numeric transformations to Qj for solving
                Qj.data *= -1
                Qj.setdiag(Qj.diagonal() + 1)
                
                # Solve for passage times using the adjusted Qj
                try:
                    solution, info = lgmres(Qj, qj.toarray().flatten())
                    if info == 0:
                        adjusted_index = np.where(mask)[0].tolist().index(orig)
                        cfpt_matrix[j, i] = solution[adjusted_index]
                    else:
                        print(f"Convergence issue with dest={dest}, orig={orig}, info={info}")
                        cfpt_matrix[j, i] = np.nan
                except Exception as e:
                    print(f"Solver failed for dest={dest}, orig={orig}: {e}")
                    cfpt_matrix[j, i] = np.nan
    return cfpt_matrix

# def CFPT(Q, R, edge_site_indices):
#     """
#     Calculate conditional first passage times for a Markov chain represented by matrix Q
#     with absorption probabilities R for each state specified in edge_site_indices.
    
#     Parameters:
#     - Q (csr_matrix): Transition probability matrix without absorption (diagonals adjusted).
#     - R (numpy.ndarray): Absorption probabilities for each state.
#     - edge_site_indices (list): Indices of states for which to calculate CFPT.
    
#     Returns:
#     - cfpt_matrix (numpy.ndarray): Matrix of conditional first passage times.
#     """

#     # create Q_scaled
#     # if more efficient can stack these
#     abs_vec = np.ravel(R)
#     # NOTE: I don't think we need the below part currently
#     # for each row (1 - abs_vec[i]) * row_values / np.sum(row_values)

#     # set diagonals
#     Q = _set_diags(Q)

#     # validate Q
#     row_sums = Q.sum(axis=1).A1
#     if not np.allclose(row_sums, 1):
#         raise ValueError("Row sums in Q should sum to 1.")

#     # in-place modifications for downstream steps
#     _opt_prep(Q)

#     # compute cfpt
#     Np = len(edge_site_indices)
#     cfpt_matrix = np.zeros((Np, Np))
    
#     for i, dest in enumerate(edge_site_indices):
#         # Mask to exclude 'dest'
#         mask = np.ones(Q.shape[0], dtype=bool)
#         mask[dest] = False

#         # Create Qj by excluding the 'dest' state
#         Qj = Q[mask, :][:, mask]

#         # Create qj as the 'dest' column, excluding itself
#         qj = Q[mask, dest]

#         for j, orig in enumerate(edge_site_indices):
#             if orig != dest:
#                 try:
#                     passage_times = spsolve(Qj, qj.toarray().flatten())
#                     print(passage_times)
#                     cfpt_matrix[j, i] = passage_times[orig]
#                 except Exception as e:
#                     print(f"Solver failed for dest={dest}, orig={orig}: {e}")
#                     cfpt_matrix[j, i] = np.nan
#     print(cfpt_matrix)
#     sys.exit()
#     return cfpt_matrix

def _set_diags(sm, offset=None):
    """
    Adjusts the diagonal elements of a sparse matrix to ensure that each
    row sums to 1, optionally incorporating an offset subtraction from each
    diagonal element.

    This function modifies the input sparse matrix in-place.

    Parameters:
    - sm (csr_matrix): The sparse matrix whose diagonals are to be adjusted.
    - offset (numpy.ndarray, optional): An array of values to subtract from
    each diagonal. If `None`, no offset is subtracted.

    Returns:
    - None: The matrix `sm` is modified in-place.
    """
    row_sums = sm.sum(axis=1).A1
    d = 1 - row_sums
    if offset is not None:
        d -= offset
    dm = diags(d, 0, shape=sm.shape, format='csr')
    sm += dm
    _validate_row_sums(sm)
    return sm


def _validate_row_sums(sm):
    """
    Validates that each row of a sparse matrix sums to 1.

    Parameters:
    - sm (csr_matrix): The sparse matrix to validate.

    Raises:
    - ValueError: If any row sum does not equal 1.
    """
    row_sums_after = sm.sum(axis=1).A1
    if not np.allclose(row_sums_after, 1):
        raise ValueError("Row sums do not sum to 1.")


def _opt_prep(sm):
    """
    Prepares a scipy.sparse matrix for optimization by negating non-zero
    elements and incrementing diagonal elements by 1, in-place.

    Parameters:
    - sm (scipy.sparse matrix): The sparse matrix to be modified in-place.

    Returns:
    - None
    """
    sm.data = -sm.data
    sm.setdiag(sm.diagonal() + 1)
