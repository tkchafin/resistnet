"""
Conditional first passage times for spatial absorbing Markov chains (SAMC).

Follows the samc R package (Marx et al. 2020; Fletcher et al. 2019): the
transient part Q of the transition matrix P = [[Q, R], [0, I]] is built from
movement weights and absorption, and conditional first passage times are
computed as in ``samc::cond_passage``.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu


def transition_matrix(W, absorption, fidelity=None):
    """
    Builds the transient transition matrix Q of an absorbing Markov chain.

    As in samc, movement probabilities are the row-normalised movement
    weights, scaled by the probability of neither being absorbed nor staying
    put: Q = diag(1 - absorption - fidelity) D^-1 W + diag(fidelity), where D
    holds the row sums of W. Each row of Q then sums to 1 - absorption.

    Args:
        W (scipy.sparse matrix): Non-negative movement weights (e.g.
            conductances); W[i, j] > 0 if a move from state i to j is
            possible. The diagonal is ignored.
        absorption (array-like): Total absorption probability of each state.
        fidelity (array-like, optional): Probability of staying in each
            state for a time step. Defaults to 0.

    Returns:
        scipy.sparse.csr_matrix: The matrix Q.
    """
    W = sp.csr_matrix(W, dtype=float)
    W.setdiag(0)
    W.eliminate_zeros()
    n = W.shape[0]
    absorption = np.asarray(absorption, dtype=float).ravel()
    if fidelity is None:
        fidelity = np.zeros(n)
    fidelity = np.asarray(fidelity, dtype=float).ravel()

    if W.data.size and W.data.min() < 0:
        raise ValueError("Movement weights must be non-negative")
    move = 1.0 - absorption - fidelity
    if np.any(absorption < 0) or np.any(fidelity < 0) or np.any(move < 0):
        raise ValueError("Absorption and fidelity must be in [0, 1] and sum "
                         "to at most 1")
    row_sums = np.asarray(W.sum(axis=1)).ravel()
    if np.any((row_sums == 0) & (move > 0)):
        raise ValueError("Some states have no possible moves")

    scale = np.divide(move, row_sums, out=np.zeros(n), where=row_sums > 0)
    Q = sp.diags(scale) @ W + sp.diags(fidelity)
    return Q.tocsr()


def cond_passage(Q, dest):
    """
    Mean first passage time to a destination state, conditional on reaching
    it (equivalent to ``samc::cond_passage(samc, dest = dest)``).

    With dest removed from Q, b = (I - Q)^-1 q is the probability of reaching
    dest from each state (q being the one-step probabilities of moving into
    dest), and the conditional passage times are t = (I - Q)^-1 b / b.

    States that are very unlikely to reach dest have tiny b, where the LU
    solves lose relative accuracy; both solves are therefore iteratively
    refined, which restores full precision since I - Q is an M-matrix.

    Args:
        Q (scipy.sparse matrix): Transient transition matrix.
        dest (int): Index of the destination state.

    Returns:
        numpy.ndarray: Conditional passage time from every state to dest
        (0 for dest itself).
    """
    Q = sp.csr_matrix(Q)
    n = Q.shape[0]
    keep = np.arange(n) != dest
    Qj = Q[keep][:, keep]
    qj = Q[keep][:, [dest]].toarray().ravel()

    A = sp.csc_matrix(sp.identity(n - 1) - Qj)
    lu = splu(A)
    b = _refine(lu, A, qj, lu.solve(qj))
    t = _refine(lu, A, b, lu.solve(b)) / b

    out = np.zeros(n)
    out[keep] = t
    return out


def _refine(lu, A, rhs, x, max_steps=10):
    """
    Iterative refinement of a solution x of A x = rhs, reusing the LU
    factorisation, until the correction is negligible relative to each
    element of x.
    """
    for _ in range(max_steps):
        dx = lu.solve(rhs - A @ x)
        x = x + dx
        if np.all(np.abs(dx) <= 1e-15 * np.abs(x)):
            break
    return x


def CFPT(Q, edge_site_indices):
    """
    Pairwise conditional first passage times among a set of states.

    Args:
        Q (scipy.sparse matrix): Transient transition matrix.
        edge_site_indices (list): State index of each site.

    Returns:
        numpy.ndarray: Matrix T where T[i, j] is the conditional passage
        time from site i to site j (0 on the diagonal).
    """
    sites = np.asarray(edge_site_indices)
    N = len(sites)
    cfpt_matrix = np.zeros((N, N))
    for j, dest in enumerate(sites):
        cfpt_matrix[:, j] = cond_passage(Q, dest)[sites]
    np.fill_diagonal(cfpt_matrix, 0.0)
    return cfpt_matrix
