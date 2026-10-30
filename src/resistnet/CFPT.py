"""
Conditional first passage times for spatial absorbing Markov chains (SAMC).

Follows the samc R package (Marx et al. 2020; Fletcher et al. 2019): the
transient part Q of the transition matrix P = [[Q, R], [0, I]] is built from
movement weights and absorption, and conditional first passage times are
computed as in ``samc::cond_passage``.
"""
import heapq

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


def cond_passage(Q, dest, absorption=None):
    """
    Mean first passage time to a destination state, conditional on reaching
    it (equivalent to ``samc::cond_passage(samc, dest = dest)``).

    With dest removed from Q, b = (I - Q)^-1 q is the probability of reaching
    dest from each state (q being the one-step probabilities of moving into
    dest), and the conditional passage times are t = (I - Q)^-1 b / b.

    States that are very unlikely to reach dest have tiny b, where the LU
    solves lose relative accuracy. Both solves are iteratively refined,
    which is enough unless some b is below 1e-15 or the refinement does not
    converge (very ill-conditioned chains); then the subtraction-free
    elimination of Grassmann, Taksar & Heyman (1985), which is accurate
    element-wise regardless of conditioning, is used instead (also when
    the LU factorisation fails as numerically singular).

    Args:
        Q (scipy.sparse matrix): Transient transition matrix.
        dest (int): Index of the destination state.
        absorption (array-like, optional): Total absorption of each state
            (1 - row sums of Q). Only needed for the fallback; if None, it is
            computed from Q (less accurately).

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
    try:
        lu = splu(A)
    except RuntimeError:
        # Numerically singular (e.g. hitting times beyond floating point
        # range): the subtraction-free elimination below
        lu = None
    if lu is not None:
        b, ok_b = _refine(lu, A, qj, lu.solve(qj))
    else:
        ok_b = False
    if ok_b and np.all(b > 1e-15):
        x, ok_x = _refine(lu, A, b, lu.solve(b))
        if ok_x:
            out = np.zeros(n)
            out[keep] = x / b
            return out

    if absorption is None:
        absorption = np.clip(1.0 - np.asarray(Q.sum(axis=1)).ravel(), 0, 1)
    return _cond_passage_gth(Q, dest, absorption)


def _refine(lu, A, rhs, x, max_steps=10, rtol=1e-13):
    """
    Iterative refinement of a solution x of A x = rhs, reusing the LU
    factorisation, until the correction is negligible relative to each
    element of x.

    Returns:
        tuple: The refined x, and whether the refinement converged.
    """
    for _ in range(max_steps):
        dx = lu.solve(rhs - A @ x)
        x = x + dx
        if np.all(np.abs(dx) <= rtol * np.abs(x)):
            return x, True
    return x, False


def _cond_passage_gth(Q, dest, absorption):
    """
    Conditional passage times to dest by subtraction-free (GTH) sparse
    elimination of I - Q (dest removed), in minimum-degree order.

    Each row's slack (its row sum: absorption plus the probability of
    moving into dest) is carried through the elimination, and pivots are
    computed as slack plus the remaining off-diagonal magnitudes, so that
    only non-negative numbers are ever added (Grassmann, Taksar & Heyman
    1985, Operations Research 33:1107-1116).
    """
    Q = sp.csr_matrix(Q)
    n = Q.shape[0]
    slack = np.asarray(absorption, dtype=float).ravel().copy()
    off = [dict() for _ in range(n)]     # off[i][j] = -A[i, j] >= 0
    q = np.zeros(n)
    for i in range(n):
        if i == dest:
            continue
        lo, hi = Q.indptr[i], Q.indptr[i + 1]
        for j, v in zip(Q.indices[lo:hi], Q.data[lo:hi]):
            if j == dest:
                q[i] = v
                slack[i] += v
            elif j != i and v > 0:
                off[i][j] = v
    into = [set() for _ in range(n)]     # rows with an entry in column j
    for i in range(n):
        for j in off[i]:
            into[j].add(i)

    alive = set(range(n)) - {dest}
    degree = {i: len(off[i]) + len(into[i]) for i in alive}
    heap = [(d, i) for i, d in degree.items()]
    heapq.heapify(heap)
    order, pivots, lower, upper = [], {}, {}, {}
    while heap:
        d, k = heapq.heappop(heap)
        if k not in alive or d != degree[k]:
            continue
        alive.remove(k)
        pivot = slack[k] + sum(off[k].values())
        mult = []
        for i in into[k]:
            if i not in alive:
                continue
            f = off[i].pop(k) / pivot
            mult.append((i, f))
            slack[i] += f * slack[k]
            for j, v in off[k].items():
                if j != i:
                    off[i][j] = off[i].get(j, 0.0) + f * v
                    into[j].add(i)
        for j in off[k]:
            into[j].discard(k)
        for i in set(off[k]) | {i for i, _ in mult}:
            if i in alive:
                degree[i] = len(off[i]) + len(into[i])
                heapq.heappush(heap, (degree[i], i))
        order.append(k)
        pivots[k], lower[k], upper[k] = pivot, mult, off[k]

    def solve(rhs):
        rhs = rhs.copy()
        for k in order:
            for i, f in lower[k]:
                rhs[i] += f * rhs[k]
        x = np.zeros(n)
        for k in reversed(order):
            x[k] = (rhs[k] + sum(v * x[j] for j, v in upper[k].items())) \
                / pivots[k]
        return x

    b = solve(q)
    x = solve(b)
    out = np.zeros(n)
    keep = np.arange(n) != dest
    out[keep] = x[keep] / b[keep]
    return out


def CFPT(Q, edge_site_indices, absorption=None):
    """
    Pairwise conditional first passage times among a set of states.

    Args:
        Q (scipy.sparse matrix): Transient transition matrix.
        edge_site_indices (list): State index of each site.
        absorption (array-like, optional): Total absorption of each state
            (see cond_passage).

    Returns:
        numpy.ndarray: Matrix T where T[i, j] is the conditional passage
        time from site i to site j (0 on the diagonal).
    """
    sites = np.asarray(edge_site_indices)
    N = len(sites)
    cfpt_matrix = np.zeros((N, N))
    for j, dest in enumerate(sites):
        cfpt_matrix[:, j] = cond_passage(Q, dest, absorption)[sites]
    np.fill_diagonal(cfpt_matrix, 0.0)
    return cfpt_matrix
