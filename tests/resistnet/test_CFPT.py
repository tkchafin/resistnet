import numpy as np
import pytest
import scipy.sparse as sp

from resistnet.CFPT import transition_matrix, cond_passage, CFPT


@pytest.fixture
def small_chain():
    """
    Small branching network with asymmetric movement weights, absorption and
    fidelity.
    """
    edges = [(0, 1), (1, 2), (2, 3), (2, 4), (4, 5), (1, 6), (6, 7)]
    w_fwd = [1.0, 2.0, 0.5, 3.0, 1.5, 0.8, 2.5]
    w_rev = [0.4, 1.0, 2.0, 0.6, 1.2, 3.0, 0.7]
    rows = [a for a, b in edges] + [b for a, b in edges]
    cols = [b for a, b in edges] + [a for a, b in edges]
    W = sp.csr_matrix((w_fwd + w_rev, (rows, cols)), shape=(8, 8))
    absorption = np.array([0.05, 0, 0, 0.1, 0, 0.02, 0, 0.2])
    fidelity = np.array([0, 0.1, 0, 0, 0.2, 0, 0, 0])
    return W, absorption, fidelity


# samc::cond_passage(samc(data = P), dest = d) for d in 1:8 (samc 3.2.1),
# with P built from small_chain; column d holds the times to destination d
SAMC_COND_PASSAGE = np.array([
    [0, 1, 3.5324366248061, 39.5843102087513, 5.88492362132889,
     8.50378591301045, 33.3579940429882, 54.9470018696965],
    [45.4682822341621, 0, 2.5324366248061, 38.5843102087513,
     4.88492362132889, 7.50378591301045, 32.3579940429882, 53.9470018696965],
    [65.1868451084136, 19.7185628742515, 0, 36.0518735839452,
     2.35248699652279, 4.97134928820434, 52.0765569172396, 73.6655647439479],
    [66.1868451084136, 20.7185628742515, 1, 0, 3.35248699652279,
     5.97134928820434, 53.0765569172396, 74.6655647439479],
    [71.6868451084136, 26.2185628742515, 6.5, 42.5518735839453, 0,
     2.61886229168155, 58.5765569172396, 80.1655647439478],
    [72.6868451084136, 27.2185628742515, 7.5, 43.5518735839453, 1, 0,
     59.5765569172396, 81.1655647439478],
    [47.6111393770192, 2.14285714285714, 4.67529376766325, 40.7271673516085,
     7.02778076418604, 9.64664305586759, 0, 21.5890078267084],
    [48.6111393770192, 3.14285714285714, 5.67529376766325, 41.7271673516085,
     8.02778076418604, 10.6466430558676, 1, 0],
])


def test_transition_matrix_rows(small_chain):
    W, absorption, fidelity = small_chain
    Q = transition_matrix(W, absorption, fidelity).toarray()
    # Rows of [Q, absorption] sum to 1, fidelity is the diagonal
    assert np.allclose(Q.sum(axis=1) + absorption, 1)
    assert np.allclose(np.diag(Q), fidelity)
    # Moves are proportional to the weights
    Wd = W.toarray()
    move = 1 - absorption - fidelity
    off = Q - np.diag(np.diag(Q))
    assert np.allclose(off, move[:, None] * Wd / Wd.sum(axis=1)[:, None])


def test_transition_matrix_rejects_invalid(small_chain):
    W, absorption, fidelity = small_chain
    with pytest.raises(ValueError):
        transition_matrix(-W, absorption)
    with pytest.raises(ValueError):
        transition_matrix(W, absorption + 0.95)
    isolated = sp.csr_matrix(W.toarray()[:7, :7])
    isolated = sp.block_diag([isolated, sp.csr_matrix((1, 1))])
    with pytest.raises(ValueError):
        transition_matrix(isolated, np.zeros(8))


def test_cond_passage_matches_samc(small_chain):
    Q = transition_matrix(*small_chain)
    T = np.column_stack([cond_passage(Q, d) for d in range(8)])
    np.testing.assert_allclose(T, SAMC_COND_PASSAGE, rtol=1e-12)


def test_CFPT_orientation(small_chain):
    # Rows are origins, columns destinations
    Q = transition_matrix(*small_chain)
    sites = [5, 0, 7]
    T = CFPT(Q, sites)
    np.testing.assert_allclose(
        T, SAMC_COND_PASSAGE[np.ix_(sites, sites)], rtol=1e-12
    )
    assert np.all(np.diag(T) == 0)


def _decimal_cond_passage(Q, dest, absorption):
    """
    60-digit reference: dense Gaussian elimination in decimal arithmetic,
    with each diagonal of I - Q from the stated absorption (so that the
    rows of [Q, absorption] sum to exactly 1).
    """
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    Q = np.asarray(Q)
    n = Q.shape[0]
    idx = [i for i in range(n) if i != dest]
    A = [[-Decimal(float(Q[i, j])) for j in idx] for i in idx]
    for r, i in enumerate(idx):
        A[r][r] = (Decimal(float(absorption[i])) + Decimal(float(Q[i, dest]))
                   + sum(Decimal(float(Q[i, k])) for k in idx if k != i))
    q = [Decimal(float(Q[i, dest])) for i in idx]

    def solve(A, rhs):
        m = len(rhs)
        A = [row[:] for row in A]
        rhs = rhs[:]
        for k in range(m):
            for i in range(k + 1, m):
                f = A[i][k] / A[k][k]
                if f:
                    for j in range(k, m):
                        A[i][j] -= f * A[k][j]
                    rhs[i] -= f * rhs[k]
        x = [Decimal(0)] * m
        for k in reversed(range(m)):
            x[k] = (rhs[k] - sum(A[k][j] * x[j] for j in range(k + 1, m))) \
                / A[k][k]
        return x

    b = solve(A, q)
    x = solve(A, b)
    out = np.zeros(n)
    out[idx] = [float(xi / bi) for xi, bi in zip(x, b)]
    return out


def _count_gth(monkeypatch):
    import resistnet.CFPT as C
    calls = []
    orig = C._cond_passage_gth

    def counted(*args):
        calls.append(1)
        return orig(*args)
    monkeypatch.setattr(C, "_cond_passage_gth", counted)
    return calls


def test_gth_matches_samc(small_chain):
    from resistnet.CFPT import _cond_passage_gth
    W, absorption, fidelity = small_chain
    Q = transition_matrix(W, absorption, fidelity)
    T = np.column_stack(
        [_cond_passage_gth(Q, d, absorption) for d in range(8)])
    np.testing.assert_allclose(T, SAMC_COND_PASSAGE, rtol=1e-12)


# A branching chain where moving towards deeper states is resistant (as for
# a directional variable that increases along the way, with a maximum
# resistance of 500); from random draws, fixed here
TREE_PARENT = [
    0, 1, 2, 3, 4, 3, 5, 5, 6, 8, 9, 11, 10, 13, 14, 14, 15, 17, 17, 17, 20,
    19, 22, 21, 24, 23, 25, 26, 28, 29, 28, 31, 32, 31, 33, 34, 34, 35, 37,
    38, 40, 41, 40, 41, 44, 43, 45, 46, 48, 48, 50, 51, 52, 53, 53, 55, 55,
    56, 56]
TREE_NOISE = [
    -0.241537, 0.023332, -0.243894, -0.1267, -0.683835, 0.221653, -0.216782,
    -0.021809, 0.13372, 0.101727, -0.39427, 0.45545, 0.227194, 0.150041,
    -0.197879, -0.481283, -0.489047, -0.062526, 0.048201, -0.211008,
    0.548834, 0.090342, 0.299556, 0.271971, -0.321141, 0.298164, 0.135061,
    -0.380589, 0.217032, 0.321593, 0.292978, 0.315556, -0.423035, -0.074249,
    -0.163715, 0.016287, 0.031875, 0.337823, 0.225645, -0.248685, 0.270439,
    0.487838, 0.091907, -0.451685, -0.193411, -0.053852, -0.487006,
    -0.041489, -0.381799, 0.237785, 0.205313, -0.022941, -0.376436,
    0.086763, -0.191241, -0.206324, 0.33301, -0.227145, -0.413043, 0.237111]
TREE_ABSORBING = [19, 35, 36, 39, 48, 55]

# Movement weights spanning 6 orders of magnitude, which create pairs of
# states that trade the walker back and forth; from random draws, fixed here
TRAP_PARENT = [0, 0, 2, 2, 1, 0, 4, 2, 7, 3, 7, 2, 7, 1, 11, 1, 12, 13, 8, 15,
               12, 19, 17, 21, 1, 0, 17, 20, 21]
TRAP_LOG10W = [
    -0.607447, 0.043095, -1.627676, 0.901165, 2.827706, -1.207721, -0.223015,
    2.349683, 0.308405, -0.471212, 1.009744, -2.815079, -2.077221, 2.22794,
    -2.050561, -2.823199, 2.558373, 1.51671, 2.864584, 1.558314, 2.53685,
    -0.933519, 1.527796, -2.01355, 1.943341, 1.393872, -1.852008, -1.896197,
    -1.987777, 1.91324, 0.985952, -1.315094, 0.926536, 1.11613, 1.133295,
    2.975801, 2.069782, -0.420087, -2.717593, -1.628428, -2.55116, 1.706704,
    1.759992, -1.714397, 2.690827, -1.061453, -0.204504, 0.4857, 2.987845,
    2.381266, -2.134106, 2.360935, -1.780855, 1.817996, 1.02883, -2.090444,
    -2.678622, 2.576578]
TRAP_ABSORPTION = {4: 0.254144, 0: 0.307516, 25: 0.258239, 9: 0.079738,
                   29: 0.153705}


def directional_tree():
    n = len(TREE_PARENT) + 1
    depth = np.zeros(n)
    for i, p in zip(range(1, n), TREE_PARENT):
        depth[i] = depth[p] + 1
    v = depth + np.array(TREE_NOISE)
    rows = list(range(1, n)) + TREE_PARENT
    cols = TREE_PARENT + list(range(1, n))
    d = v[cols] - v[rows]
    c = (d - d.min()) / (d.max() - d.min())
    W = sp.csr_matrix((1 / (1 + 499 * c), (rows, cols)), shape=(n, n))
    absorption = np.zeros(n)
    absorption[TREE_ABSORBING] = 0.005
    return transition_matrix(W, absorption), absorption, int(np.argmax(depth))


def trap_chain():
    n = len(TRAP_PARENT) + 1
    rows = list(range(1, n)) + TRAP_PARENT
    cols = TRAP_PARENT + list(range(1, n))
    W = sp.csr_matrix((10 ** np.array(TRAP_LOG10W), (rows, cols)),
                      shape=(n, n))
    absorption = np.zeros(n)
    for i, a in TRAP_ABSORPTION.items():
        absorption[i] = a
    return transition_matrix(W, absorption), absorption, 0


@pytest.mark.parametrize("case", [directional_tree, trap_chain])
def test_ill_conditioned_chains(case, monkeypatch):
    # The LU solves are inaccurate here even after refinement (by ~4e-7 and
    # ~1e-11), which does not converge, so the GTH elimination is used
    Q, absorption, dest = case()
    calls = _count_gth(monkeypatch)
    t = cond_passage(Q, dest, absorption)
    assert calls
    ref = _decimal_cond_passage(Q.toarray(), dest, absorption)
    keep = np.arange(Q.shape[0]) != dest
    np.testing.assert_allclose(t[keep], ref[keep], rtol=1e-12)
