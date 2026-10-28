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
