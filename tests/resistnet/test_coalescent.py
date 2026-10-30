import os

import numpy as np
import pandas as pd
import pytest

import resistnet.coalescent as co

REFERENCE = os.path.join(os.path.dirname(__file__), "reference")

# A small river (a tree of six demes: 0-1-2-3 main stem, 4-5 a tributary
# joining at 1), with a jump rate in each direction of every edge
# (asymmetric) and a coalescence rate per deme
TREE = [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5)]
RATES = {(0, 1): 0.7, (1, 0): 1.3, (1, 2): 0.4, (2, 1): 2.1, (2, 3): 1.0,
         (3, 2): 0.6, (1, 4): 0.9, (4, 1): 0.3, (4, 5): 1.7, (5, 4): 0.5}
GAMMA = np.array([1.0, 0.8, 1.5, 0.6, 1.2, 0.9])


def _edges():
    return list(RATES), np.array(list(RATES.values()))


def _divergences():
    """Divergences to fit: the model's coalescence times, perturbed by a
    fixed pattern (so that the residuals are not zero)."""
    edges, rates = _edges()
    C = co.coalescence_times(co.generator(rates, edges, 6), GAMMA)
    k = np.add.outer(np.arange(6), np.arange(6))
    return C * (1 + 0.03 * np.sin(k + 1))


def test_coalescence_times_equation():
    # C solves G C + C G' - diag(gamma) o C = -1 1' (their equation 7)
    edges, rates = _edges()
    G = co.generator(rates, edges, 6).toarray()
    np.testing.assert_allclose(G.sum(axis=1), 0, atol=1e-15)
    C = co.coalescence_times(G, GAMMA)
    lhs = G @ C + C @ G.T - np.diag(GAMMA * np.diag(C))
    np.testing.assert_allclose(lhs, -1, rtol=1e-12)
    np.testing.assert_allclose(C, C.T)
    # One deme: 1 / gamma
    assert co.coalescence_times(np.zeros((1, 1)), [2.0])[0, 0] == \
        pytest.approx(0.5)


def test_upper_order():
    # As R's M[upper.tri(M, diag = TRUE)]: column by column
    M = np.arange(9).reshape(3, 3)
    assert list(co.upper(M)) == [0, 1, 4, 2, 5, 8]


@pytest.mark.skipif(
    not os.path.exists(os.path.join(REFERENCE, "lundgren_ralph_H.tsv")),
    reason="reference/lundgren_ralph.R not run yet (needs their code)")
def test_matches_lundgren_ralph_code():
    # Their own HCalc and H_and_lllh (gene.flow.inference, run in R; see
    # reference/lundgren_ralph.R) on the same graph, rates, coalescence
    # rates and divergences
    edges, rates = _edges()
    H_ref = pd.read_csv(os.path.join(REFERENCE, "lundgren_ralph_H.tsv"),
                        sep="\t", header=None).to_numpy()
    C = co.coalescence_times(co.generator(rates, edges, 6), GAMMA)
    np.testing.assert_allclose(C, H_ref, rtol=1e-12)
    D = pd.read_csv(os.path.join(REFERENCE, "lundgren_ralph_inputs.tsv"),
                    sep="\t", header=None).to_numpy()
    np.testing.assert_allclose(D, _divergences(), rtol=1e-12)
    lllh = float(open(os.path.join(REFERENCE,
                                   "lundgren_ralph_lllh.txt")).read())
    lp = co.log_posterior(co.upper(D), co.upper(C), (D.mean() / 50) ** 2,
                          rates, GAMMA)
    assert lp == pytest.approx(lllh, rel=1e-12)
    model = co.Model(edges, 6, const_coal=False)
    np.testing.assert_allclose(model.hcalc(rates, GAMMA), co.upper(C))


def test_divergence():
    # Between groups: the share of loci at which two haplotypes differ,
    # averaged over all pairs; within: over pairs of distinct haplotypes
    rng = np.random.default_rng(1)
    X = (rng.random((50, 7)) < 0.4).astype(int)
    groups = np.array([0, 0, 0, 1, 1, 2, 2])
    D = co.divergence(X, groups)
    for a in range(3):
        for b in range(3):
            ia = np.flatnonzero(groups == a)
            ib = np.flatnonzero(groups == b)
            pairs = [(i, j) for i in ia for j in ib if i != j]
            ref = np.mean([np.mean(X[:, i] != X[:, j]) for i, j in pairs])
            assert D[a, b] == pytest.approx(ref)


def test_mh_recovers_rates():
    # Noise-free divergences of a two-deme graph with asymmetric rates:
    # the posterior concentrates near the truth (the start is the best of
    # prior draws; adaptive proposals)
    edges = [(0, 1), (1, 0)]
    rates = np.array([0.5, 2.0])
    C = co.coalescence_times(co.generator(rates, edges, 2), [1.0, 1.0])
    model = co.Model(edges, 2, const_coal=True)
    rng = np.random.default_rng(3)
    out = co.mh(model, C, iterations=6000, sig2=(C.mean() / 500) ** 2,
                rng=rng)
    g = out["g"][3000:].mean(axis=0)
    gam = out["gam"][3000:].mean()
    pred = co.coalescence_times(co.generator(g, edges, 2), [gam, gam])
    np.testing.assert_allclose(pred, C, rtol=0.02)


if __name__ == "__main__":
    # Writes the inputs read by reference/lundgren_ralph.R
    edges, rates = _edges()
    pd.DataFrame({"x": [e[0] for e in edges], "y": [e[1] for e in edges],
                  "rate": rates}).to_csv(
        os.path.join(REFERENCE, "lundgren_ralph_edges.tsv"), sep="\t",
        index=False)
    np.savetxt(os.path.join(REFERENCE, "lundgren_ralph_gamma.txt"), GAMMA,
               fmt="%.17g")
    np.savetxt(os.path.join(REFERENCE, "lundgren_ralph_inputs.tsv"),
               _divergences(), delimiter="\t", fmt="%.17g")
