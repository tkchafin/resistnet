"""
The coalescence-time model of Lundgren & Ralph (2019), "Are populations
like a circuit? Comparing isolation by resistance to a new coalescent-based
method", Molecular Ecology Resources 19:1388-1406 (code:
github.com/elundgre/gene-flow-inference, R package gene.flow.inference).

Lineages move back in time among demes as a continuous-time Markov chain
with generator G (G[x, y] = jump rate from x to y on each edge of a graph,
which can differ between the two directions), and two lineages in the same
deme x coalesce at rate gamma_x. The mean coalescence times C of two
lineages sampled in demes x and y satisfy (their equation 7)

    G C + C G' - diag(gamma) o C = -1 1',

solved here as they do, through the product chain (G (x) I + I (x) G).
Mean pairwise divergence D (including within-deme divergence on the
diagonal) is modelled as C plus independent Gaussian noise of variance
sig2 (their equation 8), with independent exponential priors on the jump
rates and the coalescence rates, and sampled by Metropolis-Hastings as
their findG.MH and run.mcmc.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def generator(rates, edges, n):
    """
    The generator matrix of lineage movement.

    Args:
        rates (array-like): Jump rate on each directed edge.
        edges (list): Directed edges (x, y), in the order of rates.
        n (int): Number of demes.

    Returns:
        scipy.sparse.csr_matrix: G, with diagonal -(row sums).
    """
    rates = np.asarray(rates, dtype=float)
    x = np.array([e[0] for e in edges], dtype=int)
    y = np.array([e[1] for e in edges], dtype=int)
    G = sp.csr_matrix((rates, (x, y)), shape=(n, n))
    return (G - sp.diags(np.asarray(G.sum(axis=1)).ravel())).tocsr()


def coalescence_times(G, gamma):
    """
    Mean coalescence times of pairs of lineages (their equation 7).

    Args:
        G (matrix): Generator of lineage movement (see generator()).
        gamma (float or array-like): Coalescence rate in each deme.

    Returns:
        numpy.ndarray: C[x, y], the mean coalescence time of lineages
                       sampled in demes x and y (symmetric).
    """
    G = sp.csr_matrix(G)
    n = G.shape[0]
    eye = sp.identity(n, format="csr")
    # Row-major vec(C): (G (x) I) vec(C) = vec(G C) and
    # (I (x) G) vec(C) = vec(C G')
    GG = sp.kron(G, eye) + sp.kron(eye, G)
    coal = np.zeros(n * n)
    coal[np.arange(n) * (n + 1)] = np.broadcast_to(
        np.asarray(gamma, dtype=float), (n,))
    A = (GG - sp.diags(coal)).tocsc()
    C = spl.spsolve(A, -np.ones(n * n)).reshape(n, n)
    return (C + C.T) / 2


def upper(M):
    """The upper triangle of M, diagonal included, as a vector (column by
    column, as R's M[upper.tri(M, diag = TRUE)])."""
    M = np.asarray(M)
    # (row, column) of the lower triangle, row by row, read transposed:
    # the upper triangle column by column
    r, c = np.tril_indices(M.shape[0])
    return M[c, r]


def log_posterior(h, hcalc, sig2, rates, gamma, lam_g=1.0, lam_gam=1.0):
    """
    Log posterior density, up to a constant, as their loglikelihood():
    exponential priors on the jump rates and coalescence rates, and
    Gaussian errors on the divergences (missing ones ignored).

    Args:
        h (array-like): Observed divergences (upper triangle, see upper()).
        hcalc (array-like): Coalescence times (the same entries).
        sig2 (float or array-like): Error variance.
        rates, gamma (array-like): Parameters.
        lam_g, lam_gam (float): Prior rates (1 / prior means).

    Returns:
        float: The log posterior (-inf if a parameter is not positive).
    """
    rates = np.asarray(rates, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    if not (np.all(rates > 0) and np.all(gamma > 0)):
        return -np.inf
    r = (np.asarray(hcalc) - np.asarray(h)) ** 2 / (2 * np.asarray(sig2))
    return float(-np.sum(lam_g * rates) - np.sum(lam_gam * gamma)
                 - np.nansum(r))


class Model:
    """
    The coalescence-time model on a graph: which directed edges have jump
    rates, and whether the coalescence rate is the same in every deme.

    Args:
        edges (list): Directed edges (x, y) with a free jump rate.
        n (int): Number of demes.
        const_coal (bool): One coalescence rate for all demes (True) or one
                           per deme.
    """

    def __init__(self, edges, n, const_coal=True):
        self.edges = list(edges)
        self.n = n
        self.const_coal = const_coal
        self.ng = len(self.edges)
        self.ngam = 1 if const_coal else n

    def hcalc(self, rates, gamma):
        """Coalescence times as a vector (upper triangle, see upper())."""
        gam = np.repeat(gamma, self.n) if self.const_coal else gamma
        return upper(coalescence_times(
            generator(rates, self.edges, self.n), gam))


def mh(model, H, iterations=10000, sig2=None, start=None, sd_init=1 / 200,
       fixed_sd=False, n_start=1000, lam_g=1.0, lam_gam=1.0, rng=None):
    """
    Metropolis-Hastings sampling of the jump rates and coalescence rates,
    as their findG.MH: proposals add Gaussian noise and reflect at zero
    (|x + e|); the jump rates and then the coalescence rates are updated as
    two blocks each iteration; unless fixed_sd, each block's proposal
    standard deviation is multiplied by 0.8 when more than 90 of the last
    100 proposals were rejected, and by 1.2 when fewer than 60 were.
    Without a start, the best of n_start draws from the priors is used.

    Args:
        model (Model): The graph.
        H (array-like): Observed divergences (n x n; missing values NaN).
        iterations (int): Number of iterations.
        sig2 (float): Error variance (default (mean(H) / 50) ** 2).
        start (tuple, optional): Starting (rates, gamma).
        sd_init (float or tuple): Initial proposal standard deviations
                                  (rates, gamma).
        fixed_sd (bool): Keep the proposal standard deviations fixed.
        n_start (int): Prior draws for the starting point.
        lam_g, lam_gam (float): Prior rates.
        rng (numpy.random.Generator, optional): Random numbers.

    Returns:
        dict: Samples ("g", "gam"), log posteriors ("lp"), rejections,
              proposal standard deviations ("sdG", "sdgam").
    """
    rng = np.random.default_rng() if rng is None else rng
    H = np.asarray(H, dtype=float)
    H = (H + H.T) / 2
    h = upper(H)
    if sig2 is None:
        sig2 = (np.nanmean(H) / 50) ** 2
    ng, ngam = model.ng, model.ngam

    def lp(g, gam):
        if not (np.all(g > 0) and np.all(gam > 0)):
            return -np.inf
        return log_posterior(h, model.hcalc(g, gam), sig2, g, gam, lam_g,
                             lam_gam)

    if start is None:
        gs = rng.exponential(1 / lam_g, (n_start, ng))
        gams = rng.exponential(1.0, (n_start, ngam))
        lps = [lp(gs[k], gams[k]) for k in range(n_start)]
        best = int(np.argmax(lps))
        g, gam, cur = gs[best], gams[best], lps[best]
    else:
        g = np.asarray(start[0], dtype=float).copy()
        gam = np.asarray(start[1], dtype=float).copy()
        cur = lp(g, gam)
    sdG, sdgam = (sd_init if isinstance(sd_init, tuple)
                  else (sd_init, sd_init))
    G_s = np.empty((iterations, ng))
    gam_s = np.empty((iterations, ngam))
    lps = np.empty(iterations)
    rej_g = np.zeros(iterations, dtype=int)
    rej_gam = np.zeros(iterations, dtype=int)
    sdGs = np.empty(iterations)
    sdgams = np.empty(iterations)
    G_s[0], gam_s[0], lps[0] = g, gam, cur
    sdGs[0], sdgams[0] = sdG, sdgam
    for it in range(1, iterations):
        gp = np.abs(g + rng.normal(0, sdG, ng))
        new = lp(gp, gam)
        if new > cur or rng.random() < np.exp(new - cur):
            g, cur = gp, new
        else:
            rej_g[it] = 1
        gamp = np.abs(gam + rng.normal(0, sdgam, ngam))
        new = lp(g, gamp)
        if new > cur or rng.random() < np.exp(new - cur):
            gam, cur = gamp, new
        else:
            rej_gam[it] = 1
        G_s[it], gam_s[it], lps[it] = g, gam, cur
        sdGs[it], sdgams[it] = sdG, sdgam
        if (it + 1) % 100 == 0 and it != iterations - 1 and not fixed_sd:
            window = slice(it - 99, it + 1)
            for name in ("g", "gam"):
                rej = (rej_g if name == "g" else rej_gam)[window].sum()
                factor = 0.8 if rej > 90 else 1.2 if rej < 60 else 1.0
                if name == "g":
                    sdG *= factor
                else:
                    sdgam *= factor
    return {"g": G_s, "gam": gam_s, "lp": lps, "reject_g": rej_g,
            "reject_gam": rej_gam, "sdG": sdGs, "sdgam": sdgams}


def run_mcmc(model, H, h_se, preburn=100000, burn=100000, iterations=100000,
             rng=None):
    """
    The three runs of their run.mcmc for observed divergences: divergences
    (and their standard errors) are scaled so that their mean equals the
    number of demes; a "pre-burn-in" from the best of 1000 prior draws with
    a looser likelihood (10 x the error variance), a burn-in from its last
    state with adaptive proposals, and the sampling run from the burn-in's
    last state with proposals fixed at the mean of the burn-in's second
    half.

    Args:
        model (Model): The graph.
        H (array-like): Observed divergences, n x n with the within-deme
                        divergences on the diagonal.
        h_se (float or array-like): Standard errors of the divergences
                                    (upper-triangle vector, or a scalar).
        preburn, burn, iterations (int): Lengths of the three runs.
        rng (numpy.random.Generator, optional): Random numbers.

    Returns:
        dict: The sampling run (see mh()), plus "scale" (the scaling of the
              divergences: rates are in units of the scaled divergences),
              posterior means over its second half ("g_mean", "gam_mean"),
              medians and the pre-burn-in and burn-in runs.
    """
    rng = np.random.default_rng() if rng is None else rng
    H = np.asarray(H, dtype=float)
    scale = model.n / np.nanmean(H)
    Hs = H * scale
    sig2 = (np.asarray(h_se, dtype=float) * scale) ** 2
    pre = mh(model, Hs, preburn, sig2=sig2 * 10, rng=rng)
    bur = mh(model, Hs, burn, sig2=sig2,
             start=(pre["g"][-1], pre["gam"][-1]), rng=rng)
    half = slice(int(round(burn / 2)) - 1, burn)
    sd = (float(np.mean(bur["sdG"][half])),
          float(np.mean(bur["sdgam"][half])))
    out = mh(model, Hs, iterations, sig2=sig2,
             start=(bur["g"][-1], bur["gam"][-1]), sd_init=sd,
             fixed_sd=True, rng=rng)
    last = slice(int(round(iterations / 2)) - 1, iterations)
    out.update(scale=scale, preburn=pre, burn=bur,
               g_mean=out["g"][last].mean(axis=0),
               gam_mean=out["gam"][last].mean(axis=0),
               g_median=np.median(out["g"], axis=0),
               gam_median=np.median(out["gam"], axis=0))
    return out


def divergence(haplotypes, groups):
    """
    Mean pairwise divergence between groups, and within each group (pairs
    of distinct haplotypes, sampled without replacement), per locus, as
    their genetic distance D: the mean number of differing loci between a
    haplotype from group x and one from group y, divided by the number of
    loci.

    Args:
        haplotypes (array-like): 0/1 alleles, loci x haplotypes.
        groups (array-like): Group of each haplotype (0 ... n - 1).

    Returns:
        numpy.ndarray: D (n x n).
    """
    X = np.asarray(haplotypes, dtype=float)
    groups = np.asarray(groups)
    n = int(groups.max()) + 1
    cnt = np.array([np.sum(groups == k) for k in range(n)])
    p = np.column_stack([X[:, groups == k].mean(axis=1) for k in range(n)])
    D = p.T @ (1 - p) + (1 - p).T @ p
    D /= X.shape[0]
    within = 2 * (p * (1 - p)).mean(axis=0) * cnt / (cnt - 1)
    D[np.diag_indices(n)] = within
    return D
