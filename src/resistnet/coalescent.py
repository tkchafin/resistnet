"""
Expected coalescence times and linearised FST among sites on a river
network, from Lundgren & Ralph's (2019) coalescence-time model ("Are
populations like a circuit? Comparing isolation by resistance to a new
coalescent-based method", Molecular Ecology Resources 19:1388-1406; code:
github.com/elundgre/gene-flow-inference, R package gene.flow.inference).

Lineages move back in time among demes as a continuous-time Markov chain
with generator G (G[x, y] = jump rate from x to y on each edge of a graph,
which can differ between the two directions), and two lineages in the same
deme x coalesce at rate gamma_x. The mean coalescence times C of two
lineages sampled in demes x and y satisfy (their equation 7)

    G C + C G' - diag(gamma) o C = -1 1',

solved here as they do, through the product chain (G (x) I + I (x) G).

On a river (RiverCoalescent), the demes are the points of a lattice of the
network (SegmentLattice: points at confluences, tips, the outlet and the
sites, and about cell_km apart in between). The rate of a lineage at point
x jumping to its neighbour y is the link's conductance, 1 / (its
resistance: resistance per km in series over the river it covers), times
exp(-rate x km moved upstream) with a flow term, divided by the length of
river around x (half its links): a finite volume version of a diffusion,
so that results do not depend on the lattice's spacing once it is fine.
Two lineages at the same point coalesce at rate 1 / (density x length of
river around it) (the population living there, density per km). The
genetic distance between two sites is linearised FST from the coalescence
times (Slatkin 1991), (C_ij - Cw) / Cw with Cw the mean of the two
within-site times. With a flow term, lineages traced back in time move
upstream more easily when the rate is negative, i.e. gene flow is biased
downstream.

The resistance per km of each reach is an input (1: distance alone): the
coalescent model (coalescent_network) gives it each candidate surface.
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


class SegmentLattice:
    """
    A lattice of points on a river network directed towards its outlet:
    points where reaches meet (confluences), at the network's tips, outlet
    and the given nodes (the sites), and about cell_km apart along each
    stretch of river between them, whatever the reaches' boundaries. A
    link's resistance is the sum over the reaches it overlaps of their
    resistance per km times the overlap (resistances in series), and its
    conductance the inverse; with resistance 1 per km, 1 / link length.

    Attributes:
        n (int): Number of points.
        pos (numpy.ndarray): Distance of each point to the outlet along the
                             river (km).
        a, b (numpy.ndarray): Links (both directions: a -> b).
        length (numpy.ndarray): Length of each link (km).
        upstream (numpy.ndarray): km moved upstream along each link
                                  (negative downstream).
        overlap (scipy.sparse.csr_matrix): km of each reach (columns, edge
                                           order) covered by each link
                                           (rows, as a and b).
        reach_link (numpy.ndarray): Length of the link covering each
                                    reach's middle.
        cell_km (float): Target distance between neighbouring points.
    """

    def __init__(self, K, edge_order, reachid_col, length_col, cell_km=5.0,
                 keep_nodes=()):
        """
        Args:
            K (networkx.DiGraph): The network, each edge (reach) directed
                                  downstream, with one outlet.
            edge_order (list): The reaches' IDs (their reachid_col
                               values), in the order of the resistances
                               given to weights().
            reachid_col (str): Edge attribute with the reach ID.
            length_col (str): Edge attribute with the reach length (km).
            cell_km (float): Target distance between neighbouring points.
            keep_nodes (iterable): Nodes that must be points (the sites).
        """
        if cell_km <= 0:
            raise ValueError("cell_km must be positive")
        self.cell_km = float(cell_km)
        reach_ix = {e: i for i, e in enumerate(edge_order)}
        outlets = [v for v in K.nodes if K.out_degree(v) == 0]
        if len(outlets) != 1:
            raise ValueError(f"The network has {len(outlets)} outlets; "
                             "the lattice needs one")
        keep = set(keep_nodes)
        anchors = [v for v in K.nodes
                   if K.in_degree(v) != 1 or K.out_degree(v) != 1
                   or v in keep]
        self._ix = {v: i for i, v in enumerate(anchors)}
        # Distance of every node to the outlet
        down = {u: (v, float(d[length_col])) for u, v, d in K.edges(data=True)}
        dist = {outlets[0]: 0.0}
        for v in K.nodes:
            chain = [v]
            while chain[-1] not in dist:
                chain.append(down[chain[-1]][0])
            for s in reversed(chain[:-1]):
                dist[s] = dist[down[s][0]] + down[s][1]
        P = [dist[v] for v in anchors]
        a, b, length, rows, cols, vals = [], [], [], [], [], []
        self.reach_link = np.full(len(edge_order), np.nan)
        link = 0
        for u in anchors:
            if u not in down:
                continue
            # The stretch from u down to the next anchor
            reaches, cur = [], u
            while True:
                v, L = down[cur]
                reaches.append((reach_ix[K.edges[cur, v][reachid_col]], L))
                cur = v
                if cur in self._ix:
                    break
            ends = np.cumsum([L for _, L in reaches])
            starts = ends - np.array([L for _, L in reaches])
            total = ends[-1]
            pieces = max(1, int(round(total / self.cell_km)))
            cuts = np.linspace(0.0, total, pieces + 1)
            prev = self._ix[u]
            for k in range(1, pieces + 1):
                if k < pieces:
                    P.append(dist[u] - cuts[k])
                    nxt = len(P) - 1
                else:
                    nxt = self._ix[cur]
                lo, hi = cuts[k - 1], cuts[k]
                for (e, _), s0, s1 in zip(reaches, starts, ends):
                    over = min(hi, s1) - max(lo, s0)
                    if over > 0:
                        rows.append(link)
                        cols.append(e)
                        vals.append(over)
                    mid = (s0 + s1) / 2
                    if lo <= mid <= hi:
                        self.reach_link[e] = hi - lo
                a.append(prev)
                b.append(nxt)
                length.append(hi - lo)
                prev = nxt
                link += 1
        a, b = np.array(a, dtype=int), np.array(b, dtype=int)
        length = np.array(length, dtype=float)
        cover = sp.csr_matrix((vals, (rows, cols)),
                              shape=(link, len(edge_order)))
        self.n = len(P)
        self.pos = np.array(P, dtype=float)
        self.a = np.r_[a, b]
        self.b = np.r_[b, a]
        self.length = np.r_[length, length]
        self.overlap = sp.vstack([cover, cover]).tocsr()
        self.upstream = self.pos[self.b] - self.pos[self.a]

    def states(self, nodes):
        """Point index of each of the given network nodes (which must be
        confluences, tips, the outlet or keep_nodes)."""
        return [self._ix[v] for v in nodes]

    def weights(self, reach_r=None, rate=0.0):
        """
        Movement weights between neighbouring points: 1 / the link's
        resistance (resistances per km in series over the river it covers),
        times exp(-rate x km moved upstream) with a flow term.

        Args:
            reach_r (array-like, optional): Resistance per km of each reach
                                            (default: 1, distance alone).
            rate (float): Flow rate per km (0: none).

        Returns:
            scipy.sparse.csr_matrix: Weights (row = from, column = to).
        """
        if reach_r is None:
            resist = self.length
        else:
            resist = self.overlap @ np.asarray(reach_r, dtype=float)
        w = 1.0 / resist
        if rate != 0:
            w = w * np.exp(-rate * self.upstream)
        return sp.csr_matrix((w, (self.a, self.b)), shape=(self.n, self.n))

    def reach_moves(self, rate=0.0):
        """
        Probability of moving downstream and upstream from a point inside
        each reach's stretch, away from its ends: both neighbours are one
        link (of the length covering the reach's middle) away, so only the
        flow term makes them differ: 1 / (1 + exp(-2 x rate x link length))
        downstream.

        Args:
            rate (float): Flow rate per km.

        Returns:
            numpy.ndarray: Columns (downstream, upstream), in edge order.
        """
        down = 1.0 / (1.0 + np.exp(-2.0 * rate * self.reach_link))
        return np.column_stack((down, 1.0 - down))


class RiverCoalescent:
    """
    Coalescence times and linearised FST among the sites of a river
    network (see the module's description).

    Attributes:
        lattice (SegmentLattice): The lattice the lineages move on.
        sites (list): Lattice point of each site.
        density (float): Population per km of river.
        around (numpy.ndarray): Length of river around each point (half
                                its links).
    """

    def __init__(self, K, edge_order, sites, reachid_col, length_col,
                 cell_km=5.0, density=1000.0):
        """
        Args:
            K, edge_order, reachid_col, length_col: As for SegmentLattice.
            sites (list): Network nodes of the sites, in the order of the
                          matrices returned.
            cell_km (float): Target distance between neighbouring points
                             (the coalescence times need a system the size
                             of the number of pairs of points).
            density (float): Population per km of river, which sets the
                             coalescence rates relative to the movement
                             rates.
        """
        if density <= 0:
            raise ValueError("density must be positive")
        self.density = float(density)
        self.lattice = SegmentLattice(K, edge_order, reachid_col, length_col,
                                      cell_km, keep_nodes=sites)
        self.sites = self.lattice.states(sites)
        lat = self.lattice
        self.around = np.bincount(lat.a, weights=lat.length / 2,
                                  minlength=lat.n)

    def jump_rates(self, reach_r=None, rate=0.0):
        """
        Lineage jump rates on each lattice link (lattice.a -> b).

        Args:
            reach_r (array-like, optional): Resistance per km of each reach,
                                            in edge order (default: 1,
                                            distance alone).
            rate (float): Flow rate per km (0: none).

        Returns:
            numpy.ndarray: Jump rates.
        """
        lat = self.lattice
        W = lat.weights(reach_r, rate)
        w = np.asarray(W.tocsr()[lat.a, lat.b]).ravel()
        return w / self.around[lat.a]

    def coalescence_times(self, reach_r=None, rate=0.0):
        """
        Mean coalescence times among all lattice points.

        Args:
            reach_r, rate: As for jump_rates().

        Returns:
            numpy.ndarray: C, among the lattice's points.
        """
        lat = self.lattice
        G = generator(self.jump_rates(reach_r, rate),
                      list(zip(lat.a, lat.b)), lat.n)
        gamma = 1.0 / (self.density * self.around)
        return coalescence_times(G, gamma)

    def linfst(self, reach_r=None, rate=0.0):
        """
        Linearised FST among the sites from the coalescence times.

        Args:
            reach_r, rate: As for jump_rates().

        Returns:
            numpy.ndarray: Linearised FST (sites in the order given).
        """
        C = self.coalescence_times(reach_r, rate)
        s = self.sites
        Cs = C[np.ix_(s, s)]
        w = np.diag(Cs)
        within = (w[:, None] + w[None, :]) / 2
        lin = (Cs - within) / within
        np.fill_diagonal(lin, 0.0)
        return lin
