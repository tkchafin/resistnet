"""
A river network as a lattice of points, for SAMC.

samc (Marx et al. 2020) works on raster cells of equal size: a walker
moves between neighbouring cells with conductance 1 / (mean resistance of
the two cells x distance between their centres). The network equivalent
used here: states are the nodes of the network (where reaches meet, and
the sites, which are snapped to nodes) plus points splitting each reach
into pieces of about cell_km, so that every step covers about the same
distance. Each state has a resistance per km: its reach's, for points
inside a reach; for a node, the mean over the half-pieces of river around
it (as a raster cell over a confluence), and links take samc's
conductance, 2 / ((r_i + r_j) x length).
"""
import numpy as np
import scipy.sparse as sp


class Lattice:
    """
    The point lattice of a river network directed towards its outlet.

    Attributes:
        n (int): Number of states.
        nodes (list): Network nodes, which are states 0 ... len(nodes) - 1.
        pos (numpy.ndarray): Distance of each state to the outlet along the
                             river (km).
        a, b (numpy.ndarray): Links (both directions: a -> b).
        length (numpy.ndarray): Length of each link (km).
        link_reach (numpy.ndarray): Index (in edge order) of the reach each
                                    link lies on.
        state_reach (numpy.ndarray): Reach index of each state inside a
                                     reach; -1 for nodes.
        reach_states (list): For each reach, its inner states, from upstream
                             to downstream.
        reach_link (numpy.ndarray): Length of the links of each reach.
    """

    def __init__(self, K, edge_order, reachid_col, length_col, cell_km=1.0):
        """
        Args:
            K (networkx.DiGraph): The network, edges directed downstream,
                                  with one outlet (a node with no
                                  successor).
            edge_order (list): Reach IDs, in the order of the reaches'
                               values (e.g. ResistanceNetwork._edge_order).
            reachid_col (str): Edge attribute with the reach ID.
            length_col (str): Edge attribute with the reach length (km).
            cell_km (float): Target distance between neighbouring states.
        """
        if cell_km <= 0:
            raise ValueError("cell_km must be positive")
        self.cell_km = float(cell_km)
        self.nodes = list(K.nodes)
        ix = {v: i for i, v in enumerate(self.nodes)}
        reach_ix = {e: i for i, e in enumerate(edge_order)}
        outlets = [ix[v] for v in self.nodes if K.out_degree(v) == 0]
        if len(outlets) != 1:
            raise ValueError(f"The network has {len(outlets)} outlets; "
                             "the lattice needs one")
        down = {}
        for u, v, d in K.edges(data=True):
            down[ix[u]] = (ix[v], float(d[length_col]))
        pos = np.full(len(self.nodes), np.nan)
        pos[outlets[0]] = 0.0
        for i in range(len(self.nodes)):
            chain = [i]
            while np.isnan(pos[chain[-1]]):
                chain.append(down[chain[-1]][0])
            for s in reversed(chain[:-1]):
                pos[s] = pos[down[s][0]] + down[s][1]

        P, state_reach = list(pos), [-1] * len(self.nodes)
        a, b, length, link_reach = [], [], [], []
        self.reach_states = [[] for _ in edge_order]
        self.reach_link = np.full(len(edge_order), np.nan)
        for u, v, d in K.edges(data=True):
            e = reach_ix[d[reachid_col]]
            L = float(d[length_col])
            pieces = max(1, int(round(L / self.cell_km)))
            self.reach_link[e] = L / pieces
            prev = ix[u]
            for k in range(1, pieces + 1):
                if k < pieces:
                    P.append(pos[ix[u]] - k * L / pieces)
                    state_reach.append(e)
                    cur = len(P) - 1
                    self.reach_states[e].append(cur)
                else:
                    cur = ix[v]
                a.append(prev)
                b.append(cur)
                length.append(L / pieces)
                link_reach.append(e)
                prev = cur
        a, b = np.array(a, dtype=int), np.array(b, dtype=int)
        length = np.array(length, dtype=float)
        link_reach = np.array(link_reach, dtype=int)
        self.n = len(P)
        self.pos = np.array(P, dtype=float)
        self.state_reach = np.array(state_reach, dtype=int)
        self.a = np.r_[a, b]
        self.b = np.r_[b, a]
        self.length = np.r_[length, length]
        self.link_reach = np.r_[link_reach, link_reach]
        # Upstream distance moved on each link (negative: downstream)
        self.upstream = self.pos[self.b] - self.pos[self.a]
        # Node states: half of each link touching the node, by reach
        half = self.a < len(self.nodes)
        self._node_halves = (self.a[half], self.link_reach[half],
                             self.length[half] / 2.0)
        self._ix = ix

    def states(self, nodes):
        """State index of each of the given network nodes."""
        return [self._ix[v] for v in nodes]

    def state_resistance(self, reach_r):
        """
        Resistance per km of each state: its reach's for states inside a
        reach; for nodes, the length-weighted mean over the half-links
        around them.

        Args:
            reach_r (array-like): Resistance per km of each reach (edge
                                  order).

        Returns:
            numpy.ndarray: Resistance of each state.
        """
        reach_r = np.asarray(reach_r, dtype=float)
        r = np.empty(self.n)
        inner = self.state_reach >= 0
        r[inner] = reach_r[self.state_reach[inner]]
        node, reach, w = self._node_halves
        num = np.bincount(node, weights=w * reach_r[reach],
                          minlength=len(self.nodes))
        den = np.bincount(node, weights=w, minlength=len(self.nodes))
        r[:len(self.nodes)] = num / den
        return r

    def weights(self, reach_r=None, rate=0.0):
        """
        Movement weights between neighbouring states: samc's conductance,
        2 / ((r_i + r_j) x link length), times exp(-rate x km moved
        upstream) with a flow term.

        Args:
            reach_r (array-like, optional): Resistance per km of each reach
                                            (default: 1, distance alone).
            rate (float): Flow rate per km (0: none).

        Returns:
            scipy.sparse.csr_matrix: Weights (row = from, column = to).
        """
        if reach_r is None:
            r = np.ones(self.n)
        else:
            r = self.state_resistance(reach_r)
        w = 2.0 / ((r[self.a] + r[self.b]) * self.length)
        if rate != 0:
            w = w * np.exp(-rate * self.upstream)
        return sp.csr_matrix((w, (self.a, self.b)), shape=(self.n, self.n))

    def reach_moves(self, rate=0.0):
        """
        Probability of moving downstream and upstream from a state inside
        each reach, away from its ends: both neighbours are on the reach,
        with its resistance, one link away, so only the flow term makes
        them differ: 1 / (1 + exp(-2 x rate x link length)) downstream. NaN
        for reaches of one link (no inner state).

        Args:
            rate (float): Flow rate per km.

        Returns:
            numpy.ndarray: Columns (downstream, upstream), in edge order.
        """
        down = 1.0 / (1.0 + np.exp(-2.0 * rate * self.reach_link))
        inner = np.array([len(s) > 0 for s in self.reach_states])
        down = np.where(inner, down, np.nan)
        return np.column_stack((down, 1.0 - down))


class SegmentLattice:
    """
    A coarser lattice for the coalescent model: points only where reaches
    meet (confluences), at the network's tips, outlet and sites, and about
    cell_km apart along each stretch of river between them, whatever the
    reaches' boundaries. A link's resistance is the sum over the reaches it
    overlaps of their resistance per km times the overlap (resistances in
    series), and its conductance the inverse; with resistance 1 per km,
    1 / link length, as Lattice.

    Attributes:
        As Lattice (n, pos, a, b, length, upstream, reach_link, cell_km),
        plus overlap (scipy.sparse.csr_matrix): km of each reach (columns,
        edge order) covered by each link (rows, as a and b).
    """

    def __init__(self, K, edge_order, reachid_col, length_col, cell_km=5.0,
                 keep_nodes=()):
        """
        Args:
            K, edge_order, reachid_col, length_col: As for Lattice.
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
        each reach's stretch, away from its ends (see Lattice.reach_moves),
        with the length of the link covering the reach's middle.
        """
        down = 1.0 / (1.0 + np.exp(-2.0 * rate * self.reach_link))
        return np.column_stack((down, 1.0 - down))
