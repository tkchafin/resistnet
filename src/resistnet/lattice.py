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

Directional effects multiply each move's weight by exp(-(U_to - U_from)),
the change in a potential U between the two states (samc's transition
functions can depend on the order of the two cells): U = rate x distance
to the outlet is the flow term, and a directional variable's potential
follows its values, which change along the river from state to state as
they would across raster cells (see state_values()).
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
        chains (list): For each edge of the network (a reach, or one part
                       of a reach split into several edges): its reach
                       index, upstream and downstream node states, inner
                       states from upstream to downstream, and link
                       length (km).
        span (float): Distance from the outlet to the furthest state (km).
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
        self.chains = []
        for u, v, d in K.edges(data=True):
            e = reach_ix[d[reachid_col]]
            L = float(d[length_col])
            pieces = max(1, int(round(L / self.cell_km)))
            self.reach_link[e] = L / pieces
            prev, inner = ix[u], []
            for k in range(1, pieces + 1):
                if k < pieces:
                    P.append(pos[ix[u]] - k * L / pieces)
                    state_reach.append(e)
                    cur = len(P) - 1
                    inner.append(cur)
                else:
                    cur = ix[v]
                a.append(prev)
                b.append(cur)
                length.append(L / pieces)
                link_reach.append(e)
                prev = cur
            self.reach_states[e] += inner
            self.chains.append((e, ix[u], ix[v], inner, L / pieces))
        a, b = np.array(a, dtype=int), np.array(b, dtype=int)
        length = np.array(length, dtype=float)
        link_reach = np.array(link_reach, dtype=int)
        self.n = len(P)
        self.pos = np.array(P, dtype=float)
        self.span = float(self.pos.max())
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
        self._interp = self._interpolation()

    def _interpolation(self):
        """The linear map from reach values to state values (see
        state_values()), as a sparse matrix (states x reaches)."""
        n_nodes, n_reach = len(self.nodes), len(self.reach_states)
        node, reach, w = self._node_halves
        den = np.bincount(node, weights=w, minlength=n_nodes)
        nodes = sp.csr_matrix((w / den[node], (node, reach)),
                              shape=(n_nodes, n_reach))
        rows = [nodes]
        own = sp.identity(n_reach, format="csr")
        inner = [None] * (self.n - n_nodes)
        # Each edge is its own chain of states from its upstream node to
        # its downstream node (a reach split into several edges has one
        # chain per part, each with the reach's value at its middle)
        for e, u, v, states, _ in self.chains:
            for k, s in enumerate(states):
                f = (k + 1) / (len(states) + 1)
                if f <= 0.5:
                    t = f / 0.5
                    row = (1 - t) * nodes[u] + t * own[e]
                else:
                    t = (f - 0.5) / 0.5
                    row = (1 - t) * own[e] + t * nodes[v]
                inner[s - n_nodes] = row
        rows += inner
        return sp.vstack(rows).tocsr()

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

    def state_values(self, reach_x):
        """
        A variable's value at each state, changing along the river from
        state to state as across raster cells: each reach's value at its
        middle, each node's the length-weighted mean over the half-links
        around it (as for resistance), and linear in between. The change
        between two neighbouring reaches is so spread over the river
        between their middles.

        Args:
            reach_x (array-like): Value of each reach (edge order).

        Returns:
            numpy.ndarray: Value of each state.
        """
        return self._interp @ np.asarray(reach_x, dtype=float)

    def weights(self, reach_r=None, rate=0.0, potential=None):
        """
        Movement weights between neighbouring states: samc's conductance,
        2 / ((r_i + r_j) x link length), times exp(-rate x km moved
        upstream) with a flow term, and exp(-(U_to - U_from)) with a
        directional potential U.

        Args:
            reach_r (array-like, optional): Resistance per km of each reach
                                            (default: 1, distance alone).
            rate (float): Flow rate per km (0: none).
            potential (array-like, optional): Directional potential of each
                                              state (None: none).

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
        if potential is not None:
            U = np.asarray(potential, dtype=float)
            w = w * np.exp(-(U[self.b] - U[self.a]))
        return sp.csr_matrix((w, (self.a, self.b)), shape=(self.n, self.n))

    def reach_moves(self, rate=0.0, potential=None):
        """
        Probability of moving downstream and upstream from a state inside
        each reach, away from its ends: both neighbours are on the reach,
        with its resistance, one link away, so only the directional effects
        make them differ: 1 / (1 + exp(-2 x rate x link length + U_down -
        U_up)) downstream, averaged over the reach's inner states. NaN for
        reaches of one link (no inner state).

        Args:
            rate (float): Flow rate per km.
            potential (array-like, optional): Directional potential of each
                                              state (None: none).

        Returns:
            numpy.ndarray: Columns (downstream, upstream), in edge order.
        """
        U = None if potential is None else np.asarray(potential, dtype=float)
        vals = [[] for _ in self.reach_states]
        for e, u, v, states, link in self.chains:
            if not states:
                continue
            if U is None:
                vals[e].append(1.0 / (1.0 + np.exp(-2.0 * rate * link)))
                continue
            chain = [u] + list(states) + [v]
            bias = np.array([U[chain[k + 1]] - U[chain[k - 1]]
                             for k in range(1, len(chain) - 1)])
            vals[e].extend(1.0 / (1.0 + np.exp(-2.0 * rate * link + bias)))
        down = np.array([np.mean(v) if v else np.nan for v in vals])
        return np.column_stack((down, 1.0 - down))
