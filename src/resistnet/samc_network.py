import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import scipy.sparse as sp
import momepy
import seaborn as sns
from scipy.optimize import minimize_scalar
from sortedcontainers import SortedDict

import resistnet.utils as utils
import resistnet.transform as trans
import resistnet.resist_dist as rd
import resistnet.CFPT as cfpt_samc
from resistnet.resistance_network import ResistanceNetwork

# Smallest non-zero flow rate searched (per km)
MIN_FLOW = 1e-3


class ResistanceNetworkSAMC(ResistanceNetwork):
    """
    A class to represent a resistance network and associated data for SAMC

    Attributes:
        network (NetworkX graph): The network graph.
        shapefile (str): Path to the shapefile.
        coords (list): List of coordinates.
        variables (list): List of variables.
        inmat (str): Path to input matrix file.
        pop_agg (str): Population aggregation method. Defaults to "ARITH".
        reachid_col (str): Column name for reach ID. Defaults to "HYRIV_ID".
        length_col (str): Column name for length. Defaults to "LENGTH_KM".
        output_prefix (str): Prefix for output files. Defaults to "out".
        verbose (bool): Flag to enable verbose output. Defaults to True.
        minimize (bool): Flag to minimize subgraph. Defaults to False.
        agg_opts (dict): Aggregation options.

    Methods:
        initialize_network(): Initializes the network.
        initialize_predictors(): Initializes predictors.
        build_incidence_matrix(): Builds the incidence matrix.
        parse_input_gen_mat(): Parses the input genetic distance matrix.
    """

    def __init__(self, network, shapefile, sizes, coords, variables, inmat,
                 agg_opts, pop_agg="ARITH", reachid_col="HYRIV_ID",
                 length_col="LENGTH_KM", infer_origin="NEXT_DOWN",
                 origin=None, out="out", verbose=True, minimize=False,
                 ne=100, max_resist=500, fix_resist=None, max_flow=1.0):
        """
        Constructs all the necessary attributes for the ResistanceNetwork
        object.

        Args:
            network (NetworkX graph): The network graph.
            shapefile (str): Path to the shapefile.
            sizes (str): Path to file with effective population sizes
                         (columns sample, size); if None, all sites have
                         size ne.
            coords (list): List of coordinates.
            variables (list): List of variables.
            inmat (str): Path to input matrix file.
            agg_opts (dict): Aggregation options.
            pop_agg (str): Population aggregation method. Defaults to "ARITH".
            reachid_col (str): Column name for reach ID. Defaults to "HYRIV_ID"
            length_col (str): Column name for length. Defaults to "LENGTH_KM".
            infer_origin (str): Column name used to infer origin.
                                Default "NEXT_DOWN"
            origin (int): reachid corresponding to the origin.
                            Defaults to None.
            out (str): Prefix for output files. Defaults to "out".
            verbose (bool): Flag to enable verbose output. Defaults to True.
            minimize (bool): Flag to minimize subgraph. Defaults to False.
            ne (float): Effective population size of every site, used if
                        sizes is None. Defaults to 100.
            max_resist (float): Upper limit when optimising the resistance
                        of the most resistant moves, relative to the least
                        resistant (1). Defaults to 500.
            fix_resist (float, optional): Fix that maximum resistance
                        instead of optimising it.
            max_flow (float): Upper limit of the flow term's rate (per km;
                        see flow_rate()). Defaults to 1.
        """
        self.network = network
        self.shapefile = shapefile
        self.pop_sizes = sizes
        self.ne = ne
        self.coords = coords
        self.inmat = inmat
        self.pop_agg = pop_agg
        self.reachid_col = reachid_col
        self.length_col = length_col
        self.output_prefix = out
        self.verbose = verbose
        self.minimize = minimize
        self.variables = variables
        self.agg_opts = agg_opts
        self.fitmetric = "aic"
        self.infer_origin = infer_origin
        self.origin = origin
        self.allSymmetric = False
        self.num_attributes = 5
        self.max_resist = max_resist
        self.fix_resist = fix_resist
        self.max_flow = max_flow

        # Additional attributes for internal use
        self.minimized_subgraph = None
        self.subgraph = None
        self._inc = None
        self._G = None
        self._K = None
        self._sizes = None
        self._point_coords = None
        self._points_names = None
        self._points_snapped = None
        self._points_labels = None
        self._predictors = None
        self._predictor_range = None
        self._edge_order = None
        self._gendist = None
        self._origin = None
        self._edge_origin = None
        self._reach_origin = None
        self._Kd = None
        self._edge_absorption = None
        self._R = None
        self._edge_site_indices = None
        self._flow_pos = None

        # Initialization methods
        self.initialize_network()
        self.initialize_predictors()
        self.build_incidence_matrix(self.reachid_col, out=self.output_prefix)
        if self.inmat is not None:
            self.parse_input_gen_mat(out=self.output_prefix)

        # direct graph
        self.polarise_network(self._origin)
        self._flow_pos = self.flow_positions()
        self.plot_upstream_changes(self.output_prefix)

        # calculate edge-wise absorption
        self.calculate_absorption()

    def initialize_network(self):
        """
        Initializes the network by loading and processing the network graph
        from a shapefile.

        This method involves several steps:
        - Reading the network graph from a shapefile.
        - Reading and snapping point coordinate data to the network graph.
        - Computing the subgraph based on the snapped points.
        - Annotating edges of the network for additional information.
        - Cleaning up to save memory by removing unnecessary copies of the
          graph.
        - Setting the path to the pickled network file.

        Depending on whether the network is minimized, different network files
        are set.
        """

        # Read the network graph from a shapefile
        self.read_network()

        # Read and snap point coordinate data
        points = utils.read_points_table(self.coords)

        # Read population sizes file, or use the same size for all sites
        if self.pop_sizes:
            sizes = self.read_sizes()
            missing = points.loc[
                ~points["sample"].isin(sizes["sample"]), "sample"]
            if len(missing):
                print("WARNING: No population size for samples "
                      f"{list(missing)}; these are excluded")
            # Filter to keep only entries present in both points and sizes
            points = points[points["sample"].isin(sizes["sample"])]
            sizes = sizes[sizes["sample"].isin(points["sample"])]
        else:
            sizes = pd.DataFrame(
                {"sample": points["sample"], "size": float(self.ne)})

        # Snap points to the network graph
        self.snap_points(points, self._G, self.output_prefix)

        # NOTE: If multiple pops map to same point, get the average size
        self._sizes = SortedDict()
        for key, pops in self._points_names.items():
            this_sizes = sizes[sizes['sample'].isin(pops)]
            if this_sizes.empty:
                raise ValueError(f"Population in {pops} missing size")
            self._sizes[key] = this_sizes['size'].mean()

        # Compute the subgraph based on the snapped points
        self._K = self.parse_subgraph_from_points(self.output_prefix)

        # Annotate edges of the network graph
        self._K = self.annotate_edges(self.output_prefix)

        # Identify the furthest downstream segment to anchor directional
        # calculations later on
        self.find_confluence()

        # # Modify network to contain further downstream edges
        # # This simplifies some things later on...
        new_edge = max(
            (
                data.get('EDGE_ID', -1)
                for _, _, data in self._K.edges(data=True)
            ),
            default=-1
        ) + 1
        self.extend_confluence(1, new_edge)

        # # add buffer edges around terminal sample points
        # # (makes it easier to plot later)
        samples = list(self._points_snapped.keys())
        self.buffer_points(1, points=samples)

        # Check for disconnected components
        if nx.is_connected(self._K):
            pass
        else:
            print("The network is not connected.")
            components = list(nx.connected_components(self._K))
            print(f"There are {len(components)} disconnected components.")

        # plot graph with detected origin
        self.plot_inferred_origin(self.output_prefix)

        # Clean up to save memory by removing unnecessary graph copies
        del self._G
        del self.minimized_subgraph
        del self.subgraph

        # Set the path to the pickled network file based on whether the
        # network is minimized
        if self.minimize:
            network_file_suffix = ".minimalSubgraph.net"
        else:
            network_file_suffix = ".subgraph.net"
        self.network = f"{self.output_prefix}{network_file_suffix}"

    # NOTE: Future version could add recursively so we make sure we cover
    # all upstream edges in cases of fork
    # this would make it easier to have a quick projection of resistance
    # around the whole sub-network as well
    def buffer_points(self, buffer_size=1, points=None):
        if points is None:
            # Identify all terminal nodes if none specified for extension
            points = [
                node for node in self._K.nodes() if self._K.degree(node) == 1
            ]
            if not points:
                print("No terminal points found in the subgraph.")
                return

        # Get the highest EDGE_ID currently in the network
        max_edge_id = max(
            (
                data.get('EDGE_ID', -1)
                for _, _, data in self._K.edges(data=True)
            ),
            default=-1
        ) + 1

        # Iterate over the provided terminal points
        for node in points:
            if self._K.degree(node) == 1:
                current_neighbor = list(self._K.neighbors(node))[0]

                # Identify potential new neighbors from the full graph
                neighbors = list(self._G.neighbors(node))
                new_neighbors = [
                    n for n in neighbors
                    if n not in self._K.nodes() and n != current_neighbor
                ]
                # Limit neighbors to add based on the buffer_size
                for new_neighbor in new_neighbors[:buffer_size]:
                    if self._G.has_edge(node, new_neighbor):
                        # Retrieve edge data from the full graph
                        edge_data = self._G.get_edge_data(node, new_neighbor)

                        # Add new node if not already in the subgraph
                        if not self._K.has_node(new_neighbor):
                            self._K.add_node(new_neighbor)

                        # Add new edge
                        self._K.add_edge(
                            node, new_neighbor, **edge_data,
                            EDGE_ID=max_edge_id
                        )
                max_edge_id += 1

    # scipy.sparse version
    def create_adjacency_matrix(self):
        """
        Creates an adjacency matrix representing the connections between nodes
        in the network.

        This method constructs a sparse adjacency matrix to represent the graph
        structure of the network. The matrix is filled based on the edge data,
        with entries indicating the directionality of connections (1 for
        downstream, -1 for upstream) between nodes.

        Returns:
            None: The method updates the instance attribute `_adj` with the
            constructed adjacency matrix (scipy.sparse.cnr format)
        """
        n = len(self._edge_order)
        adjacency_matrix = sp.lil_matrix((n, n), dtype=int)
        edge_to_idx = {
            edge_id: idx for idx, edge_id in enumerate(self._edge_order)}

        id = self.reachid_col

        # Helper function to determine adjacency and direction
        def add_edge_adjacency(u, v, edge_data_u_v):
            if id in edge_data_u_v:
                edge_id_u_v = edge_data_u_v[id]
                if edge_id_u_v in edge_to_idx:
                    idx_u_v = edge_to_idx[edge_id_u_v]
                    # Since the graph is oriented towards the root, successors
                    # of v are downstream, predecessors are upstream
                    for succ in self._K.successors(v):  # Downstream edges
                        edge_data_v_succ = self._K.get_edge_data(v, succ)
                        if id in edge_data_v_succ:
                            if edge_data_v_succ[id] in edge_to_idx:
                                idx_v_succ = edge_to_idx[edge_data_v_succ[id]]
                                adjacency_matrix[idx_u_v, idx_v_succ] = 1
                                adjacency_matrix[idx_v_succ, idx_u_v] = -1

        # Iterate over all edges to populate adjacency relationships
        for u, v, edge_data in self._K.edges(data=True):
            add_edge_adjacency(u, v, edge_data)

        # If the outlet could not be extended downstream, several reaches
        # flow into the origin but share no downstream reach, which would
        # split the network; link them directly instead (each is upstream of
        # the other)
        root_edges = [
            edge_to_idx[d[id]]
            for _, _, d in self._K.in_edges(self._origin, data=True)
            if d.get(id) in edge_to_idx
        ]
        if len(root_edges) > 1:
            for idx_i in root_edges:
                for idx_j in root_edges:
                    if idx_i != idx_j:
                        adjacency_matrix[idx_i, idx_j] = -1

        # # special case for edges connected to origin
        # root_edges = self._K.in_edges(self._origin, data=True)
        # root_edge_indices = []
        # for u, v, data in root_edges:
        #     if self.reachid_col in data:
        #         edge_id = data[self.reachid_col]
        #         if edge_id in edge_to_idx:
        #             root_edge_indices.append(edge_to_idx[edge_id])
        # for idx_i in root_edge_indices:
        #     for idx_j in root_edge_indices:
        #         if idx_i != idx_j:
        #             adjacency_matrix[idx_j, idx_i] = -1
        #             adjacency_matrix[idx_i, idx_j] = -1
        # self._root_edge_indices = root_edge_indices

        # If non-directional model, take absolute values of adjacency_matrix
        if self.allSymmetric:
            self._adj = abs(adjacency_matrix).tocsr()
        else:
            self._adj = adjacency_matrix.tocsr()
        return

    def polarise_network(self, origin):

        # Re-orient as a DAG
        self._K = utils.graph_to_dag_converging(self._K, self._origin)
        self.plot_dag_directionality(self.output_prefix)

        # make adjacency matrix of edges, ordered as in self._edge_order
        self.create_adjacency_matrix()
        n_comp, _ = sp.csgraph.connected_components(self._adj, directed=False)
        if n_comp > 1:
            raise ValueError(
                f"Reach adjacency is split into {n_comp} disconnected parts; "
                "SAMC needs a connected network"
            )
        dense_adj_matrix = self._adj.toarray()
        adj_out = self.output_prefix + ".adjacencyMatrix.tsv"
        np.savetxt(adj_out, dense_adj_matrix, delimiter='\t', fmt='%d')

    def flow_positions(self):
        """
        Distance along the flow from the middle of each reach to the outlet of
        the network (in the units of length_col, i.e. km).

        Differences between adjacent reaches are the distance moved upstream
        (positive) or downstream (negative) between their midpoints; up to a
        constant, this is the distance to the sea (DIST_DN_KM in
        HydroRIVERS) at each reach's midpoint.

        Returns:
            numpy.ndarray: Position of each reach, ordered as self._edge_order
        """
        eid = self.reachid_col
        length = {d[eid]: d[self.length_col]
                  for _, _, d in self._K.edges(data=True)}
        L = np.array([float(length[e]) for e in self._edge_order])
        rows, cols = self._adj.nonzero()
        down = np.asarray(self._adj[rows, cols]).ravel() == 1
        nxt = dict(zip(rows[down], cols[down]))

        # Distance from the downstream end of each reach to the outlet
        below = np.full(len(L), np.nan)
        for start in range(len(L)):
            chain = [start]
            while chain[-1] in nxt and np.isnan(below[chain[-1]]):
                chain.append(nxt[chain[-1]])
            last = chain[-1]
            if np.isnan(below[last]):
                below[last] = 0.0
            for s in reversed(chain[:-1]):
                below[s] = below[nxt[s]] + L[nxt[s]]
        return below + L / 2.0

    def calculate_absorption(self):
        # Initialize the array for absorption values, filled with zeros
        n = len(self._edge_order)
        self._edge_absorption = [0] * n

        # Mapping from edge IDs to indices in the absorption array
        edge_to_idx = {
            edge_id: idx for idx, edge_id in enumerate(self._edge_order)}

        sites = list(self._points_snapped.keys())
        tips = [None] * len(sites)
        # Iterate over edges in the graph
        for u, _, edge_data in self._K.edges(data=True):
            # Get the edge ID from edge_data using the identifier column name
            if self.reachid_col in edge_data:
                edge_id = edge_data[self.reachid_col]
                if edge_id in edge_to_idx:
                    # Determine if the 'source' node (u) is sample site
                    if u in sites:
                        # Get size
                        if u in self._sizes:
                            pop_size = self._sizes[u]
                            # Calculate the absorption value for the edge
                            # 1/2Ne (prob 2 alleles coalesce)
                            absorption_value = 1 / (2 * pop_size)
                            idx = edge_to_idx[edge_id]
                            self._edge_absorption[idx] = absorption_value
                            tips[sites.index(u)] = idx
                        else:
                            print(f"Population size missing for {u}")
        missing = [s for s, t in zip(sites, tips) if t is None]
        if missing:
            raise ValueError(
                f"No downstream reach found for sample sites {missing}"
            )
        self._R = np.array(self._edge_absorption)
        self._R = self._R.reshape(-1, 1)
        self._edge_site_indices = tips

    def read_sizes(self):
        """
        Reads the population sizes file and returns it as a DataFrame.

        Returns:
            pd.DataFrame: DataFrame containing 'sample' and 'size' columns.
        """
        if not hasattr(self, 'pop_sizes') or not self.pop_sizes:
            raise ValueError("Population sizes file path is not set.")

        # Attempt to read the file
        try:
            sizes = pd.read_csv(self.pop_sizes, sep='\t', header=0)
        except Exception as e:
            raise ValueError(f"Failed to read population sizes file: {str(e)}")

        # Ensure required columns are present
        required_columns = ['sample', 'size']
        if not all(col in sizes.columns for col in required_columns):
            raise ValueError("Missing required columns in sizes file")

        return sizes

    def find_confluence(self):
        reach_to_edge = self.create_reach_to_edge_dict()
        edges_without_next_down = []

        # Find edges without a NEXT_DOWN
        for p1, p2, dat in self.subgraph.edges(data=True):
            next_down = dat.get(self.infer_origin)
            if next_down is None or next_down not in reach_to_edge:
                edges_without_next_down.append((p1, p2))

        # Determine the origin based on the edges found
        if len(edges_without_next_down) == 1:
            edge = edges_without_next_down[0]
            if self.subgraph.degree(edge[0]) == 1:
                self._origin = edge[0]
            else:
                self._origin = edge[1]
        elif len(edges_without_next_down) > 1:
            shared_nodes = set(edges_without_next_down[0]).intersection(
                set(edges_without_next_down[1]))
            if shared_nodes:
                self._origin = shared_nodes.pop()
            else:
                raise ValueError(
                    "Could not determine a unique confluence node.")
        else:
            raise ValueError("Could not infer origin based on 'NEXT_DOWN'.")

        # Additionally, set the edge origin if not already set
        if self._origin:
            for p1, p2, dat in self.subgraph.edges(data=True):
                if self._origin in [p1, p2]:
                    self._reach_origin = dat[self.reachid_col]
                    self._edge_origin = reach_to_edge[dat[self.reachid_col]]
                    break

    def extend_confluence(self, increments=1, new_edge_id=0):
        if not hasattr(self, '_origin') or self._origin is None:
            raise ValueError("Origin node has not been set.")

        if increments <= 0:
            return  # Base case for recursion: no further increments needed

        # Determine the next downstream reach ID from the current origin
        next_down = None
        for _, _, data in self._K.edges(self._origin, data=True):
            next_down = data.get(self.infer_origin)
            if next_down:
                break

        if next_down is None:
            print("No valid NEXT_DOWN found for origin")
            return  # Base case: No further node downstream

        # Find the corresponding edge in the full graph
        next_down_edge = None
        next_down_data = None
        for u, v, data in self._G.edges(data=True):
            if data.get(self.reachid_col) == next_down:
                next_down_edge = (u, v)
                next_down_data = data
                break

        if next_down_edge is None:
            # e.g. the input network is clipped at the outlet; the reaches
            # meeting at the origin are then linked directly instead (see
            # create_adjacency_matrix)
            print(f"WARNING: Downstream reach {next_down} of the outlet is "
                  "not in the network; linking the reaches that meet at the "
                  "outlet instead")
            return

        # The new origin is the end of the downstream reach away from the
        # current origin (edge endpoints are not stored in flow order)
        u, v = next_down_edge
        if self._origin not in (u, v):
            print(f"WARNING: Downstream reach {next_down} does not meet the "
                  "outlet; linking the reaches that meet at the outlet "
                  "instead")
            return
        new_origin = v if u == self._origin else u

        if not self._K.has_node(new_origin):
            self._K.add_node(new_origin)
        self._K.add_edge(u, v, **next_down_data)

        # Increment EDGE_ID and add to new edge
        nx.set_edge_attributes(
            self._K, {(u, v): {"EDGE_ID": new_edge_id}}
        )

        # Update origin to the new downstream node
        self._origin = new_origin
        self._edge_origin = new_edge_id
        self._reach_origin = next_down_data[self.reachid_col]

        # Recursive call to process further increments
        self.extend_confluence(increments - 1, new_edge_id + 1)

    def evaluate(self, individual):
        """
        Evaluates the fitness of a given model represented by an individual.

        The model's transition matrix is used to compute conditional first
        passage times (CFPT) among sites, and the mean of the two directions
        between each pair of sites is fitted to the genetic distances with
        MLPE.

        Args:
            individual (list): Model specification, with self.num_attributes
                               values (selection, weight, transformation,
                               shape, direction) per variable, optionally
                               followed by the maximum resistance and the
                               flow rate (see flow_rate()). Direction is 0
                               (symmetric) or 1 (moving towards higher
                               transformed values is harder; a decreasing
                               transformation reverses this for the raw
                               values, see effective_directions()).

        Returns:
            tuple: A tuple containing the fitness value and the results of the
                   evaluation.
        """
        Q = self._transition_matrix(individual)
        if Q is None:
            return (float('-inf'), [np.nan, np.nan, np.nan, np.nan])

        _, res = rd.conditionalFirstPassTime(
            Q, self._edge_site_indices, self._gendist, self._R.ravel())
        if res is None:
            return (float('-inf'), [np.nan, np.nan, np.nan, np.nan])

        # Penalise AIC for the surface parameters (aic is negated, so higher
        # is better)
        k = self.surface_parameters(individual)
        res["aic"] -= 2 * k
        res["delta_aic_null"] -= 2 * k
        fitness = res[self.fitmetric].iloc[0]
        res = list(res.iloc[0])
        return fitness, res

    def surface_parameters(self, individual):
        """
        Counts the surface parameters estimated for a model (see
        ResistanceNetwork.surface_parameters), plus the direction setting of
        each selected variable unless fixAsym, the maximum resistance if it
        was optimised, and the flow rate if the model has a flow term.

        Args:
            individual (list): Model specification, as for evaluate().

        Returns:
            int: The number of surface parameters.
        """
        n_var = len(self._predictors.columns) * self.num_attributes
        n_sel = sum(
            1 for s in individual[0:n_var:self.num_attributes] if s == 1
        )
        k = super().surface_parameters(individual[:n_var])
        if not getattr(self, "fixAsym", False):
            k += n_sel
        # The maximum resistance, if it was optimised
        if len(individual) > n_var and self.fix_resist is None and n_sel:
            k += 1
        if self.flow_rate(individual) != 0:
            k += 1
        return k

    def flow_rate(self, individual):
        """
        The rate of the flow term: the value following the maximum resistance
        in the individual, or 0 (no flow term).

        Each move's weight is multiplied by exp(-rate * distance moved
        upstream in km, between reach midpoints; see flow_positions()), so a
        positive rate makes moving upstream harder and a negative rate makes
        it easier. The rate is per km, so it is not rescaled to the network.

        Args:
            individual (list): Model specification, as for evaluate().

        Returns:
            float: The flow rate (per km).
        """
        n_var = len(self._predictors.columns) * self.num_attributes
        if len(individual) > n_var + 1:
            return float(individual[n_var + 1])
        return 0.0

    def max_resistance(self, individual):
        """
        The maximum resistance for a model: the value following the variable
        attributes in the individual if present (when it is optimised),
        otherwise fix_resist, or 100 (capped at max_resist).

        Args:
            individual (list): Model specification, as for evaluate().

        Returns:
            float: The maximum resistance.
        """
        n_var = len(self._predictors.columns) * self.num_attributes
        if len(individual) > n_var:
            return float(individual[n_var])
        if self.fix_resist is not None:
            return float(self.fix_resist)
        return min(100.0, float(self.max_resist))

    def evaluate_null_model(self, oname=None):
        """
        Evaluates the reference models: distance-only and population effects
        only (see ResistanceNetwork.evaluate_null_model), and two SAMC
        models: distance alone (equal movement weights between all adjacent
        reaches; only the MLPE parameters are counted in AIC) and flow alone
        (the same, plus the flow term, with its rate fitted; see
        flow_rate()).

        Args:
            oname (str, optional): The output name for saving the results.

        Returns:
            DataFrame: Log likelihood, AIC, marginal R-squared, delta AIC
                       relative to the population-effects-only model, and
                       the fitted flow rate (per km).
        """
        table = super().evaluate_null_model()
        if table is None:
            return None

        rate, res = self._fit_flow_alone()
        rows = []
        for name, r, k, flow in [
            ("samc_distance", self._flow_alone(0.0), 0, np.nan),
            ("samc_flow", res, 1, rate),
        ]:
            r = r.iloc[0]
            rows.append(pd.DataFrame({
                "loglik": [r["loglik"]],
                "aic": [-r["aic"] + 2 * k],
                "r2m": [r["r2m"]],
                "delta_aic_null": [r["delta_aic_null"] - 2 * k],
                "flow": [flow]
            }, index=[name]))
        table = pd.concat([table.loc[["distance_only"]], *rows,
                           table.loc[["null"]]])

        if oname:
            table.to_csv(f"{oname}.Null-Model.tsv", sep="\t")
        return table

    def _flow_alone(self, rate):
        """
        MLPE results for the model with equal movement weights and a flow
        term with the given rate (None if the fit failed).
        """
        W = abs(self._adj).astype(float)
        if rate != 0:
            rows, cols = W.nonzero()
            upstream = self._flow_pos[cols] - self._flow_pos[rows]
            W = sp.csr_matrix(
                (np.exp(-rate * upstream), (rows, cols)), shape=W.shape)
        Q = cfpt_samc.transition_matrix(W, self._R.ravel())
        _, res = rd.conditionalFirstPassTime(
            Q, self._edge_site_indices, self._gendist, self._R.ravel())
        return res

    def _fit_flow_alone(self):
        """
        Fits the flow rate of the flow-alone model: a grid over both
        directions (log-spaced, up to max_flow), refined around the best.

        Returns:
            tuple: The fitted rate, and the MLPE results.
        """
        def loglik(rate):
            res = self._flow_alone(rate)
            return -np.inf if res is None else res["loglik"].iloc[0]

        mags = np.logspace(np.log10(MIN_FLOW), np.log10(self.max_flow), 13)
        grid = np.concatenate([-mags[::-1], [0.0], mags])
        lls = np.array([loglik(r) for r in grid])
        i = int(np.argmax(lls))
        lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
        best = grid[i]
        if hi > lo:
            opt = minimize_scalar(lambda r: -loglik(r), bounds=(lo, hi),
                                  method="bounded",
                                  options={"xatol": 1e-4 * self.max_flow})
            if -opt.fun > lls[i]:
                best = float(opt.x)
        return best, self._flow_alone(best)

    def model_output(self, model):
        """
        Generates the model output for a given model.

        Args:
            model (list): Model specification, as for evaluate().

        Returns:
            tuple: The CFPT matrix among sites (row = origin, column =
                   destination), and the probabilities of moving downstream
                   and upstream from each reach (columns of an array ordered
                   as self._edge_order).
        """
        Q = self._transition_matrix(model)
        if Q is None:
            return None, None
        cfpt = cfpt_samc.CFPT(Q, self._edge_site_indices, self._R.ravel())
        return cfpt, self._transition_vectors(Q)

    def _transition_vectors(self, Q):
        """
        Probability of moving downstream and upstream from each reach.

        Args:
            Q (scipy.sparse matrix): Transient transition matrix.

        Returns:
            numpy.ndarray: Array with columns (downstream, upstream); NaN
                           where a reach has no neighbour in that direction.
        """
        adj = self._adj.tocsr()
        Q = sp.csr_matrix(Q)
        down = np.asarray(Q.multiply(adj == 1).sum(axis=1)).ravel()
        up = np.asarray(Q.multiply(adj == -1).sum(axis=1)).ravel()
        down[np.asarray((adj == 1).sum(axis=1)).ravel() == 0] = np.nan
        up[np.asarray((adj == -1).sum(axis=1)).ravel() == 0] = np.nan
        return np.column_stack((down, up))

    def _transition_matrix(self, individual):
        """
        Builds the transient transition matrix Q for a model: movement
        weights from the composite surface, with absorption at the sites.

        Args:
            individual (list): Model specification, as for evaluate().

        Returns:
            scipy.sparse.csr_matrix: Q, or None if no variable is selected.
        """
        multi = self._build_composite_surface(individual)
        if multi is None:
            return None
        return cfpt_samc.transition_matrix(multi, self._R)

    def _build_composite_surface(self, individual):
        """
        Builds the composite movement weights between adjacent reaches.

        For each selected variable (transformed and rescaled to 0-1), the
        resistance of moving from reach i to j is the mean of the two reaches'
        values or, if the variable is directional, the difference (value at
        j - value at i), so that moving towards higher (transformed) values
        is harder. These are rescaled to 0-1 and summed by weight. The
        composite is rescaled to resistances from 1 to the maximum
        resistance (see max_resistance()), and inverted to give the
        movement weights (conductances). A flow term (see
        flow_rate()) then multiplies each weight by a factor depending on
        the distance moved upstream; with no variables selected, the weights
        are 1 before the flow term.

        Args:
            individual (list): Model specification, as for evaluate().

        Returns:
            scipy.sparse.csr_matrix: Movement weights, non-zero for adjacent
                                     reaches; None if no variable is selected
                                     and there is no flow term.
        """
        n_att = self.num_attributes
        rows, cols = self._adj.nonzero()
        multi = None
        for i, variable in enumerate(self._predictors.columns):
            if individual[i * n_att] != 1:
                continue
            var = self.transform(
                self._predictors[variable],
                individual[i * n_att + 2],
                individual[i * n_att + 3]
            )
            var = trans.rescaleCols(var, 0, 1).to_numpy(dtype=float)
            if individual[i * n_att + 4] == 1:
                link = var[cols] - var[rows]
            else:
                link = (var[rows] + var[cols]) / 2.0
            link = _rescale01(link) * individual[i * n_att + 1]
            multi = link if multi is None else multi + link

        rate = self.flow_rate(individual)
        if multi is None:
            if rate == 0:
                return None
            weights = np.ones(len(rows))
        else:
            max_r = self.max_resistance(individual)
            resistance = 1.0 + (max_r - 1.0) * _rescale01(multi)
            weights = 1.0 / resistance
        if rate != 0:
            upstream = self._flow_pos[cols] - self._flow_pos[rows]
            weights = weights * np.exp(-rate * upstream)
        return sp.csr_matrix((weights, (rows, cols)), shape=self._adj.shape)

    def model_spec(self, individual):
        """
        Everything needed to compute a model's movement weights for the
        reaches of another network (see resistnet.predict): for each
        selected variable, as for ResistanceNetwork.model_spec() plus its
        direction setting, the range of its link values, and its effective
        direction (see effective_directions()); and the range of the
        composite surface, the maximum resistance and the flow rate.

        Args:
            individual (list): Model specification, as for evaluate().

        Returns:
            dict: The model specification (JSON-serialisable).
        """
        n_att = self.num_attributes
        rows, cols = self._adj.nonzero()
        directions = self.effective_directions(individual)
        spec = {"type": "samc", "variables": {}}
        multi = None
        for i, v in enumerate(self._predictors.columns):
            if individual[i * n_att] != 1:
                continue
            var = self._variable_spec(v, individual[i * n_att + 2],
                                      individual[i * n_att + 3])
            var["weight"] = float(individual[i * n_att + 1])
            var["asym"] = int(individual[i * n_att + 4] == 1)
            x = self._spec_values(v, var)
            link = (x[cols] - x[rows] if var["asym"]
                    else (x[rows] + x[cols]) / 2.0)
            var["link_min"] = float(link.min())
            var["link_max"] = float(link.max())
            if v in directions:
                var["higher_harder"] = directions[v]
            spec["variables"][v] = var
            link = _rescale01(link) * var["weight"]
            multi = link if multi is None else multi + link
        if multi is not None:
            spec["composite_min"] = float(multi.min())
            spec["composite_max"] = float(multi.max())
            spec["max_resistance"] = self.max_resistance(individual)
        spec["flow"] = self.flow_rate(individual)
        return spec

    def effective_directions(self, individual):
        """
        The effective direction of each directional variable in its raw
        values: the share of moves between reaches with different values on
        which moving towards the higher value is harder (by that variable).
        The direction setting applies to the transformed values, so a
        decreasing transformation reverses it, and a humped one (Ricker)
        gives a mix.

        Args:
            individual (list): Model specification, as for evaluate().

        Returns:
            dict: Share for each selected directional variable.
        """
        n_att = self.num_attributes
        rows, cols = self._adj.nonzero()
        out = {}
        for i, v in enumerate(self._predictors.columns):
            if individual[i * n_att] != 1 or individual[i * n_att + 4] != 1:
                continue
            raw = self._predictors[v].to_numpy(dtype=float)
            x = self.transform(self._predictors[v], individual[i * n_att + 2],
                               individual[i * n_att + 3]).to_numpy(dtype=float)
            rises = np.sign(raw[cols] - raw[rows])
            moved = rises != 0
            if not moved.any():
                continue
            harder = np.sign(x[cols] - x[rows])[moved] == rises[moved]
            out[v] = float(np.mean(harder))
        return out

    def plot_inferred_origin(self, oname):
        # Extract the positions of the nodes in the graph
        pos = {n: n for n in self._K.nodes}

        # Create a new figure
        plt.figure(figsize=(8, 6))

        # Draw the edges of the network
        nx.draw_networkx_edges(self._K, pos, alpha=0.5, width=1)

        # def check_terminal(p1, p2):
        #     if self._K.degree(p1) == 1:
        #         return p1
        #     elif self._K.degree(p2) == 1:
        #         return p2
        #     else:
        #         return None

        # highlight the origin node
        nx.draw_networkx_nodes(
            self._K, pos, nodelist=[self._origin],
            node_size=50, node_color='orange')

        # Remove the axis
        plt.axis('off')

        # Save the plot to a file
        plt.title("Graph with inferred origin point")
        plt.savefig(f"{oname}.graphOrigin.pdf")
        plt.close()

    def plot_dag_directionality(self, oname):
        # Extract the positions of the nodes in the graph
        pos = {n: (n[0], n[1]) for n in self._K.nodes}

        # Create a new figure
        plt.figure(figsize=(10, 8))

        # Draw the DAG with arrows indicating direction
        nx.draw_networkx(self._K, pos, node_color='lightblue',
                         with_labels=False, arrows=True, alpha=0.6,
                         node_size=10)

        # Highlight the origin node in the DAG
        nx.draw_networkx_nodes(self._K, pos, nodelist=[self._origin],
                               node_size=30, node_color='orange')

        # Remove the axis
        plt.axis('off')
        plt.title("Directed Acyclic Graph (DAG) with Directionality")

        # Save the plot to a file
        plt.savefig(f"{oname}.graphDirectionality.pdf")
        plt.close()

    def upstream_changes(self):
        """
        How each variable changes on a step upstream: the share of upstream
        moves between adjacent reaches on which it rises, stays the same, or
        falls (moves between the reaches meeting at a clipped outlet are
        neither up nor down, so they are left out).

        A variable that changes the same way on most upstream steps (e.g.
        upstream area, discharge) gives an asymmetric effect that overlaps
        with the flow term; one that goes either way gives only local
        asymmetries.

        Returns:
            DataFrame: Columns rises, no_change, falls (shares), and steps
                       (number of upstream moves), one row per variable.
        """
        rows, cols = self._adj.nonzero()
        sign = np.asarray(self._adj[rows, cols]).ravel()
        back = np.asarray(self._adj[cols, rows]).ravel()
        up = (sign == -1) & (back != -1)
        out = {}
        for var in self._predictors.columns:
            v = self._predictors[var].to_numpy(dtype=float)
            d = v[cols][up] - v[rows][up]
            out[var] = {"rises": np.mean(d > 0), "no_change": np.mean(d == 0),
                        "falls": np.mean(d < 0), "steps": int(up.sum())}
        return pd.DataFrame(out).T

    def plot_upstream_changes(self, oname):
        """
        Writes upstream_changes() to {oname}.upstreamChanges.tsv and plots
        it to {oname}.upstreamChanges.pdf.

        Args:
            oname (str): The output prefix.
        """
        table = self.upstream_changes()
        table.to_csv(f"{oname}.upstreamChanges.tsv", sep="\t",
                     index_label="variable")
        fig, ax = plt.subplots(figsize=(8, 0.5 * len(table) + 1.5))
        left = np.zeros(len(table))
        for col, label, color in [("rises", "rises", "#2a78d6"),
                                  ("no_change", "no change", "#b0aea6"),
                                  ("falls", "falls", "#eb6834")]:
            ax.barh(table.index, table[col], left=left, label=label,
                    color=color, height=0.6, edgecolor="white")
            left += table[col].to_numpy(dtype=float)
        ax.set_xlim(0, 1)
        ax.invert_yaxis()
        ax.set_xlabel("Share of upstream steps")
        ax.set_title("How each variable changes on a step upstream")
        ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.15),
                  frameon=False)
        fig.tight_layout()
        fig.savefig(f"{oname}.upstreamChanges.pdf")
        plt.close(fig)

    def output_and_plot_model(self, oname, mat_r, edge_r, plot=True):
        """
        Outputs model results to files and optionally plots the transition
        probabilities and pairwise model.

        Args:
            oname (str): The base name for output files.
            mat_r (numpy.ndarray): CFPT matrix among sites (row = origin,
                                   column = destination).
            edge_r (numpy.ndarray): Probabilities of moving downstream and
                                    upstream from each reach (columns).
            plot (bool, optional): Boolean to control whether to generate plots
                                   Defaults to True.
        """
        if mat_r is None or edge_r is None:
            return
        edge_r = np.asarray(edge_r)
        edf = pd.DataFrame({
            "EDGE_ID": self._edge_order,
            "Downstream": edge_r[:, 0],
            "Upstream": edge_r[:, 1]
        })
        edf.to_csv(f"{oname}.TransitionEdges.tsv", sep="\t", index=False,
                   na_rep="NA")

        out2 = f"{oname}.CFPTMatrix.tsv"
        utils.write_matrix(out2, mat_r, list(self._points_labels.values()))

        if plot:
            self.plot_resistance_network(edf, oname)

            # Plot pairwise model if genetic distance data is available
            if self._gendist is not None:
                self.plot_pairwise_model(
                    (mat_r + mat_r.T) / 2.0, oname, partition=False)

    def plot_resistance_network(self, resistance_values, oname):
        """
        Plots the probabilities of moving downstream and upstream from each
        reach.

        Args:
            resistance_values (DataFrame): A DataFrame with 'EDGE_ID',
                'Downstream' and 'Upstream' columns.
            oname (str): The output name for the saved plot file.
        """
        _, edgeDF = momepy.nx_to_gdf(self._K)
        edgeDF['EDGE_ID'] = edgeDF[self.reachid_col].astype(int)
        geoDF = edgeDF.merge(resistance_values, on="EDGE_ID")

        sns.set(style="ticks")
        fig, ax = plt.subplots(1, 2, figsize=(20, 10))
        for a, col in zip(ax, ["Downstream", "Upstream"]):
            geoDF.plot(ax=a, column=col, cmap="RdYlGn", legend=True)
            a.set_title(f"Probability of moving {col.lower()}")
        fig.suptitle("Stream network transition probabilities")
        plt.savefig(f"{oname}.streamsByTransition.pdf")
        plt.clf()
        plt.close()


def _rescale01(x):
    """
    Rescales an array to 0-1 (all zeros if it is constant).
    """
    x = np.asarray(x, dtype=float)
    span = x.max() - x.min()
    if span == 0:
        return np.zeros_like(x)
    return (x - x.min()) / span


class ResistanceNetworkSAMCWorker(ResistanceNetworkSAMC):
    """
    A subclass of ResistanceNetworkSAMC that handles specific worker
    operations.

    This class is designed for internal use, to hold local data for worker
    operations.

    Args:
        network (NetworkX graph): The network graph.
        pop_agg (str): Population aggregation method.
        reachid_col (str): Column name for reach ID.
        length_col (str): Column name for length.
        variables (list): List of variables.
        agg_opts (dict): Aggregation options.
        fitmetric (str): Fitness metric.
        fixWeight (bool): Fixed weight.
        fixShape (bool): Fixed shape.
        min_weight (float): Minimum weight.
        max_shape (float): Maximum shape.
        inc (numpy.ndarray): Incidence matrix.
        point_coords (dict): Point coordinates.
        points_names (dict): Point names.
        points_snapped (dict): Snapped points.
        points_labels (dict): Point labels.
        predictors (DataFrame): Predictor variables.
        edge_order (list): Order of edges.
        gendist (numpy.ndarray): Genetic distance matrix.
        adj (scipy.sparse matrix): Reach adjacency matrix.
        origin (tuple): Origin node.
        R (numpy.ndarray): Absorption array R
        allSymmetric (bool): Ignore flow direction in the adjacency matrix.
        edge_site_indices (list): Reach index of each site.
        fixAsym (bool): Directional variables are not searched.
        max_resist (float): Upper limit of the maximum resistance.
        fix_resist (float): Fixed maximum resistance (None to optimise).
        flow_pos (numpy.ndarray): Position of each reach along the flow (see
                                  flow_positions()).
        max_flow (float): Upper limit of the flow rate.
    """
    def __init__(self, network, pop_agg, reachid_col, length_col,
                 variables, agg_opts, fitmetric, fixWeight, fixShape,
                 min_weight, max_shape, inc, point_coords, points_names,
                 points_snapped, points_labels, predictors, edge_order,
                 gendist, adj, origin, R, allSymmetric, edge_site_indices,
                 fixAsym=False, max_resist=500, fix_resist=None,
                 flow_pos=None, max_flow=1.0):

        self.network = network
        self.pop_agg = pop_agg
        self.reachid_col = reachid_col
        self.length_col = length_col
        self.output_prefix = None
        self.verbose = False
        self.variables = variables
        self.agg_opts = agg_opts

        self.fitmetric = fitmetric
        self.fixWeight = fixWeight
        self.fixShape = fixShape
        self.fixAsym = fixAsym
        self.min_weight = min_weight
        self.max_shape = max_shape
        self.allSymmetric = allSymmetric
        self.num_attributes = 5
        self.max_resist = max_resist
        self.fix_resist = fix_resist
        self.max_flow = max_flow

        self._inc = inc
        self._origin = origin
        self._G = None
        self._K = self.read_network()
        self._point_coords = point_coords
        self._points_names = points_names
        self._points_snapped = points_snapped
        self._points_labels = points_labels
        self._predictors = predictors
        self._edge_order = edge_order
        self._gendist = gendist
        self._adj = adj
        self._R = R
        self._edge_site_indices = edge_site_indices
        self._flow_pos = flow_pos
