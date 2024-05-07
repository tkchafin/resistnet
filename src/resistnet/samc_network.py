import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import scipy.sparse as sp
import momepy
import seaborn as sns
from sortedcontainers import SortedDict

import resistnet.utils as utils
import resistnet.transform as trans
import resistnet.resist_dist as rd
from resistnet.resistance_network import ResistanceNetwork


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
                 origin=None, rtol=0.00001, max_iter=1000, max_fail=1,
                 solver="iterative", out="out", verbose=True, minimize=False):
        """
        Constructs all the necessary attributes for the ResistanceNetwork
        object.

        Args:
            network (NetworkX graph): The network graph.
            shapefile (str): Path to the shapefile.
            sizes (str): Path to file with population sizes.
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
        """
        self.network = network
        self.shapefile = shapefile
        self.pop_sizes = sizes
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
        self.rtol = rtol
        self.max_iter = max_iter
        self.max_fail = max_fail
        self.solver = solver

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
        self._edge_order = None
        self._gendist = None
        self._origin = None
        self._edge_origin = None
        self._reach_origin = None
        self._Kd = None
        self._edge_absorption = None
        self._R = None
        self._edge_site_indices = None

        # Initialization methods
        self.initialize_network()
        self.initialize_predictors()
        self.build_incidence_matrix(self.reachid_col, out=self.output_prefix)
        if self.inmat is not None:
            self.parse_input_gen_mat(out=self.output_prefix)

        # direct graph
        self.polarise_network(self._origin)

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

        # Read population sizes file
        sizes = self.read_sizes()

        # Filter to keep only entries present in both points and sizes
        points = points[points["sample"].isin(sizes["sample"])]
        sizes = sizes[sizes["sample"].isin(points["sample"])]

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
        dense_adj_matrix = self._adj.toarray()
        adj_out = self.output_prefix + ".adjacencyMatrix.tsv"
        np.savetxt(adj_out, dense_adj_matrix, delimiter='\t', fmt='%d')

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
            raise ValueError("Failed to extend confluence using NEXT_DOWN.")

        if not self._K.has_node(next_down_edge[1]):
            self._K.add_node(next_down_edge[1])
        self._K.add_edge(
            next_down_edge[0], next_down_edge[1], **next_down_data)

        # Increment EDGE_ID and add to new edge
        nx.set_edge_attributes(
            self._K,
            {(next_down_edge[0], next_down_edge[1]): {"EDGE_ID": new_edge_id}}
        )

        # Update origin to the new downstream node
        self._origin = next_down_edge[1]
        self._edge_origin = new_edge_id
        self._reach_origin = next_down_data[self.reachid_col]

        # Recursive call to process further increments
        self.extend_confluence(increments - 1)

    def evaluate(self, individual):
        """
        Evaluates the fitness of a given model represented by an individual.

        This method calculates the fitness of an individual model based on its
        variable transformations. It builds a multi-surface representation of
        the variables and evaluates the model using the ResistanceNetwork's
        parsePairwise method.

        Args:
            individual (list): A list representing an individual in the genetic
                               algorithm, containing transformation and weight
                               information for variables.

        Returns:
            tuple: A tuple containing the fitness value and the results of the
                   evaluation.
        """
        multi = self._build_composite_surface(individual)
        if multi is None:
            return (float('-inf'), [np.nan, np.nan, np.nan, np.nan])

        cfpt, res = rd.conditionalFirstPassTime(
            multi, self._R, self._edge_site_indices, self._gendist,
            self.rtol, self.max_iter, self.max_fail, self.solver)

        if cfpt is not None:
            fitness = res[self.fitmetric].iloc[0]
            res = list(res.iloc[0])
            return fitness, res
        return (float('-inf'), [np.nan, np.nan, np.nan, np.nan])

    # will need to overload plotting functions for SAMC models as well
    def model_output(self, model):
        """
        Generates the model output for a given model.

        This method combines the transformed variables based on the model
        specifications to produce a multi-surface representation. It then
        calculates the effective resistance matrix for the given model.

        Args:
            model (pd.Series): A list representing an individual model in the
                          genetic algorithm, containing transformation and
                          weight information for variables.

        Returns:
            tuple: A tuple containing the effective resistance matrix and the
                   multi-surface representation for the model.
        """
        multi = self._build_composite_surface(model)
        if multi is None:
            return None, None

        # re-compute cfpt times
        cfpt, _ = rd.conditionalFirstPassTime(
            multi, self._R, self._edge_site_indices, self._gendist,
            self.rtol, self.max_iter, self.max_fail, self.solver)

        # re-compute transition rates and return as us/ds vectors
        trans = self._transition_vectors(multi)
        return cfpt, trans

    def _transition_vectors(self, mat):
        n = len(self._edge_order)
        # Initialize arrays with NaNs
        downstream_transitions = np.full(n, 0.0)
        upstream_transitions = np.full(n, 0.0)

        # Process the sparse matrix directly
        mat_coo = mat.tocoo()  # coo easier traversal
        for i, j, rate in zip(mat_coo.row, mat_coo.col, mat_coo.data):
            # Check adjacency for direction of flow
            if self._adj[i, j] == 1:  # i -> j is downstream
                downstream_transitions[i] += rate
            elif self._adj[i, j] == -1:  # i -> j is upstream for node j
                upstream_transitions[j] += rate

        # Correctly handling terminal edges w/ no neighbors
        adj_coo = self._adj.tocoo()
        has_downstream = np.zeros(n, dtype=bool)
        has_upstream = np.zeros(n, dtype=bool)

        for i, j, value in zip(adj_coo.row, adj_coo.col, adj_coo.data):
            if value == 1:
                has_downstream[i] = True
            elif value == -1:
                has_upstream[j] = True

        for i in range(n):
            if not has_downstream[i]:
                downstream_transitions[i] = np.nan
            if not has_upstream[i]:
                upstream_transitions[i] = np.nan

        return np.vstack((downstream_transitions, upstream_transitions)).T

    def _build_composite_surface(self, individual):
        """
        Builds a multi-surface representation from a model configuration.

        Args:
            config (list): A list representing an individual or model in the
                        genetic algorithm, containing transformation and weight
                        information for variables.

        Returns:
            ndarray: The multi-surface representation for the configuration.
        """
        first = True
        multi = None
        for i, variable in enumerate(self._predictors.columns):
            var_m = np.zeros_like(self._adj, dtype=float)
            if individual[0::5][i] == 1:
                var = self.transform(
                    self._predictors[variable],
                    individual[2::5][i],
                    individual[3::5][i]
                )
                var = trans.rescaleCols(var, 0, 1)
                if individual[4::5][i] == 1:
                    var_m = self._compute_transition(var, directional=True)
                    var_m = var_m.T
                else:
                    var_m = self._compute_transition(var, directional=False)
                var_m.data = utils.minmax(var_m.data)
                if first:
                    multi = var_m.copy()
                    multi.data *= individual[1::5][i]
                    first = False
                else:
                    var_m.data *= individual[1::5][i]
                    multi.data += var_m.data

        if multi is not None:
            multi.data = utils.minmax_nonzero(multi.data)
            multi.data = utils.minmax_nonzero(1 / multi.data)

        return multi

    def _compute_transition(self, var, directional=False):
        data = []
        rows = []
        cols = []
        n = self._adj.shape[0]
        var = var.reset_index(drop=True)

        if directional:
            for row, col in zip(*self._adj.nonzero()):
                new_value = var.iloc[row] - var.iloc[col]
                data.append(new_value)
                rows.append(row)
                cols.append(col)
        else:
            for row, col in zip(*self._adj.nonzero()):
                new_value = (var.iloc[row] + var.iloc[col]) / 2.0
                data.append(new_value)
                rows.append(row)
                cols.append(col)

        # Create a sparse matrix from the new values
        result_matrix = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
        return result_matrix

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

    def output_and_plot_model(self, oname, mat_r, edge_r, plot=True):
        """
        Outputs model results to files and optionally plots the resistance
        network and pairwise models.

        This method writes the edge and matrix resistance results to files and,
        if requested, plots the resistance network and pairwise models. The
        plotting is dependent on the availability of genetic distance data.

        Args:
            oname (str): The base name for output files.
            mat_r (numpy.ndarray): The average resistance matrix.
            edge_r (DataFrame): The average edge resistance values.
            plot (bool, optional): Boolean to control whether to generate plots
                                   Defaults to True.
        """
        if mat_r is None or edge_r is None:
            return
        # Write edge and matrix results to files
        out1 = f"{oname}.ResistanceEdges.tsv"
        utils.write_edges(out1, edge_r, self._edge_order)

        out2 = f"{oname}.CFPTMatrix.tsv"
        utils.write_matrix(out2, mat_r, list(self._points_labels.values()))

        # Plot resistance network and pairwise models if required
        if plot:
            edge_r = np.array(edge_r)
            if edge_r.ndim == 2 and edge_r.shape[1] == 2:
                edf = pd.DataFrame({
                    "EDGE_ID": self._edge_order,
                    "Resistance_Downstream": edge_r[:, 0],
                    "Resistance_Upstream": edge_r[:, 1]
                })
            elif edge_r.ndim == 1 or (
                    edge_r.ndim == 2 and edge_r.shape[1] == 1):
                edf = pd.DataFrame({
                    "EDGE_ID": self._edge_order,
                    "Resistance": edge_r.flatten()
                })
            else:
                raise ValueError("Unexpected shape or type for edge_r data.")
            print(edf)
            self.plot_resistance_network(edf, oname)

            # Plot pairwise model if genetic distance data is available
            if self._gendist is not None:
                self.plot_pairwise_model(mat_r, oname, partition=False)

    def plot_resistance_network(self, resistance_values, oname):
        """
        Plots the resistance network based on given resistance values.

        Args:
            resistance_values (DataFrame): A DataFrame containing resistance
                values, which should correspond to the 'EDGE_ID' in the graph.
            oname (str): The output name for the saved plot file.
        """

        # Convert the graph to a GeoDataFrame
        _, edgeDF = momepy.nx_to_gdf(self._K)
        edgeDF['EDGE_ID'] = edgeDF[self.reachid_col].astype(int)
        geoDF = edgeDF.merge(resistance_values, on="EDGE_ID")
        print(geoDF)
        # Plotting the resistance network
        sns.set(style="ticks")
        if 'Resistance_Upstream' in geoDF.columns:
            fig, ax = plt.subplots(1, 2, figsize=(20, 10))
            geoDF.plot(ax=ax[0], column='Resistance_Downstream',
                       cmap="RdYlGn_r", legend=True)
            ax[0].set_title('Downstream Resistance')
            geoDF.plot(ax=ax[1], column='Resistance_Upstream',
                       cmap="RdYlGn_r", legend=True)
            ax[1].set_title('Upstream Resistance')
            fig.suptitle("Stream Network Resistance")
        else:
            geoDF.plot(column='Resistance', cmap="RdYlGn_r", legend=True)
            plt.title("Stream network colored by resistance")

        plt.savefig(f"{oname}.streamsByResistance.pdf")
        plt.clf()
        plt.close()


class ResistanceNetworkSAMCWorker(ResistanceNetworkSAMC):
    """
    A subclass of ResistanceNetwork that handles specific worker operations.

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
        posWeight (bool): Position weight.
        fixWeight (bool): Fixed weight.
        allShapes (bool): Indicator for all shapes.
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
        adj (numpy.ndarray): Adjacency matrix.
        origin (tuple): Origin node.
        R (numpy.ndarray): Absorption array R
    """
    def __init__(self, network, pop_agg, reachid_col, length_col,
                 variables, agg_opts, fitmetric, posWeight, fixWeight,
                 allShapes, fixShape, min_weight, max_shape, inc, point_coords,
                 points_names, points_snapped, points_labels, predictors,
                 edge_order, gendist, adj, origin, R, allSymmetric,
                 edge_site_indices, rtol=0.00001, max_iter=1000, max_fail=1,
                 solver="iterative"):

        self.network = network
        self.pop_agg = pop_agg
        self.reachid_col = reachid_col
        self.length_col = length_col
        self.output_prefix = None
        self.verbose = False
        self.variables = variables
        self.agg_opts = agg_opts

        self.fitmetric = fitmetric
        self.posWeight = posWeight
        self.fixWeight = fixWeight
        self.fixShape = fixShape
        self.allShapes = allShapes
        self.min_weight = min_weight
        self.max_shape = max_shape
        self.allSymmetric = allSymmetric
        self.rtol = rtol
        self.max_iter = max_iter
        self.max_fail = max_fail
        self.solver = solver

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
