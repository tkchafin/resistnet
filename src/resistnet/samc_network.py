
import sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import scipy.sparse as sp
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
                 origin=None, out="out", verbose=True,
                 minimize=False):
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
        
        # direct and detect origin 
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

        # special case for edges connected to origin
        root_edges = self._K.in_edges(self._origin, data=True)
        root_edge_indices = []
        for u, v, data in root_edges:
            if self.reachid_col in data:
                edge_id = data[self.reachid_col]
                if edge_id in edge_to_idx:
                    root_edge_indices.append(edge_to_idx[edge_id])
        for idx_i in root_edge_indices:
            for idx_j in root_edge_indices:
                if idx_i != idx_j:
                    adjacency_matrix[idx_j, idx_i] = -1
                    adjacency_matrix[idx_i, idx_j] = -1

        # If non-directional model, take absolute values of adjacency_matrix
        if self.allSymmetric:
            self._adj = abs(adjacency_matrix).tocsr()
        else:
            self._adj = adjacency_matrix.tocsr()
        return

    # DEPRECATED: numpy ndarray version
    # def create_adjacency_matrix(self):
    #     n = len(self._edge_order)
    #     adjacency_matrix = np.zeros((n, n), dtype=int)
    #     edge_to_idx = {
    #         edge_id: idx for idx, edge_id in enumerate(self._edge_order)}

    #     id = self.reachid_col

    #     # Helper function to determine adjacency and direction
    #     def add_edge_adjacency(u, v, edge_data_u_v):
    #         if id in edge_data_u_v:
    #             edge_id_u_v = edge_data_u_v[id]
    #             if edge_id_u_v in edge_to_idx:
    #                 idx_u_v = edge_to_idx[edge_id_u_v]
    #                 # Since the graph is oriented towards the root, successors
    #                 # of v are downstream, predecessors are upstream
    #                 for succ in self._K.successors(v):  # Downstream edges
    #                     edge_data_v_succ = self._K.get_edge_data(v, succ)
    #                     if id in edge_data_v_succ:
    #                         if edge_data_v_succ[id] in edge_to_idx:
    #                             idx_v_succ = edge_to_idx[edge_data_v_succ[id]]
    #                             adjacency_matrix[idx_u_v, idx_v_succ] = 1 
    #                             adjacency_matrix[idx_v_succ, idx_u_v] = -1 

    #     # Iterate over all edges to populate adjacency relationships
    #     for u, v, edge_data in self._K.edges(data=True):
    #         add_edge_adjacency(u, v, edge_data)

    #     # special case for edges connected to origin
    #     root_edges = self._K.in_edges(self._origin, data=True)
    #     root_edge_indices = []
    #     for u, v, data in root_edges:
    #         if self.reachid_col in data:
    #             edge_id = data[self.reachid_col]
    #             if edge_id in edge_to_idx:
    #                 root_edge_indices.append(edge_to_idx[edge_id])
    #     for idx_i in root_edge_indices:
    #         for idx_j in root_edge_indices:
    #             if idx_i != idx_j:
    #                 adjacency_matrix[idx_j, idx_i] = -1
    #                 adjacency_matrix[idx_i, idx_j] = -1

    #     # TODO: If directional model, return abs(adjacency_matrix)

    #     self._adj = adjacency_matrix

    def polarise_network(self, origin):

        # Re-orient as a DAG
        self._K = utils.graph_to_dag_converging(self._K, self._origin)
        self.plot_dag_directionality(self.output_prefix)

        # make adjacency matrix of edges, ordered as in self._edge_order
        self.create_adjacency_matrix()
        dense_adj_matrix = self._adj.toarray()
        adj_out = self.output_prefix + ".adjacencyMatrix.tsv"
        np.savetxt(adj_out, dense_adj_matrix, delimiter='\t', fmt='%d')

        # work on further steps later
        pass

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
        for u, v, edge_data in self._K.edges(data=True):
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
            raise ValueError(f"Missing required columns in sizes file")

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

    # def evaluate(self, individual):
    #     """
    #     Evaluates the fitness of a given model represented by an individual.

    #     This method calculates the fitness of an individual model based on its
    #     variable transformations. It builds a multi-surface representation of
    #     the variables and evaluates the model using the ResistanceNetwork's
    #     parsePairwise method.

    #     Args:
    #         individual (list): A list representing an individual in the genetic
    #                            algorithm, containing transformation and weight
    #                            information for variables.

    #     Returns:
    #         tuple: A tuple containing the fitness value and the results of the
    #                evaluation.
    #     """
    #     first = True
    #     multi = None

    #     # Compute any transformations
    #     for i, variable in enumerate(self._predictors.columns):
    #         var_m = np.zeros_like(self._adj, dtype=float)
    #         # Perform variable transformations if the variable is selected
    #         # (indicated by a value of 1)
    #         if individual[0::5][i] == 1:
    #             var = self.transform(
    #                 self._predictors[variable],
    #                 individual[2::5][i],
    #                 individual[3::5][i]
    #             )

    #             # reScale (minmax)
    #             var = trans.rescaleCols(var, 0, 1)

    #             # Compute transition values
    #             # perform directional calculations
    #             if individual[1::5][i] == 1:
    #                 var_m = self._compute_transition(var, directional=True)
    #                 # transpose for 'backward in time' model
    #                 var_m = var_m.T
    #             else:
    #                 var_m = self._compute_transition(var, directional=False)

    #             # var_m = utils.masked_minmax(var_m, mask)
    #             var_m.data = utils.minmax(var_m.data)

    #             # sum within multivariate adjacency
    #             if first:
    #                 multi = var_m.copy()
    #                 multi.data *= individual[1::5][i]
    #                 first = False
    #             else:
    #                 var_m.data *= individual[1::5][i]
    #                 multi.data += var_m.data

    #     # If no layers are selected, return a zero fitness
    #     if multi is None:
    #         return float('-inf'), None
    #     else:
    #         # complete Q matrix
    #         # minmax scale 0-1
    #         multi.data = utils.minmax_nonzero(multi.data)

    #         # inverse to get transition rates
    #         # avoid divide-by-zero by setting zero to smallest non-zero element
    #         multi.data = utils.minmax_nonzero(1 / multi.data)

    #         # compute cfpt matrix
    #         cfpt, res = rd.conditionalFirstPassTime(
    #             multi, self._R, self._edge_site_indices, self._gendist)

    #         # plot cfpt matrix pairwise regression against genetic distance 
    #         # matrix held as self._gendist
    #         # these can be assumed to have the same order
    #         # fitness = res[self.fitmetric][0]
    #         fitness = res[self.fitmetric].iloc[0]
    #         res = list(res.iloc[0])
    #         return (fitness, res)

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
        print("EVALUATE")
        multi = self._build_composite_surface(individual)
        print(multi)
        if multi is None:
            return float('-inf'), None

        cfpt, res = rd.conditionalFirstPassTime(
            multi, self._R, self._edge_site_indices, self._gendist)
        print(cfpt)
        fitness = res[self.fitmetric].iloc[0]
        res = list(res.iloc[0])
        print(res)
        return fitness, res

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
        print("MODEL_OUTPUT")
        multi = self._build_composite_surface(model)
        print(multi)
        if multi is None:
            return None, None

        cfpt, _ = rd.conditionalFirstPassTime(
            multi, self._R, self._edge_site_indices, self._gendist)

        trans = self._transition_vectors(multi)

        print(multi)
        print(trans)
        #sys.exit()

        return cfpt, multi

    def _transition_vectors(self, mat):
        n = len(self._edge_order)
        transition_probs = np.zeros((n, 2))
        mat_csc = mat.tocsc()  # for column access

        # Downstream transitions (i -> j)
        for i in range(n):
            row_start = mat.indptr[i]
            row_end = mat.indptr[i + 1]
            if row_start != row_end:
                transition_probs[i, 0] = mat.data[row_start]

        # Upstream transitions (j -> i)
        for i in range(n):
            col_start = mat_csc.indptr[i]
            col_end = mat_csc.indptr[i + 1]
            for idx in range(col_start, col_end):
                j = mat_csc.indices[idx]
                transition_probs[j, 1] = mat_csc.data[idx]

        return transition_probs

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
                if individual[1::5][i] == 1:
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

        # Iterating through non-zero elements
        if directional:
            for row, col in zip(*self._adj.nonzero()):
                value = self._adj[row, col]
                if value == 1:
                    new_value = var.iloc[col] - var.iloc[row]
                elif value == -1:
                    new_value = var.iloc[row] - var.iloc[col]
                else:
                    continue
                data.append(new_value)
                rows.append(row)
                cols.append(col)
        else:
            for row, col in zip(*self._adj.nonzero()):
                value = self._adj[row, col]
                new_value = (var.iloc[row] + var.iloc[col]) / 2.0
                data.append(new_value)
                rows.append(row)
                cols.append(col)

            # Create a sparse matrix from the new values
            return sp.csr_matrix((data, (rows, cols)), shape=(n, n))

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

    # def output_and_plot_model(self, oname, mat_r, edge_r, plot=True):
    #     """
    #     Outputs model results to files and optionally plots the resistance
    #     network and pairwise models.

    #     This method writes the edge and matrix resistance results to files and,
    #     if requested, plots the resistance network and pairwise models. The
    #     plotting is dependent on the availability of genetic distance data.

    #     Args:
    #         oname (str): The base name for output files.
    #         mat_r (numpy.ndarray): The average resistance matrix.
    #         edge_r (DataFrame): The average edge resistance values.
    #         plot (bool, optional): Boolean to control whether to generate plots
    #                                Defaults to True.
    #     """

    #     # Write edge and matrix results to files
    #     out1 = f"{oname}.ResistanceEdges.tsv"
    #     utils.write_edges(out1, edge_r, edge_r.index)

    #     out2 = f"{oname}.ResistanceMatrix.tsv"
    #     utils.write_matrix(out2, mat_r, list(self._points_labels.values()))

    #     # Plot resistance network and pairwise models if required
    #     if plot:
    #         edf = pd.DataFrame(
    #             list(
    #                 zip(edge_r.index, edge_r)
    #             ), columns=["EDGE_ID", "Resistance"]
    #         )
    #         self.plot_resistance_network(edf, oname)

    #         # Plot pairwise model if genetic distance data is available
    #         if self._gendist is not None:
    #             self.plot_pairwise_model(mat_r, oname, partition=False)

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
                 edge_order, gendist, adj, origin, R, edge_site_indices):

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