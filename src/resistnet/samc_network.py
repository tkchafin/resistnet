
import sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

import resistnet.utils as utils
import resistnet.transform as trans
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

    def __init__(self, network, shapefile, coords, variables, inmat, agg_opts,
                 pop_agg="ARITH", reachid_col="HYRIV_ID",
                 length_col="LENGTH_KM", infer_origin="NEXT_DOWN",
                 origin=None, out="out", verbose=True,
                 minimize=False):
        """
        Constructs all the necessary attributes for the ResistanceNetwork
        object.

        Args:
            network (NetworkX graph): The network graph.
            shapefile (str): Path to the shapefile.
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

        # Additional attributes for internal use
        self.minimized_subgraph = None
        self.subgraph = None
        self._inc = None
        self._G = None
        self._K = None
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

        # Initialization methods
        self.initialize_network()
        self.initialize_predictors()
        self.build_incidence_matrix(self.reachid_col, out=self.output_prefix)
        if self.inmat is not None:
            self.parse_input_gen_mat(out=self.output_prefix)
        
        # direct and detect origin 
        print(self._reach_origin)
        self.polarise_network(self._origin)

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

        # Snap points to the network graph
        self.snap_points(points, self._G, self.output_prefix)

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

    def create_adjacency_matrix(self):
        n = len(self._edge_order)
        adjacency_matrix = np.zeros((n, n), dtype=int)
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

        # TODO: If directional model, return abs(adjacency_matrix)

        self._adj = adjacency_matrix

    def polarise_network(self, origin):

        # Re-orient as a DAG
        self._K = utils.graph_to_dag_converging(self._K, self._origin)
        self.plot_dag_directionality(self.output_prefix)
        print(self._K)

        # make adjacency matrix of edges, ordered as in self._edge_order
        self.create_adjacency_matrix()
        print(self._adj)
        adj_out = self.output_prefix + ".adjacencyMatrix.tsv"
        np.savetxt(adj_out, self._adj, delimiter='\t', fmt='%d')

        # we will work on further steps later
        pass

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
        fitness = 0
        res = None
        multi = np.zeros_like(self._adj, dtype=float)
        adj_i = np.transpose(np.nonzero(self._adj))
        first = True

        # Compute any transformations
        for i, variable in enumerate(self._predictors.columns):
            var_m = np.zeros_like(self._adj, dtype=float)
            # Perform variable transformations if the variable is selected
            # (indicated by a value of 1)
            if individual[0::5][i] == 1:
                var = self.transform(
                    self._predictors[variable],
                    individual[2::5][i],
                    individual[3::5][i]
                )
                # perform directional calculations
                if individual[1::5][i] == 1:
                    pass

                # reScale (minmax)
                var = trans.rescaleCols(var, 0, 1)
                print(var)

                # if directional, convert

                # Compute values for entries with 1
                ones_indices = adj_i[self._adj[adj_i[:, 0], adj_i[:, 1]] == 1]
                for i, j in ones_indices:
                    var_m[i, j] = (var[i] + var[j]) / 2.0

                # Compute values for entries with -1
                minus_ones_indices = adj_i[self._adj[adj_i[:, 0], adj_i[:, 1]] == -1]
                for i, j in minus_ones_indices:
                    var_m[i, j] = var[i] - var[j]

                print(var_m)
                var_m = trans.rescaleCols(var_m, 0, 1)
                print(var_m)

                # sum within multivariate adjacency
                if first:
                    multi = var * individual[1::5][i]
                    first = False
                else:
                    multi += var * individual[1::5][i]

                sys.exit()
            # inverse

        # If no layers are selected, return a zero fitness
        if first:
            fitness = float('-inf')
        else:
            pass
            # Compute P matrix 

            # fit SAMC model and compute cfpt matrix 

            # compute likelihoods with MLPE 

            # multi = trans.rescaleCols(multi, 1, 10)
            # r, res = rd.parsePairwise(
            #     self._points_snapped, self._inc, multi, self._gendist
            # )
            # # fitness = res[self.fitmetric][0]
            # fitness = res[self.fitmetric].iloc[0]
            # res = list(res.iloc[0])

        # Return fitness value and results
        return (fitness, res)

    def plot_inferred_origin(self, oname):
        # Extract the positions of the nodes in the graph
        pos = {n: n for n in self._K.nodes}

        # Create a new figure
        plt.figure(figsize=(8, 6))

        # Draw the edges of the network
        nx.draw_networkx_edges(self._K, pos, alpha=0.5, width=1)

        def check_terminal(p1, p2):
            print(self._K.degree(p1))
            print(self._K.degree(p2))
            if self._K.degree(p1) == 1:
                return p1
            elif self._K.degree(p2) == 1:
                return p2
            else:
                return None

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
    """
    def __init__(self, network, pop_agg, reachid_col, length_col,
                 variables, agg_opts, fitmetric, posWeight, fixWeight,
                 allShapes, fixShape, min_weight, max_shape, inc, point_coords,
                 points_names, points_snapped, points_labels, predictors,
                 edge_order, gendist, adj):

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