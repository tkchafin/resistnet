import os
import sys
from functools import partial
import pickle
import traceback
import itertools
import pandas as pd
import numpy as np
from sortedcontainers import SortedDict
from geopy.distance import great_circle
import seaborn as sns
import matplotlib.pyplot as plt
import pyogrio
import momepy
import networkx as nx

import resistnet.utils as utils
import resistnet.aggregators as agg
import resistnet.transform as trans
import resistnet.resist_dist as rd


class ResistanceNetwork:
    """
    A class to represent a resistance network.

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
                 length_col="LENGTH_KM", out="out", verbose=True,
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

        # Initialization methods
        self.initialize_network()
        self.initialize_predictors()
        self.build_incidence_matrix(self.reachid_col, out=self.output_prefix)
        if self.inmat is not None:
            self.parse_input_gen_mat(out=self.output_prefix)

    def initialize_predictors(self):
        """
        Initializes predictors by aggregating variables for each edge in the
        network graph.

        This method reads attributes from the network graph, aggregates chosen
        variables by edge, and rescales these variables. It also maintains the
        order of edges and handles the selection of appropriate columns for
        further analysis.

        Deprecated functionalities are commented out and may be revisited in
        later versions.
        """

        # Read attributes from the network graph and create a DataFrame
        df = utils.nx_to_df(self._K)

        # Determine the ID column based on whether the network is minimized
        id_col = "EDGE_ID" if self.minimize else self.reachid_col
        df["EDGE_ID"] = df[id_col]

        # Aggregate chosen variables by edge
        agg_funs = {}
        grouped = df.groupby('EDGE_ID')
        for v in self.variables:
            agg_funs[v] = partial(agg.aggregateDist, method=self.agg_opts[v])
        df = grouped.agg(agg_funs)

        # Save the order of edges
        self._edge_order = df.index

        # Rescale variables to a specified range (0 to 10)
        predictors = trans.rescaleCols(df[self.variables], 0, 10)

        # Deprecated functionality for edgewise genetic distance aggregation
        # if self.dist_col:
        #     agg_fun = partial(agg.aggregateDist, method=self.efit_agg)
        #     distances = grouped.agg(agg_fun)[self.dist_col]
        #     distances = distances.loc[names]
        #     distances = pd.DataFrame(
        #           list(zip(distances.index, distances)),
        #           columns=["EDGE_ID", "Edgewise Genetic Distance"]
        #      )

        # Ensure DataFrame is sorted in the same order as edge names
        predictors = predictors.loc[self._edge_order]

        # Clean and select appropriate columns for predictors, updating the
        #     variables list
        self._predictors, self.variables = utils.scrub_bad_columns(
            predictors, self.variables, verbose=True)
        # Print predictors for debugging (optional)
        # print(predictors)

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

    def read_network(self):
        """
        Reads the network from a saved file or constructs it from a shapefile.

        This method checks if a saved network file exists. If it does, the
        network is read from this file. Otherwise, the network is built from a
        provided shapefile.

        """

        # Check if a saved network file exists
        if self.network:
            # Log reading process if verbose mode is enabled
            if self.verbose:
                print("Reading network from saved file: ", self.network)

            # Read the network from the saved file
            self._G = utils.read_pickle_network(self.network)

        else:
            # Log building process if verbose mode is enabled
            if self.verbose:
                print("Building network from shapefile:", self.shapefile)
                print("WARNING: This can take a while with very large files!")

            # Read the shapefile and construct the network
            rivers = pyogrio.read_dataframe(self.shapefile)
            self._G = momepy.gdf_to_nx(
                rivers, approach="primal", directed=False, multigraph=False
            )

    def snap_points(self, points, graph, out=None):
        """
        Snaps points to the nearest node in the network graph and calculates
        distances.

        This method iterates through a set of points, snapping each to the
        closest node in the provided network graph. It optionally calculates
        the great circle distance between the original and snapped points.
        Additionally, it handles the generation of various mappings and output
        files related to these points.

        Args:
            points (DataFrame): DataFrame containing point data (like sample,
                                latitude, longitude).
            graph (NetworkX graph): The network graph to which points will be
                                    snapped.
            out (str, optional): Output prefix for saving distance and points
                                 table. Defaults to None.

        Raises:
            ValueError: If a node is not in the format of a tuple (latitude,
                        longitude).
        """

        self._point_coords = SortedDict()
        snapDists = dict()

        # Iterate through points and snap them to the closest node in the graph
        for idx, row in points.iterrows():
            name = row['sample']  # Assuming 'sample' column as the name
            original_point = (row['long'], row['lat'])  # (longitude, latitude)
            node = utils.snap_to_node(graph, original_point)

            if not isinstance(node, tuple) or len(node) != 2:
                raise ValueError(
                    "Node must be a tuple of (latitude, longitude)"
                )

            # Calculate and store the great circle distance if required
            if out:
                snapDists[name] = great_circle(
                    (row['lat'], row['long']), (node[1], node[0])
                ).km

            self._point_coords[name] = node

        # Log the number of read points if verbose mode is enabled
        if self.verbose:
            print("Read", str(len(self._point_coords.keys())), "points.")
            print()

        # Save snapping distances to a file if output is required
        if out:
            dtemp = pd.DataFrame(
                list(snapDists.items()), columns=['name', 'km']
            )
            dtout = f"{self.output_prefix}.snapDistances.txt"
            dtemp.to_csv(dtout, sep="\t", index=False)

        # make sure points are snapped to the network
        snapped = SortedDict()
        ptemp = utils.read_points_flatten(self.coords)
        for point in ptemp:
            if point not in graph.nodes():
                node = utils.snap_to_node(graph, point)
                if tuple(node) not in snapped:
                    snapped[tuple(node)] = list()
                snapped[tuple(node)] = ptemp[point]
            else:
                if tuple(point) not in snapped:
                    snapped[tuple(point)] = list()
                snapped[tuple(point)] = ptemp[point]
        self._points_names = snapped
        index = 0
        self._points_snapped = snapped.copy()
        for p in self._points_snapped.keys():
            self._points_snapped[p] = index
            index += 1
        
        # write table mapping labels to points 
        # pick first site for each in points_names as final name
        ol = list()
        self._points_labels = self._points_names.copy()
        for p in self._points_names:
            name = self._points_names[p][0]
            names = ",".join(self._points_names[p])
            self._points_labels[p] = name
            ol.append([name, p, self._points_snapped[p], names])

        # Save the table mapping labels to points if output is required
        if out:
            df = pd.DataFrame(ol, columns=['Label', 'Node', 'Index', 'Names'])
            dtout = f"{out}.pointsTable.txt"
            df.to_csv(dtout, sep="\t", index=False)

    def annotate_edges(self, out=None):
        """
        Annotates edges of the subgraph and minimized subgraph with unique
        identifiers.

        This method iterates through the edges of the minimized subgraph and
        the subgraph, assigning each edge a unique identifier. It also creates
        a mapping from reaches to edges. If an output path is provided, it
        saves the annotated subgraphs to specified files.

        Args:
            out (str, optional): Output prefix for saving annotated subgraphs.
                                 Defaults to None.

        Returns:
            NetworkX graph: Returns the minimized subgraph if minimize
                            attribute is True, otherwise returns the subgraph.
        """

        reach_to_edge = dict()
        i = 0

        # Annotate edges in the minimized subgraph
        for p1, p2, dat in self.minimized_subgraph.edges(data=True):
            reaches = dat[self.reachid_col]
            nx.set_edge_attributes(
                self.minimized_subgraph, {(p1, p2): {"EDGE_ID": int(i)}}
            )
            for r in reaches:
                reach_to_edge[r] = str(i)
            i += 1

        # Annotate edges in the subgraph
        for p1, p2, dat in self.subgraph.edges(data=True):
            edge_id = reach_to_edge[dat[self.reachid_col]]
            nx.set_edge_attributes(
                self.subgraph, {(p1, p2): {"EDGE_ID": int(edge_id)}}
            )

        # Save annotated subgraphs if an output path is provided
        if out:
            sub_out = f"{out}.subgraph.net"
            with open(sub_out, 'wb') as f:
                pickle.dump(
                    self.subgraph, f, pickle.HIGHEST_PROTOCOL
                )
            
            min_out = f"{out}.minimalSubgraph.net"
            with open(min_out, 'wb') as f:
                pickle.dump(
                    self.minimized_subgraph, f, pickle.HIGHEST_PROTOCOL
                )

        # Return the appropriate subgraph based on the minimize attribute
        return self.minimized_subgraph if self.minimize else self.subgraph

    def evaluate_null_model(self, oname=None):
        """
        Evaluates the null model by constructing a resistance vector from
        self.dist_col and aggregating it using SUM. It then uses parsePairwise
        to evaluate the model and extracts relevant metrics.

        Args:
            oname (str, optional): The output name for saving the results.
                                   Defaults to None.

        Returns:
            DataFrame: A DataFrame containing the formatted table of metrics
                       including log likelihood, Akaike Information Criterion
                       (AIC), marginal R-squared (r2m), and delta AIC null.

        Notes:
            The method extracts and aggregates self.length_col, evaluates the
            model using parsePairwise, and formats the resulting metrics into a
            table. If oname is provided, it saves the table.
        """

        # Extract and aggregate self.length_col from the network
        df = utils.nx_to_df(self._K)
        id_col = "EDGE_ID" if self.minimize else self.reachid_col
        if not self.minimize:
            df["EDGE_ID"] = df[id_col]

        # Define a function for aggregation using SUM
        def sum_aggregation(x):
            return x.sum()

        # Extract and aggregate self.dist_col from the network
        df = utils.nx_to_df(self._K)
        id_col = "EDGE_ID" if self.minimize else self.reachid_col
        if not self.minimize:
            df["EDGE_ID"] = df[id_col]

        # Aggregate using the defined sum_aggregation function for
        # self.length_col
        grouped = df.groupby('EDGE_ID')
        aggregated_dist = grouped.agg(
            {self.length_col: sum_aggregation}
        )[self.length_col]

        # Check if the aggregation result is in the correct format
        if not isinstance(aggregated_dist, pd.Series):
            print(
                "Aggregation of length_col did not result in a proper Series."
            )
            return None

        # Evaluate the model using parsePairwise and extract relevant metrics
        _, res = rd.parsePairwise(
            self._points_snapped, self._inc, aggregated_dist, self._gendist
        )
        metrics = res.copy()
        metrics["aic"] = -1 * metrics["aic"]
        metrics["aic_null"] = -1 * metrics["aic_null"]
        first_row = metrics.iloc[0]

        # Format the table of metrics
        formatted_table = pd.DataFrame({
            "loglik": [first_row["loglik"], first_row["loglik_null"]],
            "aic": [first_row["aic"], first_row["aic_null"]],
            "r2m": [first_row["r2m"], np.nan],
            "delta_aic_null": [first_row["delta_aic_null"], np.nan]
        }, index=["distance_only", "null"])

        # Save the formatted table if an output name is provided
        if oname:
            formatted_table.to_csv(f"{oname}.Null-Model.tsv", sep="\t")

        return formatted_table

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

        # Write edge and matrix results to files
        out1 = f"{oname}.ResistanceEdges.tsv"
        utils.write_edges(out1, edge_r, edge_r.index)

        out2 = f"{oname}.ResistanceMatrix.tsv"
        utils.write_matrix(out2, mat_r, list(self._points_labels.values()))

        # Plot resistance network and pairwise models if required
        if plot:
            edf = pd.DataFrame(
                list(
                    zip(edge_r.index, edge_r)
                ), columns=["EDGE_ID", "Resistance"]
            )
            self.plot_resistance_network(edf, oname)

            # Plot pairwise model if genetic distance data is available
            if self._gendist is not None:
                self.plot_pairwise_model(mat_r, oname, partition=False)

    def plot_resistance_network(self, resistance_values, oname):
        """
        Plots the resistance network based on given resistance values.

        This method converts the network graph to a GeoDataFrame and merges it
        with the provided resistance values. It then plots the resistance
        network using a color map and saves the plot to a file.

        Args:
            resistance_values (numpy.ndarray): An ndarray containing
                                               resistance values, which should
                                               correspond to the 'EDGE_ID' in
                                               the graph.
            oname (str): The output name for the saved plot file.
        """

        # Convert the graph to a GeoDataFrame
        _, edgeDF = momepy.nx_to_gdf(self._K)

        # Ensure EDGE_ID is of the correct type and merge with resistance
        edgeDF['EDGE_ID'] = edgeDF[self.reachid_col].astype(int)
        geoDF = edgeDF.merge(resistance_values, on="EDGE_ID")

        # Plotting the resistance network
        sns.set(style="ticks")
        geoDF.plot(column="Resistance", cmap="RdYlGn_r", legend=True)
        plt.title("Stream network colored by resistance")
        plt.savefig(f"{oname}.streamsByResistance.pdf")
        plt.clf()
        plt.close()

    def plot_pairwise_model(self, resistance_matrix, oname, partition=False):
        """
        Plots the pairwise model of genetic distance against resistance
        distance.

        This function visualizes the relationship between genetic distance and
        resistance distance in a pairwise manner. It supports partitioning the
        data based on population groups for more detailed analysis.

        Args:
            resistance_matrix (numpy.ndarray): A matrix containing resistance
                                               distances.
            oname (str): The output name for the saved plot file.
            partition (bool, optional): If True, partitions the data based on
                                        population groups. Defaults to False.
        """

        # Extract lower triangle of genetic and resistance distance matrices
        g = utils.get_lower_tri(self._gendist)
        r = utils.get_lower_tri(resistance_matrix)

        # Create a DataFrame for plotting
        df = pd.DataFrame(
            list(zip(g, r)),
            columns=["Genetic Distance", "Resistance Distance"]
        )

        # Set plot style
        sns.set(style="ticks")

        # Plot with or without partitioning based on the 'partition' flag
        if partition:
            npop = self._gendist.shape[0]
            ID = utils.to_from_(npop)
            df["grp"] = ID["pop1"]
            sns.lmplot(
                x="Resistance Distance",
                y="Genetic Distance",
                hue="grp",
                data=df
            )
        else:
            sns.lmplot(
                x="Resistance Distance",
                y="Genetic Distance",
                data=df
            )

        # Set plot title and save the figure
        plt.title("Pairwise-wise Resistance x Genetic Distance")
        plt.savefig(f"{oname}.PairwiseRegplot.pdf")

        # Clear the plot to free up memory
        plt.clf()
        plt.close()

    def parse_subgraph_from_points(self, out=None):
        """
        Extracts and processes subgraphs from the provided points.

        This method performs two main operations: first, it extracts a full
        subgraph from the master shapefile graph using the provided points.
        Second, it simplifies this subgraph by merging redundant paths.
        Optionally, it saves the subgraphs and a plot of the minimal subgraph.

        Args:
            out (str, optional): Output prefix for saving the subgraphs and the
                                 plot. Defaults to None.

        Returns:
            NetworkX graph: The minimal subgraph if 'minimize' is True,
                            otherwise the full subgraph.
        """

        points = self._point_coords

        # First pass: Extract full subgraph from master shapefile graph
        if self.verbose:
            print("\nExtracting full subgraph...")
        K = self._path_subgraph(
            self._G,
            points,
            self._extract_full_subgraph,
            self.reachid_col,
            self.length_col
        )

        # Save the full subgraph if an output path is provided
        if out:
            net_out = f"{out}.subgraph.net"
            nx.write_gpickle(K, net_out, pickle.HIGHEST_PROTOCOL)

        # Second pass: Simplify subgraph by merging redundant paths
        if self.verbose:
            print("\nMerging redundant paths...\n")
        Kmin = self._path_subgraph(
            K,
            points,
            self._extract_minimal_subgraph,
            self.reachid_col,
            self.length_col
        )

        # Save the minimal subgraph and plot if an output path is provided
        if out:
            net_out = f"{out}.minimalSubgraph.net"
            nx.write_gpickle(Kmin, net_out, pickle.HIGHEST_PROTOCOL)

            # Prepare node positions and colors for plotting
            pos = {n: n for n in Kmin.nodes}
            color_map = [
                "blue" if node in points.values() else "black" for node in Kmin
            ]

            # Draw and save the network plot
            nx.draw_networkx(
                Kmin, pos, with_labels=False, node_color=color_map, 
                node_size=50
            )
            network_plot = f"{out}.minimalSubgraph.pdf"
            plt.savefig(network_plot)

        # Update class attributes with the new subgraphs
        self.minimized_subgraph = Kmin
        self.subgraph = K

        # Return the appropriate subgraph based on the 'minimize' flag
        return Kmin if self.minimize else K

    def parse_input_gen_mat(self, out=None):
        """
        Parses the input genetic distance matrix.

        This method reads the input genetic distance matrix and checks its
        format. If the matrix is successfully read and formatted, it is stored
        in the class attribute. If the reading fails, the program exits with an
        error message.

        Args:
            out (str, optional): Output prefix for saving any related data.
                                 Currently not used in the method.
        """

        # Read input matrix and check its format
        gen = self._check_format_inmat(
            self.inmat, self._points_names, self.pop_agg
        )

        # Store the genetic distance matrix or exit on failure
        if gen is not None:
            self._gendist = gen
        else:
            print(
                f"Failed to read input genetic distance matrix: {self.inmat}"
            )
            sys.exit()

    def build_incidence_matrix(self, edge_id, out=None):
        """
        Builds an incidence matrix for the network based on Dijkstra's
        shortest path algorithm.

        This method constructs an incidence matrix where rows represent pairs
        of points and columns represent edges. It uses Dijkstra's algorithm to
        determine the shortest path between pairs of points and marks the edges
        involved in these paths in the incidence matrix.

        Args:
            edge_id (str): The edge ID used to identify edges in the network.
            out (str, optional): Output prefix for saving the incidence matrix.
                                 Defaults to None.
        """

        # Initialize the incidence matrix
        n_combinations = utils.nCr(len(self._points_snapped), 2)
        inc = np.zeros((n_combinations, len(self._edge_order)), dtype=int)

        # Map edge IDs to indices
        edge_to_index = {e: idx for idx, e in enumerate(self._edge_order)}

        # Function to calculate weights for Dijkstra's shortest path algorithm
        def dijkstra_weight(left, right, attributes):
            return attributes[self.length_col]

        # Compute the incidence matrix
        index = 0
        for ia, ib in itertools.combinations(
            range(len(self._points_snapped)), 2
        ):
            path = nx.bidirectional_dijkstra(
                self._K,
                list(self._points_snapped.keys())[ia],
                list(self._points_snapped.keys())[ib],
                weight=dijkstra_weight
            )
            for p1, p2, edge in self._K.edges(data=True):
                eid = edge[edge_id]
                if utils.find_pair(path[1], p1, p2):
                    inc[index, edge_to_index[eid]] = 1
                else:
                    inc[index, edge_to_index[eid]] = 0
            index += 1

        self._inc = inc

        # Save the incidence matrix to a file if an output path is provided
        if out:
            ofh = f"{out}.incidenceMatrix.txt"
            utils.write_numpy_matrix(inc, ofh)

    @staticmethod
    def _path_subgraph(graph, nodes, method, id_col, len_col):
        """
        Finds and extracts paths between points from a graph.

        This method utilizes Dijkstra's shortest path algorithm to find the
        shortest paths between a set of nodes in a graph. It builds a subgraph
        based on these paths using a specified method.

        Args:
            graph (NetworkX graph): The graph from which paths are extracted.
            nodes (dict): A dictionary of nodes for which paths are to be found
            method (function): The method used to process the paths and build
                               the subgraph.
            id_col (str): The column name for edge IDs.
            len_col (str): The column name for edge lengths.

        Returns:
            NetworkX graph: A subgraph containing the paths between the
                            specified nodes.
        """

        k = nx.Graph()

        # Function to calculate weights for Dijkstra's shortest path algorithm
        def dijkstra_weight(left, right, attributes):
            return attributes[len_col]

        # Process each pair of points
        p1 = list(nodes.values())[0]
        for p2 in list(nodes.values())[1:]:
            try:
                # Find shortest path between two points
                path = nx.bidirectional_dijkstra(
                    graph, p1, p2, weight=dijkstra_weight
                )

                # Traverse the nodes in the path to build minimal set of edges
                method(k, graph, nodes.values(), id_col, len_col, path[1])

                # Ensure the nodes are in the new subgraph
                if p1 not in k:
                    k.add_node(p1)
                if p2 not in k:
                    k.add_node(p2)

            except Exception as e:
                traceback.print_exc()
                print(f"Something unexpected happened: {e}")
                sys.exit(1)

        return k

    @staticmethod
    def _check_format_inmat(mat, points, agg_method):
        """
        Checks and formats the input genetic distance matrix.

        This method verifies if the input matrix has the correct dimensions and
        labels. If the matrix is larger than required, it aggregates distances
        based on the provided aggregation method. If the matrix has the correct
        dimensions, it formats it to match the order of points.

        Args:
            mat (str): Path to the input matrix file.
            points (dict): A dictionary mapping points to their labels.
            agg_method (str): The aggregation method to use for larger matrices

        Returns:
            numpy.ndarray: A formatted numpy array representing the genetic
                           distance matrix, or None if the file does not exist
                           or has incorrect dimensions.
        """
        order = [item for sublist in list(points.values()) for item in sublist]
        if os.path.isfile(mat):
            inmat = pd.read_csv(mat, header=0, index_col=0, sep="\t")
            inmat.index = inmat.index.astype(str)
            inmat.columns = inmat.columns.astype(str)

            if len(inmat.columns) >= len(order):
                formatted = inmat.reindex(index=order, columns=order)
                gen = np.zeros((len(points.keys()), len(points.keys())))
                gen[:] = np.nan

                if len(formatted.columns) > len(list(points.keys())):
                    for ia, ib in itertools.combinations(
                        range(len(list(points.keys()))), 2
                    ):
                        point_values_ia = points.values()[ia]
                        point_values_ib = points.values()[ib]

                        inds1 = [
                            inmat.columns.get_loc(x) for x in point_values_ia
                        ]
                        inds2 = [
                            inmat.columns.get_loc(x) for x in point_values_ib
                        ]
                        results = inmat.iloc[inds1, inds2]
                        results = agg.aggregateDist(
                            results.to_numpy(), agg_method
                        )
                        gen[ia, ib] = gen[ib, ia] = results

                    np.fill_diagonal(gen, 0.0)
                    return gen
                else:
                    return formatted.to_numpy()
            else:
                print(
                    "ERROR: Input genetic distance matrix has wrong dimensions"
                )
                sys.exit(1)
        else:
            return None

    @staticmethod
    def _extract_full_subgraph(
        subgraph, graph, nodelist, id_col, len_col, path
    ):
        """
        Extracts the full subgraph from nodes based on the provided path.

        This method traverses the path between nodes and adds the corresponding
        nodes and edges to the subgraph. It ensures that all nodes and edges
        present in the path are included in the subgraph.

        Args:
            subgraph (NetworkX graph): The subgraph to which nodes and edges
                                       will be added.
            graph (NetworkX graph): The original graph from which data will be
                                    extracted.
            nodelist (list): A list of nodes involved in the path.
            id_col (str): The column name for edge IDs. Not used in this method
            len_col (str): The column name for edge lengths. Not used in this
                           method.
            path (list): The path of nodes to be included in the subgraph.
        """

        # Iterate through pairs of nodes in the path
        for first, second in zip(path, path[1:]):
            # Add nodes to the subgraph if they are not already present
            if first not in subgraph:
                subgraph.add_node(first)
            if second not in subgraph:
                subgraph.add_node(second)

            # Retrieve edge data and add the edge to the subgraph
            dat = graph.get_edge_data(first, second)
            subgraph.add_edge(first, second, **dat)

    @staticmethod
    def _extract_minimal_subgraph(
        subgraph, graph, nodelist, id_col, len_col, path
    ):
        """
        Extracts a simplified subgraph from paths, keeping only terminal and
        junction nodes.

        This method processes a path and constructs a minimal subgraph from it.
        The subgraph includes only the terminal nodes (site nodes) and junction
        nodes. It aggregates edge attributes along the path between these nodes

        Args:
            subgraph (NetworkX graph): The subgraph to which nodes and edges
                                       will be added.
            graph (NetworkX graph): The original graph from which data will be
                                    extracted.
            nodelist (list): A list of terminal nodes (site nodes).
            id_col (str): The column name for edge IDs.
            len_col (str): The column name for edge lengths.
            path (list): The path of nodes to be included in the subgraph.
        """

        # Initialize current edge data
        curr_edge = {id_col: list(), len_col: 0.0}
        curr_start = None

        # Iterate through pairs of nodes in the path
        for first, second in zip(path, path[1:]):
            if not curr_start:
                curr_start = first
                if first in nodelist or len(graph[first]) > 2:
                    subgraph.add_node(first)

            # Add path attributes to current edge
            dat = graph.get_edge_data(first, second)

            # Determine the extension for the current edge's id_col
            if not isinstance(dat[id_col], list):
                extension = [dat[id_col]]
            else:
                extension = dat[id_col]
            
            # Extend the current edge's id_col with the determined extension
            curr_edge[id_col].extend(extension)
            curr_edge[len_col] += float(dat[len_col])

            # Check if second node is a STOP node (in nodelist or a junction)
            if second in nodelist or len(graph[second]) > 2:
                subgraph.add_node(second)
                subgraph.add_edge(curr_start, second, **curr_edge)
                curr_edge = {id_col: list(), len_col: 0}
                curr_start = second
            else:
                continue  # Continue building current edge

    def transform(self, dat, transformation, shape):
        """
        Applies a specified transformation to the data.

        This method transforms the data based on the specified transformation
        type and shape. It supports various transformation types, each
        corresponding to a different method.

        Args:
            dat (array-like): The data to be transformed.
            transformation (int): The type of transformation to apply.
            shape (float): A parameter influencing the shape of the
                           transformation.

        Returns:
            array-like: The transformed data, rescaled between 0 and 10.
        """
        d = dat

        # Apply transformation based on the specified type
        if transformation == 1:
            d = trans.ricker(dat, shape, 10)
        elif transformation == 2:
            if self.allShapes:
                d = trans.revRicker(dat, shape, 10)
            else:
                d = trans.ricker(dat, shape, 10)
        elif transformation == 3:
            if self.allShapes:
                d = trans.invRicker(dat, shape, 10)
            else:
                d = trans.revInvRicker(dat, shape, 10)
        elif transformation == 4:
            d = trans.revInvRicker(dat, shape, 10)
        elif transformation == 5:
            d = trans.monomolecular(dat, shape, 10)
        elif transformation == 6:
            if self.allShapes:
                d = trans.revMonomolecular(dat, shape, 10)
            else:
                d = trans.monomolecular(dat, shape, 10)
        elif transformation == 7:
            if self.allShapes:
                d = trans.invMonomolecular(dat, shape, 10)
            else:
                d = trans.revInvMonomolecular(dat, shape, 10)
        elif transformation == 8:
            d = trans.revInvMonomolecular(dat, shape, 10)
        elif transformation <= 0:
            # No transformation applied
            pass
        else:
            print("WARNING: Invalid transformation type.")

        return trans.rescaleCols(d, 0, 10)

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
        multi = None
        first = True

        # Build multi-surface representation of variables
        for i, variable in enumerate(self._predictors.columns):
            # Perform variable transformations if the variable is selected
            # (indicated by a value of 1)
            if individual[0::4][i] == 1:
                var = self.transform(
                    self._predictors[variable],
                    individual[2::4][i],
                    individual[3::4][i]
                )
                if first:
                    multi = var * individual[1::4][i]
                    first = False
                else:
                    multi += var * individual[1::4][i]

        # If no layers are selected, return a zero fitness
        if first:
            fitness = float('-inf')
        else:
            # Rescale multi-surface for circuitscape and evaluate using simple
            # resistance distance
            multi = trans.rescaleCols(multi, 1, 10)
            r, res = rd.parsePairwise(
                self._points_snapped, self._inc, multi, self._gendist
            )
            fitness = res[self.fitmetric][0]
            res = list(res.iloc[0])

        # Return fitness value and results
        return (fitness, res)

    def model_output(self, model):
        """
        Generates the model output for a given model.

        This method combines the transformed variables based on the model
        specifications to produce a multi-surface representation. It then
        calculates the effective resistance matrix for the given model.

        Args:
            model (list): A list representing an individual model in the
                          genetic algorithm, containing transformation and
                          weight information for variables.

        Returns:
            tuple: A tuple containing the effective resistance matrix and the
                   multi-surface representation for the model.
        """
        first = True
        multi = None

        # Combine transformed variables based on model specifications
        for i, variable in enumerate(self._predictors.columns):
            if model[0::4][i] == 1:
                var = self.transform(
                    self._predictors[variable], model[2::4][i], model[3::4][i]
                )
                if first:
                    multi = var * model[1::4][i]
                    first = False
                else:
                    multi += var * model[1::4][i]
                multi = trans.rescaleCols(multi, 1, 10)

        # Get pairwise effective resistance matrix
        r = rd.effectiveResistanceMatrix(
            self._points_snapped, self._inc, multi
        )

        return r, multi


class ResistanceNetworkWorker(ResistanceNetwork):
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
    """
    def __init__(self, network, pop_agg, reachid_col, length_col,
                 variables, agg_opts, fitmetric, posWeight, fixWeight,
                 allShapes, fixShape, min_weight, max_shape, inc, point_coords,
                 points_names, points_snapped, points_labels, predictors,
                 edge_order, gendist):

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
        self._G = None
        self._K = self.read_network()
        self._point_coords = point_coords
        self._points_names = points_names
        self._points_snapped = points_snapped
        self._points_labels = points_labels
        self._predictors = predictors
        self._edge_order = edge_order
        self._gendist = gendist


class SimResistanceNetwork(ResistanceNetwork):
    """
    A simplified subclass of ResistanceNetwork for simulating validation
    datasets.

    Args:
        network (NetworkX graph): The network graph.
        reachid_col (str): Column name for reach ID.
        length_col (str): Column name for length.
        verbose (bool, optional): Flag to enable verbose output. Defaults to
                                  False.
    """
    def __init__(self, network, reachid_col, length_col, verbose=False):
        self.network = network
        self.reachid_col = reachid_col
        self.length_col = length_col
        self.verbose = verbose
        self.minimize = False
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
        self.read_network()

    def format_sampled_points(self, samples):
        """
        Formats and stores sampled points for simulation.

        This method takes a dictionary of sampled points and formats them into
        sorted dictionaries for point coordinates, snapped points, names, and
        labels.

        Args:
            samples (dict): A dictionary of sampled points, where keys are
                            identifiers and values are coordinates.
        """
        self._point_coords = SortedDict()
        self._points_snapped = SortedDict()

        # Store points and their snapped counterparts
        for i in samples.keys():
            self._point_coords[samples[i]] = i
            self._points_snapped[i] = samples[i]

        # Copy snapped points to names and labels
        self._points_names = self._points_snapped.copy()
        self._points_labels = self._points_snapped.copy()

    def clear_sim(self):
        """
        Clears all simulation-related data from the object.

        This method resets all attributes related to simulation.
        """
        self._inc = None
        self._K = None
        self._point_coords = None
        self._points_names = None
        self._points_snapped = None
        self._points_labels = None
        self._predictors = None
        self._edge_order = None
        self._gendist = None

    def simulate_resistance(self, vars):
        """
        Simulates resistance based on various variables and their
        transformations.

        This method reads attributes from the network graph, scales the
        variables, and generates a composite resistance surface by combining
        the transformed variables with their respective weights.

        Args:
            vars (DataFrame): A DataFrame containing variables along with their
                              transformation type, shape, and weight.

        Returns:
            numpy.ndarray: The resulting composite resistance surface.
        """
        # Read attributes from the graph
        df = utils.nx_to_df(self._K)
        id_col = "EDGE_ID" if self.minimize else self.reachid_col
        if not self.minimize:
            df["EDGE_ID"] = df[id_col]

        # Save the order of edges
        self._edge_order = df.index

        # Rescale variables
        predictors = trans.rescaleCols(df[vars['VAR']], 0, 10)

        # Generate composite resistance
        multi = None
        first = True
        for i, variable in enumerate(predictors.columns):
            tr_type = vars[vars['VAR'] == variable]['TRANSFORM'].item()
            tr_shape = vars[vars['VAR'] == variable]['SHAPE'].item()
            wgt = vars[vars['VAR'] == variable]['WEIGHT'].item()
            var = self.transform(predictors[variable], tr_type, tr_shape)
            if first:
                multi = var * wgt
                first = False
            else:
                multi += var * wgt

        return multi

    def write_points(self, oname, points):
        """
        Writes the point data to a file.

        Args:
            oname (str): Output filename prefix.
            points (dict): Dictionary of point data with location and sample
                           information.
        """
        d = {"Sample": [], "Lat": [], "Lon": []}
        for p in points:
            d["Sample"].append(points[p])
            d["Lat"].append(p[1])
            d["Lon"].append(p[0])
        df = pd.DataFrame(d)
        df.to_csv(
            f"{oname}.coords", header=True, index=False, quoting=None, sep="\t"
        )

    def simulate(self, spec_file, num_reps, num_samples, out):
        """
        Simulates resistance based on the specifications file.

        Args:
            spec_file (str): Path to the specifications file.
            num_reps (int): Number of replicates to simulate.
            num_samples (int): Number of samples per replicate.
            out (str): Output filename prefix for the simulation results.
        """
        verbose = self.verbose 
        self.verbose = False

        # Read specifications file
        vars = pd.read_csv(spec_file, header=0, delimiter="\t")

        if verbose:
            print(
                f"Starting simulation with {num_reps} \
                    replicates and {num_samples} samples each."
            )

        for r in range(1, num_reps + 1):
            if verbose:
                print(f"Simulating replicate {r}/{num_reps}")

            base = f"{out}_{r}"

            # Sample random points and format them
            samples = utils.sample_nodes(self._G, num_samples)
            self.format_sampled_points(samples)

            # Compute subgraph and annotate edges
            self._K = self.parse_subgraph_from_points()
            self._K = self.annotate_edges()

            # Simulate composite resistance
            multi = self.simulate_resistance(vars)

            # Build incidence matrix and compute resistance matrix
            self.build_incidence_matrix(self.reachid_col, out=None)
            res = rd.effectiveResistanceMatrix(
                self._points_snapped, self._inc, multi
            )
            res = trans.rescaleCols(res, 0, 1)

            # Write outputs
            out_mat = f"{base}.ResistanceMatrix.tsv"
            utils.write_matrix(
                out_mat, res, list(self._points_labels.values())
            )
            self.write_points(base, self._points_labels)

            # Clear simulation data
            self.clear_sim()

            if verbose:
                print(f"Replicate {r} completed and outputs written.")

        # Restore verbosity setting
        self.verbose = verbose
