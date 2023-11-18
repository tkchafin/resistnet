import os 
import sys 

os.environ['USE_PYGEOS'] = '0'

from multiprocessing import Process, Queue
from functools import partial
import random
from datetime import datetime
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
    def __init__(self, 
                 network, 
                 shapefile, 
                 coords, 
                 variables,
                 inmat,
                 agg_opts,
                 pop_agg="ARITH",
                 reachid_col="HYRIV_ID", 
                 length_col="LENGTH_KM", 
                 out="out", 
                 verbose=True,
                 minimize=False):
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

        # Call the initialize_network method
        self.initialize_network()

        # initialize predictors
        self.initialize_predictors()

        # compute incidence matrix 
        self.build_incidence_matrix(self.reachid_col, out=self.output_prefix)

        # parse input genetic distance matrix
        if self.inmat is not None:
            self.parse_input_gen_mat(out=self.output_prefix)


    def initialize_predictors(self):

        # read attributes from graph
        df = utils.nx_to_df(self._K)
        if self.minimize:
            id_col="EDGE_ID"
        else:
            id_col=self.reachid_col
            df["EDGE_ID"] = df[id_col]

        # aggregate chosen vars by edge
        agg_funs=dict()
        grouped=df.groupby('EDGE_ID')
        for v in self.variables:
            agg_funs[v] = partial(agg.aggregateDist, method=self.agg_opts[v])
        df = grouped.agg(agg_funs)
        
        # save edge order
        self._edge_order = df.index

        # rescale vars 
        predictors=trans.rescaleCols(df[self.variables], 0, 10)

        # DEPRECATED FUNCTIONALITY -- WILL REVISIT IN LATER VERSIONS
        # if self.dist_col:
        #     agg_fun = partial(agg.aggregateDist, method=self.efit_agg)
        #     distances = grouped.agg(agg_fun)[self.dist_col]
        #     distances=distances.loc[names]
        #     distances=pd.DataFrame(list(zip(distances.index, distances)), columns=["EDGE_ID", "Edgewise Genetic Distance"])

        # make sure df is sorted the same as names
        predictors = predictors.loc[self._edge_order]
        self._predictors, self.variables = utils.scrub_bad_columns(predictors, self.variables, verbose=True)
        #print(predictors)


    def initialize_network(self):
        # get graph from shapefile
        self.read_network()

        # read and snap point coordinate data 
        points = utils.read_points_table(self.coords)

        # snap points to graph 
        self.snap_points(points, self._G, self.output_prefix)

        # compute subgraph (minimized if needed)
        self._K = self.parse_subgraph_from_points(self.output_prefix)
        self._K = self.annotate_edges(self.output_prefix)

        # remove graph copies no longer needed to save mem
        del self._G
        del self.minimized_subgraph 
        del self.subgraph

        # set path to pickled network file 
        if self.minimize:
            self.network = f"{self.output_prefix}.minimalSubgraph.net"
        else:
            self.network = f"{self.output_prefix}.subgraph.net"


    def read_network(self):
        if self.network:
            if self.verbose:
                print("Reading network from saved file: ", self.network)
            self._G = utils.read_pickle_network(self.network)
        else:
            if self.verbose:
                print("Building network from shapefile:", self.shapefile)
                print("WARNING: This can take a while with very large files!")
            rivers = pyogrio.read_dataframe(self.shapefile)
            self._G = momepy.gdf_to_nx(rivers, approach="primal", directed=False, multigraph=False)


    def snap_points(self, points, graph, out=None):
        self._point_coords = SortedDict()
        snapDists = dict()

        for idx, row in points.iterrows():
            name = row['sample']  # Assuming the first column is the name
            # Using column names for latitude and longitude
            original_point = (row['long'], row['lat'])  # (latitude, longitude)
            node = utils.snap_to_node(graph, original_point)

            # Make sure node is a tuple of (latitude, longitude)
            if not isinstance(node, tuple) or len(node) != 2:
                raise ValueError("Node must be a tuple of (latitude, longitude)")

            if out:
                # Calculate the great circle distance
                snapDists[name] = great_circle((row['lat'], row['long']), (node[1], node[0])).km
            self._point_coords[name] = node

        if self.verbose:
            print("Read", str(len(self._point_coords.keys())), "points.")
            print()

        if out:
            dtemp = pd.DataFrame(list(snapDists.items()), columns=['name', 'km'])
            dtout = str(self.output_prefix) + ".snapDistances.txt"
            dtemp.to_csv(dtout, sep="\t", index=False)
        
        #points = readPointCoords(params.coords)

        # make sure points are snapped to the network
        snapped=SortedDict()
        ptemp = utils.read_points_flatten(self.coords)
        for point in ptemp:
            if point not in graph.nodes():
                node=utils.snap_to_node(graph, point)
                if tuple(node) not in snapped:
                    snapped[tuple(node)]=list()
                snapped[tuple(node)]=ptemp[point]
            else:
                if tuple(point) not in snapped:
                    snapped[tuple(point)]=list()
                snapped[tuple(point)]=ptemp[point]
        self._points_names = snapped
        index=0
        self._points_snapped=snapped.copy()
        for p in self._points_snapped.keys():
            self._points_snapped[p]=index
            index+=1
        
        # write table mapping labels to points 
        # pick first site for each in points_names as final name
        ol=list()
        self._points_labels = self._points_names.copy()
        for p in self._points_names:
            name = self._points_names[p][0]
            names=",".join(self._points_names[p])
            self._points_labels[p] = name
            ol.append([name, p, self._points_snapped[p], names])
        if out:
            df = pd.DataFrame(ol, columns = ['Label', 'Node', 'Index', 'Names'])
            dtout = str(out) + ".pointsTable.txt"
            df.to_csv(dtout, sep="\t", index=False)


    def annotate_edges(self, out=None):
        reach_to_edge = dict()
        i=0
        edges = list()
        #print("K:",len(K.edges())
        for p1, p2, dat in self.minimized_subgraph.edges(data=True):
            reaches = dat[self.reachid_col]
            nx.set_edge_attributes(self.minimized_subgraph, {(p1, p2): {"EDGE_ID": int(i)}})
            for r in reaches:
                reach_to_edge[r] = str(i)
            i+=1

        for p1, p2, dat in self.subgraph.edges(data=True):
            #print(dat)
            edge_id = reach_to_edge[dat[self.reachid_col]]
            nx.set_edge_attributes(self.subgraph, {(p1, p2): {"EDGE_ID": int(edge_id)}})

        if out:
            sub_out = str(out) + ".subgraph.net"
            with open(sub_out, 'wb') as f:
                pickle.dump(self.subgraph, f, pickle.HIGHEST_PROTOCOL)
            min_out = str(out) + ".minimalSubgraph.net"
            with open(min_out, 'wb') as f:
                pickle.dump(self.minimized_subgraph, f, pickle.HIGHEST_PROTOCOL)
        
        if self.minimize:
            return self.minimized_subgraph
        else:
            return self.subgraph


    def evaluate_null_model(self, oname=None):
        """
        Evaluate the null model by constructing a resistance vector using the self.dist_col
        and aggregating using SUM. It then uses parsePairwise to evaluate the model and 
        extracts relevant metrics.
        """
        # Extract and aggregate self.dist_col
        df = utils.nx_to_df(self._K)
        if self.minimize:
            id_col = "EDGE_ID"
        else:
            id_col = self.reachid_col
            df["EDGE_ID"] = df[id_col]

        # Aggregate using SUM for self.dist_col
        grouped = df.groupby('EDGE_ID')
        agg_fun = lambda x: x.sum()
        aggregated_dist = grouped.agg({self.length_col: agg_fun})[self.length_col]

        # Sort the aggregated_dist in the same order as self._edge_order 
        aggregated_dist = aggregated_dist.reindex(self._edge_order)

        # Ensure the aggregated_dist is in the correct format
        if not isinstance(aggregated_dist, pd.Series):
            print("Aggregation of length_col did not result in a proper Series.")
            return None

        # Evaluate model using parsePairwise
        _, res = rd.parsePairwise(self._points_snapped, self._inc, aggregated_dist, self._gendist)

        # Extract relevant metrics
        metrics = res.copy()
        metrics["aic"] = -1 * metrics["aic"]
        metrics["aic_null"] = -1 * metrics["aic_null"]
        first_row = metrics.iloc[0]
        formatted_table = pd.DataFrame({
            "loglik": [first_row["loglik"], first_row["loglik_null"]],
            "aic": [first_row["aic"], first_row["aic_null"]],
            "r2m": [first_row["r2m"], np.nan],  
            "delta_aic_null": [first_row["delta_aic_null"], np.nan]  
        }, index=["distance_only", "null"])

        if oname:
            formatted_table.to_csv(f"{oname}.Null-Model.tsv", sep="\t")

        return formatted_table


    def output_and_plot_model(self, oname, mat_r, edge_r, plot=True):
        """
        Outputs model results to files and plots the resistance network and pairwise models.

        :param edge_avg: The average edge resistance values.
        :param matrix_avg: The average resistance matrix.
        :param oname: The base name for output files.
        :param plot: Boolean to control whether to generate plots.
        """
        # Write edge and matrix results to files
        out1 = f"{oname}.ResistanceEdges.tsv"
        utils.write_edges(out1, edge_r, edge_r.index)

        out2 = f"{oname}.ResistanceMatrix.tsv"
        utils.write_matrix(out2, mat_r, list(self._points_labels.values()))

        # Plot resistance network and pairwise models if required
        if plot:
            edf = pd.DataFrame(list(zip(edge_r.index, edge_r)), columns=["EDGE_ID", "Resistance"])
            self.plot_resistance_network(edf, oname)
            
            # Plot pairwise model
            if self._gendist:
                self.plot_pairwise_model(mat_r, oname, partition=False)


    def plot_resistance_network(self, resistance_values, oname):
        """
        Plots the resistance network.

        :param resistance_values: A DataFrame containing resistance values with an 'EDGE_ID' column.
        :param oname: The output name for the saved plot file.
        :param id_col: The column name for edge IDs. Defaults to "EDGE_ID".
        """

        # Convert the graph to a GeoDataFrame
        _, edgeDF = momepy.nx_to_gdf(self._K)

        # Ensure EDGE_ID is of the correct type and merge with resistance values
        edgeDF['EDGE_ID'] = edgeDF[self.reachid_col].astype(int)
        geoDF = edgeDF.merge(resistance_values, on="EDGE_ID")

        # Plotting
        sns.set(style="ticks")
        geoDF.plot(column="Resistance", cmap="RdYlGn_r", legend=True)
        plt.title("Stream network colored by resistance")
        plt.savefig(f"{oname}.streamsByResistance.pdf")
        plt.clf()
        plt.close()


    def plot_pairwise_model(self, resistance_matrix, oname, partition=False):
        """
        Plots the pairwise model of genetic distance against resistance distance.

        :param resistance_matrix: A matrix containing resistance distances.
        :param oname: The output name for the saved plot file.
        :param partition: Boolean indicating whether to partition the data. Defaults to False.
        """
        g = utils.get_lower_tri(self._gendist) 
        r = utils.get_lower_tri(resistance_matrix)

        df = pd.DataFrame(list(zip(g, r)), columns=["Genetic Distance", "Resistance Distance"])

        sns.set(style="ticks")
        if partition:
            npop = self._gendist.shape[0]
            ID = utils.to_from_(npop) 
            df["grp"] = ID["pop1"]
            sns.lmplot(x="Resistance Distance", y="Genetic Distance", hue="grp", data=df)
        else:
            sns.lmplot(x="Resistance Distance", y="Genetic Distance", data=df)

        plt.title("Pairwise-wise Resistance x Genetic Distance")
        plt.savefig(f"{oname}.PairwiseRegplot.pdf")
        plt.clf()
        plt.close()


    #get subgraph from inputs
    def parse_subgraph_from_points(self, out=None):
        points=self._point_coords

        #first pass grabs subgraph from master shapefile graph
        if self.verbose:
            print("\nExtracting full subgraph...")
        K=self._path_subgraph(self._G, points, self._extract_full_subgraph, self.reachid_col, self.length_col)
        if out:
            net_out=str(out) + ".subgraph.net"
            nx.write_gpickle(K, net_out, pickle.HIGHEST_PROTOCOL)

        #second pass to simplify subgraph and collapse redundant nodes
        if self.verbose:
            print("\nMerging redundant paths...\n")
        Kmin=self._path_subgraph(K, points, self._extract_minimal_subgraph, self.reachid_col, self.length_col)
        if out:
            net_out=str(out) + ".minimalSubgraph.net"
            nx.write_gpickle(Kmin, net_out, pickle.HIGHEST_PROTOCOL)
            pos=dict()
            for n in Kmin.nodes:
                pos[n] = n
            color_map = []
            for node in Kmin:
                if node in points.values():
                    color_map.append("blue")
                else:
                    color_map.append("black")
            #draw networkx
            nx.draw_networkx(Kmin, pos, with_labels=False, node_color=color_map, node_size=50)
            #nx.draw_networkx_edge_labels(Kmin, pos, font_size=6)
            network_plot=str(out) + ".minimalSubgraph.pdf"
            plt.savefig(network_plot)

        self.minimized_subgraph = Kmin
        self.subgraph = K

        if self.minimize:
            return Kmin
        else:
            return K


    def parse_input_gen_mat(self, out=None):
        #read input matrix
        gen = self._check_format_inmat(self.inmat, self._points_names, self.pop_agg)
        if gen is not None:
            self._gendist = gen
        else:
            print("Failed to read input genetic distance matrix:",self.inmat)
            sys.exit()
        

    def build_incidence_matrix(self, edge_id, out=None):
        #make matrix
        inc = np.zeros((utils.nCr(len(self._points_snapped),2), len(self._edge_order)), dtype=int)
        edge_to_index=dict()
        index=0
        for e in self._edge_order:
            edge_to_index[e]=index
            index+=1

        #function to calculate weights for Dijkstra's shortest path algorithm
        #i just invert the distance, so the shortest distance segments are favored
        def dijkstra_weight(left, right, attributes):
            #print(attributes[len_col])
            return(10000000000-attributes[self.length_col])

        #print(points)
        index=0
        for ia, ib in itertools.combinations(range(0,len(self._points_snapped)),2):
            path = nx.bidirectional_dijkstra(self._K, self._points_snapped.keys()[ia], self._points_snapped.keys()[ib], weight=dijkstra_weight)
            for p1, p2, edge in self._K.edges(data=True):
                eid=edge[edge_id]
                if utils.find_pair(path[1], p1, p2):
                    #print("yes:",edge)
                    inc[index, edge_to_index[eid]] = 1
                else:
                    #print("no")
                    inc[index, edge_to_index[eid]] = 0
            index = index+1
            #print("\n---\n")
        self._inc = inc

        if out:
            ofh=out+".incidenceMatrix.txt"
            utils.write_numpy_matrix(inc, ofh)


    @staticmethod
    #find and extract paths between points from a graph
    def _path_subgraph(graph, nodes, method, id_col, len_col):
        k=nx.Graph()

        #function to calculate weights for Dijkstra's shortest path algorithm
        #i just invert the distance, so the shortest distance segments are favored
        def dijkstra_weight(left, right, attributes):
            #print(attributes[len_col])
            return(10000000-attributes[len_col])

        p1 = list(nodes.values())[0]
        for p2 in list(nodes.values())[1:]:
        #for p1, p2 in itertools.combinations(nodes.values(),2):
            try:
                #find shortest path between the two points
                path=nx.bidirectional_dijkstra(graph, p1, p2, weight=dijkstra_weight)

                #traverse the nodes in the path to build a minimal set of edges
                method(k, graph, nodes.values(), id_col ,len_col, path[1])

                if p1 not in k:
                    k.add_node(p1)
                if p2 not in k:
                    k.add_node(p2)
            except Exception as e:
                traceback.print_exc()
                print("Something unexpected happened:",e)
                sys.exit(1)
        return(k)

    @staticmethod
    def _check_format_inmat(mat, points, agg_method):
        order = [item for sublist in list(points.values()) for item in sublist]
        if os.path.isfile(mat):
            #read and see if it has the correct dimensions
            inmat = pd.read_csv(mat, header=0, index_col=0, sep="\t")
            inmat.index = inmat.index.astype(str)
            inmat.columns = inmat.columns.astype(str)
            #if correct dimensions, check if labelled correctly
            if len(inmat.columns) >= len(order):
                formatted = inmat.reindex(index=order, columns=order)
                gen = np.zeros((len(points.keys()),len(points.keys())))
                #establish as nan
                gen[:] = np.nan
                #print(formatted.columns)
                if len(formatted.columns) > len(list(points.keys())):

                    # names=list()
                    for ia,ib in itertools.combinations(range(0,len(list(points.keys()))),2):
                        inds1 = [inmat.columns.get_loc(x) for x in points.values()[ia]]
                        inds2 = [inmat.columns.get_loc(x) for x in points.values()[ib]]
                        results = inmat.iloc[inds1, inds2]
                        results = agg.aggregateDist(results.to_numpy(), agg_method)
                        gen[ia, ib] = gen[ib, ia] = results
                    #     names.append(points.values()[ia][0])
                    #     names.append(points.values()[ib][0])
                    # names=unique(names)
                    # print("NAMES", names)
                    np.fill_diagonal(gen, 0.0)
                    #print(gen)
                    return(gen)
                else:
                    return(formatted.to_numpy())
            else:
                print("ERROR: Input genetic distance matrix has wrong dimensions")
                sys.exit(1)
        else:
            return(None)

    @staticmethod
    #extracts full subgraph from nodes
    def _extract_full_subgraph(subgraph, graph, nodelist, id_col, len_col, path):
        for first, second in zip(path, path[1:]):
            if first not in subgraph:
                subgraph.add_node(first)
            if second not in subgraph:
                subgraph.add_node(second)

            dat=graph.get_edge_data(first, second)
            subgraph.add_edge(first, second, **dat)


    @staticmethod
    #extracts a simplified subgraph from paths
    #only keeping terminal and junction nodes
    def _extract_minimal_subgraph(subgraph, graph, nodelist, id_col, len_col, path):
        curr_edge = {id_col:list(), len_col:0.0}
        curr_start=None
        #print("Path:",path)
        #print("nodelist:",nodelist)
        #for each pair of nodes in path
        for first, second in zip(path, path[1:]):
            #if first is either: 1) a site node; or 2) a junction node: add to new graph
            #if second is either:a site or junction, add edge and node to new graph
            #if not, keep edge attributes for next edge
            if not curr_start:
                curr_start=first
                if first in nodelist or len(graph[first])>2:
                    subgraph.add_node(first)
            #add path attributes to current edge
            dat=graph.get_edge_data(first, second)
            #print(dat)

            curr_edge[id_col].extend([dat[id_col]] if not isinstance(dat[id_col], list) else dat[id_col])
            curr_edge[len_col]=float(curr_edge[len_col])+float(dat[len_col])

            #if second node is a STOP node (=in nodelist or is a junction):
            if second in nodelist or len(graph[second])>2:
                #add node to subgraph
                subgraph.add_node(second)
                #link current attribute data
                subgraph.add_edge(curr_start, second, **curr_edge)
                #empty edge attributes and set current second to curr_start
                curr_edge = {id_col:list(), len_col:0}
                curr_start = second
            else:
                #otherwise continue building current edge
                continue

    def transform(self, dat, transformation, shape):
        d=dat
        if transformation <= 0:
            pass
        elif transformation == 1:
            d=trans.ricker(dat, shape, 10)
        elif transformation == 2:
            if self.allShapes:
                d=trans.revRicker(dat, shape, 10)
            else:
                d=trans.ricker(dat, shape, 10)
        elif transformation == 3:
            if self.allShapes:
                d=trans.invRicker(dat, shape, 10)
            else:
                d=trans.revInvRicker(dat, shape, 10)
        elif transformation == 4:
            d=trans.revInvRicker(dat, shape, 10)
        elif transformation == 5:
            d=trans.monomolecular(dat, shape, 10)
        elif transformation == 6:
            if self.allShapes:
                d=trans.revMonomolecular(dat, shape, 10)
            else:
                d=trans.monomolecular(dat, shape, 10)
        elif transformation == 7:
            if self.allShapes:
                d=trans.invMonomolecular(dat, shape, 10)
            else:
                d=trans.revInvMonomolecular(dat, shape, 10)
        elif transformation == 8:
            d=trans.revInvMonomolecular(dat, shape, 10)
        else:
            print("WARNING: Invalid transformation type. Returning un-transformed data.")
        return(trans.rescaleCols(d, 0, 10))

    #evaluate given a model 
    def evaluate(self, individual):
        #print("evaluate - Process",my_number)
        #vector to hold values across edges
        fitness=0
        res=None
        multi=None
        first=True

        #build multi-surface
        for i, variable in enumerate(self._predictors.columns):
            #Perform variable transformations (if desired)
            #1)Scale to 0-10; 2) Perform desired transformation; 3) Re-scale to 0-10
            #add weighted variable data to multi
            if individual[0::4][i] == 1:
                #print("Before:", predictors[variable])
                var = self.transform(self._predictors[variable], individual[2::4][i], individual[3::4][i])
                #print("After:", var)
                if first:
                    #transform(data, transformation, shape) * weight
                    multi = var*(individual[1::4][i])
                    first=False
                else:
                    multi += var*(individual[1::4][i])

        #If no layers are selected, return a zero fitness
        if first:
            fitness = float('-inf')
        else:
            #Rescale multi for circuitscape
            multi = trans.rescaleCols(multi, 1, 10)

            # evaluate using simple resistance distance
            r, res = rd.parsePairwise(self._points_snapped, self._inc, multi, self._gendist)
            fitness = res[self.fitmetric][0]
            res=list(res.iloc[0])

        #return fitness value
        return(fitness, res)

    def model_output(self, model):
        first=True
        for i, variable in enumerate(self._predictors.columns):
            if model[0::4][i] == 1:
                #print("Before:", predictors[variable])
                var = self.transform(self._predictors[variable], model[2::4][i], model[3::4][i])
                #print("Before:", var)
                if first:
                    #transform(data, transformation, shape) * weight
                    multi = var*(model[1::4][i])
                    first=False
                else:
                    multi += var*(model[1::4][i])
                multi = trans.rescaleCols(multi, 1, 10)

        # get point names
        names=list(self._points_names.values())

        # Get pairwise r
        r=rd.effectiveResistanceMatrix(self._points_snapped, self._inc, multi)

        return(r, multi)

class ResistanceNetworkWorker(ResistanceNetwork):
    def __init__(self, network, pop_agg, reachid_col, length_col, 
                 variables, agg_opts, fitmetric, posWeight, fixWeight, 
                 allShapes, fixShape, min_weight, max_shape, inc, point_coords, 
                 points_names, points_snapped, points_labels, predictors, edge_order, 
                 gendist):
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
        self._point_coords = SortedDict()
        self._points_snapped = SortedDict()
        for i in samples.keys():
            self._point_coords[samples[i]] = i
            self._points_snapped[i] = samples[i]
        self._points_names = self._points_snapped.copy()
        self._points_labels = self._points_snapped.copy()


    def clear_sim(self):
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
        # read attributes from graph
        df = utils.nx_to_df(self._K)
        if self.minimize:
            id_col="EDGE_ID"
        else:
            id_col=self.reachid_col
            df["EDGE_ID"] = df[id_col]

        # save edge order
        self._edge_order = df.index

        # rescale vars 
        predictors=trans.rescaleCols(df[vars['VAR']], 0, 10)

        # generate composite resistance
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
        
        return(multi)

    def write_points(self, oname, points):
        d = {"Sample": [], "Lat": [], "Lon": []}
        for p in points:
            d["Sample"].append(points[p])
            d["Lat"].append(p[1])
            d["Lon"].append(p[0])
        df = pd.DataFrame(d)
        df.to_csv(oname+".coords", 
        header=True, 
        index=False, 
        quoting=None,
        sep="\t")


    def simulate(self, spec_file, num_reps, num_samples, out):

        verbose = self.verbose 
        self.verbose = False
        # read specifications file
        vars = pd.read_csv(spec_file, header=0, delimiter="\t")

        if verbose:
            print(f"Starting simulation with {num_reps} replicates and {num_samples} samples each.")

        for r in range(1, num_reps + 1):
            if verbose:
                print(f"Simulating replicate {r}/{num_reps}")

            base = out + "_" + str(r)

            # sample random points
            samples = utils.sample_nodes(self._G, num_samples)
            self.format_sampled_points(samples)
            
            # compute subgraph 
            self._K = self.parse_subgraph_from_points()
            self._K = self.annotate_edges()

            # simulate composite resistance 
            multi = self.simulate_resistance(vars)

            # build incidence matrix 
            self.build_incidence_matrix(self.reachid_col, out=None)

            # compute resistance matrix 
            res = rd.effectiveResistanceMatrix(self._points_snapped, self._inc, multi)
            res = trans.rescaleCols(res, 0, 1)

            # write outputs
            out_mat = f"{base}.ResistanceMatrix.tsv"
            utils.write_matrix(out_mat, res, list(self._points_labels.values()))
            self.write_points(base, self._points_labels)

            # clear simulation 
            self.clear_sim()

            if verbose:
                print(f"Replicate {r} completed and outputs written.")

        self.verbose = verbose

