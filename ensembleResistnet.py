import sys
import os, glob
os.environ['USE_PYGEOS'] = '0'
import itertools
import math
import random
import getopt
import pickle
import momepy
import scipy as sp
import numpy as np
import traceback
import pyogrio
import pandas as pd
import geopandas as gpd
import networkx as nx
import matplotlib.pyplot as plt
import multiprocessing as mp
from datetime import datetime
from functools import partial
from collections import OrderedDict
from sortedcontainers import SortedDict
import timeit
from math import radians, degrees, sin, cos, asin, acos, sqrt

from deap import base, creator, tools, algorithms

from resistnet.params import parseArgs
import resistnet.transform as trans
import resistnet.hall_of_fame as hof
from resistnet.hall_of_fame import hallOfFame
import resistnet.hall_of_fame as hof
import resistnet.resist_dist as rd
import resistnet.aggregators as agg
import resistnet.stream_plots as splt

def main():

    params = parseArgs()

    #seed random number generator
    #random.seed(params.seed)
    if not params.seed:
        params.seed=datetime.now().timestamp()
    random.seed(params.seed)

    #########################################################
    # Step 1: Read inputs 
    #########################################################

    # check if samples should be split for each input 
    index_col=None
    if params.split_samples:
        index_col="Sample"
    
    # check if paths provided as list or directory 
    params.paths = None
    if params.input:
        params.paths = params.input
    elif params.input_list:
        params.paths = params.input_list
    
    # get graph from shapefile
    G = read_network(params.network, params.shapefile)
    if params.minimize:
        params.network = str(params.out) + ".minimalSubgraph.net"
    else:
        params.network = str(params.out) + ".subgraph.net"

    # read coordinates 
    if params.coords is None:
        coords = read_and_concat_files(params.paths, ".coords", index_column=index_col)
    else:
        coords = pd.read_csv(params.coords, sep="\t", header=0)
        coords.columns = ["Sample", "Lat", "Lon"]

    coords = coords.drop_duplicates(subset=["Lat", "Lon"], keep='first')

    coords.to_csv(params.out+".coords", 
        header=True, 
        index=False, 
        quoting=None,
        sep="\t")

    # read models 
    if params.only_best:
        local_rows=1
    else:
        local_rows=None
    hofs = read_and_concat_files(params.paths, ".HallOfFame.tsv", local_rows=local_rows)
    if (hofs['fitness'] == hofs['aic']).all():
        hofs["fitness"] = hofs["fitness"]*-1
    hofs["aic"] = hofs["aic"]*-1
    hofs.sort_values(by='fitness', ascending=False, inplace=True, ignore_index=True)
    if params.only_keep:
        hofs = hofs[hofs['keep'] == True]
    if params.hof_max is not None:
        hofs = hofs.head(params.hof_max)
    hofs["keep"] = False
    bests = hallOfFame.from_dataframe(hofs)

    # # read point coordinate data
    # points = pd.read_csv(params.coords, sep="\t", header=0)
    point_coords = snapPoints(coords, G, params.out)

    # compute subgraph (minimized if params.minimze==True)
    K, Kmin = parseSubgraphFromPoints(params, point_coords, G, out=params.out)

    annotateEdges(K, Kmin, params.reachid_col, params.out)

    #########################################################
    # Step 2: Parse models
    #########################################################
    if params.recalc_aic:
        bests.recalculate_aic()
    bests.delta_aic()
    bests.akaike_weights()
    bests.cumulative_akaike(threshold=params.awsum)
    bests.relative_variable_importance(params.only_keep)
    bests.model_average_weights()
    bests.get_best_model()
    #bests.printHOF()
    bests.printRVI()
    bests.printMAW()
    bests.printBest()
    bests.writeMAW(params.out)
    bests.writeRVI(params.out)
    bests.writeBest(params.out)
    bests.writeModelSummary(params.out)

    

    #########################################################
    # Step 3: Initialize worker processes
    #########################################################

    # set up variables and optional aggregator functions
    params.variables = bests.get_variables()
    if params.varFile is not None:
        if params.variables is not None:
            print("Warning: Variables were specified with both <-v> and <-V>... Over-riding options using file provided with <-V>")
        if params.edge_agg is None:
            params.edge_agg="ARITH"
        with open(params.varFile) as fp:
            for line in fp:
                line=line.strip()
                stuff=line.split("\t")
                if len(stuff) < 2:
                    params.agg_opts[stuff[0]]=params.edge_agg
                else:
                    params.agg_opts[stuff[0]]=stuff[1]
        params.variables=list(params.agg_opts.keys())
    else:
        for v in params.variables:
            params.agg_opts[v]=params.edge_agg

    # establish process pool 
    pool = mp.Pool(processes=params.GA_procs)

    process_list = range(1, int(params.GA_procs)+1)
    func = partial(initialize_worker, params)
    results = pool.map(func, process_list)

    #load data for master process
    global my_number
    my_number = 0

    print("Loading data...", flush=True)
    load_data(params, 0)

    #########################################################
    # Step 4: Model average resistances 
    #########################################################

    modelAverage(pool, bests.getHOF(only_keep=True), base=params.out) #set to true for production


def read_and_concat_files(paths, extension, index_column=None, local_rows=None):
    # Ensure paths is a list, even if a single string is provided
    if isinstance(paths, str):
        paths = [paths]
    
    # Initialize an empty list to store the DataFrames
    dfs = []

    # Loop through all paths in the input list
    i=0
    for path in paths:
        # Use glob to find all files with the specified extension
        if os.path.isdir(path):
            pattern = os.path.join(path, f"*{extension}")
        else:
            pattern = f"{path}*{extension}"
        for file_path in glob.glob(pattern):
            # Read the file into a DataFrame 
            df = pd.read_csv(file_path, sep="\t", header=0)

            # Add file index to the specified column if index_column is provided
            if index_column is not None and index_column in df.columns:
                df[index_column] = df[index_column].astype(str) + f"_{i}"
            
            # if local_rows set, only take top X rows from each data frame 
            if local_rows is not None:
                df = df.head(local_rows)

            # append it to list 
            dfs.append(df)

            i+=1

    # Concatenate all the DataFrames in the list
    merged_df = pd.concat(dfs, ignore_index=True)
    merged_df = merged_df.drop_duplicates()

    return merged_df



def modelOutput(stuff):
    model_num = stuff[0]
    individual = stuff[1]
    first=True
    for i, variable in enumerate(predictors.columns):
        if individual[variable] == 1:

            #print("Before:", predictors[variable])
            var = transform(predictors[variable], individual[variable+"_trans"], individual[variable+"_shape"])
            #print("After:", var)
            if first:
                #transform(data, transformation, shape) * weight
                multi = var*(individual[variable+"_weight"])
                first=False
            else:
                multi += var*(individual[variable+"_weight"])
            multi = trans.rescaleCols(multi, 1, 10)

    # OUTPUT PREFIX
    oname=str(params.out)+".Model-"+str(model_num)

    # get point names
    names=list(points_names.values())

    # Get pairwise r
    r=rd.effectiveResistanceMatrix(points, inc_matrix, multi)

    if params.report_all:
        writeMatrix(oname, r, names)

        # get edgewise R
        writeEdges(oname, multi, multi.index)

        if params.plot:
            edf=pd.DataFrame(list(zip(multi.index, multi)), columns=["EDGE_ID", "Resistance"])
            if params.minimize:
                id_col="EDGE_ID"
            else:
                id_col=params.reachid_col
            splt.plotEdgesToStreams((str(params.out)+".subgraph.net"), edf, oname, id_col)

    return(model_num, r, multi)

def modelAverage(pool, bests, base=""):
    #build model list and run Circuitscape on each model
    models=list()
    mod_num=0
    weights=dict()
    for index, row in bests.iterrows():
        n=len(predictors.columns)*4
        model_string = row.iloc[1:n+1]
        models.append([mod_num, model_string])
        weights[mod_num]=row["akaike_weight"]
        mod_num+=1

    print("Model-averaging across",mod_num,"resistance models...")

    # evaluate models
    model_results = pool.map(modelOutput, models)

    # model averaging w/ akaike_weights
    edge_avg=None
    matrix_avg=np.zeros(shape=(len(points), len(points)))

    for m in model_results:
        model_num=m[0]
        matrix_r=m[1]
        edge_r=m[2]
        oname=str(base)+".Model-"+str(model_num)

        # get weighted contributions to av model
        weight = bests.iloc[model_num]["akaike_weight"]
        if edge_avg is None:
            edge_avg=(edge_r*weight)
        else:
            edge_avg = np.add(edge_avg, (edge_r*weight))
        matrix_avg += (matrix_r*weight)


    oname=str(base) + ".Model-Average"
    writeEdges(oname, edge_avg, edge_avg.index)
    writeMatrix(oname, matrix_avg, list(points_names.values()))

    if params.plot:
        edf=pd.DataFrame(list(zip(edge_avg.index, edge_avg)), columns=["EDGE_ID", "Resistance"])
        if params.minimize:
            id_col="EDGE_ID"
        else:
            id_col=params.reachid_col
        splt.plotEdgesToStreams((str(params.out)+".subgraph.net"), edf, oname, id_col)


def writeEdges(oname, edge, ids, dist=None):
    out=(str(oname)+".ResistanceEdges.tsv")
    df=pd.DataFrame(list(zip(ids, edge)), columns=["EDGE_ID", "Resistance"])
    df.to_csv(out, sep="\t", header=True, index=False)
    return(out)


def writeMatrix(oname, mat, ids):
    out=(str(oname)+".ResistanceMatrix.tsv")
    df=pd.DataFrame(mat, columns=ids, index=ids)
    df.to_csv(out, sep="\t", header=True, index=True)
    return(out)


# DEPRECATED -- CIRCUITSCAPE VERSION
def initialize_worker(params, proc_num):
    global my_number
    my_number = proc_num
    #make new random seed, as seed+Process_number
    local_seed = int(params.seed)+my_number
    random.seed(local_seed)
    print("Worker",proc_num,"initializing...\n", flush=True)

    #print("")
    load_data(params, my_number)

    return(local_seed)

# read network
def read_network(network, shapefile):
    if network:
        print("Reading network from saved file: ", network)
        G=nx.Graph(nx.read_gpickle(network).to_undirected())
    else:
        print("Building network from shapefile:",shapefile)
        print("WARNING: This can take a while with very large files!")
        #rivers = gpd.read_file(shapefile)
        rivers = pyogrio.read_dataframe(shapefile)
        G=momepy.gdf_to_nx(rivers, approach="primal", directed=False, multigraph=False)
        #G=nx.Graph(nx.read_shp(shapefile, simplify=True, strict=True).to_undirected())
    return(G)

#returns a pandas DataFrame from points dictionary
def getPointTable(points):
    temp=list()
    for p in points:
        temp.append([p, points[p][1], points[p][0]])
    p = pd.DataFrame(temp, columns=['sample', 'lat', 'long'])
    return(p)

#get subgraph from inputs
def parseSubgraphFromPoints(params, point_coords, G, out=None):
    points=point_coords

    #first pass grabs subgraph from master shapefile graph
    print("\nExtracting full subgraph...")
    K=pathSubgraph(G, points, extractFullSubgraph, params.reachid_col, params.length_col)
    if out:
        net_out=str(out) + ".subgraph.net"
        nx.write_gpickle(K, net_out, pickle.HIGHEST_PROTOCOL)
        # pos=dict()
        # for n in K.nodes:
        #     pos[n] = n
        # color_map = []
        # for node in K:
        #     if node in points.values():
        #         color_map.append("blue")
        #     else:
        #         color_map.append("black")
        # #draw networkx
        # nx.draw_networkx(K, pos, with_labels=False, node_color=color_map, node_size=50)
        # nx.draw_networkx_edge_labels(K, pos, font_size=6)
        # network_plot=str(out) + ".subgraph.pdf"
        # plt.savefig(network_plot)
    del G

    #second pass to simplify subgraph and collapse redundant nodes
    print("\nMerging redundant paths...\n")
    Kmin=pathSubgraph(K, points, extractMinimalSubgraph, params.reachid_col, params.length_col)
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

    return(K, Kmin)

def annotateEdges(K, Kmin, reachid_col, out=None):
    reach_to_edge = dict()
    i=0
    edges = list()
    #print("K:",len(K.edges())
    for p1, p2, dat in Kmin.edges(data=True):
        reaches = dat[reachid_col]
        nx.set_edge_attributes(Kmin, {(p1, p2): {"EDGE_ID": int(i)}})
        for r in reaches:
            reach_to_edge[r] = str(i)
        i+=1

    for p1, p2, dat in K.edges(data=True):
        #print(dat)
        edge_id = reach_to_edge[dat[reachid_col]]
        nx.set_edge_attributes(K, {(p1, p2): {"EDGE_ID": int(edge_id)}})

    if out:
        kmin_out=str(out) + ".minimalSubgraph.net"
        nx.write_gpickle(Kmin, kmin_out, pickle.HIGHEST_PROTOCOL)
        k_out=str(out) + ".subgraph.net"
        nx.write_gpickle(K, k_out, pickle.HIGHEST_PROTOCOL)

    return(K, Kmin)

def load_data(p, proc_num):

    #make "local" globals (i.e. global w.r.t each process)
    global graph
    global distances
    global predictors
    global inc_matrix
    global points
    global gendist
    global edge_ids
    global params
    global id_col
    global points_names
    params = p
    my_number = proc_num
    distances = None

    #read autoStreamTree outputs
    if params.network:
        if my_number == 0:
            print("Reading network from: ", str(params.network))
        graph = readNetwork(params.network)
    else:
        sys.exit("ERROR: No network specified. Nothing to load.")

    if params.minimize:
        full_subgraph = params.network.replace("minimalSubgraph", "subgraph")
        graph = readNetwork(full_subgraph)
        df = nx_to_df(graph)
        id_col="EDGE_ID"
    else:
        df = nx_to_df(graph)
        id_col=params.reachid_col
        df["EDGE_ID"] = df[id_col]


    agg_funs=dict()
    grouped=df.groupby('EDGE_ID')

    for v in params.variables:
        agg_funs[v] = partial(agg.aggregateDist, method=params.agg_opts[v])
    df = grouped.agg(agg_funs)

    names=df.index

    predictors=trans.rescaleCols(df[params.variables], 0, 10)

    if params.dist_col:
        agg_fun = partial(agg.aggregateDist, method=params.efit_agg)
        distances = grouped.agg(agg_fun)[params.dist_col]
        distances=distances.loc[names]
        distances=pd.DataFrame(list(zip(distances.index, distances)), columns=["EDGE_ID", "Edgewise Genetic Distance"])

    # make sure df is sorted the same as names
    predictors = predictors.loc[names]

    # check if samples should be split for each input 
    index_col=None
    if params.split_samples:
        index_col="Sample"

    # read coordinates 
    if params.coords is None:
        cdf = read_and_concat_files(params.paths, ".coords", index_column=index_col)
    else:
        cdf = pd.read_csv(params.coords, sep="\t", header=0)
        cdf.columns = ["Sample", "Lat", "Lon"]

    cdf = cdf.drop_duplicates(subset=["Lat", "Lon"], keep='first')
    cdf['Lat'] = cdf['Lat'].astype(float)
    cdf['Lon'] = cdf['Lon'].astype(float)
    points = SortedDict()
    for index, row in cdf.iterrows():
        name=row['Sample']
        coords=(row['Lon'], row['Lat'])
        if coords not in points:
            points[coords]=list()
        points[coords].append(name)

    # make sure points are snapped to the network
    snapped=SortedDict()
    for point in points.keys():
        if point not in graph.nodes():
            node=snapToNode(graph, point)
            if my_number == 0:
                print("Point not found in graph, snapping to nearest node:", point, " -- ", node)
            if tuple(node) not in snapped:
                snapped[tuple(node)]=list()
            snapped[tuple(node)]=points[point]
        else:
            if tuple(point) not in snapped:
                snapped[tuple(point)]=list()
            snapped[tuple(point)]=points[point]
    points_names = snapped
    index=0
    points=snapped.copy()
    for p in points.keys():
        points[p]=index
        index+=1

    #node_point_dict=nodes_to_points(graph, points)

    # testing: generate pairwise distance matrix for specific variable directly from graph
    # for ia, ib in itertools.combinations(range(0,len(points)),2):
    # ref = getPairwisePathweights(graph, points, params.variables)
    # if my_number == 0:
    #     ofh=params.out+".rawResistMat.txt"
    #     with np.printoptions(precision=0, suppress=True):
    #         np.savetxt(ofh, ref, delimiter="\t")

    # build incidence matrix
    inc_matrix = incidenceMatrix(graph, points, params.length_col, id_col, edge_order=names)
    if my_number == 0:
        ofh=params.out+".incidenceMatrix.txt"
        with np.printoptions(precision=0, suppress=True):
            np.savetxt(ofh, inc_matrix, delimiter="\t", fmt='%i')

    # #read genetic distances
    # gendist = parseInputGenMat(graph, points_names, prefix=params.prefix, inmat=params.inmat, agg_method=params.pop_agg)
    # if my_number == 0:
    #     ofh=params.out+".genDistMat.txt"
    #     with np.printoptions(precision=0, suppress=True):
    #         np.savetxt(ofh, gendist, delimiter="\t")

    # pick first site for each in points_names as final name
    ol=list()
    for p in points_names:
        name = points_names[p][0]
        names=",".join(list(map(str,points_names[p])))
        points_names[p] = name
        ol.append([name, p, points[p], names])

    # # write table mapping final label to input pop labels and node coordinates
    # if my_number == 0:
    #     df = pd.DataFrame(ol, columns = ['Label', 'Index', 'Node', 'Names'])
    #     dtout = str(params.out) + ".pointsTable.txt"
    #     df.to_csv(dtout, sep="\t", index=False)

    #print(points_names)

def getPairwisePathweights(graph, points, attributes):

    def path_weight(G, p1, p2, attributes):
        path=nx.dijkstra_path(G, p1, p2)
        sum=0
        for att in attributes:
            sum += nx.path_weight(G, path, att)
        return(sum)

    pw=pd.DataFrame(columns=list(points.values()), index=list(points.values()))

    for ia, ib in itertools.combinations(range(0,len(points)),2):
        p1=points.keys()[ia]
        p2=points.keys()[ib]
        s=path_weight(graph, p1, p2, attributes)
        pw.loc[list(points.values())[ia], list(points.values())[ib]] = s
        pw.loc[list(points.values())[ib], list(points.values())[ia]] = s
    pw=pw.astype('float64').to_numpy()
    np.fill_diagonal(pw, 0.0)
    return(pw)

def incidenceMatrix(graph, points, len_col, id_col, edge_order):
    #make matrix
    inc = np.zeros((nCr(len(points),2), len(edge_order)), dtype=int)
    edge_to_index=dict()
    index=0
    for e in edge_order:
        edge_to_index[e]=index
        index+=1

    #function to calculate weights for Dijkstra's shortest path algorithm
    #i just invert the distance, so the shortest distance segments are favored
    def dijkstra_weight(left, right, attributes):
        #print(attributes[len_col])
        return(10000000000-attributes[len_col])

    #print(points)
    index=0
    for ia, ib in itertools.combinations(range(0,len(points)),2):
        path = nx.bidirectional_dijkstra(graph, points.keys()[ia], points.keys()[ib], weight=dijkstra_weight)
        for p1, p2, edge in graph.edges(data=True):
            eid=edge[id_col]
            if find_pair(path[1], p1, p2):
                #print("yes:",edge)
                inc[index, edge_to_index[eid]] = 1
            else:
                #print("no")
                inc[index, edge_to_index[eid]] = 0
        index = index+1
        #print("\n---\n")
    return(inc)

#utility function to test if two elements are consecutive in list (irrespective of order)
def find_pair(list, x, y):
    if x not in list or y not in list:
        return(False)
    elif abs(list.index(x)-list.index(y)) == 1:
        return(True)
    else:
        return(False)

def nx_to_df(G):
    l=list()
    for p1, p2, e in G.edges(data=True):
        e["left"] = p1
        e["right"] = p2
        l.append(e)
    df=pd.DataFrame(l)
    return(df)

def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

def checkFormatGenMat(mat, points, agg_method="ARITH"):
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
                print("\nSome populations have identical coordinates. Merging genetic distances using method", agg_method)
                
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
            return(None)
    else:
        return(None)

def parseInputGenMat(graph, points, prefix=None, inmat=None, agg_method="ARITH"):
    #read input matrix
    gen = checkFormatGenMat(inmat, points, agg_method)
    if gen is not None:
        return(gen)
    else:
        print("Failed to read input genetic distance matrix:",inmat)
        sys.exit()

# DEPRECATED
# def nodes_to_points(graph, points):
#     "node to point"
#     np=OrderedDict()
#     seen=dict()
#     node_idx=0
#     for edge in graph.edges:
#         if edge[0] not in seen.keys():
#             if edge[0] in points.keys():
#                 #print("New node:",points[edge[0]], " -- Index:",node_idx)
#                 np[node_idx] = points[edge[0]]
#             seen[edge[0]] = True
#             node_idx+=1
#         if edge[1] not in seen.keys():
#             if edge[1] in points.keys():
#                 #print("New node:",points[edge[1]], " -- Index:",node_idx)
#                 np[node_idx] = points[edge[1]]
#             seen[edge[1]] = True
#             node_idx+=1
#     return(np)


def transform(dat, transformation, shape):
    d=dat
    if transformation <= 0:
        pass
    elif transformation == 1:
        d=trans.ricker(dat, shape, 10)
    elif transformation == 2:
        if params.allShapes:
            d=trans.revRicker(dat, shape, 10)
        else:
            d=trans.ricker(dat, shape, 10)
    elif transformation == 3:
        if params.allShapes:
            d=trans.invRicker(dat, shape, 10)
        else:
            d=trans.revInvRicker(dat, shape, 10)
    elif transformation == 4:
        d=trans.revInvRicker(dat, shape, 10)
    elif transformation == 5:
        d=trans.monomolecular(dat, shape, 10)
    elif transformation == 6:
        if params.allShapes:
            d=trans.revMonomolecular(dat, shape, 10)
        else:
            d=trans.monomolecular(dat, shape, 10)
    elif transformation == 7:
        if params.allShapes:
            d=trans.invMonomolecular(dat, shape, 10)
        else:
            d=trans.revInvMonomolecular(dat, shape, 10)
    elif transformation == 8:
        d=trans.revInvMonomolecular(dat, shape, 10)
    else:
        print("WARNING: Invalid transformation type. Returning un-transformed data.")
    return(trans.rescaleCols(d, 0, 10))


def snapPoints(points, G, out=None):
    point_coords=SortedDict()
    snapDists=dict()
    verb=True
    first=True
    for idx, row in points.iterrows():
        name = row[0]
        data = tuple([row[2], row[1]])
        node = snapToNode(G, tuple([row[2], row[1]]))
        snapDists[row[0]] = great_circle(node[0], node[1], row[2], row[1])
        point_coords[name] = node

    print("Read",str(len(point_coords.keys())),"points.")
    #print(list(point_coords.keys()))
    print()

    if out:
        #plot histogram of snap distances
        dtemp = pd.DataFrame(list(snapDists.items()), columns=['name', 'km'])
        dtout = str(out) + ".snapDistances.txt"
        dtemp.to_csv(dtout, sep="\t", index=False)
    #return everything
    return(point_coords)

#TO DO: There is some redundant use of this and similar functions...
def getNodeOrder(graph, points, as_dict=False, as_index=False, as_list=True):
    node_dict=OrderedDict() #maps node (tuple) to node_index
    point_dict=OrderedDict() #maps points (a subset of nodes) to node_index
    order=list()
    node_idx=0
    #print(type(list(points.keys())[0]))
    for edge in graph.edges():
        left = edge[0]
        right = edge[1]
        #print(type(left))
        #print(type(right))
        if left not in node_dict.keys():
            #print("not in dict")
            if left in points.keys():
                node_dict[left] = node_idx
                order.append(points[left])
                point_dict[left] = points[left]
                node_idx+=1
        if right not in node_dict.keys():
            if right in points.keys():
                node_dict[right] = node_idx
                order.append(points[right])
                point_dict[right] = points[right]
                node_idx+=1
    #if as_index return ordered dict of indices
    if as_index:
        return(node_dict)
    #if as_dict return ordered dict of node names
    if as_dict:
        return(point_dict)
    if as_list:
        return(order)
    #otherwise, return list of NAMES in correct order
    return(order)

def readNetwork(network):
    graph=nx.Graph(nx.read_gpickle(network).to_undirected())
    return(graph)

def path_edge_attributes(graph, path, attribute):
    return [graph[u][v][attribute] for (u,v) in zip(path,path[1:])]

#find and extract paths between points from a graph
def pathSubgraph(graph, nodes, method, id_col, len_col):
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

#extracts full subgraph from nodes
def extractFullSubgraph(subgraph, graph, nodelist, id_col, len_col, path):
    for first, second in zip(path, path[1:]):
        if first not in subgraph:
            subgraph.add_node(first)
        if second not in subgraph:
            subgraph.add_node(second)

        dat=graph.get_edge_data(first, second)
        subgraph.add_edge(first, second, **dat)


#extracts a simplified subgraph from paths
#only keeping terminal and junction nodes
def extractMinimalSubgraph(subgraph, graph, nodelist, id_col, len_col, path):
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


#Input: Tuple of [x,y] coordinates
#output: Closest node to those coordinates
def snapToNode(graph, pos):
    #rint("closest_node call:",pos)
    nodes = np.array(graph.nodes())
    node_pos = np.argmin(np.sum((nodes - pos)**2, axis=1))
    #print(nodes)
    #print("closest to ", pos, "is",tuple(nodes[node_pos]))
    return (tuple(nodes[node_pos]))

#utility function to calculate number of combinations n choose k
def nCr(n,k):
    f = math.factorial
    return f(n) // f(k) // f(n-k)

#function to calculate great circle distances
#returns in units of KILOMETERS
def great_circle(lon1, lat1, lon2, lat2, thresh=0.0000001):
    if (abs(lon1 - lon2)) < thresh and (abs(lat1 - lat2)) < thresh:
        return(0.0)
    else:
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        return 6371 * (
            acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
        )

def readPointCoords(pfile):
    d=SortedDict()
    first=True
    with open(pfile, "r") as pfh:
        for line in pfh:
            if first:
                first=False
                continue
            stuff=line.split()
            name=stuff[0]
            coords=tuple([float(stuff[2]), float(stuff[1])])
            if coords not in d:
                d[coords]=list()
            d[coords].append(name)
    return(d)

#Object to parse command-line arguments
class parseArgs():
    def __init__(self):
        #Define options
        try:
            options, remainder = getopt.getopt(sys.argv[1:], 'ho:i:n:s:r:l:c:m:a:L:t:V:XC:', \
            ["help", "out=", "in=", "network=", "reps=", "shp=",
            "len_col=", "id_col=", "split_samples", "max_keep=", 
            "awsum=", "list=", "threads=", "edge_agg=", "varFile=",
            "allShapes", "report_all", "noPlot", "only_best", 
            "only_keep", "coords="])
        except getopt.GetoptError as err:
            print(err)
            self.display_help("\nExiting because getopt returned non-zero exit status.")
        #Default values for params
        #Input params
        self.input=None
        self.input_list=None
        self.out="sim"
        self.variables = None
        self.edge_agg="ARITH"
        self.agg_opts = dict()
        self.varFile=None
        self.network=None
        self.shapefile=None
        self.hof_max=None
        self.only_keep=False
        self.only_best=False
        self.awsum=0.95
        self.split_samples=False
        self.length_col="LENGTH_KM"
        self.dist_col=None
        self.reachid_col="HYRIV_ID"
        self.recalc_aic=False #not used currently 
        self.paths = None
        self.seed = None
        self.GA_procs = 1
        self.minimize=False
        self.allShapes=False
        self.report_all=False
        self.plot=True
        self.coords=None


        #First pass to see if help menu was called
        for o, a in options:
            if o in ("-h", "-help", "--help"):
                self.display_help("Exiting because help menu was called.")

        #Second pass to set all args.
        for opt, arg_raw in options:
            arg = arg_raw.replace(" ","")
            arg = arg.strip()
            opt = opt.replace("-","")
            #print(opt,arg)
            if opt == "h" or opt == "help":
                continue
            elif opt=="in" or opt=="i":
                self.input=arg
            elif opt=="X" or opt=="noPlot":
                self.plot=False
            elif opt=="list" or opt=="L":
                self.input_list=arg.split(",")
            elif opt=="out" or opt=="o":
                self.out=arg
            elif opt=='t' or opt=='procs':
                self.GA_procs=int(arg)
            elif opt == "network" or opt=="n":
                self.network=arg
            elif opt=='s' or opt=='shp':
                self.shapefile = arg
            elif opt == "awsum" or opt=="a":
                self.awsum=float(arg)
            elif opt == "max_keep" or opt=="m":
                self.hof_max=int(arg)
            elif opt == "split_samples":
                self.split_samples=True
            elif opt == "only_keep":
                self.only_keep=True
            elif opt == "only_best":
                self.only_best=True
            elif opt=="report_all":
                self.report_all=True
            elif opt == "C" or opt == "coords":    
                self.coords = arg
            elif opt == "id_col" or opt=="c":
                self.reachid_col=arg
            elif opt == "l" or opt=="len_col":
                self.length_col=arg
            elif opt=="edge_agg":
                self.edge_agg = arg.upper()
                if self.edge_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN", "SUM", "FIRST", "SD", "VAR", "CV"]:
                    self.display_help("Invalid option "+str(arg).upper()+" for option <--edge_agg>")
            elif opt=="V" or opt=="varFile":
                self.varFile=arg
            elif opt=="allShapes":
                self.allShapes=True
            else:
                assert False, "Unhandled option %r"%opt

        #Check manditory options are set
        if not self.input and not self.input_list:
            self.display_help("No input table provided.")
        if not self.network and not self.shapefile:
            self.display_help("No network provided.")

    def display_help(self, message=None):
        if message is not None:
            print()
            print (message)
        print ("\nensembleResistnet.py\n")
        print("Author: Tyler Chafin")
        print ("Description: Utility script for model-averaging across a collection of RestistNet outputs")
        print("""
Arguments:
-s,--shp    : Path to shapefile containing cleaned, contiguous stream reaches
-n,--network    : Input network (pickle'd networkx output; should include all sites)
-i,--in        : Directory containing resistnet outputs
-L,--list    : Alternatively, provide a comma-separated list of model output prexifes (without spaces)
-C,--coords  : Static coordinates, if used .coords files will not be read using the -L/--list prefixes

Optional Arguments:
-t,--procs    : Number of parallel processors
-a,--awsum    : Cumulative Akaike weight threshold to retain top N models [default=0.95]
--only_keep    : Only retain models where column "keep"=True
--only_best    : Only retain best model from each input
--split_samples    : Do not check for overlap in sample names, instead treating all as unique
-X,--noPlot    : Turn off plotting
-m,--max_keep    : Maximum models to keep (default = all models)
-l,--len_col    : Attribute in network corresponding to edge length (def=LENGTH_KM)
-c,--id_col    : Attribute in network corresponding to edge ID (def=EDGE_ID)
-o,--out    : Output file prefix (default=ensemble)
--report_all    : Plot per-stream resistance and generate full outputs for all retained models
--allShapes    : Allow inverse and reverse transformations
--edge_agg    : Method to use when combining variable values across segments (e.g., with --minimize)
    This can take the following options:
        ARITH        : [default] Use arithmetic mean
        MEDIAN    : Use median distance
        HARM        : Use harmonic mean
        ADJHARM    : Adjusted harmonic mean (see docs)
        GEOM        : Use geometric mean
        MIN        : Use minimum distance
        MAX        : Use maximum distance
        FIRST        : Use first value
        SD        : Use standard deviation
        VAR        : Use variance
        SUM        : Use simple sum
        CV        : Use coefficient of variation (=SD/MEAN)
-V,--varfile    : Optional file with variables provided like so:
          var1 \t <Optional aggregator function>
          var2 \t <Optional aggregator function>
          ...
          ...
""")
        print()
        sys.exit()

#Call main function
if __name__ == '__main__':
    main()
