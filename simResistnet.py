import sys
import os, glob
import geopandas as gpd
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import getopt
import traceback
import pickle
import math
from math import radians, degrees, sin, cos, asin, acos, sqrt
import momepy
import momepy
import itertools
import random
from functools import partial
from collections import OrderedDict

import resistnet.aggregators as agg
import resistnet.resist_dist as rd
import resistnet.stream_plots as splt
import resistnet.transform as trans

def main():
	params = parseArgs()

	# get graph from shapefile
	G = read_network(params.network, None)

	# read variables
	vars=pd.read_csv(params.input, header=0, delimiter="\t")

	for r in range(1, params.reps+1):

		base = params.out + "_" + str(r)

		# sample some random nodes
		samples = samplePoints(G, params.samples)
		#print(samples)

		# get sub-network
		K = parseSubgraphFromPoints(params, samples, G, out=base)
		K = annotateEdges(K)
		#write_network(K, base)

		# get edge-wise variables
		names, predictors = parsePredictors(K, vars.VAR)
		#print(predictors)

		# apply any requested transformations
		multi=None
		first=True
		for i, variable in enumerate(predictors.columns):
			#Perform variable transformations (if desired)
			#1)Scale to 0-10; 2) Perform desired transformation; 3) Re-scale to 0-10
			#add weighted variable data to multi
			#print("Before:", predictors[variable])

			tr_type=vars[vars.VAR == variable].TRANSFORM.item()
			tr_shape=vars[vars.VAR == variable].SHAPE.item()
			wgt=vars[vars.VAR == variable].WEIGHT.item()
			var = transform(predictors[variable], tr_type, tr_shape)
			#print("After:", var)
			if first:
				#transform(data, transformation, shape) * wgt
				multi = var*wgt
				first=False
			else:
				multi += var*wgt

		# get PW resistance matrix
		inc_matrix = incidenceMatrix(K, samples, params.length_col, params.id_col, edge_order=names)
		res=rd.effectiveResistanceMatrix(samples, inc_matrix, multi)
		res=trans.rescaleCols(res, 0, 1)
		
		# write outputs 
		writePoints(base, samples)
		writeMatrix(base, res, list(samples.values()))

def writePoints(oname, points):
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


def writeMatrix(oname, mat, ids):
	out=(str(oname)+".ResistanceMatrix.tsv")
	df=pd.DataFrame(mat, columns=ids, index=ids)
	df.to_csv(out, sep="\t", header=True, index=True)
	return(out)


def transform(dat, transformation, shape, allShapes=False):
	d=dat
	if transformation <= 0:
		pass
	elif transformation == 1:
		d=trans.ricker(dat, shape, 10)
	elif transformation == 2:
		if allShapes:
			d=trans.revRicker(dat, shape, 10)
		else:
			d=trans.ricker(dat, shape, 10)
	elif transformation == 3:
		if allShapes:
			d=trans.invRicker(dat, shape, 10)
		else:
			d=trans.revInvRicker(dat, shape, 10)
	elif transformation == 4:
		d=trans.revInvRicker(dat, shape, 10)
	elif transformation == 5:
		d=trans.monomolecular(dat, shape, 10)
	elif transformation == 6:
		if allShapes:
			d=trans.revMonomolecular(dat, shape, 10)
		else:
			d=trans.monomolecular(dat, shape, 10)
	elif transformation == 7:
		if allShapes:
			d=trans.invMonomolecular(dat, shape, 10)
		else:
			d=trans.revInvMonomolecular(dat, shape, 10)
	elif transformation == 8:
		d=trans.revInvMonomolecular(dat, shape, 10)
	else:
		print("WARNING: Invalid transformation type. Returning un-transformed data.")
	return(trans.rescaleCols(d, 0, 10))

def nx_to_df(G):
	l=list()
	for p1, p2, e in G.edges(data=True):
		e["left"] = p1
		e["right"] = p2
		l.append(e)
	df=pd.DataFrame(l)
	return(df)

def parsePredictors(K, variables, id_col=None, agg_opts=None):

	df = nx_to_df(K)
	if id_col is not None:
		df["EDGE_ID"] = df[id_col]
	# l = list(df.columns)
	# file = open('all_columns.txt','w')
	# for i in l:
	# 	file.write(str(i)+"\n")
	# file.close()

	agg_funs=dict()
	grouped=df.groupby('EDGE_ID')

	for v in variables:
		if agg_opts is not None:
			agg_funs[v] = partial(agg.aggregateDist, method=agg_opts[v])
		else:
			agg_funs[v] = partial(agg.aggregateDist, method="ARITH")
	df = grouped.agg(agg_funs)

	names=df.index

	predictors=trans.rescaleCols(df[variables], 0, 10)

	# make sure df is sorted the same as names
	predictors = predictors.loc[names]

	return(names, predictors)

# def write_network(K, out):
# 	k_out=str(out) + ".subgraph.net"
# 	nx.write_gpickle(K, k_out, pickle.HIGHEST_PROTOCOL)

def samplePoints(G, samples):
	ret=OrderedDict()
	nodes = random.sample(list(G.nodes), samples)
	i=1
	for s in nodes:
		ret[s] = i
		i+=1
	return(ret)

def annotateEdges(K):
	i=0
	for p1, p2, dat in K.edges(data=True):
		nx.set_edge_attributes(K, {(p1, p2): {"EDGE_ID": int(i)}})
		i+=1
	return(K)

# read network
def read_network(network, shapefile, nx3=True):
	if network:
		print("Reading network from saved file: ", network)
		G=nx.Graph(nx.read_gpickle(network).to_undirected())
	else:
		print("Building network from shapefile:",shapefile)
		print("WARNING: This can take a while with very large files!")
		rivers = gpd.read_file(shapefile)
		G=momepy.gdf_to_nx(rivers, approach="primal", directed=False, multigraph=False)
		#G=nx.Graph(nx.read_shp(shapefile, simplify=True, strict=True).to_undirected())
	return(G)

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
def great_circle(lon1, lat1, lon2, lat2):
	lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
	return 6371 * (
		acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
	)

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
		path = nx.bidirectional_dijkstra(graph, list(points.keys())[ia], list(points.keys())[ib], weight=dijkstra_weight)
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


#get subgraph from inputs
def parseSubgraphFromPoints(params, point_coords, G, out=None):
	points=point_coords

	#first pass grabs subgraph from master shapefile graph
	#print("\nExtracting full subgraph...")
	K=pathSubgraph(G, points, extractFullSubgraph, params.length_col)
	if out:
		net_out=str(out) + ".subgraph.net"
		nx.write_gpickle(K, net_out, pickle.HIGHEST_PROTOCOL)
		pos=dict()
		for n in K.nodes:
			pos[n] = n
		color_map = []
		for node in K:
			if node in points.keys():
				color_map.append("blue")
			else:
				color_map.append("black")
		#draw networkx
		nx.draw_networkx(K, pos, with_labels=False, node_color=color_map, node_size=20)
		#nx.draw_networkx_edge_labels(K, pos, font_size=6)
		network_plot=str(out) + ".subgraph.pdf"
		plt.savefig(network_plot)
	#del G
	return(K)

#extracts full subgraph from nodes
def extractFullSubgraph(subgraph, graph, nodelist, path):
	for first, second in zip(path, path[1:]):
		if first not in subgraph:
			subgraph.add_node(first)
		if second not in subgraph:
			subgraph.add_node(second)

		dat=graph.get_edge_data(first, second)
		subgraph.add_edge(first, second, **dat)

#find and extract paths between points from a graph
def pathSubgraph(graph, nodes, method, len_col):
	k=nx.Graph()

	#function to calculate weights for Dijkstra's shortest path algorithm
	#i just invert the distance, so the shortest distance segments are favored
	def dijkstra_weight(left, right, attributes):
		#print(attributes[len_col])
		return(1000000000-attributes[len_col])

	p1 = list(nodes.keys())[0]
	for p2 in list(nodes.keys())[1:]:
	#for p1, p2 in itertools.combinations(nodes.values(),2):
		try:
			#find shortest path between the two points
			path=nx.bidirectional_dijkstra(graph, p1, p2, weight=dijkstra_weight)

			#traverse the nodes in the path to build a minimal set of edges
			method(k, graph, nodes.values(), path[1])

			if p1 not in k:
				k.add_node(p1)
			if p2 not in k:
				k.add_node(p2)
		except Exception as e:
			traceback.print_exc()
			print("Something unexpected happened:",e)
			sys.exit(1)
	return(k)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'ho:i:n:s:r:l:c:', \
			["help", "out=", "in=", "network=", "reps=", "samples=",
			"len_col=", "id_col="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.input=None
		self.out="sim"
		self.network=None
		self.reps=1
		self.samples=50
		self.length_col="LENGTH_KM"
		self.id_col="EDGE_ID"


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
			elif opt=="out" or opt=="o":
				self.out=arg
			elif opt == "network" or opt=="n":
				self.network=arg
			elif opt == "samples" or opt=="s":
				self.samples=int(arg)
			elif opt == "reps" or opt=="r":
				self.reps=int(arg)
			elif opt == "id_col" or opt=="c":
				self.id_col=arg
			elif opt == "l" or opt=="len_col":
				self.length_col=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.input:
			self.display_help("No input table provided.")
		if not self.network:
			self.display_help("No network provided.")

	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nsimResistnet.py\n")
		print("Author: Tyler Chafin")
		print ("Description: Simulate data on a given network for validating resistnet")
		print("""
Arguments:
-n,--network	: Input network (pickle'd networkx output)
-i,--in		: Table giving information on which variables to use to generate resistnet input
-r,--reps	: Number of replicates
-s,--samples	: Number of random nodes to sample
-l,--len_col	: Attribute in network corresponding to edge length (def=LENGTH_KM)
-c,--id_col	: Attribute in network corresponding to edge ID (def=EDGE_ID)
-o,--out	: Output file name (default=out.fas)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
