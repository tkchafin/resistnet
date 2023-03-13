import sys
import os, glob
import itertools
import math
import random
import getopt
import pickle
import momepy
import scipy as sp
import numpy as np
import traceback
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
	# Step 1: Parsing shapefile
	#########################################################

	# get graph from shapefile
	G = read_network(params.network, params.shapefile)

	# read point coordinate data
	points = pd.read_csv(params.coords, sep="\t", header=0)
	point_coords = snapPoints(points, G, params.out)

	# compute subgraph (minimized if params.minimze==True)
	K, Kmin = parseSubgraphFromPoints(params, point_coords, G, out=params.out)

	annotateEdges(K, Kmin, params.reachid_col, params.out)

	if params.minimize:
		params.network = str(params.out) + ".minimalSubgraph.net"
	else:
		params.network = str(params.out) + ".subgraph.net"

	#########################################################
	# Step 2: Initialize worker processes
	#########################################################

	pool = mp.Pool(processes=params.GA_procs)

	process_list = range(1, int(params.GA_procs)+1)
	func = partial(initialize_worker, params)
	results = pool.map(func, process_list)

	#load data for master process
	global my_number
	my_number = 0

	print("Loading data...", flush=True)
	load_data(params, 0)

	#print(len(list(graph.edges())))
	#print(inc_matrix)
	#sys.exit()
	#print(results)
	#sys.exit()

	#########################################################
	# Step 3: Initialize GA
	#########################################################

	#initialize a single-objective GA
	creator.create("FitnessMax", base.Fitness, weights=(1.0,))
	creator.create("Individual", list, fitness=creator.FitnessMax)

	#toolbox
	global toolbox
	toolbox = base.Toolbox()

	#register GA attributes and type variables
	print("Initializing genetic algorithm parameters...\n")
	initGA(toolbox, params)

	#mp.set_start_method("spawn")

	#initialize population
	popsize=len(params.variables)*4*15
	if params.popsize:
		popsize=params.popsize
	if popsize > params.maxpopsize:
		popsize=params.maxpopsize
	print("Establishing a population of size:",str(popsize))
	pop = toolbox.population(n=popsize)
	#pop = apply_linkage(pop)


	# Evaluate the entire population
	print("\nEvaluating initial population...\n")

	fitnesses = pool.map(toolbox.evaluate, [list(i) for i in pop])

	#print(fitnesses)
	pop_list=list()
	for ind, fit in zip(pop, fitnesses):
		ind.fitness.values = fit[0],
		if fit[1] is not None:
			ind_list=list()
			ind_list.append(fit[0])
			ind_list.extend(list(ind))
			ind_list.extend(fit[1])
			pop_list.append(ind_list)
	#print(pop_list)
	bests = hof.hallOfFame(predictors.columns, params.max_hof_size, pop_list)
	#sys.exit()
	#print(pop[0])

	# CXPB  is the probability with which two individuals are crossed
	# MUTPB is the probability for mutating an individual
	cxpb, mutpb = params.cxpb, params.mutpb

	# Extracting all the fitnesses of population
	fits = [i.fitness.values[0] for i in pop]

	# Variable keeping track of the number of generations
	g = 0

	# Begin the evolution
	#NOTE: Need to implement some sort of callback for

	#########################################################
	# Step 4: Optimization
	#########################################################

	print("Starting optimization...\n")
	fails=0
	current_best=None
	logger=list()
	#while max(fits) < 5 and g < 5:
	while fails <= params.nfail and g < params.maxGens:
		# A new generation
		g = g + 1
		print("-- Generation %i --" % g)

		# Select the next generation individuals
		offspring = toolbox.select(pop, len(pop))
		# Clone the selected individuals
		offspring = list(map(toolbox.clone, offspring))

		# Apply crossover and mutation on the offspring
		for child1, child2 in zip(offspring[::2], offspring[1::2]):
			if random.random() < cxpb:
				toolbox.mate(child1, child2)
				del child1.fitness.values
				del child2.fitness.values

		for mutant in offspring:
			if random.random() < mutpb:
				toolbox.mutate(mutant)
				del mutant.fitness.values

		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]

		#offspring = apply_linkage(offspring)
		#print(offspring)
		#sys.exit()

		fitnesses = pool.map(toolbox.evaluate, [list(i) for i in invalid_ind])
		pop_list=list()
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit[0],
			if fit[1] is not None:
				ind_list=list()
				ind_list.append(fit[0])
				ind_list.extend(list(ind))
				ind_list.extend(fit[1])
				pop_list.append(ind_list)
		bests.check_population(pop_list)

		#replace population with offspring
		pop[:] = offspring

		# Gather all the fitnesses in one list and print the stats
		fits = [i.fitness.values[0] for i in pop]


		#evaluate for stopping criteria
		if g > params.burnin:

			# if params.fitmetric=="aic":
			# 	fits = [element * -1.0 for element in fits]
			# 	worst = max(fits)
			# 	best=min(fits)
			# 	if current_best is None:
			# 		current_best=best
			# 	else:
			# 		(current_best, fails) = updateFails(best, current_best, fails, params.deltaB, params.deltaB_perc, minimize=True)
			# else:
			worst=min(fits)
			best=max(fits)
			if current_best is None:
				current_best=best
			else:
				(current_best, fails) = updateFails(best, current_best, fails, params.deltaB, params.deltaB_perc, minimize=False)

			length = len(pop)
			mean = sum(fits) / length
			sum2 = sum(x*x for x in fits)
			std = abs(sum2 / length - mean**2)**0.5

			print("  Worst %s" % worst)
			print("  Best %s" % best)
			print("  Avg %s" % mean)
			print("  Std %s" % std)
			logger.append([g, worst, best, mean, std])
			print("  nFails %s" % fails)

	#best = pop[np.argmax([pool.map(toolbox.evaluate, [list(ind) for ind in pop])])]
	print("Stopping optimization after",str(g),"generations.")
	if g >= params.maxGens:
		print("Reason: Exceeded maxGens")
	elif fails > params.nfail:
		print("Reason: More than",str(fails),"generations since finding better solution.")

	#########################################################
	# Step 5: Model averaging and calculate variable importance
	#########################################################

	# if params.fitmetric == 'aic':
	# 	bests.correct_aic_fitness()
	bests.delta_aic()
	bests.akaike_weights()
	bests.cumulative_akaike(threshold=params.awsum)
	bests.relative_variable_importance(params.only_keep)
	bests.model_average_weights()
	bests.printHOF()
	bests.printRVI()
	bests.printMAW()
	bests.writeMAW(params.out)
	bests.writeRVI(params.out)
	if params.plot:
		bests.plot_ICprofile(params.out)
		bests.plotMetricPW(params.out)
		bests.plotVariableImportance(params.out)

	#write hall of fame to file
	bests.writeModelSummary(params.out)

	#write log of fitnesses
	logDF=pd.DataFrame(logger, columns=["Generation", "Worst", "Best", "Mean", "Stdev"])
	logDF.to_csv((str(params.out)+".FitnessLog.tsv"), sep="\t", header=True, index=False)
	del logger
	del logDF

	# get results for best models and model-averaged
	modelAverage(pool, bests.getHOF(only_keep=params.only_keep), base=params.out) #set to true for production
	df=pd.DataFrame(gendist, columns=list(points_names.values()), index=list(points_names.values()))
	df.to_csv((str(params.out)+".genDistMat.tsv"), sep="\t", header=True, index=True)

	#clean up temp files
	oname=".temp_"+str(os.path.basename(params.out))+"*"
	for filename in glob.glob(oname):
		os.remove(filename)

	#print("closing pool")
	#pool.close()
	pool.terminate()
	sys.exit(0)

def apply_linkage(pop):
	ret = list()
	for ind in pop:
		for i, x in enumerate(ind[0::4]):
			for j in [1,2,3]:
				ind[(i*4)+j] = ind[(i*4)+j]*x
		ret.append(ind)
	return(ret)

def updateFails(best, current_best, fails, deltB, deltB_perc, minimize=False):
	cur=current_best
	f=fails
	threshold_1=current_best
	threshold_2=current_best

	if minimize is True:
		if best < current_best:
			cur = best
		if params.deltaB is not None:
			threshold_1=current_best-params.deltaB
		if params.deltaB_perc is not None:
			threshold_2=current_best-(current_best*params.deltaB_perc)
		if best >= threshold_1 or best >= threshold_2:
			f += 1
		else:
			f = 0
	else:
		if best > current_best:
			cur = best
		if params.deltaB is not None:
			threshold_1=current_best+params.deltaB
		if params.deltaB_perc is not None:
			threshold_2=current_best+(current_best*params.deltaB_perc)
		if best <= threshold_1 or best <= threshold_2:
			f += 1
		else:
			f = 0
	return(cur, f)

def modelOutput(stuff):
	model_num = stuff[0]
	individual = stuff[1]
	first=True
	for i, variable in enumerate(predictors.columns):
		if individual[0::4][i] == 1:
			#print("Before:", predictors[variable])
			var = transform(predictors[variable], individual[2::4][i], individual[3::4][i])
			#print("Before:", var)
			if first:
				#transform(data, transformation, shape) * weight
				multi = var*(individual[1::4][i])
				first=False
			else:
				multi += var*(individual[1::4][i])
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
			if distances is not None:
				hof.plotEdgeModel(distances, edf, oname)
			hof.plotPairwiseModel(gendist, r, oname)

	return(model_num, r, multi)

def modelAverage(pool, bests, base=""):
	#build model list and run Circuitscape on each model
	models=list()
	mod_num=0
	weights=dict()
	for index, row in bests.iterrows():
		n=len(predictors.columns)*4
		model_string = row.iloc[1:n+1].to_list()
		models.append([mod_num, model_string])
		weights[mod_num]=row["akaike_weight"]
		mod_num+=1

	print("Model-averaging across",mod_num,"resistance models...")

	# evaluate models in parallel
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
		if distances is not None:
			hof.plotEdgeModel(distances, edf, oname)
		hof.plotPairwiseModel(gendist, matrix_avg, oname)


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

# def initialize_worker(params, proc_num):
# 	global my_number
# 	my_number = proc_num
# 	#make new random seed, as seed+Process_number
# 	local_seed = int(params.seed)+my_number
# 	random.seed(local_seed)
#
# 	#print("")
# 	load_data(params, my_number)
# 	return(local_seed)

# DEPRECATED -- CIRCUITSCAPE VERSION
def initialize_worker(params, proc_num):
	global my_number
	my_number = proc_num
	#make new random seed, as seed+Process_number
	local_seed = int(params.seed)+my_number
	random.seed(local_seed)
	print("Worker",proc_num,"initializing...\n", flush=True)

	# global jl
	# from julia.api import Julia
	# from julia import Base, Main
	# from julia.Main import println, redirect_stdout
	#
	# #establish connection to julia
	# print("Worker",proc_num,"connecting to Julia...\n", flush=True)
	# if params.sys_image:
	# 	jl = Julia(init_julia=True, sysimage=params.sys_image, julia=params.julia, compiled_modules=params.compiled_modules)
	# else:
	# 	print("julia", flush=True)
	# 	jl = Julia(init_julia=True, julia=params.julia, compiled_modules=params.compiled_modules)
	#
	# if my_number == 0:
	# 	print("Loading Circuitscape in Julia...\n", flush=True)
	# #jl.eval("using Pkg;")
	# jl.eval("using Circuitscape; using Suppressor;")
	# Main.eval("stdout")

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
		rivers = gpd.read_file(shapefile)
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
		# 	pos[n] = n
		# color_map = []
		# for node in K:
		# 	if node in points.values():
		# 		color_map.append("blue")
		# 	else:
		# 		color_map.append("black")
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

	points = readPointCoords(params.coords)
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
	# 	ofh=params.out+".rawResistMat.txt"
	# 	with np.printoptions(precision=0, suppress=True):
	# 		np.savetxt(ofh, ref, delimiter="\t")

	# build incidence matrix
	inc_matrix = incidenceMatrix(graph, points, params.length_col, id_col, edge_order=names)
	if my_number == 0:
		ofh=params.out+".incidenceMatrix.txt"
		with np.printoptions(precision=0, suppress=True):
			np.savetxt(ofh, inc_matrix, delimiter="\t", fmt='%i')

	#read genetic distances
	gendist = parseInputGenMat(graph, points_names, prefix=params.prefix, inmat=params.inmat, agg_method=params.pop_agg)
	if my_number == 0:
		ofh=params.out+".genDistMat.txt"
		with np.printoptions(precision=0, suppress=True):
			np.savetxt(ofh, gendist, delimiter="\t")

	# pick first site for each in points_names as final name
	ol=list()
	for p in points_names:
		name = points_names[p][0]
		names=",".join(points_names[p])
		points_names[p] = name
		ol.append([name, p, points[p], names])

	# write table mapping final label to input pop labels and node coordinates
	if my_number == 0:
		df = pd.DataFrame(ol, columns = ['Label', 'Index', 'Node', 'Names'])
		dtout = str(params.out) + ".pointsTable.txt"
		df.to_csv(dtout, sep="\t", index=False)

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
				# 	names.append(points.values()[ia][0])
				# 	names.append(points.values()[ib][0])
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
# 	"node to point"
# 	np=OrderedDict()
# 	seen=dict()
# 	node_idx=0
# 	for edge in graph.edges:
# 		if edge[0] not in seen.keys():
# 			if edge[0] in points.keys():
# 				#print("New node:",points[edge[0]], " -- Index:",node_idx)
# 				np[node_idx] = points[edge[0]]
# 			seen[edge[0]] = True
# 			node_idx+=1
# 		if edge[1] not in seen.keys():
# 			if edge[1] in points.keys():
# 				#print("New node:",points[edge[1]], " -- Index:",node_idx)
# 				np[node_idx] = points[edge[1]]
# 			seen[edge[1]] = True
# 			node_idx+=1
# 	return(np)

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


# #parallel version; actually returns a list of filenames
# #custom evaluation function
# def evaluate(individual):
# 	#vector to hold values across edges
# 	fitness=None
# 	multi=None
# 	first=True
#
# 	#build multi-surface
# 	for i, variable in enumerate(predictors.columns):
# 		#Perform variable transformations (if desired)
# 		#1)Scale to 0-10; 2) Perform desired transformation; 3) Re-scale to 0-10
# 		#add weighted variable data to multi
# 		if individual[0::2][i] == 1:
# 			if first:
# 				multi = predictors[variable]*(individual[1::2][i])
# 				first=False
# 			else:
# 				multi += predictors[variable]*(individual[1::2][i])
#
# 	#If no layers are selected, return a zero fitness
# 	if first:
# 		fitness = None
# 	else:
# 		#Rescale multi for circuitscape
# 		multi = rescaleCols(multi, 1, 10)
#
# 		#write circuitscape inputs
# 		oname=".temp"+str(params.seed)+"_"+str(random.randint(1, 100000))
# 		#oname=".temp"+str(params.seed)
# 		focal=True
# 		if params.cstype=="edgewise":
# 			focal=False
# 		cs.writeCircuitScape(oname, graph, points, multi, focalPoints=focal, fromAttribute=None)
# 		cs.writeIni(oname, cholmod=params.cholmod, parallel=int(params.CS_procs))
#
# 		fitness=oname
# 	#return fitness value
# 	return(fitness)
#
# def parallel_eval(jl, ini_list, cstype):
# 	#Call circuitscape from pyjulia
# 	results = cs.evaluateIniParallel(jl, ini_list)
#
# 	#parse circuitscape output
# 	fitnesses = list()
# 	for ini in ini_list:
# 		fitness = 0
# 		if ini is None:
# 			fitness = float('-inf')
# 		else:
# 			if cstype=="edgewise":
# 				res = cs.parseEdgewise(ini, distances)
# 				fitness = res[params.fitmetric][0]
# 			else:
# 				res = cs.parsePairwise(ini, gendist)
# 				fitness = res[params.fitmetric][0]
# 			#cs.cleanup(ini)
# 		fitnesses.append(fitness)
# 	return(fitnesses)

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

# Version for doing each individual serially
# #custom evaluation function
def evaluate(individual):
	#print("evaluate - Process",my_number)
	#vector to hold values across edges
	fitness=0
	res=None
	multi=None
	first=True

	#build multi-surface
	for i, variable in enumerate(predictors.columns):
		#Perform variable transformations (if desired)
		#1)Scale to 0-10; 2) Perform desired transformation; 3) Re-scale to 0-10
		#add weighted variable data to multi
		if individual[0::4][i] == 1:
			#print("Before:", predictors[variable])
			var = transform(predictors[variable], individual[2::4][i], individual[3::4][i])
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


		oname=".temp_"+str(os.path.basename(params.out))+"_p"+str(my_number)
		focal=True

		# DEPRECATED
		#write circuitscape inputs
		#oname=".temp"+str(params.seed)+"_"+str(mp.Process().name)
		# if params.cstype=="edgewise":
		# 	focal=False
		# cs.writeCircuitScape(oname, graph, points, multi, id_col=id_col, focalPoints=focal, fromAttribute=None)
		# cs.writeIni(oname, cholmod=params.cholmod, parallel=int(params.CS_procs))
		# #Call circuitscape from pyjulia
		# print("evaluate", flush=True)
		# cs.evaluateIni(jl, oname)

		#parse circuitscape output
		# if params.cstype=="edgewise":
		# 	res = cs.parseEdgewise(oname, distances)
		# 	fitness = res[params.fitmetric][0]
		# 	res=list(res.iloc[0])
		# else:

		# res = cs.parsePairwise(oname, gendist)
		# fitness = res[params.fitmetric][0]
		# res=list(res.iloc[0])

		# evaluate using simple resistance distance
		r, res = rd.parsePairwise(points, inc_matrix, multi, gendist)
		fitness = res[params.fitmetric][0]
		res=list(res.iloc[0])

		# oname=".rdist_"+str(os.path.basename(params.out))+"_p"+str(my_number)
		# np.savetxt(oname, r, delimiter="\t")
	#return fitness value
	return(fitness, res)

#custom mutation function
def mutate(individual, indpb):
	for (i, variable) in enumerate(predictors.columns):
		if random.random() < indpb:
			individual[0::4][i] = toolbox.feature_sel()
		if random.random() < indpb:
			individual[1::4][i] = toolbox.feature_weight()
		if random.random() < indpb:
			individual[2::4][i] = toolbox.feature_transform()
		if random.random() < indpb:
			individual[3::4][i] = toolbox.feature_shape()
	return(individual,)

def initGA(toolbox, params):
	#register attributes
	toolbox.register("feature_sel", random.randint, 0, 1)

	if params.posWeight:
		toolbox.register("feature_weight", random.uniform, params.min_weight, 1.0)
	if params.fixWeight:
		toolbox.register("feature_weight", random.uniform, 1.0, 1.0)
	if not params.fixWeight and not params.posWeight:
		toolbox.register("feature_weight", random.uniform, -1.0, 1.0)
	if not params.fixShape:
		toolbox.register("feature_transform", random.randint, 0, 8)
		toolbox.register("feature_shape", random.randint, 1, params.max_shape)
	else:
		toolbox.register("feature_transform", random.randint, 0, 0)
		toolbox.register("feature_shape", random.randint, 1, 1)

	#register type for individuals
	#these consist of chromosomes of i variables x j attributes (above)
	toolbox.register("individual", tools.initCycle, creator.Individual,(toolbox.feature_sel, toolbox.feature_weight, toolbox.feature_transform, toolbox.feature_shape), n=len(params.variables))
	#toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.feature_sel, n=len(params.variables))

	#register type for populations
	#these are just a list of individuals
	toolbox.register("population", tools.initRepeat, list, toolbox.individual)

	#register custom evaluation function
	toolbox.register("evaluate", evaluate)

	#register mating function
	toolbox.register("mate", tools.cxTwoPoint)

	#register mutation function
	toolbox.register("mutate", mutate, indpb=params.indpb) #NOTE: make indpb an argument

	#register tournament function
	toolbox.register("select", tools.selTournament, tournsize=5) #NOTE: Make tournsize an argument


#Call main function
if __name__ == '__main__':
	main()
