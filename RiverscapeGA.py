import sys
import os, glob
import itertools
import math
import random
import getopt
import scipy as sp
import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from datetime import datetime
from functools import partial
from collections import OrderedDict
from sortedcontainers import SortedDict

#julia

import timeit

#genetic algorithms
from deap import base, creator, tools, algorithms

#autoStreamTree packages
import riverscape.hall_of_fame as hof
from riverscape.acg_menu import parseArgs
import riverscape.circuitscape_runner as cs
import riverscape.transform as trans
import riverscape.stream_plots as splt


def main():
	
	global params
	params = parseArgs()
	
	#seed random number generator
	#random.seed(params.seed)
	if not params.seed:
		params.seed=int(datetime.now())
	random.seed(params.seed)
	
	pool = mp.Pool(processes=params.GA_procs)
	
	process_list = range(1, int(params.GA_procs)+1)
	func = partial(initialize_worker, params)
	results = pool.map(func, process_list)
	
	#load data for master process
	global my_number 
	my_number = 0
	load_data(params, 0)
	
	
	#print(len(list(graph.edges())))
	#print(results)
	#sys.exit()
	
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
			
			if params.fitmetric=="aic":
				fits = [element * -1.0 for element in fits]
				worst = max(fits)
				best=min(fits)
				if current_best is None:
					current_best=best
				else:
					(current_best, fails) = updateFails(best, current_best, fails, params.deltaB, params.deltaB_perc, minimize=True)
			else:
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
	
	# if params.fitmetric == 'aic':
	# 	bests.correct_aic_fitness()
	bests.delta_aic()
	bests.akaike_weights()
	bests.cumulative_akaike(threshold=params.awsum)
	bests.relative_variable_importance(params.only_keep)
	bests.printHOF()
	bests.printRVI()
	if params.plot:
		bests.plot_ICprofile(params.out)
		bests.plotMetricPW(params.out)
		bests.plotVariableImportance(params.out)
	
	#write log of fitnesses
	logDF=pd.DataFrame(logger, columns=["Generation", "Worst", "Best", "Mean", "Stdev"])
	logDF.to_csv((str(params.out)+".FitnessLog.tsv"), sep="\t", header=True, index=False)
	del logger
	del logDF
	
	#if params.modavg:
	if params.modavg:
		modelAverageCS(pool, bests.getHOF(only_keep=params.only_keep), base=params.out, plot=params.plot, report_all=params.report_all) #set to true for production
	writeMatrix((str(params.out)+".genDistMat.tsv"), gen, list(node_point_dict.values()))
	
	#write hall of fame to file
	bests.writeModelSummary(params.out)
	
	#clean up temp files
	oname=".temp_"+str(params.out)+"*"
	for filename in glob.glob(oname):
		os.remove(filename)
	
	pool.close()

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

def evaluate_ma(stuff):
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

	#write circuitscape inputs
	oname=str(params.out)+"_Model-"+str(model_num)

	cs.writeCircuitScape(oname, graph, points, multi, focalPoints=False, fromAttribute=None)
	cs.writeIni(oname, cholmod=params.cholmod, parallel=int(params.CS_procs))
	
	#Call circuitscape from pyjulia
	cs.evaluateIni(jl, oname)
	return(model_num)

def modelAverageCS(pool, bests, plot=False, base="", report_all=False):
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
	
	model_nums = pool.map(evaluate_ma, models)
	
	edge_avg=None
	matrix_avg=np.zeros(shape=(len(points), len(points)))
	
	for model_num in model_nums:
		oname=str(base)+"_Model-"+str(model_num)
		
		weight = bests.iloc[model_num]["akaike_weight"]
		
		#Extract resistance for each edge
		#print(distances)
		edge_r = cs.parseEdgewise(oname, distances, return_resistance=True)
		#print(edge_r)
		if edge_avg is None:
			edge_avg=(edge_r*weight)
		else:
			edge_avg = np.add(edge_avg, (edge_r*weight))
		
		#extract PW matrix from full result matrix
		matrix_r = cs.parsePairwiseFromAll(oname, gendist, node_point_dict, return_resistance=True)
		matrix_avg += (matrix_r*weight)
		
		if report_all:
			writeEdges(oname, edge_avg, edge_ids)
			writeMatrix(oname, matrix_avg, list(node_point_dict.values()))
			if plot:
				edf=pd.DataFrame(list(zip(edge_ids, edge_r)), columns=["EDGE_ID", "Resistance"])
				splt.plotEdgesToStreams(graph, edf, (str(params.prefix)+".streamTree.shp"), oname)
				hof.plotEdgeModel(distances, edge_r, oname)
				hof.plotPairwiseModel(gendist, matrix_r, oname)
		
	
	oname=str(base) + ".Model-Average"
	writeEdges(oname, edge_avg, edge_ids)
	writeMatrix(oname, matrix_avg, list(node_point_dict.values()))
	
	if plot:
		edf=pd.DataFrame(list(zip(edge_ids, edge_avg)), columns=["EDGE_ID", "Resistance"])
		splt.plotEdgesToStreams(graph, edf, (str(params.prefix)+".streamTree.shp"), oname)
		hof.plotEdgeModel(distances, edge_avg, oname)
		hof.plotPairwiseModel(gendist, matrix_avg, oname)
		
		#streamtree-style plot

def writeEdges(oname, edge, ids, dist=None):
	df=pd.DataFrame(list(zip(ids, edge)), columns=["EDGE_ID", "Resistance"])
	if dist is not None:
		df["FittedD"]=dist
	df.to_csv((str(oname)+".ResistanceEdges.tsv"), sep="\t", header=True, index=False)
	

def writeMatrix(oname, mat, ids):
	df=pd.DataFrame(mat, columns=ids, index=ids)
	df.to_csv((str(oname)+".ResistanceMatrix.tsv"), sep="\t", header=True, index=True)


def initialize_worker(params, proc_num):
	global my_number 
	my_number = proc_num
	#make new random seed, as seed+Process_number
	local_seed = int(params.seed)+my_number
	random.seed(local_seed)
	
	global jl
	from julia.api import Julia
	from julia import Base, Main
	from julia.Main import println, redirect_stdout
	
	#establish connection to julia
	print("Worker",proc_num,"connecting to Julia...\n")
	if params.sys_image is not None:
		jl = Julia(init_julia=False, sysimage=params.sys_image, julia=params.julia, compiled_modules=params.compiled_modules)
	else:
		jl = Julia(init_julia=False, julia=params.julia, compiled_modules=params.compiled_modules)
	
	if my_number == 0:
		print("Loading Circuitscape in Julia...\n")
	#jl.eval("using Pkg;")
	jl.eval("using Circuitscape; using Suppressor;")
	#Main.eval("stdout")
	
	load_data(params, my_number)
	
	return(local_seed)

def load_data(params, proc_num):

	#make "local" globals (i.e. global w.r.t each process)
	global graph
	global distances
	global predictors
	global inc_matrix
	global points
	global gendist
	global node_point_dict
	global edge_ids
	my_number = proc_num
	
	#read autoStreamTree outputs
	if params.network is not None:
		if my_number == 0:
			print("Reading network from: ", str(params.network))
		graph = readNetwork(params.network)
	else:
		if my_number == 0:
			print("Reading network from: ", (str(params.prefix)+".network"))
		graph = readNetwork((str(params.prefix)+".network"))
	if my_number==0:
		print("Reading autoStreamTree results from:", (str(params.prefix)+".streamTree.txt"))
	(distances, predictors, edge_ids) = readStreamTree((str(params.prefix)+".streamTree.txt"), params.variables, params.force)
	points = readPointCoords((str(params.prefix)+".pointCoords.txt"))
	
	#make sure points are snapped to the network
	snapped=SortedDict()
	for point in points.keys():
		if point not in graph.nodes():
			node=snapToNode(graph, point)
			if my_number == 0:
				print("Point not found in graph, snapping to nearest node:", point, " -- ", node)
			snapped[tuple(node)]=points[point]
		else:
			snapped[tuple(point)]=points[point]
	points = snapped
	del snapped
	
	node_point_dict=nodes_to_points(graph, points)
	
	#read genetic distances
	if params.predicted:
		if my_number == 0:
			print("Reading incidence matrix from: ", (str(params.prefix)+".incidenceMatrix.txt"))
		inc_matrix = readIncidenceMatrix((str(params.prefix)+".incidenceMatrix.txt"))
		gendist = generatePairwiseDistanceMatrix(graph, points, inc_matrix, distances)
	else:
		gendist = parseInputGenMat(graph, points, prefix=params.prefix, inmat=params.inmat)


def checkFormatGenMat(mat, order):
	if os.path.isfile(mat):
		#read and see if it has the correct dimensions
		inmat = pd.read_csv(mat, header=0, index_col=0, sep="\t")
		#if correct dimensions, check if labelled correctly
		if len(inmat.columns) >= len(order):
			if set(list(inmat.columns.values)) != set(list(inmat.columns.values)):
				#print("columns and rows don't match")
				return(None)
			# elif set(list(inmat.columns.values)) != set(order):
			# 	#print("Oh no! Input matrix columns and/ or rows don't appear to be labelled properly. Please provide an input matrix with column and row names!")
			# 	return(None)
			else:
				#this must be the one. Reorder it and return
				#print("Reading genetic distances from input matrix:",indmat)
				formatted = inmat.reindex(order)
				formatted = formatted[order]
				#print(formatted)
				gen = formatted.to_numpy()
				return(gen)
		#otherwise, skip and try the popgenmat
		else:
			#print("wrong number of columns")
			return(None)
	else:
		return(None)
		
def parseInputGenMat(graph, points, prefix=None, inmat=None):
	order = getNodeOrder(graph, points, as_list=True)
	#if no input matrix provided, infer from autoStreamTree output
	if inmat is None:
		#check if pop and ind mats both exist
		indmat = str(prefix) + ".indGenDistMat.txt"
		popmat = str(prefix) + ".popGenDistMat.txt"
		#if indmat exists
		ind = checkFormatGenMat(indmat, order)
		pop = checkFormatGenMat(popmat, order)
		#print(pop)
		#print(order)
		#print(pop.dtype)
		if pop is not None:
			return(pop)
		elif ind is not None:
			return(ind)
		else:
			print("Failed to read autoStreamTree genetic distance matrix.")
			sys.exit()
	else:
		#read input matrix instead
		gen = checkFormatGenMat(inmat)
		if gen is not None:
			return(gen)
		else:
			print("Failed to read input genetic distance matrix:",inmat)
			sys.exit()

def nodes_to_points(graph, points):
	"node to point"
	np=OrderedDict()
	seen=dict()
	node_idx=0
	for edge in graph.edges:
		if edge[0] not in seen.keys():
			if edge[0] in points.keys():
				#print("New node:",points[edge[0]], " -- Index:",node_idx)
				np[node_idx] = points[edge[0]]
			seen[edge[0]] = True
			node_idx+=1
		if edge[1] not in seen.keys():
			if edge[1] in points.keys():
				#print("New node:",points[edge[1]], " -- Index:",node_idx)
				np[node_idx] = points[edge[1]]
			seen[edge[1]] = True
			node_idx+=1
	return(np)

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

def generatePairwiseDistanceMatrix(graph, points, inc_matrix, distances):
	node_dict=getNodeOrder(graph, points, as_index=True)
	gen=np.zeros(shape=(len(points), len(points)))
	inc_row=0
	#print(node_dict.keys())
	#print(node_dict[tuple(list(points.keys())[0])])
	for ia, ib in itertools.combinations(range(0,len(points)),2):
		inc_streams=inc_matrix[inc_row,]
		#print(distances*inc_streams)
		d=sum(distances*inc_streams)
		#print(d)
		inc_row+=1
		#print(node_dict)
		#print(ia, " -- ", list(points.keys())[ia], " -- ", node_dict[list(points.keys())[ia]])
		#print(ib, " -- ", list(points.keys())[ib], " -- ", node_dict[list(points.keys())[ib]])
		gen[node_dict[list(points.keys())[ia]], node_dict[list(points.keys())[ib]]] = d
	#print(gen)
	return(gen)


def readIncidenceMatrix(inc):
	df = pd.read_csv(inc, header=None, index_col=False, sep="\t")
	return(df.to_numpy())
	
def readNetwork(network):
	graph=nx.OrderedGraph(nx.read_gpickle(network).to_undirected())
	return(graph)

#Input: Tuple of [x,y] coordinates
#output: Closest node to those coordinates
def snapToNode(graph, pos):
	#rint("closest_node call:",pos)
	nodes = np.array(graph.nodes())
	node_pos = np.argmin(np.sum((nodes - pos)**2, axis=1))
	#print(nodes)
	#print("closest to ", pos, "is",tuple(nodes[node_pos]))
	return (tuple(nodes[node_pos]))

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
			d[coords]=name
	return(d)

def readStreamTree(streamtree, variables, force=None):
	df = pd.read_csv(streamtree, header=0, index_col=False, sep="\t")
	names=sorted(set(df["EDGE_ID"].tolist()))
	
	df = df.groupby('EDGE_ID').agg('mean')
	
	#get distances (as a list, values corresponding to nodes)
	if force:
		dist=df[force].tolist()
	else:
		#get locus columns
		filter_col = [col for col in df if col.startswith('locD_')]
		data = df[filter_col]
		
		#aggregate distances
		dist = list(data.mean(axis=1))
		#print(dist)
	
	env=trans.rescaleCols(df[variables], 0, 10)
	return(dist, env, names)

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
# 		#	NOTE: Get main implementation working first
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
		d=trans.revRicker(dat, shape, 10)
	elif transformation == 3:
		d=trans.invRicker(dat, shape, 10)
	elif transformation == 4:
		d=trans.revInvRicker(dat, shape, 10)
	elif transformation == 5:
		d=trans.monomolecular(dat, shape, 10)
	elif transformation == 6:
		d=trans.revMonomolecular(dat, shape, 10)
	elif transformation == 7:
		d=trans.invMonomolecular(dat, shape, 10)
	elif transformation == 8:
		d=trans.revMonomolecular(dat, shape, 10)
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
		#	NOTE: Get main implementation working first
		#add weighted variable data to multi
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

	#If no layers are selected, return a zero fitness
	if first:
		fitness = float('-inf')
	else:
		#Rescale multi for circuitscape
		#print("Multi:",multi)
		multi = trans.rescaleCols(multi, 1, 10)

		#write circuitscape inputs
		#oname=".temp"+str(params.seed)+"_"+str(mp.Process().name)
		oname=".temp_"+str(params.out)+"_p"+str(my_number)
		#print(oname)
		focal=True
		if params.cstype=="edgewise":
			focal=False
		cs.writeCircuitScape(oname, graph, points, multi, focalPoints=focal, fromAttribute=None)
		cs.writeIni(oname, cholmod=params.cholmod, parallel=int(params.CS_procs))
		#print("evaluate")
		#Call circuitscape from pyjulia
		cs.evaluateIni(jl, oname)

		#parse circuitscape output
		if params.cstype=="edgewise":
			res = cs.parseEdgewise(oname, distances)
			fitness = res[params.fitmetric][0]
			res=list(res.iloc[0])
		else:
			res = cs.parsePairwise(oname, gendist)
			fitness = res[params.fitmetric][0]
			res=list(res.iloc[0])

	#return fitness value
	return(fitness, res)

#custom mutation function
#To decide: Should the current state inform the next state, or let it be random?
#May depend on the "gene"?
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
	toolbox.register("feature_weight", random.uniform, -1.0, 1.0)
	toolbox.register("feature_transform", random.randint, 0, 8)
	toolbox.register("feature_shape", random.randint, 1, 100)
	
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