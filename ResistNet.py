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
import timeit

from deap import base, creator, tools, algorithms

from resistnet.functions import *
from resistnet.params import parseArgs
import resistnet.hall_of_fame as hof

def main():

	global params
	params = parseArgs()

	#seed random number generator
	#random.seed(params.seed)
	if not params.seed:
		params.seed=datetime.now().timestamp()
	random.seed(params.seed)

	#IF INPUT WAS SHAPEFILE, MASTER PROC GENERATES INPUTS FOR WORKERS HERE
	# READ SHAPEFILE
	# FIND SUB_GRAPH
	# WRITE NETWORK AND TXT OUTPUT FILES

	pool = mp.Pool(processes=params.GA_procs)

	process_list = range(1, int(params.GA_procs)+1)
	func = partial(initialize_worker, params)
	results = pool.map(func, process_list)

	#load data for master process
	global my_number
	my_number = 0

	load_data(params, 0)

	sys.exit()

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

	#if params.modavg:
	if params.modavg:
		modelAverageCS(pool, bests.getHOF(only_keep=params.only_keep), base=params.out, plot=params.plot, report_all=params.report_all) #set to true for production
	writeMatrix((str(params.out)+".genDistMat.tsv"), gendist, list(node_point_dict.values()))


	#clean up temp files
	oname=".temp_"+str(params.out)+"*"
	for filename in glob.glob(oname):
		os.remove(filename)

	pool.close()


#Call main function
if __name__ == '__main__':
	main()
