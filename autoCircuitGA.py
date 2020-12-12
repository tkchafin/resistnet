import sys
import os
import itertools
import math
import random
import getopt
import scipy as sp
import numpy as np
import pandas as pd
import networkx as nx
import julia as jl
from functools import partial

from deap import base, creator, tools, algorithms

#autoStreamTree packages
import circuitscape_runner as cs
from acg_menu import parseArgs

def main():

	params = parseArgs()
	params.prefix="out3"
	params.force="fittedD"
	params.variables = ["tmp_dc_cmn", "aet_mm_cyr", "USE"]
	
	#read autoStreamTree results
	graph = readNetwork((str(params.prefix)+".network"))
	(distances, predictors) = readStreamTree((str(params.prefix)+".streamtree.txt"), params.variables, params.force)
	inc_matrix = readIncidenceMatrix((str(params.prefix)+".incidenceMatrix.txt"))
	
	# print(distances)
	# print(predictors)
	# print(inc_matrix)
	
	#Options that need to be set: 
	#1) Do we regress pairwise distance x resistance OR edge-wise distance x resistance
	#2) Use fittedD value or aggregate locD_ values?
	#3) Can do bootstraps by bootstrapping locus distances (from autoStreamTree)
	#4) Fitness calculation: AIC, loglik, or R^2? 
	
	#NOTES:
	#transformations - re-write code from ResistanceGA
	#inverse transforms aren't needed, because we allow negative weighting
	
	#Design of the genetic algorithm: 
	#There are multiple types of attributes occupying a 'chromosome':
	#1) Boolean - feature selection (decides if a variable is included or not)
	#2) Float - Provides a weight to each feature in the final additive resistence model
	#3) Categorical - specifies what type of transformation is associated with a variable
	#4) Float2 - A shape parameter associated with the transformation type
	
	#initialize a single-objective GA
	creator.create("FitnessMax", base.Fitness, weights=(1.0,))
	creator.create("Individual", list, fitness=creator.FitnessMax)
	
	#toolbox
	global toolbox
	toolbox = base.Toolbox()
	
	#register GA attributes and type variables
	initGA(toolbox, params)
	
	#initialize population
	pop = toolbox.population(n=100)
	
	# Evaluate the entire population
	fitnesses = list(map(toolbox.evaluate, pop))
	for ind, fit in zip(pop, fitnesses):
		ind.fitness.values = fit
	
	# CXPB  is the probability with which two individuals are crossed
	# MUTPB is the probability for mutating an individual
	cxpb, mutpb = 0.5, 0.2
	
	# Extracting all the fitnesses of population
	fits = [ind.fitness.values[0] for ind in pop]
	
	# Variable keeping track of the number of generations
	g = 0

	# Begin the evolution
	while max(fits) < 100 and g < 1000:
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
		
		#evaluate individuals with invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		fitnesses = map(toolbox.evaluate, invalid_ind)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit
		
		#replace population with offspring
		pop[:] = offspring
		
		# Gather all the fitnesses in one list and print the stats
		fits = [ind.fitness.values[0] for ind in pop]

		length = len(pop)
		mean = sum(fits) / length
		sum2 = sum(x*x for x in fits)
		std = abs(sum2 / length - mean**2)**0.5

		print("  Min %s" % min(fits))
		print("  Max %s" % max(fits))
		print("  Avg %s" % mean)
		print("  Std %s" % std)
	best = pop[np.argmax([toolbox.evaluate(x) for x in pop])]
	print(best)

def readIncidenceMatrix(inc):
	print("Reading incidence matrix from: ", inc)
	df = pd.read_csv(inc, header=None, index_col=False, sep="\t")
	return(df.to_numpy())
	
def readNetwork(network):
	print("Reading network from: ", network)
	graph=nx.OrderedGraph(nx.read_gpickle(network).to_undirected())
	return(graph)

def readStreamTree(streamtree, variables, force=None):
	print("Reading autoStreamTree results from:", streamtree)
	df = pd.read_csv(streamtree, header=0, index_col=False, sep="\t")
	
	df = df.groupby('EDGE_ID').agg('mean')

	#get distances (as a list, values corresponding to nodes)
	if force:
		dist=df[force].tolist()
	else:
		#get locus columns
		filter_col = [col for col in df if col.startswith('locD_')]
		data = df[filter_col]
		
		#aggregate distances
		dist = data.mean(axis=1).tolist()
	
	env=df[variables]
	env -= env.min()
	env /= env.max()
	
	return(dist, env)

#custom evaluation function
def evaluate(individual):
	fitness=sum(individual[0::2])
	fitness-=individual[1]
	fitness+=individual[3]
	fitness*=individual[-1]
	return(fitness,)

#custom mutation function
def mutate(individual, indpb):
	if random.random() < indpb:
		individual[0] = toolbox.feature_sel()
		individual[1] = toolbox.feature_weight()
		individual[2] = toolbox.feature_sel()
		individual[3] = toolbox.feature_weight()
		individual[4] = toolbox.feature_sel()
		individual[5] = toolbox.feature_weight()
	return(individual,)

def initGA(toolbox, params):
	#register attributes
	toolbox.register("feature_sel", random.randint, 0, 1)
	toolbox.register("feature_weight", random.uniform, -1.0, 1.0)
	#toolbox.register("feature_transform", random.randint, 0, 1)
	#toolbox.register("feature_scale", random.randint, 0, 1)
	
	#register type for individuals 
	#these consist of chromosomes of i variables x j attributes (above)
	toolbox.register("individual", tools.initCycle, creator.Individual,(toolbox.feature_sel, toolbox.feature_weight), n=len(params.variables))
	#toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.feature_sel, n=len(params.variables))
	
	#register type for populations
	#these are just a list of individuals
	toolbox.register("population", tools.initRepeat, list, toolbox.individual)	

	#register custom evaluation function
	toolbox.register("evaluate", evaluate)
	
	#register mating function
	toolbox.register("mate", tools.cxTwoPoint)
	
	#register mutation function
	toolbox.register("mutate", mutate, indpb=0.05) #NOTE: make indpb an argument
	
	#register tournament function
	toolbox.register("select", tools.selTournament, tournsize=3) #NOTE: Make tournsize an argument

#Call main function
if __name__ == '__main__':
	main()