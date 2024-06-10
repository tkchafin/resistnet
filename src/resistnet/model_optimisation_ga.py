# import numpy as np
# from queue import Empty
# from multiprocessing import Process, Queue
# from deap import base, creator, tools
# import matplotlib.pyplot as plt
# import seaborn as sns
# import matplotlib.patches as patches
# from matplotlib.backends.backend_pdf import PdfPages

# from resistnet.resistance_network import ResistanceNetwork
# from resistnet.resistance_network import ResistanceNetworkWorker
# from resistnet.samc_network import ResistanceNetworkSAMC
# from resistnet.samc_network import ResistanceNetworkSAMCWorker
# from resistnet.hall_of_fame import HallOfFame

# class ModelRunnerGA(ModelRunner):
#     """
#     Class to handle running the genetic algorithm for optimizing a ResistNet
#     model.

#     Attributes:
#         seed (int): Seed for random number generator.
#         resistance_network (ResistanceNetwork): An instance of
#                                                 ResistanceNetwork or subclass.
#         verbose (bool): Flag to control verbosity of output.
#         task_queue (Queue): Queue for tasks.
#         result_queue (Queue): Queue for results.
#         workers (list): List of worker processes.
#         bests (HallOfFame): Hall of fame for best models.
#         toolbox (deap.base.Toolbox): Toolbox for genetic algorithm.
#         logger (list): List to store logging information.
#         Various other GA parameters.
#     """

#     def __init__(self, resistance_network, seed=1234, verbose=True):
#         """
#         Initializes a ModelRunner with a given resistance network, seed, and
#         verbosity.

#         Args:
#             resistance_network (ResistanceNetwork or subclass): The resistance
#                                                           network to optimize.
#             seed (int): Seed for random number generator.
#             verbose (bool): Flag to control verbosity of output.
#         """
#         super().__init__(resistance_network, seed, verbose)
#         self.toolbox = None

#         self.cxpb = None
#         self.mutpb = None
#         self.indpb = None
#         self.popsize = None
#         self.maxpopsize = None
#         self.tournsize = None

#     def initialize_ga(self):
#         """
#         Initializes the genetic algorithm components and population.
#         """
#         self.toolbox = base.Toolbox()

#         # Initialize the DEAP GA components
#         creator.create("FitnessMax", base.Fitness, weights=(1.0,))
#         creator.create("Individual", list, fitness=creator.FitnessMax)

#         if self.verbose:
#             print("Initializing genetic algorithm parameters...\n")

#         self.init_ga_attributes()

#         # Initialize population
#         popsize = int(len(self.resistance_network.variables) *
#                       self.num_attributes * 15)
#         if self.popsize:
#             popsize = int(self.popsize)
#         if popsize > self.maxpopsize:
#             popsize = int(self.maxpopsize)

#         if self.verbose:
#             print("Establishing a population of size:", str(popsize))
#         self.population = self.toolbox.population(n=popsize)

#     def set_ga_parameters(self, mutpb, cxpb, indpb, popsize, maxpopsize,
#                           fixWeight, fixShape,
#                           min_weight, max_shape, max_hof_size, tournsize,
#                           fitmetric, awsum, only_keep, use_full, verbose,
#                           report_all):
#         """
#         Sets the genetic algorithm parameters.

#         Args:
#             mutpb (float): Probability of mutation per trait.
#             cxpb (float): Probability of being chosen for crossover.
#             indpb (float): Probability of mutation per individual.
#             popsize (int): Population size.
#             maxpopsize (int): Maximum population size.
#             fixWeight (bool): Constrain parameter weights to 1.0 (i.e.,
#                               unweighted).
#             fixShape (bool): Turn off feature transformation.
#             min_weight (float): Minimum allowable weight.
#             max_shape (float): Maximum shape value.
#             max_hof_size (int): Maximum Hall of Fame size.
#             tournsize (int): Tournament size.
#             fitmetric (str): Fitness metric to be used.
#             awsum (float): Cumulative Akaike weight threshold to retain top N
#                            models.
#             only_keep (bool): Only retain models where column "keep"=True.
#             verbose (bool): Flag to control verbosity of output.
#             report_all (bool): Flag to generate full outputs for all retained
#                                models.
#         """
#         self.mutpb = mutpb
#         self.cxpb = cxpb
#         self.indpb = indpb
#         self.popsize = popsize
#         self.maxpopsize = maxpopsize
#         self.fixWeight = fixWeight
#         self.fixShape = fixShape
#         self.min_weight = min_weight
#         self.max_shape = max_shape
#         self.max_hof_size = max_hof_size
#         self.tournsize = tournsize
#         self.fitmetric = fitmetric
#         self.awsum = awsum
#         self.only_keep = only_keep
#         self.use_full = use_full
#         self.verbose = verbose
#         self.report_all = report_all

#     def run_ga(self, maxgens=1, fitmetric="aic", burnin=0, deltaB=None,
#                deltaB_perc=0.001, indpb=0.5, mutpb=0.5, cxpb=0.5, nFail=50,
#                popsize=None, maxpopsize=1000, fixWeight=False,
#                fixShape=False, min_weight=0.0, max_shape=None,
#                max_hof_size=100, tournsize=10, awsum=0.95, only_keep=True,
#                use_full=False, fixed_params=None,
#                out=None, plot=True, verbose=True, report_all=False, threads=1):
#         """
#         Runs the genetic algorithm for optimizing the ResistNet model.

#         Args:
#             maxgens (int): Maximum number of generations.
#             fitmetric (str): Fitness metric to be used.
#             burnin (int): Number of initial generations to ignore for
#                           convergence check.
#             deltaB (float): Absolute threshold change in fitness for
#                             convergence.
#             deltaB_perc (float): Relative threshold change in fitness for
#                                  convergence.
#             indpb (float): Probability of mutation per individual.
#             mutpb (float): Probability of mutation per trait.
#             cxpb (float): Probability of being chosen for crossover.
#             nFail (int): Number of generations failing to improve to stop
#                          optimization.
#             popsize (int): Population size.
#             maxpopsize (int): Maximum population size.
#             fixWeight (bool): Constrain parameter weights to 1.0 (i.e.,
#                               unweighted).
#             fixShape (bool): Turn off feature transformation.
#             min_weight (float): Minimum allowable weight.
#             max_shape (float): Maximum shape value.
#             max_hof_size (int): Maximum Hall of Fame size.
#             tournsize (int): Tournament size.
#             awsum (float): Cumulative Akaike weight threshold to retain top N
#                            models.
#             only_keep (bool): Only retain models where column "keep"=True.
#             out (str): Output file prefix.
#             plot (bool): Flag to control plotting.
#             verbose (bool): Flag to control verbosity of output.
#             report_all (bool): Generate full outputs for all retained models.
#             threads (int): Number of parallel processors.

#         Raises:
#             Exception: Any exception encountered during the genetic algorithm
#                        run.
#         """
#         try:
#             # Set GA parameters
#             self.set_ga_parameters(mutpb, cxpb, indpb, popsize, maxpopsize,
#                                    fixWeight, fixShape,
#                                    min_weight, max_shape, max_hof_size,
#                                    tournsize, fitmetric, awsum, only_keep,
#                                    use_full, verbose, report_all)

#             # parse optional input fixed_params
#             if fixed_params:
#                 self.parse_fixed_params(fixed_params)

#             # write fixed params if provided
#             if out and self.fixed_params:
#                 self.write_fixed_params(str(out) + ".fixedParams.tsv")

#             # Initialize GA components and worker pool
#             self.initialize_ga()
#             self.start_workers(threads)

#             # Initialize population and set up GA run parameters
#             self.initialise_population()
#             cxpb, indpb = self.cxpb, self.indpb
#             fails, current_best = 0, None

#             # Run for maxgens generations
#             for g in range(1, maxgens + 1):
#                 if self.verbose:
#                     print("-- Generation %i --" % g)

#                 # Select the next generation individuals
#                 offspring = self.toolbox.select(
#                     self.population, len(self.population)
#                 )
#                 # Clone the selected individuals
#                 offspring = list(map(self.toolbox.clone, offspring))

#                 # Apply crossover and mutation on the offspring
#                 for child1, child2 in zip(offspring[::2], offspring[1::2]):
#                     if random.random() < self.cxpb:
#                         self.toolbox.mate(child1, child2)
#                         del child1.fitness.values
#                         del child2.fitness.values

#                 for mutant in offspring:
#                     if random.random() < self.indpb:
#                         self.toolbox.mutate(mutant)
#                         del mutant.fitness.values

#                 # Evaluate the individuals with an invalid fitness
#                 invalid_ind = [
#                     ind for ind in offspring if not ind.fitness.valid
#                 ]
#                 for ind_index, ind in enumerate(invalid_ind):
#                     self.task_queue.put(('evaluate', ind_index, list(ind)))

#                 # Collect fitness results and process them
#                 pop_list = []
#                 for _ in range(len(invalid_ind)):
#                     ind_index, fitness, res = self.result_queue.get()
#                     invalid_ind[ind_index].fitness.values = [fitness]
#                     # Process and append to pop_list
#                     ind_list = [fitness] + list(invalid_ind[ind_index])
#                     if res is not None:
#                         ind_list.extend(res)
#                     pop_list.append(ind_list)

#                 # Replace population with offspring
#                 self.population[:] = offspring

#                 # update Hall of Fame
#                 self.bests.check_population(pop_list)

#                 # Gather all the fitnesses in one list and print the stats
#                 fits = [ind.fitness.values[0] for ind in self.population]

#                 # Evaluate for stopping criteria after burn-in period
#                 if g > burnin:
#                     if fitmetric.upper() == "AIC":
#                         fits = [-element for element in fits]
#                         worst = max(fits)
#                         best = min(fits)
#                         if current_best is None:
#                             current_best = best
#                         else:
#                             current_best, fails = self.update_fails(
#                                 best, current_best, fails, deltaB, deltaB_perc,
#                                 minimize=True)
#                     else:
#                         worst = min(fits)
#                         best = max(fits)
#                         if current_best is None:
#                             current_best = best
#                         else:
#                             current_best, fails = self.update_fails(
#                                 best, current_best, fails, deltaB,
#                                 deltaB_perc, minimize=False
#                             )

#                     length = len(self.population)
#                     if length > 0:
#                         if any(math.isinf(fit) for fit in fits):
#                             # Set stats to NaN if any values are inf or -inf
#                             mean = variance = std = float('nan')
#                         else:
#                             mean = sum(fits) / length
#                             sum2 = sum(x * x for x in fits)
#                             variance = sum2 / length - mean**2

#                             std = (math.sqrt(variance)
#                                    if variance >= 0 else float('nan'))
#                     else:
#                         mean = variance = std = float('nan')

#                     if self.verbose:
#                         print("  Worst %s" % worst)
#                         print("  Best %s" % best)
#                         print("  Avg %s" % mean)
#                         print("  Std %s" % std)
#                         print("  nFails %s" % fails)
#                     self.logger.append([g, worst, best, mean, std])

#                     # Check for stopping conditions
#                     if g >= maxgens or fails > nFail:
#                         if self.verbose:
#                             print(
#                                 "Stopping optimization after",
#                                 g, "generations."
#                             )
#                             if g >= maxgens:
#                                 print("Reason: Exceeded maxGens")
#                             elif fails > nFail:
#                                 print(
#                                     "Reason: More than",
#                                     nFail,
#                                     "generations failed to improve"
#                                 )
#                         break

#             # output results
#             self.run_output(out, verbose, plot)

#         except Exception as e:
#             print(f"An error occurred: {e}")
#             traceback.print_exc()

#         finally:
#             # Terminate worker processes
#             self.terminate_workers()

#     @staticmethod
#     def update_fails(best, current_best, fails, delta_b, delta_b_perc,
#                      minimize=False):
#         """
#         Updates the fail count based on the best and current best fitness
#         scores.

#         Args:
#             best: The best fitness score.
#             current_best: The current best fitness score.
#             fails: The current fail count.
#             delta_b: Absolute improvement threshold.
#             delta_b_perc: Percentage improvement threshold.
#             minimize (bool): Flag indicating minimization or maximization.

#         Returns:
#             Tuple containing updated current best fitness score and fail count.
#         """
#         cur = current_best
#         f = fails
#         threshold_1 = current_best
#         threshold_2 = current_best

#         if minimize:
#             if best < current_best:
#                 cur = best
#             if delta_b is not None:
#                 threshold_1 = current_best - delta_b
#             if delta_b_perc is not None:
#                 threshold_2 = current_best - (current_best * delta_b_perc)
#             if best >= threshold_1 or best >= threshold_2:
#                 f += 1
#             else:
#                 f = 0
#         else:
#             if best > current_best:
#                 cur = best
#             if delta_b is not None:
#                 threshold_1 = current_best + delta_b
#             if delta_b_perc is not None:
#                 threshold_2 = current_best + (current_best * delta_b_perc)
#             if best <= threshold_1 or best <= threshold_2:
#                 f += 1
#             else:
#                 f = 0
#         return cur, f

#     def initialise_population(self):
#         """
#         Initialises the population by distributing tasks for evaluation and
#         processing the fitness results.
#         """
#         # Distribute tasks for evaluation
#         for ind_index, ind in enumerate(self.population):
#             self.task_queue.put(('evaluate', ind_index, list(ind)))

#         # Collect fitness results and process them
#         pop_list = []
#         for _ in range(len(self.population)):
#             ind_index, fitness, res = self.result_queue.get()
#             self.population[ind_index].fitness.values = [fitness]
#             ind_list = [fitness] + list(self.population[ind_index])
#             if res is not None:
#                 ind_list.extend(res)
#             pop_list.append(ind_list)

#         # Initialize Hall of Fame
#         self.bests = HallOfFame(
#             self.resistance_network._predictors.columns,
#             self.max_hof_size,
#             self.use_full,
#             pop_list)

#     def init_ga_attributes(self):
#         """
#         Initializes the genetic algorithm's attributes and registers necessary
#         functions with the toolbox, allowing dynamic specification of fixed
#         parameters per variable and attribute.
#         """

#         # Helper functions to generate attributes
#         def feature_generator(var_name, feature_type, lower, upper,
#                               is_int=False, custom_range=None):
#             if (self.fixed_params and
#                 var_name in self.fixed_params and
#                     feature_type in self.fixed_params[var_name]):
#                 # Return a lambda that always returns the fixed value
#                 return lambda: int(self.fixed_params[var_name][feature_type])
#             elif is_int:
#                 if custom_range is not None:
#                     return lambda: random.choice(custom_range)
#                 else:
#                     return lambda: random.randint(lower, upper)
#             else:
#                 return lambda: random.uniform(lower, upper)

#         # Registering attribute generators for each variable and each attribute
#         for var_name in self.resistance_network.variables:
#             self.toolbox.register(
#                 f"{var_name}_sel", feature_generator(
#                     var_name, "sel", 0, 1, is_int=True))
#             if self.fixWeight:
#                 self.toolbox.register(
#                     f"{var_name}_weight", feature_generator(
#                         var_name, "weight", 1.0, 1.0))
#             elif self.min_weight:
#                 self.toolbox.register(
#                     f"{var_name}_weight", feature_generator(
#                         var_name, "weight", self.min_weight, 1.0))
#             else:
#                 self.toolbox.register(
#                     f"{var_name}_weight", feature_generator(
#                         var_name, "weight", 0.0, 1.0))

#             # if not self.fixAsym:
#             #     self.toolbox.register(
#             #         f"{var_name}_asym", feature_generator(
#             #             var_name, "asym", 0, 1, is_int=True))
#             # else:
#             # self.toolbox.register(
#             #     f"{var_name}_asym", feature_generator(
#             #         var_name, "asym", 0, 0, is_int=True))

#             if not self.fixShape:
#                 transform_range = list(range(0, 9))
#                 self.toolbox.register(
#                     f"{var_name}_transform", feature_generator(
#                         var_name, "transform", 0, 8, is_int=True,
#                         custom_range=transform_range))
#                 self.toolbox.register(
#                     f"{var_name}_shape", feature_generator(
#                         var_name, "shape", 1, self.max_shape, is_int=True))
#             else:
#                 self.toolbox.register(
#                     f"{var_name}_transform", feature_generator(
#                         var_name, "transform", 0, 0, is_int=True))
#                 self.toolbox.register(
#                     f"{var_name}_shape", feature_generator(
#                         var_name, "shape", 0, 0, is_int=True))

#         # Define how to create an individual
#         def create_individual():
#             return creator.Individual(
#                 [getattr(self.toolbox, f"{var_name}_{attr}")()
#                     for var_name in self.resistance_network.variables
#                     for attr in [
#                         'sel', 'weight', 'transform', 'shape']]
#             )

#         # Registering individual and population creation methods
#         self.toolbox.register("individual", create_individual)
#         self.toolbox.register(
#             "population", tools.initRepeat, list, self.toolbox.individual)

#         # Evaluation, mating, mutation, and selection methods remain unchanged
#         self.toolbox.register("evaluate", lambda ind: ind)  # spoof
#         self.toolbox.register("mate", tools.cxTwoPoint)
#         self.toolbox.register("mutate", self.mutate)
#         self.toolbox.register(
#             "select", tools.selTournament, tournsize=self.tournsize)

#     def mutate(self, individual):
#         """
#         Custom mutation function for an individual that aligns with the dynamic
#         toolbox registration.
#         """
#         num_attributes = self.num_attributes
#         num_vars = len(individual) // num_attributes

#         for var_index in range(num_vars):
#             var_name = self.resistance_network.variables[var_index]
#             idx_base = var_index * num_attributes

#             if random.random() < self.mutpb:
#                 individual[idx_base] = getattr(
#                     self.toolbox, f"{var_name}_sel")()
#             if random.random() < self.mutpb:
#                 individual[idx_base + 1] = getattr(
#                     self.toolbox, f"{var_name}_weight")()
#             if random.random() < self.mutpb:
#                 individual[idx_base + 2] = getattr(
#                     self.toolbox, f"{var_name}_transform")()
#             if random.random() < self.mutpb:
#                 individual[idx_base + 3] = getattr(
#                     self.toolbox, f"{var_name}_shape")()

#         return individual,