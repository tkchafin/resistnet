import traceback
import random
import math
import pandas as pd
import numpy as np
from queue import Empty
from multiprocessing import Process, Queue
from deap import base, creator, tools

from resistnet.resistance_network import ResistanceNetwork
from resistnet.resistance_network import ResistanceNetworkWorker
from resistnet.samc_network import ResistanceNetworkSAMC
from resistnet.samc_network import ResistanceNetworkSAMCWorker
from resistnet.hall_of_fame import HallOfFame


class ModelRunner:
    """
    Class to handle running the genetic algorithm for optimizing a ResistNet
    model.

    Attributes:
        seed (int): Seed for random number generator.
        resistance_network (ResistanceNetwork): An instance of
                                                ResistanceNetwork or subclass.
        verbose (bool): Flag to control verbosity of output.
        task_queue (Queue): Queue for tasks.
        result_queue (Queue): Queue for results.
        workers (list): List of worker processes.
        bests (HallOfFame): Hall of fame for best models.
        toolbox (deap.base.Toolbox): Toolbox for genetic algorithm.
        logger (list): List to store logging information.
        Various other GA parameters.
    """

    def __init__(self, resistance_network, seed=1234, verbose=True):
        """
        Initializes a ModelRunner with a given resistance network, seed, and
        verbosity.

        Args:
            resistance_network (ResistanceNetwork or subclass): The resistance
                                                          network to optimize.
            seed (int): Seed for random number generator.
            verbose (bool): Flag to control verbosity of output.
        """
        self.seed = seed
        self.resistance_network = resistance_network
        self.verbose = verbose

        self.task_queue = Queue()
        self.result_queue = Queue()
        self.workers = []
        self.bests = None  # model hall of fame
        self.toolbox = None
        self.logger = list()

        self.cxpb = None
        self.mutpb = None
        self.indpb = None
        self.popsize = None
        self.maxpopsize = None
        self.posWeight = None
        self.fixWeight = None
        self.fixShape = None
        self.allShapes = None
        self.fixAsym = True
        self.min_weight = None
        self.max_shape = None
        self.max_hof_size = None
        self.tournsize = None
        self.fitmetric = None
        self.awsum = None
        self.only_keep = None
        self.report_all = None

    def initialize_ga(self):
        """
        Initializes the genetic algorithm components and population.
        """
        self.toolbox = base.Toolbox()

        # Initialize the DEAP GA components
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)

        if self.verbose:
            print("Initializing genetic algorithm parameters...\n")
        self.init_ga_attributes()

        # Initialize population
        popsize = len(self.resistance_network.variables) * 5 * 15
        if self.popsize:
            popsize = self.popsize
        if popsize > self.maxpopsize:
            popsize = self.maxpopsize

        if self.verbose:
            print("Establishing a population of size:", str(popsize))
        self.population = self.toolbox.population(n=popsize)

    def set_ga_parameters(self, mutpb, cxpb, indpb, popsize, maxpopsize,
                          posWeight, fixWeight, fixShape, allShapes,
                          min_weight, max_shape, max_hof_size, tournsize,
                          fitmetric, awsum, only_keep, verbose, report_all):
        """
        Sets the genetic algorithm parameters.

        Args:
            mutpb (float): Probability of mutation per trait.
            cxpb (float): Probability of being chosen for crossover.
            indpb (float): Probability of mutation per individual.
            popsize (int): Population size.
            maxpopsize (int): Maximum population size.
            posWeight (bool): Constrain parameter weights to between 0.0-1.0.
            fixWeight (bool): Constrain parameter weights to 1.0 (i.e.,
                              unweighted).
            fixShape (bool): Turn off feature transformation.
            allShapes (bool): Allow inverse and reverse transformations.
            min_weight (float): Minimum allowable weight.
            max_shape (float): Maximum shape value.
            max_hof_size (int): Maximum Hall of Fame size.
            tournsize (int): Tournament size.
            fitmetric (str): Fitness metric to be used.
            awsum (float): Cumulative Akaike weight threshold to retain top N
                           models.
            only_keep (bool): Only retain models where column "keep"=True.
            verbose (bool): Flag to control verbosity of output.
            report_all (bool): Flag to generate full outputs for all retained
                               models.
        """
        self.mutpb = mutpb
        self.cxpb = cxpb
        self.indpb = indpb
        self.popsize = popsize
        self.maxpopsize = maxpopsize
        self.posWeight = posWeight
        self.fixWeight = fixWeight
        self.fixShape = fixShape
        self.allShapes = allShapes
        self.min_weight = min_weight
        self.max_shape = max_shape
        self.max_hof_size = max_hof_size
        self.tournsize = tournsize
        self.fitmetric = fitmetric
        self.awsum = awsum
        self.only_keep = only_keep
        self.verbose = verbose
        self.report_all = report_all

    def run_ga(self, maxgens=1, fitmetric="aic", burnin=0, deltaB=None,
               deltaB_perc=0.001, indpb=0.5, mutpb=0.5, cxpb=0.5, nFail=50,
               popsize=None, maxpopsize=1000, posWeight=False, fixWeight=False,
               fixShape=False, allShapes=False, min_weight=0.0, max_shape=None,
               max_hof_size=100, tournsize=10, awsum=0.95, only_keep=True,
               out=None, plot=True, verbose=True, report_all=False, threads=1):
        """
        Runs the genetic algorithm for optimizing the ResistNet model.

        Args:
            maxgens (int): Maximum number of generations.
            fitmetric (str): Fitness metric to be used.
            burnin (int): Number of initial generations to ignore for
                          convergence check.
            deltaB (float): Absolute threshold change in fitness for
                            convergence.
            deltaB_perc (float): Relative threshold change in fitness for
                                 convergence.
            indpb (float): Probability of mutation per individual.
            mutpb (float): Probability of mutation per trait.
            cxpb (float): Probability of being chosen for crossover.
            nFail (int): Number of generations failing to improve to stop
                         optimization.
            popsize (int): Population size.
            maxpopsize (int): Maximum population size.
            posWeight (bool): Constrain parameter weights to between 0.0-1.0.
            fixWeight (bool): Constrain parameter weights to 1.0 (i.e.,
                              unweighted).
            fixShape (bool): Turn off feature transformation.
            allShapes (bool): Allow inverse and reverse transformations.
            min_weight (float): Minimum allowable weight.
            max_shape (float): Maximum shape value.
            max_hof_size (int): Maximum Hall of Fame size.
            tournsize (int): Tournament size.
            awsum (float): Cumulative Akaike weight threshold to retain top N
                           models.
            only_keep (bool): Only retain models where column "keep"=True.
            out (str): Output file prefix.
            plot (bool): Flag to control plotting.
            verbose (bool): Flag to control verbosity of output.
            report_all (bool): Generate full outputs for all retained models.
            threads (int): Number of parallel processors.

        Raises:
            Exception: Any exception encountered during the genetic algorithm
                       run.
        """
        try:
            # Set GA parameters
            self.set_ga_parameters(mutpb, cxpb, indpb, popsize, maxpopsize,
                                   posWeight, fixWeight, fixShape, allShapes,
                                   min_weight, max_shape, max_hof_size,
                                   tournsize, fitmetric, awsum, only_keep,
                                   verbose, report_all)

            # Initialize GA components and worker pool
            self.initialize_ga()
            self.start_workers(threads)

            # Initialize population and set up GA run parameters
            self.initialise_population()
            cxpb, indpb = self.cxpb, self.indpb
            fails, current_best = 0, None

            # Run for maxgens generations
            for g in range(1, maxgens + 1):
                if self.verbose:
                    print("-- Generation %i --" % g)

                # Select the next generation individuals
                offspring = self.toolbox.select(
                    self.population, len(self.population)
                )
                # Clone the selected individuals
                offspring = list(map(self.toolbox.clone, offspring))

                # Apply crossover and mutation on the offspring
                for child1, child2 in zip(offspring[::2], offspring[1::2]):
                    if random.random() < self.cxpb:
                        self.toolbox.mate(child1, child2)
                        del child1.fitness.values
                        del child2.fitness.values

                for mutant in offspring:
                    if random.random() < self.indpb:
                        self.toolbox.mutate(mutant)
                        del mutant.fitness.values

                # Evaluate the individuals with an invalid fitness
                invalid_ind = [
                    ind for ind in offspring if not ind.fitness.valid
                ]
                for ind_index, ind in enumerate(invalid_ind):
                    self.task_queue.put(('evaluate', ind_index, list(ind)))

                # Collect fitness results and process them
                pop_list = []
                for _ in range(len(invalid_ind)):
                    ind_index, fitness, res = self.result_queue.get()
                    invalid_ind[ind_index].fitness.values = [fitness]
                    # Process and append to pop_list
                    ind_list = [fitness] + list(invalid_ind[ind_index])
                    if res is not None:
                        ind_list.extend(res)
                    pop_list.append(ind_list)

                # Replace population with offspring
                self.population[:] = offspring

                # update Hall of Fame
                self.bests.check_population(pop_list)

                # Gather all the fitnesses in one list and print the stats
                fits = [ind.fitness.values[0] for ind in self.population]

                # Evaluate for stopping criteria after burn-in period
                if g > burnin:
                    if fitmetric.upper() == "AIC":
                        fits = [-element for element in fits]
                        worst = max(fits)
                        best = min(fits)
                        if current_best is None:
                            current_best = best
                        else:
                            current_best, fails = self.update_fails(
                                best, current_best, fails, deltaB, deltaB_perc,
                                minimize=True)
                    else:
                        worst = min(fits)
                        best = max(fits)
                        if current_best is None:
                            current_best = best
                        else:
                            current_best, fails = self.update_fails(
                                best, current_best, fails, deltaB,
                                deltaB_perc, minimize=False
                            )

                    length = len(self.population)
                    if length > 0:
                        if any(math.isinf(fit) for fit in fits):
                            # Set stats to NaN if any values are inf or -inf
                            mean = variance = std = float('nan')
                        else:
                            mean = sum(fits) / length
                            sum2 = sum(x * x for x in fits)
                            variance = sum2 / length - mean**2

                            std = (math.sqrt(variance)
                                   if variance >= 0 else float('nan'))
                    else:
                        mean = variance = std = float('nan')

                    if self.verbose:
                        print("  Worst %s" % worst)
                        print("  Best %s" % best)
                        print("  Avg %s" % mean)
                        print("  Std %s" % std)
                        print("  nFails %s" % fails)
                    self.logger.append([g, worst, best, mean, std])

                    # Check for stopping conditions
                    if g >= maxgens or fails > nFail:
                        if self.verbose:
                            print(
                                "Stopping optimization after",
                                g, "generations."
                            )
                            if g >= maxgens:
                                print("Reason: Exceeded maxGens")
                            elif fails > nFail:
                                print(
                                    "Reason: More than",
                                    nFail,
                                    "generations failed to improve"
                                )
                        break

            # output results
            self.run_output(out, verbose, plot)

        except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()

        finally:
            # Terminate worker processes
            self.terminate_workers()

    def build_ensemble(self, bests, awsum=0.95, only_keep=True, out=None,
                       threads=1, verbose=True):
        """
        Builds the ensemble model from the best models.

        Args:
            bests: Best models.
            awsum (float): The Akaike weight sum threshold.
            only_keep (bool): Flag to keep only certain models.
            out: Output parameter.
            threads (int): Number of parallel threads.
            verbose (bool): Verbose output flag.
        """
        try:
            self.awsum = awsum
            self.only_keep = only_keep
            self.verbose = verbose
            self.bests = bests
            self.report_all = False

            # Initialise workers
            self.start_workers(threads)

            # Compute ensemble model
            self.run_output(out, verbose, plot=True)

        except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()

        finally:
            # Terminate processes
            self.terminate_workers()

    def run_output(self, out=None, verbose=True, plot=True):
        """
        Runs the output routine for the model.

        Args:
            out: Output parameter.
            verbose (bool): Verbose output flag.
            plot (bool): Flag to enable or disable plotting.
        """
        # Various calculations and print/write operations
        self.bests.delta_aic()
        self.bests.akaike_weights()
        self.bests.cumulative_akaike(threshold=self.awsum)
        self.bests.relative_variable_importance(self.only_keep)
        self.bests.model_average_weights()

        # Printing and writing the results
        if verbose:
            self.bests.printHOF()
            self.bests.printRVI()
            self.bests.printMAW()

        if out:
            self.bests.writeMAW(out)
            self.bests.writeRVI(out)
            # Plotting, if enabled
            if plot:
                self.bests.plot_ICprofile(out)
                self.bests.plotMetricPW(out)
                self.bests.plotVariableImportance(out)

            # Writing Hall of Fame to file
            self.bests.writeModelSummary(out)

        # Writing log of fitnesses
        if out and len(self.logger) > 0:
            logDF = pd.DataFrame(self.logger,
                                 columns=[
                                     "Generation", "Worst", "Best", "Mean",
                                     "Stdev"
                                     ]
                                 )
            logDF.to_csv(
                f"{out}.FitnessLog.tsv", sep="\t", header=True, index=False
            )

        # Get results for best models and model-averaged
        self.model_average(out, plot)

        # Evaluate the null model (if genetic distances present)
        if self.resistance_network._gendist is not None:
            null = self.resistance_network.evaluate_null_model(out)
            if verbose:
                print("\nEvaluating null model...")
                with pd.option_context(
                    'display.max_rows', None, 'display.max_columns', None
                ):
                    print(null)

    def model_average(self, out="", plot=True):
        """
        Performs model averaging on a set of resistance models.

        Args:
            out (str): Output file base name.
            plot (bool): Flag to enable or disable plotting.
        """
        # Build model list and run Circuitscape on each model
        models = []
        mod_num = 0
        weights = {}
        hof = self.bests.getHOF(only_keep=self.only_keep)
        for _, row in hof.iterrows():
            n = len(self.resistance_network._predictors.columns) * 5
            model_string = row.iloc[1:n + 1].to_list()
            models.append([mod_num, model_string])
            weights[mod_num] = row["akaike_weight"]
            mod_num += 1

        if self.verbose:
            print("\nModel-averaging across", mod_num, "resistance models...")

        # Evaluate models in parallel
        for model in models:
            self.task_queue.put(('output', model[0], model[1]))

        # Retrieve model results
        model_results = []
        for _ in range(len(models)):
            mod_num, result, multi = self.result_queue.get()
            model_results.append([mod_num, result, multi])

        # Model averaging with Akaike weights
        edge_avg = None
        matrix_avg = np.zeros(
            shape=(len(self.resistance_network._points_labels),
                   len(self.resistance_network._points_labels))
            )

        for m in model_results:
            model_num = m[0]
            matrix_r = m[1]
            edge_r = m[2]

            # Get weighted contributions to average model
            weight = hof.iloc[model_num]["akaike_weight"]
            if edge_avg is None:
                edge_avg = edge_r * weight
            else:
                edge_avg = np.add(edge_avg, edge_r * weight)
            matrix_avg += matrix_r * weight

            if self.report_all:
                oname = f"{out}.Model-{model_num}"
                self.resistance_network.output_and_plot_model(
                    oname, matrix_r, edge_r
                )

        oname = f"{out}.Model-Average"
        self.resistance_network.output_and_plot_model(
            oname, matrix_avg, edge_avg
        )

    @staticmethod
    def update_fails(best, current_best, fails, delta_b, delta_b_perc,
                     minimize=False):
        """
        Updates the fail count based on the best and current best fitness
        scores.

        Args:
            best: The best fitness score.
            current_best: The current best fitness score.
            fails: The current fail count.
            delta_b: Absolute improvement threshold.
            delta_b_perc: Percentage improvement threshold.
            minimize (bool): Flag indicating minimization or maximization.

        Returns:
            Tuple containing updated current best fitness score and fail count.
        """
        cur = current_best
        f = fails
        threshold_1 = current_best
        threshold_2 = current_best

        if minimize:
            if best < current_best:
                cur = best
            if delta_b is not None:
                threshold_1 = current_best - delta_b
            if delta_b_perc is not None:
                threshold_2 = current_best - (current_best * delta_b_perc)
            if best >= threshold_1 or best >= threshold_2:
                f += 1
            else:
                f = 0
        else:
            if best > current_best:
                cur = best
            if delta_b is not None:
                threshold_1 = current_best + delta_b
            if delta_b_perc is not None:
                threshold_2 = current_best + (current_best * delta_b_perc)
            if best <= threshold_1 or best <= threshold_2:
                f += 1
            else:
                f = 0
        return cur, f

    def initialise_population(self):
        """
        Initialises the population by distributing tasks for evaluation and
        processing the fitness results.
        """
        # Distribute tasks for evaluation
        for ind_index, ind in enumerate(self.population):
            self.task_queue.put(('evaluate', ind_index, list(ind)))

        # Collect fitness results and process them
        pop_list = []
        for _ in range(len(self.population)):
            ind_index, fitness, res = self.result_queue.get()
            self.population[ind_index].fitness.values = [fitness]
            ind_list = [fitness] + list(self.population[ind_index])
            if res is not None:
                ind_list.extend(res)
            pop_list.append(ind_list)

        # Initialize Hall of Fame
        self.bests = HallOfFame(
            self.resistance_network._predictors.columns,
            self.max_hof_size,
            pop_list)

    def init_ga_attributes(self):
        """
        Initializes the genetic algorithm's attributes and registers necessary
        functions with the toolbox.
        """
        # Register attributes
        self.toolbox.register("feature_sel", random.randint, 0, 1)

        # weight
        if self.posWeight:
            self.toolbox.register(
                "feature_weight", random.uniform, self.min_weight, 1.0
            )
        if self.fixWeight:
            self.toolbox.register("feature_weight", random.uniform, 1.0, 1.0)
        if not self.fixWeight and not self.posWeight:
            self.toolbox.register("feature_weight", random.uniform, -1.0, 1.0)

        # asymmetry parameter
        if not self.fixAsym:
            self.toolbox.register("feature_asym", random.randint, 0, 1)
        else:
            self.toolbox.register("feature_asym", random.randint, 0, 0)
        
        # shape
        if not self.fixShape:
            self.toolbox.register("feature_transform", random.randint, 0, 8)
            self.toolbox.register(
                "feature_shape", random.randint, 1, self.max_shape
            )
        else:
            self.toolbox.register("feature_transform", random.randint, 0, 0)
            self.toolbox.register("feature_shape", random.randint, 1, 1)

        self.toolbox.register(
            "individual",
            tools.initCycle,
            creator.Individual,
            (self.toolbox.feature_sel,
                self.toolbox.feature_weight,
                self.toolbox.feature_transform,
                self.toolbox.feature_shape,
                self.toolbox.feature_asym),
            n=len(self.resistance_network.variables)
        )
        self.toolbox.register(
            "population", tools.initRepeat, list, self.toolbox.individual
        )
        self.toolbox.register("evaluate", lambda ind: ind)
        self.toolbox.register("mate", tools.cxTwoPoint)
        self.toolbox.register("mutate", self.mutate)
        self.toolbox.register(
            "select", tools.selTournament, tournsize=self.tournsize
        )

    def terminate_workers(self):
        """
        Terminates all worker processes.
        """
        # Send a termination signal to each worker
        for _ in self.workers:
            self.task_queue.put(("END", None, None))

        # Wait for all worker processes to complete
        for worker in self.workers:
            worker.join()

        if self.verbose:
            print("\nAll worker processes have been terminated.")

    def start_workers(self, threads):
        """
        Starts worker processes for parallel task execution.

        Args:
            threads (int): The number of worker threads to start.
        """
        if self.verbose:
            print()

        # get worker type
        # NOTE: THe order is important here, have to do subclass before
        # parent check (since both will return true)
        if isinstance(self.resistance_network, ResistanceNetworkSAMC):
            worker_type = "samc"
        elif isinstance(self.resistance_network, ResistanceNetwork):
            worker_type = "base"
        else:
            raise ValueError("Could not set worker type")

        for i in range(threads):
            worker_seed = self.seed + i
            if self.verbose:
                print(f"Starting worker {i} with seed {worker_seed}")

            # Pass only necessary and picklable information
            worker_args = {
                'worker_type': worker_type,
                'network': self.resistance_network.network,
                'pop_agg': self.resistance_network.pop_agg,
                'reachid_col': self.resistance_network.reachid_col,
                'length_col': self.resistance_network.length_col,
                'variables': self.resistance_network.variables,
                'agg_opts': self.resistance_network.agg_opts,
                'inc': self.resistance_network._inc,
                'point_coords': self.resistance_network._point_coords,
                'points_names': self.resistance_network._points_names,
                'points_snapped': self.resistance_network._points_snapped,
                'points_labels': self.resistance_network._points_labels,
                'predictors': self.resistance_network._predictors,
                'edge_order': self.resistance_network._edge_order,
                'gendist': self.resistance_network._gendist,
                "fitmetric": self.fitmetric,
                "posWeight": self.posWeight,
                "fixWeight": self.fixWeight,
                # "fixAsym": self.fixAsym,
                "allShapes": self.allShapes,
                "fixShape": self.fixShape,
                "min_weight": self.min_weight,
                "max_shape": self.max_shape
            }
            if worker_type == "samc":
                worker_args['adj'] = self.resistance_network._adj
                worker_args['origin'] = self.resistance_network._origin
                worker_args['R'] = self.resistance_network._R
                worker_args[
                    'edge_site_indices'
                    ] = self.resistance_network._edge_site_indices

            worker_process = Process(
                target=self.worker_task,
                args=(self.task_queue, self.result_queue, worker_args,
                      worker_seed)
            )
            worker_process.start()
            self.workers.append(worker_process)
        if self.verbose:
            print()

    @staticmethod
    def worker_task(task_queue, result_queue, worker_args, seed):
        """
        Static method to handle worker tasks for parallel processing.

        Args:
            task_queue (Queue): The queue from which tasks are retrieved.
            result_queue (Queue): The queue to which results are put.
            worker_args (dict): Arguments required for creating a worker
                                instance.
            seed (int): The seed value for random number generation.
        """
        random.seed(seed)

        worker_type = worker_args.get('worker_type', 'base')
        del worker_args['worker_type']
        if worker_type == "base":
            worker = ResistanceNetworkWorker(**worker_args)
        elif worker_type == "samc":
            worker = ResistanceNetworkSAMCWorker(**worker_args)
        else:
            raise ValueError("Worker type not implemented")

        while True:
            try:
                task, id, data = task_queue.get()
                if task == "END":
                    break
                elif task == "evaluate":
                    fitness, res = worker.evaluate(data)
                    if res:
                        result_queue.put((id, fitness, res[0:4]))
                    else:
                        result_queue.put((id, fitness, None))
                elif task == "output":
                    r, multi = worker.model_output(data)
                    result_queue.put((id, r, multi))
            except Empty:
                continue

    def mutate(self, individual):
        """
        Custom mutation function for an individual.

        Args:
            individual: An individual to be mutated.

        Returns:
            A tuple containing the mutated individual.
        """
        length = len(individual) // 5  # number of covariates
        for i in range(length):
            if random.random() < self.mutpb:
                individual[i * 5] = self.toolbox.feature_sel()
            if random.random() < self.mutpb:
                individual[i * 5 + 1] = self.toolbox.feature_weight()
            if random.random() < self.mutpb:
                individual[i * 5 + 2] = self.toolbox.feature_transform()
            if random.random() < self.mutpb:
                individual[i * 5 + 3] = self.toolbox.feature_shape()
            if random.random() < self.mutpb:
                individual[i * 5 + 4] = self.toolbox.feature_asym()
        return individual,
