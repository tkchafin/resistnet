import traceback
import random
import pandas as pd
import numpy as np
from queue import Empty
from hyperopt import tpe, Trials, hp
from hyperopt.base import Domain
from multiprocessing import Process, Queue
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages

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
        self.logger = list()

        self.awsum = None
        self.only_keep = None
        self.use_full = None
        self.report_all = None
        self.fitmetric = None
        self.fixed_params = None

        self.fixWeight = None
        self.fixShape = None
        self.fixAsym = False
        self.min_weight = None
        self.max_shape = None
        self.max_hof_size = None

        self.num_attributes = 4

        random.seed(self.seed)
        np.random.seed(self.seed)

    def build_ensemble(self, bests, awsum=0.95, only_keep=True, use_full=False,
                       out=None, threads=1, verbose=True):
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
            self.use_full = use_full
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
        self.bests.use_full = self.use_full
        self.bests.delta_aic()
        self.bests.akaike_weights()
        self.bests.cumulative_akaike(threshold=self.awsum)
        self.bests.relative_variable_importance(not self.only_keep)
        self.bests.model_average_weights(not self.only_keep)

        # Printing and writing the results
        if verbose:
            self.bests.printHOF()
            self.bests.printRVI()
            self.bests.printMAW()

        if out:
            self.bests.writeMAW(out)
            self.bests.writeRVI(out)
            self.bests.writeBest(out)
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
            if matrix_r is not None:
                matrix_avg += matrix_r * weight

            if model_num == 0:
                oname = f"{out}.Model-Best"
                self.resistance_network.output_and_plot_model(
                    oname, matrix_r, edge_r
                )

            if self.report_all:
                oname = f"{out}.Model-{model_num}"
                self.resistance_network.output_and_plot_model(
                    oname, matrix_r, edge_r
                )
        oname = f"{out}.Model-Average"
        self.resistance_network.output_and_plot_model(
            oname, matrix_avg, edge_avg
        )

    def parse_fixed_params(self, fixed_params):
        if isinstance(fixed_params, dict):
            self.fixed_params = fixed_params
        else:
            try:
                df = pd.read_csv(fixed_params, sep='\t')
                # Check if the DataFrame has a fitness column
                if 'fitness' in df.columns:
                    best_params = dict()
                    for var in df["variable"].unique():
                        var_df = df[df['variable'] == var]
                        best = var_df.loc[
                            var_df['fitness'] == var_df['fitness'].min()]
                        best = best.sort_values(
                            by=['shape', 'transform']).iloc[0]
                        best_params[var] = {
                            'transform': best['transform'],
                            'shape': best['shape']
                        }
                        self.fixed_params = best_params
                else:
                    # Assume file contains correct columns
                    self.fixed_params = df.set_index(
                        'variable').to_dict('index')
            except Exception as e:
                print(f"Failed to read fixed parameters: {e}")
                raise

    def write_fixed_params(self, out):
        data = []
        for var_name, params in self.fixed_params.items():
            shape = params.get("shape", np.nan)
            transform = params.get("transform", np.nan)
            data.append({
                "variable": var_name,
                "shape": shape,
                "transform": transform
            })
        df = pd.DataFrame(data)
        df.to_csv(out, index=False, sep="\t")

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
            print("\nAll worker processes have been terminated\n")

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
                "fixWeight": self.fixWeight,
                # "fixAsym": self.fixAsym,
                "fixShape": self.fixShape,
                "min_weight": self.min_weight,
                "max_shape": self.max_shape
            }
            if worker_type == "samc":
                worker_args['adj'] = self.resistance_network._adj
                worker_args['origin'] = self.resistance_network._origin
                worker_args['R'] = self.resistance_network._R
                worker_args['rtol'] = self.resistance_network.rtol
                worker_args['solver'] = self.resistance_network.solver
                worker_args['max_iter'] = self.resistance_network.max_iter
                worker_args['max_fail'] = self.resistance_network.max_fail
                worker_args[
                    'allSymmetric'
                    ] = self.resistance_network.allSymmetric
                worker_args[
                    'edge_site_indices'
                    ] = self.resistance_network._edge_site_indices
                # worker_args[
                #     'root_edge_indices'
                #     ] = self.resistance_network._root_edge_indices

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
        np.random.seed(seed)

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

    def optimise_univariate(self, fitmetric="aic", threads=1,
                            min_weight=0.0, max_shape=None,
                            out=None, plot=True, verbose=True):
        """
        Performs a grid search over transformations and shapes for each
        variable to find the optimal settings, using the given configuration
        parameters.

        Args:
            threads (int): Number of parallel processes to use.
            max_shape (int): Maximum shape value to consider in the grid search
            verbose (bool): Enable verbose output during the process.
            fitmetric (str): Fitness metric to be used for evaluating models.
            out (str or None): If provided, saves the detailed results to a
                               DataFrame with optional PDF plots

        Returns:
            dict: A dictionary containing the best transform and shape
                parameters for each variable.
        """

        try:
            self.fitmetric = fitmetric
            self.min_weight = min_weight
            self.max_shape = max_shape
            self.verbose = verbose

            if self.verbose:
                print("Starting univariate optimisation...")

            # Initialise workers
            self.start_workers(threads)

            best_params = dict()
            results_data = []

            # Define the range for transformations and shapes
            trans_range = range(0, 9)
            shape_range = range(1, max_shape+1) if max_shape else range(1, 101)

            # Iterate over each variable individually
            for var in self.resistance_network.variables:

                if self.verbose:
                    print("Performing grid search for", var)

                # Create models and distribute tasks
                for t in trans_range:
                    for s in shape_range:
                        weight = 1.0
                        model = [0] * (
                            len(self.resistance_network.variables) *
                            self.num_attributes)
                        index_var = self.resistance_network.variables.index(
                            var)
                        model[
                            index_var * self.num_attributes:index_var *
                            self.num_attributes + self.num_attributes
                        ] = [1, weight, t, s, 0]
                        model_index = f"{var};{t};{s}"
                        self.task_queue.put(('evaluate', model_index, model))

                # Collect fitness results and update best parameters
                for _ in range(len(trans_range) * len(shape_range)):
                    model_index, fitness, _ = self.result_queue.get()
                    var, t, s = model_index.split(';')
                    t, s = int(t), int(s)
                    results_data.append({
                        'variable': var,
                        'fitness': fitness,
                        'shape': s,
                        'transform': t
                    })

            # choose best parameter values, favouring lowest transformation
            df = pd.DataFrame(results_data)
            for var in self.resistance_network.variables:
                var_df = df[df['variable'] == var]
                best = var_df.loc[var_df['fitness'] == var_df['fitness'].max()]
                best = best.sort_values(by=['shape', 'transform']).iloc[0]
                best_params[var] = {
                    'transform': best['transform'],
                    'shape': best['shape']
                }

            if out:
                df.to_csv(
                    str(out)+".univariateFitness.tsv", index=False, sep="\t")

                # optional plots
                if plot:
                    # Define the PDF file path
                    pdf_path = str(out) + ".univariateFitness.pdf"
                    with PdfPages(pdf_path) as pdf:
                        for var in self.resistance_network.variables:
                            # Filter data for the current variable
                            var_data = df[df['variable'] == var]
                            # Create a pivot table for the heatmap
                            pivot_table = var_data.pivot(
                                index="shape",
                                columns="transform",
                                values="fitness"
                            )

                            fig, ax = plt.subplots(figsize=(10, 8))
                            _ = sns.heatmap(
                                pivot_table, ax=ax, annot=True, fmt=".1f",
                                cmap='viridis')
                            plt.title(f'{var}')
                            plt.xlabel('Transformation')
                            plt.ylabel('Shape')

                            # Highlight the best parameter combination
                            best_transform = best_params[var]['transform']
                            best_shape = best_params[var]['shape']
                            transform_pos = pivot_table.columns.tolist().index(
                                best_transform)
                            shape_pos = pivot_table.index.tolist().index(
                                best_shape)
                            rect = patches.Rectangle(
                                (transform_pos, shape_pos), 1, 1, linewidth=2,
                                edgecolor='red', facecolor='none')
                            ax.add_patch(rect)

                            pdf.savefig(fig)
                            plt.close(fig)

            self.fixed_params = best_params
            return best_params

        except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()

        finally:
            # Terminate processes
            self.terminate_workers()
            self.seed = random.randint(1, 10000000)


class ModelRunnerTPE(ModelRunner):
    """
    Class to handle running TPE optimization for a ResistNet model.

    Attributes:
        seed (int): Seed for random number generator.
        resistance_network (ResistanceNetwork): An instance of
                                                ResistanceNetwork.
        verbose (bool): Flag to control verbosity of output.
        task_queue (Queue): Queue for tasks.
        result_queue (Queue): Queue for results.
        workers (list): List of worker processes.
        bests (HallOfFame): Hall of fame for best models.
        Various other parameters for model configuration.
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
        super().__init__(resistance_network, seed, verbose)
        self.toolbox = None

        self.max_evals = 100
        self.num_workers = 10
        self.space = None  # This will be defined later

    def initialize_space(self):
        """Define the search space for hyperparameters."""
        self.space = {}
        for var_name in self.resistance_network.variables:
            if self.fixed_params and var_name in self.fixed_params:
                fixed_param = self.fixed_params[var_name]
                self.space[f"{var_name}_sel"] = hp.choice(
                    f"{var_name}_sel",
                    [fixed_param.get('sel', 1)]
                )
                self.space[f"{var_name}_weight"] = hp.uniform(
                    f"{var_name}_weight",
                    fixed_param.get('weight', 1.0),
                    fixed_param.get('weight', 1.0)
                )
                self.space[f"{var_name}_transform"] = hp.choice(
                    f"{var_name}_transform",
                    [fixed_param.get('transform', 0)]
                )
                self.space[f"{var_name}_shape"] = hp.uniform(
                    f"{var_name}_shape",
                    fixed_param.get('shape', 1.0),
                    fixed_param.get('shape', 1.0)
                )
            else:
                self.space[f"{var_name}_sel"] = hp.choice(
                    f"{var_name}_sel",
                    [0, 1]
                )

                if self.fixWeight:
                    self.space[f"{var_name}_weight"] = hp.uniform(
                        f"{var_name}_weight",
                        1.0,
                        1.0
                    )
                elif self.min_weight:
                    self.space[f"{var_name}_weight"] = hp.uniform(
                        f"{var_name}_weight",
                        self.min_weight,
                        1.0
                    )
                else:
                    self.space[f"{var_name}_weight"] = hp.uniform(
                        f"{var_name}_weight",
                        0.0,
                        1.0
                    )

                if not self.fixShape:
                    self.space[f"{var_name}_transform"] = hp.choice(
                        f"{var_name}_transform",
                        list(range(0, 9))
                    )
                    self.space[f"{var_name}_shape"] = hp.quniform(
                        f"{var_name}_shape",
                        1,
                        self.max_shape,
                        1.0
                    )
                else:
                    self.space[f"{var_name}_transform"] = hp.choice(
                        f"{var_name}_transform",
                        [0]
                    )
                    self.space[f"{var_name}_shape"] = hp.choice(
                        f"{var_name}_shape",
                        [0]
                    )

    def set_tpe_parameters(self, max_evals=100, num_workers=10,
                           fitmetric="likelihood", fixWeight=None,
                           fixShape=None, fixAsym=False, min_weight=None,
                           max_shape=None, max_hof_size=None, only_keep=None,
                           use_full=False, verbose=True, report_all=False,
                           awsum=0.95, fixed_params=None, out=None):
        """
        Sets the TPE optimization parameters.

        Args:
            max_evals (int): Maximum number of evaluations.
            num_workers (int): Number of worker processes.
            fitmetric (str): Fitness metric to be used.
            fixWeight (bool): Constrain parameter weights to 1.0
            fixShape (bool): Turn off feature transformation.
            fixAsym (bool): Constrain asymmetry to 0.
            min_weight (float): Minimum allowable weight.
            max_shape (float): Maximum shape value.
            max_hof_size (int): Maximum Hall of Fame size.
            only_keep (bool): Only retain models where column "keep"=True.
            use_full (bool): Flag to use full dataset.
            verbose (bool): Flag to control verbosity of output.
            report_all (bool): Flag to generate full outputs for all models.
            awsum (float): Cumulative Akaike weight threshold to retain models.
            fixed_params (dict): Fixed parameters to narrow the search space.
        """
        self.max_evals = max_evals
        self.num_workers = num_workers
        self.fitmetric = fitmetric
        self.fixWeight = fixWeight
        self.fixShape = fixShape
        self.fixAsym = fixAsym
        self.min_weight = min_weight
        self.max_shape = max_shape
        self.max_hof_size = max_hof_size
        self.only_keep = only_keep
        self.use_full = use_full
        self.verbose = verbose
        self.report_all = report_all
        self.awsum = awsum
        self.fixed_params = fixed_params
        self.out = out

    def run_tpe(self, max_evals=100, threads=4, fitmetric="loglik",
                fixWeight=None, fixShape=None, fixAsym=False, min_weight=None,
                max_shape=None, max_hof_size=None, only_keep=None,
                use_full=False, verbose=True, report_all=False, nFail=50,
                awsum=0.95, fixed_params=None, out=None, plot=True, reps=10,
                n_startup=40, n_candidates=48, gamma=0.15, pweight=0.7):
        """
        Runs the TPE optimization for optimizing the ResistNet model.

        Args:
            max_evals (int): Maximum number of evaluations.
            threads (int): Number of worker processes.
            fitmetric (str): Fitness metric to be used.
            fixWeight (bool): Constrain parameter weights to 1.0.
            fixShape (bool): Turn off feature transformation.
            fixAsym (bool): Constrain asymmetry to 0.
            min_weight (float): Minimum allowable weight.
            max_shape (float): Maximum shape value.
            max_hof_size (int): Maximum Hall of Fame size.
            only_keep (bool): Only retain models where column "keep"=True.
            use_full (bool): Flag to use full dataset.
            verbose (bool): Flag to control verbosity of output.
            report_all (bool): Flag to generate full outputs for all retained
                               models.
            nFail (int): Number of iterations to stop if fail to improve.
            awsum (float): Cumulative Akaike weight threshold to retain top N
                           models.
            fixed_params (dict): Fixed parameters to narrow the search space.
        """
        try:
            # Set TPE parameters
            self.set_tpe_parameters(max_evals, threads, fitmetric, fixWeight,
                                    fixShape, fixAsym, min_weight, max_shape,
                                    max_hof_size, only_keep, use_full, verbose,
                                    report_all, awsum, fixed_params, out)

            # Parse optional input fixed_params
            if fixed_params:
                self.parse_fixed_params(fixed_params)

            # Write fixed params if provided
            if out and self.fixed_params:
                self.write_fixed_params(str(out) + ".fixedParams.tsv")

            # Define the search space
            self.initialize_space()

            # Create the domain for hyperopt
            domain = Domain(dummy_objective, self.space)

            # Initialize Trials objects
            trials_list = [Trials() for _ in range(threads)]

            # Start workers
            self.start_workers(threads)

            best_loss = float('inf')
            global_fails = 0

            for iteration in range(max_evals):
                # Check if the global failure threshold has been reached
                if global_fails >= nFail:
                    print(f"Stopping due to {nFail} consecutive failures")
                    break

                param_sets = []
                for trials in trials_list:
                    new_ids = trials.new_trial_ids(1)
                    seed = random.randint(1, 1000000000)
                    param_set = tpe.suggest(new_ids,
                                            domain,
                                            trials,
                                            seed,
                                            prior_weight=pweight,
                                            n_startup_jobs=n_startup,
                                            n_EI_candidates=n_candidates,
                                            gamma=gamma)
                    param_sets.append(param_set[0])

                inds = {}

                # Convert and submit parameter sets to the task queue
                for idx, params in enumerate(param_sets):
                    individual = self.convert_params_to_individual(params)
                    inds[idx] = individual
                    self.task_queue.put(('evaluate', idx, individual))

                # Retrieve results from the result queue
                results = []
                for idx, params in enumerate(param_sets):
                    result_id, fitness, res = self.result_queue.get()
                    results.append((result_id, fitness, res))

                # Track if the global best loss has improved
                global_improved = False

                # Process results and update Trials objects
                for result_id, fitness, res in results:
                    trial = trials_list[result_id]
                    trial.insert_trial_doc({
                        'state': 2,  # This indicates the trial has finished
                        'tid': len(trial.trials),
                        'spec': None,
                        'result': {
                            'loss': -fitness,
                            'status': 'ok',
                            'res': res
                        },
                        'misc': {
                            'tid': len(trial.trials),
                            'cmd': 'domain_attachment',
                            'idxs': param_sets[result_id]['misc']['idxs'],
                            'vals': param_sets[result_id]['misc']['vals'],
                        },
                        'exp_key': None,
                        'owner': None,
                        'version': 0,
                        'book_time': None,
                        'refresh_time': None,
                    })
                    trial.refresh()

                    if -fitness < best_loss:
                        best_loss = -fitness
                        global_improved = True

                # Update fails count
                if global_improved:
                    global_fails = 0
                else:
                    global_fails += 1

                if verbose:
                    print(
                        f"Iteration {iteration+1}, Current Best: {best_loss}"
                    )

            if verbose:
                for i, trials in enumerate(trials_list):
                    i_best = trials.best_trial['result']['loss']
                    print(f"Worker {i + 1}: Best = {i_best}")

            # Process Hall of Fame to select the best model from each replicate
            self.process_hall_of_fame(trials_list)

            # Output results
            self.run_output(out, verbose, plot)

        except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()

        finally:
            # Terminate worker processes
            self.terminate_workers()

    def convert_params_to_individual(self, params):
        """
        Converts the hyperparameter dictionary to the correct format for
        evaluation.

        Args:
            params (dict): Hyperparameters for the model.

        Returns:
            list: Formatted individual list.
        """
        num_attributes = self.num_attributes
        num_vars = len(self.resistance_network.variables)
        individual = [0] * (num_vars * num_attributes)

        # extract parameter settings
        model = params["misc"]["vals"]

        for var_index, var_name in enumerate(
            self.resistance_network.variables
        ):
            idx_base = var_index * num_attributes
            individual[idx_base] = int(model[f"{var_name}_sel"][0])
            individual[idx_base + 1] = float(model[f"{var_name}_weight"][0])
            individual[idx_base + 2] = int(model[f"{var_name}_transform"][0])
            individual[idx_base + 3] = int(
                round(model[f"{var_name}_shape"][0])
                )
        return individual

    def process_hall_of_fame(self, trials_list):
        """
        Processes the Hall of Fame to select the best model from each replicate

        Args:
            trials_list (list): List of Trials objects, one for each replicate.
        """
        try:
            best_models = []
            for i, trials in enumerate(trials_list):
                if trials.trials:
                    # Extract the best trial
                    best_trial = min(
                        trials.trials,
                        key=lambda trial: trial['result']['loss']
                    )
                    individual = self.convert_params_to_individual(best_trial)
                    fitness = -best_trial['result']['loss']
                    res = best_trial['result'].get('res', [])
                    best_model = [fitness] + individual + res
                    best_models.append(best_model)

            # Define the column names based on the provided example structure
            num_vars = len(self.resistance_network.variables)
            num_attributes = self.num_attributes
            individual_columns = [""] * (num_vars * num_attributes)
            for var_index, var_name in enumerate(
                self.resistance_network.variables
            ):
                idx_base = var_index * num_attributes
                individual_columns[idx_base] = f"{var_name}"
                individual_columns[idx_base + 1] = f"{var_name}_weight"
                individual_columns[idx_base + 2] = f"{var_name}_trans"
                individual_columns[idx_base + 3] = f"{var_name}_shape"

            res_columns = ['loglik', 'r2m', 'aic', 'delta_aic_null']
            columns = ['fitness'] + individual_columns + res_columns

            # Create DataFrame from the list of best models
            df = pd.DataFrame(best_models, columns=columns)

            # Create the Hall of Fame
            self.bests = HallOfFame.from_dataframe(df)

        except Exception as e:
            print(f"An error occurred while processing the Hall of Fame: {e}")
            traceback.print_exc()


def weight_generator(min_weight):
    # Decide at random whether the weight will be positive or negative
    def generate():
        if random.random() < 0.5:
            # Generate a negative weight
            return random.uniform(-1.0, -min_weight)
        else:
            # Generate a positive weight
            return random.uniform(min_weight, 1.0)
    return generate


def dummy_objective(params):
    return 0
