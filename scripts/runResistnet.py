#!/usr/bin/env python3
import os
import random
from datetime import datetime

from resistnet.params import parseArgs
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunnerTPE


def main():
    """
    Main function to run the ResistNet application.
    """
    # Step 0: Read command-line arguments
    params = parseArgs()

    # Seed random number generator
    if not params.seed:
        params.seed = int(datetime.now().timestamp())
    random.seed(params.seed)

    # Step 1: Initialise network
    network = ResistanceNetwork(
        network=params.network,
        shapefile=params.shapefile,
        coords=params.coords,
        variables=params.variables,
        agg_opts=params.agg_opts,
        inmat=params.inmat,
        reachid_col=params.reachid_col,
        length_col=params.length_col,
        out=params.out,
        verbose=True
    )

    # write new output geodatabase
    if params.shapefile:
        network.write_geodataframe(params.out, params.output_driver)

    runner = ModelRunnerTPE(
        resistance_network=network,
        seed=params.seed,
        verbose=True
    )

    # Optionally optimise transformations for each parameter
    if params.gridSearch:
        runner.optimise_univariate(
            fitmetric=params.fitmetric,
            threads=params.procs,
            max_shape=params.max_shape,
            out=params.out,
            plot=True,
            verbose=True
        )

    # Step 3: Run optimisation
    runner.run_tpe(
        max_evals=params.max_iter,
        fitmetric=params.fitmetric,
        threads=params.procs,
        nFail=params.nfail,
        fixWeight=params.fixWeight,
        fixShape=params.fixShape,
        min_weight=params.min_weight,
        max_shape=params.max_shape,
        max_hof_size=params.max_hof_size,
        awsum=params.awsum,
        only_keep=True,
        use_full=params.use_full,
        fixed_params=params.transFile,
        out=params.out,
        plot=True,
        verbose=True,
        pweight=params.pweight,
        reps=params.reps,
        n_startup=params.nstart,
        n_candidates=params.ncand,
        gamma=params.gamma
    )

    # # Step 2: Initialise ModelRunner
    # runner = ModelRunnerGA(
    #     resistance_network=network,
    #     seed=params.seed,
    #     verbose=True
    # )

    # # Step 3: Run GA optimisation
    # runner.run_ga(
    #     maxgens=params.maxGens,
    #     fitmetric=params.fitmetric,
    #     threads=params.GA_procs,
    #     burnin=params.burnin,
    #     deltaB=params.deltaB,
    #     deltaB_perc=params.deltaB_perc,
    #     indpb=params.indpb,
    #     mutpb=params.mutpb,
    #     cxpb=params.cxpb,
    #     nFail=params.nfail,
    #     popsize=params.popsize,
    #     maxpopsize=params.maxpopsize,
    #     fixWeight=params.fixWeight,
    #     fixShape=params.fixShape,
    #     min_weight=params.min_weight,
    #     max_shape=params.max_shape,
    #     max_hof_size=params.max_hof_size,
    #     awsum=params.awsum,
    #     only_keep=params.only_keep,
    #     use_full=params.use_full,
    #     fixed_params=params.transFile,
    #     out=params.out,
    #     plot=True,
    #     verbose=True
    # )


# Call main function
if __name__ == '__main__':
    os.environ['USE_PYGEOS'] = '0'
    main()
