#!/usr/bin/env python3
import os
import sys
import random
from datetime import datetime

from resistnet.params import parseArgs
from resistnet.resistance_network import ResistanceNetwork
from resistnet.samc_network import ResistanceNetworkSAMC
from resistnet.model_optimisation import ModelRunner


def main():
    """
    Main function to run the ResistNet application.
    """
    # Step 0: Read command-line arguments
    params = parseArgs()

    # Seed random number generator
    if not params.seed:
        params.seed = datetime.now().timestamp()
    random.seed(params.seed)

    # Step 1: Initialise network data
    network = ResistanceNetworkSAMC(
        network=params.network,
        shapefile=params.shapefile,
        sizes=params.sizefile,
        coords=params.coords,
        variables=params.variables,
        agg_opts=params.agg_opts,
        pop_agg=params.pop_agg,
        inmat=params.inmat,
        reachid_col=params.reachid_col,
        length_col=params.length_col,
        infer_origin=params.infer_origin,
        origin=params.origin,
        rtol=params.rtol,
        solver=params.solver,
        max_iter=params.max_iter,
        max_fail=params.max_fail,
        out=params.out,
        verbose=True
    )
    # network = ResistanceNetwork(
    #     network=params.network,
    #     shapefile=params.shapefile,
    #     coords=params.coords,
    #     variables=params.variables,
    #     agg_opts=params.agg_opts,
    #     pop_agg=params.pop_agg,
    #     inmat=params.inmat,
    #     reachid_col=params.reachid_col,
    #     length_col=params.length_col,
    #     out=params.out,
    #     verbose=True
    # )

    # write new output geodatabase 
    network.write_geodataframe(params.out, params.output_driver)

    # Step 2: Initialise ModelRunner
    runner = ModelRunner(
        resistance_network=network,
        seed=params.seed,
        verbose=True
    )

    # Step 3: Run GA optimisation
    runner.run_ga(
        maxgens=params.maxGens,
        fitmetric=params.fitmetric,
        threads=params.GA_procs,
        burnin=params.burnin,
        deltaB=params.deltaB,
        deltaB_perc=params.deltaB_perc,
        indpb=params.indpb,
        mutpb=params.mutpb,
        cxpb=params.cxpb,
        nFail=params.nfail,
        popsize=params.popsize,
        maxpopsize=params.maxpopsize,
        posWeight=params.posWeight,
        fixWeight=params.fixWeight,
        fixShape=params.fixShape,
        allShapes=params.allShapes,
        min_weight=params.min_weight,
        max_shape=params.max_shape,
        max_hof_size=params.max_hof_size,
        awsum=params.awsum,
        only_keep=params.only_keep,
        out=params.out,
        plot=True,
        verbose=True
    )


# Call main function
if __name__ == '__main__':
    os.environ['USE_PYGEOS'] = '0'
    main()
