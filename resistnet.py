import sys
import os

os.environ['USE_PYGEOS'] = '0'

from datetime import datetime
import random

from resistnet.params import parseArgs
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunner

def main():

    #########################################################
    # Step 0: Read command-line arguments
    #########################################################

    params = parseArgs()

    #seed random number generator
    #random.seed(params.seed)
    if not params.seed:
        params.seed=datetime.now().timestamp()
    random.seed(params.seed)

    #########################################################
    # Step 1: Parsing shapefile
    #########################################################

    # Initialize the ResistanceNetwork object
    network = ResistanceNetwork(
        network=params.network,
        shapefile=params.shapefile,
        coords=params.coords,
        variables=params.variables,
        agg_opts = params.agg_opts,
        pop_agg = params.pop_agg,
        inmat = params.inmat,
        reachid_col=params.reachid_col,
        length_col=params.length_col,
        out = params.out, 
        verbose=True
    )


    # #########################################################
    # # Step 2: Initialise ModelRunner
    # #########################################################

    runner = ModelRunner(
        resistance_network=network,
        seed = params.seed,
        verbose=True
    )

    # #########################################################
    # # Step 3: Run GA optimisation
    # #########################################################

    runner.run_ga(
        maxgens = params.maxGens, 
        fitmetric = params.fitmetric,
        threads = params.GA_procs,
        burnin = params.burnin,
        deltaB = params.deltaB, 
        deltaB_perc = params.deltaB_perc,
        indpb = params.indpb, 
        mutpb = params.mutpb, 
        cxpb = params.cxpb,
        nFail = params.nfail,
        popsize = params.popsize, 
        maxpopsize = params.maxpopsize, 
        posWeight = params.posWeight, 
        fixWeight = params.fixWeight, 
        fixShape = params.fixShape, 
        allShapes = params.allShapes,
        min_weight = params.min_weight, 
        max_shape = params.max_shape, 
        max_hof_size = params.max_hof_size,
        awsum = params.awsum, 
        only_keep = params.only_keep,
        out = params.out, 
        plot=True, 
        verbose=True
    )

#Call main function
if __name__ == '__main__':
    main()
