#!/usr/bin/env python3
import os
import random
from datetime import datetime

from resistnet.params import parseArgs
from resistnet.resistance_network import ResistanceNetwork
from resistnet.samc_network import ResistanceNetworkSAMC
from resistnet.model_optimisation import ModelRunnerTPE
from resistnet.compare import write_comparison


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

    # Each model's outputs go to {out}/{model}/{name}.*, where name is the
    # last part of --out
    name = os.path.basename(os.path.normpath(params.out))
    runs = {}
    for model in params.models:
        prefix = os.path.join(params.out, model, name)
        os.makedirs(os.path.dirname(prefix), exist_ok=True)
        if len(params.models) > 1:
            print(f"\n=== Fitting the {model} model ===\n")

        # Step 1: Initialise network
        network = build_network(params, model, prefix)

        # write new output geodatabase
        if params.shapefile:
            network.write_geodataframe(prefix, params.output_driver)

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
                out=prefix,
                plot=True,
                verbose=True
            )

        # Step 2: Run optimisation
        runner.run_tpe(
            max_evals=params.max_iter,
            fitmetric=params.fitmetric,
            threads=params.procs,
            nFail=params.nfail,
            fixWeight=params.fixWeight,
            fixShape=params.fixShape,
            fixAsym=params.fixAsym,
            fixFlow=params.fixFlow,
            min_weight=params.min_weight,
            max_shape=params.max_shape,
            awsum=params.awsum,
            only_keep=True,
            use_full=params.use_full,
            report_all=params.report_all,
            fixed_params=params.transFile,
            out=prefix,
            plot=True,
            verbose=True,
            pweight=params.pweight,
            reps=params.reps,
            n_startup=params.nstart,
            n_candidates=params.ncand,
            gamma=params.gamma
        )
        runs[model] = (network, runner, prefix)

    # Step 3: With both models, rank their models together
    if len(runs) > 1 and params.inmat:
        write_comparison(runs, os.path.join(params.out, name))

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


def build_network(params, model, prefix):
    """
    Builds the network for a model ("resistance" or "samc").
    """
    if model == "samc":
        return ResistanceNetworkSAMC(
            network=params.network,
            shapefile=params.shapefile,
            sizes=params.sizes,
            coords=params.coords,
            variables=params.variables,
            agg_opts=params.agg_opts,
            inmat=params.inmat,
            reachid_col=params.reachid_col,
            length_col=params.length_col,
            infer_origin=params.infer_origin,
            ne=params.ne,
            max_resist=params.max_resist,
            fix_resist=params.fix_resist,
            max_flow=params.max_flow,
            out=prefix,
            verbose=True
        )
    return ResistanceNetwork(
        network=params.network,
        shapefile=params.shapefile,
        coords=params.coords,
        variables=params.variables,
        agg_opts=params.agg_opts,
        inmat=params.inmat,
        reachid_col=params.reachid_col,
        length_col=params.length_col,
        out=prefix,
        verbose=True
    )


# Call main function
if __name__ == '__main__':
    os.environ['USE_PYGEOS'] = '0'
    main()
