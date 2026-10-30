# ResistNet

`ResistNet` is a Python program for optimising environmental resistance models in dendritic spatial networks, such as riverscapes. It includes utilities for optimising models from an arbitrary number of environmental covariates using Bayesian optimisation (Tree-structured Parzen Estimators; TPE), simulating datasets (e.g., for power analyses), and model-averaging results within runs or across replicates.

## Table of Contents:
1. [Installation](#install)
    1. [Installing from GitHub](#install_dev)
2. [Usage](#usage)
    1. [Command-line Interface](#usage_cli)
    2. [Input Files](#inputs)
    3. [Outputs](#outputs)
    4. [Optimisation Parameters](#opt_params)
    5. [The Standard Model](#standard)
    6. [Directional Models (SAMC)](#samc)
    7. [Interactive Report](#html_report)
    8. [Comparing the Two Models](#compare)
3. [Tutorials](#tutorials)
    1. [Example Dataset](#example)
    2. [Utility Scripts](#scripts)
    3. [Integrating into workflows](#workflows)
4. [Contributing](#contrib)
    1. [Reporting Issues](#report)
    2. [Making Changes](#changes)
    3. [Unit Tests and CI](#testci)

## Installation <a name="install"></a>

The simplest way to install the most recent release is using conda/mamba.

With a conda installation, first create a new environment, and activate it:
```
conda create -n resistnet python=3.10
conda activate resistnet
```

Then, install resistnet (note I use mamba here, but conda will work as well):
```
mamba install -c conda-forge -c bioconda -c ecoevoinfo resistnet
```

### Installation from GitHub

To instead install the development version on GitHub, you first need to clone the repository:

```
git clone https://github.com/tkchafin/resistnet.git
cd resistnet
```

Then, you can use the `environment.yml` file to install the dependencies via conda or mamba, and then activate the environment:

```
mamba env create -f environment.yml
mamba activate resistnet
```

Then you can install `resistnet` using pip:
```
pip install -e .
```

## Usage <a name="usage"></a>

### Command-line Interface <a name="usage_cli"></a>

The core `ResistNet` functionality is accessed via the `runResistnet.py` script. This can be found in the `scripts/` directory (if cloning the repository from GitHub), or wherever your binaries are installed (if installed via conda/ pip).

All executables bundled in `ResistNet` have a help menu which can be displayed by calling them with `-h`. To view the help menu for `runResistnet.py`, simply enter:

```
runResistnet.py -h
```

This will show the currently available options:

```
usage: runResistnet.py [-h] [-g INMAT] [-s SHAPEFILE] [-c COORDS] [-n NETWORK]
                       [--seed SEED] [--reachid_col REACHID_COL]
                       [--length_col LENGTH_COL] [-t PROCS] [-o OUT]
                       [-O {SHP,GPKG,GDB}] [--gridSearch]
                       [--transFile TRANSFILE] [-v VARS] [-V VARFILE]
                       [-F NFAIL] [-i MAX_ITER] [-r REPS] [-P PWEIGHT]
                       [-S NSTART] [-C NCAND] [-G GAMMA]
                       [-f {aic,loglik,r2m,delta}] [--fixWeight] [--fixShape]
                       [--max_shape MAX_SHAPE] [--min_weight MIN_WEIGHT]
                       [--model {resistance,samc,both}] [--samc]
                       [--sizes SIZES] [--ne NE] [--infer_origin INFER_ORIGIN]
                       [--max_resist MAX_RESIST] [--fix_resist FIX_RESIST]
                       [--fixAsym] [--max_flow MAX_FLOW] [--fixFlow]
                       [--samc_absorption {none,sites}]
                       [--samc_pairs {max,mean}] [--cell_km CELL_KM]
                       [-a AWSUM] [--use_full] [--report_all]

Optimise resistance networks using Tree-structured Parzen Estimators (TPE)

options:
  -h, --help            show this help message and exit
  -g, --inmat INMAT     Input distance matrix (usually genetic distances)
                        (default: None)
  -s, --shapefile SHAPEFILE
                        Path to shapefile, geodatabase, or GPKG file (default:
                        None)
  -c, --coords COORDS   Input tsv containing sample coordinates (default:
                        None)
  -n, --network NETWORK
                        Pre-formatted pickled network file (default: None)
  --seed SEED           Random number seed (default: None)
  --reachid_col REACHID_COL
                        Reach ID in geodatabase (default: HYRIV_ID)
  --length_col LENGTH_COL
                        Length attribute in geodatabase (default: LENGTH_KM)
  -t, --procs PROCS     Number of parallel processors (default: 1)
  -o, --out OUT         Output directory and prefix: -o DIR/NAME writes each
                        model's outputs to DIR/NAME/{model}/NAME.* (and, with
                        --model both, the comparison to DIR/NAME/NAME.*)
                        (default: output)
  -O, --output_driver {SHP,GPKG,GDB}
                        Output driver for annotated geodataframe (default:
                        GPKG)
  --gridSearch          Run grid-search to find fixed transformations
                        (default: False)
  --transFile TRANSFILE
                        Read fixed transformations from file (default: None)
  -v, --vars VARS       Comma-separated list of variables to use (default:
                        None)
  -V, --varfile VARFILE
                        Optional file with variables provided as: var \t
                        <Optional aggregator function> (default: None)
  -F, --nfail NFAIL     Number of iterations without improvement to stop
                        (counted after the --nstart random start-up
                        iterations) (default: 200)
  -i, --max_iter MAX_ITER
                        Maximum number of iterations (one model per chain
                        each) (default: 500)
  -r, --reps REPS       Number of independent TPE chains (default: 10)
  -P, --pweight PWEIGHT
                        TPE prior weight (default: 0.7)
  -S, --nstart NSTART   Initial random evaluations per chain (default: 20)
  -C, --ncand NCAND     TPE expected improvement candidate points (default:
                        48)
  -G, --gamma GAMMA     TPE gamma (0-1; fraction of trials treated as "good")
                        (default: 0.15)
  -f, --fitmetric {aic,loglik,r2m,delta}
                        Fitness metric used to evaluate models (delta = change
                        in AIC relative to the null model) (default: aic)
  --fixWeight           Fix all variable weights at 1.0 (default: False)
  --fixShape            Turn off variable transformations (default: False)
  --max_shape MAX_SHAPE
                        Maximum shape parameter (shapes from 1 are searched
                        log-uniformly) (default: 15)
  --min_weight MIN_WEIGHT
                        Minimum weight parameter (default: 0.0)
  --model {resistance,samc,both}
                        Model to fit: resistance (sums of reach resistances
                        along the paths between sites; default), samc (spatial
                        absorbing Markov chain: conditional first passage
                        times among sites, allowing directional effects), or
                        both (fitted separately, then ranked together)
                        (default: None)
  --samc                Same as --model samc (default: False)
  --sizes SIZES         SAMC: tsv with columns sample, size (effective
                        population sizes; absorption at each site is 1/(2 *
                        size)) (default: None)
  --ne NE               SAMC: effective population size of every site, if no
                        --sizes (default: 100)
  --infer_origin INFER_ORIGIN
                        SAMC: attribute giving the next downstream reach ID
                        (default: NEXT_DOWN)
  --max_resist MAX_RESIST
                        Upper limit for the maximum resistance per km
                        (relative to a minimum of 1), which is optimised (as
                        ResistanceGA's max.cont) (default: 1000.0)
  --fix_resist FIX_RESIST
                        Fix the maximum resistance instead of optimising it
                        (default: None)
  --fixAsym             SAMC: do not consider directional effects of variables
                        (always the case on the lattice, where only the flow
                        term is directional; kept for compatibility) (default:
                        False)
  --max_flow MAX_FLOW   SAMC: upper limit for the rate of the flow (upstream
                        vs downstream) term, per km (default: 0.1)
  --fixFlow             SAMC: leave out the flow (upstream vs downstream) term
                        (default: False)
  --samc_absorption {none,sites}
                        SAMC: no absorption (passage times are hitting times),
                        or absorption 1/(2 Ne) at each site per visit
                        (default: none)
  --samc_pairs {max,mean}
                        SAMC: how the passage times of the two directions
                        between two sites are combined (default: max)
  --cell_km CELL_KM     SAMC: distance between neighbouring states of the
                        lattice the walker moves on (km; samc's cell size)
                        (default: 1.0)
  -a, --awsum AWSUM     Cumulative Akaike weight threshold (default: 1.0)
  --use_full            Use full rather than natural model averaging (default:
                        False)
  --report_all          Plot outputs for all retained models (default: False)
```

The use of these options are described below.

### Input Files <a name="inputs"></a>

#### Coordinates 

Coordinates for your populations or individuals should be provided as a tab-delimited table, using the `-c/ --coords` argument. An example can be found in `src/resistnet/data/nmt.pointCoords.txt`:
```
sample  lat     long
CCFE    38.37291666666618       -96.49375000000065
CCOT    38.38541666666615       -96.54791666666732
CESG    38.38124999999949       -96.164583333334
NAME    38.51041666666617       -96.30625000000067
NDUN    38.54374999999947       -96.34375000000068
...
...
```

The columns required are "sample", "lat", and "long", any other present columns will be ignored. Note that samples (or populations) defined here will require labelling which is consistent in the genetic distance matrix (described below). 

#### Genetic distance matrix 

The input genetic distance matrix, supplied via the `-g/ --inmat` argument, should be a tab-delimited file with both column and row names, matching the samples provided in the coordinates file. 

The example dataset comes with an input genetic distance matrix (pairwise linearised F<sub>ST</sub>, F<sub>ST</sub>/(1 - F<sub>ST</sub>)) for reference, at `src/resistnet/data/nmt.popGenDistMat.txt`:

```
        CCFE    CCOT    CESG    NAME    ...
CCFE    0.0     0.0059491931752027214   0.0038011061193934034   0.08831964633972295     ...
CCOT    0.0059491931752027214   0.0     0.021867693862724382    0.1290622519846013      ...
...
...
```

There are a number of packages available to generate this pairwise distance matrix; the example's was produced with [`autoStreamTree`](https://github.com/tkchafin/autostreamtree) (`-d LINFST -C loc`, Weir & Cockerham's F<sub>ST</sub> with the SNPs of each locus joined).

#### Geodatabase

The input stream network can be provided as a shapefile, geodatabase, or GPKG file, all passed using the `-s/--shapefile` option. Alternatively, a network previously built by `ResistNet` or `autoStreamTree` (e.g., `$out.subgraph.net`) can be passed with `-n/--network`, which avoids re-building the network from the shapefile. There are a number of requirements for this file in order for the result to create a valid network. I highly recommend using the existing global stream datasets provided by the [HydroLab group](https://wp.geog.mcgill.ca/hydrolab/) at McGill University, specifically the [HydroAtlas](https://www.hydrosheds.org/page/hydroatlas) or [free-flowing rivers dataset](https://wp.geog.mcgill.ca/hydrolab/free-flowing-rivers/) as these are already properly formatted for use, and contain a large number of covariates as reach-level annotations which can be analysed by `ResistNet`. You will also need to provide the reach identifier `--reachid_col` (default is "HYRIV_ID") and reach length (`--length_col`) (default "LENGTH_KM"). 

If for some reason you cannot use the HydroRIVERS dataset, you will need to do some things first before loading your shapefile into autoStreamTree. First, you will need to include two variables in the attribute table of your shapefile: 1) a reach ID must provide a unique identifier to each stream reach; and 2) a reach length variable should give the length of each segment. Next, because sometime large stream layers will have small gaps in between streams, you will need to span any small gaps between streams which should be contiguous, and also dissolve any lines that overlap with one another so that any given section of river is represented by a single line. There are some scripts in a complementary package that can help with these steps using the ArcPy API: <https://github.com/stevemussmann/StreamTree_arcpy>.

Note that a valid path is required between all sites. Thus, if you are analyzing data from multiple drainages which only share an oceanic connection (or a connection which is otherwise absent from your dataset), you will need to augment the shapefile. For example this could be accomplished by adding a vector representing the coastline to create an artificial connection among drainages.


### Outputs <a name="outputs"></a>

After parsing all of the inputs, `ResistNet` runs `-r,--reps` independent TPE chains. In each iteration, every chain proposes one model parameterization (which variables are included, and their weights, transformations, and shapes), and these are evaluated in parallel by `-t,--procs` worker processes. Each model is fitted to the genetic distances with a maximum-likelihood population effects (MLPE) mixed model (Clarke et al. 2002, as implemented in [ResistanceGA](https://github.com/wpeterman/ResistanceGA)), and scored using the fitness metric chosen with `-f,--fitmetric` (by default AIC). The first `-S,--nstart` models in each chain are random; after that, TPE uses the previous evaluations to propose new models. After each iteration, the best result so far across all chains is printed to the terminal (stdout), as a loss to be minimised (for the default `-f aic` this is the AIC itself; for the other metrics it is their negative):
```
Standard model: resistance per km, from 1 to a maximum optimised up to 1000
Building network from shapefile: ./src/resistnet/data/nmt.gpkg
WARNING: This can take a while with very large files!
Read 14 points.


Extracting full subgraph...

Merging redundant paths...


Starting worker 0 with seed 1234
Starting worker 1 with seed 1235
Starting worker 2 with seed 1236
Starting worker 3 with seed 1237

Iteration 1, Current Best: -373.28673696921066
Iteration 2, Current Best: -378.79724054203325
Iteration 3, Current Best: -378.79724054203325
Iteration 4, Current Best: -378.79724054203325
...
...
```
The optimisation stops after `-i,--max_iter` iterations, or earlier if the best result has not improved for `-F,--nfail` consecutive iterations (not counting the `-S,--nstart` random start-up iterations, which explore the space before the search proper). The best model found by each chain is then reported:
```
...
...
Iteration 50, Current Best: -451.4939698069141
Chain 1: Best = -451.43103963132376
Chain 2: Best = -411.3533917455675
Chain 3: Best = -451.4939698069141
Chain 4: Best = -427.39680589186753
...
...
```

The best model from each chain is retained in the Hall of Fame, which is used for multi-model inference (Akaike weights, relative variable importance, and model-averaged weights and resistances). A series of plots and tables will then be produced in a directory for the model: with `-o DIR/NAME`, in `DIR/NAME/{model}/`, named `NAME.*`, where `{model}` is `resistance` or `samc` (see `--model`). For example, `-o results/nmt` writes `results/nmt/resistance/nmt.HallOfFame.tsv` and so on (below, `output.*`):

| File Name                                    | Description                                                                                                                     |
|----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| output.minimalSubgraph.net                   | Compressed representation of the minimized subgraph                                                                             |
| output.HallOfFame.tsv                        | Full text specification of the retained models for the run                                                                      |
| output.bestModel.tsv                         | Weights, transformations, and shapes of the best model                                                                          |
| output.modelSpec.json                        | Each Hall of Fame model with the constants needed to apply it to another network (see `predictResistnet.py`)                    |
| output.ICprofile.pdf                         | AIC of the retained models                                                                                                      |
| output.minimalSubgraph.pdf                   | Graphical representation of the minimized subgraph                                                                              |
| output.modavgWeights.tsv                     | Table giving the model-averaged weights for each covariate                                                                      |
| output.Model-Average.PairwiseRegplot.pdf     | Regression of pairwise resistance X genetic distances                                                                           |
| output.pairPlot.pdf                          | Plots of pairwise correlation between fitness metrics across sampled models                                                     |
| output.Model-Average.ResistanceEdges.tsv     | Resistance of each reach by reach ID: per km, not multiplied by the reach's length                                                |
| output.pointsTable.txt                       | Table of points                                                                                                                 |
| output.Model-Average.ResistanceMatrix.tsv    | Pairwise effective resistance distances                                                                                         |
| output.snapDistances.txt                     | Table giving the distance (in km) between input sampling points and nearest graph nodes (useful to diagnose mis-snapped points) |
| output.Model-Average.streamsByResistance.pdf | Map of the model-averaged resistance per km of each reach in the network                                                        |
| output.Model-Best.*                          | As for the output.Model-Average.* files, for the best model                                                                     |
| output.subgraph.net                          | Compressed representation of the subgraph                                                                                       |
| output.Null-Model.tsv                        | Table giving fitness metrics for the distance-only and null (population effect only) models                                     |
| output.varImportance.pdf                     | Bar plot of relative variable importances                                                                                       |
| output.incidenceMatrix.txt                   | Incidence matrix for paths between each pair of samples (only used internally)                                                  |
| output.pairwise.tsv                          | Every pair of sites: genetic distance, river distance (km), and the best and model-averaged models' distances (for SAMC, the mean of both directions) |
| output.varImportance.tsv                     | Table giving the model-averaged relative variable importances for each covariate                                                |
| output.gpkg                                  | The input geodatabase, subset to the reaches in the subgraph (only with `-s`; format set with `-O,--output_driver`)             |
| output.univariateFitness.tsv/.pdf            | Fitness of each transformation and shape for each variable (only with `--gridSearch`)                                           |
| output.fixedParams.tsv                       | Fixed transformations and shapes used during optimisation (only with `--gridSearch` or `--transFile`)                           |

A full set of example outputs, from the commands in the [Example Dataset](#example) section, can be seen in `src/resistnet/data/example_output` (`resistance/` for the standard model, `samc/` for the SAMC model). If using the `--report_all` option, these plots will also be produced for every one of the "kept" models from the Hall of Fame, with naming as `$out.Model-#.PairwiseRegplot.pdf` (etc), where "#" is the row number from the `HallofFame.tsv` table (with 0-based indexing; i.e. 0="best" model; 1=second-best, and so on).


#### $out.HallOfFame.tsv 

The output table `$out.HallOfFame.tsv` contains textual representations of the retained models for a run, with a format like:

```
fitness run_mm_cyr      run_mm_cyr_weight       ... dor_pc_pva      dor_pc_pva_weight       dor_pc_pva_trans        dor_pc_pva_shape        max_resistance ... akaike_weight   acc_akaike_weight       keep
-451.4939698069141      0       -       ... 1       0.7540635366427507      7.0     1.0     144.60309855972295      ... 0.5078641129098956      0.5078641129098956      True
```

The table contains the best model from each of the `-r,--reps` TPE chains. Some of these columns are self-explanatory (i.e., "fitness" is the model fitness, which with the default `-f aic` is the AIC). Other summary metrics included are the Akaike weights, change in AIC relative to the null and "best" models, as well as the other fitness metrics (marginal R^2, log-likelihoods). As in ResistanceGA, the AIC counts both the parameters of the MLPE model (intercept, slope, and two variance components) and the resistance surface parameters: a transformation and shape for each included variable, their weights, and the maximum resistance per km (`max_resistance`; see [The standard model](#standard)). Only relative weights matter (the composite is a weighted mean), so a model with m included variables has m - 1 weight parameters, and 3m surface parameters in all, as ResistanceGA's equation, shape and maximum for each surface. Parameters fixed with `--fixWeight`, `--fixShape` or `--fix_resist` are not counted. The marginal R^2 is computed as in `MuMIn::r.squaredGLMM` (as used by ResistanceGA). Other values in the table include whether or not a variable is included in a model ($var, where $var is the variable name), specified as either "1" (=included) or "0" (=excluded), and the transformation type (=$var_trans), transformation shape parameter (=$var_shape), and weight of the parameter when calculating the composite resistance edges (=$var_weight). For the transformation column, the types of transformations are as follows: 0=Not transformed; 1=Ricker; 2=Reverse Ricker; 3=Inverse Ricker; 4=Reverse-Inverse Ricker (the inverse Ricker of the reversed values, as in ResistanceGA; previously this duplicated the Ricker); 5=Monomolecular; 6=Reverse Monomolecular; 7=Inverse Monomolecular; 8=Reverse-Inverse Monomolecular. In all cases, the larger the shape value, the closer each transformation gets to being linear (=essentially no transformation).

#### Covariate relative importance 

Other outputs that will typically be of interest are the `modAvgWeights.tsv` and `.varImportance.tsv` files, which give the model-averaged weights and importances for all included covariates:

```
# RVI
variable        RVI
dor_pc_pva      0.9999999990242653
tmp_dc_cyr      2.9734278463245265e-06
run_mm_cyr      0.0
sgr_dk_rav      0.0
ele_mt_cav      0.0

# MAW
variable        MAW
dor_pc_pva      0.8740561031325752
tmp_dc_cyr      0.05346755565785958
run_mm_cyr      -
sgr_dk_rav      -
ele_mt_cav      -
```

Note that relative importances only compare the retained resistance models with each other, so some variable will always appear important. To judge whether the landscape variables explain genetic distances at all, compare the AIC of the best model with that of the distance-only model in `$out.Null-Model.tsv`; the Hall of Fame gives this difference for each model as `delta_aic_distance` (positive = the model fits better than distance alone). In this example the best model (degree of regulation, `dor_pc_pva`, i.e. the effect of dams on flow, with a resistance per km from 1 to 145; AIC -451.49) fits much better than distance alone (AIC -386.51), a difference of 65.0:

```
                loglik                  aic                     r2m                     delta_aic_null
distance_only   197.25314955094646      -386.5062991018929      0.21330431591693536     28.126047275811516
null            182.1901259130407       -358.3802518260814
```

Distance is part of every model (see [The standard model](#standard)); in `$out.bestModel.tsv`, the `Contribution` column gives how much of the best model's resistance per km is due to distance (here 0.43) and to each variable (0.57 for `dor_pc_pva`). This short example run is not exhaustive, though here it finds the best single-variable model of longer runs: in three runs of 600 iterations x 8 chains, the best model was degree of regulation (inverse monomolecular, shape 1; resistance 1-145 per km; AIC -451.5) in two, and degree of regulation with runoff (`run_mm_cyr`; resistance 1-584 per km; AIC -458.7) in the third.

#### Effective resistance values 

The model-averaged resistances are computed for each reach in the subgraph (per km: from 1, distance alone, up to the models' maximum resistance), and also output for each individual model if using `--report_all`. A reach's contribution to the path sums is its resistance per km times its length (`--length_col`). 

If you would like to import these into other GIS plotting tools, you can simply perform a left join of the `.ResistanceEdges.tsv` table, using the reach IDs from `--reachid_col` (HYRIV_ID in this example; the column is named EDGE_ID):

```
EDGE_ID Resistance
70605797        88.5442340913543
70606280        90.01565746601429
70606496        91.51181008601205
70606497        91.51181008601205
...
...
```

These are also plotted for simple visualisation on the input spatial network, in the `.streamsByResistance.pdf` files. Pairwise effective values (between points) are also output as tab-delimited text (`.ResistanceMatrix.tsv`).

### Optimisation parameters <a name="opt_params"></a>

`ResistNet` searches the space of models using the TPE algorithm implemented in [`hyperopt`](https://github.com/hyperopt/hyperopt). The relevant options are:

```
  -F, --nfail NFAIL     Number of iterations without improvement to stop
                        (counted after the --nstart random start-up
                        iterations) (default: 200)
  -i, --max_iter MAX_ITER
                        Maximum number of iterations (one model per chain
                        each) (default: 500)
  -r, --reps REPS       Number of independent TPE chains (default: 10)
  -P, --pweight PWEIGHT
                        TPE prior weight (default: 0.7)
  -S, --nstart NSTART   Initial random evaluations per chain (default: 20)
  -C, --ncand NCAND     TPE expected improvement candidate points (default:
                        48)
  -G, --gamma GAMMA     TPE gamma (0-1; fraction of trials treated as "good")
                        (default: 0.15)
  --fixWeight           Fix all variable weights at 1.0 (default: False)
  --fixShape            Turn off variable transformations (default: False)
  --max_shape MAX_SHAPE
                        Maximum shape parameter (shapes from 1 are searched
                        log-uniformly) (default: 15)
  --min_weight MIN_WEIGHT
                        Minimum weight parameter (default: 0.0)
  --gridSearch          Run grid-search to find fixed transformations
                        (default: False)
  --transFile TRANSFILE
                        Read fixed transformations from file (default: None)
```

For each variable, the search space includes whether it is included in the model, its weight (between `--min_weight` and 1.0), one of the nine transformations listed above, and the transformation shape (a whole number between 1 and `--max_shape`, default 15). Shapes are drawn log-uniformly: on the 0-10 scale, small shapes give the most different curves, while large shapes all give nearly linear ones, so small shapes are tried more often. In simulations, this found a narrow optimum at a small shape (a Ricker with shape 3 on a skewed covariate) in 3 of 4 runs, against 0 of 4 when shapes were drawn uniformly from 1 to 100, and cost nothing when the truth was nearly linear (shape 50). As in ResistanceGA (shapes 0.5-14.5), the default maximum is about 15. `--fixWeight` fixes all weights at 1.0, and `--fixShape` turns off the transformations.

Each of the `-r,--reps` chains is independent, and up to `-i,--max_iter` x `-r,--reps` models will be evaluated in total. Because the best model from each chain is retained in the Hall of Fame, `-r,--reps` also sets the number of models available for model-averaging. `-t,--procs` only sets how many models are evaluated in parallel.

Optionally, transformations can be chosen for each variable separately before the TPE search. With `--gridSearch`, every combination of transformation and shape (up to `--max_shape`) is evaluated for each variable on its own (results are written to `$out.univariateFitness.tsv`, with heatmaps in `$out.univariateFitness.pdf`). The best transformation and shape for each variable are then fixed during the TPE search, which still optimises which variables are included and their weights (and the maximum resistance per km; the grid search uses 100, or `--fix_resist`). Fixed transformations can also be read from a file with `--transFile`: either a previous `$out.univariateFitness.tsv` (the best values are selected), or a tab-delimited table with columns `variable`, `transform`, and `shape` (e.g., a previous `$out.fixedParams.tsv`). Columns `sel` and `weight` can also be included, to fix whether variables are included and their weights.

Note that the "right" parameters will vary by context, and you should not necessarily just run using the defaults!


### The standard model <a name="standard"></a>

The standard model (`--model resistance`, the default) predicts genetic distances from the resistance summed along the river between sites (the effective resistance, which on a river network is the sum of the reaches' resistances along the path). Each reach's resistance is built as follows:

* Each selected variable (aggregated over the reach as given in the variables file) is rescaled to 0-10, transformed (see [$out.HallOfFame.tsv](#outputs) for the transformations), and rescaled to 0-1.
* The variables are combined as their weighted mean, the composite (0-1).
* The resistance per km is 1 + (K - 1) x composite: 1 where the composite is 0, up to the maximum resistance per km K, which is optimised (log-uniformly from 1 to `--max_resist`, default 1000; or fixed with `--fix_resist`).
* A reach's resistance is its resistance per km times its length (`--length_col`).

This is the formulation of [ResistanceGA](https://github.com/wpeterman/ResistanceGA) ([Peterman 2018](https://doi.org/10.1111/2041-210X.12984)), for reaches of different lengths. In ResistanceGA, each continuous surface is 1 + max x f(x), the surfaces are summed (and divided by their minimum), and the conductance between neighbouring cells is divided by the distance between them (`geoCorrection` in gdistance), so each step costs resistance x distance: resistances are per unit of distance. `ResistNet`'s resistance per km equals ResistanceGA's combined surface (up to a constant factor, which doesn't change the fit), with K - 1 and the weights in place of each surface's maximum (the tests check this against ResistanceGA's own functions, to 1e-12), and the AIC counts the same number of parameters (3 per variable). Distance is therefore always part of the model: K = 1 is the `distance_only` reference model (isolation by distance), and a large K means the variables dominate. How much of the resistance per km is due to distance and to each variable, averaged over the network's km (as ResistanceGA's `p.contribution`), is given in the `Contribution` column of `$out.bestModel.tsv` (with a `distance` row), and each model's AIC difference from distance alone in the Hall of Fame (`delta_aic_distance`). Because distance is already included, `LENGTH_KM` should not be used as a covariate.

In earlier versions of `ResistNet`, each reach cost its composite value (rescaled to 0-1) whatever its length, so distance could only enter as a covariate (e.g. `LENGTH_KM`). Reaches can differ in length by more than an order of magnitude (0.4-13.6 km in the example data), so that formulation counts a short reach the same as a long one. It is no longer available: it fitted worse, and it is not ResistanceGA's formulation. On the example data, full searches (600 iterations x 8 chains, three seeds) found slightly better models per km (best AIC -458.7, median -451.5) than per reach (best -455.1, median -451.2; with `LENGTH_KM` among the candidates, best -452.4, median -450.0; it was never selected, per km or per reach). In simulations with the full search, AIC preferred per km in 9 of 12 datasets simulated per km (by a median of 4.5-10 AIC), and the two tied on data simulated per reach with `LENGTH_KM` as a covariate (within 3.7 AIC); only data simulated per reach with no distance at all favoured per reach (4 of 4, by 1.2-7.9).


### Directional models (SAMC) <a name="samc"></a>

Effective resistance (circuit theory) is symmetric: moving from A to B costs the same as moving from B to A. With `--model samc` (or `--samc`), `ResistNet` instead models movement as a spatial absorbing Markov chain (SAMC; [Fletcher et al. 2019](https://doi.org/10.1111/ele.13333)), as in the [samc](https://cran.r-project.org/package=samc) R package ([Marx et al. 2020](https://doi.org/10.1111/ecog.04891)), translated from raster cells to the river network:

* samc's states are raster cells of equal size. Here they are the points of a lattice of the river: the nodes of the network (where reaches meet, and the sites, which are snapped to nodes) and points splitting each reach into pieces of about `--cell_km` (default 1 km), so that every step covers about the same distance. The network includes the reach downstream of the outlet (found using the `--infer_origin` attribute, `NEXT_DOWN` by default) and one reach beyond each terminal sample site.
* Each reach's resistance per km is the standard model's (ResistanceGA's surface: 1 + (K - 1) x the weighted mean of the transformed variables, K optimised from 1 to `--max_resist`, default 1000 as for the standard model, or fixed with `--fix_resist`; see [The standard model](#standard)). Points inside a reach take its resistance; a node takes the length-weighted mean over the half-pieces of river around it, as a raster cell over a confluence would.
* As in samc, a walker moves between neighbouring points with conductance 1 / (mean resistance of the two points x the distance between them), normalised over its neighbours and scaled by 1 - absorption. Reaches meeting at a node are all neighbours through it.
* By default nothing absorbs the walker (`--samc_absorption none`), so the passage times are hitting times. With `--samc_absorption sites`, a walker at a site is absorbed with probability 1/(2 Ne) per step, the probability that two lineages coalesce; effective population sizes can be provided per sample with `--sizes` (a tsv with columns `sample` and `size`), otherwise all sites have Ne = `--ne` (default 100).
* A *flow* term can make moving upstream harder or easier: each move's weight is multiplied by exp(-rate x km moved upstream). A positive rate makes the walker's upstream moves harder, a negative rate easier. samc's transition functions can depend on the order of the two cells, so this is a samc model with position along the river as the variable. The rate is per km and is not rescaled to the network, so it means the same in any watershed. Whether the flow term is included, its direction and its rate (from 0.001 to `--max_flow` per km, default 0.1) are optimised, unless `--fixFlow` is used (no flow term). A model can consist of the flow term alone. Variables cannot be directional (on the lattice, a reach's points share its values; `--fixAsym` is kept for compatibility).
* The distance from site i to site j is the conditional first passage time (CFPT): the expected number of steps to reach j, given that it is reached before absorption elsewhere (as `samc::cond_passage` in R). The two directions between a pair of sites are combined (`--samc_pairs`: their maximum, the default, or their mean) and fitted to the genetic distances with MLPE.

This translation is checked against samc itself: for a small river drawn on a raster, samc's transition matrix and passage times, built by samc from the raster, equal the lattice's (`tests/resistnet/test_lattice.py`).

**Reading the flow rate.** Genetic distances are symmetric, so a direction can only show through how distances vary along the river. We tested this with simulations of known gene flow (msprime stepping-stone models on the example network and the SMB Black River network, 400 SNPs, migration biased downstream or upstream by up to 0.05 per km, with and without a barrier). With the defaults (no absorption, maximum of the two directions), the fitted rate had the sign *opposite* to the bias of gene flow in 29 of 32 biased datasets, never the same sign, and found a direction in 5 of 24 datasets without bias: **the walker moves like a lineage traced back in time, so a negative rate means gene flow biased downstream.** The rate is attenuated (about a fifth to a tenth of the bias) and, with a barrier that the model does not represent, the direction was missed in half of the biased datasets. With absorption at the sites or the mean of the two directions (the formulation of earlier versions), the fitted direction did not follow the truth (directions were found in most datasets without bias). In simulations with less information (200 SNPs) the direction was mostly missed but never reversed. An AIC gain above 2 found a direction in about a fifth of the datasets without bias, where a likelihood-ratio test at that level would expect one in twenty, so pairwise MLPE overstates the evidence: a gain of 6 or more removed the false positives of the noisier simulations without losing a true one. Gradients along the river (e.g. temperature, which rises downstream) and changes in population density can mimic a direction: in the stepping-stone simulations, whose deme sizes did not follow the length of river, false directions reached gains of 25. Treat the flow term as suggestive, report it alone and with the chosen covariates, and represent known barriers in the model.

The AIC counts the MLPE parameters, the transformation and shape of each selected variable, their relative weights, the maximum resistance per km, and the flow rate (if the model has a flow term). In addition to the outputs described above, SAMC runs produce:

| File Name                                    | Description                                                                                         |
|----------------------------------------------|-----------------------------------------------------------------------------------------------------|
| output.Model-Best.CFPTMatrix.tsv             | Conditional first passage times among sites (row = origin, column = destination); not symmetric     |
| output.Model-Best.TransitionEdges.tsv        | Each reach's resistance per km, and the probabilities of moving downstream and upstream inside it (set by the flow term) |
| output.Model-Best.streamsByTransition.pdf    | Maps of the resistance per km and of the probability of moving downstream                           |
| output.adjacencyMatrix.tsv                   | Reach adjacency (1 = column reach is downstream of row reach; -1 = upstream)                        |
| output.graphOrigin.pdf, graphDirectionality.pdf | The inferred outlet, and the flow direction of each reach                                        |
| output.upstreamChanges.tsv, .pdf             | Share of upstream steps on which each variable rises, stays the same, or falls                      |
| output.flowProfile.tsv                       | AIC of flow alone at each rate tried when fitting it, in both directions (rate 0 = distance alone)  |

(and the same for `Model-Average`). The Hall of Fame includes `$var_asym` (always 0), `max_resistance` and `flow` (the rate; 0 = no flow term) columns. `Null-Model.tsv` includes two SAMC reference models: `samc_distance`, distance alone (resistance 1 per km everywhere, i.e. no landscape or flow effect), and `samc_flow`, flow alone (the same plus the flow term, with its rate fitted and given in the `flow` column). The Hall of Fame reports each model's AIC difference from them, `delta_aic_distance` and `delta_aic_flow` (positive = the model fits better, as for `delta_aic_null`), to judge whether the variables explain genetic distances beyond distance and flow alone.

If the downstream reach of the outlet is not in the input network (e.g. it was clipped), the reaches meeting at the outlet are linked to each other directly, with a warning; including it is preferable. For example, using the bundled dataset (see [Example Dataset](#example)):

```
runResistnet.py \
    -s ./src/resistnet/data/nmt.gpkg \
    -c ./src/resistnet/data/nmt.pointCoords.txt \
    -g ./src/resistnet/data/nmt.popGenDistMat.txt \
    -V ./src/resistnet/data/selected_vars.txt \
    --model samc --sizes ./src/resistnet/data/nmt.popSizes.txt \
    -i 50 -r 4 -t 4 --seed 1234 -o nmt_samc
```

In this short run, the best SAMC model is stream gradient (`sgr_dk_rav`; AIC -412.1), better than the reference model of flow alone (-0.100 per km, the limit, i.e. gene flow biased downstream; -402.8) and SAMC distance alone (-387.9). A longer run (600 iterations x 8 chains) finds degree of regulation (`dor_pc_pva`, inverse monomolecular, resistance 1-14.6 per km) without a flow term (AIC -424.0), less good than the standard model's best (-452.3 in the same search; see below). The two models fit the same genetic distances with the same MLPE, so their AICs can be compared directly.

### Interactive report <a name="html_report"></a>

Every run also writes an interactive HTML report, `DIR/NAME/NAME.report.html` (for `-o DIR/NAME`), with Plotly figures of each model's results: the best model and reference models, the Hall of Fame, the variables' importance (and, for the standard model, distance's share of the best model's resistance per km), the best model's response to each variable (as ResistanceGA's `Plot.trans`, with the network's reaches marked), the fit to the genetic distances (hover for each pair of sites), and a map of the network coloured by the model-averaged or best model's resistances per km (on a log scale; or, for SAMC, the probabilities of moving downstream and upstream), with the sampling sites; and, with `--model both`, the comparison of the two models. The report works offline (one file, about 5 MB), and is built only from the run's output files, so it can be rebuilt (e.g. after editing them) with:

```
reportResistnet.py -o DIR/NAME
```

Each figure's data is in one of the run's tables (e.g. `*.pairwise.tsv`, `*.HallOfFame.tsv`, `*.TransitionEdges.tsv`).

### Comparing the two models <a name="compare"></a>

With `--model both`, `ResistNet` fits the standard model and then the SAMC model, each with its own search, Hall of Fame, model averaging and outputs (in `DIR/NAME/resistance/` and `DIR/NAME/samc/`), and then ranks the models of both, with the reference models, by AIC. Both are fitted with the same MLPE to the same genetic distances, so their AICs are comparable; this is checked (the null and distance-only models of both must be identical). These outputs are written to `DIR/NAME/`:

| File Name                    | Description                                                                                                                                                                                              |
|------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| NAME.modelComparison.tsv     | Every Hall of Fame model of both (each once) and the reference models, ranked by AIC: family, model, number of parameters (k), log-likelihood, AIC, difference from the best, and Akaike weight over all of them |
| NAME.familySummary.tsv       | For each family (resistance, samc, reference): its best model, that model's AIC and difference from the best, and the family's summed Akaike weight                                                     |
| NAME.jointImportance.tsv     | Relative importance of each variable over both families (the summed Akaike weight of the models that include it), and the part of it from each family's models                                         |
| NAME.modelComparison.pdf     | The ranking as AIC differences from the best model, coloured by family, with each family's summed Akaike weight                                                                                        |
| NAME.jointImportance.pdf     | The joint importance of each variable, split by family                                                                                                                                                  |
| NAME.bestFits.pdf            | The genetic distances against river distance (isolation by distance) and against each family's best model, side by side                                                                               |

For example, adding `--model both` to the command in the [Example Dataset](#example) section gives `nmt/resistance/`, `nmt/samc/` and:

```
family      best_model                                                           best_aic   delta_aic  akaike_weight
resistance  dor_pc_pva (inverse monomolecular, shape 1); resistance 1-145 per km     -451.49    0.00       1.00
samc        sgr_dk_rav (reverse-inverse monomolecular, shape 1)                  -412.06    39.44      1.6e-09
reference   SAMC flow alone (-0.100/km)                                          -402.76    48.73      1.3e-11
```

This short run ranks the standard model first. With 600 iterations x 8 chains (seed 1234), the standard model's best is elevation with degree of regulation (resistance 1-995 per km; AIC -452.3; degree of regulation alone is 0.8 behind), and SAMC's best is degree of regulation without flow (-424.0), 28 AIC worse. Compare the two models only after runs long enough for both searches to settle. AIC ranks fits, not mechanisms (path resistance vs a walker's passage time), and near-ties of a few AIC units should not be over-read, since the parameter counts for the maximum resistance and flow are conventions.


## Tutorials <a name="tutorials"></a>

### Example Dataset <a name="example"></a>

The bundled example dataset (`src/resistnet/data/nmt.*`) is from the Neosho madtom (*Noturus placidus*), a small catfish endemic to the Neosho River Basin (Kansas, Oklahoma and Missouri, USA), from [Zbinden et al. (2025)](https://doi.org/10.1016/j.gecco.2025.e03783):

* `nmt.gpkg`: 176 [HydroRIVERS](https://www.hydrosheds.org/products/hydrorivers) reaches with all [HydroATLAS](https://www.hydrosheds.org/hydroatlas) covariates (Lehner & Grill 2013; Linke et al. 2019). HydroATLAS is licensed under CC-BY 4.0, and HydroRIVERS under the [HydroSHEDS license agreement](https://www.hydrosheds.org/products/hydrorivers); both require attribution. `nmt.network` is the same network, pre-built for `-n`.
* `nmt.pointCoords.txt`: the 14 sampling sites.
* `nmt.popGenDistMat.txt`: pairwise linearised F<sub>ST</sub> among the sites (from the study's SNPs, with `autoStreamTree`; see [Genetic distance matrix](#gen)).
* `nmt.popSizes.txt`: a placeholder table of effective population sizes for `--sizes` (100 for every site, as the default `--ne`).
* `selected_vars.txt`: five covariates: runoff (`run_mm_cyr`), air temperature (`tmp_dc_cyr`), stream gradient (`sgr_dk_rav`), elevation (`ele_mt_cav`) and degree of regulation (`dor_pc_pva`, the effect of dams on flow).

Zbinden et al. (2025) found that dams constrain dispersal among Neosho madtom populations. In this example, the best model is degree of regulation (`dor_pc_pva`), which fits far better than distance alone. Note that HydroATLAS's degree of regulation reflects the storage of large reservoirs: here it mainly marks the upper Neosho and the Spring River, whose sites are cut off from the others by lowhead dams near Emporia and by Grand Lake, neither of which is among the bundled covariates. You can very quickly test out if your installation is running, and explore what the outputs look like or profile runtimes, with a short run:

```
runResistnet.py \
    -s ./src/resistnet/data/nmt.gpkg \
    -c ./src/resistnet/data/nmt.pointCoords.txt \
    -g ./src/resistnet/data/nmt.popGenDistMat.txt \
    -V ./src/resistnet/data/selected_vars.txt \
    -i 50 -r 4 -t 4 --seed 1234 -o nmt
```

To fit the SAMC model as well and compare the two, add `--model both` (see [Comparing the two models](#compare)). The outputs of the command with `--model both` are in `src/resistnet/data/example_output` (`resistance/`, `samc/`, and the comparison tables). Note that this example is not necessarily exhaustive (real runs will likely need many more iterations and chains), so you can try playing around with the running parameters, as well as the covariates included (specified in the `selected_vars.txt` file; the network contains all of the HydroATLAS covariates).

> Zbinden, Z. D., Chafin, T. K., Tiemann, J. S., Edds, D. R., Martin, B. T., Hofmeier, J., Douglas, M. E., & Douglas, M. R. (2025). Historic and contemporary selection define conservation units for a short-range endemic within an anthropogenically-altered riverscape. *Global Ecology and Conservation*, 62, e03783. https://doi.org/10.1016/j.gecco.2025.e03783


### Utility Scripts <a name="scripts"></a>

The primary functions of `ResistNet` are accessed via the command-line interface `runResistnet.py`. However, there are also four other utilities that are installed alongside `runResistnet.py` (whether installing via conda or pip, or found directly in the `scripts/` directory of the GitHub repository): `simResistnet.py`, `ensembleResistnet.py`, `predictResistnet.py` and `reportResistnet.py` (see [Interactive report](#html_report)). 

#### Simulating datasets with simResistnet.py

As with `runResistnet.py`, `simResistnet.py` has a command-line interface and help menu which can be accessed by calling the script with `-h`:

```
simResistnet.py

Author: Tyler Chafin
Description: Simulate data on a given network

Arguments:
    -n, --network: Input network (pickle'd networkx output)
    -i, --in: Table giving variables to use to generate resistnet input
    -r, --reps: Number of replicates
    -s, --samples: Number of random nodes to sample
    -l, --len_col: Edge length attribute (def=LENGTH_KM)
    -c, --id_col: Reach ID attribute (def=EDGE_ID)
    -o, --out: Output file name (default=sim)
    -z, --seed: Random seed (optional)
    --max_resistance: Maximum resistance: per km (standard model; each
        reach's resistance is 1 + (K - 1) x composite per km, times its
        length) or of the composite (SAMC) (default=100)

SAMC (directional) simulations:
    --samc: Simulate conditional first passage times (mean of the two
        directions) under a SAMC model, instead of effective resistances;
        the input table may have an ASYM column (1 = directional)
    --ne: Effective population size at every site (default=100)
    --flow: Flow rate per km (positive = moving upstream harder;
        default=0, no flow term)
    --infer_origin: Next downstream reach attribute (def=NEXT_DOWN)
```

Simulations require a compressed input network (i.e., as output by `runResistnet.py` or `autoStreamTree`), and a specifications file which should be tab-delimited and provide some covariates you would like to simulate effective resistances with, along with their weights, transformations, and shapes (using the same codes as in the Hall of Fame). An optional `AGG` column gives the method used to aggregate each variable across reaches sharing an ID (default ARITH, as in `runResistnet.py`). Resistances are computed in exactly the same way as when fitting models with `runResistnet.py`, so simulated models can be recovered using the same settings: per km, with the maximum resistance per km set by `--max_resistance` (default 100). An example can be found in `src/resistnet/data/simparams.tsv`:
```
VAR     WEIGHT  TRANSFORM       SHAPE
dor_pc_pva      1.0     7       1
```

In this example, the resistance of each reach comes from its degree of regulation (`dor_pc_pva`) with an inverse monomolecular transformation (shape 1), the best model for the example dataset. The pairwise resistance matrix is then normalised between 0 and 1.

You can specify the number of samples with `-s/--samples`, which determines the number of "populations" samples (as random nodes in the graph), as well as the number of replicated datasets to produce, using `-r/--reps`. For example, to simulate 10 datasets, each with 50 random samples and the above specifications:

```
simResistnet.py -n src/resistnet/data/nmt.network -c HYRIV_ID -s 50 -r 10 -o sim -i src/resistnet/data/simparams.tsv
```

This will produce a number of outputs named `*.ResistanceMatrix.tsv` and `*.coords`, which can be supplied directly as inputs to `runResistnet.py` via `-g` and `-c` (with the same network and `--reachid_col`), for example to test if `ResistNet` model optimisation can recover the specified model in your network, given different sample sizes.

With `--samc`, the distances are instead the conditional first passage times of a SAMC model (see [Directional models](#samc)), the two directions combined as when fitting, and rescaled to 0-1; the passage times themselves (row = origin, column = destination) are written to `*.CFPTMatrix.tsv`. `--max_resistance`, `--flow`, `--ne` and `--cell_km` set the maximum resistance per km, the flow rate (per km), the effective population size at every site and the lattice's distance between states. With no variables (a header line only), the resistance is 1 per km, so e.g. `--flow 0.05` simulates flow alone. For example:

```
simResistnet.py -n src/resistnet/data/nmt.network -c HYRIV_ID -s 20 -r 10 \
    -o samc_sim -i spec.tsv --samc --max_resistance 50 --flow -0.05
```

The movement weights and passage times are computed exactly as when fitting SAMC models, so simulated models can be recovered with `runResistnet.py --samc` (with the same `--ne` and `--cell_km`). These simulations are SAMC's own model, so they test the search, not whether SAMC describes gene flow (see [Directional models](#samc)).

Note that because samples are random nodes, the network used for simulations should have a unique ID for every reach (such as `HYRIV_ID` in HydroRIVERS). If IDs are shared by several reaches (e.g., the `EDGE_ID` segments in `autoStreamTree` outputs, which are defined between the original sampling sites), new samples can fall within these segments, which cannot be represented. If the network has several unconnected parts (e.g. other basins in the reach layer, as in `nmt.network`), samples are drawn from the largest. 

#### Model-averaging across runs with ensembleResistnet.py

`runResistnet.py` by default performs model-averaging within runs, but if you want to consider outputs across replicate runs (or for example, the top or 'best' model from a series of independent runs), you can use the utility `ensembleResistnet.py`, which has the following options:

```
ensembleResistnet.py

Utility for model-averaging across ResistNet outputs

Arguments:
-s, --shp: Path to shapefile
-n, --network: Input network
-i, --in: Directory containing resistnet outputs (searched with its
          subdirectories)
-L, --list: Optional comma-separated prefixes (no spaces)
-C, --coords: Optional static coordinates file
--model: Only use runs of this model (resistance or samc), i.e. outputs
         in {out}/{model}/ directories

Optional Arguments:
-t, --procs: Number of parallel processors
--seed: RNG seed
-a, --awsum: Cumulative Akaike weight threshold [default=0.95]
--only_keep: Only retain models where column 'keep'=True
--only_best: Only retain best model from each input
--split_samples: Treat all samples as unique
-X, --noPlot: Turn off plotting
-m, --max_keep: Maximum models to keep (default = all models)
-l, --len_col: Edge length attribute (def=LENGTH_KM)
-c, --id_col: Reach ID attribute (def=HYRIV_ID)
-o, --out: Output file prefix (default=ensemble)
--report_all: Plot full outputs for all retained models
-V, --varfile: Optional file with variables provided like so:
          var1   <Optional aggregator function>
          var2   <Optional aggregator function>
          ...

SAMC runs (detected from their Hall of Fame columns):
--sizes: Effective population sizes (tsv: sample, size)
--ne: Effective population size for all sites without
      --sizes [default=100]
```

Most arguments are used in the same way as their equivalents in `runResistnet.py`, but one difference is in how a group of replicate run outputs are passed. Runs are found in the directory given and its subdirectories (e.g. the `DIR/NAME/{model}/` directories of `runResistnet.py`), and `--model resistance` or `--model samc` uses only the runs of that model. Runs of the SAMC model are recognised from their Hall of Fame columns and averaged as a SAMC model (with `--sizes` or `--ne` as for `runResistnet.py`), giving model-averaged passage times and transition probabilities; standard and SAMC runs can't be mixed in one ensemble. Standard runs of earlier versions (resistances per reach, with no `max_resistance` column) are no longer supported. The simplest way is to provide the path to an output directory, where outputs are named like `{prefix}_{rep}.*`. Only the `$prefix.HallOfFame.tsv` files are read. Example outputs, from 10 short runs of the example dataset with different seeds, can be found in `src/resistnet/data/nmt_ensemble/replicates`:
```
nmt_1.HallOfFame.tsv
nmt_10.HallOfFame.tsv
nmt_2.HallOfFame.tsv
...
...
```

Creating a model-average across the 10 replicate analyses in this directory can be done like so:

```
ensembleResistnet.py \
    -n src/resistnet/data/nmt.network \
    -i src/resistnet/data/nmt_ensemble/replicates \
    -c HYRIV_ID -C src/resistnet/data/nmt.pointCoords.txt \
    -o ensemble -t 4
```

To instead create a model-average of only the best models sampled in each run, you can add the `--only_best` argument:

```
ensembleResistnet.py \
    -n src/resistnet/data/nmt.network \
    -i src/resistnet/data/nmt_ensemble/replicates \
    -c HYRIV_ID -C src/resistnet/data/nmt.pointCoords.txt \
    -o ensemble -t 4 --only_best
```

#### Applying a model to another network with predictResistnet.py

Each run writes `$out.modelSpec.json`, giving every Hall of Fame model (best first) together with the constants needed to compute it again: for each variable, its range in the network the model was fitted to (used for the 0-10 rescaling), its transformation and shape, its range after transformation, and its weight; and the maximum resistance per km (standard models; models of earlier versions, per reach, are no longer supported; for SAMC models, also the flow rate and the lattice's `cell_km`; SAMC models of earlier versions, which moved between reaches, are no longer supported). `predictResistnet.py` applies one of these models to the reaches of another network, e.g. another watershed:

```
predictResistnet.py -s other_watershed.shp -m out/resistance/out.modelSpec.json -o other
```

Because every rescaling uses the fitted network's ranges, a reach with the same values as a reach of the fitted network gets the same result. Values outside those ranges are clamped to them (with a warning giving the share of reaches affected for each variable), so predictions far outside the fitted network's conditions are flat rather than extrapolated. The variables' combined value is a weighted mean of their clamped values, so it is always within range. If most reaches are affected, the predicted values barely differ from each other. Variables whose ranges are tied to position in the watershed (e.g. `DIST_DN_KM`, the distance to the sea) rarely overlap between watersheds, so models using them transfer poorly; the SAMC flow rate is per km and transfers as it is. Models are applied reach by reach: if the fitted network merged reaches into longer segments (e.g. the `EDGE_ID` segments of `autoStreamTree` outputs), their values were aggregated over those segments.

For a standard model, the output (`other.predicted.tsv`, and `other.predicted.gpkg` with the reaches) gives each reach's resistance per km (`resistance`, from 1 to the model's maximum resistance per km) and, if the reaches have `--length_col`, that times the reach's length (`reach_resistance`, which is what sums along paths). A SAMC model gives the same (its resistance per km is the standard model's) and the probability of moving downstream inside each reach (`p_downstream`: from a state whose two neighbours are on the reach, so it depends only on the flow term and the reach's link length, `--length_col` / its number of `cell_km` pieces; NA for reaches of one piece). `-n` selects a model by its rank (default: the best).

### Integration Into Workflows <a name="workflows"></a>

All utilities from `ResistNet` are intended to be easy to run as a sequential pipeline. An example of building these into a [`SnakeMake`](https://snakemake.github.io) workflow can be found in the [Open Science Framework](https://osf.io/4uxfj/) repository (doi: 10.17605/OSF.IO/4UXFJ).

For example, in the `run_validation.smk` file, you can find everything needed to perform simulations of varying sampling sizes across a number of specified models, across 3 example networks. 

For further examples on incorporating `ResistNet` with feature engineering (via robust PCA) and feature selection (via random forests), see the [Open Science Framework repository](https://osf.io/9zsf7/) for an upcoming paper (to be posted soon!)
## Contributing <a name="contrib"></a>

We welcome contributions in all forms, be that reporting bugs, requesting/ discussing new features, or pull requests! Please see below for some general guidelines, and please help to maintain a helpful and inclusive dialogue by following the [Contributer Covenant Code of Conduct](https://www.contributor-covenant.org/version/1/4/code-of-conduct/).

### Reporting Issues <a name="report"></a>

If you encounter a bug or any other issue, please use the [GitHub Issues](https://github.com/tkchafin/resistnet/issues) page, after first checking that someone else hasn't already reporting the same problem!

When describing your issue, please be detailed, and include any relevant system details, as well as the full command-line prompt you used, as well as any necessary data files to replicate the problem. 

You can also use this page to post any feature requests you might have. 

### Making changes <a name="changes"></a>

If you would like to make any changes or add features to `ResistNet` (always welcome!), the procedure is simple:

1. Fork the repository
2. Make changes 
3. Submit a pull request 

Note that we use `flake8` to enforce readable styling, and also have a set of integrated unit tests which will be run on submitting the pull request (see below section on continuous integration).

### Unit Tests and CI <a name="testci"></a>

We use [`flake8`](https://flake8.pycqa.org/en/latest/) to ensure that all code contributions follow a readable style, and a set of unit tests in [`pytest`](https://docs.pytest.org/en/7.4.x/) to ensure that the basic functionality of `ResistNet` is maintained before merging any changes.

Both of these are run automatically via `GitHub Actions` when you submit a pull request to the repository. If you would like to first run these locally, you can do so on the command-line, after installing both tools with conda/ mamba:

```
mamba install flake8 pytest
```

To run `flake8` linting:

```
cd resistnet
flake8 src/*.py
```

And to run the unit tests:

```
cd resistnet
pytest tests/
```

Note that these tests are NOT comprehensive! We currently have a test coverage of ~85%, and welcome any contributions to make the testing framework more robust! 
