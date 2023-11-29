# ResistNet

`ResistNet` is a Python program for optimising environmental resistance models in dendritic spatial networks, such as riverscapes. It includes utilities for optimising models from an arbitrary number of environmental covariates using a genetic algorithm, simulating datasets (e.g., for power analyses), and model-averaging results within runs or across replicates.

## Table of Contents:
1. [Installation](#install)
    1. [Installing from GitHub](#install_dev)
2. [Usage](#usage)
    1. [Command-line Interface](#usage_cli)
    2. [Input Files](#inputs)
    3. [Outputs](#outputs)
    4. [GA Parameters](#ga_params)
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
resistnet.py

Author: Tyler K Chafin
Contact: tyler.chafin@bioss.ac.uk
Description: Genetic algorithm to optimise resistance networks

Input options:
    -g, --genmat: Genetic distance matrix
    -s, --shp: Path to shapefile, geodatabase, or GPKG file
    -c, --coords: Input tsv containing sample coordinates

General options:
    --seed: Random number seed (default=taken from clock time)
    --reachid_col: Reach ID [default="REACH_ID"]
    --length_col: Length [default="LENGTH_KM"]
    -t, --procs: Number of parallel processors
    -X, --noPlot: Turn off plotting
    -o, --out: Output file prefix
    -h, --help: Displays help menu

Aggregation options:
    --edge_agg: Method for combining variables across segments
    --pop_agg: Method to combine population genetic distances
      Options: ARITH (-metic mean), MEDIAN, HARM (-monic mean),
      ADJHARM (adjusted HARM, see docs), GEOM (geometric mean),
      MIN, MAX, FIRST, SD (standard deviation), VAR (variance),
      SUM (simple sum), CV (coefficient of variation = SD/MEAN)

Genetic Algorithm Options:
    -P, --maxPop: Maximum population size [default = 100]
    -G, --maxGen: Maximum number of generations [default = 500]
    -s, --size: Manually set population size to <-p int>,
        NOTE: By default, #params * 15
    -m, --mutpb: Mutation probability per trait [default=0.5]
    --indpb: Mutation probability per individual [default=0.5]
    --cxpb: Cross-over probability [default=0.5]
    -T, --tSize: Tournament size [default=10]
    --posWeight: Constrain parameter weights to between 0.0-1.0
    --minWeight: Sets minimum allowable weight (w/--posWeight)
    --fixWeight: Constrain parameter weights to 1.0
    --fixShape: Turn off feature transformation
    --allShapes: Allow inverse and reverse transformations
    --maxShape: Maximum shape value [default=100]

Model optimization/selection options:
    -v, --vars: Comma-separated list of variables to use
    -V, --varfile: Optional file with variables provided as:
        var1 \t <Optional aggregator function>
        var2 \t <Optional aggregator function>
        ...
    -F, --nfail: Number of failed gens to stop optimization
    -d, --delt: Threshold absolute change in fitness [def.=0.0]
    -D, --deltP: Threshold as decimal percentage [def.=0.001]
    -f, --fit: Fitness metric used to evaluate models
        Options: aic (default), loglik (log-likelihood),
        r2m (marginal R^2), delta (Change in AIC vs null model)
        NOTE: Case-insensitive
    -b, --burn: Number of generations for pre-burnin [def.=0]
    --max_hof_size: Maximum models retained [default=100]

Multi-model inference options:
    -a, --awsum: Cumulative Akaike weight threshold [def.=0.95]
    --report_all: Plot outputs for all retained models
```

The use of these options are described below.

### Input Files <a name="inputs"></a>

#### Coordinates 

Coordinates for your populations or individuals should be provided as a tab-delimited table, using the `-c/ --coords` argument. An example can be found in `src/resistnet/data/test.pointCoords.txt`:
```
sample  lat     long
burk    26.927083333332924      90.39374999999947
cdikc   27.264583333332887      90.03958333333279
dakp    27.14791666666624       90.68958333333276
digl    26.893749999999564      91.75208333333276
dikc    26.881249999999547      90.26874999999947
...
...
```

The columns required are "sample", "lat", and "long", any other present columns will be ignored. Note that samples (or populations) defined here will require labelling which is consistent in the genetic distance matrix (described below). 

#### Genetic distance matrix 

The input genetic distance matrix, supplied via the `-g/ --genmat` argument, should be a tab-delimited file with both column and row names, matching the samples provided in the coordinates file. 

The example dataset comes with an input genetic distance matrix for reference, at `src/resistnet/data/test.popGenDistMat.txt`:

```
        burk    cdikc   dakp    digl    ...
burk    0.0     0.6913154781112187      0.6943568032698719      0.6975921469002094      ... 
cdikc   0.6913154781112187      0.0     0.6004024031774708      0.6024853508009713      ... 
... 
...
```

There are a number of packages available to generate this pairwise distance matrix. For convenience, all files necessary to analyse the example dataset in [`autoStreamTree`](https://github.com/tkchafin/autostreamtree) are also included in the `data/` directory. 

#### Geodatabase

The input stream network can be provided as a shapefile, geodatabase, or GPKG file, all passed uing the `-s/--shp` option. There are a number of requirements for this file in order for the result to create a valid network. I highly recommend using the existing global stream datasets provided by the [HydroLab group](https://wp.geog.mcgill.ca/hydrolab/) at McGill University, specifically the [HydroAtlas](https://www.hydrosheds.org/page/hydroatlas) or [free-flowing rivers dataset](https://wp.geog.mcgill.ca/hydrolab/free-flowing-rivers/) as these are already properly formatted for use, and contain a large number of covariates as reach-level annotations which can be analysed by `ResistNet`. You will also need to provide the reach identifier `--reachid_col` (default is "HYRIV_ID") and reach length (`--length_col`) (default "LENGTH_KM"). 

If for some reason you cannot use the HydroRIVERS dataset, you will need to do some things first before loading your shapefile into autoStreamTree. First, you will need to include two variables in the attribute table of your shapefile: 1) a reach ID must provide a unique identifier to each stream reach; and 2) a reach length variable should give the length of each segment. Next, because sometime large stream layers will have small gaps in between streams, you will need to span any small gaps between streams which should be contiguous, and also dissolve any lines that overlap with one another so that any given section of river is represented by a single line. There are some scripts in a complementary package that can help with these steps using the ArcPy API: <https://github.com/stevemussmann/StreamTree_arcpy>.

Note that a valid path is required between all sites. Thus, if you are analyzing data from multiple drainages which only share an oceanic connection (or a connection which is otherwise absent from your dataset), you will need to augment the shapefile. For example this could be accomplished by adding a vector representing the coastline to create an artificial connection among drainages.


### Outputs <a name="outputs"></a>

After parsing all of the inputs, `ResistNet` will randomly generate a population of 'individuals' (=model parameterizations), which is by default 15X the number of parameter, up to a maximum size specified by `-P,--maxPop`. Each 'generation', individuals will be selected using the desired fitness function (specified with `-f,--fit`; e.g. `-f AIC` to use AIC), and the maximum, minimum, mean, and standard deviation of population fitness values will be output to the terminal (stdout):
```
Building network from shapefile: /Users/tyler/resistnet/src/resistnet/data/test.shp
WARNING: This can take a while with very large files!
Read 37 points.


Extracting full subgraph...

Merging redundant paths...

Initializing genetic algorithm parameters...

Establishing a population of size: 10

Starting worker 0 with seed 1234
Starting worker 1 with seed 1235
Starting worker 2 with seed 1236
Starting worker 3 with seed 1237

-- Generation 1 --
  Worst -1200.4747471654575
  Best -1163.9990149524226
  Avg -1178.6391094740154
  Std 15.774111494879024
  nFails 0
...
...

```
Also included is the number of consecutive generations that the genetic algorithm has failed to find a 'better' model (defined using thresholds set with `-d,--delt` or `-D,--deltP`). After either a specified number of generations (`-G,--maxGen`) have passed, or nFail exceeds the user-specified number (`-F,--nFail`), the algorithm will stop, and report to the screen which stopping criteria was met:
```
...
...
-- Generation 5 --
  Worst -1217.011714840949
  Best -1134.9919179075587
  Avg -1162.9318785505388
  Std 32.02090507102848
  nFails 2
Stopping optimization after 5 generations.
Reason: Exceeded maxGens
...
...
```

At this time a series of plots and tables will be produced using the output path/ prefix defined with `-o,--out`):

| File Name                                    | Description                                                                                                                     |
|----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| output.FitnessLog.tsv                        | Log of fitness min, mean, max, and spread across generations                                                                    |
| output.minimalSubgraph.net                   | Compressed representation of the minimized subgraph                                                                             |
| output.HallOfFame.tsv                        | Full text specification of the retained models for the run                                                                      |
| output.minimalSubgraph.pdf                   | Graphical representation of the minimized subgraph                                                                              |
| output.modavgWeights.tsv                     | Table giving the model-averaged weights for each covariate                                                                      |
| output.Model-Average.PairwiseRegplot.pdf     | Regression of pairwise resistance X genetic distances                                                                           |
| output.pairPlot.pdf                          | Plots of pairwise correlation between fitness metrics across sampled models                                                     |
| output.Model-Average.ResistanceEdges.tsv     | Effective resistance values for each segment by reach ID                                                                        |
| output.pointsTable.txt                       | Table of points                                                                                                                 |
| output.Model-Average.ResistanceMatrix.tsv    | Pairwise effective resistance distances                                                                                         |
| output.snapDistances.txt                     | Table giving the distance (in km) between input sampling points and nearest graph nodes (useful to diagnose mis-snapped points) |
| output.Model-Average.streamsByResistance.pdf | Plot of model-averaged effective resistance values for each segment in the network                                              |
| output.subgraph.net                          | Compressed representation of the subgraph                                                                                       |
| output.Null-Model.tsv                        | Table giving fitness metrics for the distance-only and null (population effect only) models                                     |
| output.varImportance.pdf                     | Bar plot of relative variable importances                                                                                       |
| output.incidenceMatrix.txt                   | Incidence matrix for paths between each pair of samples (only used internally)                                                  |
| output.varImportance.tsv                     | Table giving the model-averaged relative variable importances for each covariate                                                |

A full set of example outputs can be seen in `src/resistnet/data/example_output`. If using the `--report_all` option, these plots will also be produced for every one of the "kept" models from the Hall of Fame, with naming as `$out.Model-#.Pairwise.pdf` (etc), where "#" is the row number from the `HallofFame.tsv` table (with 0-based indexing; i.e. 0="best" model; 1=second-best, and so on).


#### $out.HallOfFame.tsv 

The output table `$out.HallOfFame.tsv` contains textual representations of the retained models for a run, with a format like:

```
fitness run_mm_cyr      run_mm_cyr_weight       run_mm_cyr_trans        run_mm_cyr_shape ... akaike_weight   acc_akaike_weight       keep 
1165.2161145651664      1       0.5381247419039517      1   ... 0.9999999697774224      0.9999999697774224      True
```

Some of these are self-explanatory (i.e., "fitness" is the model fitness, which by default is AIC). Other summary metrics included are the Akaike weights, change in AIC relative to the null and "best" models, as well as the other fitness metrics (marginal R^2, log-likelihoods). Other values in the table include whether or not a variable is included in a model ($var, where $var is the variable name), specified as either "1" (=included) or "0" (=excluded), and the transformation type (=$var_trans), transformation shape parameter (=$var_shape), and weight of the parameter when calculating the composite resistance edges (=$var_weight). For the transformation column, the types of transformations are as follows: 0=Not transformed; 1=Ricker; 2=Reverse Ricker; 3=Inverse Ricker; 4=Reverse-Inverse Ricker; 5=Monomolecular; 6=Reverse Monomolecular; 7=Inverse Monomolecular; 8=Reverse-Inverse Monomolecular. In all cases, the larger the shape value, the closer each transformation gets to being linear (=essentially no transformation).

#### Covariate relative importance 

Other outputs that will typically be of interest are the `modAvgWeights.tsv` and `.varImportance.tsv` files, which give the model-averaged weights and importances for all included covariates:

```
# RVI
variable        RVI
run_mm_cyr      0.9999999992728014
sgr_dk_rav      0.9999999777043529
tmp_dc_cyr      0.999999974582615
dor_pc_pva      5.403280690464821e-10

# MAW
variable        MAW
sgr_dk_rav      0.9438627249077166
run_mm_cyr      0.5381247256404349
dor_pc_pva      -0.0
tmp_dc_cyr      -0.24636058036382147
```

#### Effective resistance values 

The model-averaged effective resistance values are computed for each segment in the subgraph, and also output for each individual model if using `--report-all`. 

If you would like to import these into other GIS plotting tools, you can simply perform a left join of the `.ResistanceEdges.tsv` table, using the specified `--reachid_col` (EDGE_ID in this example):

```
EDGE_ID Resistance
40830357        8.842640424633727
40832271        8.22649846960403
40833054        7.926278241834325
40833364        4.568026911233788
...
...
```

These are also plotted for simple visualisation on the input spatial network, in the `.streamsByResistance.pdf` files. Pairwise effective values (between points) are also output as tab-delimited text (`.ResistanceMatrix.tsv`).
### Genetic Algorithm parameters <a name="ga_params"></a>

ResistNet provides options for manipulating the relevant parameters of the genetic algorithm:

```
Genetic Algorithm Options:
    -P, --maxPop: Maximum population size [default = 100]
    -G, --maxGen: Maximum number of generations [default = 500]
    -s, --size: Manually set population size to <-p int>,
        NOTE: By default, #params * 15
    -m, --mutpb: Mutation probability per trait [default=0.5]
    --indpb: Mutation probability per individual [default=0.5]
    --cxpb: Cross-over probability [default=0.5]
    -T, --tSize: Tournament size [default=10]
    --posWeight: Constrain parameter weights to between 0.0-1.0
    --minWeight: Sets minimum allowable weight (w/--posWeight)
    --fixWeight: Constrain parameter weights to 1.0
    --fixShape: Turn off feature transformation
    --allShapes: Allow inverse and reverse transformations
    --maxShape: Maximum shape value [default=100]
```

The `--posWeight` and `--fixWeight` options are used to either constrain parameter weights to 0.0 - 1.0, or to fix all weights at 1.0. By default, weights can vary from -1.0 to 1.0, and all 'negatively weighted' transformations are not available (inverse and reverse). Only ricker, monomolecular, and inverse-reverse versions of both are available unless `--allShapes` is used. I don't recommend using both negative weights and `--allShapes` together, because this creates a situation where there are multiple ways to get the same possible result (i.e., ricker with -1.0 weighting has the same impact on composite surface as reverse ricker with 1.0 weighting). 

Note that the "right" parameters will vary by context, and you should not necessarily just run using the defaults!


## Tutorials <a name="tutorials"></a>

### Example Dataset <a name="example"></a>

A full example dataset is included in `src/resistnet/data`. You can very quickly test out if your installation is running, and explore what the outputs look like or profile runtimes by performing a short run using the bundled dataset:

```
runResistnet.py \
    -s ./src/resistnet/data/test.shp \
    -c ./src/resistnet/data/test.pointCoords.txt \
    -g ./src/resistnet/data/test.popGenDistMat.txt \
    -V ./src/resistnet/data/selected_vars.txt \
    -P 10 -G 5 -t 4 --seed 1234
```

Note that this example is not necessarily useful (i.e., real runs will likely need to be run for much longer, and with a larger population size), so you can try playing around with the running parameters, as well as the covariates included (specified in the `selected_vars.txt` file). This dataset contains all of the [HydroATLAS](https://www.hydrosheds.org/hydroatlas) covariates.


### Utility Scripts <a name="scripts"></a>

The primary functions of `ResistNet` are accessed via the command-line interface `runResistnet.py`. However, to there are also currently two other utilities that are installed alongside `runResistnet.py` (whether installing via conda or pip, or found directly in the `scripts/` directory of the GitHub repository): `simResistnet.py` and `ensembleResistnet.py`. 

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
```

Simulations require a compressed input network (i.e., as output by `runResistnet.py` or `autoStreamTree`), and a specifications file which should be tab-delimited and provide some covariates you would like to simulate effective resistances with. An example can be found in `src/resistnet/data/simparams.tsv`:
```
VAR     WEIGHT  TRANSFORM       SHAPE
LENGTH_KM       1.0     0       0
```

In this example, a pairwise resistance matrix will simply be the sum of pairwise segment lengths (i.e., isolation-by-distance), which will then be normalised between 0 and 1. 

You can specify the number of samples with `-s/--samples`, which determines the number of "populations" samples (as random nodes in the graph), as well as the number of replicated datasets to produce, using `-r/--reps`. For example, to simulate 10 datasets, each with 50 random samples and the above specifications:

```
simResistnet.py -n src/resistnet/data/test.network -s 50 -r 10 -o sim -i src/resistnet/data/simparams.tsv
```

This will produce a number of outputs named `*.ResistanceMatrix.tsv` and `*.coords`, which can be supplied directly as inputs to `runResistnet.py` via `-g` and `-c`, for example to test if `ResistNet` model optimisation can recover the specified model in your network, given different sample sizes. 

#### Model-averaging across runs with ensembleResistnet.py

`runResistnet.py` by default performs model-averaging within runs, but if you want to consider outputs across replicate runs (or for example, the top or 'best' model from a series of independent runs), you can use the utility `ensembleResistnet.py`, which has the following options:

```
ensembleResistnet.py

Utility for model-averaging across ResistNet outputs

Arguments:
-s, --shp: Path to shapefile
-n, --network: Input network
-i, --in: Directory containing resistnet outputs
-L, --list: Optional comma-separated prefixes (no spaces)
-C, --coords: Optional static coordinates file

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
-c, --id_col: Reach ID attribute (def=EDGE_ID)
-o, --out: Output file prefix (default=ensemble)
--report_all: Plot full outputs for all retained models
--allShapes: Allow inverse and reverse transformations
-V, --varfile: Optional file with variables provided like so:
          var1   <Optional aggregator function>
          var2   <Optional aggregator function>
          ...
```

Most arguments are used in the same way as their equivalents in `runResistnet.py`, but one difference is in how a group of replicate run outputs are passed. The simplest way is to provide the path to an output directory, where outputs are named like `{prefix}_{rep}.*`. Example outputs can be found in `src/resistnet/data/test_ensemble/replicates`:
```
test_1.out.FitnessLog.tsv
test_1.out.HallOfFame.tsv
test_1.out.ICprofile.pdf
test_1.out.Model-Average.Mantel.tsv
test_1.out.Model-Average.PairwiseRegplot.pdf
test_1.out.Model-Average.ResistanceEdges.tsv
test_1.out.Model-Average.ResistanceMatrix.tsv
test_1.out.Model-Average.streamsByResistance.pdf
test_1.out.genDistMat.tsv
test_1.out.genDistMat.txt
test_1.out.incidenceMatrix.txt
test_1.out.minimalSubgraph.net
test_1.out.minimalSubgraph.pdf
test_1.out.modavgWeights.tsv
test_1.out.pairPlot.pdf
test_1.out.pointsTable.txt
test_1.out.snapDistances.txt
test_1.out.subgraph.net
test_1.out.varImportance.pdf
test_1.out.varImportance.tsv
test_10.out.FitnessLog.tsv
test_10.out.HallOfFame.tsv
...
...
...
```

Creating a model-average across the 20 replicate analyses in this directory can be done like so:

```
ensembleResistnet.py \
    -n data/test_ensemble/test.full.network \
    -i data/test_ensemble/replicates \
    -c HYRIV_ID -C data/test_ensemble/test.pointCoords.txt \
    -o test -t 4
```

To instead create a model-average of only the best models sampled in each run, you can add the `--only_best` argument:

```
ensembleResistnet.py \
    -n data/test_ensemble/test.full.network \
    -i data/test_ensemble/replicates \
    -c HYRIV_ID -C data/test_ensemble/test.pointCoords.txt \
    -o test -t 4 --only_best
```

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
