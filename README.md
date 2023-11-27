# ResistNet
Software for examining spatial patterns of diversity and differentiation in dendritic ecological networks (=DENs). Note that the contents of this directory are a work in progress.

## Table of Contents:
1. [Installation](#install)
    1. [Installing from GitHub](#install_dev)



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



## ResistNet <a name="rscape"></a>
### Program description <a name="rscape_desc"></a>


#### Genetic algorithm background and similarities to ResistNet

All variables can vary in their model inclusion, but also in the manner by which they are transformed; transformations are the same as those used in ResistNet, and shape parameters are included for each variable in the "chromosome". Note that as shape parameters increase, the more similar the "transformed" values are to the original values, with high values for the shape parameter (e.g., 50-100) resulting in negligable transformation:

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/transforms.png)

### Usage <a name="rscape_usage"></a>

#### Options and Help Menu <a name="rscape_help"></a>

As with all of the scripts in this repository, you can view a list of the options for ResistNet by calling the help menu with <-h>:
```
resistnet.py

Author: Tyler K Chafin, Biomathematics and Statistics Scotland
Contact: tyler.chafin@bioss.ac.uk
Description: Genetic algorithm to optimize resistance models on networks

Input options:
  If using FitDistNet outputs:
	-p,--prefix	: Prefix for autoStreamTree outputs

-or-

  If manually specifying inputs:
	-g,--genmat	: Genetic distance matrix
	-n,--network	: Input graph (in pickle'd networkx format)
	<add the rest later>

General options:
	-s,--seed	: Random number seed (default=taken from clock time)
	-T,--procs	: Number of parallel processors
	-X,--noPlot	: Turn off plotting
	-o,--out	: Output file prefix
	-h,--help	: Displays help menu

Genetic Algorithm Options:
	-P,--maxPop	: Maximim population size [default = 100]
	-G,--maxGen	: Maximum number of generations [default = 500]
	-s,--size	: Manually set population size to <-p int>
			    NOTE: By default, #params * 15
	-m,--mutpb	: Probability of mutation per individual [default=0.2]
	-i,--indpb	: Probability of mutation per trait [default=0.1]
	-c,--cxpb	: Probability of being chosen for cross-over [default=0.5]
	-t,--tourn	: Tournament size [default=10]
	--posWeight	: Constrain parameter weights to between 0.0-1.0
	--fixWeight	: Constrain parameter weights to 1.0 (i.e., unweighted)
	--allShapes	: Allow inverse and reverse transformations

Model optimization/ selection options:
	-F,--nfail	: Number of generations failing to improve to stop optimization
	-d,--delt	: Threshold absolute change in fitness [default=0.0]
	-D,--deltP	: Threshold percentage change in fitness, as decimal [default=0.001]
	-f,--fit	: Fitness metric used to evaluate models <NOT IMPLEMENTED>
			    Options:
			    aic (default)
			    loglik (log-likelihood)
			    r2m (marginal R^2)
			    delta (Change in AIC versus null model)
			    NOTE: Case-insensitive
	-b,--burn	: Number of generations for pre-burnin [default=0]
	--max_hof_size	: Maximum individuals to track in the Hall of Fame [default=100]

Circuitscape options:
	--cholmod	: Turn on CHOLMOD solver
	-C,--cprocs	: Processors per Circuitscape process [default=1]
			    NOTE: Total simultaneous processes= <-T> * <-C>

Genetic distance options:
	--force		: Use XX attribute from input table as distance metric (e.g. 'fittedD')
			    NOTE: By default, the average of "locD_" columns will be taken
	--infer		: Infer pairwise distances from input table (i.e., NOT input matrix)
	-v,--vars	: Comma-separated list (no spaces) of explanatory attributes to include

Multi-model inference options:
	-A,--modavg	: Compute model-averaged resistances
			    NOTE: This involves re-running Circuitscape for each model
	-a,--awsum	: Cumulative Akaike weight threshold to retain top N models [default=0.95]
	--report_all	: Plot per-stream resistance and generate full outputs for all retained models
```

#### Input files <a name="rscape_input"></a>

For convienience, the inputs for ResistNet follow the formats of the files output by DistNet and FormatNet (not complete yet).

#### Outputs <a name="rscape_output"></a>

After parsing all of the inputs, ResistNet will randomly generate a population of 'individuals' (=model parameterizations), which is by default 15X the number of parameter, up to a maximum size specified by <-P,--maxPop>. Alternatively, you can specify a fixed size with <-s,--size>. Each 'generation', individuals will be selected using the desired fitness function (specified with <-f,--fit>; e.g. -f AIC to use AIC), and the maximum, minimum, meand, and standard deviation of population fitness values will be output to the terminal (stdout):
```
Reading network from:  ../out3.network
Reading DistNet results from: ../out3.streamTree.txt
Initializing genetic algorithm parameters...

Establishing a population of size: 50

Evaluating initial population...

Starting optimization...

-- Generation 1 --
  Worst 7466.655283180244
  Best 6644.0304460748275
  Avg 6986.83389341895
  Std 174.1935275523536
  nFails 0
-- Generation 2 --
  Worst 7499.103994724142
  Best 6637.266874648363
  Avg 6821.94687814646
  Std 149.476198344529
  nFails 0
...
...
```
Also included in this report is the number of consecutive generations that the genetic algorithm has failed to find a 'better' individual (with 'better' being defined using thresholds set with <-d,--delt> or <-D,--deltP>). After either a specified number of generations (<-G,--maxGen>) have passed, or nFail exceeds the user-specified number (-F,--nFail), the algorithm will stop, and report to the screen which stopping criteria was met:
```
-- Generation 5 --
  Worst 6854.294254647968
  Best 6596.9194357890165
  Avg 6643.346646257467
  Std 53.73712600565164
  nFails 0
Stopping optimization after 5 generations.
Reason: Exceeded maxGens
```
At this time a series of plots and tables will be produced. The fitness values through time will be reported at $out.FitnessLog.tsv (where "$out" is the output prefix defined with <-o,--out>):

| Generation | Worst             | Best              | Mean              | Stdev              |
|------------|-------------------|-------------------|-------------------|--------------------|
| 1          | 6644.389928164675 | 7466.655283180244 | 6985.505559085744 | 173.72113685174205 |
| 2          | 6639.881724568676 | 7499.539461868835 | 6820.890653185437 | 147.67937472826992 |
| 3          | 6617.603822325212 | 6872.074421364015 | 6683.739424517294 | 49.76330989463599  |
| 4          | 6610.290247013926 | 6794.057145496889 | 6642.676944232899 | 29.918610153913338 |
| 5          | 6597.775539690673 | 7018.279668244255 | 6642.276324199858 | 67.76930180708179  |

In addition to this, ResistNet maintains a "Hall of Fame", tracking the best individuals to have ever be sampled in the population, up to a maximum Hall of Fame size specified with <--max_hof_size>, which defaults to 100. This table will be sorted by fitness value, and will include calculations of the Akaike weight. Models are assigned to be retained for model averaging using the cumulative Akaike threshold set by the user <-a,--awsum>, which defaults to 0.95.
| fitness            | run_mm_cyr | run_mm_cyr_weight   | run_mm_cyr_trans | run_mm_cyr_shape | kar_pc_cse | kar_pc_cse_weight | kar_pc_cse_trans | kar_pc_cse_shape | soc_th_cav | soc_th_cav_weight | soc_th_cav_trans | soc_th_cav_shape | pst_pc_cse | pst_pc_cse_weight    | pst_pc_cse_trans | pst_pc_cse_shape | riv_tc_csu | riv_tc_csu_weight   | riv_tc_csu_trans | riv_tc_csu_shape | LENGTH_KM | LENGTH_KM_weight    | LENGTH_KM_trans | LENGTH_KM_shape | slp_dg_cav | slp_dg_cav_weight   | slp_dg_cav_trans | slp_dg_cav_shape | urb_pc_cse | urb_pc_cse_weight   | urb_pc_cse_trans | urb_pc_cse_shape | CSI | CSI_weight          | CSI_trans | CSI_shape | for_pc_cse | for_pc_cse_weight   | for_pc_cse_trans | for_pc_cse_shape | loglik              | r2m                 | aic                | delta_aic_null    | delta_aic_best     | akaike_weight          | acc_akaike_weight  | keep  |
|--------------------|------------|---------------------|------------------|------------------|------------|-------------------|------------------|------------------|------------|-------------------|------------------|------------------|------------|----------------------|------------------|------------------|------------|---------------------|------------------|------------------|-----------|---------------------|-----------------|-----------------|------------|---------------------|------------------|------------------|------------|---------------------|------------------|------------------|-----|---------------------|-----------|-----------|------------|---------------------|------------------|------------------|---------------------|---------------------|--------------------|-------------------|--------------------|------------------------|--------------------|-------|
| 6597.775539690673  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 0          | -                   | -                | -                | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | -0.9319860657625789 | 3         | 64        | 1          | -0.6901645609234912 | 2                | 89               | -3294.8877698453366 | 0.1945129136531144  | 6597.775539690673  | 632.5800876388339 | 0.0                | 0.5064833775653471     | 0.5064833775653471 | True  |
| 6598.624189654391  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 0          | -                   | -                | -                | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | -0.9319860657625789 | 3         | 64        | 1          | 0.5971939266873523  | 4                | 18               | -3295.3120948271953 | 0.1952412393909014  | 6598.624189654391  | 635.0907078516702 | 0.8486499637174347 | 0.3313471187668753     | 0.8378304963322224 | True  |
| 6600.142494390143  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 1          | 0.6205977178035142  | 2                | 35               | 1         | -0.9664157414652275 | 8               | 39              | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | -0.9319860657625789 | 3         | 64        | 1          | -0.6901645609234912 | 2                | 89               | -3296.0712471950715 | 0.18871059961215467 | 6600.142494390143  | 612.7956484064507 | 2.3669546994697157 | 0.15509132782604312    | 0.9929218241582656 | True  |
| 6608.1073678991415 | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 1          | 0.6205977178035142  | 2                | 35               | 1         | -0.9664157414652275 | 8               | 39              | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | 0.6465553956747097  | 0         | 75        | 1          | -0.6901645609234912 | 2                | 18               | -3300.0536839495708 | 0.18729419553514634 | 6608.1073678991415 | 607.3437218603012 | 10.331828208468323 | 0.0028909275456527952  | 0.9958127517039184 | False |
| 6609.615347179055  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 0          | -                   | -                | -                | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | 0.6465553956747097  | 0         | 75        | 1          | -0.6901645609234912 | 2                | 18               | -3300.8076735895274 | 0.19483502307914924 | 6609.615347179055  | 631.8005512099589 | 11.839807488381666 | 0.001360140170154527   | 0.9971728918740729 | False |
| 6610.290247013926  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 1          | 0.6205977178035142  | 2                | 35               | 1         | -0.9664157414652275 | 8               | 39              | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | 0.6465553956747097  | 0         | 75        | 1          | -0.6901645609234912 | 2                | 18               | -3301.145123506963  | 0.1852952541184347  | 6610.290247013926  | 600.6581119492184 | 12.514707323252878 | 0.0009705793121750887  | 0.9981434711862479 | False |
| 6611.183034947682  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 0          | -                   | -                | -                | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | 0.6465553956747097  | 0         | 75        | 1          | -0.6901645609234912 | 2                | 89               | -3301.591517473841  | 0.19179945737685386 | 6611.183034947682  | 621.7468721555406 | 13.407495257009032 | 0.0006211043823916513  | 0.9987645755686395 | False |
| 6611.498278064627  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 0          | -                   | -                | -                | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | 0.6465553956747097  | 0         | 75        | 1          | -0.6901645609234912 | 2                | 89               | -3301.7491390323135 | 0.19479738590965015 | 6611.498278064627  | 631.4176360681067 | 13.722738373953689 | 0.0005305305783013392  | 0.9992951061469408 | False |
| 6611.498278064627  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 0          | -                   | -                | -                | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | 0.6465553956747097  | 0         | 75        | 1          | -0.6901645609234912 | 2                | 89               | -3301.7491390323135 | 0.19479738590965015 | 6611.498278064627  | 631.4176360681067 | 13.722738373953689 | 0.0005305305783013392  | 0.9998256367252422 | False |
| 6614.424384033492  | 1          | -0.1455520926262841 | 6                | 77               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 1          | -0.09796388256653743 | 0                | 80               | 1          | 0.6205977178035142  | 2                | 35               | 1         | -0.9664157414652275 | 8               | 39              | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | -0.9319860657625789 | 3         | 64        | 1          | -0.6901645609234912 | 2                | 89               | -3303.212192016746  | 0.18791714272102264 | 6614.424384033492  | 608.633840904944  | 16.648844342818848 | 0.00012283286528141403 | 0.9999484695905235 | False |
| 6617.603822325212  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 1          | 0.07949635894691087 | 7                | 18               | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | -0.9319860657625789 | 3         | 64        | 1          | 0.5971939266873523  | 4                | 18               | -3304.801911162606  | 0.20153290148171848 | 6617.603822325212  | 653.0208456285163 | 19.828282634538482 | 2.5055803231555324e-05 | 0.9999735253937551 | False |
| 6619.827623973979  | 1          | -0.1455520926262841 | 6                | 11               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 0          | -                   | -                | -                | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | -0.9319860657625789 | 3         | 64        | 1          | 0.5971939266873523  | 4                | 18               | -3305.9138119869895 | 0.19816589207605215 | 6619.827623973979  | 641.6683238153164 | 22.052084283305703 | 8.241683590176287e-06  | 0.9999817670773452 | False |
| 6620.788954854979  | 1          | -0.1455520926262841 | 6                | 77               | 0          | -                 | -                | -                | 0          | -                 | -                | -                | 0          | -                    | -                | -                | 1          | 0.07949635894691087 | 7                | 18               | 0         | -                   | -               | -               | 1          | 0.34392503912199834 | 2                | 51               | 1          | -0.7752412475467703 | 8                | 60               | 1   | -0.9319860657625789 | 3         | 64        | 1          | -0.6901645609234912 | 2                | 89               | -3306.3944774274896 | 0.19844574255201924 | 6620.788954854979  | 642.5574955984894 | 23.013415164306025 | 5.096424430273501e-06  | 0.9999868635017755 | False |

Other values in the table include whether or not a variable is included in a model ($var, where $var is the variable name), specified as either "1" (=included) or "0" (=excluded), and the transformation type (=$var_trans), transformation shape parameter (=$var_shape), and weight of the parameter when calculating the composite resistance edges (=$var_weight). For the transformation column, the types of transformations are as follows: 0=Not transformed; 1=Ricker; 2=Reverse Ricker; 3=Inverse Ricker; 4=Reverse-Inverse Ricker; 5=Monomolecular; 6=Reverse Monomolecular; 7=Inverse Monomolecular; 8=Reverse-Inverse Monomolecular. In all cases, the larger the shape value, the closer each transformation gets to being linear (=essentially no transformation).

A plot summarizing among models called $out.ICprofile.pdf will also be produced, which shows several pieces of information: 1) How AIC supports vary among all of the Hall of Fame models (arranged from left to right = 'best' to 'worst'). Points are separated as those which were retained for model-averaging given the <-a,--awsum> cutoff, scaled by marginal R^2 values (i.e., correlation coefficient from the MLPE model), and with a red horizontal bar showing a raw delta-AIC cutoff of 2:

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/ic_profile.png)

Another plot will summarize the relationship among the various possible fitness metrics (marginal R^2, AIC, etc), as well as how these are distributed among the retained and excluded models, $out.pairPlot.pdf:

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/pairplot.png)

Relative variable importance values will be plotted as a bar plot, $out.varImportance.pdf:

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/vif.png)

If model-averaging was performed (turned on using <-A,--modavg>), ResistNet will calculate model-averaged resistance values, which will be plotted by coloring edges in the stream network, $out.Model-Average.streamsByResistance.pdf:

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/modavg_resist.png)

It will also produce plots with a simple linear regression of genetic distances against pairwise resistance ($out.Model-Average.Pairwise.pdf) and edge-wise fitted distances against effective resistances from Circuitscape ($out.Model-Average.Edgewise.pdf):

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/modavg_pw.png)

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/modavg_ew.png)

If using the <--report_all> option, these plots will also be produced for every one of the "kept" models from the Hall of Fame, with naming as $out.Model-#.Pairwise.pdf (etc), where "#" is the row number from the HallofFame.tsv table (with 0-based indexing; i.e. 0="best" model; 1=second-best, and so on).

#### Parallel execution <a name="rscape_parallel"></a>

There are two parameters for controlling parallel execution: <-T,--procs> controls the spread of individual fitness calculations in the genetic algorithm, using the multiprocessing Python library; and <-C,--cprocs> controls the number of processor cores dedicated per Circuitscape run (which is run for every individual, every generation). <-C> * <-T> should not exceed the number of CPUs on your machine.

I've found that each Circuitscape run generally doesn't take too long with moderately sized networks (more below in the [Runtimes and Benchmarking section](#rscape_benchmark)), and that the greatest gains can be gotten by maximizing the <-T> parameter. This allows not only the Circuitscape step to be parallelized, but also the generation and parsing of Circuitscape outputs, and fitting of the MLPE model. However, extremely large networks may benefit from increasing the <-C> parameter. For most users, I wouldn't recommend increasing <-C> unless your machine has twice the number of processors as there are individuals in the GA population, as <-T> will not give any benefit if increasing beyond the population size. For example, if you have a small population size of 24 and are running on a 48-core machine, increasing <-T> beyond T=24 will not increase runtimes, and you may wish to set <-T 24> and <-C 2> in order to best use the cores avaialable. However, for most use cases, the available number of CPUs is unlikely to be larger than the population size.

#### Genetic Algorithm options <a name="rscape_ga"></a>

ResistNet provides options for manipulating the relevant parameters of the genetic algorithm. If there are additional parameters you wish to control, just shoot me an email (tylerkchafin@gmail.com) or launch an Issue here and GitHub and I'll add it as soon as I can.

The parameters which can be manipulating from the command-line are as follows:
```
	Genetic Algorithm Options:
		-P,--maxPop	: Maximim population size [default = 100]
		-G,--maxGen	: Maximum number of generations [default = 500]
		-s,--size	: Manually set population size to <-p int>
				    NOTE: By default, #params * 15
		-m,--mutpb	: Probability of mutation per individual [default=0.2]
		-i,--indpb	: Probability of mutation per trait [default=0.1]
		-c,--cxpb	: Probability of being chosen for cross-over [default=0.5]
		-t,--tourn	: Tournament size [default=10]
		--posWeight	: Constrain parameter weights to between 0.0-1.0
		--fixWeight	: Constrain parameter weights to 1.0 (i.e., unweighted)
		--allShapes	: Allow inverse and reverse transformations
```

The --posWeight and --fixWeight options are used to either constrain parameter weights to 0.0 - 1.0, or to fix all weights at 1.0. By default, weights can vary from -1.0 to 1.0, and all 'negatively weighted' transformations are not available (inverse and reverse). Only ricker, monomolecular, and inverse-reverse versions of both are available unless --allShapes is used. I don't recommend using both negative weights AND --allShapes together, because this creates a situation where there are multiple ways to get the same possible result (i.e., ricker with -1.0 weighting has the same impact on composite surface as reverse ricker with 1.0 weighting).

#### Model Selection/ Optimization <a name="rscape_model"></a>

#### Circuitscape <a name="rscape_cs"></a>

#### Model-averaging and multi-model importance <a name="rscape_modavg"></a>

### Runtimes and Benchmarking <a name="rscape_benchmark"></a>
In all of these tests, I used a dataset composed of N=112 network edges (corresponding to N=2,248 contiugous segments in the input shapefile), and N=78 populations. I randomly selected 10 environmental variables, and used a population size of 50. Tests were performed on a Linux computer with 256gb of memory and 16 cores.

Increasing the <-T> parameter leads to speed increases when <-t> <= population size, although returns are diminishing at higher values, given several added costs: The need for more separate initializations, and the communication cost of sending information between the master and sub-processes. Here, for each generation, the master process has to send a model parameterization (a vector of e.g., [1, 0.89, 7, 16] x number of variables), and the sub-process returns a vector of potential fitness measures (loglik, AIC, marginal r^2, deltaAIC vs. null). Memory increases linearly as you add processors, with the slope determined by the size of the input dataset, since each sub-process maintains a deep copy of the entire graph, pairwise genetic distance matrix, etc.

In this test, I used 1-16 processors to examine N=50 individual models per generation, for 10 generations. The random number seed was fixed, to maintain comparability. Total runtimes varied from 165 minutes for <-T 1> to ~14 minutes for <-T 16>. Runtime increases linearly with number of generations, such that e.g., 200 generations with this same dataset and population size takes ~280 minutes with <-T 16>. A run of the same dataset, with a population size of 200, and 200 generations, took 840 minutes.

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/rscape_benchmark.png)

In comparing the gains from parallelization, the <-T> parameter grants the most efficient use of available CPUs:

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/t_vs_c.png)

Here, dedicating all 16 cores to <-T> parallelization, after ~3 minutes spend parsing the input files, each generation took approximately 1 minute of computation. With a population size of 50, this averages out to each CPU core evaluating 3-4 individual models, so about 15-20 seconds per individual model (that includes running Circuitscape, parsing the outputs, and fitting the MLPE model). When setting <-T 8> and <-C 2>, each model evaluation averages to 20-24 seconds each, AND each CPU had to evaluate 6-7 individuals per generation. So not only did we see loss in the number of simultaneous model evaluations at a time, Circuitscape actually got slower -- probably reflecting the cost of passing around data. With much larger networks and a larger number of pairwise evaluations (maybe several hundred samples), this relationship is likely to change. For large datasets with long computation times, I would recommend doing some short runs evaluating whether or not increasing the <-C> parameter is worthwhile for your case.

### Example Workflows <a name="rscape_workflow"></a>

## Scripts and Tools

#### Re-plotting StreamTree outputs

The default $out.streamdByFittedD plot may not be exactly what you wanted. To prevent cluttering the help menu of the main program too much, we've provided a separate script for loading up autoStreamTree outputs to re-make the plot, which has some added options for customization: scripts/plotStreamTree.py

'''
$ python3 ./scripts/plotStreamTree -h
plotStreamTree.py

Author: Tyler K Chafin, University of Arkansas
Contact: tkchafin@uark.edu
Description: Script for re-plotting StreamTree results after running DistNet

		-p,--prefix	: Prefix for autoStreamTree output
		-m,--min	: Minimum genetic distance
		-M,--max	: Maximum genetic distance
		-c,--cmap	: Colormap (any valit matplotlib cmap value)
			see: https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html
		-o,--out	: Output prefix (if not overwriting original)
		-h,--help	: Displays this help menu)
'''

For example, to re-plot values with a distance ceiling of 0.2 and a viridis color scale:
'''
python3 ./scripts/plotStreamTree.py -p out2 -m 0.0 -M 0.2 -c "viridis"
'''

#### Clustering populations using any distance matrix

For example, if you wanted to cluster individuals using their stream distances, I've provided a script called clusterPopsDB.py which will use a DBSCAN clustering algorithm to output a tab-delimited population map given any arbitrary distance matrix:
'''
$ python3 ./scripts/clusterPopsDB.py -h
'''
