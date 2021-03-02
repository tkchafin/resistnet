# DENdriscape
A collection of software tools for examining spatial patterns of diversity and differentiation in dendritic ecological networks (=DENs). Note that the contents of this directory are a work in progress. 

### Table of Contents: 
1. [Installation](#installation)
    1. [Installation with conda](#conda)
    2. [Installation with pyenv](#pyenv)
    3. [Troubleshooting PyJulia](#pyjulia)
    4. [Singlurity/ Docker](#sing)
2. [Summary of Programs](#programs)
3. [FitDistNet - Fitting distances to stream networks](#ast)
    1. [Program Description](#ast_desc)
        1. [StreamTree Background](#ast_background)
    2. [Usage](#usage)
        1. [Options and Help Menu](#ast_help)
        2. [Input File](#ast_table)
        3. [Shapefile format](#ast_input)
        4. [Output files](#ast_output)
        5. [Genetic Distance Methods](#gen)
        6. [Defining Populations](#pops)
    3. [Example Workflows](#ast_workflow)
        1. [Single-marker dataset](#ast_example1)
        2. [Large SNP datasets](#ast_example2)
        3. [Microsatellites](#ast_example3)
    4. [Runtimes and Benchmarking](#ast_benchmark)
    5. [References](#ast_refs)
4. [ResistNet - Optimizing effective resistance networks](#rscape)
    1. [Program Description](#rscape_desc)
        1. [Genetic Algorithms](#rscape_background)
        2. [Resistance Models](#rscape_background2)
    2. [Usage](#rscape_usage)
        1. [Options and Help Menu](#rscape_help)
        2. [Input Files](#rscape_inputs)
        3. [Output files](#rscape_output)
        4. [Parallel execution](#rscape_parallel)
        5. [Genetic Algorithm Options](#rscape_ga) 
        6. [Model Selection/ Optimization](#rscape_model)
        7. [Circuitscape Options](#rscape_cs)
        8. [Model-averaging and multi-model importance](#rscape_modavg)
    3. [Runtimes and Benchmarking](#rscape_benchmark)
    4. [Example Workflows](#rscape_workflow)
5. [Scripts and Other Tools](#tools)


### Installation <a name="installation"></a>

Because of the number of dependencies, I recommend setting up a virtual environment to prevent any conflicts with your system Python environment. 

If you are planning on using ResistNet, [PyJulia](https://pyjulia.readthedocs.io/en/latest/), which forms the necessary interface for Python to access Circuitscape (a Julia program), a complication is that PyJulia cannot use a Python distribution that is statically linked to libpython (such as that installed by Ubuntu or conda, or Debian-based Linux distributions). 

To check if this is a problem, you can use the following commands on Linux and Mac. If nothing prints to the screen, your Python is statically linked and will require a workaround:
```
#get full path to python3
which python3

#on linux
ldd </full/path/to/python3 | grep libpython

#on Mac
otools -L </full/path/to/python3> | grep libpython
```

If something along the lines of "libpython3.7m.so.1.0 => /usr/lib/libpython3.7m.so.1.0 (0x00007f10ef116000)" printed out, you are good to go to install the dependencies using the pip or conda instructions below (depending on your system). 

These installation instructions should work in Linux (or Windows Linux subsystem) or Mac, although Mac users can substitute homebrew commands where suggested. Note that these instructions are also comprehensive, including steps such as compiling and installing R from source -- many users will likely already have R installed on their system, and thus can skip some steps. 

#### Installation using conda <a name="conda"></a>

First, you need to set up a virtual environment:
```
conda create -n pyjulia conda-forge::python=3.6.6
```

Next, activate your environment and install dependencies:

```
#activate environment
conda activate pyjulia

#install dependencies
conda install -c conda-forge -c bioconda -c rmg julia pyjulia  numpy scipy networkx seaborn matplotlib pandas deap sortedcontainers julia geopy geopandas shapely scikit-learn r-base=4.0.3 r-essentials r-nloptr rpy2
```

Open Julia, and install the packages that we will be needing:

```
julia
julia> using Pkg
julia> Pkg.add("PyCall")
julia> Pkg.add("Circuitscape")
julia> exit()
```

If your Python is statically linked (see above ldd command), you will also need to create a custom Julia system image that you will pass to ResistNet:

```
python3 -m julia.sysimage sys.so
```


#### Installation using pyenv and manually compiled Julia/R <a name="pyenv"></a>

First, clone the repository (if on Mac, you can also use [homebrew](https://brew.sh/):
```
#cloning the repository
git clone https://github.com/pyenv/pyenv.git ~/.pyenv
#Homebrew users can install like so:
#brew install pyenv
```

Next, configure the environment:

```
#shouldn't be necessary for homebrew users
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo -e 'if command -v pyenv 1>/dev/null 2>&1; then\n eval "$(pyenv init -)"\nfi' >> ~/.bashrc
source ~/.bashrc
```

Verify the installation (this should output a list of available python versions):
```
pyenv install --list
```

Next, build python and install dependencies

```
#build python from scratch, enabling dynamic linking (needed for PyJulia to work)
PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.6.6
pyenv global 3.6.6

#verify Python is not statically linked to libpython
ldd ~/.pyenv/versions/3.6.6/bin/python3.6 | grep libpython
#if nothing prints, something didn't work

#check that the right python is in your path:
which python3 
which pip3
#output should be something like /home/USER/.pyenv/*/python3; if not, either add the correct path to .bashrc or use absolute path in the python3 and pip3 calls below

#upgrade pip and make sure venv is installed 
pip3 install --upgrade pip
```

After Python finishes building, we can set up our python virtual environment:
 
```
#cd to the directory you want to build your virtual environment in. Here I will just use the home directory 
cd ~

#make a directory to store python virtual environments in:
mkdir python_venv
cd python_venv

#set up a virtual environment called 'riverscape':
python3 -m venv riverscape

#activate the riverscape venv so we can install packaged into it:
source ./riverscape/bin/activate
```

If you do not have the [R](https://www.r-project.org/) and [Julia](https://julialang.org/) programming languages installed on you can install them directly into your new virtual environment like so. If on Mac you could use homebrew, or if on a Linux machine with root privileges you could use apt-get, but for the sake of completion I will here show how to compile them from source.

```
cd ~/python_venv/riverscape
#go to https://julialang.org/downloads/ if you want a pre-compiled binary
#NOTE: On a CentOS Linux HPC, I could only get the LTS v1.0.5 to work
git clone git://github.com/JuliaLang/julia.git
cd julia
make
ln -s ~/python_venv/riverscape/julia/bin/julia ~/python_venv/riverscape/bin/.

#follow a similar set up for R, selecting a binary or source from https://cran.r-project.org/ depending on your system
#here, I will be compiling it from source 
wget https://cran.r-project.org/src/base/R-4/R-4.0.2.tar.gz
tar -xvzf R-4.0.2.tar.gz
cd R-4.0.2
./configure --prefix=/home/tkchafin/python_venv/riverscape/bin
make
make install
```

Open Julia, and install the packages that we will be needing:

```
julia
julia> using Pkg
julia> ENV["PYTHON"] = "/home/tkchafin/python_venv/riverscape/bin/python3" #whatever your path is
julia> Pkg.add("PyCall")
julia> Pkg.add("Circuitscape")
julia> exit()
```

Finally, install the required Python dependencies:

```
#install dependencies
pip3 install --upgrade pip
pip3 install numpy scipy networkx seaborn matplotlib pandas deap sortedcontainers julia geopy geopandas shapely scikit-learn rpy2

#installing rpy2 (Python-R interface) on Mac required that I provide a path to gcc compiler:
#which gcc 
#this will output the path you need to set to 'CC' below
#env CC=/usr/bin/gcc pip3 install rpy2
```

Install the necessary R packages:

```
#check R path is the right one
which R
#should be: ~/python_venv/riverscape/bin/R
R
#install R packages
> install.packages("MuMIn")
> install.packages("Matrix")
> install.packages("lme4")
```

#### Troubleshooting PyJulia <a name="pyjulia"></a>

As noted above, PyJulia won't work with statically-linked Python. There are several workarounds in the PyJulia [documentation](https://pyjulia.readthedocs.io/en/latest/troubleshooting.html). 

One of the easiest workarounds is to turn off the Julia compiled cache, which you can do by passing ResistNet the '--no_compiled_modules' (boolean) argument. This might slow down loading Julia a little bit, but is one of the fastest ways to get up and running if you are getting PyJulia errors. 

Another option is to create a custom Julia system image, which can be passed to ResistNet using the --sys_image argument:

```
python3 -m julia.sysimage sys.so
```

#### Using Singularity/ Docker <a name="sing"></a>
Coming soon... I will at some point generate a Singularity image for the package -- if you need it in a hurry open and Issue so I will know to prioritize it

### Programs <a name="programs"></a>
- *FitDistNet.py* : This package calculates genetic and stream-distances by snapping sampling points to a network provided using an input shapefile and calculating a minimal sub-network, tests for isolation-by-distance, and fits genetic distances to stream segments using the Stream-tree algorithm by Kalinowski et al. 2008. 
- *ResistNet.py* : This package implements a genetic algorithm to optimize multi-variable resistance models on networks in Circuitscape. This is somewhat similar to the way in which ResistanceGA accomplishes this for raster datasets, but with some critical algorithmic differences (see below)
- *FormatNet.py* : (Coming soon) Integrates the network extraction functions of autoStreamTree with various pre-processing steps such as joining non-contiguous segments, merging redundant paths (e.g., for a braided stream) -- Not yet added
- *SimNet.py* : (Coming soon) Performs forward-simulations on a given stream network, given environmental 'resistance' variables (and optional weights) -- Not yet added
- *BGR_pipeline.py* : (Coming soon) Pipeline for applying BGR model given a set or covariates or resistance model
- tools/ 
	* *autoFetcher.py* : Python interface for NCBI entrez API; can be used to automatically download and parse queries
	* *clusterPopsDB.py* : Distance-based clustering of points using the DBSCAN algorithm
	* *fasta2phylip.py* : Converting between FASTA and PHYLIP formats
	* *fasta2table.py* : Converting between FASTA and autoStreamTree Table formats
	* *nremover.pl* : Filtering FASTA and PHYLIP files of concatenated SNPs
	* *plotStreamTree.py* : Re-making FitDistNet fitted distance plot with some additional options
	* *modelParser.py* : Parse model outputs of RiverscapeGA.py to generate and select transformed datasets (Coming soon)
	* [tkchafin/scripts](https://github.com/tkchafin/scripts) : For more useful file formatting and filtering scripts
	* [stevemussmann/Streamtree_arcpy](https://github.com/stevemussmann/StreamTree_arcpy) : Code for generating inputs for the original StreamTree program by Kalinowski et al.

## FitDistNet <a name="ast"></a>

### Software Description <a name="ast_desc"></a>
FitDistNet is a Python software package providing various analyses aimed at analyzing patterns of genetic differentiation among aquatic stream-dwelling organisms. The intention is to take what was previously a tedious process involving multiple discrete steps and to integrate these all in one place. 

Currently, FitDistNet provides a companion library of functions for calculating various measures of genetic distances among individuals or populations, including model-corrected p-distances (e.g. Jukes-Cantor 1969, Kimura 2-parameter, Tamura-Nei 1993) as well as those based on allele frequencies (e.g. Theta-ST, linearized Fst, Jost's D -- full list of available distance models below). It also includes integrated functions for parsing an input vector shapefile of streams (see below 'Requirements for input shapefiles') for easy calculation of pairwise stream distances between sites, as well as the ordinary or weighted least-squares fitting of reach-wise genetic distances according to the "stream tree" model of Kalinowski et al. (2008). Various plotting functions are also provided for downstream analysis, including looking at patterns of isolation-by-distance. Outputs should also be directly importable into R, with additional outputs with annotated streamtree fitted distances provided for analysis in your GIS suite of choice. 

If you use this package for analysis of fitted distances using the streamtree model, please cite the following:
* Kalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760).

#### StreamTree background <a name="ast_background"></a>

### Usage <a name="ast_usage"></a>

#### Options and Help Menu <a name="ast_help"></a>

To view all of the options for FitDistNet, call the program with the <-h> argument:
```
$ python3 FitDistNet.py -h

Exiting because help menu was called.

FitDistNet.py

Author: Tyler K Chafin, University of Arkansas
Contact: tkchafin@uark.edu
Description: Computes stream distances and genetic distances for georeferenced DNA sequences, performs tests for isolation-by-distance, and uses a least-squares method to fit distances to stream segments.

	Input file format:
		SampleName	Data	Lat	Long	seq1	[seq2]...[seqn]
		...
		...
	--NOTE: The "DATA" column can be anything- a population/ species identifier
	  (e.g. when used with --pop), or irrelevant data (e.g. GenBank accesson
	  number, if datafile produced by my autoFetcher script)

	Mandatory arguments:
		-s,--shp	: Path to shapefile containing cleaned, contiguous stream reaches
		-i,--input	: Input .tsv file containing sample coordinates and sequences

	General options:
		-o,--out	: Output prefix [default="out"]
		-n,--network	: Provide an already optimized network output from a previous run
			This will be the $out.network file written by autoStreamTree
		--overwrite	: Overwrite an input network (Only relevant with --network)
		-h,--help	: Displays help menu
		-r,--run	: Run which steps? Options: [all, gendist, ibd, streamdist, streamtree]
			ALL		: Run all steps
			GENDIST		: Only calculate genetic distance matrix
			STREAMDIST	: Only compute pairwise stream distances
			DISTANCES	: Only compute GENDIST + STREAMDIST
			IBD		: GENDIST + STREAMDIST + Mantel test
			STREAMTREE	: GENDIST + STREAMDIST + fit StreamTree model
			RUNLOCI	: Run STREAMTREE fitting on each locus
		-p,--pop		: Pool individuals based on column 2 of input file
			NOTE: The location will be taken as the centroid among individual samples
		-g,--geopop		: Pool individuals having identical coordinates
		-c,--clusterpop	: Use DBSCAN algorithm to automatically cluster populations
		--reachid_col	: Attribute name representing primary key in shapefile [default="REACH_ID"]
		--length_col	: Attribute name giving length in kilometers [default="LENGTH_KM"]

	Genetic distance options:
		-d,--dist	: Use which metric of distance? Options are:
			Substitution models (individual-based):
			  PDIST			: Uncorrected p-distances [# Differences / Length]
			  JC69 			: [default] Jukes-Cantor (1969) corrected p-distances
			  K2P			: Kimura 2-parameter distances
			  TN84			: Tajima and Nei's (1984) distance
			  TN93			: Tamura and Nei's (1993) distance
			Frequency models (when using --pop):
			  FST			: Weir and Cockerham's Fst formulation (=THETAst)
			  GST			: Hedrick's (2005) correction of Nei (1987) Gst [=G'st]
			  GSTPRIME		: Meirmans & Hedrick (2011) corrected G'st [=G''st]
			  LINFST		: [default] Rousset's (1997) Fst [=Fst/(1-Fst)]
			  JOST			: Jost's (2008) D
			  NEI72			: Nei's (1972) standard genetic distance 
			  NEI83			: Nei and Chesser (1983) Da
			  EUCLID		: Euclidean distance
			  CHORD			: Cavalli-Sforza and Edwards (1967) chord distance
			  --NOTE: Individual-based metrics can also be computed for
		  	          populations. You can set how these are aggregated w/ --pop_agg
			  --NOTE: Multiple loci for PDIST, JC69, K2P, and EUCLID distances
		  	        will be reported using the method defined in --loc_agg
			  --NOTE: TN84 will use empirical base frequencies
		-G,--genmat	: Skip calculation and use the provided labeled .tsv matrix
		--coercemat	: [Boolean] Coerce negative values in input matrix to zero
		--locmatdir	: Directory of per-locus distance matrices
		--het		: [Boolean] Count partial differences [e.g. ind1=T, ind2=W]
		--snp		: [Boolean] Data represent concatenated SNPs
		--msat		: xxx[Boolean] Data represent msat alleles [not yet implemented]
		--global_het	: Estimate Ht using global frequencies (default is averaged over pops) 
	
	DBSCAN options (only when --clusterpop):
		--min_samples	: Minimum samples per cluster [default=1]
		--epsilon		: Maximum distance (in km) within a cluster [default=20]
		
	Aggregation options: 
		-P,--pop_agg	: Define aggregator function for certain genetic distances in pop samples
		-L,--loc_agg	: Define aggregator function for aggregating locus-wise distances
			All of these can take the following options:
			  ARITH		: [default] Use arithmetic mean
			  MEDIAN	: Use median distance
			  HARM		: Use harmonic mean
			  ADJHARM	: Adjusted harmonic mean (see docs)
			  GEOM		: Use geometric mean
			  MIN		: Use minimum distance
			  MAX		: Use maximum distance

	IBD options:
		--perm		: Number of permutations for mantel test [def=1000]
		--and_log	: Also perform IBD steps with log geographic distances

	StreamTree (see Kaliowski et al. 2008) options:
		--iterative	: Prevent negative distances using the iterative approach
		-w,--weight	: Desired weighting for least-squares fitting:
			Options:
			  FM67			: Fitch and Margoliash (1967) [w = 1/D^2]
			  BEYER74		: Beyer et al. (1974) weights [w = 1/D]
			  CSE67			: [default] Cavalli-Sforza and Edwards (1967) [w = 1]
```

#### Requirements for input shapefiles <a name="ast_input"></a>

I highly recommend using the existing global stream datasets provided by the [HydroLab group](https://wp.geog.mcgill.ca/hydrolab/) at McGill University, specifically the [HydroAtlas](https://www.hydrosheds.org/page/hydroatlas) or [free-flowing rivers dataset](https://wp.geog.mcgill.ca/hydrolab/free-flowing-rivers/) as these are already properly formatted for use, and the additional variables included will make downstream analysis very easy. Because of their size, I would recommend clipping them to the relevant scale first (e.g. the drainage encompassing all of your samples).

Note that a valid path is required between all sites in order to calculate pairwise stream distances. Thus, if you are analyzing data from multiple drainages which only share an oceanic connection, you will need to augment the shapefile. For example this could be accomplished by adding a vector representing the coastline to create an artificial connection among drainages. 

If for some reason you cannot use the HydroRIVERS dataset, you will need to do some things first before loading your shapefile into autoStreamTree. First, you will need to include two variables in the attribute table of your shapefile: 1) REACH_ID (case sensitive) must provide a unique identifier to each stream reach; and 2) LENGTH_KM should give the length of each segment. Next, because sometime large stream layers will have small gaps in between streams, you will need to span any small gaps between streams which should be contiguous, and also dissolve any lines that overlap with one another so that any given section of river is represented by a single line. I will provide a tutorial for doing this in ArcMAP later, but for now there are some scripts in our complementary package that can help with these steps using the ArcPy API: https://github.com/stevemussmann/StreamTree_arcpy. Note that this package will also help you in running the original Stream Tree package on Windows, if you want to do so. 

#### Input file format <a name="ast_table"></a>

Coordinates do not need to exactly match nodes in the input shapefile, as points will be 'snapped' to network nodes after parsing. autoStreamTree will output both a table($out.snapDistances.txt) and a histogram plot ($out.snapDistances.pdf) showing distances in kilometers that samples or populations had to be snapped:

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/example.snapDistances.png)

#### FitDistNet Outputs <a name="ast_outputs"></a>

The first thing autoStreamTree will do upon reading your input shapefile is to calculate a minimally reduced sub-network which collapses the input river network into continuous reaches (="edges"), with nodes either representing sample localities or junctions. Because the full river network will likely contain many branches and contiguous reaches which do not contain samples, these are removed to speed up computation. The underlying metadata will be preserved, and the final output will consist of an annotated shapefile containing an EDGE_ID attribute which tells you how reaches were dissolved into contiguous edges in the graph, and a FittedD attribute giving the least-squares optimized distances.

The reduced sub-network will be plotted for you in a file called $OUT.subGraph.pdf:
![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/example.subGraph.png)

Here, the total cumulative stream length (in km) is plotted along edges (NOTE: Any natural curvature in the river is not preserved in this plot), with sample sites as blue dots and junctions as black dots. A geographically accurate representation, coloring individual streams to designate different dissolved edges, will be provided as $out.streamsByEdgeID.pdf: 
![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/example.networkByEdgeID.png)

After fitting genetic distances, autoStreamTree will create several other outputs. First, a table called $out.reachToEdgeTable.txt will give a tab-delimited map of how REACH_ID attributes were dissolved into contiguous edges. Second, a tabular and graphical representation of how fitted pairwise distances compare to the raw calculates (or user-provided) pairwise distances: $out.obsVersusFittedD.txt and $out.obsVersusFittedD.pdf
![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/example.obsByFittedD.png)

Finally, the fitted distances per stream edge will be output both as an added column to the original shapefile attribute table ($out.streamTree.shp and $out.streamTree.txt), and also as a plot showing how distances compare across all streams:
![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/example.networkByStreamTree.png)

#### Genetic distance models <a name="gen"></a>

The currently recommended way to run autoStreamTree is to provide a labelled matrix of pairwise genetic distances, as there are many available packages for calculating these. This input matrix is provided using the --genmat argument, and should be tab-delimited with both column and row labels matching your population or individual identifiers. 

Built-in genetic distance calculations are currently in beta, meaning they are provided conditionally and still require more extensive external testing. They appear to be by-and-large functioning as written, but I would like to have a much more comprehensive test of whether or not my implementation of each statistic produces identical (or very similar) results to those of other packages such as Arlequin or adegenet. If you use autoStreamTree for your research and have the ability to directly compare your genetic distance matrix with those of other packages, please do so and feel free to let me know how they compare. Please note that I also offer the ability to directly import a genetic distance matrix that has been externally calculated, skipping this step altogether. 

Currently I provide options for individual-based distances (p-distance and various substitution-model corrected distances), and population-based distances which use allele frequency data (such as Fst):

```
	Genetic distance options:
		-d,--dist	: Use which metric of distance? Options are:
			Substitution models (individual-based):
			  PDIST			: Uncorrected p-distances [# Differences / Length]
			  JC69 			: [default] Jukes-Cantor (1969) corrected p-distances
			  K2P			: Kimura 2-parameter distances
			  TN84			: Tajima and Nei's (1984) distance
			  TN93			: Tamura and Nei's (1993) distance
			Frequency models (when using --pop):
			  FST			: Weir and Cockerham's Fst formulation (theta)
			  GST			: Hedrick's (2005) correction of Nei's (1987) Gst [=G'st]
			  GSTPRIME		: Meirmans & Hedrick (2011) corrected G'st [=G''st]
			  LINFST		: [default] Rousset's (1997) linearized Fst [=Fst/(1-Fst)]
			  JOST			: Jost's (2008) D
			  LINJOST		: 1/1-D, where D=Jost's (2008) D
			  NEI72			: Nei's (1972) standard genetic distance 
			  NEI83			: Nei and Chesser (1983) Da
			  EUCLID		: Euclidean distance
			  CHORD			: Cavalli-Sforza and Edwards (1967) chord distance
```

Optionally, the user can also opt to aggregate individual-based distance measures (when using a priori population assignments or the --geopop option). This can be provided using the --pop_agg argument, with any of the following options available:

```
	Aggregation options: 
		--pop_agg	: Define aggregator function for certain genetic distances w/ --pops:
			All of these can take the following options:
			  ARITH		: [default] Use arithmetic mean
			  MEDIAN	: Use median distance
			  HARM		: Use harmonic mean
			  ADJHARM	: Adjusted harmonic mean (see docs)
			  GEOM		: Use geometric mean
			  MIN		: Use minimum distance
			  MAX		: Use maximum distance
```

For datasets containing multiple non-concatenated loci, note that individual-based distances (e.g. PDIST or JC69) will also need to be aggregated among loci within each pairwise calculation. Any of the above options can again be used here, provided using the --loc_agg argument. 

#### Defining populations <a name="pops"></a>
There are currently three ways in which you can define populations for population-wise analysis. The first (specified using --pop) assumes that the 2nd column in the input file contains population identifiers. These can take any form (e.g., integer or string). The second (--geopop) will group any samples into populations which "snap" to the same stream node (see below). 

A third option (--clusterpop) will automatically cluster geographically similar individuals using the DBSCAN algorithm in scikit-learn, using great-circle geographic distances (i.e., this is not informed by stream distances calculated as a part of some workflows). Two relevant options are provided for manipulating the DBSCAN results:
```
DBSCAN options (only when --clusterpop):
	--min_samples	: Minimum samples per cluster [default=1]
	--epsilon	: Maximum distance (in km) within a cluster [default=20]
```

If using population labels, whether provided in the input file (--pop/--geopop) or calculating using DBSCAN (--clusterpop), autoStreamTree will output a plot showing cluster membership in 2-D space called $OUT.clusteredPoints.pdf:

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/example.clusteredPoints.png)

In this example, DBSCAN was used (hence population IDs are formatted as "DB_"#). Population centroids, which are ultimately used to "snap" populations to the stream network are shown with an "x". Note that this means that the population will only be represented by a single point on the network! 


### Example workflows <a name="ast_workflow"></a>

#### Fitting single-marker distances <a name="ast_example1"></a>

#### Working with large SNP datasets <a name="ast_example2"></a>

#### Microsatellites <a name="ast_example3"></a>

### Runtimes and benchmarking <a name="ast_benchmark"></a>

### References <a name="ast_refs"></a>
#### Citations for FitDistNet methods 
Below is a full list of citations for the various methods used in FitDistNet. Apologies to anyone I missed - feel free to let me know if you notice any discrepancies. 
* Beyer WM, Stein M, Smith T, Ulam S. 1974. A molecular sequence metric and evolutionary trees. Mathematical Biosciences. 19: 9-25.
* Cavalli-Sforza LL, Edwards AWF. 1967. Phylogenetic analysis: model and estimation procedures. American Journal of Human Genetics. 19: 233-257.
* Ester M, Kriegel HP, Sander J, Xu X. 1996. A density-based algorithm for discovering  clusters in large spatial databases with noise. IN: Simoudis E, Han J, Fayyad UM. (eds.). Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231.
* Felsenstein J. 2004. Inferring Phylogenies: Chapter 11. Sunderland: Sinauer.
* Fitch WM, Margloiash E. 1967. Construction of phylogenetic trees. Science. 155: 279-84.
* Hagberg A, Swart P, S Chult D. 2008. Exploring network structure, dynamics, and function using NetworkX. Los Alamos National Lab.(LANL), Los Alamos, NM
* Hedrick PW. 2005. A standardized genetic differentiation measure. Evolution. 59: 1633–1638
* Jordahl K. 2014. GeoPandas: Python tools for geographic data. URL: https://github.com/geopandas/geopandas.
* Jost L. 2008. Gst and its relatives do not measure differentiation. Molecular Ecology. 17: 4015-4026.
* Jukes TH, Cantor CR. 1969. Evolution of protein molecules. New York: Academic Press.
* Kalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760)
* Kimura M. 1980. A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of Molecular Evolution. 16(2): 111-120.
* Mantel N. 1967. The detection of disease clustering and a generalized regression approach. Cancer Research 27(2): 209-220.
* Meirmans PG, Hedrick PW. 2011. Assessing population structure: Fst and related measures. Molecular Ecology Resources. 11: 5-18.
* Nei M. 1972. Genetic distance between populations. American Naturalist. 106: 283-292.
* Nei M. 1987. Molecular Evolutionary Genetics. Columbia University Press, New York
* Nei M, Chesser RK. 1983. Estimation of fixation indices and gene diversities. Annals of Human Genetics 47(3): 253-259.
* Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, Blondel M, Prettenhofer P, Weiss R, Dubourg V, Vanderplas J. 2011. Scikit-learn: Machine learning in Python. The Journal of machine Learning research. 1(12):2825-30
* Rossmann LA. DFLOW User's Manual. U.S. Environmental Protection Agency.[For description of zero-adjusted harmonic mean]
* Rousset F. 1997. Genetic differentiation and estimation of gene flow from F-statistics under isolation by distance. Genetics. 145: 1219-28.
* Tajima F, Nei M. 1984. Estimation of evolutionary distance between nucleotide sequences. Molecular Biology and Evolution 1:269-285
* Tamura K, Nei M. 1993. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Molecular Biology and Evolution. 10(3):512-526.
* Weir BS, Cockerham CC. 1984. Estimating F-statistics for the analysis of population structure. Evolution. 38: 1358-1370.

#### Other reading
Here are some recommended readings and resources:
* Comte L, Olden JD. 2018. Fish dispersal in flowing waters: A synthesis of movement- and genetic-based studies. Fish and Fisheries. 19(6): 1063-1077. 
* Comte L, Olden JD. 2018. Evidence for dispersal syndromes in freshwater fishes. Proceedings Royal Society: B. 285(1871):  
* Grill, G., Lehner, B., Thieme, M. et al. 2019. Mapping the world’s free-flowing rivers. Nature. 569:215–221.
* Linke S, Lehner B, Ouellet Dallaire C. et al. 2019. Global hydro-environmental sub-basin and river reach characteristics at high spatial resolution. Sci Data 6, 283
* Meffe GK, Vrijenhoek RC. 1988. Conservation genetics in the management of desert fishes. Conservation Biology. 2(2):157-69.
* Meirmans PG. 2012. The trouble with isolation by distance. Molecular Ecology 21(12): 2839-46.
* Sere M, Thevenon S, Belem AMG, De Meeus T. 2017. Comparison of different genetic distances to test isolation by distance between populations. 2017. 119(2):55-63.
* Thomaz AT, Christie MR, Knowles LL. 2016. The architecture of river networks can drive the evolutionary dynamics of aquatic populations. Evolution. 70(3): 731-739.
* Tonkin JD, Altermatt F, Finn DS, Heino J, Olden JD, Pauls SU, Lytle DA. 2017. The role of dispersal in river network metacommunities: Patterns, processes, and pathways. Freshwater Biology. 61(1): 141-163.
* Wright S. 1965. Isolation by distance. Genetics. 28: 114-138.

## ResistNet <a name="rscape"></a>
### Program description <a name="rscape_desc"></a>


#### Genetic algorithm background and similarities to ResistNet

All variables can vary in their model inclusion, but also in the manner by which they are transformed; transformations are the same as those used in ResistNet, and shape parameters are included for each variable in the "chromosome". Note that as shape parameters increase, the more similar the "transformed" values are to the original values, with high values for the shape parameter (e.g., 50-100) resulting in negligable transformation: 

![](https://raw.githubusercontent.com/tkchafin/DENdriscape/master/examples/plots/transforms.png)

### Usage <a name="rscape_usage"></a>

#### Options and Help Menu <a name="rscape_help"></a>

As with all of the scripts in this repository, you can view a list of the options for ResistNet by calling the help menu with <-h>:
```
$ python3 ResistNet.py -h
Exiting because help menu was called.

    o--o                                                    o-o     O  
    |   |  o                                               o       / \ 
    O-Oo      o   o  o-o  o-o  o-o   o-o   oo   o-o   o-o  |  -o  o---o
    |  \   |   \ /   |-'  |     \   |     | |   |  |  |-'  o   |  |   |
    o   o  |    o    o-o  o    o-o   o-o  o-o-  O-o   o-o   o-o   o   o
                                                |                   
                                                o                   

    Author: Tyler K Chafin, University of Arkansas
    Contact: tkchafin@uark.edu
    Description: Genetic algorithm to optimize resistance models on networks

	Input options:
	  If using autoStreamTree outputs:
		-p,--prefix	: Prefix for FitDistNet outputs
		
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

For convienience, the inputs for ResistNet follow the formats of the files output by FitDistNet and FormatNet (not complete yet). 

#### Outputs <a name="rscape_output"></a>

After parsing all of the inputs, ResistNet will randomly generate a population of 'individuals' (=model parameterizations), which is by default 15X the number of parameter, up to a maximum size specified by <-P,--maxPop>. Alternatively, you can specify a fixed size with <-s,--size>. Each 'generation', individuals will be selected using the desired fitness function (specified with <-f,--fit>; e.g. -f AIC to use AIC), and the maximum, minimum, meand, and standard deviation of population fitness values will be output to the terminal (stdout):
```
Reading network from:  ../out3.network
Reading FitDistNet results from: ../out3.streamTree.txt
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
Description: Script for re-plotting StreamTree results after running FitDistNet

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
