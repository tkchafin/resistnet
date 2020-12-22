# Riverscape Genetics
A collection of software tools for examining patterns of genetic differentiation in riverscapes 

### Table of Contents: 

### Installation

Because of the number of dependencies, I recommend setting up a virtual environment to prevent any conflicts with your system Python environment. 

If you are planning on using RiverscapeGA, [PyJulia](https://pyjulia.readthedocs.io/en/latest/), which forms the necessary interface for Python to access Circuitscape (a Julia program), a complication is that PyJulia cannot use a Python distribution that is statically linked to libpython (such as that installed by Ubuntu or conda, or Debian-based Linux distributions). 

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

#### Installation using conda

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

If your Python is statically linked (see above ldd command), you will also need to create a custom Julia system image that you will pass to RiverscapeGA:

```
python3 -m julia.sysimage sys.so
```


#### Installation using pyenv and manually compiled Julia/R

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

#check that the right python is in your path:
which python3 
which pip3
#output should be something like /home/USER/.pyenv/*/python3; if not, either add the correct path to .bashrc or use absolute path in the python3 and pip3 calls below

#upgrade pip and make sure venv is installed 
pip3 install --upgrade pip3
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

If you do not have the [R](https://www.r-project.org/) and [Julia](https://julialang.org/) programming languages installed on you systems, you can install them directly into your Python virtual environment like so. Note that you can also use pre-compiled binaries, but for the sake of completion I will here show how to compile them from source.

```
cd ~/python_venv/riverscape
#go to https://julialang.org/downloads/ if you want a pre-compiled binary
#NOTE: On a Redhat Linux HPC, I could only get the LTS v1.0.5 to work
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

#### Troubleshooting PyJulia 

As noted above, PyJulia won't work with statically-linked Python. There are several workarounds in the PyJulia [documentation](https://pyjulia.readthedocs.io/en/latest/troubleshooting.html). 

One of the easiest workarounds is to turn off the Julia compiled cache, which you can do by passing Riverscape GA the '--no_compiled_modules' (boolean) argument. This might slow down loading Julia a little bit, but is one of the fastest ways to get up and running if you are getting PyJulia errors. 

Another option is to create a custom JUlia system image, which can be passed to RiverscapeGA using the --sys_image argument:

```
python3 -m julia.sysimage sys.so
```

#### Using the pre-build Singularity image
Coming soon...

### Programs: 
- *autoStreamTree.py* : This package calculates genetic and stream-distances by snapping sampling points to a network provided using an input shapefile and calculating a minimal sub-network, tests for isolation-by-distance, and fits genetic distances to stream segments using the Stream-tree algorithm by Kalinowski et al. 2008. 
- *RiverscapeGA.py* : This package implements a genetic algorithm to optimize multi-variable resistance models on networks in Circuitscape. This is somewhat similar to the way in which ResistanceGA accomplishes this for raster datasets, but with some critical algorithmic differences (see below)
- *streamCleaner.py* : (Coming soon) Integrates the network extraction functions of autoStreamTree with various pre-processing steps such as joining non-contiguous segments, merging redundant paths (e.g., for a braided stream) -- Not yet added
- *riverscapeSim.py* : (Coming soon) Performs forward-simulations on a given stream network, given environmental 'resistance' variables (and optional weights) -- Not yet added
- *BGR_pipeline.py* : (Coming soon) Pipeline for applying BGR model given a set or covariates or resistance model
- tools/ 
	* *autoFetcher.py* : Python interface for NCBI entrez API; can be used to automatically download and parse queries
	* *clusterPopsDB.py* : Distance-based clustering of points using the DBSCAN algorithm
	* *fasta2phylip.py* : Converting between FASTA and PHYLIP formats
	* *fasta2table.py* : Converting between FASTA and autoStreamTree Table formats
	* *nremover.pl* : Filtering FASTA and PHYLIP files of concatenated SNPs
	* *plotStreamTree.py* : Re-making StreamTree fitted distance plot with some additional options
	* *modelParser.py* : Parse model outputs of RiverscapeGA.py to generate and select transformed datasets (Coming soon)
	* [tkchafin/scripts](https://github.com/tkchafin/scripts) : For more useful file formatting and filtering scripts
	* [stevemussmann/Streamtree_arcpy](https://github.com/stevemussmann/StreamTree_arcpy) : Code for generating inputs for the original StreamTree program by Kalinowski et al.

--Work in progress--

Contact: tkchafin@uark.edu 

### Quick start

## autoStreamTree

### Software Description
autoStreamTree is a Python software package providing various analyses aimed at analyzing patterns of genetic differentiation among aquatic stream-dwelling organisms. The intention is to take what was previously a tedious process involving multiple discrete steps and to integrate these all in one place. 

Currently, autoStreamTree provides a companion library of functions for calculating various measures of genetic distances among individuals or populations, including model-corrected p-distances (e.g. Jukes-Cantor 1969, Kimura 2-parameter, Tamura-Nei 1993) as well as those based on allele frequencies (e.g. Theta-ST, linearized Fst, Jost's D -- full list of available distance models below). It also includes integrated functions for parsing an input vector shapefile of streams (see below 'Requirements for input shapefiles') for easy calculation of pairwise stream distances between sites, as well as the ordinary or weighted least-squares fitting of reach-wise genetic distances according to the "stream tree" model of Kalinowski et al. (2008). Various plotting functions are also provided for downstream analysis, including looking at patterns of isolation-by-distance. Outputs should also be directly importable into R, with additional outputs with annotated streamtree fitted distances provided for analysis in your GIS suite of choice. 

If you use this package for analysis of fitted distances using the streamtree model, please cite the following:
* Kalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760).

If you find the code here useful for your research, for now please cite this repository:
* Chafin TK. autoStreamTree: Automating workflows for examining patterns of genetic differentiation in stream networks. DOI: Coming soon.

### Installation

#### Prerequisites
* Python 3
* numpy
* scipy
* Matplotlib
* Seaborn
* networkx
* pandas
* geopandas
* sortedcontainers
* shapely
* geopy
* scikit-bio
* scikit-learn

#### Installation with conda

A full conda installation will come soon. For now, you can install manually using conda like so:
```
#update conda and shared packages and set channel
conda config --add channels conda-forge
conda update --all

#create a virtual env 
conda create -n geo python=3.6

#activate virtual env
conda activate geo

#install gdal
conda install gdal

#deactivate + reactivate; sets the environmental variables
conda deactivate; conda activate geo

#install the rest of the packages
conda install scipy pandas geopandas shapely geopy sortedcontainers matplotlib networkx seaborn scikit-learn

```

Next, clone this directory and you're all set!


### Usage

#### Getting started

Fill in details later... 

Points will be 'snapped' to nodes in the stream network. autoStreamTree will output both a table($out.snapDistances.txt) and a histogram plot ($out.snapDistances.pdf) showing distances in kilometers that samples or populations had to be snapped:

![](https://raw.githubusercontent.com/tkchafin/autoStreamTree/master/examples/plots/example.snapDistances.png)

#### Genetic distance models 

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

#### Defining populations
There are currently three ways in which you can define populations for population-wise analysis. The first (specified using --pop) assumes that the 2nd column in the input file contains population identifiers. These can take any form (e.g., integer or string). The second (--geopop) will group any samples into populations which "snap" to the same stream node (see below). 

A third option (--clusterpop) will automatically cluster geographically similar individuals using the DBSCAN algorithm in scikit-learn, using great-circle geographic distances (i.e., this is not informed by stream distances calculated as a part of some workflows). Two relevant options are provided for manipulating the DBSCAN results:
```
DBSCAN options (only when --clusterpop):
	--min_samples	: Minimum samples per cluster [default=1]
	--epsilon	: Maximum distance (in km) within a cluster [default=20]
```

If using population labels, whether provided in the input file (--pop/--geopop) or calculating using DBSCAN (--clusterpop), autoStreamTree will output a plot showing cluster membership in 2-D space called $OUT.clusteredPoints.pdf:

![](https://raw.githubusercontent.com/tkchafin/autoStreamTree/master/examples/plots/example.clusteredPoints.png)

In this example, DBSCAN was used (hence population IDs are formatted as "DB_"#). Population centroids, which are ultimately used to "snap" populations to the stream network are shown with an "x". Note that this means that the population will only be represented by a single point on the network! 

#### StreamTree method

Coming soon -- some changes

The first thing autoStreamTree will do upon reading your input shapefile is to calculate a minimally reduced sub-network which collapses the input river network into continuous reaches (="edges"), with nodes either representing sample localities or junctions. Because the full river network will likely contain many branches and contiguous reaches which do not contain samples, these are removed to speed up computation. The underlying metadata will be preserved, and the final output will consist of an annotated shapefile containing an EDGE_ID attribute which tells you how reaches were dissolved into contiguous edges in the graph, and a FittedD attribute giving the least-squares optimized distances.

The reduced sub-network will be plotted for you in a file called $OUT.subGraph.pdf:
![](https://raw.githubusercontent.com/tkchafin/autoStreamTree/master/examples/plots/example.subGraph.png)

Here, the total cumulative stream length (in km) is plotted along edges (NOTE: Any natural curvature in the river is not preserved in this plot), with sample sites as blue dots and junctions as black dots. A geographically accurate representation, coloring individual streams to designate different dissolved edges, will be provided as $out.streamsByEdgeID.pdf: 
![](https://raw.githubusercontent.com/tkchafin/autoStreamTree/master/examples/plots/example.networkByEdgeID.png)

After fitting genetic distances, autoStreamTree will create several other outputs. First, a table called $out.reachToEdgeTable.txt will give a tab-delimited map of how REACH_ID attributes were dissolved into contiguous edges. Second, a tabular and graphical representation of how fitted pairwise distances compare to the raw calculates (or user-provided) pairwise distances: $out.obsVersusFittedD.txt and $out.obsVersusFittedD.pdf
![](https://raw.githubusercontent.com/tkchafin/autoStreamTree/master/examples/plots/example.obsByFittedD.png)

Finally, the fitted distances per stream edge will be output both as an added column to the original shapefile attribute table ($out.streamTree.shp and $out.streamTree.txt), and also as a plot showing how distances compare across all streams: $out.streamsByFittedD.pdf
![](https://raw.githubusercontent.com/tkchafin/autoStreamTree/master/examples/plots/example.networkByStreamTree.png)

#### Requirements for input shapefiles

I highly recommend using the existing global stream datasets provided by the [HydroLab group](https://wp.geog.mcgill.ca/hydrolab/) at McGill University, specifically the [HydroAtlas](https://www.hydrosheds.org/page/hydroatlas) or [free-flowing rivers dataset](https://wp.geog.mcgill.ca/hydrolab/free-flowing-rivers/) as these are already properly formatted for use, and the additional variables included will make downstream analysis very easy. Because of their size, I would recommend clipping them to the relevant scale first (e.g. the drainage encompassing all of your samples).

Note that a valid path is required between all sites in order to calculate pairwise stream distances. Thus, if you are analyzing data from multiple drainages which only share an oceanic connection, you will need to augment the shapefile. For example this could be accomplished by adding a vector representing the coastline to create an artificial connection among drainages. 

If for some reason you cannot use the HydroRIVERS dataset, you will need to do some things first before loading your shapefile into autoStreamTree. First, you will need to include two variables in the attribute table of your shapefile: 1) REACH_ID (case sensitive) must provide a unique identifier to each stream reach; and 2) LENGTH_KM should give the length of each segment. Next, because sometime large stream layers will have small gaps in between streams, you will need to span any small gaps between streams which should be contiguous, and also dissolve any lines that overlap with one another so that any given section of river is represented by a single line. I will provide a tutorial for doing this in ArcMAP later, but for now there are some scripts in our complementary package that can help with these steps using the ArcPy API: https://github.com/stevemussmann/StreamTree_arcpy. Note that this package will also help you in running the original Stream Tree package on Windows, if you want to do so. 

### Example workflows 

#### Re-plotting StreamTree outputs

The default $out.streamdByFittedD plot may not be exactly what you wanted. To prevent cluttering the help menu of the main program too much, we've provided a separate script for loading up autoStreamTree outputs to re-make the plot, which has some added options for customization: scripts/plotStreamTree.py

'''
$ python3 ./scripts/plotStreamTree -h
plotStreamTree.py

Author: Tyler K Chafin, University of Arkansas
Contact: tkchafin@uark.edu
Description: Script for re-plotting StreamTree results after running autoStreamTree

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

### References
#### Citations for autoStreamTree methods
Below is a full list of citations for the various methods used in autoStreamTree. Apologies to anyone I missed - feel free to let me know if you notice any discrepancies. 
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

#### Recommended reading:
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

## riverscapeGA
### Program description

### Installation 

#### Dependencies

#### Conda installation
Conda instructions coming soon... 


#### PyJulia setup

```
python3 -m pip install julia
```

You will need to open Julia and install the Julia [Circuitscape](https://github.com/Circuitscape/Circuitscape.jl) package and set it up to work with your python:
```
$ julia
julia> using Pkg; Pkg.add("Circuitscape")
julia> Pkg.test("Circuitscape")
julia> ENV["PYTHON"] = "python3"  # whatever path you have
julia> Pkg.build("PyCall")
```


Depending on your environment, you may encounter an error along the lines of ```Your Python interpreter "/Users/tyler/miniconda3/envs/geo/bin/python3" is statically linked to libpython.  Currently, PyJulia does not fully support such Python interpreter.```
If this happens, you will need to use a workaround to get Python and Julia properly working together. 

```
#install pyenv
brew install pyenv
#build python from scratch
PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.6.6
ldd ~/.pyenv/versions/3.6.6/bin/python3.6 | grep libpython
#re-install dependencies w/ new python environment
/Users/tyler/.pyenv/versions/3.6.6/bin/pip3 install numpy scipy networkx julia functools deap pandas
#open python and set up julia PyCall
/Users/tyler/.pyenv/versions/3.6.6/bin/python3 
>>> import julia
>>> julia.install()
>>> quit()
#install rpy2 on MacOS requires you point to gcc:
env CC=/usr/local/bin/gcc /Users/tyler/.pyenv/versions/3.6.6/bin/pip3 install rpy2
```


 coming soon

## Example Analysis

