# autoStreamTree
Automated pipeline for examining patterns of genetic differentiation in stream networks from single-gene, multi-gene, or SNP datasets.

--In progress--

Contact: tkchafin@uark.edu 

## Software Description
autoStreamTree is a Python software package providing various analyses aimed at analyzing patterns of genetic differentiation among aquatic stream-dwelling organisms. The intention is to take what was previously a tedious process involving multiple discrete steps and to integrate these all in one place. 

Currently, autoStreamTree provides a companion library of functions for calculating various measures of genetic distances among individuals or populations, including model-corrected p-distances (e.g. Jukes-Cantor 1969, Kimura 2-parameter, Tamura-Nei 1993) as well as those based on allele frequencies (e.g. Theta-ST, linearized Fst, Jost's D -- full list of available distance models below). It also includes integrated functions for parsing an input vector shapefile of streams (see below 'Requirements for input shapefiles') for easy calculation of pairwise stream distances between sites, as well as the ordinary or weighted least-squares fitting of reach-wise genetic distances according to the "stream tree" model of Kalinowski et al. (2008). Various plotting functions are also provided for downstream analysis, including looking at patterns of isolation-by-distance. Outputs should also be directly importable into R, with additional outputs with annotated streamtree fitted distances provided for analysis in your GIS suite of choice. 

If you use this package for analysis of fitted distances using the streamtree model, please cite the following:
* Kalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760).

If you find the code here useful for your research, for now please cite this repository:
* Chafin TK. autoStreamTree: Automating workflows for examining patterns of genetic differentiation in stream networks. DOI: Coming soon.

## Installation

### Prerequisites
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

### Installation with conda

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


## Usage

### Getting started

Fill in details later... 

Points will be 'snapped' to nodes in the stream network. autoStreamTree will output both a table($out.snapDistances.txt) and a histogram plot ($out.snapDistances.pdf) showing distances in kilometers that samples or populations had to be snapped:

![](https://raw.githubusercontent.com/tkchafin/autoStreamTree/master/examples/plots/example.snapDistances.png)

### Genetic distance models 

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

### Defining populations
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

### StreamTree method

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

### Requirements for input shapefiles

I highly recommend using the existing global stream datasets provided by the [HydroLab group](https://wp.geog.mcgill.ca/hydrolab/) at McGill University, specifically the [HydroAtlas](https://www.hydrosheds.org/page/hydroatlas) or [free-flowing rivers dataset](https://wp.geog.mcgill.ca/hydrolab/free-flowing-rivers/) as these are already properly formatted for use, and the additional variables included will make downstream analysis very easy. Because of their size, I would recommend clipping them to the relevant scale first (e.g. the drainage encompassing all of your samples).

Note that a valid path is required between all sites in order to calculate pairwise stream distances. Thus, if you are analyzing data from multiple drainages which only share an oceanic connection, you will need to augment the shapefile. For example this could be accomplished by adding a vector representing the coastline to create an artificial connection among drainages. 

If for some reason you cannot use the HydroRIVERS dataset, you will need to do some things first before loading your shapefile into autoStreamTree. First, you will need to include two variables in the attribute table of your shapefile: 1) REACH_ID (case sensitive) must provide a unique identifier to each stream reach; and 2) LENGTH_KM should give the length of each segment. Next, because sometime large stream layers will have small gaps in between streams, you will need to span any small gaps between streams which should be contiguous, and also dissolve any lines that overlap with one another so that any given section of river is represented by a single line. I will provide a tutorial for doing this in ArcMAP later, but for now there are some scripts in our complementary package that can help with these steps using the ArcPy API: https://github.com/stevemussmann/StreamTree_arcpy. Note that this package will also help you in running the original Stream Tree package on Windows, if you want to do so. 

## Example workflows 

## Accessory scripts

### Re-plotting StreamTree outputs

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

### Clustering populations using any distance matrix

For example, if you wanted to cluster individuals using their stream distances, I've provided a script called clusterPopsDB.py which will use a DBSCAN clustering algorithm to output a tab-delimited population map given any arbitrary distance matrix:
'''
$ python3 ./scripts/clusterPopsDB.py -h
'''


