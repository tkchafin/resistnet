# autoStreamTree
Automated pipeline for examining patterns of genetic differentiation in stream networks

--In progress--

## Software Description
autoStreamTree is a Python software package providing various analyses aimed at analyzing patterns of genetic differentiation among aquatic stream-dwelling organisms. The intention is to take what was previously a tedious process involving multiple discrete steps and to integrate these all in one place. 

Currently, autoStreamTree provides a companion library of functions for calculating various measures of genetic distances among individuals or populations, including model-corrected p-distances (e.g. Jukes-Cantor 1969, Kimura 2-parameter, Tamura-Nei 1993) as well as those based on allele frequencies (e.g. Theta-ST, linearized Fst, Jost's D -- full list of available distance models below). It also includes integrated functions for parsing an input vector shapefile of streams (see below 'Requirements for input shapefiles') for easy calculation of pairwise stream distances between sites, as well as the least-squares fitting of reach-wise genetic distances according to the "stream tree" model of Kalinowski et al. (2008). Various plotting functions are also provided for downstream analysis, including looking at patterns of isolation-by-distance. Outputs should also be directly importable into R, with additional outputs with annotated streamtree fitted distances provided for analysis in your GIS suite of choice. 

If you use this package for analysis of fitted distances using the streamtree model, please cite the following:
* Kalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760).

## Installation

Conda installation coming soon. 

## Usage



### Available distance models 

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


### StreamTree method

Coming soon -- some changes, including added weights (see Felsenstein 2004 chapter II)

### Requirements for input shapefiles

Coming soon

### Example workflows 

Coming soon
