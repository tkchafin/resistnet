import os
import sys

"""
Here, build an output reference string for all methods a user chose in a run.
e.g. if using StreamTree, output Kalinowski et al
For weighted LS, gen distances, etc
"""

	References:
		Beyer WM, Stein M, Smith T, Ulam S. 1974. A molecular sequence metric and
		evolutionary trees. Mathematical Biosciences. 19: 9-25.
		Cavalli-Sforza LL, Edwards AWF. 1967. Phylogenetic analysis: model and estimation
			procedures. American Journal of Human Genetics. 19: 233-257.
		Felsenstein J. 2004. Inferring Phylogenies: Chapter 11. Sunderland: Sinauer.
		Fitch WM, Margloiash E. 1967. Construction of phylogenetic trees. Science.
			155: 279-84.
		Hedrick PW. 2005. A standardized genetic differentiation measure.
			Evolution. 59: 1633–1638
		Jost L. 2008. Gst and its relatives do not measure differentiation. Molecular
			Ecology. 17: 4015-4026.
		Jukes TH, Cantor CR. 1969. Evolution of protein molecules. New York: Academic Press.
		Kimura M. 1980. A simple method for estimating evolutionary rates of base
			substitutions through comparative studies of nucleotide sequences. 
			Journal ofMolecular Evolution. 16(2): 111-120.
		Meirmans PG, Hedrick PW. 2011. Assessing population structure: Fst and related
			measures. Molecular Ecology Resources. 11: 5-18.
		Nei M (1987) Molecular Evolutionary Genetics. Columbia University Press,
			New York
		Nei M, Chesser RK. 1983. Estimation of fixation indices and gene diversities.
			Annals of Human Genetics 47(3): 253-259.
		Rossmann LA. DFLOW User's Manual. U.S. Environmental Protection Agency.
			[For description of zero-adjusted harmonic mean]
		Rousset F. 1997. Genetic differentiation and estimation of gene flow from
			F-statistics under isolation by distance. Genetics. 145: 1219-28.
		Weir BS, Cockerham CC. 1984. Estimating F-statistics for the analysis of population
			structure. Evolution. 38: 1358-1370.

	Recommended reading:
		Meirmans PG. 2012. The trouble with isolation by distance. Molecular Ecology
			21(12): 2839-46.
		Sere M, Thevenon S, Belem AMG, De Meeus T. 2017. Comparison of different genetic
			distances to test isolation by distance between populations. 2017. 119(2):55-63.
		Wright S. 1965. Isolation by distance. Genetics. 28: 114-138.