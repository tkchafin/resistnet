import sys
import os
import getopt

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hs:i:r:pd:a:lw:o:gP:L:Scn:G:', \
			["shp=", "help", "input="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.prefix="out3"
		self.force="fittedD"
		self.variables = ["tmp_dc_cmn", "aet_mm_cyr", "USE"]
		self.seed="1321"
		self.installCS=False
		self.popsize=None
		self.maxpopsize=20
		self.cstype="pairwise"
		self.fitmetric="aic"
		self.fitmetric_index=2
		self.predicted=False
		self.inmat=None
		self.cholmod=False
		self.GA_procs=3
		self.CS_procs=1
		self.deltaB=None
		self.deltaB_perc=0.01
		self.nfail=10
		self.maxGens=5
		self.tournsize=5
		self.cxpb=0.5
		self.mutpb=0.5
		self.indpb=0.1
		self.burnin=0
		self.max_hof_size=100
		self.out="output"


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg in options:
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == 's' or opt == 'shp':
				self.shapefile = arg
			elif opt == 'h' or opt == 'help':
				pass
			else:
				assert False, "Unhandled option %r"%opt


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nautoCircuitGA.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Genetic algorithm to optimize resistance models on networks")
		print("""
	Input options:
	  If using autoStreamTree outputs:
		-p,--prefix	: Prefix for autoStreamTree outputs
		-or-
	  If manually specifying inputs:
		-i,--inmat	: Genetic distance matrix
		-n,--network	: Input graph (in pickle'd networkx format)	
		
	General options:
		-s,--seed	: Random number seed (default=taken from clock time)
		
		
		-b,--boot	: Number of bootstraps <NOT IMPLEMENTED>
		-f,--fit	: Fit metric used to evaluate models <NOT IMPLEMENTED>
				    Options:
				    aic (default)
				    loglik (log-likelihood)
				    r2m (marginal R^2)
				    deltaAIC (Change in AIC versus null model)
		-p,--pair	: Compute fitness metrics using pairwise distances [default] 
		-e,--edge	: Compute fitness metrics using edgewise distances <NOT IMPLEMENTED>
		-f,--force	: Force calculation using XX distance attribute in input table (e.g., "fittedD")
		-r,--inc	: Re-create incidence matrix from .network file [default: read $out.incidenceMatrix.txt] <NOT IMPLEMENTED>
		-g,--gen	: Read genetic distances from input matrix (Disables -e/--edge option) <NOT IMPLEMENTED>
		-I,--include: Comma-separated list (no spaces) of explanatory attributes to include
		-X,--exclude: Comma-separated list (no spaces) of explanatory attributes to exclude
		-h,--help	: Displays help menu
	
	Parallel execution parameters:
		-C,--Cprocs	: Number of processors for circuitscape <NOT IMPLEMENTED>
		
	Aggregation options: 
		-L,--loc_agg	: Define aggregator function for aggregating locus-wise distances
			All of these can take the following options:
			  ARITH		: [default] Use arithmetic mean
			  MEDIAN	: Use median distance
			  HARM		: Use harmonic mean
			  ADJHARM	: Adjusted harmonic mean (see docs)
			  GEOM		: Use geometric mean
			  MIN		: Use minimum distance
			  MAX		: Use maximum distance

""")
		print()
		sys.exit()