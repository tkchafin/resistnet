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
		self.seed=None
		self.installCS=False
		self.popsize=None
		self.maxpopsize=100
		self.cstype="pairwise"
		self.fitmetric="aic"
		self.predicted=False
		self.inmat=None
		self.cholmod=False
		self.GA_procs=1
		self.CS_procs=1
		self.deltaB=None
		self.deltaB_perc=0.001
		self.nfail=50
		self.maxGens=500
		self.tournsize=10
		self.cxpb=0.5
		self.mutpb=0.2
		self.indpb=0.1
		self.burnin=0
		self.max_hof_size=100
		self.out="output"
		self.awsum=0.95
		self.modavg=False


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
		print ("\nRiverscapeGA.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Genetic algorithm to optimize resistance models on networks")
		print("""
	Input options:
	  If using autoStreamTree outputs:
		-p,--prefix	: Prefix for autoStreamTree outputs
		-or-
	  If manually specifying inputs:
		-g,--genmat	: Genetic distance matrix
		-n,--network	: Input graph (in pickle'd networkx format)	
		
	General options:
		-s,--seed	: Random number seed (default=taken from clock time)
		-P,--procs	: Number of parallel processors
		-x,--noPlots: Turn off plotting
	
	Genetic Algorithm Options:
		-P,--maxPop	: Maximim population size [default = 100]
		-M,--maxGen	: Maximum number of generations [default = 500]
		-s,--size	: Manually set population size to <-p int>
				    NOTE: By default, #params * 15
		-m,--mutpb	: Probability of mutation per individual [default=0.2]
		-i,--indpb	: Probability of mutation per trait [default=0.1]
		-c,--cxpb	: Probability of being chosen for cross-over [default=0.5]
		-t,--tourn	: Tournament size [default=10]
		
	Model optimization/ selection options:
		-N,--nfail	: Number of generations failing to improve to stop optimization
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
	
	Circuitscape options:
		-C,--cprocs	: Number of processors to use *per* Circuitscape run
		-I,--include: Comma-separated (NO SPACE) list of variables to consider
		-X,--exclude: Comma-separated (NO SPACE) list of variables to exclude
		--cholmod	: Turn on CHOLMOD solver (see Circuitscape docs)
	
	Genetic distance options:
		-f,--force	: Use XX attribute from input table as distance metric
		--infer		: Infer pairwise distances from input table (i.e., NOT pairwise matrix)
		
		-I,--include: Comma-separated list (no spaces) of explanatory attributes to include
		-X,--exclude: Comma-separated list (no spaces) of explanatory attributes to exclude
		-h,--help	: Displays help menu
	
	Multi-model inference options:
		-a,--modavg	: Compute model-averaged resistance per stream segment
			NOTE: This involves re-running Circuitscape for each model
		--avgpw	: Output model-averaged pairwise resistance distances
		-A,--awSum	: Cumulative Akaike weight threshold to retain top N models [default=0.95]
		--avgall	:Plot per-stream resistance and generate full outputs for all retained models
		--incall	: Include all models in calculation of relative variable importance
		--rvi_thresh	: Threshold 'importance' to draw on variable importance plot [default=0.8]
		--aic_thresh	: AIC difference from best model to draw in IC profile plot [default=2]

""")
		print()
		sys.exit()