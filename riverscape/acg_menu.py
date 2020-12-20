import sys
import os
import getopt

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hp:g:n:s:T:P:G:s:m:i:c:t:F:d:D:f:b:C:v:Aa:o:X', \
			["shp=", "help", "input=", "prefix=", "genmat=", "network=",
			"seed=", "procs=", "maxPop=", "maxpop=", "maxgen=", "maxGen=",
			"size=", "popsize=", "mutpb=", "indpb=", "cxpb=", "tourn=", 
			"nfail=", "nFail=", "delt=", "deltP=", "deltp=", "fit=", "metric=", "fitness=",
			"burn=", "force=", "infer", "cholmod", "cprocs=", "Cprocs=", "vars=", "modavg",
			"modAvg", "awsum=", "report_all", "noPlot", "out=", "keep_all"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.prefix="out3"
		self.force=None
		self.variables = None
		self.seed=None
		self.installCS=False
		self.popsize=None
		self.maxpopsize=100
		self.cstype="pairwise"
		self.fitmetric="aic"
		self.predicted=False
		self.inmat=None
		self.network=None
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
		self.report_all=False
		self.plot=True
		
		self.only_keep=True


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg in options:
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == 'p' or opt=='prefix':
				self.prefix = arg
			elif opt=='g' or opt=='genmat':
				self.inmat = arg
			elif opt=='n' or opt=='network':
				self.network = arg
			elif opt=='s' or opt=='seed':
				self.seed=int(arg)
			elif opt=='T' or opt=='procs':
				self.GA_procs=int(arg)
			elif opt=='P' or opt=='maxPop' or opt=='maxpop':
				self.maxpopsize=int(arg)
			elif opt=="G" or opt=="maxGen" or opt=="maxgen":
				self.maxGens=int(arg)
			elif opt=="s" or opt=="popsize" or opt=="size":
				self.popsize=int(arg)
			elif opt=="m" or opt=="mutpb":
				self.mutpb=float(arg)
			elif opt=="i" or opt=="indpb":
				self.indpb=float(arg)
			elif opt=="c" or opt=="cxpb":
				self.cxpb=float(arg)
			elif opt=="t" or opt=="tournsize":
				self.tournsize=int(arg)
			elif opt=="F" or opt=="nfail" or opt=="nFail":
				self.nfail=int(arg)
			elif opt=="d" or opt=="delt":
				self.deltaB=float(arg)
			elif opt=="D" or opt=="deltP" or opt=="deltp":
				self.deltaB_perc=float(arg)
			elif opt=="f" or opt=="fit" or opt=="metric" or opt=="fitness":
				if arg.lower() not in ["aic", "r2m", "loglik", "delta"]:
					self.diplay_help("Unrecognized fitness metric <-f, --fit>")
				else:
					self.fitmetric=arg.lower()
			elif opt=="b" or opt=="burn":
				self.burnin=int(arg)
			elif opt=="force":
				self.force=arg
			elif opt=="infer":
				self.predicted=True
			elif opt=="cholmod":
				self.cholmod=True
			elif opt=="C" or opt=="Cprocs" or opt=="cprocs":
				self.CS_procs=int(arg)
			elif opt=="v" or opt=="vars":
				self.variables=arg.split(",")
			elif opt=="A" or opt=="modavg" or opt=="modAvg":
				self.modavg=True
			elif opt=="a" or opt=="awsum":
				self.awsum=float(arg)
			elif opt=="report_all":
				self.report_all=True
			elif opt=="X" or opt=="noPlot":
				self.plot=False
			elif opt=="o" or opt=="out":
				self.out=arg
			elif opt=="avgall":
				self.only_keep=False
			elif opt == 'h' or opt == 'help':
				pass
			else:
				assert False, "Unhandled option %r"%opt
		
		if self.variables is None:
			self.display_help("No variables selected.")

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
	
	Circuitscape options:
		--cholmod	: Turn on CHOLMOD solver
		-C,--cprocs	: Processors per Circuitscape process [default=1]
				NOTE: Total simultaneous processes= <-T> * <-C>
	
	Genetic distance options:
		--force	: Use XX attribute from input table as distance metric (e.g. 'fittedD')
				NOTE: By default, the average of "locD_" columns will be taken
		--infer		: Infer pairwise distances from input table (i.e., NOT input matrix)
		-v,--vars: Comma-separated list (no spaces) of explanatory attributes to include
	
	Multi-model inference options:
		-A,--modavg	: Compute model-averaged resistances
			NOTE: This involves re-running Circuitscape for each model
		-a,--awsum	: Cumulative Akaike weight threshold to retain top N models [default=0.95]
		--report_all: Plot per-stream resistance and generate full outputs for all retained models

""")
		print()
		sys.exit()