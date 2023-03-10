import sys
import os
import getopt

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hp:g:s:s:T:P:G:s:m:i:c:t:F:d:D:f:b:C:v:V:Aa:o:Xj:n:', \
			["shp=", "help", "input=", "prefix=", "genmat=", "shapefile=",
			"seed=", "procs=", "maxPop=", "maxpop=", "maxgen=", "maxGen=",
			"size=", "popsize=", "mutpb=", "indpb=", "cxpb=", "tourn=",
			"nfail=", "nFail=", "delt=", "deltP=", "deltp=", "fit=", "metric=", "fitness=",
			"burn=", "dist_col=", "vars=", "pop_agg=", "varfile=", "maxShape=",
			"awsum=", "report_all", "noPlot", "out=", "keep_all", "minWeight=",
			"max_hof_size=", "posWeight", "fixWeight", "allShapes", "efit_agg=",
			"coords=", "length_col=", "reachid_col=", "minimize", "network=", "fixShape"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.prefix="out3"
		self.dist_col=None
		self.variables = None
		self.agg_opts = dict()
		self.varFile=None
		self.minimize=False
		self.seed=None
		self.installCS=False
		self.popsize=None
		self.length_col="LENGTH_KM"
		self.reachid_col="HYRIV_ID"
		self.pop_agg="ARITH"
		self.edge_agg="ARITH"
		self.efit_agg="SUM"
		self.maxpopsize=100
		self.cstype="pairwise"
		self.fitmetric="aic"
		self.network=None
		self.predicted=False
		self.inmat=None
		self.shapefile=None
		self.network=None
		self.coords=None
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
		self.modavg=True
		self.report_all=False
		self.plot=True

		self.posWeight=False
		self.fixWeight=False
		self.allShapes=False
		self.fixShape=False
		self.max_shape=100

		self.min_weight=0.0

		self.only_keep=True
		self.julia="julia"
		self.compiled_modules=True
		self.sys_image=None


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
			elif opt=='s' or opt=='shp':
				self.shapefile = arg
			elif opt=="c" or opt=="coords":
				self.coords = arg
			elif opt=="n" or opt=="network":
				self.network = arg
			elif opt=="minimize":
				self.minimize=True
			elif opt=='seed':
				self.seed=int(arg)
			elif opt=='t' or opt=='procs':
				self.GA_procs=int(arg)
			elif opt=='P' or opt=='maxPop' or opt=='maxpop':
				self.maxpopsize=int(arg)
			elif opt=="G" or opt=="maxGen" or opt=="maxgen":
				self.maxGens=int(arg)
			elif opt=="popsize" or opt=="size":
				self.popsize=int(arg)
			elif opt=="minWeight":
				self.min_weight=float(arg)
			elif opt=="pop_agg":
				self.pop_agg = arg.upper()
				if self.pop_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN", "SUM", "FIRST", "SD", "VAR", "CV"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--pop_agg>")
			elif opt=="edge_agg":
				self.edge_agg = arg.upper()
				if self.edge_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN", "SUM", "FIRST", "SD", "VAR", "CV"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--edge_agg>")
			elif opt=="mutpb":
				self.mutpb=float(arg)
			elif opt=="indpb":
				self.indpb=float(arg)
			elif opt=="cxpb":
				self.cxpb=float(arg)
			elif opt=="T" or opt=="tSize" or opt=="tournSize":
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
			elif opt=="dist_col":
				self.dist_col=arg
			elif opt=="infer":
				self.predicted=True
			elif opt=="reachid_col":
				self.reachid_col=arg
			elif opt=="length_col":
				self.length_col=arg
			elif opt=="cholmod":
				self.cholmod=True
			elif opt=="v" or opt=="vars":
				self.variables=arg.split(",")
			elif opt=="V" or opt=="varFile":
				self.varFile=arg
			elif opt=="A" or opt=="modavg" or opt=="modAvg":
				self.modavg=True
			elif opt=="a" or opt=="awsum":
				self.awsum=float(arg)
			elif opt=="maxShape":
				self.max_shape=float(arg)
			elif opt=="report_all":
				self.report_all=True
			elif opt=="X" or opt=="noPlot":
				self.plot=False
			elif opt=="o" or opt=="out":
				self.out=arg
			elif opt=="avgall":
				self.only_keep=False
			elif opt=="efit_agg":
				self.efit_agg = arg.upper()
				if self.efit_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN", "SUM", "FIRST", "SD", "VAR", "CV"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--efit_agg>")
			elif opt=="posWeight":
				self.posWeight=True
			elif opt=="fixWeight":
				self.fixWeight=True
			elif opt=="allShapes":
				self.allShapes=True
			elif opt=="fixShape":
				self.fixShape=True
			elif opt == 'h' or opt == 'help':
				pass
			else:
				assert False, "Unhandled option %r"%opt

		if self.posWeight and self.fixWeight:
			self.display_help("--posWeight and --fixWeight cannot be used together")

		if self.varFile is not None:
			if self.variables is not None:
				print("Warning: Variables were specified with both <-v> and <-V>... Over-riding options using file provided with <-V>")
			if self.edge_agg is None:
				self.edge_agg="ARITH"
			with open(self.varFile) as fp:
				for line in fp:
					line=line.strip()
					stuff=line.split("\t")
					if len(stuff) < 2:
						self.agg_opts[stuff[0]]=self.edge_agg
					else:
						self.agg_opts[stuff[0]]=stuff[1]
			self.variables=list(self.agg_opts.keys())
		else:
			for v in self.variables:
				self.agg_opts[v]=self.edge_agg

		if not self.variables:
			self.display_help("No variables selected.")

	def display_help(self, message=None):
		if message is not None:
			print()
			print(message)
		print ("\nresistnet.py\n")
		print("Author: Tyler K Chafin, Biomathematics and Statistics Scotland")
		print ("Contact: tyler.chafin@bioss.ac.uk")
		print ("Description: Genetic algorithm to optimize resistance models on networks")
		print("""
Input options:
	-g,--genmat	: Genetic distance matrix
	-s,--shp	: Path to shapefile containing cleaned, contiguous stream reaches
	-c,--coords	: Input .tsv file containing sample coordinates

General options:
	--seed	: Random number seed (default=taken from clock time)
	--minimize	: Minimize input graph
	--dist_col	: Optional attribute representing edge-wise distance, e.g., "fittedD"
	--reachid_col	: Attribute name representing primary key in shapefile [default="REACH_ID"]
	--length_col	: Attribute name giving length in kilometers [default="LENGTH_KM"]
	-t,--procs	: Number of parallel processors
	-X,--noPlot	: Turn off plotting
	-o,--out	: Output file prefix
	-h,--help	: Displays help menu

Aggregation options:
	--edge_agg	: Method to use when combining variable values across segments (e.g., with --minimize)
	--efit_agg	: Method to use to aggregate variable used for calculating edge-wise fit [default=SUM]
	--pop_agg	: Method to use to combine genetic distances for populations mapping to same node
		All of these can take the following options:
		  ARITH		: [default] Use arithmetic mean
		  MEDIAN	: Use median distance
		  HARM		: Use harmonic mean
		  ADJHARM	: Adjusted harmonic mean (see docs)
		  GEOM		: Use geometric mean
		  MIN		: Use minimum distance
		  MAX		: Use maximum distance
		  FIRST		: Use first value
		  SD		: Use standard deviation
		  VAR		: Use variance
		  SUM		: Use simple sum
		  CV		: Use coefficient of variation (=SD/MEAN)
		Optionally, these may be specified individuall for each variable with -V

Genetic Algorithm Options:
	-P,--maxPop	: Maximim population size [default = 100]
	-G,--maxGen	: Maximum number of generations [default = 500]
	-s,--size	: Manually set population size to <-p int>
			    NOTE: By default, #params * 15
	-m,--mutpb	: Probability of mutation per individual [default=0.2]
	--indpb	: Probability of mutation per trait [default=0.1]
	--cxpb	: Probability of being chosen for cross-over [default=0.5]
	-T,--tSize	: Tournament size [default=10]
	--posWeight	: Constrain parameter weights to between 0.0-1.0
	--minWeight : Sets a minimum allowable weight (only valid when --posWeight)
	--fixWeight	: Constrain parameter weights to 1.0 (i.e., unweighted)
	--fixShape	: Turn off feature transformation
	--allShapes	: Allow inverse and reverse transformations
	--maxShape	: Maximum shape value (all transformations approach linear as shape increases) [default=100]

Model optimization/ selection options:
	-v,--vars	: Comma-separated list (no spaces) of explanatory attributes to include
	-V,--varfile	: Optional file with variables provided like so:
		  var1 \t <Optional aggregator function>
		  var2 \t <Optional aggregator function>
		  ...
		  ...
	-F,--nfail	: Number of generations failing to improve to stop optimization
	-d,--delt	: Threshold absolute change in fitness [default=0.0]
	-D,--deltP	: Threshold percentage change in fitness, as decimal [default=0.001]
	-f,--fit	: Fitness metric used to evaluate models 
			    Options:
			    aic (default)
			    loglik (log-likelihood)
			    r2m (marginal R^2)
			    delta (Change in AIC versus null model)
			    NOTE: Case-insensitive
	-b,--burn	: Number of generations for pre-burnin [default=0]
	--max_hof_size	: Maximum individuals to track in the Hall of Fame [default=100]
	--null	: Output null (population-only) model metrics ***NOT IMPLEMENTED***

Multi-model inference options:
	-a,--awsum	: Cumulative Akaike weight threshold to retain top N models [default=0.95]
	--report_all	: Plot per-stream resistance and generate full outputs for all retained models

""")
		print()
		sys.exit()
