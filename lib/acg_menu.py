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
		self.shapefile = None
		self.locmatdir = None
		self.edgeid= "EDGE_ID"


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
	Mandatory arguments:
		-n,--network	: Input .network file from autoStreamTree
		-i,--input	: Input .streamtree.txt file from autoStreamTree

	General options:
		-o,--out	: Output prefix [default="out"]
		-b,--boot	: Number of bootstraps <NOT IMPLEMENTED>
		-f,--fit	: One of: aic, loglik, r2 [default=AIC] <NOT IMPLEMENTED>
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