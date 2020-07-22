import sys
import os
import getopt

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hs:i:r:pd:a:lw:o:gP:L:Scn:', \
			["shp=", "help", "input=", "run=", "pop", "pops","dist=", "agg_method=",
			"het", "genmat=", "snp", "snps", "msat", "msats", "log", "and_log", "iterative",
			"weight=", "out=", "method=", "plots", "plot","perm=", "phased", "median",
			"diploid", "geopop", "geopops", "global_het", "haploid", "loc_agg=", 
			"pop_agg=", "sdist_agg=", "clusterpop", "epsilon=", "min_samples=", "sclusterpop",
			"network=", "overwrite"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.shapefile = None
		self.geodb = None
		self.run = "ALL"
		self.network = None
		self.pop = False
		self.overwrite = False
		self.geopop = False
		self.clusterpop=False
		self.sclusterpop=False
		self.dist = "JC69"
		self.het = False
		self.genmat = None
		self.snps = False
		self.msat = False
		self.log = False
		self.and_log = False
		self.iterative = False
		self.weight = "CSE67"
		self.permutations = 1000
		self.method = "PEARSON"
		self.plots=False
		self.out="out"
		self.median=False
		self.ploidy=2
		self.global_het=False
		self.loc_agg = "ARITH"
		self.pop_agg = "ARITH"
		self.sdist_agg="ARITH"
		
		#dbscan Options
		self.min_samples=1
		self.epsilon=20


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
			elif opt == 'i' or opt == 'input':
				self.geodb = arg
			elif opt == 'r' or opt == 'run':
				self.run = arg.upper()
				if self.run not in ["ALL", "GENDIST", "IBD", "STREAMDIST", "STREAMTREE", "DISTANCES", "EXHAUSTIVE", "REACHFIT"]:
					self.display_help("Invalid option", arg.upper(),"for option <-r/--run>")
			elif opt == 'p' or opt == 'pop' or opt == "pops":
				self.pop = True
			elif opt == "g" or opt == "geopop" or opt == "geopops":
				self.geopop = True
			elif opt == 'd' or opt == 'dist':
				self.dist = arg.upper()
				if self.dist not in ["PDIST", "JC69", "K2P", "TN84", "TN93", "FST", "LINFST", "JOST", "NEI83", "EUCLID", "CHORD", "GST", "GSTPRIME", "LINJOST"]:
					self.display_help("Invalid option", arg.upper(),"for option <-d/--dist>")
			elif opt == "het":
				self.het = True
			elif opt == "genmat":
				self.genmat = arg
			elif opt == "snp" or opt == "snps":
				self.snps = True
			elif opt == "msat" or opt == "msats":
				self.msat = True
				print("WARNING: Option <--msat> not yet implemented")
				sys.exit(0)
			elif opt == "l" or opt == "log":
				self.log = True
			elif opt == "and_log":
				self.and_log = True
			elif opt == "iterative":
				self.iterative = True
			elif opt=="clusterpop" or opt=="c":
				self.clusterpop=True
			elif opt=="sclusterpop" or opt=="S":
				self.sclusterpop=True
			elif opt=="epsilon":
				self.epsilon = float(arg)
			elif opt=="min_samples":
				self.min_samples=int(arg)
			elif opt == "w" or opt == "weight":
				self.weight = arg.upper()
				if self.weight == "FM" or self.weight=="1/D":
					self.weight="FM67"
				if self.weight == "BEYER" or self.weight == "1/D^2":
					self.weight="BEYER74"
				if self.weight=="1" or self.weight=="CSE":
					self.weight="CSE67"
				if self.weight not in ["CSE67", "FM67", "BEYER74"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <-w/--weight>")
			elif opt == "o" or opt == "out":
				self.out = arg
			elif opt == "stream_fit":
				self.stream_fit = True
			elif opt == "overwrite":
				self.overwrite=True
			elif opt == "perm":
				self.permutations = int(arg)
			elif opt == "method":
				print("Sorry: Option --method is not yet implemented.")
				sys.exit(0)
				self.method = arg.upper()
				if self.method not in ["PEARSON", "SPEARMAN", "BOTH"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--method>")
			elif opt == "plot" or opt == "plots":
				self.plots = True
			elif opt == "phased":
				self.phased = True
				print("WARNING: Option <--snp> not yet implemented")
			elif opt=="median":
				self.median=True
			elif opt=="diploid":
				self.ploidy=2
			elif opt=="haploid":
				self.ploidy=1
			elif opt=="global_het":
				self.global_het=True
			elif opt=="network" or opt=="n":
				self.network=arg
			elif opt=="pop_agg" or opt=="P":
				self.pop_agg = arg.upper()
				if self.pop_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--pop_agg>")
			elif opt=="loc_agg" or opt=="L":
				self.loc_agg = arg.upper()
				if self.loc_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--loc_agg>")
			else:
				assert False, "Unhandled option %r"%opt

		if not self.geodb:
			self.display_help("No input provided <-i/--input>")
		if not self.shapefile and self.run != "GENDIST":
			self.display_help("No shapefile provided <-s/--shp>")
		if self.snps and self.dist in ["PDIST", "JC69", "K2P", "TN93"]:
			print("\nWARNING: Can't calculate",str(self.dist),"distances with SNP data. Will instead report proportion of SNPs which are different.")
			print("NOTE: This is equivalent to an uncorrected p-distance (PDIST) calculated on the concatenated SNP data.\n")
			self.snps = False
			self.dist = "PDIST"
		if self.ploidy > 2 or self.ploidy < 1:
			self.display_help("Ploidy of",self.ploidy,"not currently allowable. Please choose 1 (haploid) or 2 (diploid)")

		#sanity checks
		if self.dist not in ["PDIST", "TN84", "TN93", "K2P", "JC69"]:
			if not self.pop and not self.geopop:
				print("ERROR: Distance metric",self.dist,"not possible without --pop or --geopop data.")
				sys.exit(1)

		###DIE FOR OPTIONS NOT YET IMPLEMENTED
		if self.run in ["EXHAUSTIVE", "REACHFIT"]:
			print("Sorry: Option --dist",self.run," not yet implemented.")
			sys.exit(0)
		if self.sclusterpop:
			print("Sorry: Option --sclusterpop not yet implemented.")
			sys.exit(0)



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nautoStreamtree.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Computes stream distances and genetic distances for \
georeferenced DNA sequences, performs tests for isolation-by-distance, \
and uses a least-squares method to fit distances to stream segments.")
		print("""
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
			ALL			: Run all steps
			GENDIST		: Only calculate genetic distance matrix
			STREAMDIST	: Only compute pairwise stream distances
			DISTANCES	: Only compute GENDIST + STREAMDIST
			IBD		: GENDIST + STREAMDIST + Mantel test
			STREAMTREE	: GENDIST + STREAMDIST + fit StreamTree model
			xxx = NOT YET IMPLEMENTED
		-p,--pop		: Pool individuals based on column 2 of input file
			NOTE: The location will be taken as the centroid among individual samples
		-g,--geopop		: Pool individuals having identical coordinates
		-c,--clusterpop	: Use DBSCAN algorithm to automatically cluster populations

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
		--genmat	: Skip calculation and use the provided labeled .tsv matrix
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

	StreamTree (see Kaliowski et al. 2008) options:
		--iterative	: Prevent negative distances using the iterative approach
		-w,--weight	: Desired weighting for least-squares fitting:
			Options:
			  FM67			: Fitch and Margoliash (1967) [w = 1/D^2]
			  BEYER74		: Beyer et al. (1974) weights [w = 1/D]
			  CSE67			: [default] Cavalli-Sforza and Edwards (1967) [w = 1]
""")
		print()
		sys.exit()