import sys
import os, glob
import geopandas as gpd
import networkx as nx
import getopt
import pickle
import momepy

import resistnet.resist_dist as rd
import resistnet.stream_plots as splt


def main():
	params = parseArgs()

	# get graph from shapefile
	G = read_network(params.network, None)

	print(G.nodes)


# read network
def read_network(network, shapefile):
	if network:
		print("Reading network from saved file: ", network)
		G=nx.Graph(nx.read_gpickle(network).to_undirected())
	else:
		print("Building network from shapefile:",shapefile)
		print("WARNING: This can take a while with very large files!")
		rivers = gpd.read_file(shapefile)
		G=momepy.gdf_to_nx(rivers, approach="primal", directed=False, multigraph=False)
		#G=nx.Graph(nx.read_shp(shapefile, simplify=True, strict=True).to_undirected())
	return(G)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'ho:i:n:s:r:', \
			["help", "out=", "in=", "network=", "reps=", "samples="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.input=None
		self.out="sim"
		self.network=None
		self.reps=1
		self.samples=50


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "h" or opt == "help":
				continue
			elif opt=="in" or opt=="i":
				self.input=arg
			elif opt=="out" or opt=="o":
				self.out=arg
			elif opt == "network" or opt=="n":
				self.network=arg
			elif opt == "samples" or opt=="s":
				self.samples=int(arg)
			elif opt == "reps" or opt=="r":
				self.reps=int(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.input:
			self.display_help("No input table provided.")
		if not self.network:
			self.display_help("No network provided.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nsimResistnet.py\n")
		print("Author: Tyler Chafin")
		print ("Description: Simulate data on a given network for validating resistnet")
		print("""
Arguments:
-n,--network	: Input network (pickle'd networkx output)
-i,--in		: Table giving information on which variables to use to generate resistnet input
-o,--out	: Output file name (default=out.fas)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
