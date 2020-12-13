#!/usr/bin/python

import sys
import os
import getopt
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt

def main():
	params = parseArgs()
	
	s=str(params.prefix) + ".streamTree.shp"
	print("Reading shapefile:",s)
	gdf = gpd.read_file(s)
	
	if params.maxD:
		gdf.loc[gdf.fittedD > params.maxD, 'fittedD']=params.maxD
	if params.minD:
		gdf.loc[gdf.fittedD < params.minD, 'fittedD']=params.minD
		
	gdf.plot(column="fittedD", cmap = params.cmap, legend=True)
	plt.title("Stream network colored by StreamTree fitted distances")
	o=(str(params.out)+".streamsByFittedD.pdf")
	plt.savefig(o)
	print("New plot can be found at:",o)
	
	

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hp:m:M:c:o:', \
			["help", "prefix=", "cmap=", "min=", "max=", "out="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.prefix=None
		self.maxD = None
		self.minD = None
		self.cmap = "RdYlGn_r"
		self.out = None


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
			elif opt=="prefix" or opt=="p":
				self.prefix=arg
			elif opt=="c" or opt=="cmap":
				self.cmap = arg
			elif opt=="m" or opt=="min":
				self.minD = float(arg)
			elif opt=="M" or opt=="max":
				self.maxD = float(arg)
			elif opt=="o" or opt=="out":
				self.out=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.prefix:
			self.display_help("No prefix provided (-p,--prefix)")
		if not self.out:
			self.out = self.prefix



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nplotStreamTree.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Script for re-plotting StreamTree results after running autoStreamTree")
		print("""
		-p,--prefix	: Prefix for autoStreamTree output
		-m,--min	: Minimum genetic distance 
		-M,--max	: Maximum genetic distance
		-c,--cmap	: Colormap (any valit matplotlib cmap value)
			see: https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html
		-o,--out	: Output prefix (if not overwriting original)
		-h,--help	: Displays this help menu)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
