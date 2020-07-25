#!/usr/bin/python

import sys
import os
import getopt
from sklearn.cluster import DBSCAN
import pandas as pd

def main():
	params = parseArgs()
	
	
	

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hm:o:', \
			["help", "matrix=", "out="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.matrix=None


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
			elif opt=="matrix" or opt=="m":
				self.matrix=arg
			elif opt=="o" or opt=="out":
				self.out=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.matrix:
			self.display_help("No matrix provided (-m,--matrix)")
		if not self.out:
			self.out = "out"



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nclusterPopsDB.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Script for assigning samples to clusters given an arbitrary distance matrix")
		print("""
		-m,--matrix	: Input labelled matrix (column AND row labels), tab-delimited
		-o,--out	: Output prefix (if not overwriting original)
		-h,--help	: Displays this help menu)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
