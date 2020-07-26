#!/usr/bin/python

import sys
import os
import getopt
from sklearn.cluster import DBSCAN
import pandas as pd
from sortedcontainers import SortedDict

def main():
	params = parseArgs()
	
	print("\nReading in distance matrix:",params.matrix)
	mat = pd.read_csv(params.matrix, header=0, index_col=0, sep="\t")
	
	#print(mat.to_numpy())
	print("\nAssigning populations with DBSCAN using:")
	print("\tEpsilon=",str(params.epsilon))
	print("\tMin samples per cluster:",str(params.min_samples))
	print("\tAlgorithm:",params.algorithm)
	print("N processes:",int(params.procs))
	db=DBSCAN(eps=params.epsilon, min_samples=params.min_samples, algorithm=params.algorithm, n_jobs=params.procs, metric='precomputed').fit(mat.to_numpy())

	#get cluster labels
	cluster_labels=db.labels_

	#build SortedDict to return
	popmap=SortedDict()
	i=0
	for k in mat.columns.values:
		pop="DB_"+str(cluster_labels[i])
		if pop not in popmap:
			l = [k]
			popmap[pop] = l
		else:
			popmap[pop].append(k)
		i+=1

	flat = flattenPopmap(popmap)
	
	b=len(mat.columns.values)
	e=len(set(list(flat.values())))
	print("\nDone! Assigned",str(b),"samples to",e,"clusters!\n")
	
	o=str(params.out)+".popmap.txt"
	print("\nWriting output to:",o)
	temp = pd.DataFrame(list(flat.items()), columns=['IND_ID', 'POP_ID'])
	temp.to_csv(o, sep="\t", index=False)
	print()
	

#utility function, converts popmap of form key=pop; value=list(inds) to key=ind; value=pop
def flattenPopmap(popmap):
	new_popmap=dict()
	for k in popmap:
		for i in popmap[k]:
			new_popmap[i]=k
	return(new_popmap)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hm:o:M:e:a:p:', \
			["help", "matrix=", "out=", "epislon=", "min_samples=", "algorithm=", "procs="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.matrix=None
		self.out=None
		self.min_samples=1
		self.epsilon=20
		self.algorithm="auto"
		self.procs=1


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
			elif opt=="M" or opt=="min_samples":
				self.min_samples=int(arg)
			elif opt=="e" or opt=="epsilon":
				self.epsilon=float(arg)
			elif opt=="p" or opt=="procs":
				self.procs=int(arg)
			elif opt=="a" or opt=="algorithm":
				self.algorithm=arg
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
	Arguments:
		-m,--matrix	: Input labelled matrix (column AND row labels), tab-delimited
		-o,--out	: Output prefix 
		-M,--min_samples	: Minimum samples per cluster [default=1]
		-e,--epsilon	: Maximum distance (same units as dist matrix) within a cluster [default=20]
		-a,--algorithm	: auto [default], ball_tree, kd_tree, brute (see NearestNeighbors sklearn)
		-p,--procs	: Number of processors to use
		-h,--help	: Displays this help menu)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
