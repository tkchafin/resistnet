import os
import sys
import numpy as np
import scipy


def aggregateDist(method, stuff):
	if method == "HARM":
		try:
			return(scipy.stats.hmean(stuff))
		except ValueError as e:
			print(e)
			print("ERROR (DivideByZero): Harmonic mean cannot be calculated using a zero distance. Try recomputing using the \"ADJHARM\" option.")
			print("")
			sys.exit(1)
	elif method == "ARITH":
		return(np.mean(stuff))
	elif method == "GEOM":
		return(scipy.stats.mstats.gmean(stuff))
	elif method == "MEDIAN":
		return(np.median(stuff))
	elif method == "MAX":
		return(np.max(stuff))
	elif method == "MIN":
		return(np.min(stuff))
	elif method == "ADJHARM":
		return(adjustedHarmonicMean(stuff))
	elif method == "SD":
		return(np.std(stuff))
	elif method == "VAR":
		return(np.var(stuff))

#computes an harmonic mean corrected for non-positive values
def adjustedHarmonicMean(stuff):
	s=np.array(stuff)
	vals = s[s>0.0]
	bads = s[s<=0.0]
	mu = (1.0 / (np.sum([1.0/x for x in vals]) / (len(vals)-len(bads)))) * ((len(vals)-len(bads))/len(vals))
	return(mu)