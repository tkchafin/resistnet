import sys
import itertools
import numpy as np
import pandas as pd

"""
Transformations from Peterman et al. (2018) ResistanceGA package
Full Citation:
Peterman, W. E. 2018. ResistanceGA: An R package for the optimization 
of resistance surfaces using genetic algorithms. Methods in Ecology 
and Evolution doi:10.1111/2041-210X.12984
GitHub:
https://github.com/wpeterman/ResistanceGA/tree/3be50a6fa5515fe23953fca39135c5b340415607
"""

#rescales all columns in a pandas object between m and M, where M > m >= 0
def rescaleCols(df, m, M):
	df -= df.min()
	df /= df.max()
	return((df*(M-m))+m)

def ricker(dat, shape, ceiling):
	#R: parm[3]*r*exp(-1*r/parm[2])+1
	#parm[2] = shape
	#parm[3] = ceiling
	return(ceiling*dat*np.exp(-1*dat/shape)+1)

def invRicker(dat, shape, ceiling):
	return((-1*ceiling)*dat*np.exp(-1*dat/shape)-1)
	
def revInvRicker(dat, shape, ceiling):
	d = invRicker(dat, shape, ceiling)
	return(rescaleCols((-1*d), min(d), max(d)))
	#return(invRicker((-1*dat), shape, ceiling))

def revRicker(dat, shape, ceiling):
	d = rescaleCols((-1*dat), min(dat), max(dat))
	return(ricker(d, shape, ceiling))

def monomolecular(dat, shape, ceiling):
	return(ceiling*(1-np.exp(-1*dat/shape))+1)
	
def invMonomolecular(dat, shape, ceiling):
	d = ceiling*np.exp(-1*dat/shape)
	return((d-min(d))+1)

def revInvMonomolecular(dat, shape, ceiling):
	d = rescaleCols((-1*dat), min(dat), max(dat))
	return(invMonomolecular(d, shape, ceiling))

def revMonomolecular(dat, shape, ceiling):
	d = rescaleCols((-1*dat), min(dat), max(dat))
	return(monomolecular(d, shape, ceiling))

	