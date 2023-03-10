import sys
import os
import itertools
import pandas as pd
import numpy as np
import pandas as pd
from collections import OrderedDict
from io import StringIO

import resistnet.MLPE as mlpe_rga

def parsePairwise(points, inc_matrix, multi, gendist):
	#print("parsePairwise")
	#print("evaluate")
	#print(gendist, flush=True)
	r=effectiveResistanceMatrix(points, inc_matrix, multi)
	#print(r)
	res = mlpe_rga.MLPE_R(gendist, r, scale=True)
	return(r, res)

# def parsePairwise(r, gendist, return_resistance=False):
# 	res = mlpe_rga.MLPE_R(gendist, r, scale=True)
# 	return(res)

def effectiveResistanceMatrix(points, inc_matrix, edge_resistance):
	r=pd.DataFrame(columns=list(points.values()), index=list(points.values()))
	inc_row=0
	edge_resistance=np.array(edge_resistance)
	for ia, ib in itertools.combinations(range(0,len(points)),2):
		inc=np.array(inc_matrix[inc_row,])
		d=np.sum(np.multiply(edge_resistance, inc)) #effective resistance is simply a sum of serial resistances
		inc_row+=1
		r.loc[list(points.values())[ia], list(points.values())[ib]] = d
		r.loc[list(points.values())[ib], list(points.values())[ia]] = d
	r=r.astype('float64').to_numpy()
	np.fill_diagonal(r, 0.0)
	return(r)
