import sys
import os
import itertools
import pandas as pd
import numpy as np
import pandas as pd
from collections import OrderedDict
from io import StringIO 

import riverscape.MLPE as mlpe_rga


def parseEdgewise(r, edge_gendist, return_resistance=False):
	l=["from", "to", "r"]
	output_file=str(oname)+"_resistances_3columns.out"
	input=pd.read_csv((str(oname)+".graph_resistances.txt"), header=None, index_col=None, sep="\t", names=l)
	output=pd.read_csv((str(oname)+"_resistances_3columns.out"), header=None, index_col=None, sep=" ", names=l)
	output["from"] = output["from"]-1
	output["to"] = output["to"]-1
	merged=pd.merge(input, output, how="left", on=["from", "to"])
	print(merged)
	
def parsePairwise(order, points, inc_matrix, multi, gendist):
	r=effectiveResistanceMatrix(order, points, inc_matrix, multi)
	#print(rr)
	res = mlpe_rga.MLPE_R(gendist, r, scale=True)
	#print(res)
	return(res)
	
# def parsePairwise(r, gendist, return_resistance=False):
# 	res = mlpe_rga.MLPE_R(gendist, r, scale=True)
# 	return(res)

#function to write inputs for circuitscape
def effectiveResistanceMatrix(order, points, inc_matrix, distances):
	r=pd.DataFrame(columns=list(points.values()), index=list(points.values()))
	inc_row=0
	distances=np.array(distances)
	for ia, ib in itertools.combinations(range(0,len(points)),2):
		inc_streams=np.array(inc_matrix[inc_row,])
		d=np.sum(np.multiply(distances, inc_streams)) #effective resistance is simply a sum of serial resistances
		#note this won't work with a network that isn't strictly bifucating.
		inc_row+=1
		r.loc[list(points.values())[ia], list(points.values())[ib]] = d
		r.loc[list(points.values())[ib], list(points.values())[ia]] = d
	r=r.astype('float64')
	formatted=r.reindex(order)
	formatted=formatted[order]
	formatted=formatted.to_numpy()
	return(formatted)
	