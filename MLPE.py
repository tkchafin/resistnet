import sys
import itertools
import numpy as np
import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#import statsmodels.regression.mixed_linear_model as mlm
#from statsmodels.tools.sm_exceptions import ConvergenceWarning
#from statsmodels.regression.mixed_linear_model import VCSpec
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

#input should be two matrices: 
def MLPE_R(X, Y, scale=True):
	x = get_lower_tri(X)
	y = get_lower_tri(Y)
	
	#make ID table
	npop=X.shape[0]
	ID = to_from_(npop) #nrow(matrix)
	#print(ID)
	
	#make ZZ matrix
	ZZ = pd.DataFrame(ZZ_mat_(npop, ID))
	#print(ZZ)

	#center and scale y
	if scale==True:
		y = (y-np.mean(y))/np.std(y)
	
	#add x and y to data
	ID["x"] = x
	ID["y"] = y
	print(ID)
	
	r=ro.r
	r['options'](warn=-1)
	r['source']("resistanceGA_MLPE.R")
		
	mlpe=ro.globalenv['MLPE']
	with localconverter(ro.default_converter + pandas2ri.converter):
		mlpe_res=mlpe(ID, ZZ)
		
	print(mlpe_res)

def get_lower_tri(mat):
	n=mat.shape[0]
	i=np.tril_indices(n,-1)
	return(mat[i])

def testSM():
	x = np.array([2.6407333, 0.6583007, 1.9629121, 0.8529997, 2.2592001, 2.9629032, 2.0796441, 2.4179196, 0.2154603, 2.5016938])
	y = np.array([3.6407333, 1.6583007, 1.5629121, 0.4529997, 2.0592001, 2.0629032, 2.9796441, 3.1179196, 1.2154603, 1.5016938])
	
	#make ID table
	ID = to_from_(5) #nrow(matrix)
	print(ID)
	
	#make ZZ matrix
	ZZ = ZZ_mat_(5, ID)
	print(ZZ)
	
	#center and scale x and y
	x = (x-np.mean(x))/np.std(x)
	y = (y-np.mean(y))/np.std(y)
	
	#add x and y to data
	ID["x"] = x
	ID["y"] = y
	print(ID)
	
	#get matrix describing random effects structure
	vcs = getVCM(ID)
	
	#fit mixed model
	#equivalent to lme4 y ~ x + (1|pop1)
	#how do incorporate the 'ZZ' part??
	model = smf.mixedlm("y ~ x", data=ID, groups="pop1")
	#model = mlm.MixedLM(endog=ID["y"].to_numpy(), exog=ID["x"].to_numpy(), groups=ID["pop1"].to_numpy(), exog_re=ZZ)
	#model = mlm.MixedLM(endog=ID["y"].to_numpy(), exog=ID["x"].to_numpy(), groups=ID["pop1"].to_numpy(), exog_vc=vcs)
	model_fit = model.fit()
	print(model)
	print(model_fit)
	print(model_fit.summary())

def getVCM(df):
	def f(x):
		n = x.shape[0]
		g2 = x["pop2"]
		u = g2.unique()
		u.sort()
		uv = {v: k for k, v in enumerate(u)}
		mat = np.zeros((n, len(u)))
		for i in range(n):
		    mat[i, uv[g2[i]]] = 1
		colnames = ["%d" % z for z in u]
		print(mat)
		return(mat, colnames)
	vcm = df.groupby("pop1").apply(f).to_list()
	print(vcm)
	mats = [x[0] for x in vcm]
	print(mats)
	colnames = [x[1] for x in vcm]
	names = ["pop2"]
	vcs = VCSpec(names, [colnames], [mats])
	return(vcs)

def to_from_(pops):
	to = list()
	frm = list()
	for ia, ib in itertools.combinations(range(1,pops+1),2):
		to.append(ia)
		frm.append(ib)
	t = to[pops-2]
	tt = frm[pops-2]
	to[pops-2] = tt
	frm[pops-2] = t
	df = pd.DataFrame({"pop1" : to, "pop2" : frm})
	return(df)

def ZZ_mat_(pops, id):
	zz = np.zeros(shape=(pops,id.shape[0]))
	for i in range(0, pops):
		#print(i)
		for j in range(0, id.shape[0]):
			#if ID row j contains pop i+1, set zz[j, i] to 1
			if i+1 in list(id.iloc[j]):
				#print("  ",j)
				zz[i,j]=1
	return(zz)

