import os
import sys
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import warnings
warnings.simplefilter('ignore', category=UserWarning)
pd.options.mode.chained_assignment = None 

import riverscape.MLPE as MLPE

class hallOfFame():
	def __init__(self, variables, max_size, init_pop=None):
		cols=list()
		cols.append("fitness")
		for v in variables:
			cols.append(str(v))
			cols.append(str(v)+"_weight")
			cols.append(str(v)+"_trans")
			cols.append(str(v)+"_shape")
		cols.append("loglik")
		cols.append("r2m")
		cols.append("aic")
		cols.append("delta_aic_null")
		self.data = pd.DataFrame(columns=cols)
		self.variables=variables
		self.max_size=int(max_size)
		self.min_fitness=float('-inf')
		self.rvi=pd.DataFrame(columns=['variable', 'RVI'])
		self.zero_threshold=0.0000000000000001
		
		if init_pop is not None:
			self.check_population(init_pop)
			
	def check_population(self, pop):
		popDF = pd.DataFrame(pop, columns=self.data.columns)
		popDF = popDF[popDF.fitness > float('-inf')]
		if popDF.shape[0] < 1:
			return
		popDF = popDF.sort_values('fitness', ascending=False)
		popDF = popDF.drop_duplicates(keep='first', ignore_index=True)
		popDF = popDF.reset_index(drop=True)
		space = self.max_size - self.data.shape[0]
		
		if space > 0:
			#print('hall of fame not full')
			select_size=space
			if space> popDF.shape[0]:
				select_size=popDF.shape[0]
			self.data = self.data.append(popDF[:select_size], ignore_index=True)
			self.data = self.data.sort_values('fitness', ascending=False)
			self.data = self.data.drop_duplicates(keep='first', ignore_index=True)
			self.data = self.data.reset_index(drop=True)
			self.min_fitness = self.data['fitness'].min()
		else:
			if popDF['fitness'].max() > self.min_fitness:
				subset=popDF[popDF.fitness > self.min_fitness]
				self.data = self.data.append(subset, ignore_index=True)
				self.data = self.data.sort_values('fitness', ascending=False)
				self.data = self.data.drop_duplicates(keep='first', ignore_index=True)
				self.data = self.data.reset_index(drop=True)
				if self.data.shape[0] > self.max_size:
					self.data = self.data[:self.max_size]
				self.min_fitness = self.data['fitness'].min()
			else:
				return
	
	def printHOF(self, max_row=None, max_col=None):
		self.data = self.data.sort_values('fitness', ascending=False)
		self.data = self.data.reset_index(drop=True)
		with pd.option_context('display.max_rows', max_row, 'display.max_columns', max_col):  # more options can be specified also
			print(self.data)
	
	def printRVI(self, max_row=None, max_col=None):
		self.rvi = self.rvi.sort_values('RVI', ascending=False)
		self.rvi = self.rvi.reset_index(drop=True)
		with pd.option_context('display.max_rows', max_row, 'display.max_columns', max_col):  # more options can be specified also
			print(self.rvi)
	
	def delta_aic(self):
		if self.data.shape[0] <= 0:
			return
		if "delta_aic_best" in self.data.columns:
			return
		else:
			self.data["aic"] = self.data["aic"]*-1 #reverse the neg sign i added for maximizing
			best=self.data["aic"].min()
			self.data["delta_aic_best"] = self.data["aic"]-best

	
	def akaike_weights(self):
		if self.data.shape[0] <= 0:
			return
		if "delta_aic_best" not in self.data.columns:
			self.delta_aic()
		#weight(i) = e^(-1/2 delta_aic_best) / sum(e^(-1/2 delta_aic_best(k)))
		#where the denominator is summed over k models
		#delta_aic = self.data["delta_aic_best"].to_numpy()
		#sum_k = self.data["delta_aic_best"].to_numpy()
		#did a test agains MuMIn::Weights in R and this seems to be working
		self.data["akaike_weight"]=((np.exp(-0.5*self.data["delta_aic_best"])) / (sum(np.exp(-0.5*self.data["delta_aic_best"]))))
		self.data["akaike_weight"].mask(self.data["akaike_weight"] < self.zero_threshold, 0.0, inplace=True)
	
	def cumulative_akaike(self, threshold=1.0):
		if self.data.shape[0] <= 0:
			return
		self.data = self.data.reset_index(drop=True)
		threshold=float(threshold)
		if "akaike_weight" not in self.data.columns:
			self.akaike_weights()
		self.data["acc_akaike_weight"] = self.data["akaike_weight"].cumsum()
		if threshold > 0.0 and threshold < 1.0:
				if self.data["acc_akaike_weight"].max() > threshold:
					cutoff=self.data[self.data["acc_akaike_weight"].gt(threshold)].index[0]
					keep_vals = ["False"]*self.data.shape[0]
					keep_vals[:(cutoff+1)] = ["True"]*(cutoff+1)
					self.data["keep"] = keep_vals #+1 b/c above is 0-based index
				else:
					keep_vals = ["True"]*self.data.shape[0]
					self.data["keep"] = keep_vals
		else:
			keep_vals = ["True"]*self.data.shape[0]
			self.data["keep"] = keep_vals

	def relative_variable_importance(self,ignore_keep=False):
		#clear previous calculations
		self.rvi=pd.DataFrame(columns=['variable', 'RVI'])
		sub=self.data[self.data.keep=="True"]
		if ignore_keep:
			sub=self.data
		#compute sum of weights
		for v in self.variables:
			sw=(sub[v]*sub['akaike_weight']).sum()
			self.rvi.loc[len(self.rvi), :] = [v, sw]
		self.rvi = self.rvi.sort_values('RVI', ascending=False)
		self.rvi = self.rvi.reset_index(drop=True)
	
	def cleanHOF(self):
		ret = self.data.sort_values('fitness', ascending=False)
		ret = ret.reset_index(drop=True)
		for v in self.variables:
			mask=(ret[v] == 0)
			ret[(str(v)+"_weight")][mask] = "NaN"
			ret[(str(v)+"_trans")][mask] = "NaN"
			ret[(str(v)+"_shape")][mask] = "NaN"
		if ret.iloc[0]["fitness"] == (ret.iloc[0]["aic"]*-1):
			ret["fitness"] = ret["fitness"]*-1
		return(ret)
	
	def getRVI(self):
		self.rvi = self.rvi.sort_values('RVI', ascending=False)
		self.rvi = self.rvi.reset_index(drop=True)
		return(self.rvi)
	
	def getHOF(self, only_keep=False):
		self.data = self.data.sort_values('fitness', ascending=False)
		self.data = self.data.reset_index(drop=True)
		if only_keep:
			return(self.data[self.data.keep=="True"])
		else:
			return(self.data)
	
	def output(self):
		#Make sure to remove weights/ shapes where variable isn't selected
		#get absolute value of AIC (made them negative so all metrics could use maximize function)
		pass
	
	def plot_ICprofile(self, oname="out", diff=2):
		diff=int(diff)
		#X axis - order by AIC.
		dat=self.data.sort_values('aic', ascending=True)
		dat = dat.round(3)
		dat.reset_index()
		dat["model"]=dat.index + 1
		#y axis - AIC value
		sns.set(style="ticks")
		p = sns.scatterplot(data=dat, x="model", y="aic", hue="r2m", size="r2m", style="keep",  alpha=0.6)
		p.axhline((dat["aic"].min()+diff), ls="--", c="red")
		plt.title("IC Profile")
		plt.savefig((str(oname)+".ICprofile.pdf"))
		plt.clf()
		plt.close()
	
	def plotMetricPW(self, oname="out"):
		cols=["aic", "loglik", "r2m", "delta_aic_null", "keep"]
		if "akaike_weight" in self.data.columns:
			cols.append("akaike_weight")
		dat=self.data[cols]
		sns.set(style="ticks")
		sns.pairplot(dat, hue="keep", kind="scatter")
		plt.title("Pair plot")
		plt.savefig((str(oname)+".pairPlot.pdf"))
		plt.clf()
		plt.close()
	
	def plotVariableImportance(self, oname="out", cutoff=0.8):
		cutoff=float(cutoff)
		sns.set(style="ticks")
		sub=self.rvi.sort_values('RVI', ascending=False)
		p=sns.barplot(data=sub, x="RVI", y="variable")
		p.axvline(cutoff, ls="--", c="red")
		plt.title("Importance of Terms")
		plt.savefig((str(oname)+".varImportance.pdf"))
		plt.clf()
		plt.close()
	
	def writeModelSummary(self, oname):
		out_df = self.cleanHOF()
		out_df.to_csv((str(oname)+".HallOfFame.tsv"), sep="\t", index=False, na_rep="-")

def plotEdgeModel(gen, res, oname):
	sns.set(style="ticks")
	df = pd.DataFrame(list(zip(gen, res)), columns=["Fitted Genetic Distance", "Resistance"])
	#print(df)
	sns.lmplot(x="Resistance", y="Fitted Genetic Distance", data=df)
	plt.title("Edge-wise Resistance x Genetic Distance")
	plt.savefig((str(oname)+".Edgewise.pdf"))
	plt.clf()
	plt.close()

def plotPairwiseModel(gen, mat, oname, partition=False):
	g=MLPE.get_lower_tri(gen)
	r=MLPE.get_lower_tri(mat)
	#print(len(g))
	#print(len(r))
	df = pd.DataFrame(list(zip(list(g), list(r))), columns=["Genetic Distance", "Resistance Distance"])
	
	sns.set(style="ticks")
	if partition:
		npop=gen.shape[0]
		ID = MLPE.to_from_(npop)
		df["grp"] = ID["pop1"]
		sns.lmplot(x="Resistance Distance", y="Genetic Distance", hue="grp", data=df)
	else:
		sns.lmplot(x="Resistance Distance", y="Genetic Distance", data=df)
	plt.title("Pairwise-wise Resistance x Genetic Distance")
	plt.savefig((str(oname)+".Pairwise.pdf"))
	plt.clf()
	plt.close()


