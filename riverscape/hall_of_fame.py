import os
import sys
import pandas as pd 
import numpy as np

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
		self.max_size=int(max_size)
		self.min_fitness=float('-inf')
		
		if init_pop is not None:
			self.check_population(init_pop)
			
	def check_population(self, pop):
		popDF = pd.DataFrame(pop, columns=self.data.columns)
		popDF = popDF[popDF.fitness > float('-inf')]
		if popDF.shape[0] < 1:
			return
		popDF = popDF.sort_values('fitness', ascending=False)
		popDF = popDF.drop_duplicates(keep='first', ignore_index=True)
		space = self.max_size - self.data.shape[0]
		
		if space > 0:
			#print('hall of fame not full')
			select_size=space
			if space> popDF.shape[0]:
				select_size=popDF.shape[0]
			self.data = self.data.append(popDF[:select_size], ignore_index=True)
			self.min_fitness = self.data['fitness'].min()
		else:
			if popDF['fitness'].max() > self.min_fitness:
				subset=popDF[popDF.fitness > self.min_fitness]
				self.data = self.data.append(subset, ignore_index=True)
				self.data = self.data.sort_values('fitness', ascending=False)
				self.data = self.data.drop_duplicates(keep='first', ignore_index=True)
				if self.data.shape[0] > self.max_size:
					self.data = self.data.head[:self.max_size]
				self.min_fitness = self.data['fitness'].min()
			else:
				return
	
	def print(self, max_row=None, max_col=None):
		with pd.option_context('display.max_rows', max_row, 'display.max_columns', max_col):  # more options can be specified also
			print(self.data)
	
	def delta_aic(self):
		if "delta_aic_best" in self.data.columns:
			return
		else:
			self.data["aic"] = self.data["aic"]*-1 #reverse the neg sign i added for maximizing
			best=self.data["aic"].min()
			self.data["delta_aic_best"] = self.data["aic"]-best
	
	def akaike_weights(self):
		if "delta_aic_best" not in self.data.columns:
			self.delta_aic()
		#weight(i) = e^(-1/2 delta_aic_best) / sum(e^(-1/2 delta_aic_best(k)))
		#where the denominator is summed over k models
		#delta_aic = self.data["delta_aic_best"].to_numpy()
		#sum_k = self.data["delta_aic_best"].to_numpy()
		#did a test agains MuMIn::Weights in R and this seems to be working
		self.data["akaike_weight"]=((np.exp(-0.5*self.data["delta_aic_best"])) / (sum(np.exp(-0.5*self.data["delta_aic_best"]))))
	
	def cumulative_akaike(self, threshold=None):
		if "akaike_weight" not in self.data.columns:
			self.akaike_weights()
		self.data["acc_akaike_weight"] = self.data["akaike_weight"].cumsum()
		if threshold: 
			if threshold >0 and threshold < 1:
				pass

	def importance_of_terms(self):
		pass
	
	def output(self):
		#Make sure to remove weights/ shapes where variable isn't selected
		#get absolute value of AIC (made them negative so all metrics could use maximize function)
		pass

