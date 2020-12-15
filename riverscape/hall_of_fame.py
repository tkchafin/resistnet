import pandas as pd 

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
		cols.append("deltaAIC")
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
			print('hall of fame not full')
			select_size=space
			if space> popDF.shape[0]:
				select_size=popDF.shape[0]
			print(select_size)
			self.data.append(popDF.head(select_size))
			print(self.data)
			self.min_fitness = self.data['fitness'].min()
		else:
			if popDF['fitness'].max() > self.min_fitness:
				subset=popDF[popDF.fitness > self.min_fitness]
				popDF.append(subset)
				self.data = self.data.sort_values('fitness', ascending=False)
				self.data = self.data.drop_duplicates(keep='first', ignore_index=True)
				if self.data.shape[0] > self.max_size:
					self.data = self.data.head(self.max_size)
				self.min_fitness = self.data['fitness'].min()
			else:
				return
	
	def print(self, max_row=None, max_col=None):
		with pd.option_context('display.max_rows', max_row, 'display.max_columns', max_col):  # more options can be specified also
			print(self.data)
	
	def akaike_weights():
		pass

	def importance_of_terms():
		pass

