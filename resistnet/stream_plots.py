import os
import sys
import momepy
import pandas as pd
import geopandas as gpd
import numpy as np
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt

def plotEdgesToStreams(network, res, oname):
	G=nx.Graph(nx.read_gpickle(network).to_undirected())
	nodeDF, edgeDF=momepy.nx_to_gdf(G)
	geoDF = edgeDF.merge(res, on="EDGE_ID")
	sns.set(style="ticks")
	geoDF.plot(column="Resistance", cmap = "RdYlGn_r", legend=True)
	plt.title("Stream network colored by resistance")
	plt.savefig((str(oname)+".streamsByResistance.pdf"))
	plt.clf()
	plt.close()
