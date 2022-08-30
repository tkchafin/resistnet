import os
import sys
import momepy
import pandas as pd
import geopandas as gpd
import numpy as np
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt

def plotEdgesToStreams(network, res, oname, id_col="EDGE_ID"):
	G=nx.Graph(nx.read_gpickle(network).to_undirected())
	nodeDF, edgeDF=momepy.nx_to_gdf(G)
	edgeDF.EDGE_ID = edgeDF[id_col]
	edgeDF.EDGE_ID.astype(int)
	geoDF = edgeDF.merge(res, on="EDGE_ID")
	sns.set(style="ticks")
	geoDF.plot(column="Resistance", cmap = "RdYlGn_r", legend=True)
	plt.title("Stream network colored by resistance")
	plt.savefig((str(oname)+".streamsByResistance.pdf"))
	plt.clf()
	plt.close()
