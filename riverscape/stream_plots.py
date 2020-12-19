import os
import sys

import pandas as pd
import geopandas as gpd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plotEdgesToStreams(graph, res, shp_base, oname):
	geoDF = gpd.read_file(shp_base)
	geoDF = geoDF.merge(res, on="EDGE_ID")
	sns.set(style="ticks")
	geoDF.plot(column="Resistance", cmap = "RdYlGn_r", legend=True)
	plt.title("Stream network colored by resistance")
	plt.savefig((str(oname)+".streamsByResistance.pdf"))
	plt.clf()