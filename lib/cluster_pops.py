import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.cluster import DBSCAN
from geopy.distance import great_circle
from shapely.geometry import MultiPoint
from sortedcontainers import SortedDict

import matplotlib.pyplot as plt

#credit to Geoff Boeing at https://geoffboeing.com/2014/08/clustering-to-reduce-spatial-data-set-size/ for the tutorial

#here coords should be a dictionary with key=ID and value=tuple of lat,long
#epsilon=max distance (in km) between points in a cluster
def dbscan_cluster(coords, epsilon, min_samples, out=None):
	kms_per_radian = 6371.0088

	points=coordsToMatrix(coords)
	
	#run DBSCAN
	eps=epsilon/kms_per_radian
	db=DBSCAN(eps=eps, min_samples=min_samples, algorithm='ball_tree', metric='haversine').fit(np.radians(points))
	
	#save cluster labels
	cluster_labels=db.labels_
	
	#build SortedDict to return
	popmap=SortedDict()
	i=0
	for k in coords.keys():
		pop="DB_"+str(cluster_labels[i])
		if pop not in popmap:
			l = [k]
			popmap[pop] = l
		else:
			popmap[pop].append(k)
		i+=1
	return(popmap)

#function to convert SortedDict of coordinates to a data frame
def coordsToMatrix(coords):
	return(pd.DataFrame([[coords[k][0], coords[k][1]] for k in coords], columns=["long", "lat"]).to_numpy())

def coordsToDataFrame(coords):
	return(pd.DataFrame([[coords[k][0], coords[k][1]] for k in coords], columns=["long", "lat"]))
	
#function to find the centroid of a set of points
#requires a SortedDict of coordinates and a SortedDict giving population IDs
"""Coords:
	key 	value
	SampleName	Tuple(Lat, Long)
	
	popmap:
	PopulationName	list(SampleName,...)
"""
def getClusterCentroid(coords, popmap, out=None):
	centroids=SortedDict()
	ofh=None
	if out:
		ofh=out+".clusterCentroids.txt"
	log=""
	for pop in popmap.keys():
		cluster=getPopCoordsMatrix(coords, popmap[pop])
		if len(cluster)<1:
			print("ERROR: getClusterCentroid(): No coordinates in cluster:",pop)
			sys.exit(1)
		
		#add cluster to logfile (if provided)
		log=log+"Population="+pop+"\n"
		log=log+str(cluster)+"\n"
		
		#get centroid point
		centroid = (MultiPoint(cluster).centroid.x, MultiPoint(cluster).centroid.y)
		log=log+"Centroid="+str(centroid)+"\n"
		centroids[pop] = centroid
	
	#write optional logfile
	if out:
		f=open(ofh, "w")
		f.write(log)
		f.close()
	return(centroids)

#returns a matrix of coordinates from a SortedDict of sample coordinates, given a list to subset
def getPopCoordsMatrix(d, l):
	return(pd.DataFrame([[d[k][0], d[k][1]] for k in d if k in l], columns=["long", "lat"]).to_numpy())

#function plots clustered coordinates given a SortedDict of coords and a population map
def plotClusteredPoints(point_coords, popmap, out, centroids=None):
	#set output file name
	ofh=out+".clusteredPoints.pdf"
	sns.set(style="ticks")
	
	#get 1D popmap
	pmap=flattenPopmap(popmap)
	#print(pmap)
	
	df=pd.DataFrame([[ind, pmap[ind], point_coords[ind][0], point_coords[ind][1]] for ind in point_coords], columns=["sample", "pop", "long", "lat"])
	
	cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)
	ax = sns.scatterplot(x="long", y="lat", hue="pop", palette="Set2",data=df)
	
	#plot centroid positions if available
	if centroids:
		cdf=pd.DataFrame([[p, centroids[p][0], centroids[p][1]] for p in centroids], columns=["pop", "long", "lat"])
		#print(cdf)
		sns.scatterplot(x="long", y="lat", hue="pop", palette="Set2", data=cdf, legend=False, marker="X", ax=ax)
		
	plt.savefig(ofh)
	plt.clf()

#plots a simple histogram
def plotHistogram(dat, out):
	of=str(out) + ".snapDistances.pdf"
	sns.set(style="ticks")
	
	x = pd.Series(dat, name="Snap distance (km)")
	dp = sns.distplot(x, kde=True, rug=True)
	plt.savefig(of)
	plt.clf()

#utility function, converts popmap of form key=pop; value=list(inds) to key=ind; value=pop
def flattenPopmap(popmap):
	new_popmap=dict()
	for k in popmap:
		for i in popmap[k]:
			new_popmap[i]=k
	return(new_popmap)
