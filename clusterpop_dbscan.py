import os
import sys
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from geopy.distance import great_circle
from shapely.geometry import MultiPoint
from sortedcontainers import SortedDict

#credit to Geoff Boeing at https://geoffboeing.com/2014/08/clustering-to-reduce-spatial-data-set-size/ for the tutorial

#here coords should be a dictionary with key=ID and value=tuple of lat,long
#epsilon=max distance (in km) between points in a cluster
def dbscan_cluster(coords, epsilon, min_samples):
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

#function to convert SortedDict of coordinates to a numpy matrix 
def coordsToMatrix(coords):
	return(pd.DataFrame([[coords[k][0], coords[k][1]] for k in coords], columns=["long", "lat"]).to_numpy())