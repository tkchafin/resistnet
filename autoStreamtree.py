
import sys
import os
import geopy
import itertools
import traceback
import math
import getopt
import scipy
import pandas as pd
import geopandas as gpd
import numpy as np
import networkx as nx
import Mantel
import seaborn as sns
from scipy import stats
from sortedcontainers import SortedDict
from sklearn.linear_model import LinearRegression
from shapely.geometry import LineString, point, Point
from networkx import NodeNotFound
from networkx import NetworkXNoPath
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import pickle
from math import radians, degrees, sin, cos, asin, acos, sqrt

import genetic_distances as gendist
import cluster_pops as clust
from ast_menu import parseArgs

#from skbio.stats.distance import mantel
#TODO:
#--Jost D option to use global freq rather than average across pops
#--option to log stream distances
#--options for how gaps, Ns, ambigs are treated in p-distance calculation
#--add some sort of outlier detection? Maybe inspited by ABGD algorithm

#TODO: Some parsing to make sure that dist calcs are compatible with data types 
#TODO: Add msat compatibility. I might just convert them to dummy nucleotide values??

#TODO: Option to either die or remove points that can't be reached from rest of network 

def main():

	params = parseArgs()

	print("Starting...\n")

	if params.network:
		print("Reading network from saved file: ", params.network)
		G=nx.OrderedGraph(nx.read_gpickle(params.network).to_undirected())
	else:
		print("Building network from shapefile:",params.shapefile)
		print("WARNING: This can take a while with very large files! If taking too long, try clipping your shapefile to a smaller area or removing unneeded streams")
		G=nx.OrderedGraph(nx.read_shp(params.shapefile, simplify=False, strict=True).to_undirected())

	#if reading populations by 2nd column
	popmap = SortedDict()
			
	#parse dataset
	points = pd.read_csv(params.geodb, sep="\t", header="infer")
	point_coords=SortedDict()
	point_labels=dict()
	snapDists=dict()
	seqs = list()
	verb=True
	first=True
	for idx, row in points.iterrows():
		name = None
		data = None
		if params.run == "GENDIST":
			name = row[0]
			data = tuple([row[3], row[2]])
		else:
			if not params.pop and not params.clusterpop:
				#print(tuple([row[3], row[2]]))
				#--geopop and individual-level snap coordinates to nodes here
				node = snapToNode(G, tuple([row[3], row[2]]))
				snapDists[row[0]] = great_circle(node[0], node[1], row[3], row[2])
			else:
				#if pop or clusterpop, extract centroid later
				node = tuple([row[3], row[2]])
			#print(node)
			data = node
			name = row[0]
			#point_labels[node]=str(row[0])
		point_coords[name] = data
		seq_data = parseLoci(params, list(row[4:]), verbose=verb)
		verb=False
		if first:
			for i, loc in enumerate(seq_data):
				temp = dict()
				seqs.append(temp)
				seqs[i][name]=loc
			numLoci=len(seq_data)
		else:
			for i, loc in enumerate(seq_data):
				seqs[i][name] = loc
		if params.geopop:
			if point_coords[name] not in popmap:
				l = [name]
				popmap[point_coords[name]] = l
			else:
				popmap[point_coords[name]].append(row[0])
		elif params.pop:
			if row[1] not in popmap:
				l = [name]
				popmap[row[1]] = l
			else:
				popmap[row[1]].append(name)
		first=False

	print("Found",len(points.columns)-4,"loci.\n")
	#points["node"]=point_coords

	print("Read",str(len(point_coords.keys())),"individuals.")
	#print(list(point_coords.keys()))
	print()
	
	if params.pop or params.geopop:
		print("Read",str(len(popmap.keys())),"populations.")
		#print(list(popmap.keys()))
		print()
	
	"""
	For population-level analyses, generate population maps and centroids here 
	according to user-input options: --pop, --geopop, --clusterpop
	"""
	#get population centroid
	if params.pop or params.geopop or params.clusterpop:
		if params.clusterpop:
			#create population clusters using DBSCAN
			print("Running DBSCAN clustering with min_samples=",params.min_samples,"and epsilon=",params.epsilon)
			popmap=clust.dbscan_cluster(point_coords, params.epsilon, params.min_samples)
			num_clusters=len(popmap.keys())
			print("Found",str(num_clusters),"clusters!")
			print(popmap)
			print("\n")

			#calculate centroids for clusters
			pop_temp=clust.getClusterCentroid(point_coords, popmap, params.out)
			
			#now, snap pop_coords to nodes
			pop_coords = SortedDict()
			for p in pop_temp:
				node = snapToNode(G, pop_temp[p])
				snapDists[p] = great_circle(node[0], node[1], pop_temp[p][0], pop_temp[p][1])
				pop_coords[p]=node
		elif params.pop or params.geopop:
			#popmap generated earlier when parsing input file!
			#still need to calculate centroids:
			print("Calculating population centroids...")
			pop_temp=clust.getClusterCentroid(point_coords, popmap, params.out)
			#note in the case of --geopop the centroid is the joint snapped-to location
			
			#now, snap pop_coords to nodes
			pop_coords = SortedDict()
			if params.geopop:
				pop_coords = pop_temp
			else:
				for p in pop_temp:
					node = snapToNode(G, pop_temp[p])
					snapDists[p] = great_circle(node[0], node[1], pop_temp[p][0], pop_temp[p][1])
					pop_coords[p]=node
		#write popmap to file 
		flat = clust.flattenPopmap(popmap)
		temp = pd.DataFrame(list(flat.items()), columns=['IND_ID', 'POP_ID'])
		temp.to_csv((str(params.out) + ".popmap.txt"), sep="\t")
		del flat
		del temp
		
		#plot grouped samples
		#TODO: for --geopop maybe plot original coordinates with "snap" as centroid here??
		clust.plotClusteredPoints(point_coords, popmap, params.out, pop_coords)

		
	#plot histogram of snap distances
	clust.plotHistogram(list(snapDists.values()), params.out)
	dtemp = pd.DataFrame(list(snapDists.items()), columns=['name', 'km'])
	dtout = str(params.out) + ".snapDistances.txt"
	dtemp.to_csv(dtout, sep="\t")
	del dtemp
	del dtout
	del snapDists
		
		
	#traverse graph to fill: streamdists, gendists, and incidence matrix
	#calculate genetic distance matrix -- right now that is a JC69-corrected Hamming distance
	#D
	if params.run != "STREAMDIST":
		gen=None
		print("\nCalculating genetic distances...")
		if not params.genmat:
			if params.dist in ["PDIST", "TN84", "TN93", "K2P", "JC69"]:
				gen = gendist.getGenMat(params.dist, point_coords, seqs, ploidy=params.ploidy, het=params.het, loc_agg=params.loc_agg)
				print("Genetic distances:")
				np.set_printoptions(precision=3)
				print(gen, "\n")
				
				#write individual genetic distances to file
				ind_genDF = pd.DataFrame(gen, columns=list(point_coords.keys()), index=list(point_coords.keys()))
				ind_genDF.to_csv((str(params.out) + ".indGenDistMat.txt"), sep="\t", index=True)
				del ind_genDF
				
				if params.pop or params.geopop or params.clusterpop:
					print("Aggregating pairwise population genetic distances from individual distances using:",params.pop_agg)
			else:
				if not params.pop and not params.geopop:
					print("ERROR: Distance metric",params.dist,"not possible without population data.")
					sys.exit(1)
			#calculate population gendistmat
			if params.pop or params.geopop or params.clusterpop:
				pop_gen = gendist.getPopGenMat(params.dist, gen, popmap, point_coords, seqs, pop_agg=params.pop_agg, loc_agg=params.loc_agg, ploidy=params.ploidy, global_het=params.global_het)
				print("Population genetic distances:")
				np.set_printoptions(precision=3)
				print(pop_gen, "\n")
				
				#write population genetic distances to file
				#print(list(pop_coords.keys()))
				pop_genDF = pd.DataFrame(pop_gen, columns=list(pop_coords.keys()), index=list(pop_coords.keys()))
				pop_genDF.to_csv((str(params.out) + ".popGenDistMat.txt"), sep="\t", index=True)
				del pop_genDF
		else:
			print("\nReading genetic distances from provided matrix:", params.genmat)
			inmat = pd.read_csv(params.genmat, header=0, index_col=0, sep="\t")
			if (set(inmat.columns) != set(inmat.index.values)):
				print("Oh no! Input matrix columns and/ or rows don't appear to be labelled. Please provide an input matrix with column and row names!")
				sys.exit(1)
			else:
				agg=False
				#first check if it fits whatever the user input was (i.e. --pop)
				if params.pop:
					if len(inmat.columns) != len(popmap.keys()):
						print("Found",str(len(inmat.columns)), "columns in provided matrix. This doesn't match number of populations from popmap.")
						if (len(inmat.columns)) != len(point_coords):
							print("Doesn't match number of individuals either! Please check your matrix format.")
							sys.exit(1)
						else:
							print("Assuming input matrix has individual distances... Aggregating using the following method (--pop_agg):", str(params.pop_agg))
							agg=True
					else:
						#re-order using pop orders
						inmat = inmat.reindex(list(popmap.keys()))
						inmat = inmat[list(popmap.keys())]
						pop_gen = inmat.to_numpy()
						del(inmat)
				elif params.geopop or params.clusterpop:
					if (len(inmat.columns)) != len(point_coords):
						print("Found",str(len(inmat.columns)), "columns in provided matrix. This doesn't match number of individuals.")
						print("When using --geopop or --clusterpop, the provided matrix must represent individual-level distances.")
						sys.exit(1)
					else:
						#re-order using pop orders
						inmat = inmat.reindex(list(point_coords.keys()))
						inmat = inmat[list(point_coords.keys())]
						gen = inmat.to_numpy()
						agg = True
						del(inmat)
				else:
					if (len(inmat.columns)) != len(point_coords):
						print("Found",str(len(inmat.columns)), "columns in provided matrix. This doesn't match number of individuals.")
						sys.exit(1)
					else:
						#re-order using pop orders
						inmat = inmat.reindex(list(point_coords.keys()))
						inmat = inmat[list(point_coords.keys())]
						gen = inmat.to_numpy()
						del(inmat)
				#if --geopop or --clusterpop, it should be an ind matrix
				#if so, need to aggregate according to --pop_agg
				#print(pop_gen)
				if agg:
					print("Aggregating user-provided individual-level distance matrix using:",params.pop_agg)
					pop_gen = gendist.getPopGenMat("PDIST", gen, popmap, point_coords, seqs, pop_agg=params.pop_agg, loc_agg=params.loc_agg, ploidy=params.ploidy, global_het=params.global_het)
		
		if params.run == "GENDIST":
			print("Genetic distances:")
			np.set_printoptions(precision=3)
			print(gen, "\n")
			
			#write individual genetic distances to file
			ind_genDF = pd.DataFrame(gen, columns=list(point_coords.keys()), index=list(point_coords.keys()))
			ind_genDF.to_csv((str(params.out) + ".indGenDistMat.txt"), sep="\t", index=True)
			
			if params.pop or params.geopop or params.clusterpop:
				pop_gen = gendist.getPopGenMat(params.dist, gen, popmap, point_coords, seqs, pop_agg=params.pop_agg, loc_agg=params.loc_agg, ploidy=params.ploidy, global_het=params.global_het)
				print("Population genetic distances:")
				np.set_printoptions(precision=3)
				print(pop_gen, "\n")
				
				#write population genetic distances to file
				#print(list(pop_coords.keys()))
				pop_genDF = pd.DataFrame(pop_gen, columns=list(pop_coords.keys()), index=list(pop_coords.keys()))
				pop_genDF.to_csv((str(params.out) + ".popGenDistMat.txt"), sep="\t", index=True)
				del pop_genDF
			sys.exit(0)

	#for ia,ib in itertools.combinations(range(0,len(popmap)),2):
	#	print(popmap.keys()[ia])
	#	print(popmap.keys()[ib])

	#EXTRACT SUBGRAPH
	if params.run != "GENDIST":
		if params.pop or params.geopop or params.clusterpop:
			points=pop_coords
		else:
			points=point_coords
		
		#output points to table
		p = getPointTable(points)
		p.to_csv((str(params.out)+".pointCoords.txt"), sep="\t")
		del p
		
		#first pass grabs subgraph from master shapefile graph
		print("\nExtracting full subgraph...")
		ktemp=pathSubgraph(G, points, extractFullSubgraph, params.reachid_col, params.length_col)
		del G

		#second pass to simplify subgraph and collapse redundant nodes
		print("\nMerging redundant paths...\n")
		K=pathSubgraph(ktemp, points, extractMinimalSubgraph, params.reachid_col, params.length_col)
		del ktemp
		
		#grab real coordinates as node positions for plotting
		pos=dict()
		for n in K.nodes:
			pos[n] = n
		#print(pos)
		
		#make a color map to color sample points and junctions differently 
		color_map = []
		for node in K:
			if node in point_coords.values():
				color_map.append("blue")
			else:
				color_map.append("black")
		#draw networkx 
		nx.draw_networkx(K, pos, with_labels=False, node_color=color_map, node_size=50)
		
		#get LENGTH_KM attributes for labelling edges
		edge_labels = nx.get_edge_attributes(K,params.length_col)
		for e in edge_labels:
			edge_labels[e] = "{:.2f}".format(edge_labels[e])
		
		nx.draw_networkx_edge_labels(K, pos, edge_labels=edge_labels, font_size=6)
		
		#save minimized network to file (unless we already read from one)
		if not params.network:
			net_out=str(params.out) + ".network"
			nx.write_gpickle(K, net_out, pickle.HIGHEST_PROTOCOL)
		elif params.overwrite:
			net_out=str(params.out) + ".network"
			nx.write_gpickle(K, net_out, pickle.HIGHEST_PROTOCOL)
		else:
			print("NOTE: Not over-writing existing network. To change this, use --overwrite")

		network_plot=str(params.out) + ".subGraph.pdf"
		plt.savefig(network_plot)
	
	
	#sys.exit()
	
	#calculate pairwise observed stream distances and indence matrix
	#calculate incidence matrix X, which takes the form:
	#nrows = rows in column vector form of D
	#ncols = number of collapsed branches in stream network K
	if params.run in ["STREAMDIST", "DISTANCES", "STREAMTREE", "IBD", "ALL"]:
		if params.pop or params.geopop or params.clusterpop:
			points=pop_coords
			gen=pop_gen
		else:
			points=point_coords

		#calculate stream distances and incidence matrix
		(sdist, inc) = getStreamMats(points, K, params.length_col)
		print("\nStream distances:")
		print(sdist)
		sDF = pd.DataFrame(sdist, columns=list(points.keys()), index=list(points.keys()))
		sDF.to_csv((str(params.out) + ".streamDistMat.txt"), sep="\t", index=True)
		del sDF
		if params.run == "STREAMDIST":
			sys.exit(0)
		
		
		#if user requested testing isolation-by-distance:
		if params.run in ['ALL', 'IBD']:
			#HERE: Implement the IBD calculations and plots
			print("\nTesting for isolation-by-distance using Mantel test with",params.permutations,"permutations")
			testIBD(gen, sdist, params.out, params.permutations)
			
			#plot geographic x genetic DISTANCES
			gn=get_lower_tri(gen)
			go=get_lower_tri(sdist)
			sns.jointplot(gn, go, kind="reg", stat_func=r2, x="Genetic distance", y="Geographic distance (km)")
			plt.savefig((str(params.out)+".isolationByDistance.pdf"))
			# sns.jointplot(gn, np.log(go), kind="reg", stat_func=r2)
			# plt.savefig((str(params.out)+".genXlngeo.pdf"))
			
			
	
	if params.run in ["STREAMTREE", "ALL"]:
		if params.pop or params.geopop or params.clusterpop:
			gen=pop_gen
		print("\nIncidence matrix:")
		print(inc)
		ofh=params.out+".incidenceMatrix.txt"
		with np.printoptions(precision=0, suppress=True):
			np.savetxt(ofh, inc, delimiter="\t")
		print("Incidence matrix dimensions:")
		print(inc.shape)

		#fit least-squares branch lengths
		print()
		R = fitLeastSquaresDistances(gen, inc.astype(int), params.iterative, params.out,params.weight)
		print("\nFitted least-squares distances:")
		print(R)
		
		#check observed versus fitted distances:
		pred=getFittedD(points, gen, inc, R)
		print("\nComparing observed versus predicted genetic distances:")
		print(pred)
		pred.to_csv((str(params.out)+".obsVersusFittedD.txt"), sep="\t")
		sns.jointplot("observed_D", "predicted_D", data=pred, kind="reg", stat_func=r2)
		plt.savefig((str(params.out)+".obsVersusFittedD.pdf"))
		del pred
		#plt.show()
		
	
		#Now, annotate originate geoDF with dissolved reach IDs
		#also, need to collect up the stream tree fitted D to each dissolved reach
		#finally, could add residuals of fitting D vs LENGTH_KM?
		#maybe include logDxlength, DxlogLength, logDxlogLength as well?
		
		#get list of all REACHIDs to extract from geoDF
		edge_data = nx.get_edge_attributes(K,params.reachid_col)
		reach_to_edge = dict()
		i=0
		edges = list()
		for e in edge_data:
			for r in edge_data[e]:
				reach_to_edge[r] = str(i)
			edges.append(i)
			i+=1
		del edge_data
		
		#save reach_to_edge table to file
		r2eDF = pd.DataFrame(list(reach_to_edge.items()), columns=[params.reachid_col,'EDGE_ID'])
		r2eDF.to_csv((str(params.out)+".reachToEdgeTable.txt"), sep="\t")
		
		#read in original shapefile as geoDF and subset it
		print("\nExtracting attributes from original dataframe...")
		geoDF = gpd.read_file(params.shapefile)
		mask = geoDF[params.reachid_col].isin(list(reach_to_edge.keys()))
		temp = geoDF.loc[mask]
		del mask
		del reach_to_edge
		
		#join EDGE_ID to geoDF
		geoDF = temp.merge(r2eDF, on=params.reachid_col)
		del temp
		del r2eDF
		
		#plot by edge ID
		base = geoDF.plot(column="EDGE_ID", cmap = "prism")
		coords = clust.coordsToDataFrame(points)
		geo_coords = gpd.GeoDataFrame(coords, geometry=gpd.points_from_xy(coords.long, coords.lat))
		geo_coords.plot(ax=base, marker='o', color='black', markersize=10, zorder=10)
		
		plt.title("Stream network colored by EDGE_ID")
		plt.savefig((str(params.out)+".streamsByEdgeID.pdf"))
	
		#add in fitted distances & plot
		print("\nPlotting StreamTree fitted distances and writing new shapefile...")
		fittedD = pd.DataFrame({'EDGE_ID':list(edges), 'fittedD':R})
		geoDF['EDGE_ID'] = geoDF['EDGE_ID'].astype(int)
		geoDF = geoDF.merge(fittedD, on='EDGE_ID')
		geoDF.plot(column="fittedD", cmap = "RdYlGn_r", legend=True)
		plt.title("Stream network colored by StreamTree fitted distances")
		plt.savefig((str(params.out)+".streamsByFittedD.pdf"))
	
		#output a final annotated stream network layer
		geoDF.to_csv((str(params.out)+".streamTree.txt"), sep="\t")
		geoDF.to_file((str(params.out)+".streamTree.shp"))

	print("\nDone!\n")

#returns a pandas DataFrame from points dictionary
def getPointTable(points):
	temp=list()
	for p in points:
		temp.append([p, points[p][1], points[p][0]])
	p = pd.DataFrame(temp, columns=['sample', 'lat', 'long'])
	return(p)

#returns r2
def r2(x, y):
	return (stats.pearsonr(x, y)[0] ** 2)

#function takes fitted streamtree distances and calculates predicted genetic distacnes
def getFittedD(points, genmat, inc, r):
	rows=list()
	names=list(points.keys())
	inc_row=0 #tracks what ROW we are in the incidence matrix (streams are columns!)
	for ia, ib in itertools.combinations(range(0,len(points)),2):
		obs=genmat[ia,ib]
		inc_streams=inc[inc_row,]
		pred_dist=np.sum(r[inc_streams==1])
		inc_row+=1
		rows.append([names[ia], names[ib], obs, pred_dist, np.abs(obs-pred_dist)])
	D=pd.DataFrame(rows, columns=['from','to','observed_D', 'predicted_D', 'abs_diff'])
	return(D)
		

#function to calculate great circle distances
#returns in units of KILOMETERS
def great_circle(lon1, lat1, lon2, lat2):
	lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
	return 6371 * (
		acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
	)

def get_lower_tri(mat):
	n=mat.shape[0]
	i=np.tril_indices(n,-1)
	return(mat[i])

#function computes Mantel test using various transformations
def testIBD(gen, geo, out, perms):
	#get flattened lower triangle of each matrix
	gen = get_lower_tri(gen)
	geo = get_lower_tri(geo)

	#non-log pearson
	res=list(Mantel.test(gen, geo, perms=int(perms), method='pearson'))
	rows=list()
	rows.append(['genXgeo','pearson', str(perms), res[0], res[1], res[2]])

	#non-log spearman
	res=list(Mantel.test(gen, geo, perms=int(perms), method='spearman'))
	rows.append(['genXgeo','spearman', str(perms), res[0], res[1], res[2]])

	#print(rows)
	ibd=pd.DataFrame(rows,  columns=['test', 'method', 'perms', 'r', 'p' ,'z'])
	print("Mantel test results:")
	print(ibd)
	ibd.to_csv((str(out) + ".isolationByDistance.txt"))
	print()

	

#only necessary for later
#eventually will add capacity to handle phased loci and msats
#will be easier to have this in a separate function
def parseLoci(opts, data, verbose=False):
	if opts.snps:
		if "/" in data[0] and len(data[0]) > 3:
			print("ERROR: Data appear to be phased haplotypes and are incompatible with --snp option")
			sys.exit(1)
		elif "/" not in data[0] and len(data) == 1:
			if verbose:
				print("Data appears to consist of unphased concatenated SNPs...")
			return([gendist.phaseSnp(x.replace(" ","").lower()) for x in data[0]])
		elif "/" in data[0] and len(data) > 1:
			if verbose:
				print("Data appears to consist of phased un-concatenated SNPs...")
			return([str(x).replace(" ","").lower() for x in data])
		elif "/" not in data[0] and len(data) > 1:
			if verbose:
				print("Data appears to consist of unphased un-concatenated SNPs...")
			return([gendist.phaseSnp(str(x).replace(" ","").lower()) for x in data])
		else:
			print("ERROR: Unable to parse SNP input. Please check input file")
			sys.exit(1)
	else:
		#print(data)
		if "/" in data[0] and len(data[0]) > 3:
			if verbose:
				print("Data appear to consist of phased sequence data...")
			if opts.ploidy == 1:
				print("ERROR: Diplotypes were provided for haplotype data!")
				sys.exit(1)
		elif "/" in data[0] and len(data[0]) <=3:
			print("Data appear to consist of phased SNP data... Are you sure you don't need the --snp option?")
			sys.exit(1)
		else:
			if verbose:
				print("Data appear to consist of unphased sequence data... Note that autoStreamtree is unable to infer haplotypes!")
		return([str(x).replace(" ","").lower() for x in data])


#function to compute least-squares branch lengths from a vector of genetic distances D and incidence matrix X
#when iterative = True, negative distances are constrained to 0 and then recomputed
def fitLeastSquaresDistances(D, X, iterative, out, weight=None):
	num_segments = (np.size(X,1))
	#print(num_segments)
	ls = np.zeros(num_segments)
	d = vectorizeMat(D)
	
	#calculate weights matrix and write to file
	W=generateWeightsMatrix(d, weight)
	print("Weights matrix:")
	print(W)
	#ofh=out+".weightsMatrix.txt"
	#np.savetxt(ofh, W, delimiter="\t")
	
	#weighted least-squares optimization
	ls = np.matmul(np.linalg.inv(np.matmul(np.matmul(X.transpose(),W),X)), np.matmul(np.matmul(X.transpose(), W),d))
	print("\nLeast-squared optimized distances:")
	print(ls)
	#ls_ord = np.matmul(np.linalg.inv(np.matmul(X.transpose(),X)), np.matmul(X.transpose(),d))
	#print(ls_ord)
	
	#if using iterative approach
	if iterative:
		ls_old=ls
		if(np.count_nonzero(ls<0.0) > 0):
			print("\nLS-optimized distances contain negative values: Using iterative approach to re-calculate...")
		constrains = list() #save indices of all constrained values
		
		#if negative distances, use iterative procedure to re-calculate
		while (np.count_nonzero(ls<0.0) > 0):
			bad_ind = np.argmin(ls)
			constrains.append(bad_ind)
			#constrain to 0 by removing from incidence matrix
			X = np.delete(X, bad_ind, 1)
			#re-compute values
			ls = np.matmul(np.linalg.inv(np.matmul(np.matmul(X.transpose(),W),X)), np.matmul(np.matmul(X.transpose(), W),d))
		for i in reversed(constrains):
			ls=np.insert(ls, i, 0.0)
		#print(ls)
		
		#write original and constrained results to log file
		ofh=out+".leastSquaresConstrained.txt"
		df=pd.DataFrame({'LS.original':ls_old, 'LS.constrained':ls})
		df.to_csv(ofh, sep="\t")
		
		return(ls)
	else:
		return(ls)

#function generates weights matrix for least-squares method, where weights are on diagonals
def generateWeightsMatrix(d,weight):
	W=np.zeros((len(d), len(d)), dtype=float)
	row,col=np.diag_indices(W.shape[0])
	if weight.upper()=="CSE67":
		W[row,col] = np.ones(len(d))
	elif weight.upper()=="BEYER74":
		if(np.count_nonzero(d==0) > 0):
			print("WARNING: Divide-by-zero in weighted least-squares (weight=1/D).")
		W[row,col] = np.divide(1.0, d, out=np.zeros_like(d), where=d!=0)
	elif weight.upper()=="FM67":
		if(np.count_nonzero(d==0) > 0):
			print("WARNING: Divide-by-zero in weighted least-squares (weight=1/D^2).")
		W[row,col] = np.divide(1.0, np.square(d), out=np.zeros_like(d), where=d!=0)
	else:
		print("ERROR: Weight option",weight,"not recognized. Using ordinary least-squares instead.")
		W[row,col] = np.ones(len(d))
	return(W)
	

#function to convert a pairwise matrix to a 1D vector
def vectorizeMat(mat):
	size = nCr(np.size(mat,0), 2)
	vec = np.zeros((size))
	index = 0
	for ia, ib in itertools.combinations(range(0,np.size(mat,0)),2):
		vec[index] = mat[ia, ib]
		index = index+1
	#print(vec)
	return(vec)

#computes pairwise stream distances and 0/1 incidence matrix for StreamTree calculations
def getStreamMats(points, graph, len_col):
	#make matrix
	dist = np.zeros((len(points),len(points)))
	inc = np.zeros((nCr(len(points),2),len(graph.edges())),dtype=int)
	#establish as nan
	dist[:] = np.nan

	#function to calculate weights for Dijkstra's shortest path algorithm
	#i just invert the distance, so the shortest distance segments are favored
	def dijkstra_weight(attributes):
		return(attributes[len_col]*-1)

	#for each combination, get shortest path and sum the lengths
	index=0
	#print(points)
	for ia, ib in itertools.combinations(range(0,len(points)),2):
		path = nx.bidirectional_dijkstra(graph, points.values()[ia], points.values()[ib], weight=dijkstra_weight)
		if path:
			dist[ia,ib] = float(sum(path_edge_attributes(graph, path[1], len_col)))
			dist[ib,ia] = dist[ia,ib]
		#incidence matrix
		#for each edge in graph, assign 0 if not in path; 1 if in path
		#print("path:",path)
		
		for ie, edge in enumerate(graph.edges()):
			if find_pair(path[1], edge[0], edge[1]):
				#print("yes:",edge)
				inc[index, ie] = 1
			else:
				#print("no")
				inc[index, ie] = 0
		index = index+1
		#print("\n---\n")
	np.fill_diagonal(dist, 0.0)
	return((dist, inc))

#utility function to test if two elements are consecutive in list (irrespective of order)
def find_pair(list, x, y):
	if x not in list or y not in list:
		return(False)
	elif abs(list.index(x)-list.index(y)) == 1:
		return(True)
	else:
		return(False)

#utility function to calculate number of combinations n choose k
def nCr(n,k):
	f = math.factorial
	return f(n) // f(k) // f(n-k)

def path_edge_attributes(graph, path, attribute):
	return [graph[u][v][attribute] for (u,v) in zip(path,path[1:])]

#find and extract paths between points from a graph
def pathSubgraph(graph, nodes, method, id_col, len_col):
	k=nx.OrderedGraph()
	for p1, p2 in itertools.combinations(nodes.values(),2):
		try:
			#print(p1)
			#print(p2)
			
			#function to calculate weights for Dijkstra's shortest path algorithm
			#i just invert the distance, so the shortest distance segments are favored
			def dijkstra_weight(attributes):
				return(attributes[len_col]*-1)
			
			#find shortest path between the two points
			path=nx.bidirectional_dijkstra(graph, p1, p2, weight=dijkstra_weight)
			#print("path:",path)
			#traverse the nodes in the path to build a minimal set of edges
			method(k, graph, nodes.values(), id_col ,len_col, path[1])
			#calculate stream distance
			#stream_dist = sum(path_edge_attributes(graph, path[1], "LENGTH_KM")) #total length of all edges in path
			#calculate corrected genetic distance
			####later
			if p1 not in k:
				k.add_node(p1)
			if p2 not in k:
				k.add_node(p2)
		except NodeNotFound as e:
			print("Node not found:",e)
		except Exception as e:
			traceback.print_exc()
			print("Something unexpected happened:",e)
			sys.exit(1)
	return(k)

#extracts full subgraph from nodes
def extractFullSubgraph(subgraph, graph, nodelist, id_col, len_col, path):
	for first, second in zip(path, path[1:]):
		if first not in subgraph:
			subgraph.add_node(first)
		if second not in subgraph:
			subgraph.add_node(second)

		dat=graph.get_edge_data(first, second)
		subgraph.add_edge(first, second, **dat)


#extracts a simplified subgraph from paths
#only keeping terminal and junction nodes
def extractMinimalSubgraph(subgraph, graph, nodelist, id_col, len_col, path):
	curr_edge = {id_col:list(), len_col:0.0}
	curr_start=None
	#print("Path:",path)
	#print("nodelist:",nodelist)
	#for each pair of nodes in path
	for first, second in zip(path, path[1:]):
		#if first is either: 1) a site node; or 2) a junction node: add to new graph
		#if second is either:a site or junction, add edge and node to new graph
		#if not, keep edge attributes for next edge
		if not curr_start:
			curr_start=first
			if first in nodelist or len(graph[first])>2:
				subgraph.add_node(first)
		#add path attributes to current edge
		dat=graph.get_edge_data(first, second)
		curr_edge[id_col].extend([dat[id_col]] if not isinstance(dat[id_col], list) else dat[id_col])
		curr_edge[len_col]=float(curr_edge[len_col])+float(dat[len_col])
		
		#if second node is a STOP node (=in nodelist or is a junction):
		if second in nodelist or len(graph[second])>2:
			#add node to subgraph
			subgraph.add_node(second) 
			#link current attribute data
			subgraph.add_edge(curr_start, second, **curr_edge)
			#empty edge attributes and set current second to curr_start
			curr_edge = {id_col:list(), len_col:0}
			curr_start = second
		else:
			#otherwise continue building current edge
			continue


#Input: Tuple of [x,y] coordinates
#output: Closest node to those coordinates
def snapToNode(graph, pos):
	#rint("closest_node call:",pos)
	nodes = np.array(graph.nodes())
	node_pos = np.argmin(np.sum((nodes - pos)**2, axis=1))
	#print(nodes)
	#print("closest to ", pos, "is",tuple(nodes[node_pos]))
	return (tuple(nodes[node_pos]))

#Call main function
if __name__ == '__main__':
	main()
