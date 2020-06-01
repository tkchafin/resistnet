
import sys
import os
import geopy
import itertools
import traceback
import math
import getopt
import skbio
import scipy
import pandas as pd
import geopandas as gpd
import numpy as np
import networkx as nx
from sortedcontainers import SortedDict
from skbio import DistanceMatrix
from sklearn.linear_model import LinearRegression
from skbio.stats.distance import mantel
from shapely.geometry import LineString, point, Point
from networkx import NodeNotFound
from networkx import NetworkXNoPath
import matplotlib.pyplot as plt
from geopy.distance import geodesic


#TODO:
#--Jost D option to use global freq rather than average across pops
#--Add CLI and argument parsing
#--finish least-squares fitting
#--option to use Fst, or linearized Fst
#--option to log stream distances
#--options for how gaps, Ns, ambigs are treated in p-distance calculation
#--add some sort of outlier detection? Maybe inspited by ABGD algorithm
#-----How? Delete nodes w/ gen distance edges which are significantly larger than mean?
#-----Would then need to re-fit distances
#--Add method to prevent negative distances fit to segments
#----Streamtree method: Constrain most negative one to zero (remove), then re-fit.
#----Check the Felsenstein book, there might be something in there. I think he has a solution in FITCH
#----Also allow option to input this manually (some may have specialized datasets, e.g. SNPs or msats and might want to use a diff method, or calc Fst)
#--As we add in ecological information, we can consider IBE as well

#TODO: Some parsing to make sure that dist calcs are compatible with data types 
#TODO: Add msat compatibility. I might just convert them to dummy nucleotide values??

#option to apply weights, e.g. Fitch and Margouliash (see Felsenstein 2004 page 152)

def main():

	params = parseArgs()

	print("Starting...\n")
	geoDF = gpd.read_file(params.shapefile)
	#print(geoDF)
	#if params.run != "GENDIST":
	G=nx.read_shp(params.shapefile, simplify=False).to_undirected()
	#print(G.nodes)
	#nx.draw(G, with_labels=True)
	#plt.show()
	#plt.savefig("network.png")
	#sys.exit()

	#if reading populations by 2nd column
	popmap = SortedDict()
			
	#parse dataset
	points = pd.read_csv(params.geodb, sep="\t", header="infer")
	point_coords=SortedDict()
	point_labels=dict()
	numLoci=len(points.columns)-4
	seqs = list()
	verb=True
	for loc in range(0,numLoci):
		temp = dict()
		seqs.append(temp)
	for idx, row in points.iterrows():
		name = None
		data = None
		if params.run == "GENDIST":
			name = row[0]
			data = tuple([row[3], row[2]])
		else:
			#TODO: If --pop need to calculate centroids here? Can replace ind coords with centroid
			
			#print(tuple([row[3], row[2]]))
			node = snapToNode(G, tuple([row[3], row[2]]))
			#print(node)
			data = node
			name = row[0]
			#point_labels[node]=str(row[0])
		point_coords[name] = data
		seq_data = parseLoci(params, list(row[4:]), verbose=verb)
		verb=False
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

	print("Found",len(points.columns)-4,"loci.\n")
	#points["node"]=point_coords

	print("Read individuals in this order:")
	print(list(point_coords.keys()))
	print()
	
	if params.pop or params.geopop:
		print("Read populations in this order:")
		print(list(popmap.keys()))
		print()	


	if params.run != "GENDIST":
		#first pass grabs subgraph from master shapefile graph
		#print(point_coords)
		#print(G.nodes)
		print("Extracting full subgraph...")
		ktemp=pathSubgraph(G, point_coords, extractFullSubgraph)
		del G
		#sys.exit()
		#print(nx.get_node_attributes(ktemp))
		#print(ktemp.nodes)
		#second pass to simplify subgraph and collapse redundant nodes
		print("Merging redundant paths...")
		K=pathSubgraph(ktemp, point_coords, extractMinimalSubgraph)
		del ktemp
		
		#nx.relabel_nodes(K, point_labels, copy=False)
		#print("nodes:")
		#print(K.nodes)
		
		#print("pos dict")
		pos=dict()
		for n in K.nodes:
			pos[n] = n
		#print(pos)
		#pos = nx.spring_layout(K)
		nx.draw(K, pos, with_labels=False)
		edge_labels = nx.get_edge_attributes(K,'LENGTH_KM')
		#print(edge_labels)
		nx.draw_networkx_edge_labels(K, pos, edge_labels=edge_labels, font_size=6)

	network_plot=str(params.out) + ".subGraph.pdf"
	plt.savefig(network_plot)
	#sys.exit()

	#traverse graph to fill: streamdists, gendists, and incidence matrix
	#calculate genetic distance matrix -- right now that is a JC69-corrected Hamming distance
	#D
	gen = None
	pop_gen = None
	if params.dist in ["PDIST", "TN84", "TN93", "K2P", "JC69"]:
		gen = getGenMat(params, point_coords, seqs)
		print("Genetic distances:")
		np.set_printoptions(precision=3)
		print(gen, "\n")
	if params.pop or params.geopop:
		pop_gen = getPopGenMat(params, gen, popmap, point_coords, seqs)
		print("Population genetic distances:")
		np.set_printoptions(precision=3)
		print(pop_gen, "\n")
	if params.run == "GENDIST":
		sys.exit(0)
		
	#calculate pairwise observed stream distances and indence matrix
	#calculate incidence matrix X, which takes the form:
	#nrows = rows in column vector form of D
	#ncols = number of collapsed branches in stream network K
	
	if params.run in ["STREAMDIST", "DISTANCES", "STREAMTREE"]:
		(sdist, inc) = getStreamMats(point_coords, K)
		#print(len(K.edges))
		#sys.exit()
		print("Stream distances:")
		print(sdist)
		
		#TODO: If --pop or --geopop, need to summarize stream dist network here
	
	print("Stream distance dimensions:")
	print(sdist.shape)
	print("Genetic distance dimensions:")
	print(gen.shape)
	
	if params.run in ["STREAMTREE"]:
		print("Incidence matrix:")
		print(inc)
		ofh=params.out+".incidenceMatrix.txt"
		np.savetxt(ofh, inc, sep="\t")
		print("Incidence matrix dimensions:")
		print(inc.shape)

		#fit least-squares branch lengths
		R = fitLeastSquaresDistances(gen, inc, params.iterative, params.out,params.weight)
		print("Fitted least-squares distances:")
		print(R)
		
		sys.exit()

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
			return([phaseSnp(x.replace(" ","").lower()) for x in data[0]])
		elif "/" in data[0] and len(data) > 1:
			if verbose:
				print("Data appears to consist of phased un-concatenated SNPs...")
			return([str(x).replace(" ","").lower() for x in data])
		elif "/" not in data[0] and len(data) > 1:
			if verbose:
				print("Data appears to consist of unphased un-concatenated SNPs...")
			return([phaseSnp(str(x).replace(" ","").lower()) for x in data])
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

def getPopGenMat(opts, indmat, popmap, dat, seqs):
	#make matrix
	genmat = np.zeros((len(popmap),len(popmap)))
	#establish as nan
	genmat[:] = np.nan
	#for each combination, either average ind distances or calc freq dist 
	for ia,ib in itertools.combinations(range(0,len(popmap)),2):
		if opts.dist in ["JC69", "K2P", "PDIST", "TN84", "TN93"]:
			#print(popmap.keys()[ia], popmap.values()[ia])
			inds1 = [dat.index(x) for x in popmap.values()[ia]]
			inds2 = [dat.index(x) for x in popmap.values()[ib]]
			#print(inds1)
			genmat[ia,ib] = (np.mean([indmat[i, j] for i in inds1 for j in inds2]))
			genmat[ib,ia] = genmat[ia,ib]
		elif opts.dist == "JOST" or opts.dist == "LINJOST":
			results=list()
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					continue
				results.append(twoPopJostD(seqs1, seqs2, opts))
			#print(results)
			if len(results) > 1:
				if opts.dist=="LINJOST":
					D = aggregateDist(opts.loc_agg, results)
					genmat[ia,ib] = genmat[ib,ia] = (1.0/(1.0-float(D)))
				else:
					genmat[ia,ib] = genmat[ib,ia] = aggregateDist(opts.loc_agg, results)
			elif len(results) < 1:
				print("ERROR: population",popmap.values()[ia],"or",popmap.values()[ib],"lacks any data")
				sys.exit(1)
			else:
				genmat[ia,ib] = genmat[ib,ia] = results[0]
		elif opts.dist == "GST" or opts.dist == "GSTPRIME":
			HT=list()
			HS=list()
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					continue
				if opts.dist == "GST" or "GSTPRIME":
					(ht, hs) = twoPopHtHs(seqs1, seqs2, opts)
					HT.append(ht)
					HS.append(hs)
			Ht_global = np.mean(HT)
			Hs_global = np.mean(HS)
			if opts.dist == "GST":
				if Ht_global <= 0.0:
					genmat[ia,ib] = genmat[ib,ia] = 0.0
				Gst = ((Ht_global - Hs_global) / Ht_global )
				GprimeST = ((Gst * (1.0 + Hs_global)) / (1.0 - Hs_global))
				genmat[ia,ib] = genmat[ib,ia] = GprimeST
			elif opts.dist == "GSTPRIME":
				Ghedrick = ((2.0*(Ht_global - Hs_global)) / (((2.0*Ht_global) - Hs_global) * (1.0 - Hs_global)))
				genmat[ia,ib] = genmat[ib,ia] = Ghedrick
		elif opts.dist == "FST" or opts.dist == "LINFST":
			num = list() #numerator; a 
			denom = list() #denominator; a*b*c
			for loc in range(0, len(seqs)):
				seqs1 = cleanInds([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = cleanInds([seqs[loc][x] for x in popmap.values()[ib]])
				if not all("/" in x for x in seqs1) or not all("/" in x for x in seqs1):
					print("ERROR: FST estimates require phased data.")
					sys.exit(1)
				if len(seqs1) == 0 or len(seqs2) == 0:
					print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
					continue
				(n, d) = twoPopWeirCockerhamFst(seqs1, seqs2)
				num.append(n)
				denom.append(d)
			if len(num) <= 0 or len(denom) <= 0:
				print("ERROR (twoPopWeirCockerhamFst): No data for pops "+ia+" and "+ib+".")
				sys.exit(1)
			theta = np.sum(num) / np.sum(denom)
			if opts.dist == "FST":
				genmat[ia,ib] = genmat[ib,ia] = theta
			elif opts.dist == "LINFST":
				genmat[ia,ib] = genmat[ib,ia] = (theta / (1-theta))
		elif opts.dist == "NEI83":
			results = list()
			loci = 0.0
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
					continue
				loci += 1.0
				results.append(twoPopNeiDa(seqs1, seqs2))
			Da = (1.0 - (np.sum(results) / loci))
			genmat[ia,ib] = genmat[ib,ia] = Da
		elif opts.dist == "EUCLID":
			results = list()
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
					continue
				results.append(twoPopEuclidDist(seqs1, seqs2))
			euclid = np.sum(results)
			genmat[ia,ib] = genmat[ib,ia] = euclid
		elif opts.dist == "CHORD":
			num = list()
			denom = list()
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
					continue
				(n, d) = twoPopChordDist(seqs1, seqs2)
				num.append(n)
				denom.append(d)
			Dch = np.sqrt((np.sum(num))/(np.sum(denom)))
			genmat[ia,ib] = genmat[ib,ia] = Dch
	np.fill_diagonal(genmat, 0.0)
	if 0.0 in genmat:
		print("WARNING: Coercing negative distances to 0.0")
		genmat[genmat<0.0] = 0.0
	return(genmat)

#computed Nei's 1983 Da estimator 
def twoPopNeiDa(s1, s2):
	s1 = cleanList(s1, ["n", "?", "-", "N"])
	s2 = cleanList(s2, ["n", "?", "-", "N"])
	uniques = uniqAlleles(s1+s2)
	sumSqRt = 0.0
	for allele in uniques:
		if allele in ["-", "?", "n", "N"]:
			continue
		else:
			Xu = float(s1.count(allele) / len(s1))
			Yu = float(s2.count(allele) / len(s2))
			sumSqRt += np.sqrt(Xu*Yu)
	return(sumSqRt)
	
#computes euclidean distance for a single locus
def twoPopEuclidDist(s1, s2):
	uniques = uniqAlleles(s1+s2)
	s1 = cleanList(s1, ["n", "?", "-", "N"])
	s2 = cleanList(s2, ["n", "?", "-", "N"])
	sumSq = 0.0
	for allele in uniques:
		if allele in ["-", "?", "n", "N"]:
			continue
		else:
			Xu = float(s1.count(allele) / len(s1))
			Yu = float(s2.count(allele) / len(s2))
			sumSq += np.square(Xu - Yu)
	return(sumSq)

#Cavalli-Sforza and Edwards 1967 chord distance
#non-nucleotide alleles are deleted
def twoPopChordDist(s1, s2):
	s1 = cleanList(s1, ["n", "?", "-", "N"])
	s2 = cleanList(s2, ["n", "?", "-", "N"])
	uniques = uniqAlleles(s1+s2)
	sumSqRt = 0.0
	for allele in uniques:
		if allele in ["-", "?", "n", "N"]:
			continue
		else:
			Xu = float(s1.count(allele) / len(s1))
			Yu = float(s2.count(allele) / len(s2))
			sumSqRt += np.sqrt(np.square(Xu)*np.square(Yu))
	return((1.0-sumSqRt), (len(uniques)-1.0))

#two population calculator of Weir and Cockerham's THETAst Fst approximation
#returns two values: aC (a-coefficient) and (aC * bC * cC), each summed across alleles for a provided locus
#these are the numerator and denominators that must be summed across all loci
#assumes diploid
def twoPopWeirCockerhamFst(s1, s2):
	num = 0.0
	denom = 0.0
	
	#mean sample size 
	alleles1 = getAlleles(s1) #split alleles s1
	alleles2 = getAlleles(s2) #split alleles s2
	uniques = uniqAlleles(s1+s2) #list of unique alleles only
	r = 2.0 #number of pops
	n1 = float(len(s1)) #pop size of pop 1
	n2 = float(len(s2)) #pop size of pop 2
	csd = np.std([n1, n2])
	cm = np.mean([n1, n2])
	nbar = cm
	csquare = (csd*csd) / (cm*cm)
	nC   = nbar * (1.0 - (csquare/r)) #coeff of pop size variance
	for allele in uniques:
		ac1 = float(alleles1.count(allele))
		ac2 = float(alleles2.count(allele))
		p1 = ac1 / float(len(alleles1))
		p2 = ac2 / float(len(alleles2))
		h1 = getHetFromPhased(allele, s1, count=True)
		h2 = getHetFromPhased(allele, s2, count=True)
		pbar = (ac1+ac2) / (float(len(alleles1)) + float(len(alleles2)))
		ssquare = ((np.sum( [ (n1* (np.square(p1 - pbar)) ), (n2* (np.square(p2 - pbar))) ])) / ((r-1.0)*nbar))
		hbar = ((h1+h2) / (r * nbar))
		a = ((nbar/nC) * 
			(ssquare - 
			((1.0 / (nbar-1.0)) * 
			((pbar * (1.0-pbar)) - 
			((r - 1.0) * ssquare / r) - 
			(hbar / 4.0)))))
		b = ((nbar / (nbar-1.0)) * 
			((pbar * (1.0 - pbar)) - 
			((r - 1.0) * ssquare / r) -  
			(((2.0 * nbar) - 1.0) * hbar / (4.0 * nbar))))
		c = hbar/2.0
		d = a+b+c
		num += a
		denom += d
	return(num, denom)
	
#removes individuals with unknown or gap alleles
def cleanInds(inds):
	ret = list()
	for ind in inds:
		if "-" not in ind and "?" not in ind and "n" not in ind and "N" not in ind:
			ret.append(ind)
	return(ret)

#estimate Jost's D using Nei and Chessers Hs and Ht estimators
def twoPopJostD(seqs1, seqs2, opts):
	if opts.global_het:
		Ht = getGlobalHet(seqs1+seqs2)
	else:
		Ht = getAverageHet(seqs1, seqs2)
	Hs = np.mean([getGlobalHet(seqs1),getGlobalHet(seqs2)])
	harmN = scipy.stats.hmean([(len(seqs1)/opts.ploidy), (len(seqs1)/opts.ploidy)])
	Hs_est = Hs * ((2.0*harmN)/((2.0*harmN)-1.0))
	Ht_est = Ht + (Hs / harmN*2.0*2.0)
	if Ht_est == 0.0:
		return(0.0)
	D = (Ht_est - Hs_est) * 2.0 #b/c in pw estimate, N/N-1 is always 2
	if D == 2:
		return(1.0)
	elif D <= 0.0:
		return(0.0)
	return(D)

#computes Nei's Fst estimator (Gst) using Nei and Chessers Hs and Ht estimators
#also applies Hedrick's (2005) sample size corection, thus returning G'st
def twoPopHtHs(seqs1, seqs2, opts):
	if opts.global_het:
		Ht = getGlobalHet(seqs1+seqs2)
	else:
		Ht = getAverageHet(seqs1, seqs2)
	Hs = np.mean([getGlobalHet(seqs1),getGlobalHet(seqs2)])
	harmN = scipy.stats.hmean([(len(seqs1)/opts.ploidy), (len(seqs1)/opts.ploidy)])
	Hs_est = Hs * ((2.0*harmN)/((2.0*harmN)-1.0))
	Ht_est = Ht + (Hs / (harmN*2.0*2.0))
	# print("Hs_est:",Hs_est)
	# print("Ht_est:",Ht_est)
	# if Ht_est == 0.0:
	# 	return(0.0)
	# Gst = ((Ht_est - Hs_est) / Ht_est )
	# GprimeST = ((Gst * (1.0 + Hs_est)) / (1.0 - Hs_est))
	return((Ht_est, Hs_est))

def aggregateDist(method, stuff):
	if method == "HARM":
		try:
			return(scipy.stats.hmean(stuff))
		except ValueError as e:
			print(e)
			print("ERROR (DivideByZero): Harmonic mean cannot be calculated using a zero distance. Try recomputing using the \"ADJHARM\" option.")
			print("")
			sys.exit(1)
	elif method == "ARITH":
		return(np.mean(stuff))
	elif method == "GEOM":
		return(scipy.stats.mstats.gmean(stuff))
	elif method == "MEDIAN":
		return(np.median(stuff))
	elif method == "MAX":
		return(np.max(stuff))
	elif method == "MIN":
		return(np.min(stuff))
	elif method == "ADJHARM":
		return(adjustedHarmonicMean(stuff))

#computes an harmonic mean corrected for non-positive values
def adjustedHarmonicMean(stuff):
	s=np.array(stuff)
	vals = s[s>0.0]
	bads = s[s<=0.0]
	mu = (1.0 / (np.sum([1.0/x for x in vals]) / (len(vals)-len(bads)))) * ((len(vals)-len(bads))/len(vals))
	return(mu)

#returns observed heterozygosity of an allele given a list of phased 
#genotypes (e.g. allele1/allele2 for each individual)
#assumed diploid
def getHetFromPhased(allele, phasedList, count=False):
	hets = 0.0
	twoN = (len(phasedList)) * 2.0
	for genotype in phasedList:
		if "/" not in genotype:
			print("ERROR (getHetFromPhased): Phased genotypes are required.")
		gens = genotype.split("/")
		if gens[0] == allele and gens[1] != allele:
			hets += 1.0
			continue
		elif gens[1] == allele and gens[0] != allele:
			hets += 1.0
			continue
		else:
			continue
	if count==True:
		return(hets)
	else:
		return(hets/twoN)

#function to compute global expected heterozygosities (Ht) from a set of sequences
#will be returned as a list of Ht estimates per locus
def getGlobalHet(seqs):
	hom = 0.0
	#for each allele at locus
	uniq_alleles = cleanList(set(seqs), ["n", "?", "-", "N"])
	#frequency is num of allele over total size
	freqs = [np.square(float(seqs.count(x)/len(seqs))) for x in uniq_alleles]
	hom = np.sum(freqs)
	return(1.0-hom)

#function to compute mean expected heterozygosities (Ht) from populations
#will be returned as a list of Ht estimates per locus
def getAverageHet(s1, s2):
	hom = 0.0
	#for each allele at locus
	uniq_alleles = cleanList(set(s1+s2), ["n", "?", "-", "N"])
	#frequency is num of allele over total size
	freqs = [np.square(np.mean([float(s1.count(x)/len(s1)), float(s1.count(x)/len(s1))])) for x in uniq_alleles]
	hom = np.sum(freqs)
	return(1.0-hom)

def cleanList(l, bads):
	if not any(item not in bads for item in set(l)):
		return(False)
	for b in bads:
		if b in l:
			l.remove(b)
	return(l)

def getAlleles(s):
	return(sum([x.split("/") for x in s], []))

def uniqAlleles(s):
	return(set(sum([x.split("/") for x in s], [])))

#function to compute least-squares branch lengths from a vector of genetic distances D and incidence matrix X
#when iterative = True, negative distances are constrained to 0 and then recomputed
def fitLeastSquaresDistances(D, X, iterative, out, weight=None):
	num_segments = (np.size(X,1))
	print(num_segments)
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
	print("Least-squared optimized distances:")
	print(ls)
	#ls_ord = np.matmul(np.linalg.inv(np.matmul(X.transpose(),X)), np.matmul(X.transpose(),d))
	#print(ls_ord)
	
	#if using iterative approach
	if iterative:
		ls_old=ls
		if(np.count_nonzero(ls<0.0) > 0):
			print("LS-optimized distances contain negative values: Using iterative approach to re-calculate...")
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
def getStreamMats(points, graph):
	#make matrix
	dist = np.zeros((len(points),len(points)))
	inc = np.zeros((nCr(len(points),2),len(graph.edges())),dtype=int)
	#establish as nan
	dist[:] = np.nan

	#for each combination, get shortest path and sum the lengths
	index=0
	print(points)
	for ia, ib in itertools.combinations(range(0,len(points)),2):
		path = nx.bidirectional_dijkstra(graph, points.values()[ia], points.values()[ib], weight=dijkstra_weight)
		if path:
			dist[ia,ib] = float(sum(path_edge_attributes(graph, path[1], "LENGTH_KM")))
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

#function to parse alignments for base frequencies

#function computes pairwise JC69-corrected genetic distances
def getGenMat(opts, points, seqs):
	#make matrix
	genmat = np.zeros((len(points),len(points)))
	#establish as nan
	genmat[:] = np.nan

	#for models which relax equal nuc frequencies, get global frequencies for each locus
	#freqs will be a list of loci, with each locus as a dist of freqs
	if opts.dist in ["TN84", "TN93"]:
		freqs = getNucFreqs(seqs, opts.ploidy)
		index = 1
		for f in freqs:
			print("Empirical base frequencies for locus",index, end=": [ ")
			for n in f:
				print(f'{n}={f[n]:.3f} ', end="")
			print("]")
			index = index + 1

	#for each combination, calc jukes-cantor corrected distance
	for ia, ib in itertools.combinations(range(0,len(points)),2):
		results=list()
		for loc in range(0, len(seqs)):
			seq1 = seqs[loc][points.keys()[ia]]
			seq2 = seqs[loc][points.keys()[ib]]
			if "/" in seq1: 
				seq1 = DNAconsensus(seq1)
			if "/" in seq2:
				seq2 = DNAconsensus(seq2)
			if opts.dist == "JC69":
				results.append(jukes_cantor_distance(seq1, seq2, opts.het))
			elif opts.dist == "K2P":
				results.append(k2p_distance(seq1, seq2, opts.het))
			elif opts.dist == "PDIST":
				if opts.het:
					results.append(p_distance(seq1, seq2))
				else:
					results.append(hamming_distance(seq1, seq2))
			elif opts.dist == "TN84":
				results.append(tn84_distance(seq1, seq2, freqs[loc], opts.het))
			elif opts.dist == "TN93":
				results.append(tn93_distance(seq1, seq2, freqs[loc], opts.het))
		if opts.median:
			genmat[ia,ib] = np.median(results)
			genmat[ib,ia] = genmat[ia,ib]
		else:
			genmat[ia,ib] = np.mean(results)
			genmat[ib,ia] = genmat[ia,ib]
	#fill diagonals
	np.fill_diagonal(genmat, 0.0)
	return(genmat)

#function to return JC69-corrected p-distance
def jukes_cantor_distance(seq1, seq2, het=False):
	obs = 0.0
	if het:
		obs=p_distance(seq1, seq2, trans=False)
	else:
		obs=hamming_distance(seq1, seq2, trans=False)
	#print(obs)
	#print(1.0 - ((4.0*obs)/3.0))
	if obs > 0.75:
		obs = 0.75
	dist=-0.75*np.log(1.0 - ((4.0*obs)/3.0))
	#print(dist)
	if dist <= 0.0:
		return(0.0)
	return(dist)

#function to return Kimura 2-parameter distances
def k2p_distance(seq1, seq2, het=False):
	P=0.0
	Q=0.0
	if het:
		(P,Q)=p_distance(seq1, seq2, trans=True)
	else:
		(P,Q)=hamming_distance(seq1, seq2, trans=True)

	dist=-0.5*(np.log((1.0-(2.0*P)-Q) * math.sqrt(1.0-(2.0*Q))))
	#print(dist)
	if dist <= 0.0:
		return(0.0)
	return(dist)

#function to compute TN84 distances
def tn84_distance(seq1, seq2, freqs, het=False):
	D=0.0
	if het:
		D=p_distance(seq1, seq2, trans=False)
	else:
		D=hamming_distance(seq1, seq2, trans=False)

	ss=0.0
	for n in freqs:
		ss = ss + np.square(freqs[n])
	b=float(1.0-ss)
	dist=-b*np.log(1.0-((1.0/b)*D))
	#print(dist)
	if dist <= 0.0:
		return(0.0)
	return(dist)

#function to compute TN94 distances
def tn93_distance(seq1, seq2, freqs, het=False):
	D=0.0
	P=0.0
	Q=0.0
	P1=0.0
	P2=0.0
	if het:
		D=p_distance(seq1, seq2, trans=False, transSplit=False)
		(P1,P2,Q)=p_distance(seq1, seq2, trans=False, transSplit=True)
	else:
		D=hamming_distance(seq1, seq2, trans=False, transSplit=False)
		(P1,P2,Q)=hamming_distance(seq1, seq2, trans=False, transSplit=True)

	gR = float(freqs['g'] + freqs['a'])
	gY = float(freqs['c'] + freqs['t'])
	k1 = float((2.0*(freqs['a'])*(freqs['g'])) / gR)
	k2 = float((2.0*(freqs['t'])*(freqs['c'])) / gY)
	k3 = 2.0 * ((gR*gY) - ((freqs['a']*freqs['g'] * gY) / gR) - ((freqs['t']*freqs['c']*gR) / gY))
	w1 = float(1.0 - (P1/k1) - (Q/(2.0*gR)))
	w2 = float(1.0 - (P2/k3) - (Q/(2.0*gY)))
	w3 = float(1.0 - (Q/(2.0*gR*gY)))
	#s = float(-(k1*np.log(w1)) - (k2*np.log(w3)) - ((k3 - (2.0*gR*gY))*(np.log(w3))))
	#v = float(-2.0*gR*gY*np.log(w3))
	#R = float(s/v)
	#print(k1, k2, k3)
	#print(gR, gY)
	#print(P1, P2, Q)
	#print(w1, w2, w3)
	dist = float(-(k1*np.log(w1)) - (k2*np.log(w2)) - (k2*np.log(w3)))
	#print(dist)
	if dist <= 0.0:
		return(0.0)
	return(dist)

#function to compute nucleotide frequencies
#if ploidy = 1, ambiguities will be skipped
#if ploidy = 2, ambiguities will be resolved
def getNucFreqs(alns, ploidy):
	freqs = list()
	for aln in alns:
		temp = dict()
		allnucs = ""
		for samp in aln.keys():
			allnucs = allnucs + aln[samp].lower()
		badchars = ["?", "-", "N"]
		if ploidy == 1:
			badchars = ["?", "-", "n", "r", "y", "s", "w", "k", "m", "b", "d", "h", "v"]
		for char in badchars:
			allnucs = allnucs.replace(char, "")
		if ploidy == 2:
			iupacs = ["r", "y", "s", "w", "k", "m", "b", "d", "h", "v"]
			for ambig in iupacs:
				allnucs = allnucs.replace(ambig, "".join(get_iupac_caseless(ambig)))
			for nuc in ["a", "c", "t", "g"]:
				allnucs = allnucs.replace(char, char+char)
		total = len(allnucs)
		counts = {"a":0, "g":0,"c":0,"t":0}
		for c in allnucs:
			if c in counts:
				counts[c] += 1
		for nuc in counts.keys():
			counts[nuc] = float(counts[nuc]/total)
		freqs.append(counts)
	return(freqs)


#p distance = D / L (differences / length)
#L is the UNGAPPED distance #TODO: Maybe make this optional later
#ambigs are expanded
#when trans = true, returns two values: P (transitions/L) and Q (transversions/L)
def p_distance(seq1, seq2, trans=False, transSplit=False):
	L = min(len(seq1), len(seq2))
	D=0.0
	P=0.0
	Q=0.0
	P1=0.0
	P2=0.0
	for n1, n2 in zip(seq1.lower(), seq2.lower()):
		if n1 in ["?", "-", "n"] or n2 in ["?", "-", "n"] :
			L = L-1
			continue
		elif n1 == n2:
			continue
		else: #if n1 and n2 not equal and not gaps
			if n1 in ["a", "c", "g", "t"] and n2 in ["a", "c", "g", "t"]:
				if trans or transSplit:
					if n1 in ["a","g"] and n2 in ["a", "g"]:
						P = P+1
						P1 = P1+1
					elif n1 in ["c","t"] and n2 in ["c","t"]:
						P=P+1
						P2=P2+1
					else:
						Q=Q+1
				else:
					D = D + 1.0
				continue
			else:
				ex1 = get_iupac_caseless(n1)
				ex2 = get_iupac_caseless(n2)
				val=1.0
				if len(ex1) > 1 and len(ex2) > 1:
					val=0.25
				elif len(ex1) == 1 and len(ex2)== 1:
					val=1.0
				else:
					val=0.50
				for nuc1, nuc2 in zip(ex1, ex2):
					if nuc2 == nuc1:
						continue
					else:
						if trans or transSplit:
							if n1 in ["a","g"] and n2 in ["a", "g"]:
								P = P+val
								P1 = P1+val
							elif n1 in ["c","t"] and n2 in ["c","t"]:
								P = P+val
								P2 = P2+val
							else:
								Q=Q+val
						else:
							D=D+val
	if trans:
		transitions=0.0
		transversions=0.0
		if P > 0.0:
			transitions = float(P/L)
		if Q > 0.0:
			transversions = float(Q/L)
		return(transitions,transversions)
	elif transSplit:
		transition1=0.0
		transition2=0.0
		transversions=0.0
		if P1 > 0.0:
			transition1 = float(P1/L)
		if P2 > 0.0:
			transition2 = float(P2/L)
		if Q > 0.0:
			transversions = float(Q/L)
		return(transition1, transition2, transversions)
	else:
		if D <= 0.0:
			return(0.0)
		else:
			return(float(D/L))

#p distance = D / L (differences / length)
#gaps ignored
#ambigs are treated as alleles
def hamming_distance(seq1, seq2, trans=False, transSplit=False):
	L = min(len(seq1), len(seq2))
	D=0.0
	P=0.0
	Q=0.0
	P1=0.0
	P2=0.0
	for n1, n2 in zip(seq1.lower(), seq2.lower()):
		if n1 in ["?", "-", "n"] or n2 in ["?", "-", "n"]:
			L = L-1
			continue
		elif n1 == n2:
			continue
		else:
			if n1 != n2:
				if trans or transSplit:
					ex1 = get_iupac_caseless(n1)
					ex2 = get_iupac_caseless(n2)
					for nuc1, nuc2 in zip(ex1, ex2):
						if nuc2 == nuc1:
							continue
						else:
							if n1 in ["a","g"] and n2 in ["a", "g"]:
								P = P+1.0
								P1 = P1 + 1.0
							elif n1 in ["c","t"] and n2 in ["c","t"]:
								P=P+1.0
								P2 = P2 + 1.0
							else:
								Q=Q+1.0
				else:
					D = D + 1.0
				continue
	if trans:
		transitions=0.0
		transversions=0.0
		if P > 0.0:
			transitions = float(P/L)
		if Q > 0.0:
			transversions = float(Q/L)
		return(transitions,transversions)
	elif transSplit:
		transition1=0.0
		transition2=0.0
		transversions=0.0
		if P1 > 0.0:
			transition1 = float(P1/L)
		if P2 > 0.0:
			transition2 = float(P2/L)
		if Q > 0.0:
			transversions = float(Q/L)
		return(transition1, transition2, transversions)
	else:
		if D <= 0.0:
			return(0.0)
		else:
			return(float(D/L))

#function to make a consensus from alleles separated by "/"
def DNAconsensus(seq):
	alleles = seq.split("/")
	consens = ""
	if len(alleles) < 1:
		return(None)
	elif not all(len(x) == len(alleles[0]) for x in alleles):
		print("ERROR: Not all alleles are the same length:",alleles)
		sys.exit(1)
	elif len(alleles) == 1:
		return(alleles[0])
	else:
		for i in range(len(alleles[0])):
			nucs = ""
			for a in alleles:
				nucs += a[i]
			temp = listToSortUniqueString(nucs.upper())
			consens += reverse_iupac_case(temp)
	#print(consens)
	return(consens)

#Function to translate a string of bases to an iupac ambiguity code, retains case
def reverse_iupac_case(char):
	iupac = {
		'A':'A',
		'N':'N',
		'-':'-',
		'C':'C',
		'G':'G',
		'T':'T',
		'AG':'R',
		'CT':'Y',
		'AC':'M',
		'GT':'K',
		'AT':'W',
		'CG':'S',
		'CGT':'B',
		'AGT':'D',
		'ACT':'H',
		'ACG':'V',
		'ACGT':'N',
		'a':'a',
		'n':'n',
		'c':'c',
		'g':'g',
		't':'t',
		'ag':'r',
		'ct':'y',
		'ac':'m',
		'gt':'k',
		'at':'w',
		'cg':'s',
		'cgt':'b',
		'agt':'d',
		'act':'h',
		'acg':'v',
		'acgt':'n'
	}
	return iupac[char]
	
#Function to return sorted unique string from list of chars
def listToSortUniqueString(l):
	sl = sorted(set(l))
	return(str(''.join(sl)))

#Function to split character to IUPAC codes, assuing diploidy
def get_iupac_caseless(char):
	lower = False
	if char.islower():
		lower = True
		char = char.upper()
	iupac = {
		"A"	: ["A"],
		"G"	: ["G"],
		"C"	: ["C"],
		"T"	: ["T"],
		"N"	: ["A", "C", "G", "T"],
		"?"	: ["A", "C", "G", "T"],
		"-"	: ["A", "C", "G", "T", "-"],
		"R"	: ["A","G"],
		"Y"	: ["C","T"],
		"S"	: ["G","C"],
		"W"	: ["A","T"],
		"K"	: ["G","T"],
		"M"	: ["A","C"],
		"B"	: ["C","G","T"],
		"D"	: ["A","G","T"],
		"H"	: ["A","C","T"],
		"V"	: ["A","C","G"]
	}
	ret = iupac[char]
	if lower:
		ret = [c.lower() for c in ret]
	return ret

def phaseSnp(snp):
	nucs = get_iupac_caseless(snp)
	if len(nucs) > 2 or len(nucs) < 1:
		return("n/n")
	elif len(nucs) == 1:
		s = nucs[0] + '/' + nucs[0]
		return(s)
	else:
		return("/".join(nucs))

def path_edge_attributes(graph, path, attribute):
	return [graph[u][v][attribute] for (u,v) in zip(path,path[1:])]

#find and extract paths between points from a graph
def pathSubgraph(graph, nodes, method):
	k=nx.Graph()
	for p1, p2 in itertools.combinations(nodes.values(),2):
		try:
			#print(p1)
			#print(p2)
			#find shortest path between the two points
			path=nx.bidirectional_dijkstra(graph, p1, p2, weight=dijkstra_weight)
			#print("path:",path)
			#traverse the nodes in the path to build a minimal set of edges
			method(k, graph, nodes.values(), path[1])
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
def extractFullSubgraph(subgraph, graph, nodelist, path):
	for first, second in zip(path, path[1:]):
		if first not in subgraph:
			subgraph.add_node(first)
		if second not in subgraph:
			subgraph.add_node(second)

		dat=graph.get_edge_data(first, second)
		subgraph.add_edge(first, second, **dat)


#extracts a simplified subgraph from paths
#only keeping terminal and junction nodes
def extractMinimalSubgraph(subgraph, graph, nodelist, path):
	curr_edge = {"REACH_ID":list(), "LENGTH_KM":0}
	curr_start=None
	for first, second in zip(path, path[1:]):
		#if first is either: 1) a site node; or 2) a junction node: add to new graph
		#if second is either:a site or junction, add edge and node to new graph
		#if not, keep edge attributes for next edge
		if first in nodelist or len(graph[first])>2:
			#add to new graph if it doesn't exist
			if not curr_start:
				curr_start=first
			if first not in subgraph:
				subgraph.add_node(first)
				#print("adding",first)
		if second in nodelist or len(graph[second])>2 :
			dat=graph.get_edge_data(first, second)
			curr_edge["REACH_ID"].extend([dat["REACH_ID"]] if not isinstance(dat["REACH_ID"], list) else dat["REACH_ID"])
			curr_edge["LENGTH_KM"]+=dat["LENGTH_KM"]
			subgraph.add_node(second)
			subgraph.add_edge(curr_start, second, **curr_edge)
			curr_edge = {"REACH_ID":list(), "LENGTH_KM":0}
			curr_start = second


#function to calculate weights for Dijkstra's shortest path algorithm
#i just invert the distance, so the shortest distance segments are favored
def dijkstra_weight(attributes):
	return(attributes["LENGTH_KM"]*-1)

#Input: Tuple of [x,y] coordinates
#output: Closest node to those coordinates
def snapToNode(graph, pos):
	#rint("closest_node call:",pos)
	nodes = np.array(graph.nodes())
	node_pos = np.argmin(np.sum((nodes - pos)**2, axis=1))
	#print(nodes)
	#print("closest to ", pos, "is",tuple(nodes[node_pos]))
	return (tuple(nodes[node_pos]))


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hs:i:r:pd:a:lw:o:g', \
			["shp=", "help", "input=", "run=", "pop", "pops","dist=", "agg_method=",
			"het", "genmat=", "snp", "snps", "msat", "msats", "log", "and_log", "iterative",
			"weight=", "out=", "method=", "plots", "plot","perm=", "phased", "median",
			"diploid", "geopop", "geopops", "global_het", "haploid", "loc_agg=", 
			"pop_agg="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.shapefile = None
		self.geodb = None
		self.run = "ALL"
		self.pop = False
		self.geopop = False
		self.dist = "JC69"
		self.het = False
		self.genmat = None
		self.snps = False
		self.msat = False
		self.log = False
		self.and_log = False
		self.iterative = False
		self.weight = "CSE67"
		self.permutations = 1000
		self.method = "PEARSON"
		self.plots=False
		self.out="out"
		self.median=False
		self.ploidy=2
		self.global_het=False
		self.loc_agg = "HARM"
		self.pop_agg = "ARITH"


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg in options:
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == 's' or opt == 'shp':
				self.shapefile = arg
			elif opt == 'h' or opt == 'help':
				pass
			elif opt == 'i' or opt == 'input':
				self.geodb = arg
			elif opt == 'r' or opt == 'run':
				self.run = arg.upper()
				if self.run not in ["ALL", "GENDIST", "IBD", "STREAMDIST", "STREAMTREE", "DISTANCES", "EXHAUSTIVE", "REACHFIT"]:
					self.display_help("Invalid option", arg.upper(),"for option <-r/--run>")
			elif opt == 'p' or opt == 'pop' or opt == "pops":
				self.pop = True
			elif opt == "g" or opt == "geopop" or opt == "geopops":
				self.geopop = True
			elif opt == 'd' or opt == 'dist':
				self.dist = arg.upper()
				if self.dist not in ["PDIST", "JC69", "K2P", "TN84", "TN93", "FST", "LINFST", "JOST", "NEI83", "EUCLID", "CHORD", "GST", "GSTPRIME", "LINJOST"]:
					self.display_help("Invalid option", arg.upper(),"for option <-d/--dist>")
			elif opt == "het":
				self.het = True
			elif opt == "genmat":
				self.genmat = arg
			elif opt == "snp" or opt == "snps":
				self.snps = True
			elif opt == "msat" or opt == "msats":
				self.msat = True
				print("WARNING: Option <--msat> not yet implemented")
				sys.exit(0)
			elif opt == "l" or opt == "log":
				self.log = True
			elif opt == "and_log":
				self.and_log = True
			elif opt == "iterative":
				self.iterative = True
			elif opt == "w" or opt == "weight":
				self.weight = arg.upper()
				if self.weight not in ["CSE67", "FM67", "BEYER74"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <-w/--weight>")
			elif opt == "o" or opt == "out":
				self.out = arg
			elif opt == "stream_fit":
				self.stream_fit = True
			elif opt == "perm":
				self.permutations = int(arg)
			elif opt == "method":
				print("Sorry: Option --method is not yet implemented.")
				sys.exit(0)
				self.method = arg.upper()
				if self.method not in ["PEARSON", "SPEARMAN", "BOTH"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--method>")
			elif opt == "plot" or opt == "plots":
				self.plots = True
			elif opt == "phased":
				self.phased = True
				print("WARNING: Option <--snp> not yet implemented")
			elif opt=="median":
				self.median=True
			elif opt=="diploid":
				self.ploidy=2
			elif opt=="haploid":
				self.ploidy=1
			elif opt=="global_het":
				self.global_het=True
			elif opt=="pop_agg":
				self.pop_agg = arg.upper()
				if self.pop_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--pop_agg>")
			elif opt=="loc_agg":
				self.loc_agg = arg.upper()
				if self.loc_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN"]:
					self.display_help("Invalid option "+str(arg).upper()+" for option <--loc_agg>")
			else:
				assert False, "Unhandled option %r"%opt

		if not self.geodb:
			self.display_help("No input provided <-i/--input>")
		if not self.shapefile and self.run != "GENDIST":
			self.display_help("No shapefile provided <-s/--shp>")
		if self.snps and self.dist in ["PDIST", "JC69", "K2P", "TN93"]:
			print("\nWARNING: Can't calculate",str(self.dist),"distances with SNP data. Will instead report proportion of SNPs which are different.")
			print("NOTE: This is equivalent to an uncorrected p-distance (PDIST) calculated on the concatenated SNP data.\n")
			self.snps = False
			self.dist = "PDIST"
		if self.ploidy > 2 or self.ploidy < 1:
			self.display_help("Ploidy of",self.ploidy,"not currently allowable. Please choose 1 (haploid) or 2 (diploid)")


		###DIE FOR OPTIONS NOT YET IMPLEMENTED
		if self.genmat:
			print("Sorry: Option --genmat not yet implemented.")
			sys.exit(0)
		if self.dist in ["IBD", "EXHAUSTIVE", "REACHFIT"]:
			print("Sorry: Option --dist",self.dist," not yet implemented.")
			sys.exit(0)



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nautoStreamtree.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Computes stream distances and genetic distances for \
georeferenced DNA sequences, performs tests for isolation-by-distance, \
and uses a least-squares method to fit distances to stream segments.")
		print("""
	Input file format:
		SampleName	Data	Lat	Long	seq1	[seq2]...[seqn]
		...
		...
	--NOTE: The "DATA" column can be anything- a population/ species identifier
	  (e.g. when used with --pop), or irrelevant data (e.g. GenBank accesson
	  number, if datafile produced by my autoFetcher script)

	Mandatory arguments:
		-s,--shp	: Path to shapefile containing cleaned, contiguous stream reaches
		-i,--input	: Input .tsv file containing sample coordinates and sequences

	General options:
		--plots		: Output plots
		-o,--out	: Output prefix [default="out"]
		-h,--help	: Displays help menu
		-r,--run	: Run which steps? Options: [all, gendist, ibd, streamdist, streamtree]
			ALL			: Run all steps
			GENDIST		: Only calculate genetic distance matrix
			STREAMDIST	: Only compute pairwise stream distances
			DISTANCES	: Only compute GENDIST + STREAMDIST
			IBD		: xxxGENDIST + STREAMDIST + Mantel test
			STREAMTREE	: GENDIST + STREAMDIST + fit StreamTree model
			REACHFIT	: xxxSTREAMTREE + regress reach distances x length
			EXHAUSTIVE	: xxxRun all steps with all possible distances/ transformations
			xxx = NOT YET IMPLEMENTED
		-p,--pop		: Pool individuals based on column 2 of input file
			NOTE: The location will be taken as the centroid among individual samples
		-g,--geopop		: Pool individuals having identical coordinates

	Genetic distance options:
		-d,--dist	: Use which metric of distance? Options are:
			Substitution models (individual-based):
			  PDIST			: Uncorrected p-distances [# Differences / Length]
			  JC69 			: [default] Jukes-Cantor (1969) corrected p-distances
			  K2P			: Kimura 2-parameter distances
			  TN84			: Tajima and Nei's (1984) distance
			  TN93			: Tamura and Nei's (1993) distance
			Frequency models (when using --pop):
			  FST			: Weir and Cockerham's Fst formulation (theta)
			  GST			: Hedrick's (2005) correction of Nei's (1987) Gst [=G'st]
			  GSTPRIME		: Meirmans & Hedrick (2011) corrected G'st [=G''st]
			  LINFST		: [default] Rousset's (1997) linearized Fst [=Fst/(1-Fst)]
			  JOST			: Jost's (2008) D
			  LINJOST		: 1/1-D, where D=Jost's (2008) D
			  NEI72			: Nei's (1972) standard genetic distance 
			  NEI83			: Nei and Chesser (1983) Da
			  EUCLID		: Euclidean distance
			  CHORD			: Cavalli-Sforza and Edwards (1967) chord distance
			  --NOTE: Individual-based metrics can also be computed for
		  	          populations. You can set how these are aggregated w/ --pop_agg
			--NOTE: Multiple loci for PDIST, JC69, K2P, and EUCLID distances
		  	        will be reported using the method defined in --loc_agg
			--NOTE: TN84 will use empirical base frequencies from full per-locus alignments
		--genmat	: Skip calculation and use the provided labeled .tsv matrix
		--het		: [Boolean] Count partial differences [e.g. ind1=T, ind2=W] as a fraction
		--snp		: [Boolean] Bases should be considered as separate loci
		--global_het	: Estimate Ht using global frequencies (default is averaged over pops) 
	
	Aggregation options: 
		--pop_agg	: Define aggregator function for certain genetic distances w/ --pops:
		--loc_agg	: Define aggregator function for aggregating locus-wise distances:
		--dist_agg	: Define aggregator function for aggregating stream distances:
			All of these can take the following options:
			  ARITH		: [default] Use arithmetic mean
			  MEDIAN	: Use median distance
			  HARM		: Use harmonic mean
			  ADJHARM	: Adjusted harmonic mean (see docs)
			  GEOM		: Use geometric mean
			  MIN		: Use minimum distance
			  MAX		: Use maximum distance


	Stream distance/ IBD options:
		-l,--log	: Report natural log of stream distances
		--and_log	: Compute both absolute and log distances
		--perm		: Number of permutations for mantel test [def=1000]
		--centroid	: For populations, fit distances using a centroid of sample points
			NOTE: The StreamTree analysis can't be done using population-grouped 
			  coordinates, unless they are either 1) treated as a centroid; or 
			  2) have identical coordinates. 
			NOTE: When not using --centroid, population-wise distances will be 
			  aggregated using the method specified in --dist_agg
		--method	: Method used to report correlation between distance matrices
			Options:
			  PEARSON		: [default] Pearson's correlation coefficient
			  SPEARMAN		: Spearman's correlation coefficient
			  BOTH			: Run redundant analyses reporting both of the above

	StreamTree (see Kaliowski et al. 2008) options:
		--iterative	: Prevent negative distances using the iterative approach
		-w,--weight	: Desired weighting for least-squares fitting:
			Options:
			  FM67			: Fitch and Margoliash (1967) [w = 1/D^2]
			  BEYER74		: Beyer et al. (1974) weights [w = 1/D]
			  CSE67			: [default] Cavalli-Sforza and Edwards (1967) [w = 1]

	References:
		Beyer WM, Stein M, Smith T, Ulam S. 1974. A molecular sequence metric and
		evolutionary trees. Mathematical Biosciences. 19: 9-25.
		Cavalli-Sforza LL, Edwards AWF. 1967. Phylogenetic analysis: model and estimation
			procedures. American Journal of Human Genetics. 19: 233-257.
		Felsenstein J. 2004. Inferring Phylogenies: Chapter 11. Sunderland: Sinauer.
		Fitch WM, Margloiash E. 1967. Construction of phylogenetic trees. Science.
			155: 279-84.
		Hedrick PW. 2005. A standardized genetic differentiation measure.
			Evolution. 59: 1633–1638
		Jost L. 2008. Gst and its relatives do not measure differentiation. Molecular
			Ecology. 17: 4015-4026.
		Jukes TH, Cantor CR. 1969. Evolution of protein molecules. New York: Academic Press.
		Kimura M. 1980. A simple method for estimating evolutionary rates of base
			substitutions through comparative studies of nucleotide sequences. 
			Journal ofMolecular Evolution. 16(2): 111-120.
		Meirmans PG, Hedrick PW. 2011. Assessing population structure: Fst and related
			measures. Molecular Ecology Resources. 11: 5-18.
		Nei M (1987) Molecular Evolutionary Genetics. Columbia University Press,
			New York
		Nei M, Chesser RK. 1983. Estimation of fixation indices and gene diversities.
			Annals of Human Genetics 47(3): 253-259.
		Rossmann LA. DFLOW User's Manual. U.S. Environmental Protection Agency.
			[For description of zero-adjusted harmonic mean]
		Rousset F. 1997. Genetic differentiation and estimation of gene flow from
			F-statistics under isolation by distance. Genetics. 145: 1219-28.
		Weir BS, Cockerham CC. 1984. Estimating F-statistics for the analysis of population
			structure. Evolution. 38: 1358-1370.

	Recommended reading:
		Meirmans PG. 2012. The trouble with isolation by distance. Molecular Ecology
			21(12): 2839-46.
		Sere M, Thevenon S, Belem AMG, De Meeus T. 2017. Comparison of different genetic
			distances to test isolation by distance between populations. 2017. 119(2):55-63.
		Wright S. 1965. Isolation by distance. Genetics. 28: 114-138.

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
	main()
