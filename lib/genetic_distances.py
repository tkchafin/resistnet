import sys
import os
import itertools
import math
import scipy
import numpy as np
import pandas as pd
from sortedcontainers import SortedDict

import aggregators as agg

def getPopGenMat(dist, indmat, popmap, dat, seqs, pop_agg="ARITH", loc_agg="ARITH", ploidy=2, global_het=False):
	#make matrix
	genmat = np.zeros((len(popmap),len(popmap)))
	#establish as nan
	genmat[:] = np.nan
	#for each combination, either average ind distances or calc freq dist 
	for ia,ib in itertools.combinations(range(0,len(popmap)),2):
		if dist in ["JC69", "K2P", "PDIST", "TN84", "TN93"]:
			#print(popmap.keys()[ia], popmap.values()[ia])
			inds1 = [dat.index(x) for x in popmap.values()[ia]]
			inds2 = [dat.index(x) for x in popmap.values()[ib]]
			#print(inds1)
			genmat[ia,ib] = (agg.aggregateDist(pop_agg, ([indmat[i, j] for i in inds1 for j in inds2])))
			genmat[ib,ia] = genmat[ia,ib]
		elif dist == "JOST" or dist == "LINJOST":
			results=list()
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					continue
				results.append(twoPopJostD(seqs1, seqs2, ploidy, global_het))
			#print(results)
			if len(results) > 1:
				if dist=="LINJOST":
					D = agg.aggregateDist(loc_agg, results)
					genmat[ia,ib] = genmat[ib,ia] = (1.0/(1.0-float(D)))
				else:
					genmat[ia,ib] = genmat[ib,ia] = agg.aggregateDist(loc_agg, results)
			elif len(results) < 1:
				#print("ERROR: population",popmap.values()[ia],"or",popmap.values()[ib],"lacks any data")
				raise ValueError
			else:
				genmat[ia,ib] = genmat[ib,ia] = results[0]
		elif dist == "GST" or dist == "GSTPRIME":
			HT=list()
			HS=list()
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					continue
				if dist == "GST" or "GSTPRIME":
					(ht, hs) = twoPopHtHs(seqs1, seqs2, ploidy, global_het)
					HT.append(ht)
					HS.append(hs)
			Ht_global = np.mean(HT)
			Hs_global = np.mean(HS)
			if dist == "GST":
				if Ht_global <= 0.0:
					genmat[ia,ib] = genmat[ib,ia] = 0.0
				Gst = ((Ht_global - Hs_global) / Ht_global )
				GprimeST = ((Gst * (1.0 + Hs_global)) / (1.0 - Hs_global))
				genmat[ia,ib] = genmat[ib,ia] = GprimeST
			elif dist == "GSTPRIME":
				Ghedrick = ((2.0*(Ht_global - Hs_global)) / (((2.0*Ht_global) - Hs_global) * (1.0 - Hs_global)))
				genmat[ia,ib] = genmat[ib,ia] = Ghedrick
		elif dist == "FST" or dist == "LINFST":
			num = list() #numerator; a 
			denom = list() #denominator; a*b*c
			for loc in range(0, len(seqs)):
				seqs1 = cleanInds([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = cleanInds([seqs[loc][x] for x in popmap.values()[ib]])
				if not all("/" in x for x in seqs1) or not all("/" in x for x in seqs1):
					print("ERROR: FST estimates require phased data.")
					sys.exit(1)
				if len(seqs1) == 0 or len(seqs2) == 0:
					#print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
					continue
				(n, d) = twoPopWeirCockerhamFst(seqs1, seqs2)
				num.append(n)
				denom.append(d)
			if len(num) <= 0 or len(denom) <= 0:
				#print("ERROR (twoPopWeirCockerhamFst): No data for pops "+ia+" and "+ib+".")
				raise ValueError
			theta = np.sum(num) / np.sum(denom)
			if dist == "FST":
				genmat[ia,ib] = genmat[ib,ia] = theta
			elif dist == "LINFST":
				genmat[ia,ib] = genmat[ib,ia] = (theta / (1-theta))
		elif dist == "NEI83":
			results = list()
			loci = 0.0
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					#print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
					continue
				loci += 1.0
				results.append(twoPopNeiDa(seqs1, seqs2))
			Da = (1.0 - (np.sum(results) / loci))
			genmat[ia,ib] = genmat[ib,ia] = Da
		elif dist == "EUCLID":
			results = list()
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					#print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
					continue
				results.append(twoPopEuclidDist(seqs1, seqs2))
			euclid = np.sum(results)
			genmat[ia,ib] = genmat[ib,ia] = euclid
		elif dist == "CHORD":
			num = list()
			denom = list()
			for loc in range(0, len(seqs)):
				seqs1 = getAlleles([seqs[loc][x] for x in popmap.values()[ia]])
				seqs2 = getAlleles([seqs[loc][x] for x in popmap.values()[ib]])
				if not cleanList(set(seqs1), ["n", "N", "-", "?"]) or not cleanList(set(seqs2), ["n", "N", "-", "?"]):
					#print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
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

#function computes pairwise JC69-corrected genetic distances
def getGenMat(dist, points, seqs, ploidy, het, loc_agg):
	#make matrix
	genmat = np.zeros((len(points),len(points)))
	#establish as nan
	genmat[:] = np.nan

	#for models which relax equal nuc frequencies, get global frequencies for each locus
	#freqs will be a list of loci, with each locus as a dist of freqs
	if dist in ["TN84", "TN93"]:
		freqs = getNucFreqs(seqs, ploidy)
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
			if dist == "JC69":
				results.append(jukes_cantor_distance(seq1, seq2, het))
			elif dist == "K2P":
				results.append(k2p_distance(seq1, seq2, het))
			elif dist == "PDIST":
				if het:
					results.append(p_distance(seq1, seq2))
				else:
					results.append(hamming_distance(seq1, seq2))
			elif dist == "TN84":
				results.append(tn84_distance(seq1, seq2, freqs[loc], het))
			elif dist == "TN93":
				results.append(tn93_distance(seq1, seq2, freqs[loc], het))
		#aggregate results across loci
		genmat[ia,ib] = agg.aggregateDist(loc_agg, results)
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
def twoPopJostD(seqs1, seqs2, ploidy, global_het=False):
	if global_het:
		Ht = getGlobalHet(seqs1+seqs2)
	else:
		Ht = getAverageHet(seqs1, seqs2)
	Hs = np.mean([getGlobalHet(seqs1),getGlobalHet(seqs2)])
	harmN = scipy.stats.hmean([(len(seqs1)/ploidy), (len(seqs1)/ploidy)])
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
def twoPopHtHs(seqs1, seqs2, ploidy, global_het=False):
	if global_het:
		Ht = getGlobalHet(seqs1+seqs2)
	else:
		Ht = getAverageHet(seqs1, seqs2)
	Hs = np.mean([getGlobalHet(seqs1),getGlobalHet(seqs2)])
	harmN = scipy.stats.hmean([(len(seqs1)/ploidy), (len(seqs1)/ploidy)])
	Hs_est = Hs * ((2.0*harmN)/((2.0*harmN)-1.0))
	Ht_est = Ht + (Hs / (harmN*2.0*2.0))
	# print("Hs_est:",Hs_est)
	# print("Ht_est:",Ht_est)
	# if Ht_est == 0.0:
	# 	return(0.0)
	# Gst = ((Ht_est - Hs_est) / Ht_est )
	# GprimeST = ((Gst * (1.0 + Hs_est)) / (1.0 - Hs_est))
	return((Ht_est, Hs_est))

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