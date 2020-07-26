import os
import sys

"""
Here, build an output reference string for all methods a user chose in a run.
e.g. if using StreamTree, output Kalinowski et al
For weighted LS, gen distances, etc
"""

def fetch_references(params):
	refs="\n\nIf you found this software useful for your research please cite the following (a full manuscript will come at a later date):"
	refs = refs + "\nChafin TK, Mussman SM. 2020. autoStreamTree: Automated workflows for examining patterns of genetic differentiation in stream networks. github.com/tkchafin/autoStreamTree"
	
	if params.run == "STREAMTREE" or params.run == "ALL":
		refs = refs + "\n\nPlease cite the StreamTree method as: \nKalinowski ST, MH Meeuwig, SR Narum, ML Taper (2008) Stream trees: a statistical method for mapping genetic differences between populations of freshwater organisms to the sections of streams that connect them. Canadian Journal of Fisheries and Aquatic Sciences (65:2752-2760)"
		if params.weight=="FM67":
			refs = refs + "\n\nFor the LS optimization, you chose the FM67 weighting. Please also cite: \nFitch WM, Margloiash E. 1967. Construction of phylogenetic trees. Science.155: 279-84"
		elif params.weight=="CSE67":
			refs = refs + "\n\nFor the LS optimization, you chose the CSE67 weighting. Please also cite: \nCavalli-Sforza LL, Edwards AWF. 1967. Phylogenetic analysis: model and estimation procedures. American Journal of Human Genetics. 19: 233-257"
		elif params.weight=="BEYER74":
			refs = refs + "\n\nFor the LS optimization, you chose the BEYER74 weighting. Please also cite: \nBeyer WM, Stein M, Smith T, Ulam S. 1974. A molecular sequence metric and evolutionary trees. Mathematical Biosciences. 19: 9-25"
	if params.run!="STREAMDIST" and not params.genmat:
		if params.dist != "PDIST" and params.dist!="EUCLID":
			refs = refs + "\n\nPlease cite the following for your genetic distance calculations:"
			if params.dist=="JC69":
				refs = refs + "\nJukes TH, Cantor CR. 1969. Evolution of protein molecules. New York: Academic Press."
			elif params.dist=="K2P":
				refs = refs + "\nKimura M. 1980. A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of Molecular Evolution. 16(2): 111-120."
			elif params.dist=="TN84":
				refs = refs + "\nTajima F, Nei M. 1984. Estimation of evolutionary distance between nucleotide sequences. Molecular Biology and Evolution 1:269-285"
			elif params.dist=="TN93":
				refs = refs + "\nTamura K, Nei M. 1993. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Molecular Biology and Evolution. 10(3):512-526."
			elif params.dist=="FST":
				refs = refs + "\nWeir BS, Cockerham CC. 1984. Estimating F-statistics for the analysis of population structure. Evolution. 38: 1358-1370."
			elif params.dist=="GST":
				refs = refs + "\nHedrick PW. 2005. A standardized genetic differentiation measure. Evolution. 59: 1633–1638"
				refs = refs + "\nNei M. 1987. Molecular Evolutionary Genetics. Columbia University Press, New York"
			elif params.dist=="GSTPRIME":
				refs = refs + "\nMeirmans PG, Hedrick PW. 2011. Assessing population structure: Fst and related measures. Molecular Ecology Resources. 11: 5-18."
			elif params.dist=="LINFST":
				refs = refs + "\nRousset F. 1997. Genetic differentiation and estimation of gene flow from F-statistics under isolation by distance. Genetics. 145: 1219-28."
				refs = refs + "\nWeir BS, Cockerham CC. 1984. Estimating F-statistics for the analysis of population structure. Evolution. 38: 1358-1370."
			elif params.dist=="JOST":
				refs = refs + "\nJost L. 2008. Gst and its relatives do not measure differentiation. Molecular Ecology. 17: 4015-4026."
			elif params.dist=="NEI72":
				refs = refs + "\nNei M. 1972. Genetic distance between populations. American Naturalist. 106: 283-292."
			elif params.dist=="NEI83":
				refs = refs + "\nNei M, Chesser RK. 1983. Estimation of fixation indices and gene diversities. Annals of Human Genetics 47(3): 253-259."
			elif params.dist=="CHORD":
				refs = refs + "\nCavalli-Sforza LL, Edwards AWF. 1967. Phylogenetic analysis: model and estimation procedures. American Journal of Human Genetics. 19: 233-257."
	
	if params.run=="IBD":
		refs = refs + "\n\nFor the Mantel test, please cite:\nMantel N. 1967. The detection of disease clustering and a generalized regression approach. Cancer Research 27(2): 209-220."
	
	if params.clusterpop:
		refs = refs + "\n\nFor the DBSCAN approach, please cite:"
		refs = refs + "\nEster M, Kriegel HP, Sander J, Xu X. 1996. A density-based algorithm for discovering  clusters in large spatial databases with noise. IN: Simoudis E, Han J, Fayyad UM. (eds.). Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231."
		refs = refs + "\nPedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, Blondel M, Prettenhofer P, Weiss R, Dubourg V, Vanderplas J. 2011. Scikit-learn: Machine learning in Python. The Journal of machine Learning research. 1(12):2825-30"
	
	if params.pop_agg=="ADJHARM" or params.loc_agg=="ADJHARM":
		refs = refs + "\n\nFor the adjusted harmonic mean, please cite: \nRossmann LA. DFLOW User's Manual. U.S. Environmental Protection Agency."

	if params.run != "GENDIST":
		refs = refs + "\n\nFor the network methods, please cite: \nHagberg A, Swart P, S Chult D. 2008. Exploring network structure, dynamics, and function using NetworkX. Los Alamos National Lab.(LANL), Los Alamos, NM"

	refs = refs + "\n"

	return(refs)
