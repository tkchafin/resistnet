import os
import sys 

# read vcf file
def read_vcf(vcf, concat="none", popmap=None):
	bcf_in = VariantFile(vcf)

	# set up data dict
	dat=dict()
	samples = list((bcf_in.header.samples))
	for s in samples:
		if concat == "all":
			dat[s] = list()
			dat[s].append(["",""])
		else:
			dat[s] = list()

	# if popmap, make list of samples to drop that aren't in a pop
	if popmap:
		keep = list()
		for pop in popmap:
			keep.extend(popmap[pop])
		bcf_in.subset_samples(keep)

	chrom="FIRST"
	for record in bcf_in.fetch():
		for i, sample in enumerate(record.samples):
			if concat=="all":
				loc = seq.decode(record.samples[i]['GT'], record.ref, record.alts, as_list=True)
				dat[sample][-1][0]=dat[sample][-1][0]+loc[0]
				dat[sample][-1][1]=dat[sample][-1][1]+loc[1]
			elif concat=="loc":
				if record.chrom != chrom:
					dat[sample].append(["",""])
				loc = seq.decode(record.samples[i]['GT'], record.ref, record.alts, as_list=True)
				dat[sample][-1][0]=dat[sample][-1][0]+loc[0]
				dat[sample][-1][1]=dat[sample][-1][1]+loc[1]
			else:
				loc = seq.decode(record.samples[i]['GT'], record.ref, record.alts)
				dat[sample].append(loc)
		chrom=record.chrom
	if concat != "none":
		for sample in dat:
			dat[sample] = ["/".join(x) for x in dat[sample]]
	for sample in dat.keys():
		if len(dat[sample]) < 1:
			del dat[sample]
		elif len(dat[sample]) == 1 and dat[sample][0][0] == "":
			del dat[sample]
	return(dat)

# read network
def read_network(network, shapefile):
	if network:
		print("Reading network from saved file: ", network)
		G=nx.OrderedGraph(nx.read_gpickle(network).to_undirected())
	else:
		print("Building network from shapefile:",shapefile)
		print("WARNING: This can take a while with very large files!")
		G=nx.OrderedGraph(nx.read_shp(shapefile, simplify=True, strict=True).to_undirected())
	return(G)
