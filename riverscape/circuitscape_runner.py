import sys
import os
from julia.api import Julia
from julia import Base, Main
import pandas as pd
from collections import OrderedDict
from io import StringIO 

import riverscape.MLPE as mlpe_rga
#multiple mutation types: https://stackoverflow.com/questions/47720921/deap-toolbox-to-consider-different-types-and-ranges-of-genes-in-mutation-and-cr

"""
TO-DO:
Add options to only analyse sampling points as focal nodes (i.e., not junctions)
Add options to set certain nodes as 'sources'?
figure out how to suppress the output from Julia/ circuitscape
"""

def parseEdgewise(oname):
	pass

def parsePairwise(oname, gendist):
	pw=pd.read_csv((str(oname)+"_resistances.out"), header=0, index_col=0, sep=" ").to_numpy()
	res = mlpe_rga.MLPE_R(gendist, pw, scale=True)
	return(res)
	
	
def evaluateIni(jl, oname):
	ini_path = str(oname)+".ini"
	jl.eval(str("@suppress begin\ncompute(\"" + str(ini_path)+"\")\nend"))

#function to write inputs for circuitscape
def writeCircuitScape(oname, graph, points, resistance, focalPoints=False, fromAttribute=None):
	if fromAttribute is None:
		node_dict=dict()
		edge_idx=0
		node_idx=1
		net_output=""
		pts_output=""
		kept=0
		#for each edge
		#print("Number of points:",len(points))
		for edge in graph.edges:
			#get nodes on either side and index them
			#get resistance 
			#add all to output string
			#NOTE: Points should be an OrderedDict, so this should work fine
			#print(edge[0], " -- ", edge[1])
			if edge[0] not in node_dict.keys():
				node_dict[edge[0]] = node_idx
				if focalPoints and edge[0] not in points.keys():
					pass
				else:
					pts_output += str(node_idx) + ".0\n"
					kept+=1
				node_idx+=1
			if edge[1] not in node_dict.keys():
				node_dict[edge[1]] = node_idx
				if focalPoints and edge[1] not in points.keys():
					pass
				else:
					pts_output += str(node_idx) + ".0\n"
					kept+=1
				node_idx+=1
			
			net_output += str(node_dict[edge[0]]) + ".0\t" + str(node_dict[edge[1]]) + ".0\t" + str(float(resistance[edge_idx])) + "\n"
			edge_idx+=1
		#print("Number of focal points:",kept)
		with open((str(oname) + ".graph_resistances.txt"), "w") as ofh:
			ofh.write(net_output)
			ofh.close()
		with open((str(oname) + ".focal_nodes.txt"), "w") as pfh:
			pfh.write(pts_output)
			pfh.close()

def writeIni(oname, cholmod=False, parallel=1):
	with open((str(oname)+".ini"), "w") as ini:
		ini.write("[Options for advanced mode]\n")
		ini.write("ground_file_is_resistances = False\n")
		ini.write("source_file = ")
		ini.write((str(oname)+".graph_resistances.txt\n"))
		ini.write("remove_src_or_gnd = keepall\n")
		ini.write("ground_file = None\n")
		ini.write("use_unit_currents = False\n")
		ini.write("use_direct_grounds = False\n\n")
		
		ini.write("[Calculation options]\n")
		ini.write("low_memory_mode = False\n")
		if cholmod:
			ini.write("solver = cholmod\n")
		else:
			ini.write("solver = cg+amg\n")
		ini.write("print_timings = False\n")
		if parallel>1:
			ini.write("parallelize = True\n")
			ini.write(("max_parallel ="+str(parallel)+"\n"))
		else:
			ini.write("parallelize = False\n")
			ini.write("max_parallel = 0\n")
		
		ini.write("[Options for pairwise and one-to-all and all-to-one modes]\n")
		ini.write("included_pairs_file = None\n")
		ini.write("use_included_pairs = False\n")
		ini.write("point_file = ")
		ini.write((str(oname)+".focal_nodes.txt\n\n"))

		ini.write("[Output options]\n")
		ini.write("write_cum_cur_map_only = False\n")
		ini.write("log_transform_maps = False\n")
		ini.write("output_file = ")
		ini.write((str(oname)+".output.txt\n\n"))
		ini.write("write_max_cur_maps = False\n")
		ini.write("write_volt_maps = False\n")
		ini.write("set_null_currents_to_nodata = True\n")
		ini.write("set_null_voltages_to_nodata = True\n")
		ini.write("compress_grids = False\n")
		ini.write("write_cur_maps = False\n")

		ini.write("[Short circuit regions (aka polygons)]\n")
		ini.write("use_polygons = False\n")
		ini.write("polygon_file = None\n")

		ini.write("[Connection scheme for raster habitat data]\n")
		ini.write("connect_four_neighbors_only = True\n")
		ini.write("connect_using_avg_resistances = True\n")

		ini.write("[Habitat raster or graph]\n")
		ini.write("habitat_file = ")
		ini.write((str(oname)+".graph_resistances.txt\n"))
		ini.write("habitat_map_is_resistances = True\n")

		ini.write("[Options for one-to-all and all-to-one modes]\n")
		ini.write("use_variable_source_strengths = False\n")
		ini.write("variable_source_file = None\n")

		ini.write("[Version]\n")
		ini.write("version = 3.5.5\n")

		ini.write("[Mask file]\n")
		ini.write("use_mask = False\n")
		ini.write("mask_file = None\n")

		ini.write("[Circuitscape mode]\n")
		ini.write("data_type = network\n")
		ini.write("scenario = pairwise\n")
		
		ini.close()



	