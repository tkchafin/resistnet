import os 
import sys

import pandas as pd
import numpy as np
import networkx as nx
import pickle
import math
import itertools
from sortedcontainers import SortedDict

#Input: Tuple of [x,y] coordinates
#output: Closest node to those coordinates
def snap_to_node(graph, pos):
    nodes = np.array(graph.nodes())
    node_pos = np.argmin(np.sum((nodes - pos)**2, axis=1))
    return (tuple(nodes[node_pos]))


def check_dataframe_columns(df):
    for col in df.columns:
        unique_values = df[col].dropna().unique()
        
        if len(unique_values) == 0 or (len(unique_values) == 1 and np.isnan(unique_values[0])):
            print(f"Column '{col}' has all NaN values.")
            return False
        elif len(unique_values) == 1:
            print(f"Column '{col}' has all the same value: {unique_values[0]}")
            return False
    return True


def scrub_bad_columns(df, variables, verbose=False):
    bad_columns = []
    
    for col in df.columns:
        unique_values = df[col].dropna().unique()
        
        if len(unique_values) == 0 or (len(unique_values) == 1 and np.isnan(unique_values[0])):
            if verbose:
                print(f"Column '{col}' has all NaN values.")
            bad_columns.append(col)
        elif len(unique_values) == 1:
            if verbose:
                print(f"Column '{col}' has all the same value: {unique_values[0]}")
            bad_columns.append(col)
    
    # Remove bad columns from the data frame and the predictors list
    cleaned_df = df.drop(bad_columns, axis=1)
    cleaned_variables = [v for v in variables if v not in bad_columns]

    return cleaned_df, cleaned_variables


def read_pickle_network(pik):
    with open(pik, 'rb') as f:
        _G = nx.Graph(pickle.load(f)).to_undirected()
    return(_G)


#utility function to test if two elements are consecutive in list (irrespective of order)
def find_pair(list, x, y):
    if x not in list or y not in list:
        return(False)
    elif abs(list.index(x)-list.index(y)) == 1:
        return(True)
    else:
        return(False)

def nx_to_df(G):
    l=list()
    for p1, p2, e in G.edges(data=True):
        e["left"] = p1
        e["right"] = p2
        l.append(e)
    df=pd.DataFrame(l)
    return(df)


def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]


def read_points_table(pfile):
    points = pd.read_csv(pfile, sep="\t", header=0)

    # Check for required columns
    required_columns = {'sample', 'lat', 'long'}
    if not required_columns.issubset(points.columns):
        missing_cols = required_columns - set(points.columns)
        raise ValueError(f"Missing required columns in the points data: {', '.join(missing_cols)}")
    
    return(points)

def read_points_flatten(pfile):
    points = read_points_table(pfile)  # Assuming this returns a DataFrame

    d = SortedDict()
    for _, row in points.iterrows():
        name = row['sample']  
        lat = float(row['lat'])  
        long = float(row['long'])  
        coords = (long, lat)

        if coords not in d:
            d[coords] = []
        d[coords].append(name)

    return d


#utility function to calculate number of combinations n choose k
def nCr(n,k):
    f = math.factorial
    return f(n) // f(k) // f(n-k)

def write_numpy_matrix(mat, ofh):
    with np.printoptions(precision=0, suppress=True):
        np.savetxt(ofh, mat, delimiter="\t")


def path_edge_attributes(graph, path, attribute):
    return [graph[u][v][attribute] for (u,v) in zip(path,path[1:])]


def write_edges(out, edge, ids):
    df=pd.DataFrame(list(zip(ids, edge)), columns=["EDGE_ID", "Resistance"])
    df.to_csv(out, sep="\t", header=True, index=False)
    return(out)

def write_matrix(out, mat, ids):
    df=pd.DataFrame(mat, columns=ids, index=ids)
    df.to_csv(out, sep="\t", header=True, index=True)
    return(out)

def get_lower_tri(mat):
	n=mat.shape[0]
	i=np.tril_indices(n,-1)
	return(mat[i])

def to_from_(pops):
	to = list()
	frm = list()
	for ia, ib in itertools.combinations(range(1,pops+1),2):
		to.append(ia)
		frm.append(ib)
	t = to[pops-2]
	tt = frm[pops-2]
	to[pops-2] = tt
	frm[pops-2] = t
	df = pd.DataFrame({"pop1" : to, "pop2" : frm})
	return(df)
