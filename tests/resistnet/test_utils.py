import pytest
import numpy as np
import pandas as pd
import networkx as nx
from sortedcontainers import SortedDict
import resistnet.utils as utils


# Test for get_lower_tri
def test_get_lower_tri():
    mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    expected = np.array([4, 7, 8])
    assert np.array_equal(utils.get_lower_tri(mat), expected)


# Test for to_from_
def test_to_from_():
    pops = 3
    expected = pd.DataFrame({"pop1": [1, 1, 2], "pop2": [2, 3, 3]})
    result = utils.to_from_(pops)
    pd.testing.assert_frame_equal(result, expected)


# Test for sample_nodes
def test_sample_nodes():
    G = nx.path_graph(4)  # Graph with 4 nodes
    samples = 2
    result = utils.sample_nodes(G, samples)
    assert isinstance(result, SortedDict)
    assert len(result) == samples
    assert all(node in G.nodes for node in result)
    