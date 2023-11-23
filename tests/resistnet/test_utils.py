import os
import pytest
import numpy as np
import pandas as pd
import networkx as nx
import random
import tempfile
import resistnet
import resistnet.utils as utils
from sortedcontainers import SortedDict


@pytest.fixture
def points_file_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'test.pointCoords.txt')
    return file_path


@pytest.fixture
def sample_graph():
    # Create an undirected graph
    G = nx.Graph()

    # Add nodes with (long, lat) tuples as names
    nodes = [
        (random.uniform(-180, 180), random.uniform(-90, 90)) for _ in range(4)
    ]
    G.add_nodes_from(nodes)

    # Add edges with EDGE_ID and Resistance attributes
    for i in range(len(nodes) - 1):
        edge_id = i + 1
        resistance = random.uniform(1, 10)  # Example resistance value
        G.add_edge(
            nodes[i], nodes[i+1], EDGE_ID=edge_id, Resistance=resistance
        )

    return G


# Test for get_lower_tri
def test_get_lower_tri():
    mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    expected = np.array([4, 7, 8])
    assert np.array_equal(utils.get_lower_tri(mat), expected)


# Test for to_from_
def test_to_from_():
    pops = 4
    expected = pd.DataFrame(
        {"pop1": [1, 1, 4, 2, 2, 3],
         "pop2": [2, 3, 1, 3, 4, 4]}
    )
    result = utils.to_from_(pops)
    pd.testing.assert_frame_equal(result, expected)


# Test for sample_nodes using the sample_graph fixture
def test_sample_nodes(sample_graph):
    samples = 2
    result = utils.sample_nodes(sample_graph, samples)
    assert isinstance(result, SortedDict)
    assert len(result) == samples
    assert all(node in sample_graph.nodes for node in result)


# Test for write_edges
def test_write_edges():
    edge = [1.0, 2.5, 3.2]  # Example edge resistance values
    ids = [101, 102, 103]    # Example edge IDs

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        out_path = tmp.name

    # Call the function and then read the output
    utils.write_edges(out_path, edge, ids)
    df = pd.read_csv(out_path, sep="\t")

    # Clean up
    os.remove(out_path)

    # Assertions
    expected = pd.DataFrame(
        list(zip(ids, edge)), columns=["EDGE_ID", "Resistance"]
    )
    pd.testing.assert_frame_equal(df, expected)


# Test for write_matrix
def test_write_matrix():
    mat = [[1, 2], [3, 4]]   # Example matrix
    ids = ['row1', 'row2']    # Example identifiers

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        out_path = tmp.name

    # Call the function and then read the output
    utils.write_matrix(out_path, mat, ids)
    df = pd.read_csv(out_path, sep="\t", index_col=0)

    # Clean up
    os.remove(out_path)

    # Assertions
    expected = pd.DataFrame(mat, columns=ids, index=ids)
    pd.testing.assert_frame_equal(df, expected)


# test nCr function 
def test_nCr():
    assert utils.nCr(5, 2) == 10
    assert utils.nCr(10, 5) == 252
    assert utils.nCr(0, 0) == 1


# test write_numpy_matrix
def test_write_numpy_matrix():
    mat = np.array([[1.0, 0.0], [1.0, 0.0]])

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        out_path = tmp.name

    # Call the function and then read the output
    utils.write_numpy_matrix(mat, out_path)
    read_mat = np.loadtxt(out_path, delimiter="\t")

    # Clean up
    os.remove(out_path)

    # Assertions
    expected_written_mat = np.array([[1., 0.], [1., 0.]])
    assert np.array_equal(read_mat, expected_written_mat)


# test path_edge_attributes
def test_path_edge_attributes(sample_graph):
    path = list(sample_graph.nodes)[:3]
    expected_attributes = []
    for i in range(len(path) - 1):
        resistance = sample_graph[path[i]][path[i + 1]]['Resistance']
        expected_attributes.append(resistance)

    result = utils.path_edge_attributes(sample_graph, path, 'Resistance')
    assert result == expected_attributes


def test_read_points_table_success(points_file_path):
    # Test if the function reads the file correctly
    result = utils.read_points_table(points_file_path)
    assert isinstance(result, pd.DataFrame)
    assert set(result.columns) == {'sample', 'lat', 'long'}
    assert len(result) > 0  # Assuming the file is not empty


# Define different cases of invalid file contents
@pytest.mark.parametrize("file_contents", [
    ("sample\tlatitude\n1\t2.5"),  # Missing 'long' column
    ("lat\tlong\n2.5\t3.5"),       # Missing 'sample' column
    ("sample\tlong\n1\t3.5"),      # Missing 'lat' column
    ("burk\t2.5\t3.5\n"),          # no header
])
def test_read_points_table_failure(file_contents):
    # Create a temporary file with invalid contents
    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as tmp:
        tmp.write(file_contents)
        tmp.flush()

        # Test if the function raises ValueError for missing columns
        with pytest.raises(ValueError) as exc_info:
            utils.read_points_table(tmp.name)
        assert "Missing required columns" in str(exc_info.value)
