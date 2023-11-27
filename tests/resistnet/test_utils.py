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
def pickle_network_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'test.network')
    return file_path


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


def test_haversine():
    # Define two known points (Latitude, Longitude)
    coord1 = (48.8566, 2.3522)  # Paris
    coord2 = (51.5074, -0.1278)  # London

    # Known distance between Paris and London in kilometers
    known_distance = 344

    # Calculate the distance using the haversine function
    calculated_distance = utils.haversine(coord1, coord2)

    # Assert that the calculated distance is close to the known distance
    assert abs(calculated_distance - known_distance) < 1  # 1 km tolerance


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


# Define different cases fo read_points_table error checking
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


# test read_points_flatten
def test_read_points_flatten(points_file_path):
    # Read and flatten points data
    flattened_data = utils.read_points_flatten(points_file_path)

    # Verify that the result is a SortedDict
    assert isinstance(flattened_data, SortedDict)

    # Check if keys are tuples of coordinates and values are lists of samples
    for key, value in flattened_data.items():
        assert isinstance(key, tuple) and len(key) == 2
        assert all(isinstance(item, str) for item in value)


# test unique
def test_unique():
    sequence = [1, 2, 2, 3, 1, 4, 5, 4]
    expected_result = [1, 2, 3, 4, 5]

    result = utils.unique(sequence)
    assert result == expected_result

    # Test with different types of sequences (strings, etc.)
    sequence_str = ['a', 'b', 'a', 'c', 'b']
    expected_result_str = ['a', 'b', 'c']
    assert utils.unique(sequence_str) == expected_result_str

    # Test with empty sequence
    assert utils.unique([]) == []


# test nx_to_df 
def test_nx_to_df(sample_graph):
    # Convert the sample graph to DataFrame
    df = utils.nx_to_df(sample_graph)

    # Check if the result is a DataFrame
    assert isinstance(df, pd.DataFrame)

    # Verify DataFrame columns
    expected_columns = {'left', 'right', 'EDGE_ID', 'Resistance'}
    assert set(df.columns) == expected_columns


# test find_pait
def test_find_pair():
    lst = [1, 2, 3, 4, 5]

    # Check for pairs that are consecutive
    assert utils.find_pair(lst, 1, 2) is True
    assert utils.find_pair(lst, 3, 4) is True

    # Check order irrelevance
    assert utils.find_pair(lst, 4, 3) is True

    # Check for pairs that are not consecutive
    assert utils.find_pair(lst, 1, 3) is False

    # Check for elements not in list
    assert utils.find_pair(lst, 6, 7) is False


# test read_pickle_network
def test_read_pickle_network(pickle_network_path):
    # Read the network from the pickle file
    graph = utils.read_pickle_network(pickle_network_path)

    # Check if the result is a NetworkX Graph
    assert isinstance(graph, nx.Graph)

    # Verify that the graph is undirected
    assert not graph.is_directed()

    # Check that all nodes are tuples (likely representing coordinates)
    assert all(isinstance(node, tuple) for node in graph.nodes)

    # Check that edges have the required attributes
    if len(graph.edges) > 0:
        u, v = next(iter(graph.edges))
        assert 'HYRIV_ID' in graph[u][v]
        assert 'LENGTH_KM' in graph[u][v]


# Define test cases
@pytest.mark.parametrize("data, variables, expected_result, \
                         expected_vars, expected_exception", [
    # Test with columns containing NaNs and single unique values
    (
        pd.DataFrame({
            'A': [1, 1, 1], 'B': [np.nan, np.nan, np.nan], 'C': [1, 2, 3]
        }),
        ['A', 'B', 'C'],
        pd.DataFrame({'C': [1, 2, 3]}),
        ['C'],
        None
    ),
    # Test with valid columns
    (
        pd.DataFrame({'A': [1, 2, 1], 'B': [3, 4, 5]}),
        ['A', 'B'],
        pd.DataFrame({'A': [1, 2, 1], 'B': [3, 4, 5]}),
        ['A', 'B'],
        None
    ),
    # Valid columns but invalid variables
    (
        pd.DataFrame({'A': [1, 2, 1], 'B': [3, 4, 5]}),
        ['D', 'E'],
        None,
        None,
        ValueError
    ),
    # Valid variables but all columns fail
    (
        pd.DataFrame({'A': [1, 1, 1], 'B': [1, 1, 1]}),
        ['A', 'B'],
        None,
        None,
        ValueError
    ),
    # Test with empty DataFrame (should raise ValueError)
    (
        pd.DataFrame(),
        [],
        None,
        None,
        ValueError
    )
])
def test_scrub_bad_columns(
        data, variables, expected_result, expected_vars, expected_exception):
    if expected_exception:
        # Check if the function raises the expected exception
        with pytest.raises(expected_exception):
            utils.scrub_bad_columns(data, variables)
    else:
        # Call the scrub_bad_columns function
        cleaned_df, cleaned_variables = utils.scrub_bad_columns(
            data, variables
        )

        # Assert the cleaned DataFrame and variables list
        pd.testing.assert_frame_equal(cleaned_df, expected_result)
        assert cleaned_variables == expected_vars


# Define test cases
@pytest.mark.parametrize("data, expected_result", [
    # DataFrame with valid columns
    (
        pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}),
        True
    ),
    # DataFrame with a column containing only NaN values
    (
        pd.DataFrame({'A': [1, 2, 3], 'B': [np.nan, np.nan, np.nan]}),
        False
    ),
    # DataFrame with a column containing a single unique value
    (
        pd.DataFrame({'A': [1, 1, 1], 'B': [2, 3, 4]}),
        False
    ),
    # Test with empty DataFrame (should raise ValueError)
    (
        pd.DataFrame(),
        ValueError
    )
])
def test_check_dataframe_columns(data, expected_result):
    if expected_result is ValueError:
        # Check if the function raises the expected exception
        with pytest.raises(ValueError):
            utils.check_dataframe_columns(data)
    else:
        # Call the check_dataframe_columns function
        result = utils.check_dataframe_columns(data)

        # Assert the result
        assert result == expected_result


def test_snap_to_node(sample_graph):
    # Test with positions that exactly match node coordinates
    for node in sample_graph.nodes:
        closest_node = utils.snap_to_node(sample_graph, node)
        assert closest_node == node

    # Test with slightly shifted positions
    for node in sample_graph.nodes:
        shifted_pos = tuple(np.array(node) + np.random.normal(0, 0.01, 2))
        closest_node = utils.snap_to_node(sample_graph, shifted_pos)
        assert closest_node == node

    # Test with empty graph (should raise ValueError)
    empty_graph = nx.Graph()
    with pytest.raises(ValueError):
        utils.snap_to_node(empty_graph, (0, 0))

    # Test with invalid position (should raise ValueError)
    with pytest.raises(ValueError):
        utils.snap_to_node(sample_graph, "invalid_pos")