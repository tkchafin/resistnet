import os
import pytest
import tempfile
import resistnet
import networkx as nx
import pandas as pd
import numpy as np
from unittest.mock import patch
from resistnet.resistance_network import ResistanceNetwork
from resistnet.resistance_network import ResistanceNetworkWorker
from resistnet.resistance_network import SimResistanceNetwork
import resistnet.utils as utils


@pytest.fixture
def fake_spec_file():
    spec_content = pd.DataFrame({
        "VAR": ["LENGTH_KM"],
        "WEIGHT": [1.0],
        "TRANSFORM": [0],
        "SHAPE": [0]
    })
    with tempfile.NamedTemporaryFile(
            mode='w+', suffix='.tsv', delete=False) as tmpfile:
        spec_content.to_csv(tmpfile.name, sep='\t', index=False)
        yield tmpfile.name
        os.remove(tmpfile.name)


@pytest.fixture
def network_graph_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'test.network')
    return file_path


@pytest.fixture
def shapefile_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'test.shp')
    return file_path


@pytest.fixture
def coords_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'test.pointCoords.txt')
    return file_path


@pytest.fixture
def inmat_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'test.popGenDistMat.txt')
    return file_path


@pytest.fixture
def sim_resistance_network_params(network_graph_path):
    return {
        "network": network_graph_path,
        "reachid_col": "EDGE_ID",
        "length_col": "LENGTH_KM",
        "verbose": True
    }


@pytest.fixture
def resistance_network_fixture(
        network_graph_path, shapefile_path, coords_path, inmat_path):
    # Use temporary directory for output
    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "output")

        network = ResistanceNetwork(
            network=network_graph_path,
            shapefile=shapefile_path,
            coords=coords_path,
            variables=["run_mm_cyr"],
            inmat=inmat_path,
            agg_opts={"run_mm_cyr": "ARITH"},
            pop_agg="ARITH",
            reachid_col="EDGE_ID",
            length_col="LENGTH_KM",
            out=output_prefix,
            verbose=False
        )

        yield network


def test_initialize_predictors(resistance_network_fixture):
    resistance_network_fixture.initialize_predictors()

    # Assertions
    # Ensure that predictors are initialized and variables are updated
    assert resistance_network_fixture._predictors is not None
    assert len(resistance_network_fixture.variables) > 0


def test_read_network_with_existing_file(
        resistance_network_fixture, monkeypatch):
    # Mock the utility function that reads the pickle file
    mock_network = nx.Graph()
    monkeypatch.setattr(
        utils, 'read_pickle_network', lambda file: mock_network)

    # Setup: Assume network file exists
    resistance_network_fixture.network = 'existing_network_file'

    resistance_network_fixture.read_network()

    # Assertions
    assert resistance_network_fixture._G == mock_network


def test_evaluate_null_model(resistance_network_fixture):
    mock_df = pd.DataFrame({
        'EDGE_ID': [1, 2, 3],
        resistance_network_fixture.length_col: [10, 20, 30]
    })

    # Create a mock result DataFrame to be returned by parsePairwise
    mock_parse_result = pd.DataFrame({
        "aic": [-100],
        "aic_null": [-200],
        "loglik_null": [-300],
        "loglik": [-300],
        "delta_aic_null": [100],
        "r2m": [0.5]
    }, index=[1])

    # Mocking external functions and methods used in evaluate_null_model
    with patch('resistnet.utils.nx_to_df', return_value=mock_df), \
         patch('resistnet.resist_dist.parsePairwise',
               return_value=(None, mock_parse_result)):
        # Call evaluate_null_model
        result = resistance_network_fixture.evaluate_null_model()

    assert isinstance(result, pd.DataFrame)
    assert 'loglik' in result.columns
    assert 'aic' in result.columns
    assert 'r2m' in result.columns
    assert 'delta_aic_null' in result.columns
    assert result.loc['distance_only', 'loglik'] == -300
    assert result.loc['distance_only', 'aic'] == 100


def test_model_output_mock(resistance_network_fixture):
    # Mock model representing an individual in the genetic algorithm
    model = [1, 1, 1, 1, 1]

    # Mock effective resistance matrix
    mock_effective_resistance_matrix = np.array([[0.5, 1.0], [1.0, 0.5]])

    # Mock the effectiveResistanceMatrix function
    with patch(
        'resistnet.resist_dist.effectiveResistanceMatrix',
        return_value=mock_effective_resistance_matrix) \
            as mock_effectiveResistanceMatrix:
        r, multi = resistance_network_fixture.model_output(model)

        mock_effectiveResistanceMatrix.assert_called_once()


def test_model_output(resistance_network_fixture):
    # Mock model representing an individual in the genetic algorithm
    model = [1, 1, 1, 1, 1]

    try:
        r, multi = resistance_network_fixture.model_output(model)
        exception_raised = False
    except Exception as e:
        exception_raised = True
        # Optionally, you can print the exception or take some action
        print(f"Exception occurred: {e}")

    # Assert no exception was raised
    assert not exception_raised, \
        "An exception was raised during the model_output execution"

    # Assert that r is a numpy.ndarray
    assert isinstance(r, np.ndarray), \
        "Expected r to be a numpy.ndarray"

    # Assert that multi is a pandas.Series
    assert isinstance(multi, pd.Series), \
        "Expected multi to be a pandas.Series"


def test_resistance_network_worker_creation(resistance_network_fixture):
    # Prepare additional arguments with defaults or extracted from the fixture
    worker_args = {
        'network': resistance_network_fixture.network,
        'pop_agg': resistance_network_fixture.pop_agg,
        'reachid_col': resistance_network_fixture.reachid_col,
        'length_col': resistance_network_fixture.length_col,
        'variables': resistance_network_fixture.variables,
        'agg_opts': resistance_network_fixture.agg_opts,
        'inc': resistance_network_fixture._inc,
        'point_coords': resistance_network_fixture._point_coords,
        'points_names': resistance_network_fixture._points_names,
        'points_snapped': resistance_network_fixture._points_snapped,
        'points_labels': resistance_network_fixture._points_labels,
        'predictors': resistance_network_fixture._predictors,
        'edge_order': resistance_network_fixture._edge_order,
        'gendist': resistance_network_fixture._gendist,
        # Use defaults or specific values for the additional attributes
        'fitmetric': 'aic',
        'posWeight': False,
        'fixWeight': False,
        'allShapes': False,
        'fixShape': False,
        'min_weight': 0.0,
        'max_shape': 100.0
    }

    # Instantiate the ResistanceNetworkWorker
    worker = ResistanceNetworkWorker(**worker_args)

    # Assertions to verify that the worker is correctly instantiated
    assert worker.pop_agg == resistance_network_fixture.pop_agg
    assert worker.reachid_col == resistance_network_fixture.reachid_col
    assert worker.length_col == resistance_network_fixture.length_col
    assert worker.variables == resistance_network_fixture.variables
    assert worker.agg_opts == resistance_network_fixture.agg_opts

    # Validate unique attributes for ResistanceNetworkWorker
    assert worker.fitmetric == 'aic'
    assert worker.posWeight is False
    assert worker.fixWeight is False
    assert worker.allShapes is False
    assert worker.fixShape is False
    assert worker.min_weight == 0.0
    assert worker.max_shape == 100.0


def test_initialization(sim_resistance_network_params):
    sim_engine = SimResistanceNetwork(**sim_resistance_network_params)
    assert sim_engine.reachid_col == \
        sim_resistance_network_params["reachid_col"]
    assert sim_engine.length_col == \
        sim_resistance_network_params["length_col"]
    assert sim_engine.verbose == \
        sim_resistance_network_params["verbose"]


@pytest.mark.parametrize("num_reps,num_samples", [(1, 10), (2, 5), (3, 3)])
def test_simulate_with_different_parameters(
        sim_resistance_network_params, fake_spec_file, num_reps, num_samples):
    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "simulation_output")

        # Initialize the SimResistanceNetwork with mocked components
        # with patch("resistnet.utils") as mock_utils:
        #     mock_utils.sample_nodes.return_value = {
        #         i: (0.0, 0.0) for i in range(num_samples)
        #     }
        sim_engine = SimResistanceNetwork(**sim_resistance_network_params)

        # Simulate the resistance network
        sim_engine.simulate(
            spec_file=fake_spec_file,
            num_reps=num_reps,
            num_samples=num_samples,
            out=output_prefix
        )

        # Check for the correct number of output files and their contents
        for r in range(1, num_reps + 1):
            resistance_matrix_file = \
                f"{output_prefix}_{r}.ResistanceMatrix.tsv"
            coords_file = f"{output_prefix}_{r}.coords"

            # Check if files exist
            assert os.path.exists(
                resistance_matrix_file), \
                f"{resistance_matrix_file} missing"
            assert os.path.exists(coords_file), f"{coords_file} missing"

            # Check the number of rows in the coords file
            with open(coords_file, 'r') as file:
                num_lines = len(file.readlines())
                expected_lines = num_samples + 1  # +1 for the header
                assert num_lines == expected_lines, \
                    f"{coords_file} should have {expected_lines} \
                        lines but has {num_lines}"
            # Check the dimensions of the resistance matrix
            resistance_matrix_file = \
                f"{output_prefix}_{r}.ResistanceMatrix.tsv"
            df = pd.read_csv(
                resistance_matrix_file, sep='\t', header=0, index_col=0)
            expected_shape = (num_samples, num_samples)
            assert df.shape == expected_shape, \
                f"{resistance_matrix_file} should have shape {expected_shape}"
