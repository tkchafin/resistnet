import os
import pytest
import tempfile
import resistnet
import networkx as nx
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock
from resistnet.resistance_network import ResistanceNetwork
import resistnet.utils as utils


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

    print(result)
    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert 'loglik' in result.columns
    assert 'aic' in result.columns
    assert 'r2m' in result.columns
    assert 'delta_aic_null' in result.columns
    assert result.loc['distance_only', 'loglik'] == -300
    assert result.loc['distance_only', 'aic'] == 100


def test_output_and_plot_model(resistance_network_fixture):
    # Mock data for mat_r and edge_r
    mat_r = np.array([[1, 2], [3, 4]])
    edge_r = pd.DataFrame({'Resistance': [1, 2, 3]}, index=[1, 2, 3])

    # Using a temporary directory for output
    with tempfile.TemporaryDirectory() as tmpdirname:
        oname = f"{tmpdirname}/test_output"

        # Mock the utility functions and plotting methods
        with patch('resistnet.utils.write_edges') as mock_write_edges, \
             patch('resistnet.utils.write_matrix') as mock_write_matrix, \
             patch.object(
                ResistanceNetwork,
                'plot_resistance_network'
                ) as mock_plot_network, \
             patch.object(
                ResistanceNetwork,
                'plot_pairwise_model') as mock_plot_pairwise:

            # Call output_and_plot_model
            resistance_network_fixture.output_and_plot_model(
                oname, mat_r, edge_r, plot=True
            )

            # Assertions for file writing
            mock_write_edges.assert_called_with(
                f"{oname}.ResistanceEdges.tsv", edge_r, edge_r.index
            )
            mock_write_matrix.assert_called_with(
                f"{oname}.ResistanceMatrix.tsv",
                mat_r, list(resistance_network_fixture._points_labels.values())
            )

            # Assertions for plotting
            mock_plot_network.assert_called_with(edge_r, oname)
            if resistance_network_fixture._gendist is not None:
                mock_plot_pairwise.assert_called_with(
                    mat_r, oname, partition=False
                )

            # Call output_and_plot_model with plot=False
            resistance_network_fixture.output_and_plot_model(
                oname, mat_r, edge_r, plot=False
            )

            # Assertions for no plotting
            mock_plot_network.assert_not_called()
            mock_plot_pairwise.assert_not_called()