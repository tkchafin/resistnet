import os
import pytest
import tempfile
import resistnet
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunnerTPE


@pytest.fixture
def tpe_parameters_fixture():
    return {
        'fixWeight': False,
        'fixShape': True,
        'min_weight': 0.5,
        'max_shape': 2.0,
        'max_hof_size': 100,
        'fitmetric': 'aic',
        'awsum': 0.75,
        'only_keep': True,
        'verbose': False,
        'use_full': False,
        'report_all': False
    }


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


@pytest.fixture
def model_runner_fixture(resistance_network_fixture):
    seed = 1234
    verbose = True

    model_runner = ModelRunnerTPE(
        resistance_network=resistance_network_fixture,
        seed=seed, verbose=verbose)

    yield model_runner


def test_model_runner_initialization(resistance_network_fixture):
    seed = 1234
    verbose = True

    model_runner = ModelRunnerTPE(
        resistance_network=resistance_network_fixture,
        seed=seed, verbose=verbose)

    # Check if the passed parameters are correctly assigned
    assert model_runner.seed == seed
    assert model_runner.resistance_network == resistance_network_fixture
    assert model_runner.verbose == verbose

    # Check if the default values are correctly initialized
    assert isinstance(model_runner.workers, list)
    assert model_runner.bests is None
    assert isinstance(model_runner.logger, list)

    # Check the default values of GA parameters
    assert model_runner.fixWeight is None
    assert model_runner.fixShape is None
    assert model_runner.min_weight is None
    assert model_runner.max_shape is None
    assert model_runner.max_hof_size is None
    assert model_runner.fitmetric is None
    assert model_runner.awsum is None
    assert model_runner.only_keep is None
    assert model_runner.report_all is None


def test_set_tpe_parameters(
        model_runner_fixture, tpe_parameters_fixture):
    # Call set_ga_parameters on the fixture instance with the test parameters
    model_runner_fixture.set_tpe_parameters(
        **tpe_parameters_fixture
    )

    # Assert that each parameter is set correctly
    for param, value in tpe_parameters_fixture.items():
        assert getattr(model_runner_fixture, param) == \
            value, f"Parameter {param} not set correctly"


def test_run_tpe(model_runner_fixture, tpe_parameters_fixture):
    with tempfile.TemporaryDirectory() as temp_dir:
        out_prefix = os.path.join(temp_dir, "tpe_output")
        tpe_parameters_fixture["verbose"] = False
        tpe_parameters_fixture["out"] = out_prefix
        tpe_parameters_fixture["threads"] = 1
        tpe_parameters_fixture["max_evals"] = 1

        # Run the genetic algorithm with the modified parameters
        try:
            model_runner_fixture.run_tpe(**tpe_parameters_fixture)
            error_occurred = False
        except Exception as e:
            error_occurred = True
            print(f"Error during TPE run: {e}")

        # Assert that no errors occurred during the run
        assert not error_occurred, "An error occurred during the TPE run"

        # Assert that the best models are identified (bests is not None)
        assert model_runner_fixture.bests is not None, \
            "ModelRunner.bests missing"
