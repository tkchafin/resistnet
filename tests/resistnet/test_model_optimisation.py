import os
import pytest
import tempfile
import resistnet
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunner


@pytest.fixture
def ga_parameters_fixture():
    return {
        'mutpb': 0.1,
        'cxpb': 0.8,
        'indpb': 0.05,
        'popsize': 10,
        'maxpopsize': 10,
        'fixWeight': False,
        'fixShape': True,
        'min_weight': 0.5,
        'max_shape': 2.0,
        'max_hof_size': 100,
        'tournsize': 3,
        'fitmetric': 'aic',
        'awsum': 0.75,
        'only_keep': True,
        'verbose': False,
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

    model_runner = ModelRunner(
        resistance_network=resistance_network_fixture,
        seed=seed, verbose=verbose)

    yield model_runner


def test_model_runner_initialization(resistance_network_fixture):
    seed = 1234
    verbose = True

    model_runner = ModelRunner(
        resistance_network=resistance_network_fixture,
        seed=seed, verbose=verbose)

    # Check if the passed parameters are correctly assigned
    assert model_runner.seed == seed
    assert model_runner.resistance_network == resistance_network_fixture
    assert model_runner.verbose == verbose

    # Check if the default values are correctly initialized
    assert isinstance(model_runner.workers, list)
    assert model_runner.bests is None
    assert model_runner.toolbox is None
    assert isinstance(model_runner.logger, list)

    # Check the default values of GA parameters
    assert model_runner.cxpb is None
    assert model_runner.mutpb is None
    assert model_runner.indpb is None
    assert model_runner.popsize is None
    assert model_runner.maxpopsize is None
    assert model_runner.fixWeight is None
    assert model_runner.fixShape is None
    assert model_runner.min_weight is None
    assert model_runner.max_shape is None
    assert model_runner.max_hof_size is None
    assert model_runner.tournsize is None
    assert model_runner.fitmetric is None
    assert model_runner.awsum is None
    assert model_runner.only_keep is None
    assert model_runner.report_all is None


def test_set_ga_parameters(
        model_runner_fixture, ga_parameters_fixture):
    # Call set_ga_parameters on the fixture instance with the test parameters
    model_runner_fixture.set_ga_parameters(
        **ga_parameters_fixture
    )

    # Assert that each parameter is set correctly
    for param, value in ga_parameters_fixture.items():
        assert getattr(model_runner_fixture, param) == \
            value, f"Parameter {param} not set correctly"


def test_run_ga(model_runner_fixture, ga_parameters_fixture):
    with tempfile.TemporaryDirectory() as temp_dir:
        out_prefix = os.path.join(temp_dir, "ga_output")
        ga_parameters_fixture["verbose"] = False
        ga_parameters_fixture["out"] = out_prefix
        ga_parameters_fixture["threads"] = 1
        ga_parameters_fixture["maxgens"] = 1

        # Run the genetic algorithm with the modified parameters
        try:
            model_runner_fixture.run_ga(**ga_parameters_fixture)
            error_occurred = False
        except Exception as e:
            error_occurred = True
            print(f"Error during GA run: {e}")

        # Assert that no errors occurred during the run
        assert not error_occurred, "An error occurred during the GA run"

        # Assert that the best models are identified (bests is not None)
        assert model_runner_fixture.bests is not None, \
            "ModelRunner.bests missing"

        # Assert that the expected output files are generated
        expected_files = [
            f"{out_prefix}.varImportance.tsv",
            f"{out_prefix}.HallOfFame.tsv",
            f"{out_prefix}.FitnessLog.tsv",
            f"{out_prefix}.Model-Average.streamsByResistance.pdf"
        ]

        for file in expected_files:
            assert os.path.isfile(file), f"Expected output file not found: \
                {file}"

