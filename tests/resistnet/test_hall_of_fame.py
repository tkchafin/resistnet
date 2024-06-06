import os
import io
import tempfile
import contextlib
import pytest
import pandas as pd
import numpy as np
import resistnet
from resistnet.hall_of_fame import HallOfFame


def create_fake_population(num_models, num_variables):
    """Creates a fake population for testing."""
    population = []
    for _ in range(num_models):
        model = {'fitness': np.random.rand()}
        for var in range(num_variables):
            model[f'var{var+1}'] = np.random.rand()
            model[f'var{var+1}_weight'] = 1.0
            model[f'var{var+1}_trans'] = 0
            model[f'var{var+1}_shape'] = 1.0
        model.update(
            {'loglik': np.random.rand(), 'r2m': np.random.rand(),
             'aic': -1 * np.random.rand(), 'delta_aic_null': np.random.rand()}
        )
        population.append(model)
    return population


@pytest.fixture
def hall_of_fame_df_fixture():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(
        base_path, 'data', 'test_ensemble', 'replicates',
        'test_18.out.HallOfFame.tsv'
    )

    # Read the TSV file into a DataFrame
    df = pd.read_csv(file_path, sep='\t', header=0)
    return df


@pytest.fixture
def hall_of_fame_fixture():
    variables = ['var1', 'var2']
    max_size = 10
    num_models = 5

    fake_population = create_fake_population(num_models, len(variables))

    hall_of_fame = HallOfFame(variables, max_size, init_pop=fake_population)
    return hall_of_fame


@pytest.fixture
def hall_of_fame_fixture_complete(hall_of_fame_fixture):
    hall_of_fame_fixture.delta_aic()
    hall_of_fame_fixture.akaike_weights()
    hall_of_fame_fixture.cumulative_akaike(0.95)
    hall_of_fame_fixture.model_average_weights(ignore_keep=False)
    hall_of_fame_fixture.relative_variable_importance(ignore_keep=False)
    return hall_of_fame_fixture


def test_hall_of_fame_initialization():
    variables = ['var1', 'var2']
    max_size = 10
    init_pop = None

    # Initialize the HallOfFame
    hall_of_fame = HallOfFame(variables, max_size, init_pop)

    # Check the DataFrame structure
    expected_columns = ["fitness"]
    for v in variables:
        expected_columns.extend(
            [str(v), f"{v}_weight", f"{v}_trans", f"{v}_shape"]
        )
    expected_columns.extend(["loglik", "r2m", "aic", "delta_aic_null"])

    assert list(hall_of_fame.data.columns) == expected_columns, \
        "Data columns do not match expected columns"

    # Check other attributes
    assert hall_of_fame.variables == variables, \
        "Variables attribute not set correctly"
    assert hall_of_fame.max_size == max_size, \
        "Max size attribute not set correctly"
    assert hall_of_fame.min_fitness == float("-inf"), \
        "Min fitness attribute not set correctly"
    assert hall_of_fame.rvi is None, \
        "RVI attribute should be None initially"
    assert hall_of_fame.maw is None, \
        "MAW attribute should be None initially"
    assert hall_of_fame.best is None, \
        "Best attribute should be None initially"
    assert hall_of_fame.zero_threshold == 1e-17, \
        "Zero threshold attribute not set correctly"


def test_hall_of_fame_from_dataframe(hall_of_fame_df_fixture):
    max_size = 200  # Example max_size, adjust as needed

    # Create a HallOfFame instance from the DataFrame
    hall_of_fame_instance = HallOfFame.from_dataframe(hall_of_fame_df_fixture,
                                                      max_size)

    # Asserts to verify the HallOfFame instance
    assert len(hall_of_fame_instance.data) == len(hall_of_fame_df_fixture), \
        "Data length in HallOfFame instance does not match input DataFrame"


def test_check_population():
    variables = ['var1', 'var2']
    max_size = 10
    num_models = 5

    # Create a fake population
    fake_population = create_fake_population(num_models, len(variables))

    # Initialize HallOfFame
    hall_of_fame = HallOfFame(variables, max_size)

    # Check the population
    hall_of_fame.check_population(fake_population)

    # Assertions to verify the HallOfFame data is updated correctly
    assert len(hall_of_fame.data) <= max_size, \
        "Hall of fame size exceeds maximum size"
    assert all(
        col in hall_of_fame.data.columns for col in ['fitness', 'var1', 'var2']
    ), "Hall of fame data does not contain expected columns"


def test_print_hof(hall_of_fame_fixture):
    # Redirect the stdout to a string buffer
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        hall_of_fame_fixture.printHOF()

    # Get the content from the buffer
    output = buffer.getvalue()

    # Asserts to verify the output
    assert "fitness" in output, "Output should contain 'fitness'"
    assert "var1" in output, "Output should contain 'var1'"
    assert "var2" in output, "Output should contain 'var2'"


def test_calculate_bic(hall_of_fame_fixture):
    n = 10

    hall_of_fame_fixture.calculate_bic(n)

    assert "bic" in hall_of_fame_fixture.data.columns, \
        "BIC column not added to data"


def test_delta_aic(hall_of_fame_fixture):
    hall_of_fame_fixture.delta_aic()

    assert "delta_aic_best" in hall_of_fame_fixture.data.columns, \
        "delta_aic_best column not added to data"


def test_akaike_weights(hall_of_fame_fixture):
    hall_of_fame_fixture.akaike_weights()

    assert "akaike_weight" in hall_of_fame_fixture.data.columns, \
        "akaike_weight column not added to data"


def test_cumulative_akaike(hall_of_fame_fixture):
    threshold = 0.8
    hall_of_fame_fixture.cumulative_akaike(threshold)

    assert "acc_akaike_weight" in hall_of_fame_fixture.data.columns, \
        "acc_akaike_weight column not added to data"
    assert "keep" in hall_of_fame_fixture.data.columns, \
        "keep column not added to data"


def test_relative_variable_importance(hall_of_fame_fixture):
    hall_of_fame_fixture.delta_aic()
    hall_of_fame_fixture.akaike_weights()
    hall_of_fame_fixture.cumulative_akaike(0.95)
    hall_of_fame_fixture.relative_variable_importance(ignore_keep=False)

    assert not hall_of_fame_fixture.rvi.empty, \
        "RVI DataFrame should not be empty"
    assert "variable" in hall_of_fame_fixture.rvi.columns, \
        "RVI DataFrame should contain 'variable' column"
    assert "RVI" in hall_of_fame_fixture.rvi.columns, \
        "RVI DataFrame should contain 'RVI' column"


def test_get_best_model(hall_of_fame_fixture):
    best_model_df = hall_of_fame_fixture.get_best_model()

    assert best_model_df is not None, "Best model DataFrame should not be None"
    assert not best_model_df.empty, "Best model DataFrame should not be empty"
    # Check if the returned DataFrame contains expected columns
    for var in hall_of_fame_fixture.variables:
        assert any(best_model_df["Variable"] == var), \
            f"Best model should contain variable {var}"


def test_model_average_weights(hall_of_fame_fixture):
    hall_of_fame_fixture.delta_aic()
    hall_of_fame_fixture.akaike_weights()
    hall_of_fame_fixture.cumulative_akaike(0.95)
    hall_of_fame_fixture.model_average_weights(ignore_keep=False)

    assert not hall_of_fame_fixture.maw.empty, \
        "MAW DataFrame should not be empty"
    assert "variable" in hall_of_fame_fixture.maw.columns, \
        "MAW DataFrame should contain 'variable' column"
    assert "MAW" in hall_of_fame_fixture.maw.columns, \
        "MAW DataFrame should contain 'MAW' column"


def test_get_variables(hall_of_fame_fixture):
    variables = hall_of_fame_fixture.get_variables()

    assert isinstance(variables, list), "get_variables should return a list"
    assert variables == hall_of_fame_fixture.variables, \
        "Returned variables do not match the expected list"


def test_get_best(hall_of_fame_fixture):
    best_model = hall_of_fame_fixture.getBest()

    assert isinstance(best_model, pd.DataFrame), \
        "getBest should return a DataFrame"
    assert not best_model.empty, "Best model DataFrame should not be empty"


def test_get_rvi(hall_of_fame_fixture):
    hall_of_fame_fixture.delta_aic()
    hall_of_fame_fixture.akaike_weights()
    hall_of_fame_fixture.cumulative_akaike(0.95)
    hall_of_fame_fixture.relative_variable_importance(ignore_keep=False)
    rvi_data = hall_of_fame_fixture.getRVI()

    assert isinstance(rvi_data, pd.DataFrame), \
        "getRVI should return a DataFrame"
    assert not rvi_data.empty, "RVI DataFrame should not be empty"


def test_get_maw(hall_of_fame_fixture):
    hall_of_fame_fixture.delta_aic()
    hall_of_fame_fixture.akaike_weights()
    hall_of_fame_fixture.cumulative_akaike(0.95)
    hall_of_fame_fixture.model_average_weights(ignore_keep=False)
    maw_data = hall_of_fame_fixture.getMAW()

    assert isinstance(maw_data, pd.DataFrame), \
        "getMAW should return a DataFrame"
    assert not maw_data.empty, "MAW DataFrame should not be empty"


def test_get_hof(hall_of_fame_fixture):
    hof_data = hall_of_fame_fixture.getHOF(only_keep=False)

    assert isinstance(hof_data, pd.DataFrame), \
        "getHOF should return a DataFrame"
    assert not hof_data.empty, "Hall of Fame DataFrame should not be empty"


def test_plot_icprofile(hall_of_fame_fixture_complete):
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_ICprofile")
        hall_of_fame_fixture_complete.plot_ICprofile(oname=filename)

        assert os.path.isfile(f"{filename}.ICprofile.pdf"), \
            "IC profile plot file was not created"
        assert os.path.getsize(f"{filename}.ICprofile.pdf") > 0, \
            "IC profile plot file is empty"


def test_plot_metric_pw(hall_of_fame_fixture_complete):
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_metricPW")
        hall_of_fame_fixture_complete.plotMetricPW(oname=filename)

        assert os.path.isfile(f"{filename}.pairPlot.pdf"), \
            "Metric PW plot file was not created"
        assert os.path.getsize(f"{filename}.pairPlot.pdf") > 0, \
            "Metric PW plot file is empty"


def test_plot_variable_importance(hall_of_fame_fixture_complete):
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_varImportance")
        hall_of_fame_fixture_complete.plotVariableImportance(oname=filename)

        assert os.path.isfile(f"{filename}.varImportance.pdf"), \
            "Variable importance plot file was not created"
        assert os.path.getsize(f"{filename}.varImportance.pdf") > 0, \
            "Variable importance plot file is empty"


def test_plot_model_averaged_weights(hall_of_fame_fixture_complete):
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_modavgWeights")
        hall_of_fame_fixture_complete.plotModelAveragedWeights(oname=filename)

        assert os.path.isfile(f"{filename}.modavgWeights.pdf"), \
            "Model-averaged weights plot file was not created"
        assert os.path.getsize(f"{filename}.modavgWeights.pdf") > 0, \
            "Model-averaged weights plot file is empty"


def test_write_model_summary(hall_of_fame_fixture_complete):
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_HallOfFame")
        hall_of_fame_fixture_complete.writeModelSummary(filename)

        assert os.path.isfile(f"{filename}.HallOfFame.tsv"), \
            "Model summary TSV file was not created"
        assert os.path.getsize(f"{filename}.HallOfFame.tsv") > 0, \
            "Model summary TSV file is empty"


def test_write_maw(hall_of_fame_fixture_complete):
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_modavgWeights")
        hall_of_fame_fixture_complete.writeMAW(filename)

        assert os.path.isfile(f"{filename}.modavgWeights.tsv"), \
            "MAW TSV file was not created"
        assert os.path.getsize(f"{filename}.modavgWeights.tsv") > 0, \
            "MAW TSV file is empty"


def test_write_maw_empty(hall_of_fame_fixture_complete):
    hall_of_fame_fixture_complete.maw = None
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_modavgWeights")
        hall_of_fame_fixture_complete.writeMAW(filename)

        assert os.path.isfile(f"{filename}.modavgWeights.tsv"), \
            "MAW TSV file was not created"
        assert os.path.getsize(f"{filename}.modavgWeights.tsv") > 0, \
            "MAW TSV file is empty"


def test_write_rvi(hall_of_fame_fixture_complete):
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_varImportance")
        hall_of_fame_fixture_complete.writeRVI(filename)

        assert os.path.isfile(f"{filename}.varImportance.tsv"), \
            "RVI TSV file was not created"
        assert os.path.getsize(f"{filename}.varImportance.tsv") > 0, \
            "RVI TSV file is empty"


def test_write_rvi_empty(hall_of_fame_fixture_complete):
    hall_of_fame_fixture_complete.rvi = None
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_varImportance")
        hall_of_fame_fixture_complete.writeRVI(filename)

        assert os.path.isfile(f"{filename}.varImportance.tsv"), \
            "RVI TSV file was not created"
        assert os.path.getsize(f"{filename}.varImportance.tsv") > 0, \
            "RVI TSV file is empty"


def test_write_best(hall_of_fame_fixture_complete):
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, "test_bestModel")
        hall_of_fame_fixture_complete.writeBest(filename)

        assert os.path.isfile(f"{filename}.bestModel.tsv"), \
            "Best model TSV file was not created"
        assert os.path.getsize(f"{filename}.bestModel.tsv") > 0, \
            "Best model TSV file is empty"
