import os
import io
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
