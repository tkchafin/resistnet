import pytest
import numpy as np
import pandas as pd
import resistnet.transform as tr


@pytest.fixture
def example_data():
    dat = pd.Series(np.random.uniform(0, 5, size=(10)))
    return dat


@pytest.mark.parametrize("data, m, M, expected_exception", [
    # Test with valid DataFrame and range
    (pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}), 0, 10, None),
    # Test with valid Series and range
    (pd.Series([1, 2, 3]), 0, 10, None),
    # Test with invalid range (M <= m)
    (pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}), 5, 5, ValueError),
    # Test with non-DataFrame input
    ([1, 2, 3], 0, 10, TypeError),
    # Test with empty DataFrame
    (pd.DataFrame(), 0, 10, ValueError)
])
def test_rescaleCols(data, m, M, expected_exception):
    if expected_exception:
        with pytest.raises(expected_exception):
            tr.rescaleCols(data, m, M)
    else:
        rescaled_df = tr.rescaleCols(data, m, M)
        assert rescaled_df.min().min() >= m
        assert rescaled_df.max().max() <= M
        assert rescaled_df.shape == data.shape


# Parameterized test for transformation functions
@pytest.mark.parametrize("function", [
    tr.ricker,
    tr.invRicker,
    tr.revInvRicker,
    tr.revRicker,
    tr.monomolecular,
    tr.invMonomolecular,
    tr.revInvMonomolecular,
    tr.revMonomolecular
])
@pytest.mark.parametrize("shape, ceiling, expected_exception", [
    (2, 10, None),  # valid inputs
    (0, 10, ValueError),  # invalid shape (zero)
    ('invalid', 10, ValueError),  # invalid type for shape
])
def test_transformations(function, shape, ceiling, expected_exception,
                         example_data):
    if expected_exception:
        with pytest.raises(expected_exception):
            function(example_data, shape, ceiling)
    else:
        transformed_data = function(example_data, shape, ceiling)
        assert isinstance(transformed_data, (np.ndarray, pd.Series))
