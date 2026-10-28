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


def test_reverse_inverse_ricker():
    # As Inv.Rev.Ricker in ResistanceGA: the inverse Ricker of the reversed
    # data (it was the inverse of the inverse Ricker, i.e. the Ricker)
    x = pd.Series(np.linspace(0, 10, 51))
    np.testing.assert_allclose(tr.revInvRicker(x, 3, 1),
                               tr.invRicker(10 - x, 3, 1))


def test_transformations_distinct():
    # After rescaling to 0-1, the eight transformations all differ
    x = pd.Series(np.linspace(0, 10, 51))
    curves = []
    for t in range(1, 9):
        y = tr.apply_transformation(x, t, 3)
        curves.append(tr.rescaleCols(y, 0, 1).to_numpy())
    for i in range(8):
        for j in range(i + 1, 8):
            assert np.abs(curves[i] - curves[j]).max() > 0.05, (i + 1, j + 1)
    # No transformation (code 0, or shape 0) returns the data unchanged
    assert tr.apply_transformation(x, 0, 3) is x
    assert tr.apply_transformation(x, 5, 0) is x


def test_numpy_shape_accepted():
    # Shapes read back from a DataFrame row are numpy scalars
    x = pd.Series(np.linspace(0, 10, 11))
    for shape in [np.int64(3), np.float64(3.0)]:
        np.testing.assert_allclose(tr.monomolecular(x, shape, 1),
                                   tr.monomolecular(x, 3, 1))
