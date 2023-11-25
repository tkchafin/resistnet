from unittest.mock import patch, MagicMock
import numpy as np
import pandas as pd
import resistnet.MLPE as mlpe


# Assuming MLPE_R is a function within resistnet.MLPE
@patch('resistnet.MLPE.ro.globalenv')
def test_MLPE_R_with_mocked_R(mock_globalenv):
    X_mock = np.random.rand(10, 10)
    Y_mock = np.random.rand(10, 10)

    # Mock the MLPE function in R
    mock_mlpe = MagicMock()
    mock_mlpe.return_value = "Expected Result"
    mock_globalenv.__getitem__.return_value = mock_mlpe

    # Call your function
    result = mlpe.MLPE_R(X_mock, Y_mock, scale=True)

    # Assert that the result is as expected
    assert result == "Expected Result"
    mock_mlpe.assert_called_once()


def test_ZZ_mat_valid_input():
    # Create a mock ID dataframe
    id_data = pd.DataFrame({
        'pop1': [1, 2, 3],
        'pop2': [2, 3, 1]
    })

    # Define number of populations
    pops = 3

    # Call the function
    zz_matrix = mlpe.ZZ_mat_(pops, id_data)
    print(zz_matrix)
    # Assertions to check if the output is as expected
    assert zz_matrix.shape == (pops, id_data.shape[0])
    # Check if the values in the matrix are correctly set
    expected_matrix = np.array([
        [1, 0, 1],
        [1, 1, 0],
        [0, 1, 1]
    ])
    assert np.array_equal(zz_matrix, expected_matrix)


def test_ZZ_mat_empty_dataframe():
    # Test with an empty dataframe
    id_data = pd.DataFrame()
    pops = 3
    zz_matrix = mlpe.ZZ_mat_(pops, id_data)

    # Check if the output is an empty matrix with correct dimensions
    assert zz_matrix.shape == (pops, 0)
    assert np.array_equal(zz_matrix, np.zeros((pops, 0)))
