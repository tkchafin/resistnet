import numpy as np
import pandas as pd

"""
Transformations from Peterman et al. (2018) ResistanceGA package
Full Citation:
Peterman, W. E. 2018. ResistanceGA: An R package for the optimization
of resistance surfaces using genetic algorithms. Methods in Ecology
and Evolution doi:10.1111/2041-210X.12984
GitHub: https://github.com/wpeterman/ResistanceGA
"""


def rescaleCols(df, m, M):
    if not isinstance(df, (pd.DataFrame, pd.Series)):
        raise TypeError("Input 'df' must be a pandas DataFrame or Series.")

    if not (isinstance(m, (int, float)) and isinstance(M, (int, float))):
        raise ValueError("'m' and 'M' must be numeric values.")

    if M <= m:
        raise ValueError(
            "The maximum value M must be greater than the minimum value m."
            )

    # Convert Series to DataFrame for uniform processing
    series = False
    if isinstance(df, pd.Series):
        series = True
        df = df.to_frame()

    # Check for empty input
    if df.empty or df.size == 0:
        raise ValueError("Input is empty.")

    # Perform rescaling
    rescaled_df = df.copy()
    rescaled_df -= rescaled_df.min()
    rescaled_df /= rescaled_df.max()
    rescaled_df = (rescaled_df * (M - m)) + m

    # Convert back to Series if the original input was a Series
    if series:
        return rescaled_df.iloc[:, 0]
    return rescaled_df


def ricker(dat, shape, ceiling):
    """
    Apply the Ricker model transformation to data.

    This function transforms the input data using the Ricker model. Before
    applying the transformation, it validates the inputs to ensure they meet
    the expected criteria.

    Args:
        dat: The input data, expected to be a numeric array or similar
             iterable.
        shape: The shape parameter of the Ricker model, must be a numeric value
               and non-zero.
        ceiling: The ceiling parameter of the Ricker model, must be a numeric
                 value.

    Returns:
        ndarray: The transformed data.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    validate_arguments(dat, shape, ceiling)
    return ceiling * dat * np.exp(-1 * dat / shape) + 1


def invRicker(dat, shape, ceiling):
    """
    Apply the inverse Ricker model transformation to data.

    Args:
        dat: The input data.
        shape: The shape parameter of the Ricker model.
        ceiling: The ceiling parameter of the Ricker model.

    Returns:
        ndarray: The transformed data.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    validate_arguments(dat, shape, ceiling)
    return (-1 * ceiling) * dat * np.exp(-1 * dat / shape) - 1


def revInvRicker(dat, shape, ceiling):
    """
    Apply the reverse inverse Ricker model transformation to data.

    Args:
        dat: The input data.
        shape: The shape parameter of the Ricker model.
        ceiling: The ceiling parameter of the Ricker model.

    Returns:
        pandas.DataFrame: The transformed and rescaled DataFrame.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    validate_arguments(dat, shape, ceiling)
    d = invRicker(dat, shape, ceiling)
    return rescaleCols((-1 * d), min(d), max(d))


def revRicker(dat, shape, ceiling):
    """
    Apply the reverse Ricker model transformation to data.

    Args:
        dat: The input data.
        shape: The shape parameter of the Ricker model.
        ceiling: The ceiling parameter of the Ricker model.

    Returns:
        ndarray: The transformed data.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    validate_arguments(dat, shape, ceiling)
    d = rescaleCols((-1 * dat), min(dat), max(dat))
    return ricker(d, shape, ceiling)


def monomolecular(dat, shape, ceiling):
    """
    Apply the monomolecular model transformation to data.

    Args:
        dat: The input data.
        shape: The shape parameter of the monomolecular model.
        ceiling: The ceiling parameter of the monomolecular model.

    Returns:
        ndarray: The transformed data.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    validate_arguments(dat, shape, ceiling)
    return ceiling * (1 - np.exp(-1 * dat / shape)) + 1


def invMonomolecular(dat, shape, ceiling):
    """
    Apply the inverse monomolecular model transformation to data.

    Args:
        dat: The input data.
        shape: The shape parameter of the monomolecular model.
        ceiling: The ceiling parameter of the monomolecular model.

    Returns:
        ndarray: The transformed data.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    validate_arguments(dat, shape, ceiling)
    d = ceiling * np.exp(-1 * dat / shape)
    return (d - min(d)) + 1


def revInvMonomolecular(dat, shape, ceiling):
    """
    Apply the reverse inverse monomolecular model transformation to data.

    Args:
        dat: The input data.
        shape: The shape parameter of the monomolecular model.
        ceiling: The ceiling parameter of the monomolecular model.

    Returns:
        pandas.DataFrame: The transformed and rescaled DataFrame.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    validate_arguments(dat, shape, ceiling)
    d = rescaleCols((-1 * dat), min(dat), max(dat))
    return invMonomolecular(d, shape, ceiling)


def revMonomolecular(dat, shape, ceiling):
    """
    Apply the reverse monomolecular model transformation to data.

    Args:
        dat: The input data.
        shape: The shape parameter of the monomolecular model.
        ceiling: The ceiling parameter of the monomolecular model.

    Returns:
        ndarray: The transformed data.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    validate_arguments(dat, shape, ceiling)
    d = rescaleCols((-1 * dat), min(dat), max(dat))
    return monomolecular(d, shape, ceiling)


def validate_arguments(dat, shape, ceiling):
    """
    Validate the input arguments for transformation functions.

    Args:
        dat: The input data, expected to be a numeric array or similar
             iterable.
        shape: A numeric value representing a parameter of the transformation.
        ceiling: A numeric value representing another parameter of the
                 transformation.

    Raises:
        TypeError: If 'dat' is not a numeric array.
        ValueError: If 'shape' or 'ceiling' is not a numeric value, or if
                    'shape' is zero.
    """
    if not isinstance(dat, (pd.Series)):
        raise TypeError(
            "Input 'dat' must be a pd.Series"
        )

    if not (isinstance(shape, (int, float)) and
            isinstance(ceiling, (int, float))):
        raise ValueError("'shape' and 'ceiling' must be numeric values.")

    if shape == 0:
        raise ValueError("'shape' must not be zero.")
