import numpy as np

"""
Transformations from Peterman et al. (2018) ResistanceGA package
Full Citation:
Peterman, W. E. 2018. ResistanceGA: An R package for the optimization
of resistance surfaces using genetic algorithms. Methods in Ecology
and Evolution doi:10.1111/2041-210X.12984
GitHub: https://github.com/wpeterman/ResistanceGA
"""


def rescaleCols(df, m, M):
    """
    Rescale all columns in a pandas DataFrame between m and M, where M > m >= 0

    Args:
        df: A pandas DataFrame.
        m: The minimum value of the range to scale to.
        M: The maximum value of the range to scale to.

    Returns:
        pandas.DataFrame: The rescaled DataFrame.
    """
    df -= df.min()
    df /= df.max()
    return (df * (M - m)) + m


def ricker(dat, shape, ceiling):
    """
    Apply the Ricker model transformation to data.

    Args:
        dat: The input data.
        shape: The shape parameter of the Ricker model.
        ceiling: The ceiling parameter of the Ricker model.

    Returns:
        ndarray: The transformed data.
    """
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
    """
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
    """
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
    """
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
    """
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
    """
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
    """
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
    """
    d = rescaleCols((-1 * dat), min(dat), max(dat))
    return monomolecular(d, shape, ceiling)
