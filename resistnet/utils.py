import itertools
import math
import pickle
import random

import networkx as nx
import numpy as np
import pandas as pd
from sortedcontainers import SortedDict


def snap_to_node(graph, pos):
    """
    Find the closest node in a graph to the given coordinates.

    Args:
        graph: A graph object containing nodes.
        pos: A tuple of (x, y) coordinates.

    Returns:
        A tuple representing the coordinates of the closest node in the graph.

    Raises:
        ValueError: If 'graph' or 'pos' is not in the expected format or empty.
    """
    if not graph.nodes() or not pos:
        raise ValueError("Invalid graph or position data.")

    nodes = np.array(graph.nodes())
    node_pos = np.argmin(np.sum((nodes - pos) ** 2, axis=1))
    return tuple(nodes[node_pos])


def check_dataframe_columns(df):
    """
    Check if any columns in the dataframe contain only NaN values or a single
    unique value.

    This function iterates through each column in the dataframe and checks for
    two conditions:

    - If the column contains only NaN values.
    - If the column contains only a single unique value (other than NaN).

    Args:
        df: A pandas DataFrame.

    Returns:
        bool: True if all columns have more than one unique value
              (other than NaN), otherwise False.

    Raises:
        ValueError: If 'df' is empty or not a DataFrame.
    """
    if df.empty or not isinstance(df, pd.DataFrame):
        raise ValueError("Invalid DataFrame.")

    for col in df.columns:
        unique_values = df[col].dropna().unique()

        if len(unique_values) == 0:
            print(f"Column '{col}' has all NaN values.")
            return False
        elif len(unique_values) == 1:
            print(f"Column '{col}' has all the same value: {unique_values[0]}")
            return False

    return True


def scrub_bad_columns(df, variables, verbose=False):
    """
    Remove columns from a DataFrame that contain only NaN values or a single
    unique value.

    This function iterates through each column in the DataFrame, identifying
    and removing

    any that meet the following criteria:
    - Contains only NaN values.
    - Contains only a single unique value.

    Args:
        df: A pandas DataFrame to be scrubbed.
        variables: A list of variables (columns) that are of interest.
        verbose: A boolean flag to enable printing of information about removed
                 columns.

    Returns:
        tuple: A tuple containing the cleaned DataFrame and a list of remaining
               variables.

    Raises:
        ValueError: If 'df' or 'variables' is not in the expected format or
                    empty.
    """
    if df.empty or not isinstance(df, pd.DataFrame):
        raise ValueError("Invalid DataFrame.")
    if not variables or not isinstance(variables, list):
        raise ValueError("Invalid variables list.")

    bad_columns = []

    # Iterate through columns to identify bad columns
    for col in df.columns:
        uniq_vals = df[col].dropna().unique()

        if (len(uniq_vals) == 0 or
                (len(uniq_vals) == 1 and np.isnan(uniq_vals[0]))):
            if verbose:
                print(f"Column '{col}' has all NaN values.")
            bad_columns.append(col)
        elif len(uniq_vals) == 1:
            if verbose:
                print(
                    f"Column '{col}' has all the same value:" f"{uniq_vals[0]}"
                    )
            bad_columns.append(col)

    # Remove bad columns from the DataFrame and the variables list
    cleaned_df = df.drop(bad_columns, axis=1)
    cleaned_variables = [v for v in variables if v not in bad_columns]

    return cleaned_df, cleaned_variables


def read_pickle_network(pik):
    """
    Read a network graph from a pickle file and convert it to an undirected
    graph.

    Args:
        pik: The path to the pickle file containing the network graph.

    Returns:
        A networkx.Graph object representing the undirected graph.

    Raises:
        FileNotFoundError: If the pickle file is not found.
        IOError: If there's an error reading from the file.
    """
    try:
        with open(pik, 'rb') as f:
            _G = nx.Graph(pickle.load(f)).to_undirected()
        return _G
    except FileNotFoundError:
        raise FileNotFoundError(f"Pickle file {pik} not found.")
    except IOError:
        raise IOError("Error occurred while reading from the file.")


def find_pair(lst, x, y):
    """
    Check if two elements are consecutive in a list, irrespective of their
    order.

    Args:
        lst: The list in which to check for consecutive elements.
        x: The first element to check.
        y: The second element to check.

    Returns:
        bool: True if x and y are consecutive in lst, False otherwise.
    """
    if x not in lst or y not in lst:
        return False
    elif abs(lst.index(x) - lst.index(y)) == 1:
        return True
    else:
        return False


def nx_to_df(G):
    """
    Convert a networkx Graph object to a pandas DataFrame.

    The function iterates over the edges of the graph, extracting edge data and
    the connected nodes, then stores this information in a DataFrame.

    Args:
        G: A networkx Graph object.

    Returns:
        A pandas DataFrame representing the edges of the graph.

    Raises:
        TypeError: If 'G' is not a networkx Graph object.
    """
    if not isinstance(G, nx.Graph):
        raise TypeError("The input must be a networkx Graph object.")

    edges_data = []
    for p1, p2, e in G.edges(data=True):
        e["left"] = p1
        e["right"] = p2
        edges_data.append(e)
    return pd.DataFrame(edges_data)


def unique(sequence):
    """
    Return a list of unique elements from a sequence, preserving the order.

    Args:
        sequence: A sequence (like a list or tuple) from which unique elements
                  are extracted.

    Returns:
        list: A list containing unique elements of the given sequence.
    """
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]


def read_points_table(pfile):
    """
    Read a tab-separated values (TSV) file into a pandas DataFrame.

    The function checks for the presence of required columns
    ('sample', 'lat', 'long').

    If any required columns are missing, a ValueError is raised.

    Args:
        pfile: Path to the TSV file.

    Returns:
        pandas.DataFrame: A DataFrame containing the data from the file.

    Raises:
        ValueError: If any required columns are missing in the file.
    """
    points = pd.read_csv(pfile, sep="\t", header=0)

    required_columns = {'sample', 'lat', 'long'}
    if not required_columns.issubset(points.columns):
        miss = required_columns - set(points.columns)
        error_message = (
            f"Missing required columns in the points data: {', '.join(miss)}"
        )
        raise ValueError(error_message)

    return points


def read_points_flatten(pfile):
    """
    Read a points file, flatten the data, and store it in a SortedDict.

    The function reads a file into a DataFrame, then iterates over its rows to
    create a SortedDict. The keys are tuples of ('long', 'lat'), and the values
    are lists of 'sample' names.

    Args:
        pfile: Path to the file containing points data.

    Returns:
        SortedDict: A dictionary with coordinates as keys and lists of sample
                    names as values.
    """
    points = read_points_table(pfile)  # Assuming this returns a DataFrame

    flattened_data = SortedDict()
    for _, row in points.iterrows():
        name = row['sample']
        lat = float(row['lat'])
        long = float(row['long'])
        coords = (long, lat)

        if coords not in flattened_data:
            flattened_data[coords] = []
        flattened_data[coords].append(name)

    return flattened_data


def nCr(n, k):
    """
    Calculate the number of combinations for choosing k items from n items.

    Utilizes the mathematical formula for combinations, n choose k, which is
    defined as n! / (k!(n-k)!).

    Args:
        n: The total number of items.
        k: The number of items to choose.

    Returns:
        int: The number of possible combinations.
    """
    f = math.factorial
    return f(n) // f(k) // f(n - k)


def write_numpy_matrix(mat, ofh):
    """
    Write a NumPy matrix to a file with specified formatting.

    This function saves a NumPy matrix to a file, using tab-delimited format
    and suppressing scientific notation.

    Args:
        mat: A NumPy matrix to be written to the file.
        ofh: File handle or file path to which the matrix will be written.

    Raises:
        IOError: If there is an error in writing to the file.
    """
    try:
        with np.printoptions(precision=0, suppress=True):
            np.savetxt(ofh, mat, delimiter="\t")
    except IOError:
        raise IOError("Error occurred while writing the matrix to the file.")


def path_edge_attributes(graph, path, attribute):
    """
    Extract the values of a specified attribute for each edge in a path on a
    graph.

    This function takes a path within a graph and returns the values of a
    specified attribute for each edge along the path.

    Args:
        graph: A networkx Graph object.
        path: A list of nodes representing a path in the graph.
        attribute: The edge attribute whose values are to be returned.

    Returns:
        list: A list containing the values of the specified attribute for each
              edge in the path.
    """
    return [graph[u][v][attribute] for (u, v) in zip(path, path[1:])]


def write_edges(out, edge, ids):
    """
    Write edge data to a CSV file.

    This function creates a DataFrame from edge IDs and their corresponding
    resistance values, and writes this data to a CSV file.

    Args:
        out: The file path where the CSV will be written.
        edge: A list of edge resistance values.
        ids: A list of edge IDs.

    Returns:
        str: The output file path.
    """
    df = pd.DataFrame(list(zip(ids, edge)), columns=["EDGE_ID", "Resistance"])
    df.to_csv(out, sep="\t", header=True, index=False)
    return out


def write_matrix(out, mat, ids):
    """
    Write a matrix to a CSV file.

    This function creates a DataFrame from a matrix and writes it to a CSV file
    using the provided ids for both row and column headers.

    Args:
        out: The file path where the CSV will be written.
        mat: The matrix to be written.
        ids: A list of identifiers for the matrix rows and columns.

    Returns:
        str: The output file path.
    """
    df = pd.DataFrame(mat, columns=ids, index=ids)
    df.to_csv(out, sep="\t", header=True, index=True)
    return out


def get_lower_tri(mat):
    """
    Extract the lower triangular part of a matrix.

    Args:
        mat: A NumPy matrix.

    Returns:
        ndarray: An array containing the elements of the lower triangular part
                 of the matrix.
    """
    n = mat.shape[0]
    i = np.tril_indices(n, -1)
    return mat[i]


def to_from_(pops):
    """
    Generate a DataFrame of combinations of population pairs.

    This function creates a DataFrame containing all combinations of population
    pairs from a given number of populations.

    Args:
        pops: The total number of populations.

    Returns:
        pandas.DataFrame: A DataFrame with two columns 'pop1' and 'pop2'
                          representing population pairs.
    """
    to, frm = list(), list()
    for ia, ib in itertools.combinations(range(1, pops + 1), 2):
        to.append(ia)
        frm.append(ib)
    to[-1], frm[-1] = frm[-1], to[-1]  # Swap the last elements
    return pd.DataFrame({"pop1": to, "pop2": frm})


def sample_nodes(G, samples):
    """
    Randomly sample a specified number of nodes from a graph.

    Args:
        G: A networkx Graph object.
        samples: The number of nodes to sample.

    Returns:
        SortedDict: A dictionary with nodes as keys and a sequential integer as
                    values.
    """
    ret = SortedDict()
    nodes = random.sample(list(G.nodes), samples)
    for i, s in enumerate(nodes, 1):
        ret[s] = i
    return ret
