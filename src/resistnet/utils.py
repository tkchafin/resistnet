import contextlib
import io
import itertools
import random
from math import radians, cos, sin, asin, sqrt
import networkx as nx
import numpy as np
import pandas as pd
from sortedcontainers import SortedDict
import warnings
from autostreamtree import network as astnet
# Shared with autoStreamTree
from autostreamtree.network import find_pair, nCr  # noqa: F401


def graph_to_dag(K, origin):
    """
    Convert an undirected graph to a Directed Acyclic Graph (DAG) where all
    edges are directed toward the tips from the root/ source node.

    Args:
        K (networkx.Graph): The undirected graph to be converted.
        origin (node): The node in the graph that will act as the source point
                       or root.

    Returns:
        networkx.DiGraph: A directed version of the input graph with all paths
                          oriented away from the specified root node.
    """
    # Create an empty Directed Graph
    D = nx.DiGraph()

    # Start BFS from the root node, adding edges to D
    for parent, child in nx.bfs_edges(K, source=origin):
        D.add_edge(parent, child)  # direction is root -> tips
        # copy attributes
        if K.has_edge(parent, child):
            for key, value in K[parent][child].items():
                D[child][parent][key] = value
    return D


def graph_to_dag_converging(K, origin):
    """
    Convert an undirected graph to a Directed Acyclic Graph (DAG) where all
    edges are directed to converge on a specified root node.

    Args:
        K (networkx.Graph): The undirected graph to be converted.
        origin (node): The node in the graph that will act as the convergence
                       point or root.

    Returns:
        networkx.DiGraph: A directed version of the input graph with all paths
                          oriented to converge on the specified root node.
    """
    # Create an empty Directed Graph
    D = nx.DiGraph()

    # Start BFS from the root node, adding edges w direction towards the root
    for parent, child in nx.bfs_edges(K, source=origin):
        D.add_edge(child, parent)  # Reverse direction

        # Copy attributes to the new directed graph
        if K.has_edge(parent, child):
            for key, value in K[parent][child].items():
                D[child][parent][key] = value
    return D


def minmax(X):
    """
    Perform min-max scaling on a NumPy array, scaling all values to be between
    0 and 1.

    Args:
        X (numpy.ndarray): The input array to be scaled.

    Returns:
        numpy.ndarray: The scaled array, with all values between 0 and 1.
    """
    X_min = X.min()
    X_max = X.max()
    X_scaled = (X - X_min) / (X_max - X_min)
    return X_scaled


def minmax_nonzero(X):
    """
    Perform min-max scaling on a NumPy array, initially scaling all values to
    be between 0 and 1, and then adjusting the scale so that the smallest
    non-zero value becomes 1 and other values are adjusted accordingly.

    Args:
        X (numpy.ndarray): The input array to be scaled.

    Returns:
        numpy.ndarray: The scaled array, with values adjusted as described.
    """
    X_min = X.min()
    X_max = X.max()
    X_scaled = (X - X_min) / (X_max - X_min)
    smallest_nonzero = np.min(X_scaled[X_scaled > 0])
    X_scaled[X_scaled == 0] = smallest_nonzero
    return X_scaled


def haversine(coord1, coord2):
    """
    Calculate the great circle distance in kilometers between two points
    on the earth (specified as tuples of (latitude, longitude))

    Args:
    coord1 (tuple): (Latitude, Longitude) of point 1.
    coord2 (tuple): (Latitude, Longitude) of point 2.

    Returns:
    float: Distance between the two points in kilometers.
    """

    lat1, lon1 = coord1
    lat2, lon2 = coord2

    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of Earth in kilometers
    return c * r


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
    # Check if 'graph' is a NetworkX graph and not empty
    if not isinstance(graph, nx.Graph) or not graph:
        raise ValueError("Invalid graph or graph is empty.")

    # Check if 'pos' is a tuple with two numeric values
    if not isinstance(pos, tuple) or \
            len(pos) != 2 or not \
            all(isinstance(coord, (int, float)) for coord in pos):
        raise ValueError("Position must be a tuple of two numeric values.")

    return astnet.snap_to_node(graph, pos)


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
        if (len(unique_values) == 0 or
                (len(unique_values) == 1 and np.isnan(unique_values[0]))):
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
        ValueError: If 'df', 'variables', or outputs not in the expected
                    format or empty.
    """
    if not isinstance(variables, list):
        raise ValueError("Invalid variables list.")
    if df.empty or not isinstance(df, pd.DataFrame):
        raise ValueError("Invalid DataFrame.")

    # check variables are present
    variables = [v for v in variables if v in df.columns]
    if not variables:
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

    # Check if the cleaned DataFrame is empty
    if cleaned_df.empty:
        raise ValueError(
            "All columns were removed; resulting DataFrame is empty."
        )

    # Check if the cleaned variables list is empty
    if not cleaned_variables:
        raise ValueError(
            "All specified variables were removed; \
                resulting variables list is empty."
        )

    return cleaned_df, cleaned_variables


def read_pickle_network(pik):
    """
    Read a saved network graph (pickle) as an undirected graph, with
    autoStreamTree's reader.

    Args:
        pik: The path to the pickle file containing the network graph.

    Returns:
        A networkx.Graph object representing the undirected graph.
    """
    with warnings.catch_warnings(), \
            contextlib.redirect_stdout(io.StringIO()):
        warnings.filterwarnings(
            "ignore",
            message=r"Unpickling a shapely <2\.0 geometry object.*",
            category=UserWarning
        )
        return astnet.read_network(pik, None)


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
    # Define the data types for the columns
    dtype_dict = {'sample': str, 'lat': float, 'long': float}

    # Read the TSV file, specifying the data types
    points = pd.read_csv(pfile, sep="\t", header=0, dtype=dtype_dict)

    required_columns = {'sample', 'lat', 'long'}
    if not required_columns.issubset(points.columns):
        miss = required_columns - set(points.columns)
        error_message = (
            f"Missing required columns in the points data: {', '.join(miss)}"
        )
        raise ValueError(error_message)

    return points


def vcf_distances(vcf, popmap, coords, out, dist="LINFST", concat="none"):
    """
    Computes the genetic distances among the sites from a VCF with
    autoStreamTree (as its --run GENDIST with -p/--pop), and writes them to
    {out}.popGenDistMat.txt, in the order of the sites in coords.

    Args:
        vcf: Path to the VCF file.
        popmap: Path to a file with an individual and its site on each line
            (whitespace-separated; a header line is ignored). Individuals of
            other sites, or not in the VCF, are left out.
        coords: Path to the sites' coordinates table (as -c).
        out: Output prefix.
        dist: Genetic distance (autoStreamTree's -d/--dist).
        concat: How SNPs are combined (autoStreamTree's -C/--concat).

    Returns:
        str: The path of the matrix written.

    Raises:
        ValueError: If a site has no individuals in the VCF.
    """
    # Imported here: autoStreamTree's genetic distances need pysam
    from autostreamtree.functions import population_distances, read_popmap
    sites = list(read_points_table(coords)["sample"])
    inds = {i: s for i, s in read_popmap(popmap).items() if s in sites}
    mat = population_distances(vcf, inds, dist=dist, concat=concat)
    missing = [s for s in sites if s not in mat.index]
    if missing:
        raise ValueError("No individuals in the VCF for sites: "
                         + ", ".join(missing))
    mat = mat.loc[sites, sites]
    fname = f"{out}.popGenDistMat.txt"
    mat.to_csv(fname, sep="\t")
    return fname


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


def write_edges(out, edge, ids):
    """
    Write edge data to a CSV file.

    This function creates a DataFrame from edge IDs and their corresponding
    resistance values, and writes this data to a CSV file. If the edge array
    contains two columns, it treats them as downstream and upstream
    resistances.

    Args:
        out: The file path where the CSV will be written.
        edge: A numpy array of edge resistance values, possibly with two
              columns.
        ids: A list of edge IDs.

    Returns:
        str: The output file path.
    """
    edge = np.array(edge)
    if edge.ndim == 1:
        df = pd.DataFrame({
            "EDGE_ID": ids,
            "Resistance": edge
        })
    elif edge.ndim == 2 and edge.shape[1] == 2:
        df = pd.DataFrame({
            "EDGE_ID": ids,
            "Resistance_Downstream": edge[:, 0],
            "Resistance_Upstream": edge[:, 1]
        })
    else:
        raise ValueError("Unexpected shape of edge data.")

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

    Elements are returned in column-major order (as R's lower.tri), so that
    they line up with the population pairs given by to_from_().

    Args:
        mat: A NumPy matrix.

    Returns:
        ndarray: An array containing the elements of the lower triangular part
                 of the matrix.
    """
    n = mat.shape[0]
    col, row = np.triu_indices(n, 1)
    return mat[row, col]


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
    to = list()
    frm = list()

    for ia, ib in itertools.combinations(range(1, pops+1), 2):
        to.append(ia)
        frm.append(ib)

    t = to[pops-2]
    tt = frm[pops-2]
    to[pops-2] = tt
    frm[pops-2] = t

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
