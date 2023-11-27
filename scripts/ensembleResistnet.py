#!/usr/bin/env python3
import sys
import os
import random
from datetime import datetime
import getopt
import glob
import pandas as pd
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunner
from resistnet.hall_of_fame import HallOfFame

# Set environment variable
os.environ['USE_PYGEOS'] = '0'


def main():
    """
    Main function to run the resistnet application.
    """
    # Step 0: Read/format input args
    params = ParseArgs()

    # Seed random number generator
    if not params.seed:
        params.seed = datetime.now().timestamp()
    random.seed(params.seed)

    # Check if paths provided as list or directory
    params.paths = None
    if params.input:
        params.paths = params.input
    elif params.input_list:
        params.paths = params.input_list

    # Step 1: Collect models
    # Read models
    local_rows = 1 if params.only_best else None
    hofs = read_and_concat_files(
        params.paths, ".HallOfFame.tsv", local_rows=local_rows
    )

    if (hofs['fitness'] == hofs['aic']).all():
        hofs["fitness"] = hofs["fitness"] * -1
    hofs["aic"] = hofs["aic"] * -1
    hofs.sort_values(
        by='fitness', ascending=False, inplace=True, ignore_index=True
    )

    if params.only_keep:
        hofs = hofs[hofs['keep'] is True]
    if params.hof_max is not None:
        hofs = hofs.head(params.hof_max)

    hofs["keep"] = False
    bests = HallOfFame.from_dataframe(hofs)

    # Get list of variables
    variables = bests.get_variables()
    agg_opts = get_agg_opts(variables, params.varFile, params.edge_agg)

    # Step 2: Build network data
    # Check if samples should be split for each input
    index_col = "sample" if params.split_samples else None

    # Read coordinates
    # If coordinates not provided, collate those from the runs
    if not params.coords:
        coords = read_and_concat_files(
            params.paths, ".pointCoords.txt", index_column=index_col
        )
        coords.columns = ["sample", "lat", "long"]
        coords = coords.drop_duplicates(subset=["lat", "long"], keep='first')
        coords.to_csv(params.out + ".coords",
                      header=True,
                      index=False,
                      quoting=None,
                      sep="\t")
        params.coords = params.out + ".coords"

    # Get network
    network = ResistanceNetwork(
        network=params.network,
        shapefile=params.shapefile,
        coords=params.coords,
        variables=variables,
        agg_opts=agg_opts,
        inmat=None,
        reachid_col=params.reachid_col,
        length_col=params.length_col,
        out=params.out,
        verbose=True
    )

    # Step 3: Set up ModelRunner
    runner = ModelRunner(
        resistance_network=network,
        seed=params.seed,
        verbose=True
    )

    # Step 4: Build ensemble model
    runner.build_ensemble(
        bests=bests,
        awsum=params.awsum,
        only_keep=params.only_keep,
        out=params.out,
        threads=params.GA_procs,
        verbose=True
    )


def read_and_concat_files(
        paths, extension, index_column=None, local_rows=None):
    """
    Reads and concatenates files from given paths with a specific extension.

    Args:
        paths (str or list): A string or a list of paths to search for files.
        extension (str): The file extension to search for.
        index_column (str, optional): Column name to append file index to.
                                      Defaults to None.
        local_rows (int, optional): Number of rows to read from each file.
                                    Defaults to None.

    Returns:
        pandas.DataFrame: A DataFrame concatenated from all read files.
    """
    # Ensure paths is a list
    if isinstance(paths, str):
        paths = [paths]

    dfs = []  # List to store DataFrames

    i = 0
    for path in paths:
        # Use glob to find files with specified extension
        if os.path.isdir(path):
            pattern = os.path.join(path, f"*{extension}")
        else:
            pattern = f"{path}*{extension}"
        for file_path in glob.glob(pattern):
            df = pd.read_csv(file_path, sep="\t", header=0)

            # Append file index to specified column if provided
            if index_column is not None and index_column in df.columns:
                df[index_column] = df[index_column].astype(str) + f"_{i}"

            # Take top rows if local_rows is set
            if local_rows is not None:
                df = df.head(local_rows)

            dfs.append(df)
            i += 1

    # Concatenate all DataFrames
    merged_df = pd.concat(dfs, ignore_index=True)
    merged_df = merged_df.drop_duplicates()

    return merged_df


def get_agg_opts(variables, varFile=None, edge_agg="ARITH"):
    """
    Generates a dictionary of aggregation options for given variables.

    Args:
        variables (list): List of variable names.
        varFile (str, optional): Path to a file containing variable and
                                 aggregation method pairs. Defaults to None.
        edge_agg (str, optional): Default aggregation method. Defaults to
                                  "ARITH".

    Returns:
        dict: A dictionary mapping each variable to its aggregation method.
    """
    agg_opts = {}

    if varFile is not None:
        with open(varFile) as fp:
            for line in fp:
                line = line.strip()
                parts = line.split("\t")
                agg_opts[parts[0]] = parts[1] if len(parts) > 1 else edge_agg
    else:
        for v in variables:
            agg_opts[v] = edge_agg

    return agg_opts


class ParseArgs:
    """
    Class to parse command-line arguments for ensembleResistnet application.

    Attributes are set based on the provided command-line arguments.
    """

    def __init__(self):
        """
        Initializes the ParseArgs object by parsing command-line arguments.
        """
        try:
            options, _ = getopt.getopt(
                sys.argv[1:], 'ho:i:n:s:r:l:c:m:a:L:t:V:XC:',
                ["help", "out=", "in=", "network=", "reps=", "shp=",
                 "len_col=", "id_col=", "split_samples", "max_keep=",
                 "awsum=", "list=", "threads=", "edge_agg=", "varFile=",
                 "allShapes", "report_all", "noPlot", "only_best",
                 "only_keep", "coords=", "seed="]
            )

        except getopt.GetoptError as err:
            print(err)
            self.display_help(
                "\nExiting because getopt returned non-zero exit status."
            )

        # set defaults
        self.set_default_values()

        # parse arguments
        self.set_arguments(options)

        # check mandatory settings
        self.check_mandatory_options()

    def set_default_values(self):
        """ Sets the default values for all the parameters. """
        self.input = None
        self.input_list = None
        self.out = "sim"
        self.variables = None
        self.edge_agg = "ARITH"
        self.agg_opts = dict()
        self.varFile = None
        self.network = None
        self.shapefile = None
        self.hof_max = None
        self.only_keep = False
        self.only_best = False
        self.awsum = 0.95
        self.split_samples = False
        self.length_col = "LENGTH_KM"
        self.dist_col = None
        self.reachid_col = "HYRIV_ID"
        self.seed = None
        self.paths = None
        self.seed = None
        self.GA_procs = 1
        self.minimize = False
        self.allShapes = False
        self.report_all = False
        self.plot = True
        self.coords = None

    def check_mandatory_options(self):
        """ Checks if mandatory options are set """
        if not self.input and not self.input_list:
            self.display_help("No input table provided.")
        if not self.network and not self.shapefile:
            self.display_help("No network provided.")

    def set_arguments(self, options):
        """
        Sets the arguments for the class based on provided options from command
        line.

        Args:
            options (list of tuples): Command line options and their values.
        """

        # First pass to see if help menu was called
        for opt, _ in options:
            if opt in ("-h", "--help"):
                self.display_help("Exiting because help menu was called.")

        for opt, arg_raw in options:
            arg = arg_raw.strip()
            opt = opt.replace("-", "")

            if opt in ("h", "help"):
                continue
            elif opt in ("in", "i"):
                self.input = arg
            elif opt in ("X", "noPlot"):
                self.plot = False
            elif opt in ("list", "L"):
                self.input_list = arg.split(",")
            elif opt in ("out", "o"):
                self.out = arg
            elif opt in ("t", "procs"):
                self.GA_procs = int(arg)
            elif opt in ("seed"):
                self.seed = int(arg)
            elif opt in ("network", "n"):
                self.network = arg
            elif opt in ("s", "shp"):
                self.shapefile = arg
            elif opt in ("awsum", "a"):
                self.awsum = float(arg)
            elif opt in ("max_keep", "m"):
                self.hof_max = int(arg)
            elif opt == "split_samples":
                self.split_samples = True
            elif opt == "only_keep":
                self.only_keep = True
            elif opt == "only_best":
                self.only_best = True
            elif opt == "report_all":
                self.report_all = True
            elif opt in ("C", "coords"):
                self.coords = arg
            elif opt in ("id_col", "c"):
                self.reachid_col = arg
            elif opt in ("l", "len_col"):
                self.length_col = arg
            elif opt == "edge_agg":
                self.edge_agg = arg.upper()
                if self.edge_agg not in (
                    "HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX",
                    "MIN", "SUM", "FIRST", "SD", "VAR", "CV"
                ):
                    self.display_help(
                        "Invalid option " + str(arg).upper() +
                        " for option <--edge_agg>"
                    )
            elif opt in ("V", "varFile"):
                self.varFile = arg
            elif opt == "allShapes":
                self.allShapes = True
            else:
                assert False, f"Unhandled option {opt!r}"

    def display_help(self, message=None):
        """
        Displays the help message for the command-line interface.

        Args:
            message (str, optional): An additional message to print before the
                                     help message.
        """
        if message is not None:
            print()
            print(message)

        print("\nensembleResistnet.py\n")
        print("Utility for model-averaging across ResistNet outputs")

        print(
            "\nArguments:\n"
            "-s, --shp: Path to shapefile\n"
            "-n, --network: Input network\n"
            "-i, --in: Directory containing resistnet outputs\n"
            "-L, --list: Optional comma-separated prefixes (no spaces)\n"
            "-C, --coords: Optional static coordinates file\n\n"

            "Optional Arguments:\n"
            "-t, --procs: Number of parallel processors\n"
            "--seed: RNG seed\n"
            "-a, --awsum: Cumulative Akaike weight threshold [default=0.95]\n"
            "--only_keep: Only retain models where column 'keep'=True\n"
            "--only_best: Only retain best model from each input\n"
            "--split_samples: Treat all samples as unique\n"
            "-X, --noPlot: Turn off plotting\n"
            "-m, --max_keep: Maximum models to keep (default = all models)\n"
            "-l, --len_col: Edge length attribute (def=LENGTH_KM)\n"
            "-c, --id_col: Reach ID attribute (def=EDGE_ID)\n"
            "-o, --out: Output file prefix (default=ensemble)\n"
            "--report_all: Plot full outputs for all retained models\n"
            "--allShapes: Allow inverse and reverse transformations\n"
            "-V, --varfile: Optional file with variables provided like so:\n"
            "          var1 \t <Optional aggregator function>\n"
            "          var2 \t <Optional aggregator function>\n"
            "          ...\n"
        )

        print()
        sys.exit()


# Main script execution
if __name__ == '__main__':
    main()
