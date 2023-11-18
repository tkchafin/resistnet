import sys
import os

os.environ['USE_PYGEOS'] = '0'

from datetime import datetime
import random
import getopt
import glob
import pandas as pd

from resistnet.params import parseArgs
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunner
from resistnet.hall_of_fame import HallOfFame

def main():

    #########################################################
    # Step 0: Read/ format input args
    #########################################################

    params = parseArgs()

    #seed random number generator
    #random.seed(params.seed)
    if not params.seed:
        params.seed=datetime.now().timestamp()
    random.seed(params.seed)

    # check if paths provided as list or directory 
    params.paths = None
    if params.input:
        params.paths = params.input
    elif params.input_list:
        params.paths = params.input_list
    
    print(params.paths)

    #########################################################
    # Step 1: Collect models 
    #########################################################

    # read models 
    if params.only_best:
        local_rows=1
    else:
        local_rows=None
    hofs = read_and_concat_files(params.paths, ".HallOfFame.tsv", local_rows=local_rows)
    if (hofs['fitness'] == hofs['aic']).all():
        hofs["fitness"] = hofs["fitness"]*-1
    hofs["aic"] = hofs["aic"]*-1
    hofs.sort_values(by='fitness', ascending=False, inplace=True, ignore_index=True)
    if params.only_keep:
        hofs = hofs[hofs['keep'] == True]
    if params.hof_max is not None:
        hofs = hofs.head(params.hof_max)
    hofs["keep"] = False
    bests = HallOfFame.from_dataframe(hofs)

    # get list of variables 
    variables = bests.get_variables()
    agg_opts = get_agg_opts(variables, params.varFile, params.edge_agg)

    #########################################################
    # Step 2: Build network data 
    #########################################################

    # check if samples should be split for each input 
    index_col=None
    if params.split_samples:
        index_col="sample"

    # read coordinates 
    # if coordinates not provided, collate those from the runs 
    if not params.coords:
        coords = read_and_concat_files(params.paths, ".pointCoords.txt", index_column=index_col)
        coords.columns = ["sample", "lat", "long"]
        coords = coords.drop_duplicates(subset=["lat", "long"], keep='first')
        coords.to_csv(params.out+".coords", 
            header=True, 
            index=False, 
            quoting=None,
            sep="\t")
        params.coords=params.out+".coords"


    # get network 
    network = ResistanceNetwork(
        network = params.network,
        shapefile = params.shapefile,
        coords = params.coords,
        variables = variables, 
        agg_opts = agg_opts,
        inmat = None,
        reachid_col=params.reachid_col,
        length_col=params.length_col,
        out = params.out, 
        verbose=True
    )


    # #########################################################
    # # Step 3: Set up ModelRunner 
    # #########################################################

    runner = ModelRunner(
        resistance_network=network,
        seed = params.seed,
        verbose=True
    )


    # #########################################################
    # # Step 4: Build ensemble model 
    # #########################################################

    runner.build_ensemble(
        bests = bests,
        awsum = params.awsum, 
        only_keep = params.only_keep,
        out = params.out,
        threads = params.GA_procs, 
        verbose = True)


def read_and_concat_files(paths, extension, index_column=None, local_rows=None):
    # Ensure paths is a list, even if a single string is provided
    if isinstance(paths, str):
        paths = [paths]
    
    # Initialize an empty list to store the DataFrames
    dfs = []

    # Loop through all paths in the input list
    i=0
    for path in paths:
        # Use glob to find all files with the specified extension
        if os.path.isdir(path):
            pattern = os.path.join(path, f"*{extension}")
        else:
            pattern = f"{path}*{extension}"
        for file_path in glob.glob(pattern):
            # Read the file into a DataFrame 
            df = pd.read_csv(file_path, sep="\t", header=0)

            # Add file index to the specified column if index_column is provided
            if index_column is not None and index_column in df.columns:
                df[index_column] = df[index_column].astype(str) + f"_{i}"
            
            # if local_rows set, only take top X rows from each data frame 
            if local_rows is not None:
                df = df.head(local_rows)

            # append it to list 
            dfs.append(df)

            i+=1

    # Concatenate all the DataFrames in the list
    merged_df = pd.concat(dfs, ignore_index=True)
    merged_df = merged_df.drop_duplicates()

    return merged_df


def get_agg_opts(variables, varFile=None, edge_agg="ARITH"):
    agg_opts = dict()
    if varFile is not None:
        with open(varFile) as fp:
            for line in fp:
                line=line.strip()
                stuff=line.split("\t")
                if len(stuff) < 2:
                    agg_opts[stuff[0]]=edge_agg
                else:
                    agg_opts[stuff[0]]=stuff[1]
    else:
        for v in variables:
            agg_opts[v]=edge_agg
    return agg_opts 


#Object to parse command-line arguments
class parseArgs():
    def __init__(self):
        #Define options
        try:
            options, remainder = getopt.getopt(sys.argv[1:], 'ho:i:n:s:r:l:c:m:a:L:t:V:XC:', \
            ["help", "out=", "in=", "network=", "reps=", "shp=",
            "len_col=", "id_col=", "split_samples", "max_keep=", 
            "awsum=", "list=", "threads=", "edge_agg=", "varFile=",
            "allShapes", "report_all", "noPlot", "only_best", 
            "only_keep", "coords="])
        except getopt.GetoptError as err:
            print(err)
            self.display_help("\nExiting because getopt returned non-zero exit status.")
        #Default values for params
        #Input params
        self.input=None
        self.input_list=None
        self.out="sim"
        self.variables = None
        self.edge_agg="ARITH"
        self.agg_opts = dict()
        self.varFile=None
        self.network=None
        self.shapefile=None
        self.hof_max=None
        self.only_keep=False
        self.only_best=False
        self.awsum=0.95
        self.split_samples=False
        self.length_col="LENGTH_KM"
        self.dist_col=None
        self.reachid_col="HYRIV_ID"
        self.recalc_aic=False #not used currently 
        self.paths = None
        self.seed = None
        self.GA_procs = 1
        self.minimize=False
        self.allShapes=False
        self.report_all=False
        self.plot=True
        self.coords=None


        #First pass to see if help menu was called
        for o, a in options:
            if o in ("-h", "-help", "--help"):
                self.display_help("Exiting because help menu was called.")

        #Second pass to set all args.
        for opt, arg_raw in options:
            arg = arg_raw.replace(" ","")
            arg = arg.strip()
            opt = opt.replace("-","")
            #print(opt,arg)
            if opt == "h" or opt == "help":
                continue
            elif opt=="in" or opt=="i":
                self.input=arg
            elif opt=="X" or opt=="noPlot":
                self.plot=False
            elif opt=="list" or opt=="L":
                self.input_list=arg.split(",")
            elif opt=="out" or opt=="o":
                self.out=arg
            elif opt=='t' or opt=='procs':
                self.GA_procs=int(arg)
            elif opt == "network" or opt=="n":
                self.network=arg
            elif opt=='s' or opt=='shp':
                self.shapefile = arg
            elif opt == "awsum" or opt=="a":
                self.awsum=float(arg)
            elif opt == "max_keep" or opt=="m":
                self.hof_max=int(arg)
            elif opt == "split_samples":
                self.split_samples=True
            elif opt == "only_keep":
                self.only_keep=True
            elif opt == "only_best":
                self.only_best=True
            elif opt=="report_all":
                self.report_all=True
            elif opt == "C" or opt == "coords":    
                self.coords = arg
            elif opt == "id_col" or opt=="c":
                self.reachid_col=arg
            elif opt == "l" or opt=="len_col":
                self.length_col=arg
            elif opt=="edge_agg":
                self.edge_agg = arg.upper()
                if self.edge_agg not in ["HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN", "SUM", "FIRST", "SD", "VAR", "CV"]:
                    self.display_help("Invalid option "+str(arg).upper()+" for option <--edge_agg>")
            elif opt=="V" or opt=="varFile":
                self.varFile=arg
            elif opt=="allShapes":
                self.allShapes=True
            else:
                assert False, "Unhandled option %r"%opt

        #Check manditory options are set
        if not self.input and not self.input_list:
            self.display_help("No input table provided.")
        if not self.network and not self.shapefile:
            self.display_help("No network provided.")

    def display_help(self, message=None):
        if message is not None:
            print()
            print (message)
        print ("\nensembleResistnet.py\n")
        print("Author: Tyler Chafin")
        print ("Description: Utility script for model-averaging across a collection of RestistNet outputs")
        print("""
Arguments:
-s,--shp    : Path to shapefile containing cleaned, contiguous stream reaches
-n,--network    : Input network (pickle'd networkx output; should include all sites)
-i,--in        : Directory containing resistnet outputs
-L,--list    : Alternatively, provide a comma-separated list of model output prexifes (without spaces)
-C,--coords  : Static coordinates, if used .coords files will not be read using the -L/--list prefixes

Optional Arguments:
-t,--procs    : Number of parallel processors
-a,--awsum    : Cumulative Akaike weight threshold to retain top N models [default=0.95]
--only_keep    : Only retain models where column "keep"=True
--only_best    : Only retain best model from each input
--split_samples    : Do not check for overlap in sample names, instead treating all as unique
-X,--noPlot    : Turn off plotting
-m,--max_keep    : Maximum models to keep (default = all models)
-l,--len_col    : Attribute in network corresponding to edge length (def=LENGTH_KM)
-c,--id_col    : Attribute in network corresponding to edge ID (def=EDGE_ID)
-o,--out    : Output file prefix (default=ensemble)
--report_all    : Plot per-stream resistance and generate full outputs for all retained models
--allShapes    : Allow inverse and reverse transformations
-V,--varfile    : Optional file with variables provided like so:
          var1 \t <Optional aggregator function>
          var2 \t <Optional aggregator function>
          ...
          ...
""")
        print()
        sys.exit()

#Call main function
if __name__ == '__main__':
    main()
