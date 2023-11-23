import sys
import getopt


class parseArgs():
    """
    A class to parse command-line arguments for ResistNet

    Attributes:
        Various attributes for command-line options and their default values.

    Methods:
        display_help(message=None): Prints the help message and exits
    """
    def __init__(self):
        # Define options
        try:
            options, r = getopt.getopt(
                sys.argv[1:],
                'hp:g:s:s:T:P:G:s:m:i:c:t:F:d:D:f:b:C:v:V:Aa:o:Xj:n:',
                ["shp=", "help", "input=", "prefix=", "genmat=", "shapefile=",
                 "seed=", "procs=", "maxPop=", "maxpop=", "maxgen=", "maxGen=",
                 "size=", "popsize=", "mutpb=", "indpb=", "cxpb=", "tourn=",
                 "nfail=", "nFail=", "delt=", "deltP=", "deltp=", "fit=",
                 "metric=", "fitness=", "burn=", "dist_col=", "vars=",
                 "pop_agg=", "varfile=", "maxShape=", "awsum=", "report_all",
                 "noPlot", "out=", "keep_all", "minWeight=", "max_hof_size=",
                 "posWeight", "fixWeight", "allShapes", "efit_agg=", "coords=",
                 "length_col=", "reachid_col=", "minimize", "network=",
                 "fixShape", "max_gens=", "max_gen=", "max_pop="
                 ]
            )
        except getopt.GetoptError as err:
            print(err)
            self.display_help(
                "\nExiting because getopt returned non-zero exit status."
            )

        # Default values for params
        self.prefix = "out"
        self.dist_col = None
        self.variables = None
        self.agg_opts = dict()
        self.varFile = None
        self.minimize = False
        self.seed = None
        self.installCS = False
        self.popsize = None
        self.length_col = "LENGTH_KM"
        self.reachid_col = "HYRIV_ID"
        self.pop_agg = "ARITH"
        self.edge_agg = "ARITH"
        self.efit_agg = "SUM"
        self.maxpopsize = 1000
        self.cstype = "pairwise"
        self.fitmetric = "aic"
        self.network = None
        self.predicted = False
        self.inmat = None
        self.shapefile = None
        self.network = None
        self.coords = None
        self.cholmod = False
        self.GA_procs = 1
        self.CS_procs = 1
        self.deltaB = None
        self.deltaB_perc = 0.001
        self.nfail = 50
        self.maxGens = 500
        self.tournsize = 10
        self.cxpb = 0.5
        self.mutpb = 0.5
        self.indpb = 0.5
        self.burnin = 0
        self.max_hof_size = 100
        self.out = "output"
        self.awsum = 0.95
        self.modavg = True
        self.report_all = False
        self.plot = True

        self.posWeight = False
        self.fixWeight = False
        self.allShapes = False
        self.fixShape = False
        self.max_shape = 100

        self.min_weight = 0.0

        self.only_keep = True
        self.julia = "julia"
        self.compiled_modules = True
        self.sys_image = None

        # First pass to see if help menu was called
        for o, a in options:
            if o in ("-h", "-help", "--help"):
                self.display_help("Exiting because help menu was called.")

        # Second pass to set all args.
        for opt, arg in options:
            arg = arg.strip()
            opt = opt.replace("-", "")
            if opt in ('p', 'prefix'):
                self.prefix = arg
            elif opt in ('g', 'genmat'):
                self.inmat = arg
            elif opt in ('s', 'shp'):
                self.shapefile = arg
            elif opt in ("c", "coords"):
                self.coords = arg
            elif opt in ("n", "network"):
                self.network = arg
            elif opt == "minimize":
                self.minimize = True
            elif opt == 'seed':
                self.seed = int(arg)
            elif opt == 't' or opt == 'procs':
                self.GA_procs = int(arg)
            elif opt in ('P', 'maxPop', 'maxpop', "max_pop"):
                self.maxpopsize = int(arg)
            elif opt in ("G", "maxGen", "maxgen", "max_gen", "max_gens"):
                self.maxGens = int(arg)
            elif opt in ("popsize", "size"):
                self.popsize = int(arg)
            elif opt == "minWeight":
                self.min_weight = float(arg)
            elif opt == "pop_agg":
                self.pop_agg = arg.upper()
                if self.pop_agg not in [
                    "HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX",
                    "MIN", "SUM", "FIRST", "SD", "VAR", "CV"
                ]:
                    self.display_help(
                        "Invalid option " + str(arg).upper() +
                        " for option <--pop_agg>"
                    )
            elif opt == "edge_agg":
                self.edge_agg = arg.upper()
                if self.edge_agg not in [
                    "HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX",
                    "MIN", "SUM", "FIRST", "SD", "VAR", "CV"
                ]:
                    self.display_help(
                        "Invalid option " + str(arg).upper() +
                        " for option <--edge_agg>"
                    )
            elif opt == "mutpb":
                self.mutpb = float(arg)
            elif opt == "indpb":
                self.indpb = float(arg)
            elif opt == "cxpb":
                self.cxpb = float(arg)
            elif opt in ('T', 'tSize', 'tournSize'):
                self.tournsize = int(arg)
            elif opt in ('F', 'nfail', 'nFail'):
                self.nfail = int(arg)
            elif opt in ('d', 'delt'):
                self.deltaB = float(arg)
            elif opt in ('D', 'deltP', 'deltp'):
                self.deltaB_perc = float(arg)
            elif opt in ('f', 'fit', 'metric', 'fitness'):
                if arg.lower() not in ["aic", "r2m", "loglik", "delta"]:
                    self.diplay_help("Unrecognized fitness metric <-f, --fit>")
                else:
                    self.fitmetric = arg.lower()
            elif opt in ('b', 'burn'):
                self.burnin = int(arg)
            elif opt == "dist_col":
                self.dist_col = arg
            elif opt == "infer":
                self.predicted = True
            elif opt == "reachid_col":
                self.reachid_col = arg
            elif opt == "length_col":
                self.length_col = arg
            elif opt == "cholmod":
                self.cholmod = True
            elif opt in ('v', 'vars'):
                self.variables = arg.split(",")
            elif opt in ('V', 'varFile'):
                self.varFile = arg
            elif opt in ('A', 'modavg', 'modAvg'):
                self.modavg = True
            elif opt in ('a', 'awsum'):
                self.awsum = float(arg)
            elif opt == "maxShape":
                self.max_shape = float(arg)
            elif opt == "report_all":
                self.report_all = True
            elif opt in ('X', 'noPlot'):
                self.plot = False
            elif opt in ('o', 'out'):
                self.out = arg
            elif opt == "avgall":
                self.only_keep = False
            elif opt == "efit_agg":
                self.efit_agg = arg.upper()
                if self.efit_agg not in [
                    "HARM", "ADJHARM", "ARITH", "GEOM", "MEDIAN", "MAX", "MIN",
                    "SUM", "FIRST", "SD", "VAR", "CV"
                ]:
                    self.display_help(
                        "Invalid option " + str(arg).upper() +
                        " for option <--efit_agg>"
                    )
            elif opt == "posWeight":
                self.posWeight = True
            elif opt == "fixWeight":
                self.fixWeight = True
            elif opt == "allShapes":
                self.allShapes = True
            elif opt == "fixShape":
                self.fixShape = True
            elif opt in ('h', 'help'):
                pass
            else:
                assert False, f"Unhandled option {opt!r}"

        if self.posWeight and self.fixWeight:
            self.display_help(
                "--posWeight and --fixWeight cannot be used together"
            )

        if self.varFile is not None:
            if self.variables is not None:
                print(
                    "Warning: Variables were specified with both <-v> and "
                    "<-V>... Over-riding options using file provided "
                    "with <-V>"
                )
            if self.edge_agg is None:
                self.edge_agg = "ARITH"
            with open(self.varFile) as fp:
                for line in fp:
                    line = line.strip()
                    stuff = line.split("\t")
                    if len(stuff) < 2:
                        self.agg_opts[stuff[0]] = self.edge_agg
                    else:
                        self.agg_opts[stuff[0]] = stuff[1]
            self.variables = list(self.agg_opts.keys())
        else:
            for v in self.variables:
                self.agg_opts[v] = self.edge_agg

        if not self.variables:
            self.display_help("No variables selected.")

    def display_help(self, message=None):
        """
        Displays the help message and exits the program.

        Args:
            message (str, optional): An additional message to print before the
                                     help message.
        """
        if message is not None:
            print()
            print(message)
        print("\nresistnet.py\n")
        print("Author: Tyler K Chafin")
        print("Contact: tyler.chafin@bioss.ac.uk")
        print("Description: Genetic algorithm to optimise resistance networks")
        print(
            "\nInput options:\n"
            "    -g, --genmat: Genetic distance matrix\n"
            "    -s, --shp: Path to shapefile, geodatabase, or GPKG file\n"
            "    -c, --coords: Input tsv containing sample coordinates\n\n"

            "General options:\n"
            "    --seed: Random number seed (default=taken from clock time)\n"
            "    --reachid_col: Reach ID [default=\"REACH_ID\"]\n"
            "    --length_col: Length [default=\"LENGTH_KM\"]\n"
            "    -t, --procs: Number of parallel processors\n"
            "    -X, --noPlot: Turn off plotting\n"
            "    -o, --out: Output file prefix\n"
            "    -h, --help: Displays help menu\n\n"

            "Aggregation options:\n"
            "    --edge_agg: Method for combining variables across segments\n"
            "    --pop_agg: Method to combine population genetic distances\n"
            "      Options: ARITH (-metic mean), MEDIAN, HARM (-monic mean),\n"
            "      ADJHARM (adjusted HARM, see docs), GEOM (geometric mean),\n"
            "      MIN, MAX, FIRST, SD (standard deviation), VAR (variance),\n"
            "      SUM (simple sum), CV (coefficient of variation = SD/MEAN)\n"

            "\nGenetic Algorithm Options:\n"
            "    -P, --maxPop: Maximum population size [default = 100]\n"
            "    -G, --maxGen: Maximum number of generations [default = 500]\n"
            "    -s, --size: Manually set population size to <-p int>,\n"
            "        NOTE: By default, #params * 15\n"
            "    -m, --mutpb: Mutation probability per trait [default=0.5]\n"
            "    --indpb: Mutation probability per individual [default=0.5]\n"
            "    --cxpb: Cross-over probability [default=0.5]\n"
            "    -T, --tSize: Tournament size [default=10]\n"
            "    --posWeight: Constrain parameter weights to between 0.0-1.0\n"
            "    --minWeight: Sets minimum allowable weight (w/--posWeight)\n"
            "    --fixWeight: Constrain parameter weights to 1.0\n"
            "    --fixShape: Turn off feature transformation\n"
            "    --allShapes: Allow inverse and reverse transformations\n"
            "    --maxShape: Maximum shape value [default=100]\n\n"

            "Model optimization/selection options:\n"
            "    -v, --vars: Comma-separated list of variables to use\n"
            "    -V, --varfile: Optional file with variables provided as:\n"
            "        var1 \\t <Optional aggregator function>\n"
            "        var2 \\t <Optional aggregator function>\n"
            "        ...\n"
            "    -F, --nfail: Number of failed gens to stop optimization\n"
            "    -d, --delt: Threshold absolute change in fitness [def.=0.0]\n"
            "    -D, --deltP: Threshold as decimal percentage [def.=0.001]\n"
            "    -f, --fit: Fitness metric used to evaluate models\n"
            "        Options: aic (default), loglik (log-likelihood),\n"
            "        r2m (marginal R^2), delta (Change in AIC vs null model)\n"
            "        NOTE: Case-insensitive\n"
            "    -b, --burn: Number of generations for pre-burnin [def.=0]\n"
            "    --max_hof_size: Maximum models retained [default=100]\n\n"

            "Multi-model inference options:\n"
            "    -a, --awsum: Cumulative Akaike weight threshold [def.=0.95]\n"
            "    --report_all: Plot outputs for all retained models\n\n"
        )
        print()
        sys.exit()
