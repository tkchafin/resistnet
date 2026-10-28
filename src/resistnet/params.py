import sys
import argparse


class parseArgs:
    def __init__(self):
        parser = argparse.ArgumentParser(
            description=(
                'Optimise resistance networks using Tree-structured Parzen '
                'Estimators (TPE)'
            ),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        # Input options
        parser.add_argument(
            '-g', '--inmat', type=str,
            help='Input distance matrix (usually genetic distances)'
        )
        parser.add_argument(
            '-s', '--shapefile', type=str,
            help='Path to shapefile, geodatabase, or GPKG file'
        )
        parser.add_argument(
            '-c', '--coords', type=str,
            help='Input tsv containing sample coordinates'
        )
        parser.add_argument(
            '-n', '--network', type=str,
            help='Pre-formatted pickled network file'
        )

        # General options
        parser.add_argument(
            '--seed', type=int,
            help='Random number seed'
        )
        parser.add_argument(
            '--reachid_col', type=str, default='HYRIV_ID',
            help='Reach ID in geodatabase'
        )
        parser.add_argument(
            '--length_col', type=str, default='LENGTH_KM',
            help='Length attribute in geodatabase'
        )
        parser.add_argument(
            '-t', '--procs', type=int, default=1,
            help='Number of parallel processors'
        )
        parser.add_argument(
            '-o', '--out', type=str, default='output',
            help='Output file prefix'
        )
        parser.add_argument(
            '-O', '--output_driver', type=str, choices=['SHP', 'GPKG', 'GDB'],
            default='GPKG',
            help='Output driver for annotated geodataframe'
        )

        # Univariate parameter settings
        parser.add_argument(
            '--gridSearch', action='store_true',
            help='Run grid-search to find fixed transformations'
        )
        parser.add_argument(
            '--transFile', type=str,
            help='Read fixed transformations from file'
        )

        # Model optimization/selection options
        parser.add_argument(
            '-v', '--vars', type=str, required=False,
            help='Comma-separated list of variables to use'
        )
        parser.add_argument(
            '-V', '--varfile', type=str, required=False,
            help=(
                'Optional file with variables provided as: '
                'var \\t <Optional aggregator function>'
            )
        )
        parser.add_argument(
            '-F', '--nfail', type=int, default=50, required=False,
            help='Number of iterations without improvement to stop'
        )
        parser.add_argument(
            '-i', '--max_iter', type=int, default=500, required=False,
            help='Maximum number of iterations (one model per chain each)'
        )
        parser.add_argument(
            '-r', '--reps', type=int, default=10, required=False,
            help='Number of independent TPE chains'
        )
        parser.add_argument(
            '-P', '--pweight', type=float, default=0.7, required=False,
            help='TPE prior weight'
        )
        parser.add_argument(
            '-S', '--nstart', type=int, default=20, required=False,
            help='Initial random evaluations per chain'
        )
        parser.add_argument(
            '-C', '--ncand', type=int, default=48, required=False,
            help='TPE expected improvement candidate points'
        )
        parser.add_argument(
            '-G', '--gamma', type=float, default=0.15, required=False,
            help='TPE gamma (0-1; fraction of trials treated as "good")'
        )
        parser.add_argument(
            '-f', '--fitmetric', type=str, choices=['aic', 'loglik', 'r2m',
                                                    'delta'],
            default='aic', required=False,
            help=(
                'Fitness metric used to evaluate models (delta = change in '
                'AIC relative to the null model)'
            )
        )
        parser.add_argument(
            '--fixWeight', action='store_true',
            help='Fix all variable weights at 1.0'
        )
        parser.add_argument(
            '--fixShape', action='store_true',
            help='Turn off variable transformations'
        )
        parser.add_argument(
            '--max_shape', type=int, default=100, required=False,
            help='Maximum shape parameter'
        )
        parser.add_argument(
            '--min_weight', type=float, default=0.0, required=False,
            help='Minimum weight parameter'
        )

        # SAMC model options
        parser.add_argument(
            '--samc', action='store_true',
            help=(
                'Use the spatial absorbing Markov chain (SAMC) model: '
                'conditional first passage times among sites, allowing '
                'directional effects (requires --sizes)'
            )
        )
        parser.add_argument(
            '--sizes', type=str,
            help=(
                'SAMC: tsv with columns sample, size (population sizes; '
                'absorption at each site is 1/(2 * size))'
            )
        )
        parser.add_argument(
            '--infer_origin', type=str, default='NEXT_DOWN',
            help='SAMC: attribute giving the next downstream reach ID'
        )
        parser.add_argument(
            '--fixAsym', action='store_true',
            help='SAMC: do not consider directional effects of variables'
        )

        # Multi-model inference options
        parser.add_argument(
            '-a', '--awsum', type=float, default=1.0, required=False,
            help='Cumulative Akaike weight threshold'
        )
        parser.add_argument(
            '--use_full', action='store_true',
            help='Use full rather than natural model averaging'
        )
        parser.add_argument(
            '--report_all', action='store_true',
            help='Plot outputs for all retained models'
        )

        # Parse arguments and set as attributes
        args = parser.parse_args()
        for key, value in vars(args).items():
            setattr(self, key, value)

        if self.samc and not self.sizes:
            parser.error("--samc requires --sizes")

        # 'delta' refers to the MLPE output column delta_aic_null
        if self.fitmetric == 'delta':
            self.fitmetric = 'delta_aic_null'

        # Set variables based on vars or varfile
        if self.vars:
            self.variables = self.vars.split(',')
        elif self.varfile:
            self.variables = []
            self.agg_opts = {}
            with open(self.varfile, 'r') as file:
                for line in file:
                    fields = line.strip().split('\t')
                    if not fields[0]:
                        continue
                    self.variables.append(fields[0])
                    # Optional aggregator function (default ARITH)
                    self.agg_opts[fields[0]] = (
                        fields[1].strip() if len(fields) > 1 else "ARITH"
                    )
        else:
            print("Error: Either --vars or --varfile must be specified.")
            sys.exit(1)
        if self.vars:
            self.agg_opts = {var: "ARITH" for var in self.variables}
