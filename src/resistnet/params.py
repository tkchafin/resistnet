import argparse


class parseArgs:
    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Genetic algorithm to optimize resistance networks',
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
            '-t', '--procs', type=int,
            help='Number of parallel processors'
        )
        parser.add_argument(
            '-o', '--out', type=str, default='output',
            help='Output file prefix'
        )
        parser.add_argument(
            '-O', '--gdf_out', type=str, choices=['SHP', 'GPKG', 'GDB'],
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
            '-v', '--vars', type=str,
            help='Comma-separated list of variables to use'
        )
        parser.add_argument(
            '-V', '--varfile', type=str,
            help=(
                'Optional file with variables provided as: '
                'var \\t <Optional aggregator function>'
            )
        )
        parser.add_argument(
            '-F', '--nfail', type=int, default=50,
            help='Number of failed gens to stop optimization'
        )
        parser.add_argument(
            '-i', '--max_iter', type=int, default=500,
            help='Maximum number of generations'
        )
        parser.add_argument(
            '-r', '--reps', type=int, default=10,
            help='Number of independent TPE chains'
        )
        parser.add_argument(
            '-P', '--pweight', type=float, default=0.7,
            help='Prior weight'
        )
        parser.add_argument(
            '-S', '--nstart', type=int, default=20,
            help='Initial random evaluations'
        )
        parser.add_argument(
            '-C', '--ncand', type=int, default=48,
            help='EI candidate points'
        )
        parser.add_argument(
            '-G', '--gamma', type=float, default=0.15,
            help='Exploration factor'
        )
        parser.add_argument(
            '-f', '--fitmetric', type=str, choices=['aic', 'loglik', 'r2m',
                                                    'delta'],
            default='loglik',
            help='Fitness metric used to evaluate models'
        )
        parser.add_argument(
            '--max_hof_size', type=int, default=100,
            help='Maximum models retained'
        )
        parser.add_argument(
            '--fixWeight', action='store_true',
            help='Fix the weight parameter during optimization'
        )
        parser.add_argument(
            '--fixShape', action='store_true',
            help='Fix the shape parameter during optimization'
        )
        parser.add_argument(
            '--max_shape', type=int, default=100,
            help='Maximum shape parameter'
        )
        parser.add_argument(
            '--min_weight', type=float, default=0.0,
            help='Minimum weight parameter'
        )

        # Multi-model inference options
        parser.add_argument(
            '-a', '--awsum', type=float, default=1.0,
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

        # Set variables based on vars or varfile
        if self.vars:
            self.variables = self.vars.split(',')
        elif self.varfile:
            self.variables = []
            with open(self.varfile, 'r') as file:
                for line in file:
                    var = line.strip().split('\t')[0]
                    self.variables.append(var)
        else:
            self.variables = None
        self.agg_opts = {var: "ARITH" for var in self.variables}
