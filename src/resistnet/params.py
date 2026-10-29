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

        # Genetic distances from a VCF (with autoStreamTree), instead of -g
        parser.add_argument(
            '--vcf', type=str,
            help=(
                'VCF of the individuals, to compute the genetic distances '
                'among sites with autoStreamTree (instead of -g; needs '
                '--popmap)'
            )
        )
        parser.add_argument(
            '--popmap', type=str,
            help=(
                'With --vcf: tsv of individual and site (the site labels '
                'of -c) on each line'
            )
        )
        parser.add_argument(
            '--gendist', type=str, default='LINFST',
            choices=['LINFST', 'FST', 'JOST', 'HARMD', 'CHORD', 'PDIST'],
            help=(
                "With --vcf: genetic distance (autoStreamTree's -d/--dist): "
                "LINFST = FST / (1 - FST), FST (Weir & Cockerham's theta), "
                "JOST (Jost's D), HARMD (harmonic mean of Jost's D over "
                "loci), CHORD (Cavalli-Sforza & Edwards), PDIST (mean "
                "p-distance between individuals)"
            )
        )
        parser.add_argument(
            '--concat', type=str, default='none',
            choices=['none', 'loc', 'all'],
            help=(
                "With --vcf: combine SNPs (autoStreamTree's -C/--concat): "
                "none (each SNP a locus), loc (the SNPs of each CHROM joined "
                "into one locus, e.g. RAD loci), or all"
            )
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
            help=(
                'Output directory and prefix: -o DIR/NAME writes each '
                "model's outputs to DIR/NAME/{model}/NAME.* (and, with "
                '--model both, the comparison to DIR/NAME/NAME.*)'
            )
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
            '-F', '--nfail', type=int, default=200, required=False,
            help=(
                'Number of iterations without improvement to stop (counted '
                'after the --nstart random start-up iterations)'
            )
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
            '--max_shape', type=int, default=15, required=False,
            help=('Maximum shape parameter (shapes from 1 are searched '
                  'log-uniformly)')
        )
        parser.add_argument(
            '--min_weight', type=float, default=0.0, required=False,
            help='Minimum weight parameter'
        )

        # Model and SAMC model options
        parser.add_argument(
            '--model', type=str, choices=['resistance', 'samc', 'both'],
            default=None,
            help=(
                'Model to fit: resistance (sums of reach resistances along '
                'the paths between sites; default), samc (spatial absorbing '
                'Markov chain: conditional first passage times among sites, '
                'allowing directional effects), or both (fitted separately, '
                'then ranked together)'
            )
        )
        parser.add_argument(
            '--samc', action='store_true',
            help='Same as --model samc'
        )
        parser.add_argument(
            '--per_reach', action='store_true',
            help=(
                'Standard model: each reach costs its composite value '
                'whatever its length (the formulation of earlier versions), '
                'instead of a resistance per km times its length'
            )
        )
        parser.add_argument(
            '--sizes', type=str,
            help=(
                'SAMC: tsv with columns sample, size (effective population '
                'sizes; absorption at each site is 1/(2 * size))'
            )
        )
        parser.add_argument(
            '--ne', type=float, default=100,
            help='SAMC: effective population size of every site, if no --sizes'
        )
        parser.add_argument(
            '--infer_origin', type=str, default='NEXT_DOWN',
            help='SAMC: attribute giving the next downstream reach ID'
        )
        parser.add_argument(
            '--max_resist', type=float, default=None,
            help=(
                'Upper limit for the maximum resistance (relative to a '
                'minimum of 1), which is optimised: per km for the standard '
                "model (default 1000, as ResistanceGA's max.cont), and of "
                'the composite surface for SAMC (default 500)'
            )
        )
        parser.add_argument(
            '--fix_resist', type=float,
            help='Fix the maximum resistance instead of optimising it'
        )
        parser.add_argument(
            '--fixAsym', action='store_true',
            help='SAMC: do not consider directional effects of variables'
        )
        parser.add_argument(
            '--max_flow', type=float, default=1.0,
            help=(
                'SAMC: upper limit for the rate of the flow (upstream vs '
                'downstream) term, per km'
            )
        )
        parser.add_argument(
            '--fixFlow', action='store_true',
            help='SAMC: leave out the flow (upstream vs downstream) term'
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

        if self.model is None:
            self.model = "samc" if self.samc else "resistance"
        elif self.samc and self.model != "samc":
            parser.error(f"--samc conflicts with --model {self.model}")
        self.models = (["resistance", "samc"] if self.model == "both"
                       else [self.model])
        if self.vcf and self.inmat:
            parser.error("--vcf and -g/--inmat can't both be used")
        if self.vcf and not self.popmap:
            parser.error("--vcf needs --popmap")
        if self.max_resist is not None and self.max_resist <= 1:
            parser.error("--max_resist must be greater than 1")
        if self.fix_resist is not None and self.fix_resist < 1:
            parser.error("--fix_resist must be at least 1")
        # The smallest flow rate searched (samc_network.MIN_FLOW)
        if self.max_flow <= 0.001:
            parser.error("--max_flow must be greater than 0.001")

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

    def max_resist_for(self, model):
        """
        The upper limit of the maximum resistance for a model: --max_resist,
        or by default 1000 for the standard model (per km, as ResistanceGA's
        max.cont) and 500 for SAMC (of the composite surface).

        Args:
            model (str): "resistance" or "samc".

        Returns:
            float: The upper limit.
        """
        if self.max_resist is not None:
            return float(self.max_resist)
        return 500.0 if model == "samc" else 1000.0
