import sys
import getopt
from resistnet.resistance_network import SimResistanceNetwork


def main():
    params = ParseArgs()

    # Step 1: Read network
    sim_engine = SimResistanceNetwork(
        network=params.network,
        reachid_col=params.id_col,
        length_col=params.length_col,
        verbose=True
    )

    # Step 2: Write simulated outputs
    sim_engine.simulate(
        spec_file=params.input,
        num_reps=params.reps,
        num_samples=params.samples,
        out=params.out
    )


class ParseArgs:
    def __init__(self):
        # Define options
        try:
            options, _ = getopt.getopt(
                sys.argv[1:], 'ho:i:n:s:r:l:c:',
                ["help", "out=", "in=", "network=", "reps=", "samples=",
                 "len_col=", "id_col="]
            )
        except getopt.GetoptError as err:
            print(err)
            self.display_help(
                "\nExiting because getopt returned non-zero exit status."
            )

        # Default values for params
        self.input = None
        self.out = "sim"
        self.network = None
        self.reps = 1
        self.samples = 50
        self.length_col = "LENGTH_KM"
        self.id_col = "EDGE_ID"

        # First pass to see if help menu was called
        for o, _ in options:
            if o in ("-h", "--help"):
                self.display_help("Exiting because help menu was called.")

        # Second pass to set all args.
        for opt, arg_raw in options:
            arg = arg_raw.strip()
            opt = opt.lstrip('-')
            if opt in ("h", "help"):
                continue
            elif opt in ("in", "i"):
                self.input = arg
            elif opt in ("out", "o"):
                self.out = arg
            elif opt in ("network", "n"):
                self.network = arg
            elif opt in ("samples", "s"):
                self.samples = int(arg)
            elif opt in ("reps", "r"):
                self.reps = int(arg)
            elif opt in ("id_col", "c"):
                self.id_col = arg
            elif opt in ("len_col", "l"):
                self.length_col = arg
            else:
                assert False, f"Unhandled option {opt!r}"

        # Check mandatory options are set
        if not self.input:
            self.display_help("No input table provided.")
        if not self.network:
            self.display_help("No network provided.")

    def display_help(self, message=None):
        if message is not None:
            print()
            print(message)
        print("\nsimResistnet.py\n")
        print("Author: Tyler Chafin")
        print("Description: Simulate data on a given network")
        print("""
Arguments:
    -n, --network: Input network (pickle'd networkx output)
    -i, --in: Table giving variables to use to generate resistnet input
    -r, --reps: Number of replicates
    -s, --samples: Number of random nodes to sample
    -l, --len_col: Edge length attribute (def=LENGTH_KM)
    -c, --id_col: Reach ID attribute (def=EDGE_ID)
    -o, --out: Output file name (default=sim)
""")
        sys.exit()


# Call main function
if __name__ == '__main__':
    main()