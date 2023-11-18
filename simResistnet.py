import sys
import os

import getopt

from resistnet.resistance_network import SimResistanceNetwork

def main():
    params = parseArgs()

    #########################################################
    # Step 1: Read network
    #########################################################

    sim_engine = SimResistanceNetwork(
        network=params.network,
        reachid_col=params.id_col,
        length_col=params.length_col,
        verbose = True
    )

    #########################################################
    # Step 2: Write simulated outputs
    #########################################################

    sim_engine.simulate(
        spec_file = params.input,
        num_reps = params.reps,
        num_samples = params.samples,
        out = params.out
    )


#Object to parse command-line arguments
class parseArgs():
    def __init__(self):
        #Define options
        try:
            options, remainder = getopt.getopt(sys.argv[1:], 'ho:i:n:s:r:l:c:', \
            ["help", "out=", "in=", "network=", "reps=", "samples=",
            "len_col=", "id_col="])
        except getopt.GetoptError as err:
            print(err)
            self.display_help("\nExiting because getopt returned non-zero exit status.")
        #Default values for params
        #Input params
        self.input=None
        self.out="sim"
        self.network=None
        self.reps=1
        self.samples=50
        self.length_col="LENGTH_KM"
        self.id_col="EDGE_ID"


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
            elif opt=="out" or opt=="o":
                self.out=arg
            elif opt == "network" or opt=="n":
                self.network=arg
            elif opt == "samples" or opt=="s":
                self.samples=int(arg)
            elif opt == "reps" or opt=="r":
                self.reps=int(arg)
            elif opt == "id_col" or opt=="c":
                self.id_col=arg
            elif opt == "l" or opt=="len_col":
                self.length_col=arg
            else:
                assert False, "Unhandled option %r"%opt

        #Check manditory options are set
        if not self.input:
            self.display_help("No input table provided.")
        if not self.network:
            self.display_help("No network provided.")

    def display_help(self, message=None):
        if message is not None:
            print()
            print (message)
        print ("\nsimResistnet.py\n")
        print("Author: Tyler Chafin")
        print ("Description: Simulate data on a given network for validating resistnet")
        print("""
Arguments:
-n,--network    : Input network (pickle'd networkx output)
-i,--in        : Table giving information on which variables to use to generate resistnet input
-r,--reps    : Number of replicates
-s,--samples    : Number of random nodes to sample
-l,--len_col    : Attribute in network corresponding to edge length (def=LENGTH_KM)
-c,--id_col    : Attribute in network corresponding to edge ID (def=EDGE_ID)
-o,--out    : Output file name (default=sim)
""")
        print()
        sys.exit()

#Call main function
if __name__ == '__main__':
    main()
