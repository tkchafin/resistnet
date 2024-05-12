import sys
import os
import tempfile
from unittest.mock import patch
from resistnet.params import parseArgs


def test_parse_args_basic():
    test_args = [
        "script_name", "-o", "my_prefix", "-g", "genmat_file",
        "--shp", "shapefile.shp", "-v", "example1"
    ]
    sys.argv = test_args
    args = parseArgs()
    assert args.out == "my_prefix"
    assert args.inmat == "genmat_file"
    assert args.shapefile == "shapefile.shp"


def test_display_help_called():
    test_args = ["script_name", "-h"]
    sys.argv = test_args

    with patch('sys.exit') as mock_exit, patch('builtins.print') as mock_print:
        _ = parseArgs()
        mock_exit.assert_called()
        mock_print.assert_called()


def test_parse_args_conflicting_weight_options():
    test_args = [
        "script_name", "-o", "my_prefix", "-g", "genmat_file",
        "--shp", "shapefile.shp", "--fixWeight",
        "-v", "example1,example2,example3"]
    sys.argv = test_args

    with patch('sys.exit') as mock_exit, patch('builtins.print') as mock_print:
        _ = parseArgs()
        mock_exit.assert_called()
        mock_print.assert_called()


def test_parse_args_variable_file_input():
    # Create a temporary file and write some data to it
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file.write("var1\tSUM\nvar2")
        temp_file_name = temp_file.name

    # Set sys.argv to simulate command line arguments
    test_args = [
        "script_name", "-o", "my_prefix", "-g", "genmat_file",
        "--shp", "shapefile.shp", "--varFile", temp_file_name
    ]
    sys.argv = test_args

    # Run the function under test
    args = parseArgs()

    # Assertions
    assert args.varFile == temp_file_name
    assert args.agg_opts == {"var1": "SUM", "var2": "ARITH"}

    # Clean up the temporary file
    os.remove(temp_file_name)


def test_parse_args_defaults():
    test_args = [
        "script_name", "-g", "genmat_file",
        "--shp", "shapefile.shp", "-v", "example1"
    ]
    sys.argv = test_args

    # Create an instance of parseArgs
    args = parseArgs()

    # Assertions to check if defaults are set correctly
    assert args.GA_procs == 1
    assert args.deltaB is None
    assert args.deltaB_perc == 0.001
    assert args.nfail == 50
    assert args.maxGens == 500
    assert args.tournsize == 10
    assert args.cxpb == 0.5
    assert args.mutpb == 0.5
    assert args.indpb == 0.5
    assert args.burnin == 0
    assert args.maxpopsize == 1000
    assert args.max_hof_size == 100
    assert args.out == "output"
    assert args.awsum == 0.95
    assert args.modavg is True
    assert args.report_all is False
    assert args.plot is True
    assert args.fixWeight is False
    assert args.fixShape is False
    assert args.max_shape == 100
    assert args.min_weight == 0.0


def test_parse_args_nondefaults():
    test_args = [
        "script_name", "-o", "my_prefix", "-g", "genmat_file",
        "--shp", "shapefile.shp", "-v", "example1", "-P", "123",
        "-G", "123", "--threads", "4", "-d", "0.001", "-D", "0.01",
        "-F", "42", "-T", "42", "--cxpb", "0.42", "--indpb", "0.42",
        "--mutpb", "0.42", "--burn", "2", "--report_all",
        "-o", "example", "--awsum", "0.42", "--max_hof_size", "123",
        "--min_weight", "0.042", "--max_shape", "42"
    ]
    sys.argv = test_args

    # Create an instance of parseArgs
    args = parseArgs()

    # Assertions to check if defaults are set correctly
    assert args.GA_procs == 4
    assert args.deltaB == 0.001
    assert args.deltaB_perc == 0.01
    assert args.nfail == 42
    assert args.maxGens == 123
    assert args.tournsize == 42
    assert args.maxpopsize == 123
    assert args.cxpb == 0.42
    assert args.mutpb == 0.42
    assert args.indpb == 0.42
    assert args.burnin == 2
    assert args.max_hof_size == 123
    assert args.out == "example"
    assert args.awsum == 0.42
    assert args.modavg is True
    assert args.report_all is True
    assert args.plot is True
    assert args.fixWeight is False
    assert args.fixShape is False
    assert args.max_shape == 42
    assert args.min_weight == 0.042
