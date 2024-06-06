import pytest
from unittest.mock import patch
import sys
from resistnet.params import parseArgs


def test_parseArgs_with_args():
    testargs = ["prog", "--seed", "42", "--inmat", "matrix.txt", "--shapefile", "file.shp", "--coords", "coords.tsv",
                "--network", "network.pkl", "--procs", "4", "--out", "results", "--gdf_out", "SHP",
                "--gridSearch", "--transFile", "trans.txt", "--vars", "var1,var2", "--varfile", "vars.txt",
                "--nfail", "60", "--max_iter", "600", "--reps", "20", "--pweight", "0.8", "--nstart", "30",
                "--ncand", "50", "--gamma", "0.2", "--fitmetric", "r2m", "--max_hof_size", "150",
                "--fixWeight", "--fixShape", "--max_shape", "200", "--min_weight", "0.1", "--awsum", "2.0",
                "--use_full", "--report_all"]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.seed == 42
        assert args.inmat == "matrix.txt"
        assert args.shapefile == "file.shp"
        assert args.coords == "coords.tsv"
        assert args.network == "network.pkl"
        assert args.procs == 4
        assert args.out == "results"
        assert args.gdf_out == "SHP"
        assert args.gridSearch
        assert args.transFile == "trans.txt"
        assert args.vars == "var1,var2"
        assert args.varfile == "vars.txt"
        assert args.nfail == 60
        assert args.max_iter == 600
        assert args.reps == 20
        assert args.pweight == 0.8
        assert args.nstart == 30
        assert args.ncand == 50
        assert args.gamma == 0.2
        assert args.fitmetric == "r2m"
        assert args.max_hof_size == 150
        assert args.fixWeight
        assert args.fixShape
        assert args.max_shape == 200
        assert args.min_weight == 0.1
        assert args.awsum == 2.0
        assert args.use_full
        assert args.report_all


def test_variables_from_vars():
    testargs = ["prog", "--vars", "var1,var2,var3"]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.variables == ["var1", "var2", "var3"]


def test_variables_from_varfile(tmp_path):
    varfile = tmp_path / "vars.txt"
    varfile.write_text("var1\nvar2\nvar3\n")
    testargs = ["prog", "--varfile", str(varfile)]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.variables == ["var1", "var2", "var3"]


def test_no_vars_or_varfile():
    testargs = ["prog"]
    with patch.object(sys, 'argv', testargs):
        with pytest.raises(SystemExit):
            parseArgs()


if __name__ == "__main__":
    pytest.main()
