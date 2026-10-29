import pytest
from unittest.mock import patch
import sys
from resistnet.params import parseArgs


def test_parseArgs_with_args():
    testargs = ["prog", "--seed", "42", "--inmat", "matrix.txt", "--shapefile", "file.shp", "--coords", "coords.tsv",
                "--network", "network.pkl", "--procs", "4", "--out", "results",
                "--gridSearch", "--transFile", "trans.txt", "--vars", "var1,var2", "--varfile", "vars.txt",
                "--nfail", "60", "--max_iter", "600", "--reps", "20", "--pweight", "0.8", "--nstart", "30",
                "--ncand", "50", "--gamma", "0.2", "--fitmetric", "r2m",
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


def test_varfile_aggregators(tmp_path):
    varfile = tmp_path / "vars.txt"
    varfile.write_text("var1\tSUM\nvar2\nvar3\tMAX\n")
    testargs = ["prog", "--varfile", str(varfile)]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.variables == ["var1", "var2", "var3"]
        assert args.agg_opts == {
            "var1": "SUM", "var2": "ARITH", "var3": "MAX"
        }


def test_fitmetric_default_aic():
    testargs = ["prog", "--vars", "var1"]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.fitmetric == "aic"
        assert args.nfail == 200


def test_samc_defaults():
    testargs = ["prog", "--vars", "var1", "--samc"]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.samc
        assert args.sizes is None
        assert args.ne == 100
        assert args.max_resist is None
        assert args.max_resist_for("samc") == 500
        assert args.max_resist_for("resistance") == 1000
        assert args.fix_resist is None
        assert args.max_flow == 1.0
        assert not args.fixFlow
        assert not args.fixAsym


@pytest.mark.parametrize("opt", [["--max_resist", "1"],
                                 ["--fix_resist", "0.5"],
                                 ["--max_flow", "0.001"]])
def test_samc_resistance_limits(opt):
    testargs = ["prog", "--vars", "var1", "--samc"] + opt
    with patch.object(sys, 'argv', testargs):
        with pytest.raises(SystemExit):
            parseArgs()


def test_fitmetric_delta():
    testargs = ["prog", "--vars", "var1", "-f", "delta"]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.fitmetric == "delta_aic_null"
        assert args.procs == 1


def test_no_vars_or_varfile():
    testargs = ["prog"]
    with patch.object(sys, 'argv', testargs):
        with pytest.raises(SystemExit):
            parseArgs()


if __name__ == "__main__":
    pytest.main()


@pytest.mark.parametrize("opts,model,models", [
    ([], "resistance", ["resistance"]),
    (["--samc"], "samc", ["samc"]),
    (["--model", "samc"], "samc", ["samc"]),
    (["--model", "both"], "both", ["resistance", "samc"]),
    (["--model", "samc", "--samc"], "samc", ["samc"]),
])
def test_model_option(opts, model, models):
    testargs = ["prog", "--vars", "var1"] + opts
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.model == model
        assert args.models == models


def test_model_conflicts_with_samc():
    testargs = ["prog", "--vars", "var1", "--samc", "--model", "both"]
    with patch.object(sys, 'argv', testargs):
        with pytest.raises(SystemExit):
            parseArgs()


def test_units_and_resistance_limits():
    # Standard model: resistance per km by default (--per_reach: earlier
    # formulation); --max_resist applies to both models when given
    testargs = ["prog", "--vars", "var1"]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert not args.per_reach
    testargs = ["prog", "--vars", "var1", "--per_reach", "--max_resist",
                "50", "--model", "both"]
    with patch.object(sys, 'argv', testargs):
        args = parseArgs()
        assert args.per_reach
        assert args.max_resist_for("resistance") == 50
        assert args.max_resist_for("samc") == 50
