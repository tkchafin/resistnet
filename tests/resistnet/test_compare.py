import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

import resistnet
import resistnet.MLPE as mlpe
from resistnet.compare import comparison_tables, write_comparison
from resistnet.model_optimisation import ModelRunnerTPE
from resistnet.resistance_network import ResistanceNetwork
from resistnet.samc_network import ResistanceNetworkSAMC

DATA = os.path.join(os.path.dirname(resistnet.__file__), "data")
VARIABLES = ["sgr_dk_rav", "ele_mt_cav", "dor_pc_pva"]


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    # Both model families fitted (briefly) to the NMT genetic distances
    tmp = tmp_path_factory.mktemp("compare")
    kw = dict(network=None, shapefile=os.path.join(DATA, "nmt.gpkg"),
              coords=os.path.join(DATA, "nmt.pointCoords.txt"),
              variables=VARIABLES,
              inmat=os.path.join(DATA, "nmt.popGenDistMat.txt"),
              agg_opts={v: "ARITH" for v in VARIABLES}, verbose=False)
    out = {}
    for family in ("resistance", "samc"):
        prefix = str(tmp / family / "nmt")
        os.makedirs(os.path.dirname(prefix))
        net = (ResistanceNetworkSAMC(sizes=None, out=prefix, **kw)
               if family == "samc" else ResistanceNetwork(out=prefix, **kw))
        runner = ModelRunnerTPE(net, seed=3, verbose=False)
        runner.run_tpe(max_evals=6, threads=1, reps=3, max_shape=10,
                       n_startup=3, out=prefix, verbose=False, plot=False)
        out[family] = (net, runner, prefix)
    return out


def test_comparison_tables(runs):
    table, summary, importance = comparison_tables(runs)
    # One ranking by AIC, with joint Akaike weights
    assert np.all(np.diff(table["aic"]) >= 0)
    assert list(table["rank"]) == list(range(1, len(table) + 1))
    assert table["akaike_weight"].sum() == pytest.approx(1.0)
    np.testing.assert_allclose(table["delta_aic"],
                               table["aic"] - table["aic"].min())
    # Each reference model once, with its number of parameters
    refs = table[table["family"] == "reference"].set_index("model")["k"]
    assert len(refs) == 4
    assert refs["population effects only (null)"] == 3
    assert refs["distance only (km; standard model)"] == 4
    assert refs["SAMC distance alone"] == 4
    # Each model once per family
    assert not table.duplicated(["family", "model", "aic"]).any()
    # k is the MLPE parameters plus the model's surface parameters
    for family, (net, runner, prefix) in runs.items():
        hof = pd.read_csv(f"{prefix}.HallOfFame.tsv", sep="\t",
                          na_values="-")
        ks = set(table.loc[table["family"] == family, "k"])
        expected = {4 + net.surface_parameters(runner._row_individual(r))
                    for _, r in hof.iterrows()}
        assert ks == expected
    # The summary sums each family's weights
    np.testing.assert_allclose(summary["akaike_weight"].sum(), 1.0)
    for _, row in summary.iterrows():
        part = table[table["family"] == row["family"]]
        assert row["best_aic"] == part["aic"].min()
        assert row["akaike_weight"] == pytest.approx(
            part["akaike_weight"].sum())
    assert set(importance["variable"]) <= set(VARIABLES)
    assert (importance["RVI"] <= 1 + 1e-12).all()


def test_comparison_requires_same_data(runs, tmp_path):
    # Families fitted to different genetic distances can't be compared
    net, runner, prefix = runs["samc"]
    moved = str(tmp_path / "nmt")
    for ext in ("HallOfFame.tsv", "Null-Model.tsv"):
        table = open(f"{prefix}.{ext}").read()
        if ext == "Null-Model.tsv":
            table = table.replace("-386.0", "-380.0")
        with open(f"{moved}.{ext}", "w") as fh:
            fh.write(table)
    with pytest.raises(ValueError, match="same genetic distances"):
        comparison_tables({"resistance": runs["resistance"],
                           "samc": (net, runner, moved)})


def test_write_comparison(runs, tmp_path):
    prefix = str(tmp_path / "nmt")
    write_comparison(runs, prefix, verbose=False)
    for ext in ("modelComparison", "familySummary", "jointImportance"):
        assert os.path.exists(f"{prefix}.{ext}.tsv")
    for ext in ("modelComparison", "jointImportance", "bestFits"):
        assert os.path.getsize(f"{prefix}.{ext}.pdf") > 0
    # Each variable's importance is the sum of its parts from each family
    imp = pd.read_csv(f"{prefix}.jointImportance.tsv", sep="\t")
    np.testing.assert_allclose(imp["RVI"], imp["RVI_resistance"]
                               + imp["RVI_samc"])


def test_run_both_models_layout(tmp_path):
    # runResistnet.py --model both writes each model's outputs to its own
    # subdirectory, and the comparison above them
    script = os.path.join(os.path.dirname(__file__), "..", "..", "scripts",
                          "runResistnet.py")
    out = str(tmp_path / "res" / "nmt")
    p = subprocess.run(
        [sys.executable, script, "-n", os.path.join(DATA, "nmt.network"),
         "-c", os.path.join(DATA, "nmt.pointCoords.txt"),
         "-g", os.path.join(DATA, "nmt.popGenDistMat.txt"),
         "--vars", "dor_pc_pva", "--model", "both", "-i", "4", "-r", "2",
         "-S", "2", "-t", "1", "--seed", "1", "-o", out],
        capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    for model in ("resistance", "samc"):
        assert os.path.exists(
            os.path.join(out, model, "nmt.HallOfFame.tsv"))
    assert os.path.exists(
        os.path.join(out, "samc", "nmt.Model-Best.CFPTMatrix.tsv"))
    for ext in ("modelComparison", "familySummary", "jointImportance"):
        assert os.path.exists(os.path.join(out, f"nmt.{ext}.tsv"))

    # The pairwise tables: genetic distances as in the input, and model
    # distances as in the written matrices (SAMC: mean of both directions)
    gen = pd.read_csv(os.path.join(DATA, "nmt.popGenDistMat.txt"),
                      sep="\t", index_col=0)
    for model, kind in (("resistance", "ResistanceMatrix"),
                        ("samc", "CFPTMatrix")):
        prefix = os.path.join(out, model, "nmt")
        pw = pd.read_csv(f"{prefix}.pairwise.tsv", sep="\t")
        assert len(pw) == 14 * 13 // 2
        np.testing.assert_allclose(
            pw["genetic_distance"],
            [gen.loc[a, b] for a, b in zip(pw["site1"], pw["site2"])],
            rtol=1e-12)
        M = pd.read_csv(f"{prefix}.Model-Best.{kind}.tsv", sep="\t",
                        index_col=0)
        M.columns = M.columns.astype(str)
        np.testing.assert_allclose(
            pw["best_model_distance"],
            [(M.loc[a, b] + M.loc[b, a]) / 2
             for a, b in zip(pw["site1"], pw["site2"])], rtol=1e-12)
        # River distance is the distance-only reference model's distance:
        # MLPE on it gives that model's log-likelihood
        n = 14
        iu = np.triu_indices(n, 1)
        sites = list(dict.fromkeys(list(pw["site1"]) + list(pw["site2"])))
        pos = {s: i for i, s in enumerate(sites)}
        X, Y = np.zeros((n, n)), np.zeros((n, n))
        for a, b, km, g in zip(pw["site1"], pw["site2"], pw["river_km"],
                               pw["genetic_distance"]):
            i, j = pos[a], pos[b]
            X[i, j] = X[j, i] = km
            Y[i, j] = Y[j, i] = g
        assert np.all(X[iu] > 0)
        null = pd.read_csv(f"{prefix}.Null-Model.tsv", sep="\t",
                           index_col=0, keep_default_na=False,
                           na_values=[""])
        assert mlpe.MLPE(Y, X, scale=True)["loglik"].iloc[0] == \
            pytest.approx(float(null.loc["distance_only", "loglik"]),
                          abs=1e-6)

    # The interactive report, with both models and their comparison
    with open(os.path.join(out, "nmt.report.html")) as fh:
        page = fh.read()
    for section in ('id="comparison"', 'id="resistance"', 'id="samc"'):
        assert section in page
    assert page.count('class="plotly-graph-div"') >= 10

    # ensembleResistnet.py finds runs in the model subdirectories; with
    # both models there, --model picks one
    ens = os.path.join(os.path.dirname(script), "ensembleResistnet.py")
    args = [sys.executable, ens, "-n", os.path.join(DATA, "nmt.network"),
            "-i", str(tmp_path / "res"), "-c", "HYRIV_ID",
            "-C", os.path.join(DATA, "nmt.pointCoords.txt"), "-t", "1", "-X"]
    p = subprocess.run(args + ["-o", str(tmp_path / "mix")],
                       capture_output=True, text=True)
    assert p.returncode != 0
    assert "--model" in p.stderr + p.stdout
    p = subprocess.run(args + ["-o", str(tmp_path / "ens"), "--model",
                               "samc"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    assert os.path.exists(str(tmp_path / "ens.Model-Best.CFPTMatrix.tsv"))
