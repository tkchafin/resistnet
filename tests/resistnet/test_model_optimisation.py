import os

import numpy as np
import pandas as pd
import pytest
import tempfile
import resistnet
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunnerTPE


@pytest.fixture
def tpe_parameters_fixture():
    return {
        'fixWeight': False,
        'fixShape': True,
        'min_weight': 0.5,
        'max_shape': 2.0,
        'fitmetric': 'aic',
        'awsum': 0.75,
        'only_keep': True,
        'verbose': False,
        'use_full': False,
        'report_all': False
    }


@pytest.fixture
def network_graph_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'nmt.network')
    return file_path


@pytest.fixture
def shapefile_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'nmt.gpkg')
    return file_path


@pytest.fixture
def coords_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'nmt.pointCoords.txt')
    return file_path


@pytest.fixture
def inmat_path():
    base_path = os.path.dirname(resistnet.__file__)
    file_path = os.path.join(base_path, 'data', 'nmt.popGenDistMat.txt')
    return file_path


@pytest.fixture
def resistance_network_fixture(
        network_graph_path, shapefile_path, coords_path, inmat_path):
    # Use temporary directory for output
    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "output")

        network = ResistanceNetwork(
            network=network_graph_path,
            shapefile=shapefile_path,
            coords=coords_path,
            variables=["run_mm_cyr"],
            inmat=inmat_path,
            agg_opts={"run_mm_cyr": "ARITH"},
            pop_agg="ARITH",
            reachid_col="HYRIV_ID",
            length_col="LENGTH_KM",
            out=output_prefix,
            verbose=False
        )

        yield network


@pytest.fixture
def model_runner_fixture(resistance_network_fixture):
    seed = 1234
    verbose = True

    model_runner = ModelRunnerTPE(
        resistance_network=resistance_network_fixture,
        seed=seed, verbose=verbose)

    yield model_runner


def test_model_runner_initialization(resistance_network_fixture):
    seed = 1234
    verbose = True

    model_runner = ModelRunnerTPE(
        resistance_network=resistance_network_fixture,
        seed=seed, verbose=verbose)

    # Check if the passed parameters are correctly assigned
    assert model_runner.seed == seed
    assert model_runner.resistance_network == resistance_network_fixture
    assert model_runner.verbose == verbose

    # Check if the default values are correctly initialized
    assert isinstance(model_runner.workers, list)
    assert model_runner.bests is None
    assert isinstance(model_runner.logger, list)

    # Check the default values of GA parameters
    assert model_runner.fixWeight is None
    assert model_runner.fixShape is None
    assert model_runner.min_weight is None
    assert model_runner.max_shape is None
    assert model_runner.fitmetric is None
    assert model_runner.awsum is None
    assert model_runner.only_keep is None
    assert model_runner.report_all is None


def test_set_tpe_parameters(
        model_runner_fixture, tpe_parameters_fixture):
    # Call set_ga_parameters on the fixture instance with the test parameters
    model_runner_fixture.set_tpe_parameters(
        **tpe_parameters_fixture
    )

    # Assert that each parameter is set correctly
    for param, value in tpe_parameters_fixture.items():
        assert getattr(model_runner_fixture, param) == \
            value, f"Parameter {param} not set correctly"


def test_run_tpe(model_runner_fixture, tpe_parameters_fixture):
    with tempfile.TemporaryDirectory() as temp_dir:
        out_prefix = os.path.join(temp_dir, "tpe_output")
        tpe_parameters_fixture["verbose"] = False
        tpe_parameters_fixture["out"] = out_prefix
        tpe_parameters_fixture["threads"] = 1
        tpe_parameters_fixture["max_evals"] = 5

        # Run the genetic algorithm with the modified parameters
        try:
            model_runner_fixture.run_tpe(**tpe_parameters_fixture)
            error_occurred = False
        except Exception as e:
            error_occurred = True
            print(f"Error during TPE run: {e}")

        # Assert that no errors occurred during the run
        assert not error_occurred, "An error occurred during the TPE run"

        # Assert that the best models are identified (bests is not None)
        assert model_runner_fixture.bests is not None, \
            "ModelRunner.bests missing"


def test_nfail_counts_after_startup(model_runner_fixture, monkeypatch):
    # The random start-up iterations never count towards nFail, so at least
    # n_startup + 1 iterations run; the workers' argument file is removed
    seen = {}
    orig = ModelRunnerTPE.process_hall_of_fame

    def keep(self, trials_list):
        seen["n"] = [len(t.trials) for t in trials_list]
        seen["args_file"] = self._worker_args_file
        seen["exists"] = os.path.exists(self._worker_args_file)
        return orig(self, trials_list)
    monkeypatch.setattr(ModelRunnerTPE, "process_hall_of_fame", keep)
    model_runner_fixture.run_tpe(
        max_evals=12, threads=1, reps=1, n_startup=6, nFail=1,
        fixShape=True, max_shape=2, verbose=False, out=None)
    assert seen["n"][0] >= 7
    assert seen["exists"]
    assert not os.path.exists(seen["args_file"])
    assert model_runner_fixture._worker_args_file is None


RUN_SCRIPT = """
import os, sys
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunnerTPE

if __name__ == "__main__":
    data, out = sys.argv[1], sys.argv[2]
    net = ResistanceNetwork(
        network=os.path.join(data, "nmt.network"), shapefile=None,
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=["run_mm_cyr"],
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={"run_mm_cyr": "ARITH"}, reachid_col="HYRIV_ID",
        out=out, verbose=False)
    if sys.argv[3] == "missing_network":
        net.network = out + ".missing.net"
    runner = ModelRunnerTPE(net, seed=1, verbose=False)
    try:
        runner.run_tpe(max_evals=4, threads=2, reps=2, max_shape=10,
                       verbose=False, out=None)
        print("FINISHED")
    except ValueError as e:
        print("ERROR", e)
"""


def _run_in_subprocess(tmp_path, mode):
    import subprocess
    import sys
    script = tmp_path / "run.py"
    script.write_text(RUN_SCRIPT)
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    # A timeout turns a hang into a failure
    p = subprocess.run(
        [sys.executable, str(script), data, str(tmp_path / "net"), mode],
        cwd=cwd, capture_output=True, text=True, timeout=300)
    return p.stdout, sorted(os.listdir(cwd))


def test_worker_start_failure_does_not_hang(tmp_path):
    stdout, _ = _run_in_subprocess(tmp_path, "missing_network")
    assert "ERROR No valid models found" in stdout


def test_run_without_output_writes_no_files(tmp_path):
    stdout, files = _run_in_subprocess(tmp_path, "ok")
    assert "FINISHED" in stdout
    assert files == []


def test_shape_search_space(model_runner_fixture):
    # Shapes are integers from 1 to max_shape, drawn log-uniformly (most
    # often small, where the transformations differ most)
    import hyperopt.pyll.stochastic as stoch
    runner = model_runner_fixture
    runner.set_tpe_parameters(max_shape=15)
    runner.initialize_space()
    name = [k for k in runner.space if k.endswith("_shape")][0]
    rng = np.random.default_rng(1)
    draws = np.array([stoch.sample(runner.space[name], rng=rng)
                      for _ in range(2000)])
    assert np.all(draws == np.round(draws))
    assert draws.min() >= 1 and draws.max() <= 15
    assert np.median(draws) <= 5


def test_max_resistance_per_km(model_runner_fixture, tmp_path):
    # Standard model per km: the maximum resistance per km is searched
    # log-uniformly from 1 to max_resist, kept in the Hall of Fame (whose
    # AIC is the model's, re-evaluated) and the best model's table
    import hyperopt.pyll.stochastic as stoch
    runner = model_runner_fixture
    net = runner.resistance_network
    assert runner.per_km and net.max_resist == 1000
    runner.set_tpe_parameters(max_shape=15)
    runner.initialize_space()
    rng = np.random.default_rng(1)
    draws = np.array([stoch.sample(runner.space["max_resistance"], rng=rng)
                      for _ in range(2000)])
    assert draws.min() >= 1 and draws.max() <= 1000
    assert np.median(draws) == pytest.approx(np.sqrt(1000), rel=0.15)

    out = str(tmp_path / "km")
    runner.run_tpe(max_evals=6, threads=1, reps=2, max_shape=15,
                   n_startup=3, out=out, verbose=False, plot=False)
    hof = pd.read_csv(f"{out}.HallOfFame.tsv", sep="\t", na_values="-")
    assert "max_resistance" in hof.columns
    assert hof["max_resistance"].between(1, 1000).all()
    best = runner.bests.getHOF().iloc[0]
    ind = runner._row_individual(best)
    assert ind[-1] == best["max_resistance"]
    assert -net.evaluate(ind)[1][2] == pytest.approx(float(best["aic"]),
                                                     abs=1e-9)
    table = pd.read_csv(f"{out}.bestModel.tsv", sep="\t", na_values="-")
    rows = table.set_index("Variable")
    assert rows.loc["max_resistance", "Value"] == pytest.approx(
        best["max_resistance"])
    assert rows["Contribution"].sum() == pytest.approx(1.0, abs=1e-9)
