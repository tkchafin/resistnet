import json
import os

import numpy as np
import pandas as pd
import pytest

import resistnet
import resistnet.predict as pr
from resistnet.resistance_network import ResistanceNetwork
from resistnet.model_optimisation import ModelRunnerTPE

DATA = os.path.join(os.path.dirname(resistnet.__file__), "data")
VARIABLES = ["run_mm_cyr", "tmp_dc_cyr", "sgr_dk_rav"]


@pytest.fixture(scope="module", params=["km", "reach"])
def net(request, tmp_path_factory):
    # Both formulations of the standard model: resistance per km (default)
    # and per reach (earlier versions)
    out = tmp_path_factory.mktemp("predict")
    return ResistanceNetwork(
        network=os.path.join(DATA, "nmt.network"), shapefile=None,
        coords=os.path.join(DATA, "nmt.pointCoords.txt"),
        variables=VARIABLES,
        inmat=os.path.join(DATA, "nmt.popGenDistMat.txt"),
        agg_opts={v: "ARITH" for v in VARIABLES}, reachid_col="HYRIV_ID",
        out=str(out / "net"), verbose=False, units=request.param)


def model(spec, net=None, max_resistance=37.5):
    """spec: {variable: (weight, transformation, shape)}; per km, with the
    maximum resistance per km after the variables"""
    ind = []
    for v in VARIABLES:
        ind += [1, *spec[v]] if v in spec else [0, 1.0, 0, 1]
    if net is not None and net.units == "km":
        ind.append(max_resistance)
    return ind


def raw_values(net):
    """Raw (aggregated) values, from the 0-10 predictors and ranges."""
    raw = net._predictors.copy()
    for v, (lo, hi) in net._predictor_range.items():
        raw[v] = lo + raw[v] / 10 * (hi - lo)
    return raw


@pytest.mark.parametrize("t,s", [(0, 1), (1, 3), (2, 5), (3, 4), (4, 6),
                                 (5, 2), (6, 8), (7, 3), (8, 2)])
def test_spec_reproduces_resistances(net, t, s):
    # The saved specification, applied to the training edges' raw values,
    # gives exactly the fitted model's edge resistances
    ind = model({"run_mm_cyr": (0.7, t, s), "sgr_dk_rav": (0.3, 5, 9)},
                net)
    spec = net.model_spec(ind)
    assert spec["type"] == "standard"
    assert spec["units"] == net.units
    expected = net.resistance_surface(ind).to_numpy(dtype=float)
    np.testing.assert_allclose(
        pr.standard_resistance(raw_values(net), spec), expected,
        rtol=1e-12, atol=1e-15)
    if net.units == "km":
        # From 1 (distance alone) to the maximum resistance per km
        assert spec["max_resistance"] == 37.5
        assert expected.min() >= 1.0 and expected.max() <= 37.5
        # Times the reaches' lengths: the fitted model's reach resistances
        reaches = raw_values(net)
        reaches.insert(0, "HYRIV_ID", reaches.index.to_numpy())
        reaches["LENGTH_KM"] = net._lengths
        table, _, _ = pr.predict_reaches(
            reaches, {"reachid_col": "HYRIV_ID", "models": [spec]})
        np.testing.assert_allclose(
            table["reach_resistance"],
            net.edge_resistances(ind).to_numpy(dtype=float),
            rtol=1e-12, atol=1e-15)


def test_spec_file_and_reach_prediction(net, tmp_path):
    # A run writes the specifications; applying the best one to a table of
    # reaches clamps values outside the training range
    runner = ModelRunnerTPE(net, seed=2, verbose=False)
    out = str(tmp_path / "run")
    runner.run_tpe(max_evals=3, threads=1, reps=2, max_shape=10,
                   n_startup=2, out=out, verbose=False, plot=False)
    spec = pr.load_spec(f"{out}.modelSpec.json")
    assert spec["reachid_col"] == "HYRIV_ID"
    best = spec["models"][0]
    assert best["rank"] == 1
    raw = raw_values(net)
    reaches = raw.copy()
    reaches.insert(0, "HYRIV_ID", raw.index.to_numpy())
    table, moves, clamped = pr.predict_reaches(reaches, spec)
    assert moves is None
    assert all(v == 0 for v in clamped.values())
    hof = pd.read_csv(f"{out}.HallOfFame.tsv", sep="\t", na_values="-")
    ind = runner._row_individual(runner.bests.getHOF().iloc[0])
    np.testing.assert_allclose(
        table["resistance"],
        net.resistance_surface(ind).to_numpy(dtype=float),
        rtol=1e-12, atol=1e-15)
    assert len(spec["models"]) == len(hof)
    # Out-of-range reaches get the value at the edge of the range
    for v, var in best["variables"].items():
        far = reaches.copy()
        far[v] = np.where(far[v] == far[v].max(), far[v] + 100, far[v])
        far_table, _, far_clamped = pr.predict_reaches(far, spec)
        np.testing.assert_allclose(far_table["resistance"],
                                   table["resistance"])
        assert far_clamped[v] > 0
    json.dumps(spec)


def test_composite_clamped_and_reported():
    # Each variable is within its training range, but their combination
    # can still fall outside the training composite range (the variables'
    # minima didn't occur together): it is clamped, and reported
    var = {"min": 0.0, "max": 10.0, "transform": 0, "shape": 0.0,
           "tmin": 0.0, "tmax": 10.0, "weight": 1.0}
    model = {"type": "standard", "variables": {"a": dict(var),
                                               "b": dict(var)},
             "composite_min": 0.5, "composite_max": 1.5}
    spec = {"reachid_col": "HYRIV_ID", "models": [model]}
    reaches = pd.DataFrame({"HYRIV_ID": [1, 2, 3, 4],
                            "a": [0.0, 1.0, 10.0, 5.0],
                            "b": [0.0, 1.0, 5.0, 5.0]})
    table, moves, clamped = pr.predict_reaches(reaches, spec)
    assert moves is None
    assert clamped["a"] == 0 and clamped["b"] == 0
    # Combined values 0, 0.2, 1.5, 1.0 -> (x - 0.5) / 1, clamped to 0-1
    np.testing.assert_allclose(table["resistance"], [0, 0, 1, 0.5])
    assert clamped["composite"] == 0.5

    # SAMC: the same for each move (here two reaches, one link each way)
    samc = dict(model, type="samc", max_resistance=10.0, flow=0.0)
    for v in samc["variables"].values():
        v.update(asym=0, link_min=0.0, link_max=1.0)
    reaches = pd.DataFrame({"HYRIV_ID": [1, 2], "NEXT_DOWN": [2, 0],
                            "LENGTH_KM": [1.0, 1.0],
                            "a": [0.0, 1.0], "b": [0.0, 1.0]})
    spec = {"reachid_col": "HYRIV_ID", "length_col": "LENGTH_KM",
            "models": [samc]}
    table, moves, clamped = pr.predict_reaches(reaches, spec)
    # Links (mean of 0 and 0.1 on the 0-1 scale) sum to 0.1 < 0.5: clamped
    assert clamped["composite"] == 1.0
    np.testing.assert_allclose(moves["weight"], 1.0)
