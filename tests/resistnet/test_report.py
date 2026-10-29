import os

import pandas as pd

from resistnet.report import RunOutputs, describe


def test_describe():
    row = pd.Series({"a": 1, "a_weight": 0.7, "a_trans": 5, "a_shape": 3.0,
                     "a_asym": 1, "b": 1, "b_weight": 0.3, "b_trans": 0,
                     "b_shape": 1.0, "b_asym": 0, "c": 0, "flow": -0.05})
    assert describe(row, ["a", "b", "c"], samc=True) == (
        "a (monomolecular, shape 3, directional, weight 0.70) + "
        "b (linear, symmetric, weight 0.30) + flow -0.050/km")
    assert describe(row, ["b"], samc=False) == "b (linear)"


def test_sites_with_numpy_scalars(tmp_path):
    # Site nodes can be written with numpy's scalar reprs
    prefix = str(tmp_path / "x")
    pd.DataFrame({"fitness": [1.0], "a": [1], "a_weight": [1.0],
                  "a_trans": [0], "a_shape": [1.0], "aic": [-1.0],
                  "r2m": [0.5], "akaike_weight": [1.0], "keep": [True]}
                 ).to_csv(f"{prefix}.HallOfFame.tsv", sep="\t", index=False)
    pd.DataFrame({
        "Label": ["s1", "s2"],
        "Node": ["(-96.5, 38.25)",
                 "(np.float64(-95.03541666666736), np.float64(36.98125))"],
        "Index": [0, 1], "Names": ["s1", "s2"]}
    ).to_csv(f"{prefix}.pointsTable.txt", sep="\t", index=False)
    sites = RunOutputs(prefix, "resistance").sites()
    assert list(sites["lon"]) == [-96.5, -95.03541666666736]
    assert list(sites["lat"]) == [38.25, 36.98125]
    assert os.path.exists(f"{prefix}.HallOfFame.tsv")


def test_describe_per_km():
    # Standard model per km: the range of the resistance per km
    row = pd.Series({"a": 1, "a_weight": 1.0, "a_trans": 7, "a_shape": 1.0,
                     "max_resistance": 144.6})
    assert describe(row, ["a"], samc=False) == (
        "a (inverse monomolecular, shape 1); resistance 1-145 per km")


def test_response_curves_match_fitted_model():
    # The best model's response to its variable, from modelSpec.json: on
    # the bundled example (one variable), each reach's point is its fitted
    # resistance per km (Model-Best.ResistanceEdges.tsv)
    import numpy as np
    import resistnet
    from resistnet.report import fig_response
    prefix = os.path.join(os.path.dirname(resistnet.__file__), "data",
                          "example_output", "resistance", "nmt")
    run = RunOutputs(prefix, "resistance")
    assert run.per_km
    fig = fig_response(run)
    curve, dots = fig.data
    assert curve.y.min() >= 1.0
    edges = pd.read_csv(f"{prefix}.Model-Best.ResistanceEdges.tsv",
                        sep="\t", index_col=0)["Resistance"]
    fitted = edges.loc[[int(r) for r in dots.customdata]].to_numpy()
    np.testing.assert_allclose(dots.y, fitted, rtol=1e-9)


def test_log_ticks():
    from resistnet.report import _log_ticks
    assert _log_ticks(1.0, 145.0) == [1, 2, 5, 10, 20, 50, 100]
    assert len(_log_ticks(74.0, 78.0)) >= 2


def test_contributions_figure():
    # Per km: distance's and each variable's share of the best model's
    # resistance per km (bestModel.tsv), which add up to 1
    import resistnet
    from resistnet.plotting import COLOURS
    from resistnet.report import fig_contributions
    prefix = os.path.join(os.path.dirname(resistnet.__file__), "data",
                          "example_output", "resistance", "nmt")
    fig = fig_contributions(RunOutputs(prefix, "resistance"))
    bar = fig.data[0]
    assert "distance" in bar.y
    assert abs(sum(bar.x) - 1) < 1e-9
    assert bar.marker.color[list(bar.y).index("distance")] == \
        COLOURS["reference"]


def test_duplicates_per_km_merged(tmp_path):
    # Chains converging on the same model per km find slightly different
    # maximum resistances: listed once, with the chains counted and their
    # Akaike weights summed; a different model stays separate
    prefix = str(tmp_path / "x")
    pd.DataFrame({
        "fitness": [-503.57177, -503.57172, -495.6],
        "a": [1, 1, 1], "a_weight": [1.0, 1.0, 1.0], "a_trans": [7, 7, 7],
        "a_shape": [1.0, 1.0, 2.0],
        "max_resistance": [122.4, 123.6, 35.7],
        "loglik": [258.79, 258.79, 254.8],
        "aic": [-503.57177, -503.57172, -495.6], "r2m": [0.77] * 3,
        "akaike_weight": [0.49, 0.49, 0.02], "keep": [True] * 3}
    ).to_csv(f"{prefix}.HallOfFame.tsv", sep="\t", index=False)
    models = RunOutputs(prefix, "resistance").models
    assert len(models) == 2
    assert models.iloc[0]["chains"] == 2
    assert abs(models.iloc[0]["akaike_weight"] - 0.98) < 1e-12
