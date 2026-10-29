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
