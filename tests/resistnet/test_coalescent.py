# The coalescent model (Lundgren & Ralph's coalescence times, computed by
# autostreamtree.coalescent, which tests the computation itself) on the
# NMT network
import os

import numpy as np
import pandas as pd
import pytest


def _coal_net(tmp_path, cell_km, density=1000.0, variables=("dor_pc_pva",)):
    import resistnet
    from resistnet.coalescent_network import ResistanceNetworkCoalescent
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    v = list(variables)
    return ResistanceNetworkCoalescent(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"), sizes=None,
        coords=os.path.join(data, "nmt.pointCoords.txt"), variables=v,
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={x: "ARITH" for x in v},
        out=str(tmp_path / f"coal{cell_km:g}"), verbose=False,
        cell_km=cell_km, density=density)


def test_coalescent_network_symmetric_is_river_distance(tmp_path):
    # With symmetric movement (no flow), resistance 1 per km and high
    # gene flow, linearised FST from the coalescence times is river
    # distance up to a shift and scale (coalescence and commute times are
    # equal when hitting times are symmetric; Lundgren & Ralph 2019), so
    # MLPE fits it as well as the standard model's distance alone
    net = _coal_net(tmp_path, 2.0, density=1e5)
    lin = net._lattice.linfst()
    from resistnet.resistance_network import ResistanceNetwork
    import resistnet
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    std = ResistanceNetwork(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"),
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=["dor_pc_pva"],
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={"dor_pc_pva": "ARITH"}, out=str(tmp_path / "std"),
        verbose=False)
    null = std.evaluate_null_model()
    res = net._flow_alone(0.0)
    assert float(res["loglik"].iloc[0]) == pytest.approx(
        float(null.loc["distance_only", "loglik"]), abs=0.05)
    assert lin.min() >= 0
    assert np.allclose(lin, lin.T)


def test_coalescent_network_spacing_and_direction(tmp_path):
    # The lattice's spacing barely changes the fit (finite-volume
    # movement and coalescence rates); a flow term changes it, and the
    # model's specification records the coalescent settings
    fits = {}
    for cell in (2.0, 5.0):
        net = _coal_net(tmp_path, cell)
        fits[cell] = [float(net._flow_alone(r)["loglik"].iloc[0])
                      for r in (0.0, -0.02)]
    np.testing.assert_allclose(fits[2.0], fits[5.0], atol=1.5)
    assert fits[2.0][1] != pytest.approx(fits[2.0][0])
    spec = net.model_spec([1, 1.0, 7, 1.0, 0, 50.0, -0.02])
    assert spec["type"] == "coalescent" and spec["density"] == 1000.0
    assert "absorption" not in spec


def test_run_coalescent_model(tmp_path):
    # runResistnet.py --model coalescent: its outputs, and the coalescent
    # reference models (distance alone, flow alone)
    import subprocess
    import sys
    import resistnet
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    script = os.path.join(os.path.dirname(__file__), "..", "..", "scripts",
                          "runResistnet.py")
    out = str(tmp_path / "run" / "nmt")
    p = subprocess.run(
        [sys.executable, script, "-n", os.path.join(data, "nmt.network"),
         "-c", os.path.join(data, "nmt.pointCoords.txt"),
         "-g", os.path.join(data, "nmt.popGenDistMat.txt"),
         "--vars", "dor_pc_pva", "--model", "coalescent", "-i", "3", "-r",
         "1", "-S", "2", "-t", "1", "--seed", "1", "-o", out],
        capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    prefix = os.path.join(out, "coalescent", "nmt")
    for ext in ("HallOfFame.tsv", "Model-Best.CoalescenceMatrix.tsv",
                "Model-Average.TransitionEdges.tsv", "modelSpec.json",
                "pairwise.tsv"):
        assert os.path.exists(f"{prefix}.{ext}"), ext
    null = pd.read_csv(f"{prefix}.Null-Model.tsv", sep="\t", index_col=0)
    assert {"coal_distance", "coal_flow"} <= set(null.index)
    pw = pd.read_csv(f"{prefix}.pairwise.tsv", sep="\t")
    assert "best_model_distance" in pw.columns
