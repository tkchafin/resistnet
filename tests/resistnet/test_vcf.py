import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

import autostreamtree
from autostreamtree.functions import population_distances, read_popmap

import resistnet.utils as utils
from resistnet.params import parseArgs

# autoStreamTree's test data: a VCF, a population map, the individuals'
# coordinates and a HydroRIVERS network with HydroATLAS covariates
AST = os.path.join(os.path.dirname(autostreamtree.__file__), "data")
VCF = os.path.join(AST, "test_sub100.vcf.gz")
POPMAP = os.path.join(AST, "test.popmap")
SITES = ["sing", "pipp", "dort", "kame", "pach", "thew"]


@pytest.fixture
def site_coords(tmp_path):
    # Sites at their individuals' mean coordinates, in a non-sorted order
    inds = pd.read_csv(os.path.join(AST, "test.coords"), sep="\t")
    pops = read_popmap(POPMAP)
    inds["site"] = inds["sample"].map(pops)
    xy = inds[inds["site"].isin(SITES)].groupby("site")[["lat", "long"]]
    xy = xy.mean().loc[SITES].reset_index().rename(columns={"site": "sample"})
    path = tmp_path / "sites.tsv"
    xy.to_csv(path, sep="\t", index=False)
    return str(path)


def test_vcf_distances_matches_autostreamtree(site_coords, tmp_path):
    # The matrix is autoStreamTree's for the sites' individuals, in the
    # order of the coordinates table
    out = str(tmp_path / "run")
    fname = utils.vcf_distances(VCF, POPMAP, site_coords, out,
                                dist="LINFST", concat="none")
    assert fname == f"{out}.popGenDistMat.txt"
    mat = pd.read_csv(fname, sep="\t", index_col=0)
    assert list(mat.index) == SITES and list(mat.columns) == SITES
    pops = {i: s for i, s in read_popmap(POPMAP).items() if s in SITES}
    ref = population_distances(VCF, pops, dist="LINFST", concat="none")
    np.testing.assert_allclose(mat.to_numpy(),
                               ref.loc[SITES, SITES].to_numpy(),
                               rtol=0, atol=1e-12)
    # Only the sites' individuals count: the same with all populations
    full = population_distances(VCF, POPMAP, dist="LINFST", concat="none")
    np.testing.assert_allclose(mat.to_numpy(),
                               full.loc[SITES, SITES].to_numpy(),
                               rtol=0, atol=1e-12)


def test_vcf_distances_site_without_individuals(site_coords, tmp_path):
    xy = pd.read_csv(site_coords, sep="\t")
    xy.loc[len(xy)] = ["nowhere", xy["lat"].iloc[0], xy["long"].iloc[0]]
    path = tmp_path / "more_sites.tsv"
    xy.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="nowhere"):
        utils.vcf_distances(VCF, POPMAP, str(path), str(tmp_path / "x"))


@pytest.mark.parametrize("args, message", [
    (["--vcf", VCF, "--popmap", POPMAP, "-g", "gen.tsv"], "can't both"),
    (["--vcf", VCF], "needs --popmap"),
])
def test_vcf_options_checked(args, message, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["runResistnet.py", "-c", "c.tsv",
                                      "-s", "net.shp", "-v", "x"] + args)
    with pytest.raises(SystemExit):
        parseArgs()
    assert message in capsys.readouterr().err


def test_run_with_vcf(site_coords, tmp_path):
    # runResistnet.py --vcf: the matrix is written to the output directory
    # and fitted, as if given with -g
    script = os.path.join(os.path.dirname(__file__), "..", "..", "scripts",
                          "runResistnet.py")
    out = str(tmp_path / "res" / "test")
    p = subprocess.run(
        [sys.executable, script, "-s", os.path.join(AST, "test.shp"),
         "-c", site_coords, "--vcf", VCF, "--popmap", POPMAP,
         "--gendist", "FST", "--vars", "dor_pc_pva", "-i", "3", "-r", "1",
         "-S", "2", "-t", "1", "--seed", "1", "-o", out],
        capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    fname = os.path.join(out, "test.popGenDistMat.txt")
    mat = pd.read_csv(fname, sep="\t", index_col=0)
    ref = population_distances(VCF, POPMAP, dist="FST")
    np.testing.assert_allclose(mat.to_numpy(),
                               ref.loc[SITES, SITES].to_numpy(),
                               rtol=0, atol=1e-12)
    pw = pd.read_csv(os.path.join(out, "resistance", "test.pairwise.tsv"),
                     sep="\t")
    assert len(pw) == len(SITES) * (len(SITES) - 1) // 2
    np.testing.assert_allclose(
        pw["genetic_distance"],
        [mat.loc[a, b] for a, b in zip(pw["site1"], pw["site2"])],
        rtol=1e-12)
