import itertools
import os
import pickle
import pytest
import tempfile
import resistnet
import networkx as nx
import pandas as pd
import numpy as np
from unittest.mock import patch
from resistnet.resistance_network import ResistanceNetwork
from resistnet.resistance_network import ResistanceNetworkWorker
from resistnet.resistance_network import SimResistanceNetwork
import resistnet.utils as utils


@pytest.fixture
def fake_spec_file():
    spec_content = pd.DataFrame({
        "VAR": ["LENGTH_KM"],
        "WEIGHT": [1.0],
        "TRANSFORM": [0],
        "SHAPE": [0]
    })
    with tempfile.NamedTemporaryFile(
            mode='w+', suffix='.tsv', delete=False) as tmpfile:
        spec_content.to_csv(tmpfile.name, sep='\t', index=False)
        yield tmpfile.name
        os.remove(tmpfile.name)


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
def sim_resistance_network_params(network_graph_path):
    return {
        "network": network_graph_path,
        "reachid_col": "HYRIV_ID",
        "length_col": "LENGTH_KM",
        "verbose": True
    }


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


def test_initialize_predictors(resistance_network_fixture):
    resistance_network_fixture.initialize_predictors()

    # Assertions
    # Ensure that predictors are initialized and variables are updated
    assert resistance_network_fixture._predictors is not None
    assert len(resistance_network_fixture.variables) > 0


def test_read_network_with_existing_file(
        resistance_network_fixture, monkeypatch):
    # Mock the utility function that reads the pickle file
    mock_network = nx.Graph()
    monkeypatch.setattr(
        utils, 'read_pickle_network', lambda file: mock_network)

    # Setup: Assume network file exists
    resistance_network_fixture.network = 'existing_network_file'

    resistance_network_fixture.read_network()

    # Assertions
    assert resistance_network_fixture._G == mock_network


def test_evaluate_null_model(resistance_network_fixture):
    mock_df = pd.DataFrame({
        resistance_network_fixture.reachid_col: [1, 2, 3],
        resistance_network_fixture.length_col: [10, 20, 30]
    })

    # Create a mock result DataFrame to be returned by parsePairwise
    mock_parse_result = pd.DataFrame({
        "aic": [-100],
        "aic_null": [-200],
        "loglik_null": [-300],
        "loglik": [-300],
        "delta_aic_null": [100],
        "r2m": [0.5]
    }, index=[1])

    # Mocking external functions and methods used in evaluate_null_model
    with patch('resistnet.utils.nx_to_df', return_value=mock_df), \
         patch('resistnet.resist_dist.parsePairwise',
               return_value=(None, mock_parse_result)):
        # Call evaluate_null_model
        result = resistance_network_fixture.evaluate_null_model()

    assert isinstance(result, pd.DataFrame)
    assert 'loglik' in result.columns
    assert 'aic' in result.columns
    assert 'r2m' in result.columns
    assert 'delta_aic_null' in result.columns
    assert result.loc['distance_only', 'loglik'] == -300
    assert result.loc['distance_only', 'aic'] == 100


def test_model_output_mock(resistance_network_fixture):
    # Mock model representing an individual in the genetic algorithm
    model = [1, 1, 1, 1, 1]

    # Mock effective resistance matrix
    mock_effective_resistance_matrix = np.array([[0.5, 1.0], [1.0, 0.5]])

    # Mock the effectiveResistanceMatrix function
    with patch(
        'resistnet.resist_dist.effectiveResistanceMatrix',
        return_value=mock_effective_resistance_matrix) \
            as mock_effectiveResistanceMatrix:
        r, multi = resistance_network_fixture.model_output(model)

        mock_effectiveResistanceMatrix.assert_called_once()


def test_model_output(resistance_network_fixture):
    # Mock model representing an individual in the genetic algorithm
    model = [1, 1, 1, 1, 1]

    try:
        r, multi = resistance_network_fixture.model_output(model)
        exception_raised = False
    except Exception as e:
        exception_raised = True
        # Optionally, you can print the exception or take some action
        print(f"Exception occurred: {e}")

    # Assert no exception was raised
    assert not exception_raised, \
        "An exception was raised during the model_output execution"

    # Assert that r is a numpy.ndarray
    assert isinstance(r, np.ndarray), \
        "Expected r to be a numpy.ndarray"

    # Assert that multi is a pandas.Series
    assert isinstance(multi, pd.Series), \
        "Expected multi to be a pandas.Series"


def test_resistance_network_worker_creation(resistance_network_fixture):
    # Prepare additional arguments with defaults or extracted from the fixture
    worker_args = {
        'network': resistance_network_fixture.network,
        'pop_agg': resistance_network_fixture.pop_agg,
        'reachid_col': resistance_network_fixture.reachid_col,
        'length_col': resistance_network_fixture.length_col,
        'variables': resistance_network_fixture.variables,
        'agg_opts': resistance_network_fixture.agg_opts,
        'inc': resistance_network_fixture._inc,
        'point_coords': resistance_network_fixture._point_coords,
        'points_names': resistance_network_fixture._points_names,
        'points_snapped': resistance_network_fixture._points_snapped,
        'points_labels': resistance_network_fixture._points_labels,
        'predictors': resistance_network_fixture._predictors,
        'edge_order': resistance_network_fixture._edge_order,
        'gendist': resistance_network_fixture._gendist,
        # Use defaults or specific values for the additional attributes
        'fitmetric': 'aic',
        'fixWeight': False,
        'fixShape': False,
        'min_weight': 0.0,
        'max_shape': 100.0
    }

    # Instantiate the ResistanceNetworkWorker
    worker = ResistanceNetworkWorker(**worker_args)

    # Assertions to verify that the worker is correctly instantiated
    assert worker.pop_agg == resistance_network_fixture.pop_agg
    assert worker.reachid_col == resistance_network_fixture.reachid_col
    assert worker.length_col == resistance_network_fixture.length_col
    assert worker.variables == resistance_network_fixture.variables
    assert worker.agg_opts == resistance_network_fixture.agg_opts

    # Validate unique attributes for ResistanceNetworkWorker
    assert worker.fitmetric == 'aic'
    assert worker.fixWeight is False
    assert worker.fixShape is False
    assert worker.min_weight == 0.0
    assert worker.max_shape == 100.0


def _three_variables(units):
    """A network of three variables, for counting parameters only."""
    net = object.__new__(ResistanceNetwork)
    net._predictors = pd.DataFrame(columns=["a", "b", "c"])
    net.num_attributes = 4
    net.units = units
    net.fix_resist = None
    return net


def test_surface_parameters():
    # sel, weight, trans, shape for three variables
    one = [1, 0.5, 1, 3, 0, 0.2, 0, 1, 0, 0.9, 5, 10]
    three = [1, 0.5, 1, 3, 1, 0.2, 0, 1, 1, 0.9, 5, 10]
    none = [0, 0.5, 1, 3, 0, 0.2, 0, 1, 0, 0.9, 5, 10]

    # Per reach: transform + shape per variable, plus m - 1 relative
    # weights
    net = _three_variables("reach")
    assert net.surface_parameters(one) == 2
    assert net.surface_parameters(three) == 8
    assert net.surface_parameters(none) == 0

    net.fixWeight = True
    assert net.surface_parameters(three) == 6
    net.fixShape = True
    assert net.surface_parameters(three) == 0
    net.fixWeight = False
    assert net.surface_parameters(three) == 2

    # Per km: plus the maximum resistance per km if it was optimised (3
    # per variable in all, as ResistanceGA's equation, shape and maximum)
    net = _three_variables("km")
    assert net.surface_parameters(one + [50.0]) == 3
    assert net.surface_parameters(three + [50.0]) == 9
    # K is not read as a selected variable (even if 1), nor counted
    # without variables
    assert net.surface_parameters(none + [1.0]) == 0
    # Not optimised: absent (e.g. the grid search), or fixed
    assert net.surface_parameters(three) == 8
    net.fix_resist = 50.0
    assert net.surface_parameters(three + [50.0]) == 8


def test_evaluate_penalises_aic(resistance_network_fixture):
    fitness, res = resistance_network_fixture.evaluate([1, 1.0, 1, 3])
    loglik, _, aic, delta_aic_null, _, aic_null = res
    # negated AIC: MLPE's k = 4 (2 fixed effects, 2 variance parameters)
    # plus 2 surface parameters
    assert aic == pytest.approx(2 * loglik - 2 * (4 + 2))
    assert delta_aic_null == pytest.approx(aic - aic_null)
    assert fitness == aic


def _random_coords(network_path, n, seed, path):
    # Write random network nodes (of its largest connected part) as a
    # coordinates file
    G = utils.read_pickle_network(network_path)
    rng = np.random.default_rng(seed)
    nodes = sorted(max(nx.connected_components(G), key=len))
    pick = rng.choice(len(nodes), size=n, replace=False)
    pd.DataFrame({
        "sample": [f"s{i}" for i in range(n)],
        "lat": [nodes[i][1] for i in pick],
        "long": [nodes[i][0] for i in pick],
    }).to_csv(path, sep="\t", index=False)


def test_incidence_matches_paths(network_graph_path, tmp_path):
    # Several reaches can share an ID (e.g. the EDGE_ID segments of
    # autoStreamTree outputs), so an ID can be partly on a path; it must be
    # counted if any of its reaches are. Here pairs of reaches share an ID.
    G = utils.read_pickle_network(network_graph_path)
    reaches = sorted({d["HYRIV_ID"] for _, _, d in G.edges(data=True)})
    seg = {r: i // 2 for i, r in enumerate(reaches)}
    nx.set_edge_attributes(G, {(a, b): seg[d["HYRIV_ID"]]
                               for a, b, d in G.edges(data=True)}, "SEG_ID")
    network = tmp_path / "shared.network"
    with open(network, "wb") as fh:
        pickle.dump(G, fh)
    coords = tmp_path / "coords.txt"
    _random_coords(network_graph_path, 12, 1, coords)
    net = ResistanceNetwork(
        network=str(network), shapefile=None, coords=str(coords),
        variables=["LENGTH_KM"], inmat=None,
        agg_opts={"LENGTH_KM": "ARITH"}, reachid_col="SEG_ID",
        length_col="LENGTH_KM", out=str(tmp_path / "out"), verbose=False)
    edge_to_index = {e: i for i, e in enumerate(net._edge_order)}
    nodes = list(net._points_snapped.keys())
    expected = np.zeros_like(net._inc)
    partly = 0
    for row, (i, j) in enumerate(itertools.combinations(range(len(nodes)), 2)):
        path = nx.shortest_path(net._K, nodes[i], nodes[j])
        on = {net._K[a][b]["SEG_ID"] for a, b in zip(path, path[1:])}
        for s in on:
            expected[row, edge_to_index[s]] = 1
        # IDs with a reach on the path and another one off it
        partly += sum(
            1 for a, b, d in net._K.edges(data=True)
            if d["SEG_ID"] in on and not utils.find_pair(path, a, b))
    assert partly > 0
    assert np.array_equal(net._inc, expected)


def test_snap_points_merges_samples(
        network_graph_path, coords_path, tmp_path):
    # Two different coordinates snapping to the same node keep both samples
    pts = pd.read_csv(coords_path, sep="\t")
    first = pts.iloc[0]
    extra = pd.DataFrame([{"sample": "extra", "lat": first.lat + 1e-6,
                           "long": first.long + 1e-6}])
    coords = tmp_path / "coords.txt"
    pd.concat([pts, extra]).to_csv(coords, sep="\t", index=False)
    net = ResistanceNetwork(
        network=network_graph_path, shapefile=None, coords=str(coords),
        variables=["LENGTH_KM"], inmat=None,
        agg_opts={"LENGTH_KM": "ARITH"}, reachid_col="HYRIV_ID",
        length_col="LENGTH_KM", out=str(tmp_path / "out"), verbose=False)
    names = [n for v in net._points_names.values() for n in v]
    assert sorted(names) == sorted(list(pts["sample"]) + ["extra"])


@pytest.mark.parametrize("units", ["km", "reach"])
def test_simulate_matches_fitted_model(
        sim_resistance_network_params, network_graph_path, tmp_path, units):
    # Simulated resistances must equal those of the same model as fitted
    # (per km, with the same maximum resistance per km)
    spec = pd.DataFrame({
        "VAR": ["LENGTH_KM", "UPLAND_SKM"],
        "WEIGHT": [1.0, 0.4],
        "TRANSFORM": [0, 5],
        "SHAPE": [0, 20],
        "AGG": ["SUM", "MAX"],
    })
    spec_file = tmp_path / "spec.tsv"
    spec.to_csv(spec_file, sep="\t", index=False)
    out = str(tmp_path / "sim")
    sim_engine = SimResistanceNetwork(**sim_resistance_network_params,
                                      seed=5, units=units,
                                      max_resistance=37.0)
    sim_engine.simulate(str(spec_file), 1, 12, out)

    sim = pd.read_csv(f"{out}_1.ResistanceMatrix.tsv", sep="\t", index_col=0)
    sim.columns = sim.columns.astype(int)
    assert np.array_equal(sim.to_numpy(), sim.to_numpy().T)

    net = ResistanceNetwork(
        network=network_graph_path, shapefile=None,
        coords=f"{out}_1.coords", variables=list(spec["VAR"]), inmat=None,
        agg_opts=dict(zip(spec["VAR"], spec["AGG"])),
        reachid_col="HYRIV_ID", length_col="LENGTH_KM",
        out=str(tmp_path / "fit"), verbose=False, units=units)
    model = []
    for _, s in spec.set_index("VAR").loc[net._predictors.columns].iterrows():
        model.extend([1, s["WEIGHT"], int(s["TRANSFORM"]), int(s["SHAPE"])])
    if units == "km":
        model.append(37.0)
    fitted, _ = net.model_output(model)
    labels = [int(x) for x in net._points_labels.values()]
    np.testing.assert_allclose(
        sim.loc[labels, labels].to_numpy(), fitted / fitted.max(),
        rtol=0, atol=1e-12)


def test_initialization(sim_resistance_network_params):
    sim_engine = SimResistanceNetwork(**sim_resistance_network_params)
    assert sim_engine.reachid_col == \
        sim_resistance_network_params["reachid_col"]
    assert sim_engine.length_col == \
        sim_resistance_network_params["length_col"]
    assert sim_engine.verbose == \
        sim_resistance_network_params["verbose"]


@pytest.mark.parametrize("num_reps,num_samples", [(1, 10), (2, 5), (3, 3)])
def test_simulate_with_different_parameters(
        sim_resistance_network_params, fake_spec_file, num_reps, num_samples):
    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "simulation_output")

        # Initialize the SimResistanceNetwork with mocked components
        # with patch("resistnet.utils") as mock_utils:
        #     mock_utils.sample_nodes.return_value = {
        #         i: (0.0, 0.0) for i in range(num_samples)
        #     }
        sim_engine = SimResistanceNetwork(**sim_resistance_network_params)

        # Simulate the resistance network
        sim_engine.simulate(
            spec_file=fake_spec_file,
            num_reps=num_reps,
            num_samples=num_samples,
            out=output_prefix
        )

        # Check for the correct number of output files and their contents
        for r in range(1, num_reps + 1):
            resistance_matrix_file = \
                f"{output_prefix}_{r}.ResistanceMatrix.tsv"
            coords_file = f"{output_prefix}_{r}.coords"

            # Check if files exist
            assert os.path.exists(
                resistance_matrix_file), \
                f"{resistance_matrix_file} missing"
            assert os.path.exists(coords_file), f"{coords_file} missing"

            # Check the number of rows in the coords file
            with open(coords_file, 'r') as file:
                num_lines = len(file.readlines())
                expected_lines = num_samples + 1  # +1 for the header
                assert num_lines == expected_lines, \
                    f"{coords_file} should have {expected_lines} \
                        lines but has {num_lines}"
            # Check the dimensions of the resistance matrix
            resistance_matrix_file = \
                f"{output_prefix}_{r}.ResistanceMatrix.tsv"
            df = pd.read_csv(
                resistance_matrix_file, sep='\t', header=0, index_col=0)
            expected_shape = (num_samples, num_samples)
            assert df.shape == expected_shape, \
                f"{resistance_matrix_file} should have shape {expected_shape}"


def test_simulate_samples_connected_part(
        sim_resistance_network_params, network_graph_path, tmp_path):
    # nmt.network has several unconnected parts (other reaches in the
    # layer): samples come from the largest, so every pair has a path
    G = utils.read_pickle_network(network_graph_path)
    parts = list(nx.connected_components(G))
    assert len(parts) > 1
    largest = max(parts, key=len)
    spec = tmp_path / "spec.tsv"
    pd.DataFrame({"VAR": ["LENGTH_KM"], "WEIGHT": [1.0], "TRANSFORM": [0],
                  "SHAPE": [0]}).to_csv(spec, sep="\t", index=False)
    sim = SimResistanceNetwork(**sim_resistance_network_params, seed=3)
    sim.simulate(str(spec), 3, 10, str(tmp_path / "sim"))
    nodes = np.array(list(G.nodes()))
    in_largest = np.array([tuple(n) in largest for n in map(tuple, nodes)])
    for r in range(1, 4):
        pts = pd.read_csv(tmp_path / f"sim_{r}.coords", sep="\t")
        assert len(pts) == 10
        # Each written point is a network node (to the precision written)
        for lat, lon in zip(pts["lat"], pts["long"]):
            d = np.hypot(nodes[:, 0] - lon, nodes[:, 1] - lat)
            assert d.min() < 1e-9
            assert in_largest[np.argmin(d)]


def test_mlpe_reference_inputs_nmt(tmp_path):
    # The genetic distances and the standard-model distances behind the
    # lme4 references (test_MLPE.py) are what ResistNet computes for NMT:
    # river distance in km, and degree of regulation (inverse
    # monomolecular, shape 1)
    import resistnet.resist_dist as rd
    ref = os.path.join(os.path.dirname(__file__), "reference")
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    net = ResistanceNetwork(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"),
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=["dor_pc_pva"],
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={"dor_pc_pva": "ARITH"}, out=str(tmp_path / "net"),
        verbose=False)
    Y = np.loadtxt(os.path.join(ref, "nmt_mlpe_Y.tsv"))
    # (input parsing can differ in the last digit between pandas versions)
    np.testing.assert_allclose(net._gendist, Y, rtol=1e-12, atol=0)
    length = {d[net.reachid_col]: d[net.length_col]
              for _, _, d in net._K.edges(data=True)}
    km = np.array([float(length[e]) for e in net._edge_order])
    np.testing.assert_allclose(
        rd.effectiveResistanceMatrix(net._points_snapped, net._inc, km),
        np.loadtxt(os.path.join(ref, "nmt_mlpe_X_streamkm.tsv")),
        rtol=1e-12)
    np.testing.assert_allclose(
        rd.effectiveResistanceMatrix(
            net._points_snapped, net._inc,
            net._build_composite_surface([1, 1.0, 7, 1])),
        np.loadtxt(os.path.join(ref, "nmt_mlpe_X_dor_inv_mono.tsv")),
        rtol=1e-12)


def test_per_km_surface_matches_resistancega(tmp_path):
    # ResistanceGA's composite surface (R/Combine_Surfaces.R: each
    # continuous surface 1 + max * A * f, summed, divided by its minimum),
    # computed in R with its own functions for the NMT reaches
    # (reference/resistancega_surfaces.R), equals the resistance per km,
    # 1 + (K - 1) * weighted mean of f, divided by its minimum, with
    # weights max_i * A_i and K - 1 = sum(max_i * A_i) / n
    ref_dir = os.path.join(os.path.dirname(__file__), "reference")
    surfaces = pd.read_csv(os.path.join(ref_dir, "resistancega_surfaces.tsv"),
                           sep="\t")
    constants = pd.read_csv(
        os.path.join(ref_dir, "resistancega_surface_constants.tsv"),
        sep="\t")
    code = {"Ricker": 1, "Rev.Ricker": 2, "Inv.Ricker": 3,
            "Inv.Rev.Ricker": 4, "Monomolecular": 5, "Rev.Monomolecular": 6,
            "Inv.Monomolecular": 7, "Inv.Rev.Monomolecular": 8}
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    variables = ["dor_pc_pva", "sgr_dk_rav", "ele_mt_cav"]
    net = ResistanceNetwork(
        network=os.path.join(data, "nmt.network"), shapefile=None,
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=variables, inmat=None,
        agg_opts={v: "ARITH" for v in variables}, reachid_col="HYRIV_ID",
        out=str(tmp_path / "nmt"), verbose=False, units="km")
    for s, c in constants.groupby("set"):
        c = c.set_index("variable")
        ind = []
        for v in variables:
            if v in c.index:
                ind += [1, c.loc[v, "max"] * c.loc[v, "A"],
                        code[c.loc[v, "transformation"]], c.loc[v, "shape"]]
            else:
                ind += [0, 1.0, 0, 1]
        ind.append(1.0 + (c["max"] * c["A"]).sum() / len(c))
        rho = net.resistance_surface(ind)
        ref = surfaces[surfaces["set"] == s].set_index("HYRIV_ID")["value"]
        ref = ref.loc[list(net._edge_order)].to_numpy()
        np.testing.assert_allclose((rho / rho.min()).to_numpy(), ref,
                                   rtol=1e-12, err_msg=f"set {s}")
        # Each reach's resistance is its resistance per km times its
        # length (ResistanceGA's geoCorrection: per distance between cells)
        np.testing.assert_allclose(net.edge_resistances(ind).to_numpy(),
                                   rho.to_numpy() * net._lengths, rtol=1e-15)


def test_per_km_contributions(tmp_path):
    # Shares of the resistance per km due to distance (the baseline of 1)
    # and to each variable, averaged over the network's km: they sum to 1,
    # and distance is everything when K = 1
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    variables = ["dor_pc_pva", "sgr_dk_rav"]
    net = ResistanceNetwork(
        network=os.path.join(data, "nmt.network"), shapefile=None,
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=variables, inmat=None,
        agg_opts={v: "ARITH" for v in variables}, reachid_col="HYRIV_ID",
        out=str(tmp_path / "nmt"), verbose=False, units="km")
    ind = [1, 0.6, 7, 1, 1, 0.2, 5, 3, 40.0]
    shares = net.contributions(ind)
    assert set(shares) == {"distance", *variables}
    assert sum(shares.values()) == pytest.approx(1.0, abs=1e-12)
    rho = net.resistance_surface(ind).to_numpy()
    km = net._lengths
    assert shares["distance"] == pytest.approx(
        np.sum(km / rho) / km.sum(), rel=1e-12)
    f = net.transform(net._predictors["dor_pc_pva"], 7, 1).to_numpy()
    assert shares["dor_pc_pva"] == pytest.approx(
        np.sum(km * 39.0 * 0.75 * f / rho) / km.sum(), rel=1e-12)
    assert net.contributions(ind[:-1] + [1.0])["distance"] == \
        pytest.approx(1.0, abs=1e-12)
    # Per reach: not defined
    net.units = "reach"
    assert net.contributions(ind) is None


def test_paths_follow_river_length(tmp_path):
    # With a loop (two routes between the sites: 1 km in two reaches, or
    # one 10 km reach), the subgraph, incidence matrix and river distance
    # follow the shorter route (autoStreamTree's paths by river length)
    G = nx.Graph()
    a, x, b = (0.0, 0.0), (0.5, 0.5), (1.0, 0.0)
    G.add_edge(a, x, HYRIV_ID=1, LENGTH_KM=0.5, var=1.0)
    G.add_edge(x, b, HYRIV_ID=2, LENGTH_KM=0.5, var=2.0)
    G.add_edge(a, b, HYRIV_ID=3, LENGTH_KM=10.0, var=3.0)
    with open(tmp_path / "loop.network", "wb") as fh:
        pickle.dump(G, fh)
    pd.DataFrame({"sample": ["s1", "s2"], "lat": [0.0, 0.0],
                  "long": [0.0, 1.0]}).to_csv(tmp_path / "coords.txt",
                                              sep="\t", index=False)
    net = ResistanceNetwork(
        network=str(tmp_path / "loop.network"), shapefile=None,
        coords=str(tmp_path / "coords.txt"), variables=["var"], inmat=None,
        agg_opts={"var": "ARITH"}, reachid_col="HYRIV_ID",
        length_col="LENGTH_KM", out=str(tmp_path / "out"), verbose=False)
    assert sorted(net._edge_order) == [1, 2]
    assert net._inc.tolist() == [[1, 1]]
    from resistnet.report import river_distance_matrix
    assert river_distance_matrix(net)[0, 1] == pytest.approx(1.0)
