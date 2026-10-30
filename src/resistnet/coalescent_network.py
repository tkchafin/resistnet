"""
The coalescent model on a river network: Lundgren & Ralph's (2019)
coalescence times (resistnet.coalescent) on a lattice of the river
(resistnet.lattice.SegmentLattice: points at confluences, tips, the outlet
and the sites, and about cell_km apart in between), with movement from the
standard model's resistance surface.

Lineages move back in time between neighbouring points of the lattice.
The rate of a lineage at point x jumping to its neighbour y is the link's
conductance, 1 / (its resistance: resistance per km in series over the
river it covers; r: the standard model's surface, 1 + (K - 1) x the
weighted mean of the transformed variables), times exp(-rate x km moved
upstream) with a flow term, divided by the length of river around x (half
its links): a finite volume version of a diffusion, so that results do not
depend on the lattice's spacing once it is fine. Two lineages at the same
point coalesce at rate 1 / (density x length of river around it) (the
population living there, density per km). The mean coalescence times C
solve Lundgren & Ralph's equation 7, and the genetic distance between two
sites is
linearised FST from them (Slatkin 1991), (C_ij - Cw) / Cw with Cw the mean
of the two within-site times, fitted with MLPE as the other models.

With a flow term, lineages traced back in time move upstream more easily
when the rate is negative, i.e. gene flow is biased downstream.
"""
import numpy as np

import resistnet.MLPE as mlpe
import resistnet.utils as utils
import resistnet.coalescent as co
from resistnet.lattice import SegmentLattice
from resistnet.samc_network import ResistanceNetworkSAMC
from resistnet.samc_network import ResistanceNetworkSAMCWorker


class CoalescentMethods:
    """The coalescent model's computations, shared by the network and its
    worker (see ResistanceNetworkCoalescent)."""

    null_prefix = "coal"

    def _river_around(self):
        """Length of river around each lattice point (half its links)."""
        lat = self._lattice
        return np.bincount(lat.a, weights=lat.length / 2, minlength=lat.n)

    def _jump_rates(self, individual):
        """Lineage jump rates on each lattice link (lattice.a -> b), or None
        if no variable is selected and there is no flow term."""
        W = self._movement_weights(individual)
        if W is None:
            return None
        lat = self._lattice
        w = np.asarray(W.tocsr()[lat.a, lat.b]).ravel()
        return w / self._river_around()[lat.a]

    def _distance_rates(self, rate=0.0):
        """Jump rates of distance alone (resistance 1 per km), with a flow
        term of the given rate."""
        lat = self._lattice
        W = lat.weights(rate=rate)
        w = np.asarray(W.tocsr()[lat.a, lat.b]).ravel()
        return w / self._river_around()[lat.a]

    def coalescence_times(self, rates):
        """Mean coalescence times among all lattice points."""
        lat = self._lattice
        G = co.generator(rates, list(zip(lat.a, lat.b)), lat.n)
        gamma = 1.0 / (self.density * self._river_around())
        return co.coalescence_times(G, gamma)

    def linfst(self, rates):
        """Linearised FST among the sites from the coalescence times."""
        C = self.coalescence_times(rates)
        s = self._site_states
        Cs = C[np.ix_(s, s)]
        w = np.diag(Cs)
        within = (w[:, None] + w[None, :]) / 2
        lin = (Cs - within) / within
        np.fill_diagonal(lin, 0.0)
        return lin

    def _fit(self, rates):
        """MLPE results for a matrix of linearised FST (None if not
        finite)."""
        lin = self.linfst(rates)
        if not np.all(np.isfinite(lin)):
            return lin, None
        return lin, mlpe.MLPE(self._gendist, lin, scale=True)

    def evaluate(self, individual):
        """
        Fitness of a model (as ResistanceNetworkSAMC.evaluate()): the
        linearised FST predicted from the coalescence times among sites,
        fitted to the genetic distances with MLPE.
        """
        rates = self._jump_rates(individual)
        if rates is None:
            return (float('-inf'), [np.nan, np.nan, np.nan, np.nan])
        _, res = self._fit(rates)
        if res is None:
            return (float('-inf'), [np.nan, np.nan, np.nan, np.nan])
        k = self.surface_parameters(individual)
        res["aic"] -= 2 * k
        res["delta_aic_null"] -= 2 * k
        fitness = res[self.fitmetric].iloc[0]
        return fitness, list(res.iloc[0])

    def _flow_alone(self, rate):
        """MLPE results for distance alone (resistance 1 per km) with a
        flow term of the given rate."""
        return self._fit(self._distance_rates(rate))[1]

    def model_output(self, model):
        """
        The predicted linearised FST among sites, and for each reach its
        resistance per km and the probabilities of a lineage moving
        downstream and upstream inside it (see Lattice.reach_moves()).
        """
        rates = self._jump_rates(model)
        if rates is None:
            return None, None
        lin, _ = self._fit(rates)
        rho = self.resistance_surface(model)
        rho = (np.ones(len(self._edge_order)) if rho is None
               else rho.to_numpy(dtype=float))
        moves = self._lattice.reach_moves(self.flow_rate(model))
        return lin, np.column_stack((rho, moves))

    def pair_distance(self, M):
        """The predicted genetic distances are already symmetric."""
        return np.asarray(M, dtype=float)


class ResistanceNetworkCoalescent(CoalescentMethods, ResistanceNetworkSAMC):
    """
    The coalescent model on a river network (see the module's description).

    Args:
        As ResistanceNetworkSAMC (absorption and pairs are not used), plus
        density (float): Population per km of river, which sets the
                         coalescence rates relative to the movement rates.
                         Defaults to 1000; on the example data, fits do not
                         change from 100 up (high gene flow, linearised FST
                         of 0.4 and less on average), and are poorer at 10.
        cell_km defaults to 5 (the coalescence times need a system the
        size of the number of pairs of lattice points).
    """

    def __init__(self, *args, density=1000.0, cell_km=5.0, **kwargs):
        if density <= 0:
            raise ValueError("density must be positive")
        self.density = float(density)
        super().__init__(*args, cell_km=cell_km, **kwargs)

    def build_lattice(self):
        """
        The segment lattice of the network (see
        resistnet.lattice.SegmentLattice), with the sites as points (no
        absorption: the coalescence rates are the model's).
        """
        sites = list(self._points_snapped.keys())
        self._lattice = SegmentLattice(self._K, list(self._edge_order),
                                       self.reachid_col, self.length_col,
                                       self.cell_km, keep_nodes=sites)
        self._site_states = self._lattice.states(sites)
        self._R_lat = np.zeros(self._lattice.n)

    def model_spec(self, individual):
        """As ResistanceNetworkSAMC.model_spec(), for the coalescent model
        (type "coalescent", with its density)."""
        spec = super().model_spec(individual)
        spec["type"] = "coalescent"
        spec["density"] = self.density
        for key in ("absorption", "pairs"):
            spec.pop(key, None)
        return spec

    def output_and_plot_model(self, oname, mat_r, edge_r, plot=True):
        """
        Writes the predicted linearised FST among sites
        ({oname}.CoalescenceMatrix.tsv) and the reaches' resistances and
        lineage moves ({oname}.TransitionEdges.tsv), and plots them.
        """
        if mat_r is None or edge_r is None:
            return
        import pandas as pd
        edge_r = np.asarray(edge_r)
        edf = pd.DataFrame({
            "EDGE_ID": self._edge_order,
            "Resistance": edge_r[:, 0],
            "Downstream": edge_r[:, 1],
            "Upstream": edge_r[:, 2]
        })
        edf.to_csv(f"{oname}.TransitionEdges.tsv", sep="\t", index=False,
                   na_rep="NA")
        utils.write_matrix(f"{oname}.CoalescenceMatrix.tsv", mat_r,
                           list(self._points_labels.values()))
        if plot:
            self.plot_resistance_network(edf, oname)
            if self._gendist is not None:
                self.plot_pairwise_model(mat_r, oname, partition=False)


class ResistanceNetworkCoalescentWorker(CoalescentMethods,
                                        ResistanceNetworkSAMCWorker):
    """Worker for the coalescent model (as ResistanceNetworkSAMCWorker,
    plus the density)."""

    def __init__(self, *args, density=1000.0, **kwargs):
        self.density = float(density)
        super().__init__(*args, **kwargs)
