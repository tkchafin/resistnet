"""
The coalescent model on a river network: Lundgren & Ralph's (2019)
coalescence times on a lattice of the river
(resistnet.coalescent.RiverCoalescent), with movement from the
standard model's resistance surface.

Lineages move back in time between neighbouring points of the lattice
(points at confluences, tips, the outlet and the sites, and about cell_km
apart in between). The rate of a lineage at point x jumping to its
neighbour y is the link's conductance, 1 / (its resistance: resistance per
km in series over the river it covers; r: the standard model's surface,
1 + (K - 1) x the weighted mean of the transformed variables), times
exp(-rate x km moved upstream) with a flow term, divided by the length of
river around x (half its links): a finite volume version of a diffusion,
so that results do not depend on the lattice's spacing once it is fine.
Two lineages at the same point coalesce at rate 1 / (density x length of
river around it) (the population living there, density per km). The mean
coalescence times C solve Lundgren & Ralph's equation 7, and the genetic
distance between two sites is linearised FST from them (Slatkin 1991),
(C_ij - Cw) / Cw with Cw the mean of the two within-site times, fitted
with MLPE as the other models.

With a flow term, lineages traced back in time move upstream more easily
when the rate is negative, i.e. gene flow is biased downstream.
"""
import numpy as np

import resistnet.MLPE as mlpe
import resistnet.utils as utils
from resistnet.coalescent import RiverCoalescent
from resistnet.samc_network import ResistanceNetworkSAMC
from resistnet.samc_network import ResistanceNetworkSAMCWorker


class CoalescentMethods:
    """The coalescent model's computations, shared by the network and its
    worker (see ResistanceNetworkCoalescent)."""

    null_prefix = "coal"

    def _fit(self, reach_r, rate):
        """MLPE results for the linearised FST predicted with the given
        resistance per km of each reach (None: 1, distance alone) and flow
        rate (None if not finite)."""
        lin = self._lattice.linfst(reach_r, rate)
        if not np.all(np.isfinite(lin)):
            return lin, None
        return lin, mlpe.MLPE(self._gendist, lin, scale=True)

    def evaluate(self, individual):
        """
        Fitness of a model (as ResistanceNetworkSAMC.evaluate()): the
        linearised FST predicted from the coalescence times among sites,
        fitted to the genetic distances with MLPE.
        """
        movement = self._surface_and_flow(individual)
        if movement is None:
            return (float('-inf'), [np.nan, np.nan, np.nan, np.nan])
        _, res = self._fit(*movement)
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
        return self._fit(None, rate)[1]

    def model_output(self, model):
        """
        The predicted linearised FST among sites, and for each reach its
        resistance per km and the probabilities of a lineage moving
        downstream and upstream inside it (see
        resistnet.coalescent.SegmentLattice.reach_moves()).
        """
        movement = self._surface_and_flow(model)
        if movement is None:
            return None, None
        rho, rate = movement
        lin, _ = self._fit(rho, rate)
        if rho is None:
            rho = np.ones(len(self._edge_order))
        moves = self._lattice.lattice.reach_moves(rate)
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
        The coalescent model of the network (see
        resistnet.coalescent.RiverCoalescent: its lattice, with the
        sites as points, and the density; no absorption: the coalescence
        rates are the model's).
        """
        sites = list(self._points_snapped.keys())
        self._lattice = RiverCoalescent(
            self._K, list(self._edge_order), sites, self.reachid_col,
            self.length_col, cell_km=self.cell_km, density=self.density)
        self._site_states = self._lattice.sites
        self._R_lat = np.zeros(self._lattice.lattice.n)

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
    """Worker for the coalescent model (as ResistanceNetworkSAMCWorker; its
    lattice, a RiverCoalescent, holds the density)."""
