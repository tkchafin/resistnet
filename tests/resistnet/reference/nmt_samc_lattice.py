"""Exports the lattice SAMC models of test_samc_network.R_SAMC_MODELS for
the NMT example (their transition matrices with an absorbing state, as
samc(data = P) takes them, and the site states), for
nmt_samc_cond_passage.R; and the symmetrised passage times of the
sgr_flow model, the SAMC model distance of the lme4 references
(nmt_mlpe_X_samc_sgr_flow.tsv). Run in this directory:
python nmt_samc_lattice.py <scratch dir>"""
import os
import sys
import warnings

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                ".."))
from test_samc_network import R_SAMC_MODELS, _ref_net, _ref_transition_matrix  # noqa: E402,E501
from resistnet.CFPT import CFPT  # noqa: E402

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
net = _ref_net(out)
R = net._R_lat
for name in R_SAMC_MODELS:
    Q = _ref_transition_matrix(net, name)
    n = Q.shape[0]
    P = sp.bmat([[Q, sp.csr_matrix(R.reshape(-1, 1))],
                 [sp.csr_matrix((1, n)), sp.csr_matrix([[1.0]])]]).tocsc()
    sio.mmwrite(os.path.join(out, f"{name}.P.mtx"), P, precision=17)
    np.savetxt(os.path.join(out, f"{name}.dests.txt"),
               np.array(net._site_states) + 1, fmt="%d")
    if name == "sgr_flow":
        T = CFPT(Q, net._site_states, R)
        np.savetxt("nmt_mlpe_X_samc_sgr_flow.tsv", (T + T.T) / 2,
                   delimiter="\t", fmt="%.17g")
print(f"{len(R_SAMC_MODELS)} models, {net._lattice.n} states")
