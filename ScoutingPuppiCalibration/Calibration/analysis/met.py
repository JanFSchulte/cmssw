"""MET performance: bias/resolution of the x/y components of PF MET vs
GenMET, binned in pileup -- the standard, boson-independent MET performance
metric (unlike hadronic recoil, which needs a real leptonically-decaying
boson to define a recoil axis against -- not available yet, see
results/met_analysis notes).

response is measured as PF MET's cartesian components minus GenMET's, not a
pt ratio: GenMET in QCD (and even DY, absent real neutrinos beyond the
occasional semileptonic decay) runs close to zero, so a ratio metric's
denominator blows up in exactly the region with the most events -- the
standard CMS "u1/u2" or dx/dy component decomposition avoids that by
construction.
"""

import numpy as np


def met_xy(pt, phi):
    return pt * np.cos(phi), pt * np.sin(phi)


def component_resolution_vs_pu(dx, dy, nTrueInt, pu_bins=None):
    """Median + IQR/2 resolution of dx=METx-GenMETx and dy=METy-GenMETy,
    binned in true pileup. Returns (centers, bias_x, res_x, bias_y, res_y, counts).
    """
    if pu_bins is None:
        pu_bins = np.arange(0, 80, 5)
    idx = np.digitize(nTrueInt, pu_bins) - 1
    centers, bias_x, res_x, bias_y, res_y, counts = [], [], [], [], [], []
    for b in range(len(pu_bins) - 1):
        sel = idx == b
        n = np.count_nonzero(sel)
        if n < 20:
            continue
        qx16, qx50, qx84 = np.percentile(dx[sel], [16, 50, 84])
        qy16, qy50, qy84 = np.percentile(dy[sel], [16, 50, 84])
        centers.append(0.5 * (pu_bins[b] + pu_bins[b + 1]))
        bias_x.append(qx50)
        res_x.append(0.5 * (qx84 - qx16))
        bias_y.append(qy50)
        res_y.append(0.5 * (qy84 - qy16))
        counts.append(n)
    return (np.array(centers), np.array(bias_x), np.array(res_x),
            np.array(bias_y), np.array(res_y), np.array(counts))
