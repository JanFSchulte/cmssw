"""Hadronic-recoil-vs-Z-pT helpers: dimuon Z reconstruction (from the new
scoutingPuppiCalibMuonTable) and the standard u_parallel/u_perp recoil
decomposition, u_vec = -(MET_vec + Z_pT_vec) -- at truth level (MET=0, no
real invisible particles in DY) this exactly balances the Z: u_vec =
-Z_pT_vec, so u_parallel (projected onto -Z_hat) responds around +Z_pT and
u_perp scatters around 0; both are standard MET/recoil calibration
diagnostics (CMS Z+jets recoil corrections use exactly this decomposition).
"""

import awkward as ak
import numpy as np

Z_MASS = 91.1876


def select_z_candidates(pt, eta, phi, mass, charge, mass_window=20.0):
    """Best (mass-closest-to-mZ) opposite-charge dimuon pair per event.

    Returns (has_z, z_pt, z_phi, z_mass) as flat-per-event numpy arrays
    (has_z=False events get NaN for the others).
    """
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e = np.sqrt(px * px + py * py + pz * pz + mass * mass)

    pairs = ak.combinations(ak.zip({"px": px, "py": py, "pz": pz, "e": e, "charge": charge}), 2,
                             fields=["a", "b"])
    os_pairs = pairs[pairs.a.charge != pairs.b.charge]

    sum_px = os_pairs.a.px + os_pairs.b.px
    sum_py = os_pairs.a.py + os_pairs.b.py
    sum_pz = os_pairs.a.pz + os_pairs.b.pz
    sum_e = os_pairs.a.e + os_pairs.b.e
    m2 = sum_e * sum_e - sum_px * sum_px - sum_py * sum_py - sum_pz * sum_pz
    pair_mass = np.sqrt(ak.where(m2 > 0, m2, 0.0))
    pair_pt = np.sqrt(sum_px * sum_px + sum_py * sum_py)
    pair_phi = np.arctan2(sum_py, sum_px)

    dmass = np.abs(pair_mass - Z_MASS)
    in_window = dmass < mass_window
    dmass_masked = ak.where(in_window, dmass, np.inf)

    has_z = ak.any(in_window, axis=1)
    best = ak.argmin(dmass_masked, axis=1, keepdims=True)

    z_pt = ak.firsts(pair_pt[best])
    z_phi = ak.firsts(pair_phi[best])
    z_mass = ak.firsts(pair_mass[best])

    has_z_np = np.asarray(has_z)
    z_pt_np = np.where(has_z_np, np.asarray(ak.fill_none(z_pt, np.nan)), np.nan)
    z_phi_np = np.where(has_z_np, np.asarray(ak.fill_none(z_phi, np.nan)), np.nan)
    z_mass_np = np.where(has_z_np, np.asarray(ak.fill_none(z_mass, np.nan)), np.nan)
    return has_z_np, z_pt_np, z_phi_np, z_mass_np


def puppi_met_from_candidates(cand_pt, cand_phi, weight):
    """METx/METy (negative vector sum of weighted candidates), per event."""
    wpx = cand_pt * weight * np.cos(cand_phi)
    wpy = cand_pt * weight * np.sin(cand_phi)
    metx = -np.asarray(ak.sum(wpx, axis=1))
    mety = -np.asarray(ak.sum(wpy, axis=1))
    return metx, mety


def recoil_components(met_x, met_y, z_pt, z_phi):
    """u_parallel (along -Z direction, so it responds around +z_pt) and
    u_perp (perpendicular), from u_vec = -(MET_vec + Z_pT_vec).
    """
    z_x = z_pt * np.cos(z_phi)
    z_y = z_pt * np.sin(z_phi)
    u_x = -(met_x + z_x)
    u_y = -(met_y + z_y)
    # project onto -z_hat = (-cos(z_phi), -sin(z_phi))
    u_par = -(u_x * np.cos(z_phi) + u_y * np.sin(z_phi))
    u_perp = -u_x * np.sin(z_phi) + u_y * np.cos(z_phi)
    return u_par, u_perp
