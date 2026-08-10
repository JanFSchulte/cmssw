"""Standalone JEC application for the offline FastJet-reclustered jets in
puppi_refit.py, using the same AK4PFHLT payload (L1FastJet, L2Relative,
L3Absolute, L2L3Residual) that production's patJetCorrFactors applies (see
scoutingPFJetReclusterCorrFactors in scoutingToMiniAODDerivedCollections_cff.py).

The payload text files (data/jec/150X_mcRun3_2024_realistic_v2_*_AK4PFHLT.txt)
were extracted from the GlobalTag with CondTools/JetMET's JetCorrectorDBReader
-- see test/dump_jec_payload_cfg.py (run once with cmsRun; the GlobalTag/
conddb machinery isn't available outside CMSSW, so these must be regenerated
if the GlobalTag ever changes).

Needs PyROOT + CondFormats/JetMETObjects (libCondFormatsJetMETObjects.so from
the CMSSW release), not a pip package -- run inside a `cmsenv`.
"""

import os

import ROOT

_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "jec")
_GT = "150X_mcRun3_2024_realistic_v2"
_LEVELS = ["L1FastJet", "L2Relative", "L3Absolute", "L2L3Residual"]

_corrector = None


def get_corrector():
    """Lazily build and cache the FactorizedJetCorrector (AK4PFHLT, same
    levels/order as production's patJetCorrFactors)."""
    global _corrector
    if _corrector is not None:
        return _corrector
    ROOT.gSystem.Load("libCondFormatsJetMETObjects")
    vpar = ROOT.std.vector("JetCorrectorParameters")()
    for level in _LEVELS:
        path = os.path.join(_DATA_DIR, "%s_%s_AK4PFHLT.txt" % (_GT, level))
        if not os.path.exists(path):
            raise FileNotFoundError(
                "JEC payload not found: %s -- regenerate with "
                "`cmsRun test/dump_jec_payload_cfg.py` (see module docstring)" % path)
        vpar.push_back(ROOT.JetCorrectorParameters(path))
    _corrector = ROOT.FactorizedJetCorrector(vpar)
    return _corrector


def correction_factor(pt, eta, phi, area, rho):
    """Multiplicative JEC factor for one jet (all levels combined).
    Corrected pt = pt * correction_factor(...)."""
    corr = get_corrector()
    corr.setJetPt(float(pt))
    corr.setJetEta(float(eta))
    corr.setJetPhi(float(phi))
    corr.setJetA(float(area))
    corr.setRho(float(rho))
    return corr.getCorrection()
