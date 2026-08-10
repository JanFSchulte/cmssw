"""CRAB submission for the scouting-PUPPI calibration study.

Patterned on ScoutingNanoProduction/submit_scoutingNano.py's MC CRAB block,
but: (a) only a representative subset of the QCD HT bins (low/mid/high,
enough to span the pt/pileup regimes needed for response/resolution tuning,
not the full bin list), and (b) capped statistics per bin (Data.totalUnits)
since this is a calibration study, not a full-dataset production -- O(10-50k)
events per bin is enough.

Dataset names copied locally from ScoutingNanoProduction/submit_scoutingNano.py
(not imported cross-package, to keep this package genuinely standalone).

Usage:
    python submit_scoutingPuppiCalib.py --submit
    python submit_scoutingPuppiCalib.py --resubmit
    python submit_scoutingPuppiCalib.py --report
"""

import argparse
import glob
import os

from CRABClient.UserUtilities import config
from CRABAPI.RawCommand import crabCommand

parser = argparse.ArgumentParser()
parser.add_argument('--submit', '-s', action='store_true', help='submit new CRAB tasks')
parser.add_argument('--resubmit', '-r', action='store_true', help='resubmit failed jobs in existing tasks')
parser.add_argument('--report', action='store_true', help='report events processed for existing tasks')
args = parser.parse_args()

# Representative low/mid/high-HT subset of the QCD_HT-binned MINIAODSIM
# samples (full bin list in ScoutingNanoProduction/submit_scoutingNano.py) --
# spans the pt/pileup regimes response/resolution tuning needs without
# processing every bin.
mc_samples = {
    'QCD_HT-100to200': '/QCD-4Jets_Bin-HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM',
    'QCD_HT-400to600': '/QCD-4Jets_Bin-HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM',
    'QCD_HT-1000to1200': '/QCD-4Jets_Bin-HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM',
    'QCD_HT-2000': '/QCD-4Jets_Bin-HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM',
    # Z(mumu)+jets -- QCD has no real leptonically-decaying boson, so it
    # can't support the standard hadronic-recoil-vs-boson-pT calibration
    # technique (recoil = -(MET + Z pT vector), binned in true Z pT) no
    # matter what columns get added to the QCD production; this sample is
    # what actually makes that analysis possible. Run through the same
    # calibration nano config as the QCD bins (same jet/candidate/MET
    # columns), so recoil can be studied against the same PUPPI variants.
    'DY_MLL-50': '/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v6/MINIAODSIM',
    # WWW (triple-W, 4-flavor) -- genuinely boosted hadronic W decays (unlike
    # QCD/DY, where an AK8 jet is mostly wide-angle QCD radiation, not a
    # single boosted resonance), needed to test the impact of the PUPPI
    # recalibration on softdrop mass (_addSoftDropTable, scoutingPuppiCalibration_cff.py)
    # against a real W->qq mass peak rather than a falling QCD spectrum.
    'WWW_4F': '/WWW-4F_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM',
}

# Calibration-scale statistics per bin: FileBased splitting with a small
# unitsPerJob, capped by totalUnits (~O(10-50k) events/bin, not the full
# dataset). Tune totalUnits after checking events/file for these samples.
UNITS_PER_JOB = 2
TOTAL_UNITS = 40  # -> up to UNITS_PER_JOB*TOTAL_UNITS files processed per sample

# v1 (crab_projects/crab_scoutingPuppiCalib_*, no trailing _v2) predates the
# vtxAssocOff variant (useVertexAssociation=False + UseFromPVLooseTight=True,
# see variants.py). v2 predates vtxAssocImproved (useImprovedVertexAssociation
# in Run3ScoutingParticleToPackedCandidateProducer.cc -- redoes the LV/PU
# split from a real per-vertex dz scan instead of trusting particle.vertex(),
# which is unusable for charged candidates -- see the long comment on
# vtxAssocImproved in variants.py). v3 predates vtxAssocTrackMatched
# (useTrackMatchedVertexAssociation -- refines vtxAssocImproved's linear-
# z-shift dz approximation into a genuine track.dz(Point) using the same
# kinematically-matched reco::Track already used to embed track details) and
# the PuppiRawAlphaToValueMapProducer diagnostic columns (puppiRawAlpha/
# puppiAlphaMed/puppiAlphaRms on the reference variant's candidate table, for
# offline MedEtaSF/RMSEtaSF/MinNeutralPt/MinNeutralPtSlope recalibration --
# see variants.py and _addPuppiAlphaDiagnostics in
# scoutingPuppiCalibration_cff.py). v4's diagnostic columns turned out to be
# corrupted by a real bug in PuppiRawAlphaToValueMapProducer (PuppiAlphasMed/
# PuppiAlphasRms are per-candidate arrays, not [algoBlock][candidate]-flattened
# like PuppiRawAlphas genuinely is -- the plugin was broadcasting one arbitrary
# candidate's value across the whole event). v5 fixes that, adds the new rho
# table (fixedGridRhoFastjetAll, needed for JEC in the offline FastJet
# recalibration pipeline -- see analysis/puppi_refit.py + analysis/jec.py),
# and is the first version whose output can actually support a real offline
# PUPPI-parameter recalibration scan (the whole pipeline -- weight formula,
# reclustering, JEC, gen-matching -- was validated to reproduce production's
# actual nominal PUPPI jets essentially exactly, see the scouting-vertex-
# association-bug memory for the full derivation). Bump to v5 so this
# resubmission gets fresh request names/output areas instead of colliding
# with the completed v1/v2/v3/v4 tasks, whose results are kept for reference.
# v6 adds: (a) a fix to PuppiRawAlphaToValueMapProducer's raw-alpha block
# selection -- it used to hardcode the central block for every candidate,
# silently wrong for the ~14% of candidates in the forward (|eta|>=2.5)
# region, needed to calibrate MedEtaSF/RMSEtaSF there too, not just central;
# (b) AK8 validation -- plain AK8 (customizeForScoutingAK8ReclusteredJets,
# already in PhysicsTools/PatFromScouting but never wired into this package
# before) plus a new PUPPI-weighted AK8 collection per variant, reusing the
# same per-candidate PUPPI weights AK4 already computes; (c) PF MET + GenMET
# NanoAOD tables (previously entirely absent from this package's output),
# for MET/hadronic-recoil validation; (d) the new DY_MLL-50 sample above,
# since hadronic recoil needs a real leptonically-decaying boson QCD can't
# provide. See _addAK8Jets/_addMETTables in scoutingPuppiCalibration_cff.py
# and the etaBoundaries comment in PuppiRawAlphaToValueMapProducer.cc.
# v7 adds: (a) genJetAK8Idx/GenJetAK8Table (fatJetMCTable, gated at pt>100
# to match GenJetAK8Table's own cut) replacing an AK8 MC-table bug where
# v6 silently reused AK4's jetMCTable/pt>10 convention with no GenJetAK8
# table behind it -- v6's AK8 gen-matching indices are unusable, this is
# what actually makes AK8 response/resolution/efficiency analysis possible;
# (b) the "optimized" variant (central_MedEtaSF=0.5, central_MinNeutralPtSlope
# =0.01, the recalibrated operating point validated at full v5/v6 statistics
# this session) is now the DEFAULT referenceVariant, i.e. what
# scoutingPuppiCalib_mc_cfg.py builds unless overridden -- "nominal" (true
# PF/stock defaults) and "vtxAssocTrackMatched" are kept alongside it for
# continued baseline comparison, see variants.py; (c) a minimal muon table
# (_addMuonTable, off slimmedMuons) for dimuon Z reconstruction, needed for
# the hadronic-recoil-vs-boson-pT technique the DY sample was added for in
# v6 but couldn't yet use; (d) soft-drop-groomed mass companion collections
# per AK8 flavor (_addSoftDropTable, standard zcut=0.1/beta=0.0/R0=0.8 CMS
# Run3 working point, own pt/eta/phi/mass table, match to the ungroomed AK8
# jets by deltaR offline -- see its own comment for why not index-embedded
# at production time); (e) the new WWW_4F sample above, for a real boosted
# hadronic-W mass peak to validate softdrop mass against (QCD/DY's AK8 jets
# are mostly wide-angle radiation, not a single boosted resonance).
# v8 ports PUPPI v15's (CMS DP-2021/001) charged-particle handling into the
# scouting setup. Root cause (traced this session, cross-checking DP-2021/001
# against CommonTools/PileupAlgos/plugins/PuppiProducer.cc in this release):
# the release's stock `puppi` config (Puppi_cff.py) already IS v15 -- neutral
# high-pT weight floor (PtMaxNeutralsStartSlope/PtMaxNeutrals), FromPV2Recovery,
# etc. are already its defaults -- but every scouting variant here sets
# useVertexAssociation=True (required, see vtxAssocOff's comment in
# variants.py for why the alternative is broken for scouting), which routes
# PuppiProducer.cc through a *different*, much simpler charged-particle-
# categorization branch (~L200-224 pre-patch) that implements NONE of v15's
# charged-particle protections (PtMaxCharged, UseFromPV2Recovery/
# PtMinForFromPV2Recovery) -- those lived only in the offline dz-based branch,
# structurally unreachable when useVertexAssociation=True. Measured impact on
# the v7 DY sample: 65% of ALL charged scouting candidates (still 13-16%
# above 20-50 GeV) have pvAssocQuality=NotReconstructedPrimary (no vertex
# association at all, mostly a scouting HLT-tracking-coarseness artifact, not
# genuine pileup) and were falling straight through to id=0 (neutral),
# getting PUPPI-suppressed (mean weight 0.43-0.56) with zero recovery.
# PuppiProducer.cc's fUseVertexAssociation branch was patched this session to
# apply the same two v15 protections there, reusing the existing PtMaxCharged/
# UseFromPV2Recovery/PtMinForFromPV2Recovery knobs (no new C++ parameters):
# unconditionally keep high-pT (PtMaxCharged) charged candidates as LV
# regardless of PU/unassociated categorization, and apply a softer pT floor
# (UseFromPV2Recovery, already True/4.GeV by default on _stockPuppi) to
# unassociated candidates specifically. UseFromPV2Recovery's floor is now
# live for EVERY variant below once rebuilt (it was always inherited from
# Puppi_cff.py, just previously dead code on this branch); the new
# "chargedV15" variant additionally sets ptMaxCharged=20 on top of
# "optimized"'s MedEtaSF/slope, to isolate that knob's incremental effect.
# See scoutingPuppiCalibration_cff.py's puppiClone PtMaxCharged comment and
# variants.py's chargedV15 comment for the full derivation.
VERSION = 'v8'

if not args.resubmit and not args.report:
    for sample, dataset in mc_samples.items():
        task_name = f'scoutingPuppiCalib_{sample}_{VERSION}'
        print(f'Configuring task: {task_name}  ({dataset})')

        cfg = config()
        cfg.General.requestName = task_name
        cfg.General.workArea = 'crab_projects'
        cfg.General.transferLogs = True

        cfg.JobType.pluginName = 'Analysis'
        cfg.JobType.psetName = '../test/scoutingPuppiCalib_mc_cfg.py'
        cfg.JobType.allowUndistributedCMSSW = True
        cfg.JobType.maxMemoryMB = 5000
        cfg.JobType.numCores = 2
        cfg.JobType.maxJobRuntimeMin = 2750

        cfg.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

        cfg.Data.inputDataset = dataset
        cfg.Data.outputDatasetTag = f'ScoutingPuppiCalib_{sample}_{VERSION}'
        cfg.Data.outLFNDirBase = '/store/user/jschulte/ScoutingPuppiCalibration/'
        cfg.Data.splitting = 'FileBased'
        cfg.Data.unitsPerJob = UNITS_PER_JOB
        cfg.Data.totalUnits = TOTAL_UNITS
        cfg.Data.ignoreLocality = True
        cfg.Data.publication = False

        cfg.Site.storageSite = 'T2_US_Purdue'
        cfg.Site.whitelist = ['T2_*']
        cfg.Site.blacklist = ['T2_BR_UERJ', 'T2_US_Florida', 'T2_US_Wisconsin', 'T2_US_Caltech', 'T2_US_Nebraska']

        if args.submit:
            try:
                crabCommand('submit', config=cfg)
            except Exception as exc:
                print(f"  can't submit ({exc}); task may already exist")

if args.resubmit:
    project_dirs = sorted(glob.glob(os.path.join('crab_projects', 'crab_scoutingPuppiCalib_*')))
    if not project_dirs:
        print('No existing CRAB projects found to resubmit.')
    for d in project_dirs:
        print(f'Resubmitting {d} ...')
        try:
            crabCommand('resubmit', dir=d, siteblacklist='T2_BR_UERJ,T2_US_Florida,T2_US_Wisconsin,T2_US_Caltech,T2_US_Nebraska')
        except Exception:
            print('  failed to resubmit, most likely there are no failed jobs')

if args.report:
    print('Reporting events processed for scouting-PUPPI-calibration CRAB tasks:')
    total_events = 0
    for sample in mc_samples:
        task_name = f'scoutingPuppiCalib_{sample}_{VERSION}'
        d = os.path.join('crab_projects', f'crab_{task_name}')
        if not os.path.isdir(d):
            print(f'  [{sample}] no CRAB project directory found ({d})')
            continue
        try:
            res = crabCommand('report', dir=d)
            events = res.get('eventsRead', 0)
            total_events += events
            print(f'  [{sample}] {events} events read')
        except Exception as exc:
            print(f'  [{sample}] report failed: {exc}')
    print(f'Total events read (all samples): {total_events}')
