import FWCore.ParameterSet.Config as cms

from ScoutingPuppiCalibration.Calibration.variants import DEFAULT_VARIANTS

# See ScoutingPuppiCalibration/python/variants.py for the variant grid
# itself (kept in a separate, CMSSW-import-free module so analysis scripts
# can read the variant names/params without needing a full CMSSW
# environment). Briefly: "nominal" is the unmodified stock puppi
# configuration (baseline); coneSmall/coneLarge/rmsPtMinLoose probe the
# per-particle alpha-discrimination axis; neutralPtLoose/neutralPtTight probe
# the jet-level neutral-pt-threshold axis. Scouting jets/taggers elsewhere in
# this workflow are restricted to |eta|<=2.5 (entirely the "central" PUPPI
# region), so the forward-region knobs may turn out not to matter for this
# final state -- check that on the first pass before reading much into the
# forward-region comparisons.


def _applyVariantParams(mod, params):
    """Mutate a puppi.clone() in place according to one DEFAULT_VARIANTS entry.

    mod.algos[0] is the central-eta PSet (wraps puppiCentral), mod.algos[1] is
    the forward-eta PSet (wraps puppiForward) -- see CommonTools/PileupAlgos/
    python/Puppi_cff.py. Safe to mutate per-variant because .clone() deep-copies
    the whole algos VPSet.
    """
    if "central_cone" in params:
        mod.algos[0].puppiAlgos[0].cone = cms.double(params["central_cone"])
    if "forward_cone" in params:
        mod.algos[1].puppiAlgos[0].cone = cms.double(params["forward_cone"])
    if "central_rmsPtMin" in params:
        mod.algos[0].puppiAlgos[0].rmsPtMin = cms.double(params["central_rmsPtMin"])
    if "forward_rmsPtMin" in params:
        mod.algos[1].puppiAlgos[0].rmsPtMin = cms.double(params["forward_rmsPtMin"])
    # central_* below mutate algos[0]'s MedEtaSF/RMSEtaSF/MinNeutralPt/
    # MinNeutralPtSlope, which the offline recalibration scans (this
    # session's puppi_recompute.py/run_puppi_recalibration.py) tune as
    # medEtaSF/rmsEtaSF/minNeutralPt/minNeutralPtSlope -- not explicitly set
    # by Puppi_cff.py's algos[0] override (it only sets puppiAlgos=
    # puppiCentral), so these come from PuppiProducer_cfi's registered
    # defaults (1.0/1.0/0.2/0.015) until overridden here, exactly like
    # fwd_MedEtaSF etc. override algos[1]'s cfi-registered forward defaults.
    if "central_MedEtaSF" in params:
        mod.algos[0].MedEtaSF = cms.vdouble(params["central_MedEtaSF"])
    if "central_RMSEtaSF" in params:
        mod.algos[0].RMSEtaSF = cms.vdouble(params["central_RMSEtaSF"])
    if "central_MinNeutralPt" in params:
        mod.algos[0].MinNeutralPt = cms.vdouble(params["central_MinNeutralPt"])
    if "central_MinNeutralPtSlope" in params:
        mod.algos[0].MinNeutralPtSlope = cms.vdouble(params["central_MinNeutralPtSlope"])
    if "fwd_MinNeutralPt" in params:
        mod.algos[1].MinNeutralPt = cms.vdouble(*params["fwd_MinNeutralPt"])
    if "fwd_MinNeutralPtSlope" in params:
        mod.algos[1].MinNeutralPtSlope = cms.vdouble(*params["fwd_MinNeutralPtSlope"])
    if "fwd_RMSEtaSF" in params:
        mod.algos[1].RMSEtaSF = cms.vdouble(*params["fwd_RMSEtaSF"])
    if "fwd_MedEtaSF" in params:
        mod.algos[1].MedEtaSF = cms.vdouble(*params["fwd_MedEtaSF"])


def customiseForScoutingPuppiCalibration(process, pName, variants=None, referenceVariant="optimized"):
    """Build one PUPPI-weighted AK4 jet collection per parameter variant.

    Mirrors the CHS block in PhysicsTools.PatFromScouting.
    scoutingToMiniAODDerivedCollections_cff.customizeForScoutingAK4ReclusteredJets:
    kinematics + JEC + charge + gen-matching only, no tagger rerun (the UParT/
    ParticleNet taggers are trained on un-weighted, non-CHS candidates, so
    PUPPI-reweighted constituents would be yet another train/inference
    mismatch). Requires customizeForScoutingAK4ReclusteredJets(process, pName)
    (i.e. customiseScoutingNanoDerived) to already have been run, since it
    reuses packedPFCandidates, offlineSlimmedPrimaryVertices, scoutingTracks,
    slimmedGenJets, prunedGenParticles and patJetPartonsNano from that step.
    """
    if variants is None:
        variants = DEFAULT_VARIANTS

    from CommonTools.PileupAlgos.Puppi_cff import puppi as _stockPuppi
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
    from RecoJets.JetAssociationProducers.ak4JTA_cff import ak4JetTracksAssociatorAtVertex
    from PhysicsTools.PatAlgos.recoLayer0.jetTracksCharge_cff import patJetCharge
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
    from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import patJetGenJetMatch, patJetPartonMatch
    from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import patJetFlavourAssociation
    from PhysicsTools.NanoAOD.jetMC_cff import patJetPartonsNano, jetMCTable
    from PhysicsTools.NanoAOD.common_cff import Var, P4Vars
    from PhysicsTools.NanoAOD.simplePATJetFlatTableProducer_cfi import simplePATJetFlatTableProducer

    if not hasattr(process, "patJetPartonsNano"):
        process.patJetPartonsNano = patJetPartonsNano

    # Same kinematic/energy-fraction/multiplicity variables used for the
    # plain and CHS jet tables in scoutingToMiniAODDerivedCollections_cff.py,
    # redefined here since that module keeps its copy as a local variable.
    PFJetVariables = cms.PSet(
        P4Vars,
        area = Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
        chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision=10),
        neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision=10),
        chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision=10),
        neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision=10),
        hfHEF = Var("HFHadronEnergyFraction()", float, doc="hadronic Energy Fraction in HF", precision=10),
        hfEmEF = Var("HFEMEnergyFraction()", float, doc="electromagnetic Energy Fraction in HF", precision=10),
        muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision=10),
        chHadMultiplicity = Var("chargedHadronMultiplicity()", "int16", doc="number of charged hadrons in the jet"),
        neHadMultiplicity = Var("neutralHadronMultiplicity()", int, doc="number of neutral hadrons in the jet"),
        hfHadMultiplicity = Var("HFHadronMultiplicity()", int, doc="number of HF hadrons in the jet"),
        hfEMMultiplicity = Var("HFEMMultiplicity()", int, doc="number of HF EMs in the jet"),
        muMultiplicity = Var("muonMultiplicity()", int, doc="number of muons in the jet"),
        elMultiplicity = Var("electronMultiplicity()", int, doc="number of electrons in the jet"),
        phMultiplicity = Var("photonMultiplicity()", int, doc="number of photons in the jet"),
        nConstituents = Var("numberOfDaughters()", int, doc="number of particles in the jet"),
    )

    variantTasks = []
    variantMCTasks = []
    variantTableTasks = []
    variantMCTableTasks = []

    # vtxAssocImproved (see variants.py) needs its own packedPFCandidates
    # instance, built with useImprovedVertexAssociation=True, since that flag
    # changes the recoCands/vtxass products themselves rather than just how
    # PuppiProducer reads them. Built once, shared by any variant that asks
    # for it, and only if such a variant is actually present (this must run
    # after customiseScoutingNanoDerived, which is what defines
    # process.packedPFCandidates in the first place).
    if any(params.get("improvedVertexAssociation") and not params.get("trackMatchedVertexAssociation")
           for params in variants.values()):
        process.packedPFCandidatesImprovedVtxAssoc = process.packedPFCandidates.clone(
            useImprovedVertexAssociation = cms.bool(True),
        )
        process.scoutingPuppiCalibImprovedVtxAssocTask = cms.Task(process.packedPFCandidatesImprovedVtxAssoc)
        process.scoutingNanoSequence.associate(process.scoutingPuppiCalibImprovedVtxAssocTask)

    # vtxAssocTrackMatched (see variants.py) needs its own packedPFCandidates
    # instance too, distinct from packedPFCandidatesImprovedVtxAssoc: both set
    # useImprovedVertexAssociation=True, but only this one also sets
    # useTrackMatchedVertexAssociation=True, so its recoCands/vtxass are
    # generally different (genuine per-vertex dz where a confident track
    # match exists, vtxAssocImproved's linear-shift dz otherwise).
    if any(params.get("trackMatchedVertexAssociation") for params in variants.values()):
        process.packedPFCandidatesTrackMatchedVtxAssoc = process.packedPFCandidates.clone(
            useImprovedVertexAssociation = cms.bool(True),
            useTrackMatchedVertexAssociation = cms.bool(True),
        )
        process.scoutingPuppiCalibTrackMatchedVtxAssocTask = cms.Task(process.packedPFCandidatesTrackMatchedVtxAssoc)
        process.scoutingNanoSequence.associate(process.scoutingPuppiCalibTrackMatchedVtxAssocTask)

    for vname, params in variants.items():
        if params.get("trackMatchedVertexAssociation"):
            candSrcMod = "packedPFCandidatesTrackMatchedVtxAssoc"
        elif params.get("improvedVertexAssociation"):
            candSrcMod = "packedPFCandidatesImprovedVtxAssoc"
        else:
            candSrcMod = "packedPFCandidates"
        # EDM module labels must not contain underscores: outputCommands
        # keep/drop rules split on '_' expecting exactly 4 fields
        # (type_label_instance_process), and a label with embedded
        # underscores breaks that parse (confirmed at cmsRun time -- "keep
        # nanoaodFlatTable_scoutingPFJetReclusterPUPPI_nominal_Table_*_*" was
        # rejected as having too many fields). Use camelCase instead.
        vcap = vname[0].upper() + vname[1:]
        puppiMod = "scoutingPuppi%s" % vcap
        recoJetMod = "recoScoutingPFJetReclusterPUPPI%s" % vcap
        corrMod = "scoutingPFJetReclusterPUPPI%sCorrFactors" % vcap
        tavMod = "scoutingPFJetReclusterPUPPI%sTracksAssociatorAtVertex" % vcap
        chgMod = "scoutingPFJetReclusterPUPPI%sCharge" % vcap
        patMod = "patScoutingPFJetReclusterPUPPI%s" % vcap
        slimMod = "slimmedJetsPUPPI%s" % vcap
        gjmMod = "scoutingPFJetReclusterPUPPI%sGenJetMatch" % vcap
        gpmMod = "scoutingPFJetReclusterPUPPI%sGenPartonMatch" % vcap
        flavMod = "scoutingPFJetReclusterPUPPI%sFlavourAssociation" % vcap
        tblMod = "scoutingPFJetReclusterPUPPI%sTable" % vcap
        mcTblMod = "scoutingPFJetReclusterPUPPI%sMCTable" % vcap

        # --- PuppiProducer instance for this variant ---
        # useVertexAssociation=True + vertexAssociation="packedPFCandidates:vtxass"
        # reuses the Association<VertexCollection> + parallel ValueMap<int> of
        # pat::PackedCandidate::PVAssociationQuality that
        # Run3ScoutingParticleToPackedCandidateProducer already emits under
        # that instance label -- no new C++ needed.
        # useVertexAssociation defaults to True (trust the scouting candidate's
        # precomputed vtxass quality flag); the "vtxAssocOff" variant flips it
        # to False so PuppiProducer instead derives its own dz-threshold-based
        # categorization from the packed candidate's real dz()/fromPV() --
        # see the comment on vtxAssocOff in variants.py. vertexAssociation is
        # harmless to leave set either way: PuppiProducer only consumes()/
        # reads it when useVertexAssociation=True (PuppiProducer.cc L117-121).
        puppiClone = _stockPuppi.clone(
            candName = cms.InputTag(candSrcMod, "recoCands", pName),
            vertexName = cms.InputTag("offlineSlimmedPrimaryVertices", "", pName),
            useVertexAssociation = params.get("useVertexAssociation", True),
            # required alongside useVertexAssociation=False -- see the long
            # comment on vtxAssocOff in variants.py for why (without it,
            # charged pileup rejection is effectively disabled for scouting's
            # |eta|<2.5, mostly-soft candidates).
            UseFromPVLooseTight = params.get("useFromPVLooseTight", False),
            vertexAssociation = cms.InputTag(candSrcMod, "vtxass", pName),
            clonePackedCands = False,
            puppiDiagnostics = (vname == referenceVariant),
            # PtMaxCharged: stock default is -1 (disabled). PuppiProducer.cc's
            # fUseVertexAssociation branch (the one every variant except
            # vtxAssocOff runs) was patched this session to read this same
            # top-level v15 knob -- see the "v15-style protection" comment
            # there -- so setting it here unconditionally protects high-pT
            # charged candidates the scouting vertex fit assigned to PU or
            # failed to associate at all, mirroring PUPPI v15 (DP-2021/001).
            # UseFromPV2Recovery/PtMinForFromPV2Recovery need no override:
            # Puppi_cff.py already sets them True/4. on _stockPuppi, so the
            # softer pT-floor half of that same patch is already active for
            # every variant below without any params.py change.
            PtMaxCharged = params.get("ptMaxCharged", -1.),
        )
        _applyVariantParams(puppiClone, params)
        setattr(process, puppiMod, puppiClone)

        # --- reclustered AK4 jets: full (non-CHS) candidates, PUPPI-weighted ---
        setattr(process, recoJetMod, ak4PFJets.clone(
            src = (candSrcMod, "recoCands", pName),
            applyWeight = True,
            srcWeights = cms.InputTag(puppiMod),
            jetPtMin = 20,
        ))

        # --- JEC: reuse AK4PFHLT, same approximation the CHS jets already use ---
        setattr(process, corrMod, patJetCorrFactors.clone(
            src = recoJetMod,
            levels = cms.vstring(
                "L1FastJet",
                "L2Relative",
                "L3Absolute",
                "L2L3Residual"),
            payload = cms.string("AK4PFHLT"),
            primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices", "", pName),
        ))

        setattr(process, tavMod, ak4JetTracksAssociatorAtVertex.clone(
            jets = cms.InputTag(recoJetMod),
            coneSize = cms.double(0.4),
            tracks = cms.InputTag("scoutingTracks"),
            pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices", "", pName),
        ))
        setattr(process, chgMod, patJetCharge.clone(
            src = cms.InputTag(tavMod),
        ))

        # --- PAT-ify: kinematics/JEC/charge/gen-matching only, no tagger rerun ---
        setattr(process, patMod, _patJets.clone(
            jetSource = recoJetMod,
            addJetCorrFactors = True,
            jetCorrFactorsSource = [corrMod],
            addBTagInfo = False,
            addDiscriminators = False,
            addAssociatedTracks = False,
            addJetCharge = True,
            jetChargeSource = chgMod,
            addGenPartonMatch = True,
            embedGenPartonMatch = True,
            genPartonMatch = cms.InputTag(gpmMod),
            addGenJetMatch = True,
            embedGenJetMatch = True,
            genJetMatch = cms.InputTag(gjmMod),
            getJetMCFlavour = True,
            useLegacyJetMCFlavour = False,
            addJetFlavourInfo = True,
            JetFlavourInfoSource = cms.InputTag(flavMod),
        ))

        # --- slim: PUPPI does not drop constituents -> rekey against the
        # plain (non-CHS) packedPFCandidates, not packedPFCandidatesCHS ---
        setattr(process, slimMod, cms.EDProducer("PATJetSlimmer",
            src = cms.InputTag(patMod),
            packedPFCandidates = cms.InputTag(candSrcMod, "", pName),
            dropJetVars = cms.string("1"),
            dropDaughters = cms.string("0"),
            rekeyDaughters = cms.string("1"),
            dropTrackRefs = cms.string("1"),
            dropSpecific = cms.string("0"),
            dropTagInfos = cms.string("1"),
            modifyJets = cms.bool(True),
            mixedDaughters = cms.bool(False),
            modifierConfig = cms.PSet(modifications = cms.VPSet())
        ))

        # --- gen matching, same pattern as the plain/CHS blocks ---
        setattr(process, gjmMod, patJetGenJetMatch.clone(
            src = cms.InputTag(recoJetMod),
            matched = cms.InputTag("slimmedGenJets"),
            resolveByMatchQuality = cms.bool(True)
        ))
        setattr(process, gpmMod, patJetPartonMatch.clone(
            src = cms.InputTag(recoJetMod),
            matched = cms.InputTag("prunedGenParticles"),
        ))
        setattr(process, flavMod, patJetFlavourAssociation.clone(
            jets = cms.InputTag(recoJetMod),
            # required whenever the jet collection was clustered with
            # applyWeight/srcWeights (PUPPI): JetFlavourClustering needs the
            # same per-constituent weights for its internal ghost clustering,
            # otherwise it throws "No weights (e.g. PUPPI) given for weighted
            # jet collection" -- same pattern as stock
            # PhysicsTools/PatAlgos/python/tools/puppiJetMETReclusteringTools.py
            # (patJetFlavourAssociationPuppi.weights = cms.InputTag(puppiLabel)).
            weights = cms.InputTag(puppiMod),
            bHadrons = cms.InputTag("patJetPartonsNano", "bHadrons"),
            cHadrons = cms.InputTag("patJetPartonsNano", "cHadrons"),
            partons = cms.InputTag("patJetPartonsNano", "physicsPartons"),
            leptons = cms.InputTag("patJetPartonsNano", "leptons"),
        ))

        # --- NanoAOD tables ---
        setattr(process, tblMod, simplePATJetFlatTableProducer.clone(
            src = cms.InputTag(patMod),
            name = cms.string("ScoutingPFJetReclusterPUPPI_%s" % vname),
            doc = cms.string("AK4 scouting jets reclustered with PUPPI weighting, variant=%s "
                              "(JEC payload AK4PFHLT is an approximation, no dedicated "
                              "scouting-PUPPI payload exists)" % vname),
            cut = cms.string(""),
            variables = cms.PSet(
                PFJetVariables,
                rawFactor = Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=10),
                charge = Var("jetCharge()", float, doc="charge", precision=10),
            ),
        ))
        setattr(process, mcTblMod, jetMCTable.clone(
            src = cms.InputTag(patMod),
            name = getattr(process, tblMod).name,
            cut = getattr(process, tblMod).cut,
        ))

        variantTasks.append(cms.Task(
            getattr(process, puppiMod),
            getattr(process, recoJetMod),
            getattr(process, corrMod),
            getattr(process, tavMod),
            getattr(process, chgMod),
            getattr(process, patMod),
            getattr(process, slimMod),
        ))
        variantMCTasks.append(cms.Task(
            getattr(process, gjmMod),
            getattr(process, gpmMod),
            getattr(process, flavMod),
        ))
        variantTableTasks.append(cms.Task(getattr(process, tblMod)))
        variantMCTableTasks.append(cms.Task(getattr(process, mcTblMod)))

    process.scoutingPuppiCalibJetTask = cms.Task(*variantTasks)
    process.scoutingPuppiCalibTableTask = cms.Task(*variantTableTasks)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibJetTask)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibTableTask)

    runOnMC = hasattr(process, "NANOEDMAODSIMoutput") or hasattr(process, "NANOAODSIMoutput")
    if runOnMC:
        process.scoutingPuppiCalibJetMCTask = cms.Task(*variantMCTasks)
        process.scoutingPuppiCalibMCTableTask = cms.Task(*variantMCTableTasks)
        process.scoutingNanoSequence.associate(process.scoutingPuppiCalibJetMCTask)
        process.scoutingNanoSequence.associate(process.scoutingPuppiCalibMCTableTask)

    _addCandidateTable(process, pName, variants)
    for vname, params in variants.items():
        if params.get("trackMatchedVertexAssociation"):
            _addTrackMatchedCandidateTable(process, pName, vname)
        elif params.get("improvedVertexAssociation"):
            _addImprovedCandidateTable(process, pName, vname)
    _addPuppiAlphaDiagnostics(process, pName, referenceVariant, variants)
    _addPileupTruthTable(process, pName)
    _addGenJetTable(process)
    _addRhoTable(process)

    return process


def _addCandidateTable(process, pName, variants):
    """Shared per-candidate table: kinematics + the (online-proxy) vertex-
    association quality + every variant's PUPPI weight as sibling columns.

    All three are keyed to the same packedPFCandidates:recoCands collection/
    index order (candidate kinematics directly, pvAssocQuality via the
    "vtxass" ValueMap<int>, each PUPPI weight via PuppiProducer's own
    ValueMap<float> output), so this is a single flat table with no
    Association-hopping between reco::PFCandidate and pat::PackedCandidate
    needed. SimplePFCandidateFlatTableProducer (PhysicsTools/NanoAOD/plugins/
    SimpleFlatTableProducerPlugins.cc) is templated on reco::PFCandidate and
    so reads packedPFCandidates:recoCands directly.
    """
    from PhysicsTools.NanoAOD.common_cff import CandVars, ExtVar

    externalVariables = dict(
        pvAssocQuality = ExtVar(
            cms.InputTag("packedPFCandidates", "vtxass", pName), int,
            doc="pat::PackedCandidate::PVAssociationQuality, derived from the HLT "
                "nearest-pixel-vertex match (particle.vertex()); an ONLINE proxy for "
                "pileup origin, NOT MC truth. 0=NotReconstructedPrimary, "
                "5=CompatibilityDz (associated to a pileup vertex), "
                "7=UsedInFitTight (associated to the leading/PV vertex)"),
    )
    for vname, params in variants.items():
        if params.get("improvedVertexAssociation") or params.get("trackMatchedVertexAssociation"):
            # its PuppiProducer weight ValueMap is keyed to
            # packedPFCandidatesImprovedVtxAssoc:recoCands or
            # packedPFCandidatesTrackMatchedVtxAssoc:recoCands, a *different*
            # product than this table's src (packedPFCandidates:recoCands) --
            # attaching it here would throw "ValueMap: no associated value
            # for given product and index" at runtime. It gets its own table
            # instead, see _addImprovedCandidateTable/_addTrackMatchedCandidateTable.
            continue
        vcap = vname[0].upper() + vname[1:]
        # NOTE: "puppiWeight_<variant>" here is a PSet member name / NanoAOD
        # column name, not an EDM module label, so the underscore is fine --
        # only the InputTag's module-label part (built the same way as
        # puppiMod above) must be underscore-free.
        externalVariables["puppiWeight_%s" % vname] = ExtVar(
            cms.InputTag("scoutingPuppi%s" % vcap), float,
            doc="PUPPI weight, variant=%s" % vname, precision=12)

    process.scoutingPuppiCalibCandTable = cms.EDProducer("SimplePFCandidateFlatTableProducer",
        src = cms.InputTag("packedPFCandidates", "recoCands", pName),
        cut = cms.string(""),
        name = cms.string("ScoutingPuppiCalibCand"),
        doc = cms.string("Scouting PF candidates used as PUPPI input for the calibration study; "
                          "pvAssocQuality is an ONLINE HLT-vertex-matching proxy, NOT MC truth "
                          "(no TrackingParticle/TrackingVertex truth link exists for scouting)"),
        singleton = cms.bool(False),
        extension = cms.bool(False),
        variables = cms.PSet(CandVars),
        externalVariables = cms.PSet(**externalVariables),
    )
    process.scoutingPuppiCalibCandTask = cms.Task(process.scoutingPuppiCalibCandTable)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibCandTask)


def _addImprovedCandidateTable(process, pName, vname):
    """Second, small candidate table for variants with improvedVertexAssociation=True.

    packedPFCandidatesImprovedVtxAssoc:recoCands is a *different* product
    instance than packedPFCandidates:recoCands (same underlying particles,
    freshly reprocessed by a second Run3ScoutingParticleToPackedCandidate
    Producer with useImprovedVertexAssociation=True), so its vtxass ValueMap
    and this variant's PuppiProducer weight ValueMap can't be attached as
    ExternalVariables on the main scoutingPuppiCalibCandTable (src=
    packedPFCandidates:recoCands) -- an EDM Association/ValueMap lookup
    requires the value map's product ID to match the src collection's,
    and these come from a different producer instance. Hence a standalone
    table, keyed to packedPFCandidatesImprovedVtxAssoc:recoCands directly,
    so the vtxAssocOff-style before/after pvAssocQuality-vs-puppiWeight
    sanity check (analysis/candidates.py) can be rerun against it.
    """
    from PhysicsTools.NanoAOD.common_cff import CandVars, ExtVar

    vcap = vname[0].upper() + vname[1:]
    puppiMod = "scoutingPuppi%s" % vcap
    tblName = "ScoutingPuppiCalibImprovedVtxAssocCand"

    process.scoutingPuppiCalibImprovedVtxAssocCandTable = cms.EDProducer("SimplePFCandidateFlatTableProducer",
        src = cms.InputTag("packedPFCandidatesImprovedVtxAssoc", "recoCands", pName),
        cut = cms.string(""),
        name = cms.string(tblName),
        doc = cms.string("Scouting PF candidates as re-associated by useImprovedVertexAssociation=True "
                          "(variant=%s), for direct before/after comparison against ScoutingPuppiCalibCand's "
                          "pvAssocQuality" % vname),
        singleton = cms.bool(False),
        extension = cms.bool(False),
        variables = cms.PSet(CandVars),
        externalVariables = cms.PSet(
            pvAssocQuality = ExtVar(
                cms.InputTag("packedPFCandidatesImprovedVtxAssoc", "vtxass", pName), int,
                doc="pat::PackedCandidate::PVAssociationQuality from the improved, dz-scan-based "
                    "association -- see useImprovedVertexAssociation in Run3ScoutingParticleTo"
                    "PackedCandidateProducer.cc"),
            **{"puppiWeight_%s" % vname: ExtVar(
                cms.InputTag(puppiMod), float,
                doc="PUPPI weight, variant=%s" % vname, precision=12)},
        ),
    )
    process.scoutingPuppiCalibImprovedVtxAssocCandTask = cms.Task(process.scoutingPuppiCalibImprovedVtxAssocCandTable)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibImprovedVtxAssocCandTask)


def _addTrackMatchedCandidateTable(process, pName, vname):
    """Third candidate table, for variants with trackMatchedVertexAssociation=True.

    Same rationale as _addImprovedCandidateTable: packedPFCandidatesTrackMatched
    VtxAssoc:recoCands is yet another distinct product instance, so its vtxass
    ValueMap and this variant's PuppiProducer weight ValueMap need their own
    table, keyed to that instance directly.
    """
    from PhysicsTools.NanoAOD.common_cff import CandVars, ExtVar

    vcap = vname[0].upper() + vname[1:]
    puppiMod = "scoutingPuppi%s" % vcap
    tblName = "ScoutingPuppiCalibTrackMatchedVtxAssocCand"

    process.scoutingPuppiCalibTrackMatchedVtxAssocCandTable = cms.EDProducer("SimplePFCandidateFlatTableProducer",
        src = cms.InputTag("packedPFCandidatesTrackMatchedVtxAssoc", "recoCands", pName),
        cut = cms.string(""),
        name = cms.string(tblName),
        doc = cms.string("Scouting PF candidates as re-associated by useImprovedVertexAssociation=True "
                          "+ useTrackMatchedVertexAssociation=True (variant=%s), for direct before/after "
                          "comparison against ScoutingPuppiCalibCand's/ScoutingPuppiCalibImprovedVtxAssoc"
                          "Cand's pvAssocQuality" % vname),
        singleton = cms.bool(False),
        extension = cms.bool(False),
        variables = cms.PSet(CandVars),
        externalVariables = cms.PSet(
            pvAssocQuality = ExtVar(
                cms.InputTag("packedPFCandidatesTrackMatchedVtxAssoc", "vtxass", pName), int,
                doc="pat::PackedCandidate::PVAssociationQuality from the track-matched, genuine-dz "
                    "association -- see useTrackMatchedVertexAssociation in Run3ScoutingParticleTo"
                    "PackedCandidateProducer.cc"),
            **{"puppiWeight_%s" % vname: ExtVar(
                cms.InputTag(puppiMod), float,
                doc="PUPPI weight, variant=%s" % vname, precision=12)},
        ),
    )
    process.scoutingPuppiCalibTrackMatchedVtxAssocCandTask = cms.Task(process.scoutingPuppiCalibTrackMatchedVtxAssocCandTable)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibTrackMatchedVtxAssocCandTask)


def _addPuppiAlphaDiagnostics(process, pName, referenceVariant, variants):
    """Expose PuppiProducer's puppiDiagnostics=True output (enabled only for
    referenceVariant -- see puppiDiagnostics=(vname==referenceVariant) above)
    as per-candidate columns on referenceVariant's own candidate table:
    puppiRawAlpha (this candidate's own alpha, the value PuppiAlgo::compute()
    compares to the event's median/RMS to derive the weight -- not the
    weight itself) and puppiAlphaMed/puppiAlphaRms (that event's per-algo
    median/RMS, broadcast to every candidate for convenience).
    PuppiRawAlphaToValueMapProducer (ScoutingPuppiCalibration/Calibration/
    plugins) does the actual vector->ValueMap conversion PuppiProducer's raw
    diagnostics products need, since they're bare event-level vector<double>
    products, not ValueMaps keyed to the candidate collection.

    With these three numbers plus a candidate's own pt, PuppiAlgo::compute()'s
    chi2 = (alpha-med)*|alpha-med|/rms^2, weight = chisquared_cdf(chi2, ndof)
    can be recomputed offline for any hypothetical MedEtaSF/RMSEtaSF/
    MinNeutralPt/MinNeutralPtSlope -- see the long comment in
    PuppiRawAlphaToValueMapProducer.cc for why those four parameters (unlike
    cone/rmsPtMin, which affect the alpha computation itself) never require
    rerunning PuppiProducer to retune.
    """
    from PhysicsTools.NanoAOD.common_cff import ExtVar

    if referenceVariant not in variants:
        return
    refParams = variants[referenceVariant]
    if refParams.get("trackMatchedVertexAssociation"):
        candSrcMod = "packedPFCandidatesTrackMatchedVtxAssoc"
        tableAttr = "scoutingPuppiCalibTrackMatchedVtxAssocCandTable"
    elif refParams.get("improvedVertexAssociation"):
        candSrcMod = "packedPFCandidatesImprovedVtxAssoc"
        tableAttr = "scoutingPuppiCalibImprovedVtxAssocCandTable"
    else:
        candSrcMod = "packedPFCandidates"
        tableAttr = "scoutingPuppiCalibCandTable"

    vcap = referenceVariant[0].upper() + referenceVariant[1:]
    puppiMod = "scoutingPuppi%s" % vcap

    process.puppiAlphaDiagnostics = cms.EDProducer("PuppiRawAlphaToValueMapProducer",
        candidates = cms.InputTag(candSrcMod, "recoCands", pName),
        rawAlphas = cms.InputTag(puppiMod, "PuppiRawAlphas"),
        alphasMed = cms.InputTag(puppiMod, "PuppiAlphasMed"),
        alphasRms = cms.InputTag(puppiMod, "PuppiAlphasRms"),
        # ascending |eta| upper edges of every PuppiRawAlphas block but the
        # last -- {2.5} matches Puppi_cff.py's production "algos" VPSet
        # (algos[0]=central |eta|<2.5, algos[1]=forward |eta|>=2.5, covering
        # both its 2.5-3.0/3.0-10.0 sub-bins with one shared raw-alpha block
        # since they share cone/rmsPtMin -- see PuppiRawAlphaToValueMapProducer.cc).
        etaBoundaries = cms.vdouble(2.5),
    )
    process.scoutingPuppiCalibAlphaDiagnosticsTask = cms.Task(process.puppiAlphaDiagnostics)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibAlphaDiagnosticsTask)

    table = getattr(process, tableAttr)
    table.externalVariables.puppiRawAlpha = ExtVar(
        cms.InputTag("puppiAlphaDiagnostics", "rawAlpha"), float,
        doc="PUPPI's internal per-candidate shape variable (alpha) for variant=%s, the value "
            "PuppiAlgo::compute() compares to the event's median/RMS (puppiAlphaMed/puppiAlphaRms) "
            "to derive the weight -- NOT the weight itself. Correctly picks the central "
            "(|eta|<2.5, cone=0.4/rmsPtMin=0.1) or forward (|eta|>=2.5, cone=0.4/rmsPtMin=0.5) "
            "raw-alpha block per candidate's own eta." % referenceVariant, precision=12)
    table.externalVariables.puppiAlphaMed = ExtVar(
        cms.InputTag("puppiAlphaDiagnostics", "alphaMed"), float,
        doc="That event's median alpha for this candidate's own eta-region algo (already "
            "MedEtaSF-scaled: central and forward candidates get different values, and the two "
            "forward MedEtaSF sub-bins, 0.90 for 2.5-3.0 / 0.75 for 3.0-10.0, are both already "
            "correctly applied here even though rawAlpha's forward block is shared); divide out "
            "the configured MedEtaSF to recover the raw median.", precision=12)
    table.externalVariables.puppiAlphaRms = ExtVar(
        cms.InputTag("puppiAlphaDiagnostics", "alphaRms"), float,
        doc="Same as puppiAlphaMed but the RMS (RMSEtaSF-scaled); divide out the configured "
            "RMSEtaSF to recover the raw RMS.",
        precision=12)


def _addPileupTruthTable(process, pName):
    """Wire in Pileup_nTrueInt (not present in the stock scouting NanoAOD
    sequence) so jet response/resolution can be binned vs true pileup.
    """
    if hasattr(process, "puTable"):
        return
    from PhysicsTools.NanoAOD.globals_cff import puTable
    process.puTable = puTable.clone(
        pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices", "", pName),
    )
    process.scoutingPuppiCalibPUTask = cms.Task(process.puTable)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibPUTask)


def _addRhoTable(process):
    """Wire in the standard Rho NanoAOD table (fixedGridRhoFastjetAll etc),
    not present in the stock scouting NanoAOD sequence. Needed so an offline
    FastJet reclustering of the raw candidates (analysis/puppi_refit.py) can
    apply the same L1FastJet JEC term (which needs rho) production's
    patJetCorrFactors uses (default rho source, see PhysicsTools/PatAlgos/
    python/recoLayer0/jetCorrFactors_cfi.py's useRho/rho='fixedGridRhoFastjetAll',
    unmodified by scoutingToMiniAODDerivedCollections_cff.py) -- fixedGridRhoFastjetAll
    itself is a standard MiniAOD product (not scouting-specific), already
    present in the input file; only the NanoAOD table exposing it was missing.
    """
    if hasattr(process, "rhoTable"):
        return
    from PhysicsTools.NanoAOD.globals_cff import rhoTable
    process.rhoTable = rhoTable.clone()
    process.scoutingPuppiCalibRhoTask = cms.Task(process.rhoTable)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibRhoTask)


def _addGenJetTable(process):
    """Wire in the standard GenJet NanoAOD table (src=slimmedGenJets).

    Nothing in the stock scouting NanoAOD workflow (custom_run3scouting_cff)
    or in scoutingToMiniAODDerivedCollections_cff.py adds this table -- the
    latter only uses slimmedGenJets as an EDM product (patJetGenJetMatch's
    "matched" collection), never exposes its kinematics as a branch. Without
    it, the genJetIdx column on each jet table has nothing to index into
    from a flat-table analysis, so this is required for the response/
    resolution comparisons, not optional.
    """
    if hasattr(process, "genJetTable"):
        return
    from PhysicsTools.NanoAOD.jetMC_cff import genJetTable
    process.genJetTable = genJetTable.clone()
    process.scoutingPuppiCalibGenJetTask = cms.Task(process.genJetTable)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibGenJetTask)


def _addAK8Jets(process, pName, variants, referenceVariant):
    """AK8 validation for the PUPPI recalibration: wires in the existing
    plain-AK8 reclustering (PhysicsTools.PatFromScouting.
    scoutingToMiniAODDerivedCollections_cff.customizeForScoutingAK8Reclustered
    Jets, already in the codebase but never called by customiseScoutingNano
    Derived or exposed as a NanoAOD table -- confirmed by grepping this repo
    for "FatJet"/"fatJet", no hits) and adds one PUPPI-weighted AK8
    collection per variant alongside it.

    The PUPPI AK8 collections reuse the SAME scoutingPuppi<Variant>
    PuppiProducer instances customiseForScoutingPuppiCalibration already
    built for AK4 -- PUPPI weights are a per-candidate quantity, independent
    of which jet clustering radius consumes them, so no new PuppiProducer is
    needed, only a differently-configured jet clustering step. Must run
    AFTER customiseForScoutingPuppiCalibration for that reason.

    jetPtMin is left at the existing plain-AK8 default (170 GeV, matching
    the boosted-object trigger thresholds these jets are for) rather than
    lowered for broader validation coverage -- these jets are only needed
    for boosted particle decays, not a general AK8 validation sample.
    """
    from PhysicsTools.PatFromScouting.scoutingToMiniAODDerivedCollections_cff import (
        customizeForScoutingAK8ReclusteredJets,
    )
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
    from RecoJets.JetAssociationProducers.ak4JTA_cff import ak4JetTracksAssociatorAtVertex
    from PhysicsTools.PatAlgos.recoLayer0.jetTracksCharge_cff import patJetCharge
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
    from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import patJetGenJetMatch
    from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import patJetFlavourAssociation
    from PhysicsTools.NanoAOD.jetMC_cff import patJetPartonsNano, fatJetMCTable, genJetAK8Table
    from PhysicsTools.NanoAOD.common_cff import Var, P4Vars
    from PhysicsTools.NanoAOD.simplePATJetFlatTableProducer_cfi import simplePATJetFlatTableProducer

    if variants is None:
        variants = DEFAULT_VARIANTS
    if not hasattr(process, "patJetPartonsNano"):
        process.patJetPartonsNano = patJetPartonsNano

    # AK8 gen truth: NOT the same table/branch as AK4's GenJet -- fatJetMCTable's
    # genJetAK8Idx (below) is gated at pt>100 to match this table's own cut
    # (both hardcoded in PhysicsTools/NanoAOD/python/jetMC_cff.py), unlike
    # AK4's jetMCTable/GenJetTable pair which gate at pt>10. Using the AK4
    # jetMCTable for AK8 collections (an earlier version of this function did)
    # would silently produce a genJetIdx branch pointing past the end of a
    # GenJetAK8 table that doesn't exist -- must use the matched pair.
    if not hasattr(process, "genJetAK8Table"):
        process.genJetAK8Table = genJetAK8Table.clone()
        process.scoutingPuppiCalibGenJetAK8Task = cms.Task(process.genJetAK8Table)
        process.scoutingNanoSequence.associate(process.scoutingPuppiCalibGenJetAK8Task)

    AK8JetVariables = cms.PSet(
        P4Vars,
        area = Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
        chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision=10),
        neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision=10),
        chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision=10),
        neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision=10),
        nConstituents = Var("numberOfDaughters()", int, doc="number of particles in the jet"),
    )

    # --- plain AK8: wire in the existing (never-called) customization ---
    process = customizeForScoutingAK8ReclusteredJets(process, pName)
    process.scoutingNanoSequence.associate(process.scoutingFatPFJetRecluster2Task)
    process.scoutingNanoSequence.associate(process.scoutingFatPFJetRecluster2MCTask)

    process.scoutingFatPFJetRecluster2Table = simplePATJetFlatTableProducer.clone(
        src = cms.InputTag("slimmedJetsAK8"),
        name = cms.string("ScoutingFatPFJetRecluster2"),
        doc = cms.string("AK8 scouting jets reclustered from plain (non-CHS, non-PUPPI) candidates, "
                          "for boosted-decay validation"),
        cut = cms.string(""),
        variables = cms.PSet(
            AK8JetVariables,
            rawFactor = Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=10),
            charge = Var("jetCharge()", float, doc="charge", precision=10),
        ),
    )
    process.scoutingFatPFJetRecluster2MCTable = fatJetMCTable.clone(
        src = cms.InputTag("slimmedJetsAK8"),
        name = process.scoutingFatPFJetRecluster2Table.name,
        cut = process.scoutingFatPFJetRecluster2Table.cut,
    )
    process.scoutingPuppiCalibAK8TableTask = cms.Task(
        process.scoutingFatPFJetRecluster2Table, process.scoutingFatPFJetRecluster2MCTable,
    )
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibAK8TableTask)

    # --- softdrop-groomed companion collections (for msoftdrop, boosted-W/H
    # decay validation e.g. WWW). Recipe copied from RecoJets/JetProducers/
    # python/ak8PFJets_cfi.py's ak8PFJetsCHSSoftDrop/ak8PFJetsPuppiSoftDrop
    # (same underlying FastjetJetProducer as our own AK8 clustering, just
    # useSoftDrop=True + the standard zcut=0.1/beta=0.0/R0=0.8 CMS Run3
    # working point), WITHOUT that recipe's writeCompound=True/SubJets
    # output -- we only want the scalar groomed jet mass, not per-subjet
    # substructure, so plain (non-compound) output already has the groomed
    # mass directly in .mass(), no addSubjets()/groomedMass() indirection
    # needed. A companion table, own pt/eta/phi/mass only (no JEC/tagging) --
    # deliberately NOT index- or dR-matched back onto the ungroomed AK8
    # table here (would need a dedicated reco-to-reco matcher, since
    # PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi's GenJetMatcher only
    # accepts a GenJetCollection on the "matched" side); do that offline via
    # a simple deltaR nearest-match instead, analysis-side, the same way
    # puppi_refit.py's _match_to_gen already does for gen matching.
    def _addSoftDropTable(label, src, applyWeight, srcWeights, doc):
        recoMod = "reco%sSoftDrop" % label
        patMod = "pat%sSoftDrop" % label
        tblMod = "scouting%sSoftDropTable" % label
        setattr(process, recoMod, ak4PFJets.clone(
            src = src,
            applyWeight = applyWeight,
            srcWeights = srcWeights if srcWeights else cms.InputTag(""),
            rParam = 0.8,
            jetPtMin = 170.0,
            useSoftDrop = cms.bool(True),
            zcut = cms.double(0.1),
            beta = cms.double(0.0),
            R0 = cms.double(0.8),
        ))
        setattr(process, patMod, _patJets.clone(
            jetSource = recoMod,
            addJetCorrFactors = False,
            addBTagInfo = False,
            addDiscriminators = False,
            addAssociatedTracks = False,
            addJetCharge = False,
            addGenPartonMatch = False,
            addGenJetMatch = False,
            getJetMCFlavour = False,
            addJetFlavourInfo = False,
        ))
        setattr(process, tblMod, simplePATJetFlatTableProducer.clone(
            src = cms.InputTag(patMod),
            name = cms.string("Scouting%sSoftDrop" % label),
            doc = cms.string(doc),
            cut = cms.string(""),
            variables = cms.PSet(P4Vars),
        ))
        task = cms.Task(getattr(process, recoMod), getattr(process, patMod), getattr(process, tblMod))
        setattr(process, "scoutingPuppiCalib%sSoftDropTask" % label, task)
        process.scoutingNanoSequence.associate(task)

    _addSoftDropTable(
        "FatPFJetRecluster2", cms.InputTag("packedPFCandidates", "recoCands", pName),
        False, None,
        "Soft-drop-groomed mass (zcut=0.1, beta=0.0, R0=0.8) of the plain AK8 collection, "
        "own pt/eta/phi/mass only -- match to ScoutingFatPFJetRecluster2 by deltaR offline "
        "(not index-aligned by construction).",
    )

    # --- PUPPI AK8, one per variant, reusing the AK4 block's PuppiProducer ---
    variantTasks, variantMCTasks, variantTableTasks, variantMCTableTasks = [], [], [], []
    for vname, params in variants.items():
        if params.get("trackMatchedVertexAssociation"):
            candSrcMod = "packedPFCandidatesTrackMatchedVtxAssoc"
        elif params.get("improvedVertexAssociation"):
            candSrcMod = "packedPFCandidatesImprovedVtxAssoc"
        else:
            candSrcMod = "packedPFCandidates"
        vcap = vname[0].upper() + vname[1:]
        puppiMod = "scoutingPuppi%s" % vcap  # already built by customiseForScoutingPuppiCalibration

        _addSoftDropTable(
            "FatPFJetReclusterPUPPI%s" % vcap, cms.InputTag(candSrcMod, "recoCands", pName),
            True, cms.InputTag(puppiMod),
            "Soft-drop-groomed mass (zcut=0.1, beta=0.0, R0=0.8) of the PUPPI-weighted AK8 "
            "collection, variant=%s, own pt/eta/phi/mass only -- match to "
            "ScoutingFatPFJetReclusterPUPPI_%s by deltaR offline (not index-aligned by "
            "construction)." % (vname, vname),
        )

        recoJetMod = "recoScoutingFatPFJetReclusterPUPPI%s" % vcap
        corrMod = "scoutingFatPFJetReclusterPUPPI%sCorrFactors" % vcap
        tavMod = "scoutingFatPFJetReclusterPUPPI%sTracksAssociatorAtVertex" % vcap
        chgMod = "scoutingFatPFJetReclusterPUPPI%sCharge" % vcap
        patMod = "patScoutingFatPFJetReclusterPUPPI%s" % vcap
        gjmMod = "scoutingFatPFJetReclusterPUPPI%sGenJetMatch" % vcap
        flavMod = "scoutingFatPFJetReclusterPUPPI%sFlavourAssociation" % vcap
        tblMod = "scoutingFatPFJetReclusterPUPPI%sTable" % vcap
        mcTblMod = "scoutingFatPFJetReclusterPUPPI%sMCTable" % vcap

        setattr(process, recoJetMod, ak4PFJets.clone(
            src = (candSrcMod, "recoCands", pName),
            applyWeight = True,
            srcWeights = cms.InputTag(puppiMod),
            rParam = 0.8,
            jetPtMin = 170.0,
        ))
        setattr(process, corrMod, patJetCorrFactors.clone(
            src = recoJetMod,
            levels = cms.vstring("L1FastJet", "L2Relative", "L3Absolute", "L2L3Residual"),
            payload = cms.string("AK8PFHLT"),
            primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices", "", pName),
        ))
        # coneSize=0.4 (not 0.8) matches the plain-AK8 tracksAssociator's own
        # existing choice in customizeForScoutingAK8ReclusteredJets -- kept
        # consistent with that precedent rather than deviating for this
        # parallel PUPPI collection.
        setattr(process, tavMod, ak4JetTracksAssociatorAtVertex.clone(
            jets = cms.InputTag(recoJetMod),
            coneSize = cms.double(0.4),
            tracks = cms.InputTag("scoutingTracks"),
            pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices", "", pName),
        ))
        setattr(process, chgMod, patJetCharge.clone(src = cms.InputTag(tavMod)))

        setattr(process, patMod, _patJets.clone(
            jetSource = recoJetMod,
            addJetCorrFactors = True,
            jetCorrFactorsSource = [corrMod],
            addBTagInfo = False,
            addDiscriminators = False,
            addAssociatedTracks = False,
            addJetCharge = True,
            jetChargeSource = chgMod,
            addGenPartonMatch = False,
            embedGenPartonMatch = False,
            addGenJetMatch = True,
            embedGenJetMatch = True,
            genJetMatch = cms.InputTag(gjmMod),
            getJetMCFlavour = True,
            useLegacyJetMCFlavour = False,
            addJetFlavourInfo = True,
            JetFlavourInfoSource = cms.InputTag(flavMod),
        ))

        setattr(process, gjmMod, patJetGenJetMatch.clone(
            src = cms.InputTag(recoJetMod),
            matched = cms.InputTag("slimmedGenJetsAK8"),
            resolveByMatchQuality = cms.bool(True),
        ))
        setattr(process, flavMod, patJetFlavourAssociation.clone(
            jets = cms.InputTag(recoJetMod),
            rParam = cms.double(0.8),
            weights = cms.InputTag(puppiMod),
            bHadrons = cms.InputTag("patJetPartonsNano", "bHadrons"),
            cHadrons = cms.InputTag("patJetPartonsNano", "cHadrons"),
            partons = cms.InputTag("patJetPartonsNano", "physicsPartons"),
            leptons = cms.InputTag("patJetPartonsNano", "leptons"),
        ))

        setattr(process, tblMod, simplePATJetFlatTableProducer.clone(
            src = cms.InputTag(patMod),
            name = cms.string("ScoutingFatPFJetReclusterPUPPI_%s" % vname),
            doc = cms.string("AK8 scouting jets reclustered with PUPPI weighting, variant=%s, for "
                              "boosted-decay validation (JEC payload AK8PFHLT is an approximation, "
                              "no dedicated scouting-PUPPI-AK8 payload exists)" % vname),
            cut = cms.string(""),
            variables = cms.PSet(
                AK8JetVariables,
                rawFactor = Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=10),
                charge = Var("jetCharge()", float, doc="charge", precision=10),
            ),
        ))
        setattr(process, mcTblMod, fatJetMCTable.clone(
            src = cms.InputTag(patMod),
            name = getattr(process, tblMod).name,
            cut = getattr(process, tblMod).cut,
        ))

        variantTasks.append(cms.Task(
            getattr(process, recoJetMod), getattr(process, corrMod),
            getattr(process, tavMod), getattr(process, chgMod), getattr(process, patMod),
        ))
        variantMCTasks.append(cms.Task(getattr(process, gjmMod), getattr(process, flavMod)))
        variantTableTasks.append(cms.Task(getattr(process, tblMod)))
        variantMCTableTasks.append(cms.Task(getattr(process, mcTblMod)))

    process.scoutingPuppiCalibAK8JetTask = cms.Task(*variantTasks)
    process.scoutingPuppiCalibAK8JetMCTask = cms.Task(*variantMCTasks)
    process.scoutingPuppiCalibAK8PuppiTableTask = cms.Task(*variantTableTasks)
    process.scoutingPuppiCalibAK8PuppiMCTableTask = cms.Task(*variantMCTableTasks)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibAK8JetTask)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibAK8JetMCTask)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibAK8PuppiTableTask)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibAK8PuppiMCTableTask)

    return process


def _addMETTables(process, pName):
    """Wire in PF MET + GenMET NanoAOD tables, for MET/hadronic-recoil
    validation. Neither is present in the stock scouting NanoAOD workflow or
    in scoutingToMiniAODDerivedCollections_cff.py.

    slimmedMETs already exists as a process product from the base scouting->
    MiniAOD customization (PhysicsTools.PatFromScouting.scoutingToMiniAOD_cff,
    Run3ScoutingMETProducer reading hltScoutingPFPacker's raw pfMetPt/pfMetPhi
    -- confirmed the only genuine MET this HLT scouting chain has, there's no
    separate Type-1-corrected version) -- customiseScoutingNano already runs
    before this, so it's just never been exposed as a NanoAOD table.

    GenMET needs more: PhysicsTools.NanoAOD.met_cff's metMCTable reads
    pat::MET::genMET(), which Run3ScoutingMETProducer never sets (confirmed
    0 GenMET branches in existing v5 output) -- and the standard genMetTrue
    recipe (RecoMET/Configuration/python/GenMETParticles_cff.py) needs
    AOD-level genParticles, not available in a MiniAOD-derived workflow like
    this one (only prunedGenParticles/packedGenParticles survive slimming).
    So this rebuilds the same InputGenJetsParticleSelector->GenMETProducer
    chain pointed at packedGenParticles instead (the standard MiniAOD-level
    substitute), then sets it on slimmedMETs via Run3ScoutingMETProducer's
    new opt-in genMET parameter (default off, so every other user of that
    shared plugin is unaffected -- see its own comment).

    A PUPPI-weighted MET is deliberately NOT added here: it's fully
    derivable offline from the existing per-candidate puppiWeight_<variant>
    + kinematics columns (vector-sum of -weight*p_T over all candidates), no
    new production information needed for that piece -- see analysis-side
    code instead.
    """
    from PhysicsTools.NanoAOD.met_cff import pfmetTable, rawMetTable, metMCTable

    if not hasattr(process, "pfMetTable"):
        process.pfMetTable = pfmetTable.clone(src = cms.InputTag("slimmedMETs", "", pName))
        process.rawPFMetTable = rawMetTable.clone(src = process.pfMetTable.src)
        process.scoutingPuppiCalibMETTask = cms.Task(process.pfMetTable, process.rawPFMetTable)
        process.scoutingNanoSequence.associate(process.scoutingPuppiCalibMETTask)

    runOnMC = hasattr(process, "NANOEDMAODSIMoutput") or hasattr(process, "NANOAODSIMoutput")
    if runOnMC and not hasattr(process, "genMetTable"):
        # MiniAOD-level substitute for RecoMET/Configuration/python/
        # GenMETParticles_cff.py's genParticlesForMETAllVisible (same
        # ignoreParticleIDs list -- neutrinos/LSP-like invisibles skipped in
        # the negative-vector-sum), src redirected to packedGenParticles.
        process.genParticlesForMETAllVisible = cms.EDProducer("InputGenJetsParticleSelector",
            src = cms.InputTag("packedGenParticles"),
            partonicFinalState = cms.bool(False),
            excludeResonances = cms.bool(False),
            excludeFromResonancePids = cms.vuint32(),
            tausAsJets = cms.bool(False),
            ignoreParticleIDs = cms.vuint32(
                1000022, 1000012, 1000014, 1000016, 2000012, 2000014, 2000016,
                1000039, 5100039, 4000012, 4000014, 4000016,
                9900012, 9900014, 9900016, 39, 12, 14, 16,
            ),
        )
        process.genMetTrue = cms.EDProducer("GenMETProducer",
            src = cms.InputTag("genParticlesForMETAllVisible"),
            alias = cms.string("genMetTrue"),
            onlyFiducialParticles = cms.bool(False),
            globalThreshold = cms.double(0.0),
            usePt = cms.bool(True),
            applyFiducialThresholdForFractions = cms.bool(False),
        )
        process.slimmedMETs.genMET = cms.InputTag("genMetTrue")

        process.genMetTable = metMCTable.clone(src = process.pfMetTable.src)
        process.scoutingPuppiCalibGenMETTask = cms.Task(
            process.genParticlesForMETAllVisible, process.genMetTrue, process.genMetTable,
        )
        process.scoutingNanoSequence.associate(process.scoutingPuppiCalibGenMETTask)


def _addMuonTable(process, pName):
    """Wire in a minimal muon NanoAOD table, for hadronic-recoil validation
    (reconstructing a Z->mumu candidate to define the recoil axis against --
    see met.py / the MET analysis notes on why QCD alone can't support that
    technique). slimmedMuons already exists as a process product from the
    base scouting->MiniAOD customization (PhysicsTools.PatFromScouting.
    scoutingToMiniAOD_cff, PatFromScoutingMuonProducer) -- same situation
    slimmedMETs was in before this session's MET table addition -- just
    never exposed as a NanoAOD table.

    Deliberately minimal (kinematics + charge + dz/dxy/pfRelIso03, all
    direct pat::Muon accessors needing no extra producers): the standard
    stock muonTable additionally needs slimmedMuonsUpdated/isoForMu/
    linkedObjects (miniIso, jetIdx, svIdx, jet-based isolation) built by the
    full PAT cross-linking chain, which doesn't exist in this workflow and
    isn't needed just to build a dimuon Z candidate.

    Named scoutingPuppiCalibMuonTable / "ScoutingPuppiCalibMuon", NOT
    scoutingMuonTable / "ScoutingMuon" -- the stock scouting NanoAOD
    workflow (PhysicsTools.NanoAOD.custom_run3scouting_cff, loaded before
    this runs) ALREADY defines its own process.scoutingMuonTable, a
    completely different producer (SimpleRun3ScoutingMuonCollectionFlatTable
    Producer, built from the raw HLT Run3ScoutingMuon collection, table
    name="ScoutingMuon") -- confirmed via edmConfigDump after this module's
    first version silently did nothing (its hasattr(process,
    "scoutingMuonTable") guard, meant to prevent double-registration on a
    second call, instead saw the STOCK module and returned immediately
    without ever creating this one). Reusing either the module label or the
    table name would collide with that pre-existing, unrelated table.
    """
    from PhysicsTools.NanoAOD.simplePATMuonFlatTableProducer_cfi import simplePATMuonFlatTableProducer
    from PhysicsTools.NanoAOD.common_cff import CandVars, Var

    if hasattr(process, "scoutingPuppiCalibMuonTable"):
        return
    process.scoutingPuppiCalibMuonTable = simplePATMuonFlatTableProducer.clone(
        src = cms.InputTag("slimmedMuons", "", pName),
        name = cms.string("ScoutingPuppiCalibMuon"),
        doc = cms.string("Scouting muons (slimmedMuons), minimal variable set for dimuon "
                          "(Z candidate / hadronic recoil) reconstruction"),
        variables = cms.PSet(
            CandVars,
            dz = Var("dB('PVDZ')", float, doc="dz (with sign) wrt first PV, in cm", precision=10),
            dxy = Var("dB('PV2D')", float, doc="dxy (with sign) wrt first PV, in cm", precision=10),
            pfRelIso03_all = Var(
                "(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt "
                "+ pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",
                float, doc="PF relative isolation dR=0.3, total (deltaBeta corrections)", precision=8),
        ),
    )
    process.scoutingPuppiCalibMuonTask = cms.Task(process.scoutingPuppiCalibMuonTable)
    process.scoutingNanoSequence.associate(process.scoutingPuppiCalibMuonTask)


def customiseScoutingPuppiCalibrationNano(process, pName="NANO", variants=None, referenceVariant="optimized"):
    """Single entry point for the standalone calibration driver cfg: reuses
    the existing scouting NanoAOD customization unchanged, then adds the
    PUPPI variant collections and calibration tables on top.
    """
    from PhysicsTools.NanoAOD.custom_run3scouting_cff import customiseScoutingNano
    from PhysicsTools.PatFromScouting.scoutingToMiniAODDerivedCollections_cff import customiseScoutingNanoDerived

    process = customiseScoutingNano(process)
    process = customiseScoutingNanoDerived(process, pName)
    process = customiseForScoutingPuppiCalibration(process, pName, variants, referenceVariant)
    process = _addAK8Jets(process, pName, variants, referenceVariant)
    _addMETTables(process, pName)
    _addMuonTable(process, pName)
    return process
