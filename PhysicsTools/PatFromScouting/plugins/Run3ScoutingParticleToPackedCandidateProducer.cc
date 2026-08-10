/**
 * Run3ScoutingParticleToPackedCandidateProducer
 *
 * Converts Run3ScoutingParticle to pat::PackedCandidate for MiniAOD compatibility.
 * Matches charged candidates to reco::Tracks to embed track details
 * (hasTrackDetails() == true, dxyError/dzError/normalizedChi2/hit counts available).
 *
 * Requires vertices and scoutingTracks to be produced first.
 */

#include <memory>
#include <cmath>
#include <limits>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

class Run3ScoutingParticleToPackedCandidateProducer : public edm::stream::EDProducer<> {
public:
  typedef edm::Association<reco::VertexCollection> CandToVertex;

  explicit Run3ScoutingParticleToPackedCandidateProducer(const edm::ParameterSet&);
  ~Run3ScoutingParticleToPackedCandidateProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> particleToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> pdtToken_;
  const bool useCHS_;
  const int covarianceVersion_;
  const int covarianceSchema_;
  const bool useImprovedVertexAssociation_;
  const double maxDzForPrimaryAssignment_;
  const double maxDzSigForPrimaryAssignment_;
  const bool useTrackMatchedVertexAssociation_;
};

Run3ScoutingParticleToPackedCandidateProducer::Run3ScoutingParticleToPackedCandidateProducer(
    const edm::ParameterSet& iConfig)
    : particleToken_(consumes<std::vector<Run3ScoutingParticle>>(iConfig.getParameter<edm::InputTag>("src"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      pdtToken_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>()),
      useCHS_(iConfig.getParameter<bool>("CHS")),
      covarianceVersion_(iConfig.getParameter<int>("covarianceVersion")),
      covarianceSchema_(iConfig.getParameter<int>("covarianceSchema")),
      useImprovedVertexAssociation_(iConfig.getParameter<bool>("useImprovedVertexAssociation")),
      maxDzForPrimaryAssignment_(iConfig.getParameter<double>("maxDzForPrimaryAssignment")),
      maxDzSigForPrimaryAssignment_(iConfig.getParameter<double>("maxDzSigForPrimaryAssignment")),
      useTrackMatchedVertexAssociation_(iConfig.getParameter<bool>("useTrackMatchedVertexAssociation")) {
  produces<reco::PFCandidateCollection>("recoCands");
  produces<pat::PackedCandidateCollection>();
  produces<edm::Association<pat::PackedCandidateCollection>>();
  produces<edm::ValueMap<int>>("quality");
  produces<edm::ValueMap<int>>("vtxass");
  produces<CandToVertex>("vtxass");
}

void Run3ScoutingParticleToPackedCandidateProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto outputReco = std::make_unique<reco::PFCandidateCollection>();
  auto output = std::make_unique<pat::PackedCandidateCollection>();

  const auto& pdt = iSetup.getData(pdtToken_);
  const auto& particles = iEvent.get(particleToken_);

  auto verticesHandle = iEvent.getHandle(vertexToken_);
  const auto& vertices = *verticesHandle;
  reco::VertexRefProd vertexRefProd(verticesHandle);

  edm::Handle<reco::TrackCollection> trackHandle;
  iEvent.getByToken(trackToken_, trackHandle);
  const auto& tracks = *trackHandle.product();
  std::vector<int> mappingTk(tracks.size(), -1);

  // Build a "used" flag so each reco::Track is matched at most once
  std::vector<bool> trackUsed(tracks.size(), false);

  output->reserve(particles.size());
  std::vector<int> mapping;
  mapping.reserve(particles.size());
  std::vector<int> vtx_ass;
  std::vector<int> qual;
  std::vector<int> pfToPVVector;

  for (unsigned int ic = 0, nc = particles.size(); ic < nc; ++ic) {
	  const auto& particle = particles[ic];
  //for (const auto& particle : particles) {
    if (useCHS_ && particle.vertex() > 0) {
      continue;
    }

    int pdgId = particle.pdgId();
    int charge = (pdgId == 22 || pdgId == 130 || pdgId == 1 || pdgId == 2 || pdgId == 0) ? 0 : (pdgId > 0) - (pdgId < 0);

    const HepPDT::ParticleData* pdtData = pdt.particle(HepPDT::ParticleID(particle.pdgId()));
    if (!pdtData) {
      continue;
    }
    float mass = pdtData->mass().value();

    // Index this candidate will occupy in the output collections. Particles can
    // be skipped (CHS pileup drop, unknown PDG id) so this differs from ic.
    const int outIdx = static_cast<int>(output->size());

    reco::PFCandidate::ParticleType particleType = reco::PFCandidate::ParticleType::X;
    switch(std::abs(pdgId)) {
      case 211: // charge hadron
        particleType = reco::PFCandidate::ParticleType::h;
        break;
      case 11: // electron
        particleType = reco::PFCandidate::ParticleType::e;
        break;
      case 13: // muon
        particleType = reco::PFCandidate::ParticleType::mu;
        break;
      case 22: // gamma
        particleType = reco::PFCandidate::ParticleType::gamma;
        break;
      case 130: // neutral hadron
        particleType = reco::PFCandidate::ParticleType::h0;
        break;
      case 1: // HF hadron
        particleType = reco::PFCandidate::ParticleType::h_HF;
        break;
      case 2: // HF em
        particleType = reco::PFCandidate::ParticleType::egamma_HF;
        break;
      case 0:
      default:
        particleType = reco::PFCandidate::ParticleType::X;
        break;
    }



    float pt = particle.pt();
    float eta = particle.eta();
    float phi = particle.phi();
    float px = pt * std::cos(phi);
    float py = pt * std::sin(phi);
    float pz = pt * std::sinh(eta);
    float energy = std::sqrt(px * px + py * py + pz * pz + mass * mass);
    math::XYZTLorentzVector p4(px, py, pz, energy);

    bool relativeTrackVars = particle.relative_trk_vars();
    float trkPt = relativeTrackVars ? particle.trk_pt() + particle.pt() : particle.trk_pt();
    float trkEta = relativeTrackVars ? particle.trk_eta() + particle.eta() : particle.trk_eta();
    float trkPhi = relativeTrackVars ? particle.trk_phi() + particle.phi() : particle.trk_phi();

    // Find a kinematically-matched reco::Track early, before deciding vtxIdx,
    // so a confident match (if any) can be used below for genuine per-vertex
    // dz/dxy via tracks[matchedTrackIdx].dz(Point)/.dxy(Point) -- exact
    // analytic closest-approach formulas (DataFormats/TrackReco/TrackBase.h),
    // depending on the vertex's full (x,y,z) position -- instead of
    // approximating from particle.dz() by a z-only shift. The track-detail
    // embedding block further down reuses this same search result rather
    // than repeating it.
    bool eligibleForTrackMatch =
        (pdgId != 22 && pdgId != 130 && pdgId != 2 && pdgId != 1 && trkPt > 0);
    int matchedTrackIdx = -1;
    float matchedTrackMetric = 999.f;
    if (eligibleForTrackMatch) {
      for (size_t iTk = 0; iTk < tracks.size(); ++iTk) {
        if (trackUsed[iTk])
          continue;
        const auto& tk = tracks[iTk];
        float dEta = trkEta - tk.eta();
        float dPhi = reco::deltaPhi(trkPhi, tk.phi());
        float dR2 = dEta * dEta + dPhi * dPhi;
        float dPtRel = std::abs(trkPt - tk.pt()) / trkPt;
        float metric = dR2 + dPtRel * dPtRel;
        if (metric < matchedTrackMetric) {
          matchedTrackMetric = metric;
          matchedTrackIdx = static_cast<int>(iTk);
        }
      }
    }
    constexpr float kTrackMatchMaxMetric = 0.01f;
    bool hasConfidentTrackMatch = (matchedTrackIdx >= 0 && matchedTrackMetric < kTrackMatchMaxMetric);

    int vtxIdx = particle.vertex();
    float dxy = particle.dxy();
    float dz = particle.dz();

    if (useImprovedVertexAssociation_ && charge != 0 && trkPt > 0 && !vertices.empty()) {
      // particle.vertex() is essentially unusable for charged candidates: at HLT
      // (HLTScoutingPFProducer.cc) it comes from a 100-micron 3D-position match
      // between the PF candidate's *stored* vertex point and the scouting vertex
      // positions. For charged candidates that stored point is the track's helix
      // reference point near the beamline (set in PFAlgo.cc), which essentially
      // never coincides with a reconstructed vertex position -- so the large
      // majority of charged candidates end up with vertex()<0 ("unassociated")
      // regardless of their true origin, not because of any real dz-based
      // incompatibility with the leading vertex. particle.dz()/dzsig(), by
      // contrast, ARE a genuine trk->dz() against vertex 0 (see
      // HLTScoutingPFProducer.cc), so redo a proper nearest-vertex-in-dz search
      // here instead of trusting vertex().
      if (useTrackMatchedVertexAssociation_ && hasConfidentTrackMatch) {
        // tracks[matchedTrackIdx] is a full reco::Track (pat::makeRecoTrack,
        // built by Run3ScoutingTrackToRecoTrackProducer from the same
        // underlying track, with its own reference point/momentum/
        // covariance), so .dz(Point)/.dxy(Point) give a genuine per-vertex
        // impact parameter honoring the vertex's x/y position too, not just
        // its z -- unlike the linear z-shift approximation below, and the
        // same calculation offline PrimaryVertexAssignment-style code
        // performs, just reusing the persisted track covariance for the
        // error rather than a fresh refit.
        const auto& matchedTrack = tracks[matchedTrackIdx];
        double trackDzError = matchedTrack.dzError();
        if (trackDzError > 0.0 && std::isfinite(trackDzError)) {
          double bestDist = std::numeric_limits<double>::max();
          int bestVtx = -1;
          double bestDz = 0.0;
          double bestDzE = 0.0;
          for (size_t iv = 0; iv < vertices.size(); ++iv) {
            double dzI = matchedTrack.dz(vertices[iv].position());
            double dzE = std::hypot(trackDzError, static_cast<double>(vertices[iv].zError()));
            double dist = (dzI / dzE) * (dzI / dzE);
            if (dist < bestDist) {
              bestDist = dist;
              bestVtx = static_cast<int>(iv);
              bestDz = dzI;
              bestDzE = dzE;
            }
          }
          if (bestVtx >= 0 && std::abs(bestDz) < maxDzForPrimaryAssignment_ &&
              std::abs(bestDz) / bestDzE < maxDzSigForPrimaryAssignment_) {
            vtxIdx = bestVtx;
            dz = static_cast<float>(bestDz);
            dxy = static_cast<float>(matchedTrack.dxy(vertices[bestVtx].position()));
          } else {
            vtxIdx = -1;
          }
        } else {
          vtxIdx = -1;
        }
      } else {
        // Fallback used whenever useTrackMatchedVertexAssociation_ is off, or
        // no confident track match was found for this particle: approximate
        // dz to a vertex other than 0 the same way pat::PackedCandidate::
        // dz(ipv) itself approximates dz to a vertex it wasn't packed
        // against -- shifting by the Delta_z between vertex positions --
        // mirroring CommonTools/RecoAlgos/PrimaryVertexAssignment's dz+dzSig
        // window test.
        float dzSig0 = particle.dzsig();
        double dzError0 = (dzSig0 != 0.f) ? std::abs(static_cast<double>(dz) / static_cast<double>(dzSig0)) : 0.0;
        if (dzError0 > 0.0) {
          double bestDist = std::numeric_limits<double>::max();
          int bestVtx = -1;
          double bestDz = 0.0;
          double bestDzE = 0.0;
          for (size_t iv = 0; iv < vertices.size(); ++iv) {
            double dzI = static_cast<double>(dz) - (vertices[iv].position().z() - vertices[0].position().z());
            double dzE = std::hypot(dzError0, static_cast<double>(vertices[iv].zError()));
            double dist = (dzI / dzE) * (dzI / dzE);
            if (dist < bestDist) {
              bestDist = dist;
              bestVtx = static_cast<int>(iv);
              bestDz = dzI;
              bestDzE = dzE;
            }
          }
          if (bestVtx >= 0 && std::abs(bestDz) < maxDzForPrimaryAssignment_ &&
              std::abs(bestDz) / bestDzE < maxDzSigForPrimaryAssignment_) {
            vtxIdx = bestVtx;
            dz = static_cast<float>(bestDz);
          } else {
            vtxIdx = -1;
          }
        }
      }
    }
    pfToPVVector.push_back(vtxIdx);

    reco::VertexRef::key_type pvKey = 0;
    if (vtxIdx >= 0 && static_cast<size_t>(vtxIdx) < vertices.size()) {
      pvKey = static_cast<reco::VertexRef::key_type>(vtxIdx);
    }

    math::XYZPoint pvPos(0, 0, 0);
    if (pvKey < vertices.size()) {
      pvPos = vertices[pvKey].position();
    }

    float sinPhi = std::sin(phi);
    float cosPhi = std::cos(phi);

    math::XYZPoint vtxPos(pvPos.X() - dxy * sinPhi, pvPos.Y() + dxy * cosPhi, pvPos.Z() + dz);

    reco::PFCandidate pfCand(charge, p4, particleType);
    pfCand.setVertex(vtxPos);

    pat::PackedCandidate cand(p4, vtxPos, trkPt, trkEta, trkPhi, particle.pdgId(), vertexRefProd, pvKey);

    // Set lost inner hits
    pat::PackedCandidate::LostInnerHits lostHits = pat::PackedCandidate::noLostInnerHits;
    uint8_t scoutingLostHits = particle.lostInnerHits();
    if (scoutingLostHits == 0) {
      lostHits = pat::PackedCandidate::noLostInnerHits;
    } else if (scoutingLostHits == 1) {
      lostHits = pat::PackedCandidate::oneLostInnerHit;
    } else if (scoutingLostHits >= 2) {
      lostHits = pat::PackedCandidate::moreLostInnerHits;
    }
    cand.setLostInnerHits(lostHits);

    // Set track quality
    int quality = particle.quality();
    bool highPurity = (quality & 4);
    cand.setTrackHighPurity(highPurity);

    // Set PV association quality
    if (vtxIdx == 0) {
      cand.setAssociationQuality(pat::PackedCandidate::UsedInFitTight);
    } else if (vtxIdx > 0) {
      cand.setAssociationQuality(pat::PackedCandidate::CompatibilityDz);
    } else {
      cand.setAssociationQuality(pat::PackedCandidate::NotReconstructedPrimary);
    }


    // Match charged candidates to reco::Tracks and embed track details.
    // Reuses the search already performed above (matchedTrackIdx/
    // matchedTrackMetric/hasConfidentTrackMatch), rather than repeating it.
    if (eligibleForTrackMatch) {
      if (hasConfidentTrackMatch) {
        reco::TrackRef trackRef(trackHandle, matchedTrackIdx);
        pfCand.setTrackRef(trackRef);
        cand.setTrackProperties(tracks[matchedTrackIdx], covarianceSchema_, covarianceVersion_);
        trackUsed[matchedTrackIdx] = true;
      }

      if (pfCand.trackRef().isNonnull() && pfCand.trackRef().id() == trackHandle.id()) {
          mappingTk[pfCand.trackRef().key()] = outIdx;
      }

    }

    mapping.push_back(outIdx);
    outputReco->push_back(pfCand);
    output->push_back(cand);
    qual.push_back(cand.hasTrackDetails() ? cand.pseudoTrack().qualityMask() : (1 << reco::TrackBase::loose));
    vtx_ass.push_back(cand.pvAssociationQuality());
  }

  auto pfHandle = iEvent.put(std::move(outputReco), "recoCands");
  assert(mapping.size() == pfHandle->size());
  auto oh = iEvent.put(std::move(output));
  auto pf2pc = std::make_unique<edm::Association<pat::PackedCandidateCollection>>(oh);
  edm::Association<pat::PackedCandidateCollection>::Filler pf2pcFiller(*pf2pc);
  pf2pcFiller.insert(pfHandle, mapping.begin(), mapping.end());
  pf2pcFiller.insert(trackHandle, mappingTk.begin(), mappingTk.end());

  pf2pcFiller.fill();
  iEvent.put(std::move(pf2pc));

  std::unique_ptr<edm::ValueMap<int>> quality_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_quality(*quality_VM);
  filler_quality.insert(pfHandle, qual.begin(), qual.end());
  filler_quality.fill();
  iEvent.put(std::move(quality_VM), "quality");  

  std::unique_ptr<edm::ValueMap<int>> vtx_ass_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_vtx_ass(*vtx_ass_VM);
  filler_vtx_ass.insert(pfHandle, vtx_ass.begin(), vtx_ass.end());
  filler_vtx_ass.fill();
  iEvent.put(std::move(vtx_ass_VM), "vtxass");  

  std::unique_ptr<CandToVertex> pfCandToOriginalVertexOutput(new CandToVertex(vertexRefProd));
  CandToVertex::Filler cand2VertexFiller(*pfCandToOriginalVertexOutput);
  cand2VertexFiller.insert(pfHandle, pfToPVVector.begin(), pfToPVVector.end());
  cand2VertexFiller.fill();
  iEvent.put(std::move(pfCandToOriginalVertexOutput), "vtxass");
}

void Run3ScoutingParticleToPackedCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltScoutingPFPacker"))
      ->setComment("Input scouting particle collection");
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"))
      ->setComment("Input vertex collection for vertex references");
  desc.add<edm::InputTag>("tracks", edm::InputTag("scoutingTracks"))
      ->setComment("Input reco::Track collection for embedding track details");
  desc.add<bool>("CHS", false)->setComment("Apply Charged Hadron Subtraction (skip vtx > 0)");
  desc.add<int>("covarianceVersion", 1)->setComment("Covariance parameterization version (0=Phase0, 1=Phase1)");
  desc.add<int>("covarianceSchema", 520)->setComment("Covariance packing schema");
  desc.add<bool>("useImprovedVertexAssociation", false)
      ->setComment(
          "Default off, existing consumers unaffected. When true, ignore particle.vertex() for charged "
          "candidates (see comment at its use in produce()) and instead redo a nearest-vertex-in-dz search "
          "from particle.dz()/dzsig(), mirroring CommonTools/RecoAlgos/PrimaryVertexAssignment's dz+dzSig "
          "window test.");
  desc.add<double>("maxDzForPrimaryAssignment", 0.1)
      ->setComment("cm; only used if useImprovedVertexAssociation=True. Matches PrimaryVertexAssignment default.");
  desc.add<double>("maxDzSigForPrimaryAssignment", 5.0)
      ->setComment("only used if useImprovedVertexAssociation=True. Matches PrimaryVertexAssignment default.");
  desc.add<bool>("useTrackMatchedVertexAssociation", false)
      ->setComment(
          "Default off. Only has an effect when useImprovedVertexAssociation=True. When true, for "
          "charged candidates with a confident kinematic match to a 'tracks' entry, use that "
          "reco::Track's own dz(Point)/dxy(Point)/dzError() for the nearest-vertex-in-dz search "
          "(genuine, x/y/z-dependent impact parameter) instead of approximating dz to a vertex "
          "other than 0 with a z-only linear shift. Falls back to the linear-shift approximation "
          "when no confident match exists.");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(Run3ScoutingParticleToPackedCandidateProducer);
