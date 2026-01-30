#include "RecoVertex/AdaptiveVertexFinder/plugins/InclusiveVertexFinder.h"
#include "RecoVertex/AdaptiveVertexFinder/plugins/InclusiveVertexFinderScouting.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/View.h"

typedef TemplatedInclusiveVertexFinder<reco::TrackCollection, reco::Vertex> InclusiveVertexFinder;
typedef TemplatedInclusiveVertexFinder<edm::View<reco::Candidate>, reco::VertexCompositePtrCandidate>
    InclusiveCandidateVertexFinder;
typedef TemplatedInclusiveVertexFinderScouting<reco::TrackCollection, reco::Vertex> InclusiveVertexFinderScouting;
typedef TemplatedInclusiveVertexFinderScouting<edm::View<reco::Candidate>, reco::VertexCompositePtrCandidate>
    InclusiveCandidateVertexFinderScouting;

DEFINE_FWK_MODULE(InclusiveVertexFinder);
DEFINE_FWK_MODULE(InclusiveCandidateVertexFinder);
DEFINE_FWK_MODULE(InclusiveVertexFinderScouting);
DEFINE_FWK_MODULE(InclusiveCandidateVertexFinderScouting);
