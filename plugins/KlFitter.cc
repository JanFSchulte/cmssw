#include "KlFitter.h"

KlFitter::KlFitter(std::vector<reco::TransientTrack> &vecttrk) {
  KalmanVertexFitter vtxFitter(true);
  dimuvtx = vtxFitter.vertex(vecttrk);
  if (!dimuvtx.isValid())
    status_ = false;
  else {
    refited = dimuvtx.refittedTracks();
    prob_ = ChiSquaredProbability(dimuvtx.totalChiSquared(), dimuvtx.degreesOfFreedom());
    normalchi2_ = dimuvtx.totalChiSquared() / dimuvtx.degreesOfFreedom();
  }
}

KlFitter::~KlFitter(){};

void KlFitter::fillNtuple(NtupleContent &nt, bool isTuneP) {
  if (isTuneP) {
    if (status_) {
      nt.branches["pair_tuneP_svprob"] = (float)prob_;
      nt.branches["pair_tuneP_fit_mass"] = (float)DimuonMass(refited[0].track().pt(),
                                          refited[0].track().eta(),
                                          refited[0].track().phi(),
                                          refited[1].track().pt(),
                                          refited[1].track().eta(),
                                          refited[1].track().phi());
      nt.branches["pair_tuneP_normalchi2"] = (float)normalchi2_;
    } else {
      nt.branches["pair_tuneP_svprob"] = (float)-1;
      nt.branches["pair_tuneP_fit_mass"] = (float)-1;
      nt.branches["pair_tuneP_normalchi2"] = (float)-1;
    }
  } else {
    if (status_) {
      nt.branches["pair_svprob"] = (float)prob_;
      nt.branches["pair_fit_mass"] = (float)DimuonMass(refited[0].track().pt(),
                                    refited[0].track().eta(),
                                    refited[0].track().phi(),
                                    refited[1].track().pt(),
                                    refited[1].track().eta(),
                                    refited[1].track().phi());
      nt.branches["pair_normalchi2"] = (float)normalchi2_;
    } else {
      nt.branches["pair_svprob"] = (float)-1;
      nt.branches["pair_fit_mass"] = (float)-1;
      nt.branches["pair_normalchi2"] = (float)-1;
    }
  }
}
