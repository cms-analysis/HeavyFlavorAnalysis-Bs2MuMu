#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"

#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

using reco::Vertex;
using reco::Track;
using std::vector;
using std::cout;
using std::endl;

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFKalmanVertexFit::HFKalmanVertexFit(const TransientTrackBuilder *TTB, 
				     Vertex  &PV,
				     int type,
				     int verbose) {
  
  fpTTB                 = TTB;
  fPV                   = PV;
  fType                 = type; 
  fVerbose              = verbose;
  
  setNoCuts();
}


// ----------------------------------------------------------------------
HFKalmanVertexFit::~HFKalmanVertexFit() {
}


// ----------------------------------------------------------------------
void HFKalmanVertexFit::setNoCuts() {
  fMaxDoca  = 9999.;
  fVtxChi2  = 9999.;
  fVtxSigXY = -99.;
  fVtxSig3d = -99.;
  fCosAngle = -99.;
  fPtCand   = -99.;
}

// ----------------------------------------------------------------------
void HFKalmanVertexFit::doFit(vector<Track>  &trackList, 
			      vector<int>    &trackIndices,
			      vector<double> &trackMasses, 
			      int type, 
			      int ntracks
			      ) {
  
  if (fVerbose > 5) cout << "==>HFKalmanVertexFit> doFit()" << endl;

  vector<reco::TransientTrack> RecoTransientTrack;
  RecoTransientTrack.clear();
  
  if (ntracks < 0) {
    for (unsigned int i = 0; i < trackList.size(); ++i) {
      RecoTransientTrack.push_back(fpTTB->build(trackList[i]));
    }
  } else {
    for (int i = 0; i < ntracks; ++i) {
      RecoTransientTrack.push_back(fpTTB->build(trackList[i]));
    }
  }


  // -- Compute minimal and maximal closest approach for any combination of tracks
  TwoTrackMinimumDistance md;
  double dist(99.); 
  double minDist(99.), mini, minj; 
  double maxDist(-99.), maxi, maxj; 
  for (unsigned int i = 0; i < RecoTransientTrack.size(); ++i) {
    for (unsigned int j = i+1; j < RecoTransientTrack.size(); ++j) {
      md.calculate(RecoTransientTrack[i].initialFreeState(), RecoTransientTrack[j].initialFreeState());
      if (md.status()) {
	dist = md.distance();
	if (dist < minDist) {
	  mini = i; 
	  minj = j;
	  minDist = dist;
	}

	if (dist > maxDist) {
	  maxi = i; 
	  maxj = j;
	  maxDist = dist;
	}
      }
    }
  }

  if (maxDist > fMaxDoca) {
    if (fVerbose > 1) cout << "maxDist = " << maxDist << " > " << fMaxDoca << ", not fitting, return ..." << endl;
    return;
  }

  // -- Do the vertexing
  KalmanVertexFitter theFitter(true);
  TransientVertex TransSecVtx = theFitter.vertex(RecoTransientTrack); 
  if (TransSecVtx.isValid() ) {
    if (isnan(TransSecVtx.position().x()) 
	|| isnan(TransSecVtx.position().y()) 
	|| isnan(TransSecVtx.position().z()) ) {
      if (fVerbose > 0) cout << "==>doKalmanVertexFit> Something went wrong! SecVtx nan - continue ... " << endl;
      return; 
    }
  } else {
    if (fVerbose > 0) {
      cout << "==>doKalmanVertexFit> KVF " << type << " failed for tracks:";
      for (unsigned int i = 0; i < trackIndices.size(); ++i) cout << " " << trackIndices[i];
      cout << ", continue ..." << endl;
    }
    return; 
  }

  if (TransSecVtx.totalChiSquared() < 0) return; 
  
  // -- Get refitted tracks
  vector<reco::TransientTrack> refTT = TransSecVtx.refittedTracks();
  vector<Track> refT; refT.clear(); 
  for(vector<reco::TransientTrack>::const_iterator i = refTT.begin(); i != refTT.end(); i++) {
    const Track &ftt = i->track();
    refT.push_back(ftt);
  }
  
  // -- Build composite
  TLorentzVector comp, m1, comp2;
  if (ntracks < 0) {
    for (unsigned int i = 0; i < trackList.size(); ++i) {
      m1.SetXYZM(refT[i].px(), refT[i].py(), refT[i].pz(), TMath::Abs(trackMasses[i])); 
      comp += m1; 
      if (trackMasses[i] > 0) comp2 += m1; 
    }
  } else {
    for (int i = 0; i < ntracks; ++i) {
      m1.SetXYZM(refT[i].px(), refT[i].py(), refT[i].pz(), TMath::Abs(trackMasses[i])); 
      comp += m1; 
      if (trackMasses[i] > 0) comp2 += m1; 
    }
    for (unsigned int i = ntracks; i < trackList.size(); ++i) {
      m1.SetXYZM(trackList[i].px(), trackList[i].py(), trackList[i].pz(), TMath::Abs(trackMasses[i])); 
      comp += m1; 
      if (trackMasses[i] > 0) comp2 += m1; 
    }

  }


//   // -- Build a kinematic particle to determine the best PV (?)
//   float sigma = 0.00001*mass;
//   kinParticles.push_back(pFactory.particle(refT[0],trackMasses[0],0.0f,0.0f,sigma)); 
//   kinParticles.push_back(pFactory.particle(refT[1],trackMasses[1],0.0f,0.0f,sigma)); 
		
//   md.calculate(RecoTransientTrack[i].initialFreeState(), RecoTransientTrack[j].initialFreeState());
//   if (md.status()) {
//     dist = md.distance();
//   }


  // -- Build "rootio" vertex 
  TAnaVertex anaVtx;
  ChiSquared chi(TransSecVtx.totalChiSquared(), TransSecVtx.degreesOfFreedom());
  anaVtx.setInfo(chi.value(), chi.degreesOfFreedom(), chi.probability(), 0, 0);
  anaVtx.fPoint.SetXYZ(TransSecVtx.position().x(), 
		       TransSecVtx.position().y(), 
		       TransSecVtx.position().z());
  
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    anaVtx.addTrack(trackIndices[i]);
  }

  VertexDistanceXY axy;
  anaVtx.fDxy     = axy.distance(fPV, TransSecVtx).value();
  anaVtx.fDxyE    = axy.distance(fPV, TransSecVtx).error();

  VertexDistance3D a3d;
  anaVtx.fD3d     = a3d.distance(fPV, TransSecVtx).value();
  anaVtx.fD3dE    = a3d.distance(fPV, TransSecVtx).error();


  // -- Get covariance matrix
  
  double cov[9];
  cov[0] = TransSecVtx.positionError().cxx();
  cov[1] = TransSecVtx.positionError().cyx();
  cov[2] = TransSecVtx.positionError().czx();
  cov[3] = TransSecVtx.positionError().cyx();
  cov[4] = TransSecVtx.positionError().cyy();
  cov[5] = TransSecVtx.positionError().czy();
  cov[6] = TransSecVtx.positionError().czx();
  cov[7] = TransSecVtx.positionError().czy();
  cov[8] = TransSecVtx.positionError().czz();
  anaVtx.setCovXX(cov);

  // -- Apply cuts
  if (comp.Vect().Perp() < fPtCand) {
    if (fVerbose > 1) cout << "pT(Cand) = " << comp.Vect().Perp() << " < " << fPtCand << ", not filling, return ..." << endl;
    return;
  }

  if (anaVtx.fDxy/anaVtx.fDxyE < fVtxSigXY) {
    if (fVerbose > 1) cout << "axyV/axyE = " << anaVtx.fDxy/anaVtx.fDxyE  << " < " << fVtxSigXY << ", not filling, return ..." << endl;
    return;
  }

  if (anaVtx.fD3d/anaVtx.fD3dE < fVtxSig3d) {
    if (fVerbose > 1) cout << "a3dV/a3dE = " << anaVtx.fD3d/anaVtx.fD3dE  << " < " << fVtxSig3d << ", not filling, return ..." << endl;
    return;
  }

  if (anaVtx.fChi2 > fVtxChi2) {
    if (fVerbose > 1) cout << "vertexChi2 = " << anaVtx.fChi2  << " < " << fVtxChi2 << ", not filling, return ..." << endl;
    return;
  }

  // -- Compute 2d (??) pointing angle
  TVector2 tPV = TVector2(fPV.position().x(), fPV.position().y());
  TVector2 tSV = TVector2(anaVtx.fPoint.X(), anaVtx.fPoint.Y());
  TVector2 tl  = tSV - tPV;
  TVector2 tB  = comp.Vect().XYvector();
  double CosAngle = TMath::Cos(tB.DeltaPhi(tl));

  if (CosAngle < fCosAngle) {
    if (fVerbose > 1) cout << "CosAngle = " << CosAngle  << " < " << fCosAngle << ", not filling, return ..." << endl;
    return;
  }

 
  // -- fill candidate
  if (fVerbose > 1) cout << "Filling this one" << std::endl;
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fPlab    = comp.Vect();
  pCand->fMass    = comp.M();
  pCand->fVtx     = anaVtx;    
  pCand->fType    = (type != 0? type: fType);
  pCand->fDau1    = -1;
  pCand->fDau2    = -1;
  pCand->fSig1    = gHFEvent->nSigTracks();
  pCand->fSig2    = pCand->fSig1 + trackList.size() - 1;

  pCand->fMinDoca = minDist;
  pCand->fMaxDoca = maxDist;
  pCand->fVar1    = comp2.M();
  
  // -- fill refitted sig tracks
  Track trk; 
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    TAnaTrack *pTrack = gHFEvent->addSigTrack();
    //  pTrack->fMCID     = trackList[i].charge()*13;  //??? FIXME ???
    pTrack->fMCID     = gHFEvent->getRecTrack(trackIndices[i])->fMCID;
    pTrack->fMuID     = gHFEvent->getRecTrack(trackIndices[i])->fMuID;    
    pTrack->fGenIndex = gHFEvent->getRecTrack(trackIndices[i])->fGenIndex; 
    pTrack->fQ        = trackList[i].charge();

    if (ntracks < 0) {
      trk = refT[i];
    } else {
      if (i < static_cast<unsigned int>(ntracks)) {
	trk = refT[i];
      } else {
	trk = trackList[i]; 
      }
    }

    pTrack->fPlab.SetXYZ(trk.px(),
			 trk.py(),
			 trk.pz()
			 ); 
    pTrack->fIndex    = trackIndices[i];
  }

}


// ----------------------------------------------------------------------
void HFKalmanVertexFit::doNotFit(vector<Track>  &trackList, 
				 vector<int>    &trackIndices,
				 vector<double> &trackMasses, 
				 int type
				 ) {

  // -- Build composite
  TLorentzVector comp, m1, comp2;
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    m1.SetXYZM(trackList[i].px(), trackList[i].py(), trackList[i].pz(), TMath::Abs(trackMasses[i])); 
    comp += m1; 
    if (trackMasses[i] > 0) comp2 += m1; 
  }
    
  // -- Build dummy vertex for ntuple
  TAnaVertex anaVtx;
  anaVtx.setInfo(-99., -99., -99., -99, -99);
  anaVtx.fPoint.SetXYZ(-9999., -9999., -9999.);
  
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    anaVtx.addTrack(trackIndices[i]);
  }
  
  anaVtx.fDxy     = -99.;
  anaVtx.fDxyE    = -99.;
  anaVtx.fD3d     = -99.;
  anaVtx.fD3dE    = -99.;

  // -- Build transient tracks so that the distance can be calculated
  double dist(99.); 
  double minDist(99.), mini, minj; 
  double maxDist(-99.), maxi, maxj; 
  if (fpTTB) {
    vector<reco::TransientTrack> RecoTransientTrack;
    RecoTransientTrack.clear();
    
    for (unsigned int i = 0; i < trackList.size(); ++i) {
      RecoTransientTrack.push_back(fpTTB->build(trackList[i]));
    }
    
    // -- Compute minimal and maximal closest approach for any combination of tracks
    TwoTrackMinimumDistance md;
    for (unsigned int i = 0; i < RecoTransientTrack.size(); ++i) {
      for (unsigned int j = i+1; j < RecoTransientTrack.size(); ++j) {
	md.calculate(RecoTransientTrack[i].initialFreeState(), RecoTransientTrack[j].initialFreeState());
	if (md.status()) {
	  dist = md.distance();
	  if (dist < minDist) {
	    mini = i; 
	    minj = j;
	    minDist = dist;
	  }
	  
	  if (dist > maxDist) {
	    maxi = i; 
	    maxj = j;
	    maxDist = dist;
	  }
	}
      }
    }
  } else {
    minDist = 99.;
    maxDist = 99.;
  }
        
  // -- fill candidate
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fPlab    = comp.Vect();
  pCand->fMass    = comp.M();
  pCand->fVtx     = anaVtx;    
  pCand->fType    = (type != 0? type: fType);
  pCand->fDau1    = -1;
  pCand->fDau2    = -1;
  pCand->fSig1    = gHFEvent->nSigTracks();
  pCand->fSig2    = pCand->fSig1 + trackList.size() - 1;

  pCand->fMinDoca = minDist;
  pCand->fMaxDoca = maxDist;
  pCand->fVar1    = comp2.M();
  
  // -- fill original tracks for sig tracks
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    TAnaTrack *pTrack = gHFEvent->addSigTrack();
    //  pTrack->fMCID     = trackList[i].charge()*13;  //??? FIXME ???
    pTrack->fMCID     = gHFEvent->getRecTrack(trackIndices[i])->fMCID;
    pTrack->fMuID     = gHFEvent->getRecTrack(trackIndices[i])->fMuID;
    pTrack->fGenIndex = gHFEvent->getRecTrack(trackIndices[i])->fGenIndex; 
    pTrack->fQ        = trackList[i].charge();
    pTrack->fPlab.SetXYZ(trackList[i].px(),
			 trackList[i].py(),
			 trackList[i].pz()
			 ); 
    pTrack->fIndex    = trackIndices[i];
  }

}


