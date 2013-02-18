#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

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
  
//   for (unsigned int i = 0; i < trackList.size(); ++i) {
//     anaVtx.addTrack(trackIndices[i]);
//   }
  
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
  
  // -- fill original tracks for sig tracks
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    TAnaTrack *pTrack = gHFEvent->addSigTrack();
    int mcid          =  gHFEvent->getSimpleTrackMCID(trackIndices[i]);
    TAnaMuon *pM      =  gHFEvent->getSimpleTrackMuon(trackIndices[i]);
    pTrack->fMCID     = mcid;
    pTrack->fMuID     = (pM == 0? 0 : pM->fMuID);
    pTrack->fGenIndex = gHFEvent->getSimpleTrack(trackIndices[i])->getGenIndex(); 
    pTrack->fQ        = trackList[i].charge();
    pTrack->fPlab.SetXYZ(trackList[i].px(),
			 trackList[i].py(),
			 trackList[i].pz()
			 ); 
    pTrack->fIndex    = trackIndices[i];
  }

}


