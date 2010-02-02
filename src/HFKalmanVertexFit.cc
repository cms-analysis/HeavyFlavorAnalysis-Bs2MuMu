#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFKalmanVertexFit::HFKalmanVertexFit(const TransientTrackBuilder *TTB, 
				     Vertex  &PV,
				     int type,
				     int verbose) {
  
  fpTTB     = TTB;
  fPV       = PV;
  fType     = type; 
  fVerbose  = verbose;
}

// ----------------------------------------------------------------------
HFKalmanVertexFit::~HFKalmanVertexFit() {
}


// ----------------------------------------------------------------------
void HFKalmanVertexFit::doFit(vector<Track>  &trackList, 
			      vector<int>    &trackIndices,
			      vector<double> &trackMasses, 
			      int type
			      ) {
  
  vector<reco::TransientTrack> RecoTransientTrack;
  RecoTransientTrack.clear();
  
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    RecoTransientTrack.push_back(fpTTB->build(trackList[i]));
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
    if (fVerbose > 0) cout << "==>doKalmanVertexFit> KVF failed! continue ..." << endl;
    return; 
  }
  
  // -- Get refitted tracks
  vector<TransientTrack> refTT = TransSecVtx.refittedTracks();
  vector<Track> refT; refT.clear(); 
  for(vector<TransientTrack>::const_iterator i = refTT.begin(); i != refTT.end(); i++) {
    const Track &ftt = i->track();
    refT.push_back(ftt);
  }
  
  // -- Build composite
  TLorentzVector comp, m1;
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    m1.SetXYZM(refT[i].px(), refT[i].py(), refT[i].pz(), trackMasses[i]); 
    comp += m1; 
  }
  
  
  // -- Build vertex for ntuple
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
        
  // -- fill candidate
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fPlab = comp.Vect();
  pCand->fMass = comp.M();
  pCand->fVtx  = anaVtx;    
  pCand->fType = (type != 0? type: fType);
  pCand->fDau1 = -1;
  pCand->fDau2 = -1;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + trackList.size() - 1;
  
  // -- fill refitted sig tracks
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    TAnaTrack *pTrack = gHFEvent->addSigTrack();
    pTrack->fMCID     = trackList[i].charge()*13;  //??? FIXME ???
    pTrack->fGenIndex = gHFEvent->getRecTrack(trackIndices[i])->fGenIndex; 
    pTrack->fQ        = trackList[i].charge();
    pTrack->fPlab.SetXYZ(refT[i].px(),
			 refT[i].py(),
			 refT[i].pz()
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
  TLorentzVector comp, m1;
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    m1.SetXYZM(trackList[i].px(), trackList[i].py(), trackList[i].pz(), trackMasses[i]); 
    comp += m1; 
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
        
  // -- fill candidate
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fPlab = comp.Vect();
  pCand->fMass = comp.M();
  pCand->fVtx  = anaVtx;    
  pCand->fType = (type != 0? type: fType);
  pCand->fDau1 = -1;
  pCand->fDau2 = -1;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + trackList.size() - 1;
  
  // -- fill original tracks for sig tracks
  for (unsigned int i = 0; i < trackList.size(); ++i) {
    TAnaTrack *pTrack = gHFEvent->addSigTrack();
    pTrack->fMCID     = trackList[i].charge()*13;  //??? FIXME ???
    pTrack->fGenIndex = gHFEvent->getRecTrack(trackIndices[i])->fGenIndex; 
    pTrack->fQ        = trackList[i].charge();
    pTrack->fPlab.SetXYZ(trackList[i].px(),
			 trackList[i].py(),
			 trackList[i].pz()
			 ); 
    pTrack->fIndex    = trackIndices[i];
  }

}


