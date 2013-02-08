#include <iostream>

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" 
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h" 

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpUtilities.hh"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

using namespace std;
using namespace edm;
using namespace reco;

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
void fillSimpleTrack(TSimpleTrack *pTrack, const reco::Track &trackView, 
		     int tidx, int mid, int gidx, const reco::VertexCollection *vc) {

  // -- index and momentum
  pTrack->setIndex(tidx);
  pTrack->setP(trackView.px(), trackView.py(), trackView.pz());
  
  // -- high purity flag
  if (trackView.quality(reco::TrackBase::qualityByName("highPurity"))) {
    pTrack->setHighPurity(1); 
  } else {
    pTrack->setHighPurity(0); 
  }
  
  // -- charge
  pTrack->setCharge(trackView.charge());

  // -- muon index
  pTrack->setMuonID(mid);

  // -- PV index
  if (vc) {
    int pvidx = getPv(tidx, vc); 
    pTrack->setPvIndex(pvidx);
    //    cout << "sT " << tidx << " setting pv index: " << pvidx << "/"; 
    pvidx = pTrack->getPvIndex();
    //    cout << pvidx ;
  } else {
    pTrack->setPvIndex(0xFFF); 
  }


  // -- truth matching with deltaR comparision
  pTrack->setGenIndex(gidx);
  
  //  cout << "/" << pTrack->getPvIndex() << "  (genIdx = " << gidx << ")" << endl;
}


// ----------------------------------------------------------------------
void fillAnaTrack(TAnaTrack *pTrack, const reco::Track &trackView, int tidx, int gidx,
		  const reco::VertexCollection *vc, const reco::MuonCollection *mc, const reco::BeamSpot *bs) {
  
  math::XYZPoint refPt; 
  double xE(0.), yE(0), zE(0.);
  int pvIdx = getPv(tidx, vc);
  //  cout << "aT " << tidx << " setting pv index: " << pvIdx << endl;
  if (pvIdx > -1) {
    refPt = vc->at(pvIdx).position();
    xE = vc->at(pvIdx).xError();
    yE = vc->at(pvIdx).yError();
    zE = vc->at(pvIdx).zError();
  } else {
    if (bs) {
      refPt = bs->position(); 
      xE = bs->BeamWidthX();
      yE = bs->BeamWidthY();
      zE = bs->sigmaZ();
    } else {
      xE = 99.;
      yE = 99.;
      zE = 99.;
    } 
  }
  
  pTrack->fIndex = tidx;
  pTrack->fPlab.SetPtEtaPhi(trackView.pt(),
			    trackView.eta(),
			    trackView.phi()
			    );
  pTrack->fPtE   = trackView.ptError();
  pTrack->fPhiE  = trackView.phiError();
  pTrack->fEtaE  = trackView.etaError();
  
  pTrack->fTip   = trackView.dxy(refPt);
  pTrack->fTipE  = sqrt(trackView.d0Error()*trackView.d0Error() + 0.5*xE*xE + 0.5*yE*yE); 
  
  pTrack->fLip   = trackView.dz(refPt); 
  pTrack->fLipE  = sqrt(trackView.dzError()*trackView.dzError() + zE*zE); 
  pTrack->fPvIdx = pvIdx;
  
  if (bs) {
    math::XYZPoint bsPt(bs->x0(), bs->y0(), bs->z0());
    pTrack->fBsTip  = trackView.dxy(bsPt);
    pTrack->fBsTipE = sqrt(trackView.d0Error()*trackView.d0Error() 
			   + 0.5*bs->BeamWidthX()*bs->BeamWidthX() 
			   + 0.5*bs->BeamWidthY()*bs->BeamWidthY()); 
    pTrack->fBsLip  = trackView.dz(bsPt);
    pTrack->fBsLipE = sqrt(trackView.dzError()*trackView.dzError() + bs->sigmaZ()*bs->sigmaZ());
  } else {
    pTrack->fBsTip  = -99.;
    pTrack->fBsTipE = -99.;		      		      
    pTrack->fBsLip  = -99.;
    pTrack->fBsLipE = -99.;
  }
  pTrack->fdxy  = trackView.dxy();
  pTrack->fdxyE = trackView.dxyError();
  pTrack->fd0   = trackView.d0();
  pTrack->fd0E  = trackView.d0Error();
  pTrack->fdz   = trackView.dz();
  pTrack->fdzE  = trackView.dzError();
  pTrack->fdsz  = trackView.dsz();
  pTrack->fdszE = trackView.dszError();
  
  pTrack->fQ = trackView.charge();
  pTrack->fChi2 = trackView.chi2();
  pTrack->fDof = int(trackView.ndof());
  pTrack->fValidHits = trackView.numberOfValidHits();  
  pTrack->fAlgorithm = trackView.algo(); 
  
  // -- see https://indico.cern.ch/getFile.py/access?contribId=2&resId=1&materialId=slides&confId=123067
  pTrack->fValidHitFraction = trackView.validFraction();
  pTrack->fLayersWithHits   = trackView.hitPattern().trackerLayersWithMeasurement();
  
  
  // -- from: RecoBTag/TrackProbability/src/TrackClassFilter.cc
  reco::TrackBase::TrackQuality trackQualityUndef               =  reco::TrackBase::qualityByName("undefQuality");
  reco::TrackBase::TrackQuality trackQualityLoose               =  reco::TrackBase::qualityByName("loose");
  reco::TrackBase::TrackQuality trackQualityTight               =  reco::TrackBase::qualityByName("tight");
  reco::TrackBase::TrackQuality trackQualityhighPur             =  reco::TrackBase::qualityByName("highPurity");
  reco::TrackBase::TrackQuality trackQualityConfirmed           =  reco::TrackBase::qualityByName("confirmed");
  reco::TrackBase::TrackQuality trackQualityGoodIterative       =  reco::TrackBase::qualityByName("goodIterative");
  reco::TrackBase::TrackQuality trackQualityLooseSetWithPV      =  reco::TrackBase::qualityByName("looseSetWithPV");
  reco::TrackBase::TrackQuality trackQualityHighPuritySetWithPV =  reco::TrackBase::qualityByName("highPuritySetWithPV");
  
  int trakQuality  = 0;
  if (trackView.quality(trackQualityUndef))               trakQuality |= 0x1<<10;
  if (trackView.quality(trackQualityLoose))               trakQuality |= 0x1<<0;
  if (trackView.quality(trackQualityTight))               trakQuality |= 0x1<<1;
  if (trackView.quality(trackQualityhighPur))             trakQuality |= 0x1<<2;
  if (trackView.quality(trackQualityConfirmed))           trakQuality |= 0x1<<3;
  if (trackView.quality(trackQualityGoodIterative))       trakQuality |= 0x1<<4;
  if (trackView.quality(trackQualityLooseSetWithPV))      trakQuality |= 0x1<<5;
  if (trackView.quality(trackQualityHighPuritySetWithPV)) trakQuality |= 0x1<<6;
  pTrack->fTrackQuality = trakQuality; 
  
  // -- Muon ID
  pTrack->fMuIndex = -1; 
  pTrack->fMuID    = 0; 
  if (mc) {
    for (unsigned int i = 0; i < mc->size(); ++i) {
      const reco::Muon muon = mc->at(i); 
      TrackRef track = muon.innerTrack();
      if (static_cast<unsigned int>(tidx) == track.index()) {
	pTrack->fMuIndex = i; 
	pTrack->fMuID    = muonID(muon);
	break;
      }
    }
  }
  
  // -- hits of the track
  const reco::HitPattern& p = trackView.hitPattern();
  for (int i=0; i<p.numberOfHits(); i++) {
    uint32_t hit = p.getHitPattern(i);
    if (i < 20) pTrack->fHitPattern[i] = hit; 
  }
  
  // -- truth-matching information
  pTrack->fGenIndex = gidx;
  pTrack->fMCID     = ((gidx > -1 && gidx < gHFEvent->nGenCands()) ? gHFEvent->getGenCand(gidx)->fID : 0);

}


// ----------------------------------------------------------------------
int getPv(int tidx, const reco::VertexCollection *vc) {
  if (vc) {
    for (unsigned int i = 0; i < vc->size(); ++i) {
      Vertex::trackRef_iterator v1TrackIter;
      Vertex::trackRef_iterator v1TrackBegin = vc->at(i).tracks_begin();
      Vertex::trackRef_iterator v1TrackEnd   = vc->at(i).tracks_end();
      for (v1TrackIter = v1TrackBegin; v1TrackIter != v1TrackEnd; v1TrackIter++) {
	if (static_cast<unsigned int>(tidx) == v1TrackIter->key()) return i;
      }
    }
  }
  return -1; 
}


// ----------------------------------------------------------------------
int muonID(const Muon &rm) {
  int MuID(0); 
  if (muon::isGoodMuon(rm, muon::AllStandAloneMuons))               MuID |= 0x1<<0; 
  if (muon::isGoodMuon(rm, muon::AllGlobalMuons))                   MuID |= 0x1<<1; 
  if (muon::isGoodMuon(rm, muon::AllTrackerMuons))                  MuID |= 0x1<<2; 
  if (muon::isGoodMuon(rm, muon::TrackerMuonArbitrated))            MuID |= 0x1<<4; 
  if (muon::isGoodMuon(rm, muon::GlobalMuonPromptTight))            MuID |= 0x1<<6; 
  if (muon::isGoodMuon(rm, muon::TMLastStationLoose))               MuID |= 0x1<<7; 
  if (muon::isGoodMuon(rm, muon::TMLastStationTight))               MuID |= 0x1<<8; 
  if (muon::isGoodMuon(rm, muon::TM2DCompatibilityLoose))           MuID |= 0x1<<9; 
  if (muon::isGoodMuon(rm, muon::TM2DCompatibilityTight))           MuID |= 0x1<<10; 
  if (muon::isGoodMuon(rm, muon::TMOneStationLoose))                MuID |= 0x1<<11; 
  if (muon::isGoodMuon(rm, muon::TMOneStationTight))                MuID |= 0x1<<12; 
  if (muon::isGoodMuon(rm, muon::TMLastStationOptimizedLowPtLoose)) MuID |= 0x1<<13; 
  if (muon::isGoodMuon(rm, muon::TMLastStationOptimizedLowPtTight)) MuID |= 0x1<<14;
  //MuID |= 0x1<<15; // used in sigtrack for tightmu selection w.r.t. selected PV
  //MuID |= 0x1<<16; // used in sigtrack for tightmu selection w.r.t. SV
  return MuID; 
}
