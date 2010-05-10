#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpMuons.h"

#include <iostream>


#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h" 

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaMuon.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TTrgObj.hh"

#include <TFile.h>
#include <TH1.h>

// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;


// ----------------------------------------------------------------------
HFDumpMuons::HFDumpMuons(const edm::ParameterSet& iConfig):
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fCaloMuonsLabel(iConfig.getUntrackedParameter<InputTag>("calomuonsLabel")),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)) {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpMuons constructor" << endl;
  cout << "---  fMuonsLabel:             " << fMuonsLabel.encode() << endl;
  cout << "---  fCaloMuonsLabel:         " << fCaloMuonsLabel.encode() << endl;
  cout << "---  fDoTruthMatching:        " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS
  cout << "---  fVerbose:                " << fVerbose << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
HFDumpMuons::~HFDumpMuons() {

}


// ----------------------------------------------------------------------
void HFDumpMuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // -- global muons
  Handle<MuonCollection> hMuons;
  if (fVerbose > 0) cout << "==> HFDumpMuons> " << fMuonsLabel << endl;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  
  // -- calo muons
  Handle<CaloMuonCollection> cMuons;
  if (fVerbose > 0) cout << "==> HFDumpMuons> " << fCaloMuonsLabel << endl;
  iEvent.getByLabel(fCaloMuonsLabel, cMuons);
  
  

  int im(0);
  for (reco::MuonCollection::const_iterator iMuon = hMuons->begin(); iMuon != hMuons->end();  iMuon++) {
    fillMuon(*iMuon, im); 
    ++im;
  }

  for (reco::CaloMuonCollection::const_iterator cMuon = cMuons->begin(); cMuon != cMuons->end();  cMuon++) {
    fillCaloMuon(*cMuon, im); 
    ++im;
  }


  if (fVerbose > 0) {
    for (int im = 0; im < gHFEvent->nMuons(); ++im) {
      gHFEvent->getMuon(im)->dump();
    }
  }

}



// ----------------------------------------------------------------------
void HFDumpMuons::fillMuon(const reco::Muon& rm, int im) {

  TAnaMuon *pM = gHFEvent->addMuon();    
  if (rm.innerTrack().isNonnull()) {
    pM->fIndex = rm.innerTrack().index();
  } else {
    pM->fIndex = -23;
  }
  pM->fMuIndex = im;
  pM->fMuID    = muonID(rm);
  pM->fQ       = rm.charge();

  pM->fTimeInOut  = rm.time().timeAtIpInOut; 
  pM->fTimeInOutE = rm.time().timeAtIpInOutErr; 
  pM->fTimeOutIn  = rm.time().timeAtIpOutIn; 
  pM->fTimeOutInE = rm.time().timeAtIpOutInErr; 
  pM->fTimeNdof   = rm.time().nDof;

  TrackRef gTrack = rm.globalTrack();
  TrackRef iTrack = rm.innerTrack();
  TrackRef oTrack = rm.outerTrack();

  if (gTrack.isNonnull()) {
    Track trk(*gTrack);
    pM->fGlobalPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());

    vector<unsigned int> hits = muonStatHits(trk);
    pM->fNhitsDT  = hits.at(0); 
    pM->fNhitsCSC = hits.at(1); 
    pM->fNhitsRPC = hits.at(2); 

  }

  if (iTrack.isNonnull()) {
    Track trk(*iTrack);
    pM->fInnerPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
  }

  if (oTrack.isNonnull()) {
    Track trk(*oTrack);
    pM->fOuterPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
  }

}

// ----------------------------------------------------------------------
void HFDumpMuons::fillCaloMuon(const reco::CaloMuon& rm, int im) {

  TAnaMuon *pM = gHFEvent->addMuon();    
  if (rm.innerTrack().isNonnull()) {
    pM->fIndex = rm.innerTrack().index();
  } else {
    pM->fIndex = -23;
  }
  pM->fMuIndex = im; 
  pM->fMuID   |= 0x1<<15;

  TrackRef iTrack = rm.innerTrack();

  if (iTrack.isNonnull()) {
    Track trk(*iTrack);
    pM->fInnerPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
  }

}



// ----------------------------------------------------------------------
int HFDumpMuons::muonID(const Muon &rm) {
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
  return MuID; 
}

// ----------------------------------------------------------------------
vector<unsigned int> HFDumpMuons::muonStatHits(const reco::Track& tr) {
  vector<unsigned int> theMuonHits;
  unsigned int nRecHitDT(0), nRecHitCSC(0), nRecHitRPC(0);

  for (trackingRecHit_iterator recHit = tr.recHitsBegin(); recHit != tr.recHitsEnd(); ++recHit){
    DetId detIdHit = (*recHit)->geographicalId();
    if (detIdHit.det() == DetId::Muon ){
      if (detIdHit.subdetId() == MuonSubdetId::DT ) nRecHitDT++;
      else if (detIdHit.subdetId() == MuonSubdetId::CSC ) nRecHitCSC++;
      else if (detIdHit.subdetId() == MuonSubdetId::RPC ) nRecHitRPC++;
    }
  }
  
  theMuonHits.push_back(nRecHitDT); 
  theMuonHits.push_back(nRecHitCSC);
  theMuonHits.push_back(nRecHitRPC);
  return theMuonHits; 
}


// ------------ method called once each job just before starting event loop  ------------
void  HFDumpMuons::beginJob() {

  gHFFile->cd();
  // H1D *h1 = new TH1D("h2", "h2", 20, 0., 20.);

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpMuons::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpMuons);
