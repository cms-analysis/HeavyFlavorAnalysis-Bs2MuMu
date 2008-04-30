#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpMuons.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"


// -- Yikes!
extern TAna00Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace reco;
using namespace edm;

// ----------------------------------------------------------------------
BmmDumpMuons::BmmDumpMuons(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fMuonsLabel1(iConfig.getUntrackedParameter<InputTag>("muonsLabel1")),
  fMuonsLabel2(iConfig.getUntrackedParameter<InputTag>("muonsLabel2")),
  fTracksLabel1(iConfig.getUntrackedParameter<InputTag>("tracksLabel1")),
  fTracksLabel2(iConfig.getUntrackedParameter<InputTag>("tracksLabel2")),
  fTracksLabel3(iConfig.getUntrackedParameter<InputTag>("tracksLabel3")),
  fTracksLabel4(iConfig.getUntrackedParameter<InputTag>("tracksLabel4")),
  fTracksLabel5(iConfig.getUntrackedParameter<InputTag>("tracksLabel5")),
  fTracksLabel6(iConfig.getUntrackedParameter<InputTag>("tracksLabel6")),
  fTracksLabel7(iConfig.getUntrackedParameter<InputTag>("tracksLabel7"))

{
  using namespace std;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- BmmDumpMuons constructor" << endl;
  cout << "--- Verbose         : " << fVerbose << endl;
  cout << "--- muonsLabel1      : " << fMuonsLabel1 << endl;
  cout << "--- muonsLabel2      : " << fMuonsLabel2 << endl;
  cout << "--- tracksLabel1     : " << fTracksLabel1 << endl;
  cout << "--- tracksLabel2     : " << fTracksLabel2 << endl;
  cout << "--- tracksLabel3     : " << fTracksLabel3 << endl;
  cout << "--- tracksLabel4     : " << fTracksLabel4 << endl;
  cout << "--- tracksLabel5     : " << fTracksLabel5 << endl;
  cout << "--- tracksLabel6     : " << fTracksLabel6 << endl;
  cout << "--- tracksLabel7     : " << fTracksLabel7 << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fNevt = 0;

}


// ----------------------------------------------------------------------
BmmDumpMuons::~BmmDumpMuons() {
  
}


// ----------------------------------------------------------------------
void BmmDumpMuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  fNevt++;

  using reco::TrackCollection;
  using reco::MuonCollection;

  TAnaTrack *pTrack;

  // -- Look at muons
  const reco::Track* atg = 0;
  const reco::Track* ttg = 0;
  const reco::Track* ctg = 0;

  int mcnt(0), type(0), idrec(0);

  // ===============================================================================================
  // -- Muon Collection: Global Muons
  // ===============================================================================================

  Handle<reco::MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel1, hMuons);

  if (fVerbose > 0) cout << "==>BmmDumpMuons> nGlobalMuons = " << hMuons->size() << ", event: " << fNevt << endl;

  for (reco::MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
  
    mcnt++;

    idrec = (muon->track()).index();

    // -- standalone muon
    type = 31;
    TrackRef stagTrack = muon->standAloneMuon();
    atg = &(*stagTrack);

    pTrack   = gHFEvent->addSigTrack(); 
    pTrack->fMuType   = type;
    pTrack->fMuID     = mcnt;
    pTrack->fIndex    = idrec;
    pTrack->fGenIndex = -1; 
    pTrack->fMCID     = atg->charge()*-13; 
    pTrack->fQ        = atg->charge();

    pTrack->fPlab.SetPtEtaPhi(atg->pt(),
			      atg->eta(),
			      atg->phi()
			      );

    pTrack->fTip  = atg->d0();
    pTrack->fLip  = atg->dz();
    pTrack->fChi2 = atg->chi2();
    pTrack->fDof  = int(atg->ndof());
    pTrack->fHits = atg->numberOfValidHits();  


    if (fVerbose > 0) {

      cout << "--------------------------------------------------------------------" << endl; 
      cout << "-- Muon " <<fMuonsLabel1 << ", #" << mcnt <<  ",  track index (" 
	   << (muon->track()).index() << ")" << endl;
  
      cout << "%%> Type-" << type << ": " ; 
      pTrack->dump(); 
    }
    
    // -- tracker muon
    type = 32;
    TrackRef trkgTrack = muon->track();
    ttg = &(*trkgTrack);

    pTrack   = gHFEvent->addSigTrack(); 
    pTrack->fMuType   = type;
    pTrack->fMuID     = mcnt;
    pTrack->fIndex    = idrec;
    pTrack->fGenIndex = -1; 
    pTrack->fMCID     = ttg->charge()*-13; 
    pTrack->fQ        = ttg->charge();

    pTrack->fPlab.SetPtEtaPhi(ttg->pt(),
			      ttg->eta(),
			      ttg->phi()
			      );

    pTrack->fTip  = ttg->d0();
    pTrack->fLip  = ttg->dz();
    pTrack->fChi2 = ttg->chi2();
    pTrack->fDof  = int(ttg->ndof());
    pTrack->fHits = ttg->numberOfValidHits();  

    if (fVerbose > 0) {
      cout << "%%> Type-" << type <<  ": "; 
      pTrack->dump(); 
    }

    // -- combined muon
    type = 33;
    TrackRef comgTrack = muon->combinedMuon();
    ctg = &(*comgTrack);

    pTrack   = gHFEvent->addSigTrack(); 
    pTrack->fMuType   = type;
    pTrack->fMuID     = mcnt;
    pTrack->fIndex    = idrec;
    pTrack->fGenIndex = -1; 
    pTrack->fMCID     = ctg->charge()*-13; 
    pTrack->fQ        = ctg->charge();

    pTrack->fPlab.SetPtEtaPhi(ctg->pt(),
			      ctg->eta(),
			      ctg->phi()
			      );

    pTrack->fTip  = ctg->d0();
    pTrack->fLip  = ctg->dz();
    pTrack->fChi2 = ctg->chi2();
    pTrack->fDof  = int(ctg->ndof());
    pTrack->fHits = ctg->numberOfValidHits();  

    if (fVerbose > 0) {
      cout << "%%> Type-" << type << ": "; 
      pTrack->dump(); 
    }

    if (fVerbose > 0) {
      cout << "--------------------------------------------------------------------" << endl;; 
    }
  }

  if  (fNevt < 10000) {  
    
    TH2D *h2 = (TH2D*)gHFFile->Get("nmuons");
    h2->SetBinContent(fNevt - 100*int(fNevt/100), int(fNevt/100), mcnt);
  }  



  // -- Look at muons
  const reco::Track* atm = 0;
  const reco::Track* ttm = 0;
  const reco::Track* ctm = 0;

  int tcnt(0);

  // ===============================================================================================
  // -- Muon Collection: Tracker Muons
  // ===============================================================================================

  Handle<reco::MuonCollection> tkMuons;
  iEvent.getByLabel(fMuonsLabel2, tkMuons);

  if (fVerbose > 0) cout << "==>BmmDumpMuons> nTrackerMuons = " << tkMuons->size() << ", event: " << fNevt << endl;

  for (reco::MuonCollection::const_iterator tkmu = tkMuons->begin(); tkmu != tkMuons->end(); ++tkmu) {
    
    tcnt++;

    idrec = (tkmu->track()).index();

    // -- standalone muon (no links are filled in CMSSW_1_6_x)
//     type = 41;
//     TrackRef statTrack = tkmu->standAloneMuon();
//     atm = &(*statTrack);

//     pTrack   = gHFEvent->addSigTrack(); 
//     pTrack->fMuType   = type;
//     pTrack->fMuID     = mcnt;
//     pTrack->fIndex    = idrec;
//     pTrack->fGenIndex = -1; 
//     pTrack->fMCID     = atm->charge()*-13; 
//     pTrack->fQ        = atm->charge();

//     pTrack->fPlab.SetPtEtaPhi(atm->pt(),
// 			      atm->eta(),
// 			      atm->phi()
// 			      );

//     pTrack->fTip  = atm->d0();
//     pTrack->fLip  = atm->dz();
//     pTrack->fChi2 = atm->chi2();
//     pTrack->fDof  = int(atm->ndof());
//     pTrack->fHits = atm->numberOfValidHits();  

    if (fVerbose > 0) {

      cout << "--------------------------------------------------------------------" << endl; 
      cout << "-- Muon " <<fMuonsLabel2 << ", #" << tcnt <<  ",  track index (" 
	   << (tkmu->track()).index() << ")" << endl;
    }
   
//     if ( fVerbose > 0 )
//       cout << "%%> Type-" << type << ": "; 
//       pTrack->dump(); 
//     }
    
    // -- tracker muon
    type = 42;
    TrackRef trktTrack = tkmu->track();
    ttm = &(*trktTrack);

    pTrack   = gHFEvent->addSigTrack(); 
    pTrack->fMuType   = type;
    pTrack->fMuID     = mcnt;
    pTrack->fIndex    = idrec;
    pTrack->fGenIndex = -1; 
    pTrack->fMCID     = ttm->charge()*-13; 
    pTrack->fQ        = ttm->charge();

    pTrack->fPlab.SetPtEtaPhi(ttm->pt(),
			      ttm->eta(),
			      ttm->phi()
			      );

    pTrack->fTip  = ttm->d0();
    pTrack->fLip  = ttm->dz();
    pTrack->fChi2 = ttm->chi2();
    pTrack->fDof  = int(ttm->ndof());
    pTrack->fHits = ttm->numberOfValidHits();  

    if (fVerbose > 0) {
      cout << "%%> Type-" << type <<  ": "; 
      pTrack->dump(); 
    }

    // -- combined muon (no links are filled in CMSSW_1_6_x)
//     type = 43;
//     TrackRef comtTrack = tkmu->combinedMuon();
//     ctm = &(*comtTrack);

//     pTrack   = gHFEvent->addSigTrack(); 
//     pTrack->fMuType   = type;
//     pTrack->fMuID     = mcnt;
//     pTrack->fIndex    = idrec;
//     pTrack->fGenIndex = -1; 
//     pTrack->fMCID     = ctm->charge()*-13; 
//     pTrack->fQ        = ctm->charge();

//     pTrack->fPlab.SetPtEtaPhi(ctm->pt(),
// 			      ctm->eta(),
// 			      ctm->phi()
// 			      );

//     pTrack->fTip  = ctm->d0();
//     pTrack->fLip  = ctm->dz();
//     pTrack->fChi2 = ctm->chi2();
//     pTrack->fDof  = int(ctm->ndof());
//     pTrack->fHits = ctm->numberOfValidHits();  

//     if (fVerbose > 0) {
//       cout << "%%> Type-" << type  << ": "; 
//       pTrack->dump(); 
//     }

    if (fVerbose > 0) {
      cout << "--------------------------------------------------------------------" << endl;; 
    }
  }


  // ===============================================================================================
  // -- Track Collections
  // ===============================================================================================

  // -- globalMuons

  type = 100;

  edm::Handle<TrackCollection> globalMuons;
  iEvent.getByLabel(fTracksLabel1, globalMuons);    

  if (fVerbose > 0) {
    cout << "--------------------------------------------------------------------" << endl; 
    cout << "-- Muon " <<fTracksLabel1 << endl;
  }

  for (unsigned int i = 0; i < globalMuons->size(); ++i){
    
    TrackRef rTrack(globalMuons, i);
    Track track(*rTrack);    

    pTrack = gHFEvent->addSigTrack();
    pTrack->fMuType   = type;
    pTrack->fMuID     = i;
    pTrack->fMCID     = track.charge()*-13; 
    pTrack->fGenIndex = -1;
    pTrack->fIndex    = -1;
    //    pTrack->fIndex = rTrack.index();
    pTrack->fPlab.SetPtEtaPhi(track.pt(),
			      track.eta(),
			      track.phi()
			      );
    pTrack->fTip = track.d0();
    pTrack->fLip = track.dz();
    pTrack->fQ = track.charge();
    pTrack->fChi2 = track.chi2();
    pTrack->fDof = int(track.ndof());
    pTrack->fHits = track.numberOfValidHits();  
    pTrack->fMuID = 1.; 
     
    if (fVerbose > 0) {
      cout << "%%> Type-" << type << ": "; 
      pTrack->dump(); 
    }
  }


  // -- standAloneMuons

  edm::Handle<TrackCollection> staMuons;
  iEvent.getByLabel(fTracksLabel2, staMuons);    

  type = 101;

  if (fVerbose > 0) {
    cout << "--------------------------------------------------------------------" << endl; 
    cout << "-- Muon " <<fTracksLabel2 << endl;
  }

  for (unsigned int i = 0; i < staMuons->size(); ++i){
    
    TrackRef rTrack(staMuons, i);
    Track track(*rTrack);    

    pTrack = gHFEvent->addSigTrack();
    pTrack->fMuType   = type;
    pTrack->fMuID     = i;
    pTrack->fMCID     = track.charge()*-13; 
    pTrack->fGenIndex = -1;
    pTrack->fIndex    = -1;
    //    pTrack->fIndex = rTrack.index();
    pTrack->fPlab.SetPtEtaPhi(track.pt(),
			      track.eta(),
			      track.phi()
			      );
    pTrack->fTip = track.d0();
    pTrack->fLip = track.dz();
    pTrack->fQ = track.charge();
    pTrack->fChi2 = track.chi2();
    pTrack->fDof = int(track.ndof());
    pTrack->fHits = track.numberOfValidHits();  
    pTrack->fMuID = 1.;  
     
    if (fVerbose > 0) {
      cout << "%%> Type-" << type << ": "; 
      pTrack->dump(); 
    }
  }


  // -- standAloneMuons updated at Vtx

  edm::Handle<TrackCollection> staVtxMuons;
  iEvent.getByLabel(fTracksLabel3, staVtxMuons);    

  type = 102;

  if (fVerbose > 0) {
    cout << "--------------------------------------------------------------------" << endl; 
    cout << "-- Muon " <<fTracksLabel3 << endl;
  }

  for (unsigned int i = 0; i < staVtxMuons->size(); ++i){
    
    TrackRef rTrack(staVtxMuons, i);
    Track track(*rTrack);    

    pTrack = gHFEvent->addSigTrack();
    pTrack->fMuType   = type;
    pTrack->fMuID     = i;
    pTrack->fMCID     = track.charge()*-13; 
    pTrack->fGenIndex = -1;
    pTrack->fIndex    = -1;
    //    pTrack->fIndex = rTrack.index();
    pTrack->fPlab.SetPtEtaPhi(track.pt(),
			      track.eta(),
			      track.phi()
			      );
    pTrack->fTip = track.d0();
    pTrack->fLip = track.dz();
    pTrack->fQ = track.charge();
    pTrack->fChi2 = track.chi2();
    pTrack->fDof = int(track.ndof());
    pTrack->fHits = track.numberOfValidHits();  
    pTrack->fMuID = 1.;  
     
    if (fVerbose > 0) {
      cout << "%%> Type-" << type << ": "; 
      pTrack->dump(); 
    }
  }

  // -- hlt2

  edm::Handle<TrackCollection> hlt2Muons;
  iEvent.getByLabel(fTracksLabel4, hlt2Muons);    

  type = 103;

  if (fVerbose > 0) {
    cout << "--------------------------------------------------------------------" << endl; 
    cout << "-- Muon " <<fTracksLabel4 << endl;
  }

  for (unsigned int i = 0; i < hlt2Muons->size(); ++i){
    
    TrackRef rTrack(hlt2Muons, i);
    Track track(*rTrack);    

    pTrack = gHFEvent->addSigTrack();
    pTrack->fMuType   = type;
    pTrack->fMuID     = i;
    pTrack->fMCID     = track.charge()*-13; 
    pTrack->fGenIndex = -1;
    pTrack->fIndex    = -1;
    //    pTrack->fIndex = rTrack.index();
    pTrack->fPlab.SetPtEtaPhi(track.pt(),
			      track.eta(),
			      track.phi()
			      );
    pTrack->fTip = track.d0();
    pTrack->fLip = track.dz();
    pTrack->fQ = track.charge();
    pTrack->fChi2 = track.chi2();
    pTrack->fDof = int(track.ndof());
    pTrack->fHits = track.numberOfValidHits();  
    pTrack->fMuID = 1.;  
     
    if (fVerbose > 0) {
      cout << "%%> Type-" << type << ": "; 
      pTrack->dump(); 
    }
  }  


  // -- hlt2 updated at Vtx

  edm::Handle<TrackCollection> hlt2VtxMuons;
  iEvent.getByLabel(fTracksLabel5, hlt2VtxMuons);    

  type = 104;

  if (fVerbose > 0) {
    cout << "--------------------------------------------------------------------" << endl; 
    cout << "-- Muon " <<fTracksLabel5 << endl;
  }

  for (unsigned int i = 0; i < hlt2VtxMuons->size(); ++i){
    
    TrackRef rTrack(hlt2VtxMuons, i);
    Track track(*rTrack);    

    pTrack = gHFEvent->addSigTrack();
    pTrack->fMuType   = type;
    pTrack->fMuID     = i;
    pTrack->fMCID     = track.charge()*-13; 
    pTrack->fGenIndex = -1;
    pTrack->fIndex    = -1;
    //    pTrack->fIndex = rTrack.index();
    pTrack->fPlab.SetPtEtaPhi(track.pt(),
			      track.eta(),
			      track.phi()
			      );
    pTrack->fTip = track.d0();
    pTrack->fLip = track.dz();
    pTrack->fQ = track.charge();
    pTrack->fChi2 = track.chi2();
    pTrack->fDof = int(track.ndof());
    pTrack->fHits = track.numberOfValidHits();  
    pTrack->fMuID = 1.;  
     
    if (fVerbose > 0) {
      cout << "%%> Type-" << type << ": "; 
      pTrack->dump(); 
    }
  }


  // -- hlt3

  edm::Handle<TrackCollection> hlt3Muons;
  iEvent.getByLabel(fTracksLabel6, hlt3Muons);    

  type = 105;

  if (fVerbose > 0) {
    cout << "--------------------------------------------------------------------" << endl; 
    cout << "-- Muon " <<fTracksLabel6 << endl;
  }

  for (unsigned int i = 0; i < hlt3Muons->size(); ++i){
    
    TrackRef rTrack(hlt3Muons, i);
    Track track(*rTrack);    

    pTrack = gHFEvent->addSigTrack();
    pTrack->fMuType   = type;
    pTrack->fMuID     = i;
    pTrack->fMCID     = track.charge()*-13; 
    pTrack->fGenIndex = -1;
    pTrack->fIndex    = -1;
    //    pTrack->fIndex = rTrack.index();
    pTrack->fPlab.SetPtEtaPhi(track.pt(),
			      track.eta(),
			      track.phi()
			      );
    pTrack->fTip = track.d0();
    pTrack->fLip = track.dz();
    pTrack->fQ = track.charge();
    pTrack->fChi2 = track.chi2();
    pTrack->fDof = int(track.ndof());
    pTrack->fHits = track.numberOfValidHits();  
    pTrack->fMuID = 1.;  
     
    if (fVerbose > 0) {
      cout << "%%> Type-" << type << ": "; 
      pTrack->dump(); 
    }
  }  


  // -- hlt3 seeded

  edm::Handle<TrackCollection> hlt3seededMuons;
  iEvent.getByLabel(fTracksLabel7, hlt3seededMuons);    

  type = 106;

  if (fVerbose > 0) {
    cout << "--------------------------------------------------------------------" << endl; 
    cout << "-- Muon " <<fTracksLabel7 << endl;
  }

  for (unsigned int i = 0; i < hlt3seededMuons->size(); ++i){
    
    TrackRef rTrack(hlt3seededMuons, i);
    Track track(*rTrack);    

    pTrack = gHFEvent->addSigTrack();
    pTrack->fMuType   = type;
    pTrack->fMuID     = i;
    pTrack->fMCID     = track.charge()*-13; 
    pTrack->fGenIndex = -1;
    pTrack->fIndex = -1;
    //    pTrack->fIndex = rTrack.index();
    pTrack->fPlab.SetPtEtaPhi(track.pt(),
			      track.eta(),
			      track.phi()
			      );
    pTrack->fTip = track.d0();
    pTrack->fLip = track.dz();
    pTrack->fQ = track.charge();
    pTrack->fChi2 = track.chi2();
    pTrack->fDof = int(track.ndof());
    pTrack->fHits = track.numberOfValidHits();  
    pTrack->fMuID = 1.;  
     
    if (fVerbose > 0) {
      cout << "%%> Type-" << type << ": "; 
      pTrack->dump(); 
    }
  }




}

// ------------ method called once each job just before starting event loop  ------------
void  BmmDumpMuons::beginJob(const edm::EventSetup& setup) {
}

// ------------ method called once each job just after ending the event loop  ------------
void  BmmDumpMuons::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(BmmDumpMuons);
