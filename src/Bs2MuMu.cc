// -*- C++ -*-
//
// Package:    Bs2MuMu
// Class:      Bs2MuMu
// 
/**\class Bs2MuMu Bs2MuMu.cc HeavyFlavorAnalysis/Bs2MuMu/src/Bs2MuMu.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Christina Eggel
//         Created:  Mon Oct 23 15:14:30 CEST 2006
// $Id: Bs2MuMu.cc,v 1.8 2007/09/28 16:33:47 ceggel Exp $
//
//

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/Bs2MuMu.h"


#include <memory>
#include <iostream>
#include <string>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TParameter.h>
#include <TH1.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>

// ----------------------------------------------------------------------

using namespace reco;
using namespace std;
using namespace edm;
using namespace HepMC;


class TrackAssociatorByHits;
class TrackAssociatorByChi2; 
class VertexAssociatorByTrack;
class TrackerHitAssociator;

// ----------------------------------------------------------------------
struct anaStuff {

  reco::VertexCollection *theVtxCollection;

  const reco::MuonCollection       *theMuonCollection;
  const reco::TrackCollection      *theTkCollection;
  const TrackingParticleCollection *theTPCollection;

  reco::RecoToSimCollection        *recSimCollection;
  //  reco::VertexRecoToSimCollection  *recSimCollectionVertex;

  // ----------Tracks--------------------------     
       
  // -- Muons (TM=2 or GM=3) 
  std::vector<const reco::Track*> MuonRecTracks; 
  std::vector<int> MuonRecTracksIndex;

  // -- Muons (full TM=1)
  std::vector<const reco::Track*> BmmRecTracks;
  std::vector<int> BmmRecTracksIndex;          
  std::vector<int> BmmRecTracksB;                 // Index of B-mother

  std::vector<const reco::Track*> JpsiRecTracks; 
  std::vector<int> JpsiRecTracksIndex;
  std::vector<int> JpsiRecTracksB;                // Index of B-grandmother

  // -- Kaons (full TM=1)
  std::vector<const reco::Track*> KaonRecTracks; 
  std::vector<int> KaonRecTracksIndex;
  std::vector<int> KaonRecTracksB;                // Index of B-mother


  // --------Candidates------------------------

  // -- Muons signal
  std::vector< std::pair<const reco::Track*, const reco::Track*> > BmmPairTracks;
  std::vector< std::pair<int, int> > BmmPairTracksIndex;

  // -- Muons norm
  std::vector< std::pair<const reco::Track*, const reco::Track*> > JpsiPairTracks;
  std::vector< std::pair<int, int> > JpsiPairTracksIndex;

  // -- Kaons norm
  std::vector<const reco::Track*> KaonTrack; 
  std::vector<int> KaonTrackIndex;
  // ------------------------------------------
  // -- Vertex tracks
  std::vector<const reco::Track*> RecTracks; 
  std::vector<int> RecTracksIndex;
  std::vector<reco::Track> RefittedTracks;
  // ------------------------------------------ 


  int fBmmSel;

  reco::Vertex primaryVertex;
  reco::Vertex primaryVertex2;

  std::vector<double> SecVtxChi2;
  std::vector<double> InvMass;

};

// ----------------------------------------------------------------------
// ======================================================================

Bs2MuMu::Bs2MuMu(const edm::ParameterSet& iConfig) {

  cout << endl << "===>> Bs2MuMu >>> ctor, instantiating histogramms, etc." << endl;

  // -- Setup "class variables"
  fStuff = new anaStuff;

  // -- Counters
  fNevt = 0; 
  fNgen = 0;
  fNrec = 0;


  // -- Config. File input
  //  fLabel          = iConfig.getUntrackedParameter("moduleLabel",std::string("source"));
  fSourceLabel    = iConfig.getParameter<string>("HepMC");
  fTracksLabel    = iConfig.getParameter<string>("tracks");
  fAssocLabel     = iConfig.getParameter<string>("associator");
  fMuonLabel      = iConfig.getParameter<string>("Muons");

  fStuff->fBmmSel = iConfig.getParameter<int>("bmmsel");

  fVerbose     = iConfig.getParameter<int>("Verbose");
  fGenVerbose     = iConfig.getParameter<int>("GenVerbose");
  fSimVerbose     = iConfig.getParameter<int>("SimVerbose");
  fRecVerbose     = iConfig.getParameter<int>("RecVerbose");
  fGlbVerbose     = iConfig.getParameter<int>("GlbVerbose");
  fR2SVerbose     = iConfig.getParameter<int>("Rec2SimVerbose");

  if ( strcmp( (fAssocLabel.c_str()), "TrackAssociatorByHits") && 
       strcmp( (fAssocLabel.c_str()), "TrackAssociatorByChi2") )  {
  
    cout << " Please set your track associator option to either \"TrackAssociatorByHits\"                                                              or \"TrackAssociatorByChi2\" .... Abort!" << endl;
    return;
  }


  // -- ROOT output
  fChannel = iConfig.getParameter<string>("channel");
  fFile = new TFile(iConfig.getParameter<string>("fileName").c_str(), "RECREATE");
  fTree = new TTree("T1","CMSSW Bs -> mu+mu- tree");
  fEvent = new TAna00Event(0);
  fTree->Branch("TAna00Event", "TAna00Event", &fEvent, 256000/8, 1);

  // -- Decay mode
  decayChannel(fChannel.c_str());

  // -- Troubleshoot histogramm
  fEff = new TH1D("eff", "Efficiencies", 100, 0., 1000. );

  // -- Invariant mass & control histograms

  for (int i = 0; i < 3; i++) {

    fM000[i]  = new TH1D(Form("m000_%i", i+1), Form("inv. Mass Cand0 (all kaon cand.), sel = %i", i+1), 1000, 0., 10. );
    fM100[i]  = new TH1D(Form("m531_%i", i+1), Form("inv. Mass Cand1 (bmm), sel = %i", i+1), 1000, 0., 10. );
    fM200[i]  = new TH1D(Form("m443_%i", i+1), Form("inv. Mass Cand2 (jpsi), sel = %i", i+1), 1000, 0., 10. );
    fM300[i]  = new TH1D(Form("m521_%i", i+1), Form("inv. Mass Cand3 (bjk), sel = %i", i+1), 1000, 0., 10. );
  }

  // -- Kaon Selection
  fK100  = new TH2D("k100", "inv. Mass Cand2 vs. Chi2 (candidates)", 500, 0., 50., 500, 0., 100. );
  fK200  = new TH2D("k200", "inv. Mass Cand2 vs. Chi2 (selected cand.)", 500, 0., 50., 500, 0., 100. );

  // -- pT Resolution (filled in printReco2Sim)
  if ( fR2SVerbose == 1) {

    fPT300 = new TH1D("p300", "pT Resoultion (hits)", 400, -2., 2. );
    fPT310 = new TH1D("p310", "pT Resoultion (chi2)", 400, -2., 2. );
    fPT320 = new TH1D("p320", "pT Resoultion (hits)", 400, -2., 2. );
  }
}

// ======================================================================
// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// ======================================================================

Bs2MuMu::~Bs2MuMu() {

  cout << "===>> Bs2MuMu >>> dtor, writing histogramms to file" << endl;

  // -- Save output
  fFile->cd();
 
  fEff->Write();

  for (int i = 0; i < 3; i++) {
   
    fM000[i]->Write(); 
    fM100[i]->Write(); 
    fM200[i]->Write();
    fM300[i]->Write();
  }

  fK100->Write(); 
  fK200->Write();

  if ( fR2SVerbose == 1) {
   
    fPT300->Write(); 
    fPT310->Write(); 
    fPT320->Write();
  } 

  fTree->Write();

  fFile->Write();
  fFile->Close();
  delete fFile;

}

// ======================================================================
// ----------------------------------------------------------------------



// ------------ method called to for each event  ------------
void Bs2MuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  ++fNevt; 

  cout << endl << "*********************************************************************" << endl;
  cout << "===>> Bs2MuMu >>> Start with event: " << fNevt << endl;

  using reco::TrackCollection;
  using reco::MuonCollection;

  // === Initialize event record ===

  fEvent->Clear();

  fEvent->fRunNumber   = iEvent.id().run();
  fEvent->fEventNumber = iEvent.id().event();

  fStuff->theVtxCollection = 0;
  fStuff->theMuonCollection = 0;
  fStuff->theTkCollection = 0;
  fStuff->theTPCollection = 0;
  fStuff->recSimCollection = 0;

  clearTracks();
  clearCandidateTracks();


  fStuff->RecTracks.clear();
  fStuff->RecTracksIndex.clear();

//   fStuff->primaryVertex  = 0;
//   fStuff->primaryVertex2 = 0;

  // === Fill event record ===

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  fStuff->theTkCollection  = tracks.product();
  // fStuff->theTkCollection  = new reco::TrackCollection(*(tracks.product()));

  // -- get the collection of TrackingParticles 
  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
  iEvent.getByLabel("trackingtruth","TrackTruth",TPCollectionH);
  fStuff->theTPCollection  = TPCollectionH.product();  
  // fStuff->theTPCollection  = new TrackingParticleCollection(*(TPCollectionH.product()));  
   
  // -- get the collection of MuonTracks 
  edm::Handle<reco::MuonCollection> MuCollection;
  iEvent.getByLabel(fMuonLabel.c_str(), MuCollection);
  fStuff->theMuonCollection   = MuCollection.product();
  // fStuff->theMuonCollection   = new reco::MuonCollection(*(MuCollection.product()));  

  // -- get the collection of sim. Tracks
//   Handle<SimTrackContainer> simTrackCollection;
//   iEvent.getByLabel("g4SimHits", simTrackCollection);
//   const SimTrackContainer simTC = *(simTrackCollection.product());

  // -- get the collection of sim. Vertices
//   Handle<SimVertexContainer> simVertexCollection;
//   iEvent.getByLabel("g4SimHits", simVertexCollection);
//   const SimVertexContainer simVC = *(simVertexCollection.product());

  // -- get the collection of TrackingVertices
  edm::Handle<TrackingVertexCollection>  TVCollectionH ;
  iEvent.getByLabel("trackingtruth","VertexTruth",TVCollectionH);
  const TrackingVertexCollection tVC   = *(TVCollectionH.product());

  // -- get the collection of primary Vertices
  edm::Handle<reco::VertexCollection>  primaryVertexH ;
  iEvent.getByLabel("offlinePrimaryVerticesFromCTFTracks","",primaryVertexH);
  const reco::VertexCollection primaryVertexCollection   = *(primaryVertexH.product());

  // -- perform association by either hits or chi2 (in config-file)
  if ( !strcmp( (fAssocLabel.c_str()), "TrackAssociatorByHits") ) {
      
    if (fVerbose) cout << "==> -- Track Associator by hits --" << endl;  
    fStuff->recSimCollection = new 
      reco::RecoToSimCollection(associatorByHits->associateRecoToSim(tracks, TPCollectionH, &iEvent));
  }
  if ( !strcmp( (fAssocLabel.c_str()), "TrackAssociatorByChi2") ) {
      
    if (fVerbose) cout << "==> -- Track Associator by chi2 --" << endl;  
    fStuff->recSimCollection = new   
      reco::RecoToSimCollection(associatorByChi2->associateRecoToSim(tracks, TPCollectionH, &iEvent));
  }
    

//**  if (fVerbose) cout << "==> -- Vertex Associator by tracks --" << endl; 
//**   fStuff->recSimCollectionVertex = new
//**     reco::VertexRecoToSimCollection(associatorByTracks->associateRecoToSim(primaryVertexH, TVCollectionH, iEvent
//**  									   , (*fStuff->recSimCollection)));
  // === Start analysis ===

  // -- Generator level
  fillGeneratorBlock(iEvent);

  // -- Print Gen-/Sim-/Rec-Block
  if ( fGenVerbose == 1 ) {
    
    printGenTracks(iEvent);
  }
  if ( fSimVerbose == 1 ) {
    
    printSimTracks(iEvent);
  }
  if ( fRecVerbose == 1 ) {
    
    printRecTracks(iEvent);
  }
  if ( fGlbVerbose == 1 ) {
    
    printMuonTracks(iEvent);
  }

  // -- Reconstructed Tracks
  fillRecTracks(iEvent);

  // -- Signal Tracks (Global Muons)
  trimBmmTracks(iEvent);

  // -- Primary Vertex
  int npv = primaryVertex(iEvent);
  if (!npv) {
    cout << "*********** No primary vertex found !!!" << endl;  
  }

  // -- (Rec.) Signal Tracks
  if (fStuff->fBmmSel == 1) {

    bmmTracks1(iEvent);
    truthCandTracks(iEvent, iSetup);
    secondaryVertex(iEvent, iSetup);

  } else if (fStuff->fBmmSel == 2) {

    bmmTracks2(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

  } else if (fStuff->fBmmSel == 3) {

    bmmTracks3(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);
  } else {

    fStuff->fBmmSel = 1;
    bmmTracks1(iEvent);
    truthCandTracks(iEvent, iSetup);
    secondaryVertex(iEvent, iSetup);

    fStuff->fBmmSel = 2;
    bmmTracks2(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

    fStuff->fBmmSel = 3;
    bmmTracks3(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

    fStuff->fBmmSel = -1;
  }

  // -- Dump tree
  fTree->Fill();

  cout << endl << "===>> Bs2MuMu >>> Done with event: " << fNevt << endl;
  cout << "*********************************************************************" << endl;
    
}



// ----------------------------------------------------------------------
// ======================================================================

void Bs2MuMu::fillGeneratorBlock(const edm::Event &iEvent) {

  fNgen++;

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>fillGeneratorBlock> Starting to fill generator block, event: " << fNevt << endl;

  TGenCand *pCand;

  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel(fSourceLabel.c_str(), evt);
  const HepMC::GenEvent *genEvent = evt->GetEvent();
  

  int gcnt = genEvent->particles_size();
  if (fVerbose) cout << "Counted " << gcnt  << " genTracks in generator block" << endl;
  gcnt = 0;

  int aid(0); 
  int tmpDau1(-1);
  int tmpDau2(-1);
  int scr0(0), scr1(0), scr2(0);
  for (HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin();
       p != genEvent->particles_end();
       ++p) {
   
    
    pCand = fEvent->addGenCand();
    pCand->fNumber = gcnt;

    pCand->fID     = (*p)->ParticleID();
    pCand->fStatus = (*p)->StatusCode();
    pCand->fMass   = (*p)->Mass();
    pCand->fMom1   = (*p)->Mother()-1;
    pCand->fMom2   = (*p)->SecondMother()-1;


//     pCand->fDau1   = (*p)->BeginDaughters()-1;
//     pCand->fDau2   = (*p)->EndDaughters()-1; // this is screwed!?

    tmpDau1   = (*p)->BeginDaughters()-1;
    tmpDau2   = (*p)->EndDaughters()-1; // this is screwed!?

    pCand->fDau1 = -1;
    pCand->fDau2 = -1;

    std::vector<HepMC::GenParticle*> children = (*p)->listChildren();
    std::vector<HepMC::GenParticle*>::const_iterator aDaughter;

    if ((*p)->numChildren() > 1) {
      
      int first(1); 
      for (aDaughter = children.begin(); aDaughter != children.end(); aDaughter++) {
	if (first == 1) {
	  pCand->fDau1 = (*p)->BeginDaughters()-1;
	  first = 0;
	}
	pCand->fDau2   = (*aDaughter)->barcode()-1;
      } 
    } else if ((*p)->numChildren() > 0 ) {

      pCand->fDau1 = (*p)->BeginDaughters()-1;
      pCand->fDau2 = -1;

    } else {

            pCand->fDau1   = -1;
      pCand->fDau2   = -1;
    }

    // -- Sometimes the pointers are swapped?!
    if (pCand->fDau1 > pCand->fDau2) {
      aid          = pCand->fDau1; 
      pCand->fDau1 = pCand->fDau2; 
      pCand->fDau2 = aid;
    } 

    // -- Print-out of Bs-Muon-Indeces in Generator Block
     if ( pCand->fDau1 != tmpDau1 && pCand->fDau2 == tmpDau2 ) {
//        if (fVerbose) cout << " SCREWED INDEED --> 1: DAU1 "
//  	   << tmpDau1 << " != " << pCand->fDau1 << " !!!!!!!!!!" << endl;
       scr0++;
     }
     if ( pCand->fDau1 == tmpDau1 && pCand->fDau2 != tmpDau2 ) {
//        if (fVerbose) cout << " SCREWED INDEED --> 1: DAU2 "
//  	   << tmpDau2 << " != " << pCand->fDau2 << " !!!!!!!!!!" << endl;
       scr1++;
     }
     if ( pCand->fDau1 != tmpDau1 && pCand->fDau2 != tmpDau2 ) {
//        if (fVerbose) cout << " SCREWED INDEED --> 2: DAU1 "
//  	   << tmpDau1 << " != " << pCand->fDau1 << " and DAU2 "
//  	   << tmpDau2 << " != " << pCand->fDau2 << " !!!!!!!!!!" << endl;
       scr2++;
     }

    if (abs(pCand->fID) == 531) {

      if (fVerbose) cout << "Daughters: #" << (*p)->beginDaughters()->barcode() - 1
	   << " - #" << (*p)->endDaughters()->barcode() - 1 << endl;
      if (fVerbose) cout << "1. Muon (fDau1)    " << pCand->fDau1 << endl;
      if (fVerbose) cout << "2. Muon (fDau2)    " << pCand->fDau2 << endl;
    }


    pCand->fP.SetPxPyPzE((*p)->Momentum().px(),
  			 (*p)->Momentum().py(),
  			 (*p)->Momentum().pz(),
  			 (*p)->Momentum().e()
  			 );
    
    pCand->fV.SetXYZ((*p)->DecayVertex().x(),
  		     (*p)->DecayVertex().y(),
  		     (*p)->DecayVertex().z()
  		     );

    ++gcnt;



    aid = abs( (*p)->pdg_id() );
    //     if (5 == aid) {
    //       if (fVerbose) cout << "  ----> Found a b quark ..." << endl;
    //     }

    if (531 == aid) {
      if (fVerbose) cout << "  ----> Found a Bs -- event " << fNgen << " should be kept!! pT = " << (*p)->Momentum().perp() << endl;
    }
  }

  fEff->AddBinContent(10, gcnt);
  fEff->AddBinContent(11, scr0);
  fEff->AddBinContent(12, scr1);
  fEff->AddBinContent(13, scr2);

//   if (fVerbose) cout << "============================================================================================" << endl;
//   fEvent->dumpGenBlock();
//   if (fVerbose) cout << "============================================================================================" << endl;
//   if (fVerbose) cout << "==>Bs2MuMu> Done with generator block." << endl;

}


// ----------------------------------------------------------------------
void Bs2MuMu::fillRecTracks(const edm::Event &iEvent) {
  
  fNrec++;

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>fillRecTracks> Starting to fill reconstructed tracks, event: " << fNevt << endl;

  TAnaTrack *pTrack;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product());

  int rcnt = tracks->size();
  if (fVerbose) cout << "Counted " << rcnt  << " reconstructed tracks." << endl;
  rcnt = 0;

  // -- get the collection of TrackingParticles 
  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
  iEvent.getByLabel("trackingtruth","TrackTruth",TPCollectionH);
  const TrackingParticleCollection tPC   = *(TPCollectionH.product()); 

  int scnt = TPCollectionH->size();
  if (fVerbose) cout << "Counted " << scnt  << " simulated tracks." << endl;
  scnt = 0;

  // -- Print Reco-To-Sim-Levels
  if ( fR2SVerbose == 1) {

    printReco2Sim(iEvent, "input");
  }
  else if ( fR2SVerbose == 2) {

    printReco2Sim(iEvent, "hits");
    printReco2Sim(iEvent, "chi2");
  }

  // -- Fill Event
  for(TrackCollection::size_type i=0; i<recTC.size(); ++i) {

    TrackRef track(tracks, i);

    reco::Track tt(*track);

    pTrack = fEvent->addRecTrack();
    pTrack->fIndex = track.index();

    fillTrack(iEvent, pTrack, &tt, i, 0);

  }

  for (MuonCollection::const_iterator glbMuon = (*fStuff->theMuonCollection).begin(); 
       glbMuon != (*fStuff->theMuonCollection).end(); 
       ++glbMuon){

    TrackRef glbTrack = glbMuon->combinedMuon();
    const reco::Track* tt = &(*glbTrack);

    int idrec = idRecTrack(tt);

    // -- trouble histo for muon2reco-----------------
    fEff->Fill(50.1);
    if ( idrec > -1 ) {
      fEff->Fill(51.1);
    } else {
      fEff->Fill(52.1);
    }
    // -----------------------------------------------
  }
}


// ----------------------------------------------------------------------
void Bs2MuMu::bmmTracks1(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>bmmTracks1> Matching particle, truth-matched to be particle PDG #" 
       << fTruthMC_I << " and a decay products of particle PDG #" << fTruthMC_mom << ", event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product()); 
 
 // -- track association
  reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection); 

  // -- Clear tracks
  clearTracks();
  clearCandidateTracks();

 // -- perform association by either hits or chi2
//   reco::RecoToSimCollection recSimColl; 
//   if ( !strcmp( (fAssocLabel.c_str()), "TrackAssociatorByHits") ) {

//     if (fVerbose) cout << "==>bmmTracks1> -- Associator by hits --" << endl;  
//     recSimColl = associatorByHits->associateRecoToSim (tracks, TPCollectionH, &iEvent);
//   }
//   if ( !strcmp( (fAssocLabel.c_str()), "TrackAssociatorByChi2") ) {

//     if (fVerbose) cout << "==>bmmTracks1> -- Associator by chi2 --" << endl;  
//     recSimColl = associatorByChi2->associateRecoToSim (tracks, TPCollectionH, &iEvent );
//   }

  const reco::Track* tt = 0;
  const HepMC::GenParticle* genGmo = 0;
  const HepMC::GenParticle* genMom = 0;
  const HepMC::GenParticle* genPar = 0; 

  int mcand1(0), mcand2(0), kcand(0);
  int gmo_pdg_id(-99999), gmo_id(-99999);
  int mom_pdg_id(-99999), mom_id(-99999);
  int gen_pdg_id(-99999), gen_id(-99999);
  int gen_cnt(0);

  for(TrackCollection::size_type i=0; i<recTC.size(); ++i) {

    TrackRef track(tracks, i);

    try{ 

      std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
//       for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
// 	     it != tp.end(); ++it) {

      TrackingParticleRef tpr = tp.begin()->first;  // ??? This takes the first associated simulated track ???	    
      double assocChi2 = tp.begin()->second;

      if (fVerbose) cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-." <<endl;
      
      gen_cnt = 0; gen_pdg_id = -99999; gen_id = -99999;

      for ( TrackingParticle::genp_iterator pit=tpr->genParticle_begin(); pit!=tpr->genParticle_end(); pit++ ){
	
	genPar=pit->get();
	
	if (fVerbose) {

	  cout << ".. Rec. Track #" << track.index() << " pT = " << track->pt() << endl;
	  cout << ".. Sim. Particle #" << tpr.index() << " pT = " << tpr->pt()  
	       << " PDG ID: " << tpr->pdgId() << " (NShared: "  << assocChi2 << ")" << endl;
	}

	gen_pdg_id = genPar->ParticleID();
	gen_id     = genPar->barcode()-1;
	
	if (fVerbose) cout << ".. Gen. Particle #" << gen_id << " pT = " << genPar->Momentum().perp() << " ---> PDG ID: " << gen_pdg_id << endl;
	


	// -- Two body decays
	if ( (abs(gen_pdg_id) == fTruthMC_I) ||  
	     (abs(gen_pdg_id) == fTruthMC_II) ) { 
	  
	  genMom = genPar->mother();
	  mom_pdg_id = genMom->ParticleID();
	  mom_id     = genMom->barcode()-1;
	  
	  
	  if (fVerbose) cout << endl << ".. and Mother Particle #" << mom_id << " pT = " << genMom->Momentum().perp() 
			     << " ---> PDG ID: " << mom_pdg_id;
	  
	  // -- Cand 1 (Bs)
	  if ( abs(mom_pdg_id) == fTruthMC_mom ) { 
	    
	    mcand1++;
	    
	    tt = &(*track);
	    
	    fStuff->BmmRecTracks.push_back(tt);
	    fStuff->BmmRecTracksIndex.push_back(track.index());
	    fStuff->BmmRecTracksB.push_back(mom_id);

	    fStuff->MuonRecTracks.push_back(tt);
	    fStuff->MuonRecTracksIndex.push_back(track.index());
	    
	    if (fVerbose) cout << "    *** " << fPrintChannel << " (" << mcand1 << ") *** ";
	  }	  
	  
	  // -- Cand 2 (J/Psi)
	  if ( abs(mom_pdg_id) == fTruthMC2_mom ) { 
	    
	    genGmo = genMom->mother();
	    gmo_pdg_id = genGmo->ParticleID();
	    gmo_id     = genGmo->barcode()-1;
	    
	    if (fVerbose) cout << endl << ".. and Grandmother Particle #" << gmo_id << " pT = " 
			       << genGmo->Momentum().perp() << " ---> PDG ID: " << gmo_pdg_id;
	    
	    if ( abs(gmo_pdg_id) == fTruthMC2_gmo ) { 
	      
	      mcand2++;
	      
	      tt = &(*track);
	      
	      fStuff->JpsiRecTracks.push_back(tt);
	      fStuff->JpsiRecTracksIndex.push_back(track.index());
	      fStuff->JpsiRecTracksB.push_back(gmo_id);

	      fStuff->MuonRecTracks.push_back(tt);
	      fStuff->MuonRecTracksIndex.push_back(track.index());
	      
	      if (fVerbose) cout << "    *** " << fPrintChannel2 << " (" << mcand2 << ". mu) *** ";
	    }
	    
	  }
	  
	  if (fVerbose) cout << endl; 
	}
	
	// -- Cand 3 (only Kaon)
	if ( abs(gen_pdg_id) == fTruthMC2  ) { 
	  
	  genMom = genPar->mother();
	  mom_pdg_id = genMom->ParticleID();
	  mom_id     = genMom->barcode()-1;
	  
	  
	  if (fVerbose) cout << endl << ".. and Mother Particle #" << mom_id << " pT = " << genMom->Momentum().perp() 
			     << " ---> PDG ID: " << mom_pdg_id;
	  
	  if ( abs(mom_pdg_id) == fTruthMC2_gmo ) { 
	    
	    kcand++;
	    
	    tt = &(*track);
	    
	    fStuff->KaonRecTracks.push_back(tt);
	    fStuff->KaonRecTracksIndex.push_back(track.index());
	    fStuff->KaonRecTracksB.push_back(mom_id);
	    
	    if (fVerbose) cout << "    *** " << fPrintChannel2 << " (" << kcand << ". K) *** ";
	  }
	  
	  if (fVerbose) cout << endl; 
	}
	     
	gen_cnt++;	 
      }


	
      if ( gen_cnt == 0 ) {
	  
	if (fVerbose) cout << "%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
      }

	if ( gen_cnt > 1 ) {
	
	if (fVerbose) cout << " ==> !!! More than one gen. particle for sim. particle #" << tpr.index() << " !!!" << endl;
      }
    } catch (Exception event) {

      if (fVerbose) cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	   << setprecision(2) << setw(6) << track->pt() 
	   <<  " matched to 0 MC Tracks" << endl;
    } 
    
  } 
}


// ----------------------------------------------------------------------
void Bs2MuMu::bmmTracks2(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>bmmTracks2> Matching particle, truth matched to be particle PDG #" 
       << fTruthMC_I << " in the generator block, event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product());

 // -- track association
  reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection);

  // -- Clear tracks
  clearTracks();
  clearCandidateTracks();

  const reco::Track* tt = 0;
  const HepMC::GenParticle* genGmo = 0;
  const HepMC::GenParticle* genMom = 0;
  const HepMC::GenParticle* genPar = 0; 

  int mcand(0);
  int gmo_pdg_id(-99999), gmo_id(-99999);
  int mom_pdg_id(-99999), mom_id(-99999);
  int gen_pdg_id(-99999), gen_id(-99999);
  int gen_cnt(0);

  for(TrackCollection::size_type i=0; i<recTC.size(); ++i) {

    TrackRef track(tracks, i);

    try{ 

      std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
//       for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
// 	     it != tp.end(); ++it) {

      TrackingParticleRef tpr = tp.begin()->first;  // ??? This takes the first associated simulated track ???	    
      double assocChi2 = tp.begin()->second;

      if (fVerbose) cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-." <<endl;
      
      gen_cnt = 0; gen_pdg_id = -99999; gen_id = -99999;

      for ( TrackingParticle::genp_iterator pit=tpr->genParticle_begin(); pit!=tpr->genParticle_end(); pit++ ){
	
	genPar=pit->get();
		
	gen_pdg_id = genPar->ParticleID();
	gen_id     = genPar->barcode()-1;

	if ( (abs(gen_pdg_id) == fTruthMC_I) || 
	     (abs(gen_pdg_id) == fTruthMC_II) ) { 
	  
	  mcand++;
	  tt = &(*track);

	  if (fVerbose) {
	    
	    cout << ".. Rec. Track #" << track.index() << " pT = " << track->pt() << endl;
	    cout << ".. Sim. Particle #" << tpr.index() << " pT = " << tpr->pt() 
		 << " PDG ID: " << tpr->pdgId() << " (NShared: "  << assocChi2 << ")" << endl;
	    cout << ".. Gen. Particle #" << gen_id << " pT = " << genPar->Momentum().perp() 
		 << " ---> PDG ID: " << gen_pdg_id << endl;
	  
	  }

	  genMom = genPar->mother();
	  mom_pdg_id = genMom->ParticleID();
	  mom_id     = genMom->barcode()-1;
     

	  if (fVerbose) cout << endl << ".. Mother Particle #" << mom_id << " pT = " << genMom->Momentum().perp() 
	               << " ---> PDG ID: " << mom_pdg_id;

	  genGmo = genMom->mother();
	  gmo_pdg_id = genGmo->ParticleID();
	  gmo_id     = genGmo->barcode()-1;
	  
	  if (fVerbose) cout << endl << ".. and Grandmother Particle #" << gmo_id << " pT = " 
			     << genGmo->Momentum().perp() << " ---> PDG ID: " << gmo_pdg_id;

	    
	  fStuff->MuonRecTracks.push_back(tt);
	  fStuff->MuonRecTracksIndex.push_back(track.index());
 
	  if (fVerbose) cout << "    *** " << fPrintChannel << " (" << mcand << ") *** " << endl;
	}

	gen_cnt++;	 
      }

      // -- eff. histogram ---------------------------

      if ( fStuff->MuonRecTracks.size() >= 0 ) {

	if ( fStuff->MuonRecTracks.size() < 9 ) {
	  
	  fEff->Fill(60.1 + fStuff->MuonRecTracks.size());

	} else {
	  
	  fEff->Fill(69.1);
	}
      }

      // ----------------------------------------------

      if ( gen_cnt == 0 ) {
	
	if (fVerbose) cout << "%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
      }

      if ( gen_cnt > 1 ) {
	
	if (fVerbose) cout << " ==> !!! More than one gen. particle for sim. particle #" << tpr.index() << " !!!" << endl;
      }
    } catch (Exception event) {

      if (fVerbose) cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	   << setprecision(2) << setw(6) << track->pt() 
	   <<  " matched to 0 MC Tracks" << endl;
    }
  }
}


// ----------------------------------------------------------------------
void Bs2MuMu::bmmTracks3(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>bmmTracks3> Matching muons, rec PID as muons, event: " << fNevt << endl;

  // -- Clear tracks
  clearTracks();
  clearCandidateTracks();

  if (0 == fStuff->theMuonCollection ) {
    
    if (fVerbose) cout << "==>bmmTrack3> getting the muon collection ..." << endl;
    
    edm::Handle<reco::MuonCollection> MuCollection;
    iEvent.getByLabel(fMuonLabel.c_str(), MuCollection);
    fStuff->theMuonCollection   = new reco::MuonCollection(*(MuCollection.product()));
  }

  int nMuon  = fStuff->theMuonCollection->size();
  int nTrack = fStuff->theTkCollection->size();

  if (fVerbose) cout << "==>bmmTracks3>  tracks found: " << nTrack << endl
       << "==>bmmTracks3>  muons  found: " << nMuon  << endl;  

//   reco::MuonRefVector glbMuons1, glbMuons2;
  int mcnt(0);

  const reco::Track* tt = 0;

  for (MuonCollection::const_iterator glbMuon = (*fStuff->theMuonCollection).begin(); 
       glbMuon != (*fStuff->theMuonCollection).end(); 
       ++glbMuon){

    TrackRef glbTrack = glbMuon->combinedMuon();
    tt = &(*glbTrack);

    int idrec = idRecTrack(tt);

    if ( idrec > -1 ) {

	    
      fStuff->MuonRecTracks.push_back(tt);
      fStuff->MuonRecTracksIndex.push_back(idrec);

    } else {
      if (fVerbose) cout << "==>bmmTracks3> Coulnd't find rec. track for muon #" << mcnt << endl;  
    }
    
    mcnt++;
//     glbMuons1.push_back(MuonRef(MuCollection,mcnt-1));
//     glbMuons2.push_back(MuonRef(MuCollection,mcnt-1));
  }

  // -- eff. histogram ---------------------------
  
  if ( fStuff->MuonRecTracks.size() >= 0 ) {
    
    if ( fStuff->MuonRecTracks.size() < 9 ) {
      
      fEff->Fill(70.1 + fStuff->MuonRecTracks.size());
      
    } else {
      
      fEff->Fill(79.1);
    }
  }
  // ----------------------------------------------
}


// ----------------------------------------------------------------------
void Bs2MuMu::truthCandTracks(const edm::Event &iEvent, const edm::EventSetup& iSetup) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>truthCandTracks> Sort muon and kaon candidates depending on ID of B-mother, event: " 
		     << fNevt << endl;

  if ( fStuff->BmmRecTracks.size() > 1 ) {

    for ( unsigned int i=0; i<fStuff->BmmRecTracks.size(); i++ ) { 

      for ( unsigned int j=i; j<fStuff->BmmRecTracks.size(); j++ ) { 

	if (i == j) { continue; }

	if (fStuff->BmmRecTracksB[i] == fStuff->BmmRecTracksB[j]) {
	  
	  if ((*fStuff->BmmRecTracks[j]).pt() > (*fStuff->BmmRecTracks[i]).pt() ) { 
	    fStuff->BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (fStuff->BmmRecTracks[j],fStuff->BmmRecTracks[i]));
	    fStuff->BmmPairTracksIndex.push_back(pair<int, int>
						 (fStuff->BmmRecTracksIndex[j],fStuff->BmmRecTracksIndex[i]));
	  } else {
	    fStuff->BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (fStuff->BmmRecTracks[i],fStuff->BmmRecTracks[j]));
	    fStuff->BmmPairTracksIndex.push_back(pair<int, int>
						 (fStuff->BmmRecTracksIndex[i],fStuff->BmmRecTracksIndex[j]));
	  }

	}
      }
    }
  }  

  
  if ( fStuff->JpsiRecTracks.size() > 1 ) {
    
    for ( unsigned int i=0; i<fStuff->JpsiRecTracks.size(); i++ ) { 
      
      for ( unsigned int j=i; j<fStuff->JpsiRecTracks.size(); j++ ) { 
	
	if (i == j) { continue; }
	
	if (fStuff->JpsiRecTracksB[i] == fStuff->JpsiRecTracksB[j]) {
	  

	  if ( fStuff->KaonRecTracks.size() > 0 ) {

	    for ( unsigned int k=0; k<fStuff->KaonRecTracks.size(); k++ ) {
	      
	      if (fStuff->JpsiRecTracksB[i] == fStuff->KaonRecTracksB[k]) {
		
		if ((*fStuff->JpsiRecTracks[j]).pt() > (*fStuff->JpsiRecTracks[i]).pt() ) { 
		  fStuff->JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
						   (fStuff->JpsiRecTracks[j],fStuff->JpsiRecTracks[i]));
		  fStuff->JpsiPairTracksIndex.push_back(pair<int, int>
							(fStuff->JpsiRecTracksIndex[j],fStuff->JpsiRecTracksIndex[i]));
		} else {
		  fStuff->JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
						   (fStuff->JpsiRecTracks[i],fStuff->JpsiRecTracks[j]));
		  fStuff->JpsiPairTracksIndex.push_back(pair<int, int>
							(fStuff->JpsiRecTracksIndex[i],fStuff->JpsiRecTracksIndex[j]));
		}
		
		fStuff->KaonTrack.push_back(fStuff->KaonRecTracks[k]);
		fStuff->KaonTrackIndex.push_back(fStuff->KaonRecTracksIndex[k]);
	      }
	    }
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
void Bs2MuMu::muonCandTracks(const edm::Event &iEvent, const edm::EventSetup& iSetup, double m_cand1, double m_cand2) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>muonCandTracks> Sort opposite charge muons into pairs depending on inv. mass ("
		     << m_cand1 << " GeV) or (" << m_cand2 << " GeV), event: " << fNevt << endl;

  int type(0);

  if ( fStuff->MuonRecTracks.size() > 1 ) {

    for ( unsigned int i=0; i<fStuff->MuonRecTracks.size(); i++ ) { 

      for ( unsigned int j=i; j<fStuff->MuonRecTracks.size(); j++ ) {  // < ---- ??????

	if (i == j) { continue; }

	fStuff->RecTracks.clear();
	fStuff->RecTracksIndex.clear();

	fStuff->RecTracks.push_back(fStuff->MuonRecTracks[i]);
	fStuff->RecTracks.push_back(fStuff->MuonRecTracks[j]);

	fStuff->RecTracksIndex.push_back(fStuff->MuonRecTracksIndex[i]);
	fStuff->RecTracksIndex.push_back(fStuff->MuonRecTracksIndex[j]);

	type = massMuonCand(iEvent, iSetup, m_cand1, m_cand2);   //  < ------- here

	if ( type == 1 ) {

	  if ((*fStuff->MuonRecTracks[j]).pt() > (*fStuff->MuonRecTracks[i]).pt() ) { 
	    fStuff->BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (fStuff->MuonRecTracks[j],fStuff->MuonRecTracks[i]));
	    fStuff->BmmPairTracksIndex.push_back(pair<int, int>
						 (fStuff->MuonRecTracksIndex[j],fStuff->MuonRecTracksIndex[i]));
	  } else {
	    fStuff->BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (fStuff->MuonRecTracks[i],fStuff->MuonRecTracks[j]));
	    fStuff->BmmPairTracksIndex.push_back(pair<int, int>
						 (fStuff->MuonRecTracksIndex[i],fStuff->MuonRecTracksIndex[j]));
	  }
	}

	    
	if ( type == 2 ) {
	      
	  if ((*fStuff->MuonRecTracks[j]).pt() > (*fStuff->MuonRecTracks[i]).pt() ) { 
	    fStuff->JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					     (fStuff->MuonRecTracks[j],fStuff->MuonRecTracks[i]));
	    fStuff->JpsiPairTracksIndex.push_back(pair<int, int>
						  (fStuff->MuonRecTracksIndex[j],fStuff->MuonRecTracksIndex[i]));
	  } else {
	    fStuff->JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					     (fStuff->MuonRecTracks[i],fStuff->MuonRecTracks[j]));
	    fStuff->JpsiPairTracksIndex.push_back(pair<int, int>
						  (fStuff->MuonRecTracksIndex[i],fStuff->MuonRecTracksIndex[j]));
	  }
	}
      }
    }
  }
}


// ----------------------------------------------------------------------
void Bs2MuMu::kaonCandTracks(const edm::Event &iEvent, const edm::EventSetup& iSetup, double cone) { 
   
  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>kaonCandTracks> Find " << fStuff->JpsiPairTracks.size() 
		     << " kaon candidate(s) with best chi2, from all tracks in cone "
		     << cone << " around reconstructed J/Psi , event: " << fNevt << endl;

  const reco::Track* tt = 0;

  int index(0), trk(-1), ncand(0), cand(-1);
  double chi2_min(99999.), chi2(-1.);

  if ( fStuff->JpsiPairTracks.size() > 0 ) {
    
    for ( unsigned int i=0; i<fStuff->JpsiPairTracks.size(); i++ ) {

      // -- Muons from J/Psi
      fStuff->RecTracks.clear();
      fStuff->RecTracksIndex.clear();

      fStuff->RecTracks.push_back(fStuff->JpsiPairTracks[i].first);
      fStuff->RecTracksIndex.push_back(fStuff->JpsiPairTracksIndex[i].first);

      fStuff->RecTracks.push_back(fStuff->JpsiPairTracks[i].second); 
      fStuff->RecTracksIndex.push_back(fStuff->JpsiPairTracksIndex[i].second);  

      // ------------------------------------------------------------------------
      // -- Find kaon for given muon pair (if no kaon cand. found, index remains -1)
      fStuff->KaonTrack.push_back(tt);
      fStuff->KaonTrackIndex.push_back(-1);

      index = 0; trk = -1; ncand = 0; cand = -1;

      fStuff->SecVtxChi2.clear();
      fStuff->InvMass.clear();

      for (TrackCollection::const_iterator it = (*fStuff->theTkCollection).begin(); 
	   it != (*fStuff->theTkCollection).end(); 
	   ++it){
	
	if (index == fStuff->RecTracksIndex[0]) { index++;  continue; }
	if (index == fStuff->RecTracksIndex[1]) { index++;  continue; }
	
	tt = &(*it);
	
	fStuff->RecTracks.push_back(tt);
	fStuff->RecTracksIndex.push_back(index);
     
	chi2 =  rmmKaonCand(iEvent, iSetup, cone);   //  < ------- here: chi2 = -1, if outside cone!

	fStuff->RecTracks.pop_back();
	fStuff->RecTracksIndex.pop_back();
	
	if (chi2 > 0) {

	  if (chi2 < chi2_min) {
	    
	    chi2_min  = chi2;
	    
	    fStuff->KaonTrack.pop_back();
	    fStuff->KaonTrackIndex.pop_back();
	    
	    fStuff->KaonTrack.push_back(tt);
	    fStuff->KaonTrackIndex.push_back(index);
	    
	    trk  = index;
	    cand = ncand;
	  }

	  ncand++;
	}
	
	index++;	
      }

      // ------------------------------------------------------------------------  

      if (ncand)  fK200->Fill(fStuff->SecVtxChi2[cand], fStuff->InvMass[cand]);
      
      if (fVerbose) {

	if (ncand) { 
	  
	  cout << "==>kaonCandTracks> Fount " << ncand << " kaon candidate tracks for J/Psi candidate #"
	       << i << ". Selected track #" << trk << ", with chi2 = " << chi2_min << "   (chi2=" 
	       << fStuff->SecVtxChi2[cand] << ", m=" << fStuff->InvMass[cand] << ")." << endl; 

	} else {	
	  
	  cout << "==>kaonCandTracks> No kaon candidates found for J/Psi candidate #" 
	       << i << "." << endl;
	}
      }
 
      // -- eff. histogram ---------------------------

      if ( ncand < 9 ) {
	
	fEff->Fill(80.1 + ncand );
	
      } else {
	
	fEff->Fill(89.1);
      }

      // ----------------------------------------------  
 
    }
  }
}
// ----------------------------------------------------------------------
int Bs2MuMu::massMuonCand(const edm::Event &iEvent, const edm::EventSetup& iSetup, double m_cand1, double m_cand2) {  

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>massMuonCand> Looking for two muons in mass-windows around "
		     << m_cand1 << " and " << m_cand2 << ", event: " << fNevt << endl;

  TLorentzVector m0, m1, cand; 
  double mass(-1.);
  double dmass_1(0.6);
  double dmass_2(0.2);
  
  m0.SetXYZM((*fStuff->RecTracks[0]).px(),
	     (*fStuff->RecTracks[0]).py(),
	     (*fStuff->RecTracks[0]).pz(),
	     0.1056583);
      
  m1.SetXYZM((*fStuff->RecTracks[1]).px(),
	     (*fStuff->RecTracks[1]).py(),
	     (*fStuff->RecTracks[1]).pz(),
	     0.1056583);
    
  cand  = m0 + m1;
  mass  = cand.M();

  if ( mass > (m_cand1 - dmass_1)  && mass < (m_cand1 + dmass_1) &&
       fStuff->RecTracks[0]->charge() != fStuff->RecTracks[1]->charge() ) { 


    if (fVerbose) cout << endl << "==>massCandidate> : Found B candidates with muons \t #"
		       << fStuff->RecTracksIndex[0] << " \t #" <<  fStuff->RecTracksIndex[1]
		       << "\t --->  Inv. Mass: "  << mass
		       << "( q1 = " << fStuff->RecTracks[0]->charge()
		       << ", q2 = " << fStuff->RecTracks[1]->charge()
		       << ")" << endl;
    return 1;
	  
  } else if ( mass > (m_cand2 - dmass_2)  && mass < (m_cand2 + dmass_2) &&
	      fStuff->RecTracks[0]->charge() != fStuff->RecTracks[1]->charge() ) { 

    if (fVerbose) cout << endl << "==>massCandidate> : Found J/Psi candidates with muons \t #"
		       << fStuff->RecTracksIndex[0] << " \t #" <<  fStuff->RecTracksIndex[1]
		       << "\t --->  Inv. Mass: "  << mass
		       << "( q1 = " << fStuff->RecTracks[0]->charge()
		       << ", q2 = " << fStuff->RecTracks[1]->charge()
		       << ")" << endl;
    return 2;

  } else {

    if (fVerbose) cout << endl << "==>massCandidate> : Muons not in any of the mass windows \t #"
		       << fStuff->RecTracksIndex[0] << " \t #" <<  fStuff->RecTracksIndex[1]
		       << "\t --->  Inv. Mass: "  << mass
		       << "( q1 = " << fStuff->RecTracks[0]->charge()
		       << ", q2 = " << fStuff->RecTracks[1]->charge()
		       << ")" << endl;

    return -1;
  }
}


// ----------------------------------------------------------------------
double Bs2MuMu::rmmKaonCand(const edm::Event &iEvent, const edm::EventSetup& iSetup, double cone) {  

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>rmmKaonCand> Looking for kaon in cone of "
		     << cone << " around J/Psi, event: " << fNevt << endl;

  double chi2 = kalmanVertexFit(iEvent, iSetup, 0, 3);  //  Vertex TYPE = 0, B+ -> mu mu K (all K-candidates)
  
  TLorentzVector m0, m1, kp;
  TLorentzVector jpsi, bs;
  
  m0.SetXYZM(fStuff->RefittedTracks[0].px(),
	     fStuff->RefittedTracks[0].py(),
	     fStuff->RefittedTracks[0].pz(),
	     0.1056583);
  
  m1.SetXYZM(fStuff->RefittedTracks[1].px(),
	     fStuff->RefittedTracks[1].py(),
	     fStuff->RefittedTracks[1].pz(),
	     0.1056583);
  
  kp.SetXYZM(fStuff->RefittedTracks[2].px(),
	     fStuff->RefittedTracks[2].py(),
	     fStuff->RefittedTracks[2].pz(),
	     0.493677);
  
  jpsi = m0 + m1;
  bs   = m0 + m1 + kp;
  
  TVector3 jpsi_plab = jpsi.Vect();
  TVector3 kp_plab   = kp.Vect();
  
  double dphi = jpsi_plab.DeltaPhi(kp_plab);
  double deta = kp_plab.Eta() - jpsi_plab.Eta();
  double dr   = TMath::Sqrt(dphi*dphi + deta*deta);
  
  if ( dr < cone ) {
    
    if (fVerbose) cout << "==>rmmKaonCand> Fount kaon candidate, corresponding to track #" 
		       << fStuff->RecTracksIndex[2] 	   << " with dr = " << dr << ", B mass " 
		       << bs.M() << " and chi2 " << chi2 << endl;
    
    
    fStuff->InvMass.push_back(bs.M());
    fStuff->SecVtxChi2.push_back(chi2);

    fK100->Fill(chi2, bs.M());

    return chi2;
  
  } else {
    
    if (fVerbose) cout << "==>rmmKaonCand> !!! Kaon candidate rejected, ouside cone:  dr = " << dr
		       << ", B mass " << bs.M() << " and chi2 " << chi2 << " !!!" <<endl;

    return -1.;
  }
}

// ----------------------------------------------------------------------
void Bs2MuMu::trimBmmTracks(const edm::Event &iEvent) {  

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>trimBmmTracks> Trimming global muon vector of size " << fStuff->MuonRecTracks.size() 
       << ", event: " << fNevt << endl;

  TAnaTrack *pTrack;

  // -- get the two highest-pT muons 
  int m0(-1), m1(-1);
  double ptMu(0.), pT0(0.), pT1(0.);
  for (unsigned int i = 0; i < fStuff->MuonRecTracks.size(); ++i) {
    
    ptMu = (*fStuff->MuonRecTracks[i]).pt();
    if (fVerbose) cout << "  " << i << " with pT = " << ptMu << endl;
    if (ptMu > pT0) {
      if (m0 > -1) {
 	m1 = m0;
 	pT1 = pT0;
      }
      
      pT0 = ptMu;
      m0 = i;
    }  else {
      if (ptMu > pT1) {
 	pT1 = ptMu;
 	m1 = i;
      }	  
    }
  }
  

  if ((m0 > -1) && (m1 > -1)) {
    
    // -- RecTrack vector
    const reco::Track* tM0 = fStuff->MuonRecTracks[m0];
    const reco::Track* tM1 = fStuff->MuonRecTracks[m1];
    
    fStuff->RecTracks.clear();
    
    fStuff->RecTracks.push_back(tM0);
    fStuff->RecTracks.push_back(tM1);
    
    // -- Index into bmm-track vector 
    m0 = fStuff->MuonRecTracksIndex[m0];
    m1 = fStuff->MuonRecTracksIndex[m1];
    
    fStuff->RecTracksIndex.clear();
    
    fStuff->RecTracksIndex.push_back(m0);
    fStuff->RecTracksIndex.push_back(m1);

    for (unsigned int i = 0; i < fStuff->RecTracks.size(); ++i) {


      // -- Add to signal tracks block in ntuple
      reco::Track tt(*fStuff->RecTracks[i]);

      if (fVerbose) cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-" << endl;
      if (fVerbose) cout << "==>trimBmmTracks>  Adding to signal tracks, track # " << fStuff->RecTracksIndex[i]
			 << " (index in RecTracks), pT = " << tt.pt()
			 << endl;

      pTrack = fEvent->addSigTrack();
      pTrack->fIndex = fStuff->RecTracksIndex[i];

      fillTrack(iEvent, pTrack, &tt, fStuff->RecTracksIndex[i], 1);
    }

  } else {

    if (fVerbose) cout << "========> Not enough tracks for signal block" << endl;  
      
    for (unsigned int i = 0; i < fStuff->RecTracks.size(); ++i) {

      const reco::Track tt(*fStuff->RecTracks[i]);

      if (fVerbose) cout << " NOT added to signal tracks, track #" << fStuff->RecTracksIndex[i]
			 << " (index in RecTracks), pT = " << tt.pt()
			 << endl;
    }
  }
}


// ----------------------------------------------------------------------
int Bs2MuMu::primaryVertex(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>primaryVertex> List of primary vertices, event: " << fNevt << endl;

  int nvtx(0);
  const reco::Vertex* pV = 0;

  edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel("offlinePrimaryVerticesFromCTFTracks",recoPrimaryVertexCollection);
  const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());

  for(reco::VertexCollection::const_iterator v=recoPrimaryVertexCollection->begin(); 
      v!=recoPrimaryVertexCollection->end(); 
      ++v){


    nvtx++;

    printf ("%i. Primary Vertex (x, y, z) = ( %5.4f, %5.4f, %5.4f)\n", nvtx, v->x(), v->y(), v->z());

    pV = &(*v);   // ??? This takes the last Vertex - Selection ???

//     pV.SetX(v->x()*10);   //*10 to get mm (same unit as gen info)
//     pV.SetY(v->y()*10);
//     pV.SetZ(v->z()*10);
    
  }

  if (fVerbose) cout << endl << "====> Primary Vertices from CTF Tracks: " << nvtx << " vertices found." << endl;

  if ( nvtx < 9 )  { fEff->Fill(20.1 + nvtx); } else { fEff->Fill(29.1); }


  if (nvtx > 0) {

    //    printf ("Taking Primary Vertex (x, y, z) = ( %5.4f, %5.4f, %5.4f)\n", pV->x(), pV->y(), pV->z());
    printf ("Taking Primary Vertex (x, y, z) = ( %5.4f, %5.4f, %5.4f)\n", 10*pV->x(), 10*pV->y(), 10*pV->z());

    ChiSquared chi2(pV->chi2(),pV->ndof());

    TAnaVertex *pVtx;
    pVtx = new TAnaVertex();

    //    pVtx->setInfo(chi2.value(), chi2.degreesOfFreedom(), chi2.probability(), 0, 2); ??? change ndof to double ???
    pVtx->setInfo(chi2.value(), int(chi2.degreesOfFreedom()), chi2.probability(), 0, 2);
    pVtx->fPoint.SetXYZ(pV->position().x(),
		        pV->position().y(),
		        pV->position().z());

    fStuff->primaryVertex2 = *pV;
    fEvent->fPrimaryVertex2 = *pVtx;

  } else {
    
    TAnaVertex *pVtx;
    pVtx = new TAnaVertex();

    //    pVtx->setInfo(chi2.value(), chi2.degreesOfFreedom(), chi2.probability(), 0, 2); ??? change ndof to double ???
    pVtx->setInfo(-1., -1, -1., -1, -1);
    pVtx->fPoint.SetXYZ(-100., -100., -100.);

    fEvent->fPrimaryVertex2 = *pVtx;
  }

  return nvtx;
}

// ----------------------------------------------------------------------
void Bs2MuMu::secondaryVertex(const edm::Event &iEvent, const edm::EventSetup& iSetup) { 
  
  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>secondaryVertex> Fitting secondary vertices, event: " << fNevt << endl;

  if (fVerbose) {

    cout << "==>secondaryVertex> " << fStuff->BmmPairTracks.size() << " Bs -> mu+ mu- pair(s), event: " 
	 << fNevt << ", sel = " << fStuff->fBmmSel << endl; 
    cout << "==>secondaryVertex> " << fStuff->JpsiPairTracks.size() << " J/Psi -> mu+ mu- pair(s), event: " 
	 << fNevt << ", sel = " << fStuff->fBmmSel << endl; 
    cout << "==>secondaryVertex> " << fStuff->KaonTrack.size() << " Kaon(s) for B+ -> J/Psi K+, event: " 
	 << fNevt << ", sel = " << fStuff->fBmmSel << endl; 
  }

  double chi2(-1.);

  int nbs0  = fStuff->BmmPairTracks.size();
  int njpsi = fStuff->JpsiPairTracks.size();
  int nkaon = fStuff->KaonTrack.size();

  if ( fStuff->BmmPairTracks.size() > 0 ) {
    
    for ( unsigned int i=0; i<fStuff->BmmPairTracks.size(); i++ ) { 
      
      fStuff->RecTracks.clear();
      fStuff->RecTracksIndex.clear();

      fStuff->RecTracks.push_back(fStuff->BmmPairTracks[i].first);
      fStuff->RecTracks.push_back(fStuff->BmmPairTracks[i].second);

      fStuff->RecTracksIndex.push_back(fStuff->BmmPairTracksIndex[i].first);
      fStuff->RecTracksIndex.push_back(fStuff->BmmPairTracksIndex[i].second);

      chi2 = kalmanVertexFit(iEvent, iSetup, 5310+fStuff->fBmmSel, 2);  // Vertex TYPE = 531x  
                                                                        // -> mu mu or BG from B, 2 TRACKS => Cand 1
    }
  }
 
  if ( fStuff->JpsiPairTracks.size() > 0 ) {
    
    if ( fStuff->JpsiPairTracks.size() == fStuff->KaonTrack.size() ) {

      for ( unsigned int i=0; i<fStuff->JpsiPairTracks.size(); i++ ) { 
	
	if ( fStuff->KaonTrackIndex[i] > 0 ) {
	  
	  fStuff->RecTracks.clear();
	  fStuff->RecTracksIndex.clear();
	  
	  fStuff->RecTracks.push_back(fStuff->JpsiPairTracks[i].first);
	  fStuff->RecTracksIndex.push_back(fStuff->JpsiPairTracksIndex[i].first);
	  
	  fStuff->RecTracks.push_back(fStuff->JpsiPairTracks[i].second);
	  fStuff->RecTracksIndex.push_back(fStuff->JpsiPairTracksIndex[i].second);
	  
	  chi2 = kalmanVertexFit(iEvent, iSetup, 4430+fStuff->fBmmSel , 2);  // Vertex TYPE = 443x 
	                                                                     // -> mu mu from J/Psi, 2 TRACKS => Cand 2
	  
	  fStuff->RecTracks.push_back(fStuff->KaonTrack[i]);
	  fStuff->RecTracksIndex.push_back(fStuff->KaonTrackIndex[i]);
	  
	  chi2 = kalmanVertexFit(iEvent, iSetup, 5210+fStuff->fBmmSel, 3);  // Vertex TYPE = 521x  
	                                                                    // -> mu mu K, 3 TRACKS => Cand 3
	}
      }
    } else {

      cout << "****ERROR: Number of kaon candidates: " << fStuff->KaonTrack.size()
	   << ", is not equal to number of J/Psi candidates: " <<  fStuff->JpsiPairTracks.size() << "!!!! ****" <<endl;
      nkaon = 9;
      njpsi = 9;
    }
  }

  if (nbs0  >= 0 && nbs0  < 10) fEff->Fill(100.1 + (fStuff->fBmmSel-1)*100 + nbs0);
  if (njpsi >= 0 && njpsi < 10) fEff->Fill(120.1 + (fStuff->fBmmSel-1)*100 + njpsi);
  if (nkaon >= 0 && nkaon < 10) fEff->Fill(140.1 + (fStuff->fBmmSel-1)*100 + nkaon);
}


// ----------------------------------------------------------------------
double Bs2MuMu::kalmanVertexFit(const edm::Event &iEvent, const edm::EventSetup& iSetup
				, int type, unsigned int ntracks) {

  // TYPE =    -1: fillVertex without filling vertex/candidate in ntuple
  //     
  //      =     0: B->mu mu K (all candidates)
  //      =  531x: B->mu mu                    ----> where x = bmmsel mode
  //      =  443x: J/Psi->mu mu
  //      =  521x: B->mu mu K

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>kalmanVertexFit> Starting with secondary vertex, event: " << fNevt << endl;

  // -- MC truth vertices
//   Handle<SimVertexContainer> simVertexCollection;
//   iEvent.getByLabel("g4SimHits", simVertexCollection);
//   const SimVertexContainer simVC = *(simVertexCollection.product());

  // -- Kalman Vertex Fit
  std::vector<reco::TransientTrack> RecoTransientTrack;
  RecoTransientTrack.clear();
  
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  if (fVerbose) cout << "==>kalmanVertexFit> Vertexing with " << fStuff->RecTracks.size() << " reco. tracks." << endl;
   
  //get transient track 
  for ( unsigned int i=0; i<fStuff->RecTracks.size(); i++ ) { 

    RecoTransientTrack.push_back( (theB->build( &(*fStuff->RecTracks[i]) )) );
  }

  //get secondary vertex 
  if (fVerbose) cout << "==>kalmanVertexFit> Vertexing with " << RecoTransientTrack.size() << " trans. tracks " << endl;
  KalmanVertexFitter theFitter(true);
  
  if (fVerbose) cout << "==>kalmanVertexFit> Starting fitter TransSecVtx" << endl;

  // -- Does not work for CMSSW_1_2_0
  TransientVertex TransSecVtx = theFitter.vertex(RecoTransientTrack); 
  
  if ( TransSecVtx.isValid() ) {
    cout << "==>kalmanVertexFit> KVF successful!" << endl; 

    if ( isnan(TransSecVtx.position().x()) || isnan(TransSecVtx.position().y()) || isnan(TransSecVtx.position().z()) ) {

      cout << "==>kalmanVertexFit> Something went wrong!" << endl;
      cout << "==>kalmanVertexFit> SecVtx nan - Aborting... !" << endl;
      return -2.;
    }
  } else {
    cout << "==>kalmanVertexFit> KVF failed!" << endl;
    cout << "==>kalmanVertexFit> Aborting... !" << endl;
    return -1.;
  }

  if (fVerbose) cout << "==>kalmanVertexFit> Filling vector SecVtx" << endl;
  
  TVector3 SecVtx; 

  SecVtx.SetX(TransSecVtx.position().x()*10);  //*10 to get mm (same unit as gen info)
  SecVtx.SetY(TransSecVtx.position().y()*10);
  SecVtx.SetZ(TransSecVtx.position().z()*10);
  
  printf ("RECO SecVtx (x, y, z) = ( %5.4f, %5.4f, %5.4f)\n", SecVtx.X(), SecVtx.Y(), SecVtx.Z());

  ChiSquared chi(TransSecVtx.totalChiSquared(), TransSecVtx.degreesOfFreedom());

  if (fVerbose) cout << "Chi2 of SecVtx-Fit: " << chi.value() << endl;

  fillVertex(iEvent, iSetup, &TransSecVtx, type, ntracks);

  return chi.value();
}


// ----------------------------------------------------------------------
void Bs2MuMu::fillTrack(const edm::Event &iEvent, TAnaTrack *pTrack, reco::Track *it, int idx, int verb) {

  if ( verb == 1 ) {
    if (fVerbose) cout << "==>fillTrack> Filling track #" << idx << ", event: " << fNevt << endl;
  }

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  


  // -- Track parameters

  pTrack->fPlab.SetPtEtaPhi(it->pt(),
			    it->eta(),
			    it->phi()
			    );

  pTrack->fTip = it->d0();
  pTrack->fLip = it->dz();
  pTrack->fQ = it->charge();
  pTrack->fChi2 = it->chi2();
  //   pTrack->fDof = it->ndof();       ??? change ndof to double ???
  pTrack->fDof = int(it->ndof());
  pTrack->fHits = it->numberOfValidHits();  


  // -- Track matching to MC Track
  int index(-99999), type(-99999), mu(0);
  reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection); 

  std::vector<reco::Track> trackvector;
  //  TrackRef track = it;
  TrackRef track(tracks, idx);
  
  // -- trouble histo for reco2sim -----------------
  fEff->Fill(40.1);
  // -----------------------------------------------

  try{ 

    std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
    //       for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
    // 	     it != tp.end(); ++it) {

    TrackingParticleRef tpr = tp.begin()->first;  // ??? This takes the first associated simulated track ???	    
    double assocChi2 = tp.begin()->second;

    type = (*tpr).pdgId();
    index = fEvent->getGenIndex((*tpr).momentum().x(), (*tpr).momentum().y(), (*tpr).momentum().z(), type);     
    
    // -- Compare with CMSSW match
    const HepMC::GenParticle* genPar = 0; 
    int gen_pdg_id(-99999), gen_id(-99999); 
    int gen_cnt(0);

    for ( TrackingParticle::genp_iterator pit=tpr->genParticle_begin(); pit!=tpr->genParticle_end(); pit++ ){

      genPar=pit->get();

      gen_pdg_id = (*genPar).ParticleID();
      gen_id     = (*genPar).barcode()-1;

      gen_cnt++;

    }

    // -- trouble histo for sim2gen -----------------
    fEff->Fill(0.1);
    if ( abs(type) != abs(gen_pdg_id) ) {
      fEff->Fill(1.1);
    }
    if ( gen_cnt > 1 ) {
      fEff->Fill(2.1);
    }
    if ( gen_cnt == 0 ) {
      fEff->Fill(3.1);
    }
    if ( index < 0 ) {
      fEff->Fill(4.1);
    }
    if ( index != gen_id ) {
      fEff->Fill(5.1);
    }
    // -- trouble histo for reco2sim -----------------
    fEff->Fill(41.1);
    // -----------------------------------------------

    if ( (gen_cnt > 0) && (abs(type) == abs(gen_pdg_id)) ) {
      index = gen_id;
    // -- trouble histo for sim2gen -----------------
      fEff->Fill(6.1);
    }
    else {
      fEff->Fill(7.1);
    // -----------------------------------------------
    }

    if ( (index > 0) && (verb == 1) ) {
      
      if (fVerbose) cout << "--> Track #" << idx 
	   << " (pT = " << it->pt()
	   << ", eta = " << it->eta() << ", " << it->charge() << ")"
	   << " matched to Gen. Part. #" << index 
	   << " (" << tpr->pt()
	   << ", " << tpr->eta() << ")"
	   << ", PDG: " << type << ", NShared: " << assocChi2 << endl;
    }
    else if (verb == 1) {
      
      if (fVerbose) cout << "\t%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
    }
  } catch (Exception event) {

    if (verb == 1) {

      if (fVerbose) cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	   << setprecision(2) << setw(6) << track->pt() 
	   <<  " matched to 0 MC Tracks" << endl;
    }

    // -- trouble histo for reco2sim ---------------
    fEff->Fill(42.1);
    // ----------------------------------------------
  }

   pTrack->fGenIndex = index;
   pTrack->fMCID     = type;


  // -- Track matching to Muon Track 
  int idrec(-99999);
  int mcnt(0);
  pTrack->fMuID = 0.0;

  const reco::Track* tt = 0;

  for (MuonCollection::const_iterator glbMuon = (*fStuff->theMuonCollection).begin(); 
       glbMuon != (*fStuff->theMuonCollection).end(); 
       ++glbMuon){

    TrackRef glbTrack = glbMuon->combinedMuon();
    tt = &(*glbTrack);

    idrec = idRecTrack(tt);

    if ( idrec == idx ) {

      mu = 1;
      pTrack->fMuID = 1.0;

      if (verb == 1) {
	if (fVerbose) cout << " --> Track #" << idx 
	     << " (pT = " << it->pt()
	     << ", eta = " << it->eta() << ", " << it->charge() << ")"
	     << " matched to Global Muon #" << mcnt
	     << " (" << tt->pt()
	     << ", " << tt->eta() 
	     << ", " << tt->charge()<< ")" << endl;
      } 
      break;
    }
    mcnt++;
  }

  // -- trouble histo for muon2reco ---------------
  if (TMath::Abs(type) == 13) {
    if (mu) {
      fEff->Fill(43.1);
    }
    else {
      fEff->Fill(44.1);
    }
  }
  else {
    if (!mu) {
      fEff->Fill(45.1);
    }
    else {
      fEff->Fill(46.1);
    }
  }
  // ----------------------------------------------
}


// ----------------------------------------------------------------------
void Bs2MuMu::fillVertex(const edm::Event &iEvent, const edm::EventSetup& iSetup
			 , TransientVertex *v, int type, unsigned int ntracks) {

  // TYPE =    -1: without filling vertex/candidate in ntuple
  //     
  //      =     0: B->mu mu K (all candidates)
  //      =  531x: B->mu mu                    ----> where x = bmmsel mode
  //      =  443x: J/Psi->mu mu
  //      =  521x: B->mu mu K

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>fillVertex> Filling vertex, event: " << fNevt << endl;
  
  using namespace edm;
  using namespace reco;

  std::vector<reco::TransientTrack> refTT = v->refittedTracks();
  std::vector<reco::TransientTrack> orgTT = v->originalTracks();

  std::vector<reco::Track> refitted;
  std::vector<reco::Track> original;

  fStuff->RefittedTracks.clear();

  for(vector<reco::TransientTrack>::const_iterator i = refTT.begin(); 
      i != refTT.end(); i++) {

    const Track & ftt = i->track();
    // if (fVerbose) cout << "mmmmmmm : " << ftt->innermomentum() << endl;
    refitted.push_back(ftt);
    fStuff->RefittedTracks.push_back(ftt);
  }

  for(vector<reco::TransientTrack>::const_iterator i = orgTT.begin(); 
      i != orgTT.end(); i++) {

    const Track & ftt = i->track();
    original.push_back(ftt);
  }
  
  TLorentzVector m0, m1, kp;
  TLorentzVector bs;

  int i1(0), i2(0);
  double mass(-1.);

  fEff->Fill(30.1);

  // ==============================================================================================
  // ----------------------------- REFITTED TRACKS ------------------------------------------------

  if ( refitted.size() == ntracks ) {

    fEff->Fill(31.1);

    if ( ntracks == 2 ) {
      
      if (fVerbose) cout << "--- Vertex with number of refitted tracks: " << refitted.size() << " of 2 ---" << endl;
      
      i1 = fStuff->RecTracksIndex[0];
      i2 = fStuff->RecTracksIndex[1];

      m0.SetXYZM(refitted[0].px(),
		 refitted[0].py(),
		 refitted[0].pz(),
		 0.1056583);
    
      m1.SetXYZM(refitted[1].px(),
		 refitted[1].py(),
		 refitted[1].pz(),
		 0.1056583);
  
      bs = m0 + m1;
   
      mass = bs.M();
    }
      
    if ( ntracks == 3 ) {
      
      if (fVerbose) cout << "--- Vertex with number of refitted tracks: " << refitted.size() << " of 3 ---" << endl;
      
      i1 = fStuff->RecTracksIndex[0];
      i2 = fStuff->RecTracksIndex[2];  // This should be the K+

      m0.SetXYZM(refitted[0].px(),
		 refitted[0].py(),
		 refitted[0].pz(),
		 0.1056583);
      
      m1.SetXYZM(refitted[1].px(),
		 refitted[1].py(),
		 refitted[1].pz(),
		 0.1056583);
      
      kp.SetXYZM(refitted[2].px(),
		 refitted[2].py(),
		 refitted[2].pz(),
		 0.493677);
      
      bs = m0 + m1 + kp;
      
      mass = bs.M();     
    }
  }
  
  // ----------------------------- ORIGINAL TRACKS --------------------------------------------------
  
  else if ( original.size() == ntracks ) {

    if (fVerbose) cout << "--- Not enough refitted tracks for this Vertex: " << refitted.size()  << " ---" << endl;
    if (fVerbose) cout << "--- Use original tracks instead ---" << endl;
    
    fEff->Fill(32.1);

    if ( ntracks == 2 ) {
      
      if (fVerbose) cout << "--- Vertex with number of original tracks: " << original.size() << " of 2 ---" << endl;
      
      i1 = fStuff->RecTracksIndex[0];
      i2 = fStuff->RecTracksIndex[1];
      
      m0.SetXYZM(original[0].px(),
		 original[0].py(),
		 original[0].pz(),
		 0.1056583);
      
      m1.SetXYZM(original[1].px(),
		 original[1].py(),
		 original[1].pz(),
		 0.1056583);
      
      bs = m0 + m1;

      mass = bs.M();
    }
      
    if ( ntracks == 3 ) {
      
      if (fVerbose) cout << "--- Vertex with number of original tracks: " << original.size() << "of 3 ---" << endl;
      
      i1 = fStuff->RecTracksIndex[0];
      i2 = fStuff->RecTracksIndex[2];  // This should be the K+

      m0.SetXYZM(original[0].px(),
		 original[0].py(),
		 original[0].pz(),
		 0.1056583);
      
      m1.SetXYZM(original[1].px(),
		 original[1].py(),
		 original[1].pz(),
		 0.1056583);
      
      kp.SetXYZM(original[2].px(),
		 original[2].py(),
		 original[2].pz(),
		 0.493677);
      
      bs = m0 + m1 + kp;

      mass = bs.M();
    }
  }

  // ----------------------------- NO TRACKS ?! --------------------------------------------------

  else {

    fEff->Fill(33.1);
    if (fVerbose) cout << "!!! Not enough tracks for this Vertex: " << refitted.size() 
		       << " refitted and " << original.size() << "original tracks !!!" << endl;
  }
  
  // ---------------------------------------------------------------------------------------------
  // ==============================================================================================


  if (type == 0)            fM000[fStuff->fBmmSel-1]->Fill(mass);
  if (int(type/10) == 531)  fM100[fStuff->fBmmSel-1]->Fill(mass);
  if (int(type/10) == 443)  fM200[fStuff->fBmmSel-1]->Fill(mass);
  if (int(type/10) == 521)  fM300[fStuff->fBmmSel-1]->Fill(mass);


  ChiSquared chi(v->totalChiSquared(), v->degreesOfFreedom());

  VertexDistanceXY axy;
  double dXY      = axy.distance(fStuff->primaryVertex2, *v).value();
  double dXYE     = axy.distance(fStuff->primaryVertex2, *v).error();
  double compXY   = axy.compatibility(fStuff->primaryVertex2, *v);

  VertexDistance3D a3d;
  double d3d      = a3d.distance(fStuff->primaryVertex2, *v).value();
  double d3dE     = a3d.distance(fStuff->primaryVertex2, *v).error();
  double comp3d   = a3d.compatibility(fStuff->primaryVertex2, *v);


  if (fVerbose) {

    cout << "   KVF Vertex/cand of type " << type << " points to signal tracks: " << i1 << " and " << i2 << endl;
    cout << "      mass: " << mass << endl;
    cout << "      PVF xy: " << dXY << " +/- " << dXYE << "   comp: " << compXY << endl;
    cout << "      PVF 3d: " << d3d << " +/- " << d3dE << "   comp: " << comp3d << endl;
  }

  if ( type > -1 ) {  


    // -- Adding vertex to ntuple ...
    TAnaVertex *pVtx = new TAnaVertex();
    
    //    pVtx->setInfo(chi2.value(), chi2.degreesOfFreedom(), chi2.probability(), 0, 2); ??? change ndof to double ???
    pVtx->setInfo(chi.value(), int(chi.degreesOfFreedom()), chi.probability(), 1, type);
    pVtx->fPoint.SetXYZ(v->position().x(), v->position().y(), v->position().z());
    
    pVtx->addTrack(i1);
    pVtx->addTrack(i2);
    
    pVtx->fDxy  = dXY; 
    pVtx->fDxyE = dXYE; 
    pVtx->fCxy  = compXY; 
    
    pVtx->fD3d  = d3d; 
    pVtx->fD3dE = d3dE; 
    pVtx->fC3d  = comp3d; 
    
    // -- Adding Candidate to ntuple.
    TAnaCand  *pCand = fEvent->addCand();
    
    pCand->fPlab = bs.Vect();
    pCand->fMass = bs.M();
    
    pCand->fSig1 = i1;  
    pCand->fSig2 = i2;  
    pCand->fType = type;
    
    pCand->fVtx  = *pVtx;    

    if (fVerbose) cout << "fillVertex> Vertex/cand Added to ntuple ... " << endl; 
  }    
  


  // -----------------------------------------------------------------------------------
  //** -- Vertex reco2sim

//**   reco::VertexRecoToSimCollection recSimCollVertex = (*fStuff->recSimCollectionVertex); 


  // -------------------------------not used---------------------------------------------------
//   reco::RecoToSimCollection recSimColl             = (*fStuff->recSimCollection); 

  // -- get the collection of TrackingVertices
//   edm::Handle<TrackingVertexCollection>  TVCollectionH ;
//   iEvent.getByLabel("trackingtruth","VertexTruth",TVCollectionH);
//   const TrackingVertexCollection tVC   = *(TVCollectionH.product());

  // -- get the collection of primary Vertices
//   edm::Handle<reco::VertexCollection>  primaryVertexH ;
//   iEvent.getByLabel("offlinePrimaryVerticesFromCTFTracks","",primaryVertexH);
//   const reco::VertexCollection primaryVertexCollection   = *(primaryVertexH.product());

//   reco::VertexRecoToSimCollection vR2S = associatorByTracks ->
//     associateRecoToSim(primaryVertexH,TVCollectionH,iEvent,recSimColl);
  // ------------------------------not used-----------------------------------------------------

//**   for (reco::VertexRecoToSimCollection::const_iterator iR2S = recSimCollVertex.begin();
//**        iR2S != recSimCollVertex.end(); ++iR2S) {

//**     math::XYZPoint recoPos = (iR2S->key)->position();

//**     double nreco = (iR2S->key)->tracksSize();

//**     std::vector<std::pair<TrackingVertexRef, double> > vVR = iR2S -> val;

//**     for (std::vector<std::pair<TrackingVertexRef, double> >::const_iterator
//** 	   iMatch = vVR.begin(); iMatch != vVR.end(); ++iMatch) {

//**       TrackingVertexRef trueV =  iMatch->first;
//**       HepLorentzVector simVec = (iMatch->first)->position();

//**       double ntrue = trueV->daughterTracks().size();

//**       math::XYZPoint simPos = math::XYZPoint(simVec.x(),simVec.y(),simVec.z());

//**       double qual  = iMatch->second;
//**       double xmiss = simPos.X() - recoPos.X();
//**       double ymiss = simPos.Y() - recoPos.Y();
//**       double zmiss = simPos.Z() - recoPos.Z();
//**       double rmiss = sqrt(xmiss*xmiss+ymiss*ymiss+zmiss*zmiss);

//**       if (fVerbose) cout << " Qual. = " << qual
//** 	   << " xmiss = " << xmiss
//** 	   << " ymiss = " << ymiss
//** 	   << " zmiss = " << zmiss
//** 	   << " rmiss = " << rmiss
//** 	   << endl;
//**     }
//**   }

//   if (simv) {
//     pVtx->fSimPoint.SetXYZ(simv->position().x(),
// 			   simv->position().y(),
// 			   simv->position().z()
// 			   );
//   } else {
//     pVtx->fSimPoint.SetXYZ(-99.,
// 			   -99.,
// 			   -99.
// 			   );
//   }

//   if (fVerbose) cout << "      sim position: pos  = " 
//        << " " << pVtx->fSimPoint.X()
//        << " " << pVtx->fSimPoint.Y()
//        << " " << pVtx->fSimPoint.Z()
//        << endl;//
//

}

// ----------------------------------------------------------------------
int Bs2MuMu::idRecTrack(const reco::Track *track) {

  int found(-1), index(0);

  double ept(0.2), ephi(0.01), eeta(0.01);
  double mdpt(9999.), mdphi(9999.), mdeta(9999.);
  double dpt(0.),  dphi(0.),  deta(0.);

  //  double dhits(0.), mdhits(9999.)

  double pt  = track->pt();
  double phi = track->phi();
  double eta = track->eta();

  for (TrackCollection::const_iterator it = (*fStuff->theTkCollection).begin(); 
       it != (*fStuff->theTkCollection).end(); 
       ++it){


    dpt  = fabs(pt - it->pt());
    dphi = phi - it->phi();
    while (dphi >= M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    dphi = fabs(dphi);
    deta = fabs(eta - it->eta());

    if ((dpt < mdpt)
	&& (dphi < mdphi)
	&& (deta < mdeta)
	) {
      mdpt = dpt;
      mdphi = dphi;
      mdeta = deta;
      found = index;
    }

    ++index;
  }

  // if (fVerbose) cout << mdpt << " " << mdphi << " " << mdeta << " " << found << endl;

  if ((mdpt < ept)
      && (mdphi < ephi)
      && (mdeta < eeta)
      ) {
    return found;
  } else {
    return -1;
  }

}

// ----------------------------------------------------------------------
void Bs2MuMu::decayChannel(const char *fileName) {

  fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = -1;
  fTruthMC_mom = 531;                 fTruthMC2_mom = -1;
  fTruthMC_gmo = -1;                  fTruthMC2_gmo = -1;
  fMass = 5.367;                      fMass2 = -1;

  fPrintChannel = "N/D";              fPrintChannel2 = "N/D";

  //  if ( !strcmp("bd2pi", fileName.Data()) ) {
  if ( !strcmp("Bs2KpKm", fileName) ) {

    fTruthMC_I = 321;  fTruthMC_II = -1;
    fTruthMC_mom = 531;
    fMass = 5.367;

    fPrintChannel = "Bs -> K+ K-";
    if (fVerbose) cout << "Selected decay mode: Bs (" << fTruthMC_mom << ") to K+ K- (" << fTruthMC_I << ")" << endl;

  } else if ( !strcmp("Bd2PiPi", fileName) ) {

    fTruthMC_I = 211;  fTruthMC_II = -1;
    fTruthMC_mom = 511;
    fMass = 5.279;

    fPrintChannel = "Bd -> pi+ pi-";
    if (fVerbose) cout << "Selected decay mode: Bd (" << fTruthMC_mom << ") to pi+ pi- (" << fTruthMC_I << ")" << endl;

  } else if ( !strcmp("Lb2PKm", fileName) ) {

    fTruthMC_I = 321;  fTruthMC_II = 2212;
    fTruthMC_mom = 5122; 
    fMass = 5.624;                 

    fPrintChannel = "Lb -> p+ K-";
    if (fVerbose) cout << "Selected decay mode: Lb (" << fTruthMC_mom << ") to p+ (" 	 << fTruthMC_II 
		       << ") K- (" << fTruthMC_I << ")" << endl;

  } else {
    
    fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = 321;
    fTruthMC_mom = 531;                 fTruthMC2_mom = 443;
    fTruthMC_gmo = -1;                  fTruthMC2_gmo = 521;
    fMass = 5.367;                      fMass2 = 3.096;

    fPrintChannel = "Bs -> mu+ mu-";
    fPrintChannel2 = "B+ -> J/Psi K+";

    if (fVerbose) cout << "Selected decay modes: Bs (" << fTruthMC_mom << ") to mu+ mu- (" << fTruthMC_I << ")" << endl
		       << "     and B+ (" << fTruthMC2_gmo << ") to J/Psi (" << fTruthMC2_mom << ") K+ (" << fTruthMC2
		       << ") where J/Psi to mu+ mu- (" << fTruthMC_I  << ")" << endl;
  }          

  if (fVerbose) cout << "   ---> Channels (" << fChannel << "): " << fPrintChannel << " " << fPrintChannel2 << endl;

//   } else if ( !strcmp("bsmumug", fileName.Data()) || 
// 	    !strcmp("bsmumup0", fileName.Data()) ||
// 	    !strcmp("bu3munu", fileName.Data()) ||
// 	    !strcmp("bc3munu", fileName.Data()) ||
// 	    !strcmp("bcjpmunu", fileName.Data()) ) {
    
//     fTruthMC = 13;  fTruthMC2 = 13;
//     fChannel = TString("2mu_rare");

}

// ----------------------------------------------------------------------
void Bs2MuMu::clearTracks() {

  fStuff->BmmRecTracks.clear();
  fStuff->BmmRecTracksIndex.clear();
  fStuff->BmmRecTracksB.clear();

  fStuff->JpsiRecTracks.clear();
  fStuff->JpsiRecTracksIndex.clear();
  fStuff->JpsiRecTracksB.clear();

  fStuff->KaonRecTracks.clear();
  fStuff->KaonRecTracksIndex.clear();
  fStuff->KaonRecTracksB.clear();

  fStuff->MuonRecTracks.clear();
  fStuff->MuonRecTracksIndex.clear(); 
}

// ----------------------------------------------------------------------
void Bs2MuMu::clearCandidateTracks() {


  fStuff->BmmPairTracks.clear();
  fStuff->JpsiPairTracks.clear();
  fStuff->KaonTrack.clear();

  fStuff->BmmPairTracksIndex.clear();
  fStuff->JpsiPairTracksIndex.clear();
  fStuff->KaonTrackIndex.clear();
}

// ==============================================================================
//                           PRINT OUT
// ==============================================================================

void Bs2MuMu::printGenTracks(const edm::Event &iEvent) {
  
  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>printGenTracks> Starting to print generator block, event: " << fNevt << endl;


  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel(fSourceLabel.c_str(), evt);  
  const HepMC::GenEvent *genEvent = evt->GetEvent();

  int gcnt = genEvent->particles_size();
  if (fVerbose) cout << "Counted " << gcnt  << " genTracks in generator block" << endl;
  gcnt = 0;
  
  for (HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin();
       p != genEvent->particles_end();
       ++p) {
      
    if (fVerbose) cout << "GenTrack " <<  gcnt
	 << "   Particle ID: " << (*p)->pdg_id()
	 << ", Status "        << (*p)->StatusCode()
	 << ", pT "            << (*p)->Momentum().perp()
	 << ", eta "           << (*p)->Momentum().eta() << endl;
    
    gcnt++;
  }
}

// ----------------------------------------------------------------------
void Bs2MuMu::printSimTracks(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>printSimTracks> Starting to print simulated tracks, event: " << fNevt << endl;


  int scnt = (*fStuff->theTPCollection).size();
  if (fVerbose) cout << "Counted " << scnt  << " simulated tracks." << endl;
  scnt = 0;

  const HepMC::GenParticle* genPar = 0; 
  int gen_pdg_id(-99999), gen_id(-99999); 
  int gen_cnt(0);

  for (TrackingParticleCollection::const_iterator tpr = (*fStuff->theTPCollection).begin(); 
       tpr != (*fStuff->theTPCollection).end(); 
       ++tpr){

    if (fVerbose) cout << "SimTrack #" <<  scnt
	 << ": pT "     <<  (*tpr).pt()
	 << ", eta "    <<  (*tpr).eta()
	 << ", charge " <<  (*tpr).charge()
	 << ", PDG ID " <<  (*tpr).pdgId() << endl;
    

    // -- Using fEvent->getGenIndex to find matching track in generator block
    int index = fEvent->getGenIndex((*tpr).momentum().x(), (*tpr).momentum().y(), (*tpr).momentum().z(), 
				    (*tpr).pdgId());
    if ( index < 0 ) {
      
      if (fVerbose) cout << "%%> Did not find a match in generator block (fEvent->getGenIndex)" << endl;
    }
    else {
      if (fVerbose) cout << "==> matches to generator block particle #" << index << " (fEvent->getGenIndex)" <<endl;
    }


    // -- Using CMSSW's genp_iterator to find matching track in generator block
    gen_cnt = 0; gen_pdg_id = -99999; gen_id = -99999;
    for ( TrackingParticle::genp_iterator pit=tpr->genParticle_begin(); pit!=tpr->genParticle_end(); pit++ ){

      genPar=pit->get();

      gen_pdg_id = (*genPar).ParticleID();
      gen_id     = (*genPar).barcode()-1;

      if ( gen_id == index ) {

	if (fVerbose) cout << "==> matches to generator block particle #" << gen_id
	     << " with pT = " << (*genPar).Momentum().perp() 
	     << " and (gen) PDG ID: " << gen_pdg_id << endl;
      }
      else {
	if (fVerbose) cout << "==> matches to generator block particle #" << gen_id << " (disagreement!)"
	     << " with pT = " << (*genPar).Momentum().perp() 
	     << " and (gen) PDG ID: " << gen_pdg_id << endl;
      }

      gen_cnt++;

    }
    

    if ( gen_cnt == 0 ) {
	
      if (fVerbose) cout << "%%> no match in gen. block for sim. particle #" << scnt << endl;
    }
    else if ( gen_cnt > 1 ) {

      if (fVerbose) cout << " ==> !!! More than one gen. particle for sim. particle #" << scnt << " !!!" << endl;
    }

    if (fVerbose) cout << endl;
    scnt++;
  }
}


// ----------------------------------------------------------------------
void Bs2MuMu::printRecTracks(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>printRecTracks> Starting to print reconstructed tracks, event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );

 // -- track association
  reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection); 

  int index(-99999), type(-99999);

  int rcnt = (*fStuff->theTkCollection).size();
  if (fVerbose) cout << "Counted " << rcnt  << " reconstructed tracks." << endl;
  rcnt = 0;  

  for (TrackCollection::const_iterator itTrack = (*fStuff->theTkCollection).begin(); 
       itTrack != (*fStuff->theTkCollection).end(); 
       ++itTrack){

    // -- Track parameters
    if (fVerbose) cout << "RecTrack #"  <<  rcnt
	 << "\t : momentum " <<  (*itTrack).momentum() 
	 << ", pT "        <<  (*itTrack).pt()
	 << ", eta "       <<  (*itTrack).eta() << " +- " << (*itTrack).etaError() 
	 << ", phi "       <<  (*itTrack).phi() << " +- " << (*itTrack).phiError() << endl
	 << "\t\t   chi2 " <<  (*itTrack).chi2()
	 << ", TIP "       <<  (*itTrack).d0() << " +- " << (*itTrack).d0Error()
	 << ", LIP "       <<  (*itTrack).dz() << " +- " << (*itTrack).dzError()
	 << ", charge "    <<  (*itTrack).charge() 
	 << ", ndof "      <<  (*itTrack).ndof()
	 << ", hits "      <<  (*itTrack).numberOfValidHits()  << endl;
    
    // -- Track matching to MC Track

    index = -99999;
    type  = -99999;
    TrackRef track(tracks, rcnt);
    
    try{ 

      std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
      
      for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	   it != tp.end(); ++it) {
	
	TrackingParticleRef tpr = it->first;
	double assocChi2 = it->second;

	type  = (*tpr).pdgId();
	index = fEvent->getGenIndex((*tpr).momentum().x(), (*tpr).momentum().y(), (*tpr).momentum().z(), type);

	if ( index > 0 ) {

	  if (fVerbose) cout << "\t ----> Found gen. MC Track #" << index << " with PDG ID " << type 
	       << " (NShared: " << assocChi2 << ")" << endl;
	}
	else {

	  if (fVerbose) cout << "\t%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
	}
      }   

    } catch (Exception event) {

      if (fVerbose) cout << "\t%%> no MC match found for rec. Track #" << rcnt << endl;
    }

    rcnt++;
  }    
}

// ----------------------------------------------------------------------
void Bs2MuMu::printMuonTracks(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>printMuonTracks> Starting to print tracks reconstructed as muons, event: " << fNevt << endl;


  int mcnt = (*fStuff->theMuonCollection).size();
  if (fVerbose) cout << "Counted " << mcnt  << " reconstructed tracks." << endl;
  mcnt = 0;  

  for (MuonCollection::const_iterator glbMuon = (*fStuff->theMuonCollection).begin(); 
       glbMuon != (*fStuff->theMuonCollection).end(); 
       ++glbMuon){

    TrackRef glbTrack = glbMuon->combinedMuon();

    if (fVerbose) cout << "MuonTrack #" <<  mcnt
	 << ": pT "     <<  (*glbTrack).pt()
	 << ", eta "    <<  (*glbTrack).eta()
	 << ", charge " <<  (*glbTrack).charge()
	 << " matched to track #" << idRecTrack(&*glbTrack)
	 << endl; 

    mcnt++;
  }
}

// ----------------------------------------------------------------------
void Bs2MuMu::printReco2Sim(const edm::Event &iEvent, const char *option) {
 

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>printReco2Sim> Associating recontructed tracks to simulated tracks, event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  &recTC = *(tracks.product());

  // -- get the collection of TrackingParticles 
  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
  iEvent.getByLabel("trackingtruth","TrackTruth",TPCollectionH);
  const TrackingParticleCollection tPC   = *(TPCollectionH.product()); 


  // -- perform association by hits
  if ( !strcmp(option,"hits") ) {

    if (fVerbose) cout << "-- Associator by hits --" << endl;  
    reco::RecoToSimCollection recSimColl = 
      associatorByHits->associateRecoToSim (tracks, TPCollectionH, &iEvent);
    
    double ptRes(1000.);
    for(TrackCollection::size_type i=0; i<recTC.size(); ++i) {
      TrackRef track(tracks, i);
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
	if (fVerbose) cout << "Reco Track " << setw(2) << track.index() << " pT: "  << setw(6) << track->pt() 
	     <<  " matched to " << tp.size() << " MC Tracks" << std::endl;

	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  double assocChi2 = it->second;
	  if (fVerbose) cout << "\t\tMCTrack " << setw(2) << tpr.index() << " pdgId: " << setw(2) << tpr->pdgId() 
	       << " pT: " << setw(6) << tpr->pt() 
	       << " NShared: " << assocChi2 << endl;
	  
 	  ptRes=track->pt()-sqrt(tpr->momentum().Perp2());
 	  if (fVerbose) cout << "  res = " << ptRes << endl;
 	  fPT310->Fill(ptRes);
	}
      } catch (Exception event) {

	if (fVerbose) cout << "%%> Rec. Track #" << setw(2) << track.index() << " pT: " 
	     << setprecision(2) << setw(6) << track->pt() 
	     <<  " matched to 0 MC Tracks" << endl;
      }
    }
  }
  
  // -- perform association by chi2
  else if ( !strcmp(option,"chi2") ) {

    if (fVerbose) cout << "-- Associator by chi2 --" << endl;  
    reco::RecoToSimCollection recSimColl = 
      associatorByChi2->associateRecoToSim (tracks, TPCollectionH, &iEvent );

    double ptRes(1000.);
    for(TrackCollection::size_type i=0; i<recTC.size(); ++i) {
      TrackRef track(tracks, i);
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
	if (fVerbose) cout << "Reco Track " << setw(2) << track.index() << " pT: "  << setw(6) << track->pt() 
	     <<  " matched to " << tp.size() << " MC Tracks" << std::endl;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  double assocChi2 = it->second;
	  if (fVerbose) cout << "\t\tMCTrack " << setw(2) << tpr.index() << " pdgId: " << setw(2) << tpr->pdgId() 
	       << " pT: " << setw(6) << tpr->pt() 
	       << " NShared: " << assocChi2 << endl;
	  
	  ptRes=track->pt()-sqrt(tpr->momentum().Perp2());
	  if (fVerbose) cout << "  res = " << ptRes << endl;
	  fPT320->Fill(ptRes);
	}
      } catch (Exception event) {

	if (fVerbose) cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	     << setprecision(2) << setw(6) << track->pt() 
	     <<  " matched to 0 MC Tracks" << endl;
      }
    }
  }
  else {

    if (fVerbose) cout << "-- " << fAssocLabel.c_str() << " --" << endl;  

    reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection);

    double ptRes(1000.);

    for(TrackCollection::size_type i=0; i<(*fStuff->theTkCollection).size(); ++i) {

      TrackRef track(tracks, i);
      
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
	if (fVerbose) cout << "Reco Track " << setw(2) << track.index() << " pT: "  << setw(6) << track->pt() 
	     <<  " matched to " << tp.size() << " MC Tracks" << std::endl;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  double assocChi2 = it->second;
	  if (fVerbose) cout << "\t\tMCTrack " << setw(2) << tpr.index() << " pdgId: " << setw(2) << tpr->pdgId() 
	       << " pT: " << setw(6) << tpr->pt() 
	       << " NShared: " << assocChi2 << endl;
	  
	  ptRes=track->pt()-sqrt(tpr->momentum().Perp2());
	  if (fVerbose) cout << "  res = " << ptRes << endl;
	  fPT300->Fill(ptRes);
	}
      } catch (Exception event) {

	if (fVerbose) cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	     << setprecision(2) << setw(6) << track->pt() 
	     <<  " matched to 0 MC Tracks" << endl;
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void Bs2MuMu::beginJob(const edm::EventSetup& setup) {  
 
   edm::ESHandle<MagneticField> theMF;
   setup.get<IdealMagneticFieldRecord>().get(theMF);

   edm::ESHandle<TrackAssociatorBase> theChiAssociator;
   setup.get<TrackAssociatorRecord>().get("TrackAssociatorByChi2",theChiAssociator);
   associatorByChi2 = (TrackAssociatorBase *) theChiAssociator.product();

   edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
   setup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
   associatorByHits = (TrackAssociatorBase *) theHitsAssociator.product();

//**    edm::ESHandle<VertexAssociatorBase> theTracksAssociator;
//**    setup.get<VertexAssociatorRecord>().get("VertexAssociatorByTracks",theTracksAssociator);
//**    associatorByTracks = (VertexAssociatorBase *) theTracksAssociator.product();

}

// ------------ method called once each job just after ending the event loop  ------------
void Bs2MuMu::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(Bs2MuMu);
