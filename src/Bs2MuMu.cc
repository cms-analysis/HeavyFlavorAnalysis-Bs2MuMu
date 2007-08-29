// -*- C++ -*-
//
// Package:    Bs2MuMu
// Class:      Bs2MuMu
// 
/**\class Bs2MuMu Bs2MuMu.cc analysis/Bs2MuMu/src/Bs2MuMu.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Urs Langenegger
//         Created:  Mon Oct 23 15:14:30 CEST 2006
// $Id: Bs2MuMu.cc,v 1.28 2007/07/13 15:43:02 eggel Exp $
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

  TH1D *fH1;

  reco::VertexCollection *theVtxCollection;

  const reco::MuonCollection       *theMuonCollection;
  const reco::TrackCollection      *theTkCollection;
  const TrackingParticleCollection *theTPCollection;

  reco::RecoToSimCollection        *recSimCollection;
  //  reco::VertexRecoToSimCollection  *recSimCollectionVertex;

  std::vector<const reco::Vertex*> BmmPVTracks;

  std::vector<const reco::Track*> BmmRecTracks; 
  std::vector<int> BmmRecTracksIndex;

  double fPtMin;
  int fBmmSel;

  reco::Vertex primaryVertex;
  reco::Vertex primaryVertex2;

  TVector3 secondaryVertex;
};



// ----------------------------------------------------------------------
struct candStuff {

  std::vector<const reco::Track*> BmmRecTracks; 
  std::vector<int> BmmRecTracksIndex;

  std::vector<int> parPdgId;
  std::vector<int> parIndex;

  std::vector<int> momPdgId;
  std::vector<int> momIndex;

  std::vector<int> gmoPdgId;
  std::vector<int> gmoIndex;

  std::vector<double> SecVtxChi2;
  std::vector<double> InvMass;
};


// ----------------------------------------------------------------------
// ======================================================================

Bs2MuMu::Bs2MuMu(const edm::ParameterSet& iConfig) {

  cout << endl << "===>> Bs2MuMu >>> ctor, instantiating histogramms, etc." << endl;

  // -- Setup "class variables"
  fStuff = new anaStuff;
  fCand  = new candStuff;

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
  fFile = new TFile(iConfig.getParameter<string>( "fileName" ).c_str(), "RECREATE");
  fTree = new TTree("T1","CMSSW Bs -> mu+mu- tree");
  fEvent = new TAna00Event(0);
  fTree->Branch("TAna00Event", "TAna00Event", &fEvent, 256000/8, 1);


  // -- Troubleshoot histogramm
  fEff = new TH1D("eff", "Efficiencies", 100, 0., 100. );

  // -- Invariant mass & control histograms
  fM100  = new TH1D("m100", "inv. Mass Cand1", 1000, 0., 10. );
  fM200  = new TH1D("m200", "inv. Mass Cand2", 1000, 0., 10. );

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

  fM100->Write(); 
  fM200->Write();

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

  fStuff->BmmRecTracks.clear();
  fStuff->BmmPVTracks.clear();
  fStuff->BmmRecTracksIndex.clear();

//   fStuff->primaryVertex  = 0;
//   fStuff->primaryVertex2 = 0;


   fCand->BmmRecTracksIndex.clear();
   fCand->BmmRecTracks.clear();

   fCand->SecVtxChi2.clear();
   fCand->InvMass.clear();

   fCand->parIndex.clear();
   fCand->parPdgId.clear();
   fCand->momIndex.clear();
   fCand->momPdgId.clear();
   fCand->gmoIndex.clear();
   fCand->gmoPdgId.clear();

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
      
    cout << "==> -- Track Associator by hits --" << endl;  
    fStuff->recSimCollection = new 
      reco::RecoToSimCollection(associatorByHits->associateRecoToSim(tracks, TPCollectionH, &iEvent));
  }
  if ( !strcmp( (fAssocLabel.c_str()), "TrackAssociatorByChi2") ) {
      
    cout << "==> -- Track Associator by chi2 --" << endl;  
    fStuff->recSimCollection = new   
      reco::RecoToSimCollection(associatorByChi2->associateRecoToSim(tracks, TPCollectionH, &iEvent));
  }
    

//**  cout << "==> -- Vertex Associator by tracks --" << endl; 
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

  // -- (Rec.) Signal Tracks
  if (fStuff->fBmmSel == 1) {
    bmmTracks1(iEvent);
  } else if (fStuff->fBmmSel == 2) {
    bmmTracks2(iEvent);
  } else if (fStuff->fBmmSel == 3) {
    bmmTracks3(iEvent);
  } else if (fStuff->fBmmSel == 100) {
    bmmkTracks1(iEvent, iSetup);
  } else if (fStuff->fBmmSel == 300) {
    bmmkTracks3(iEvent, iSetup);
  }

  if ( fStuff->fBmmSel/100. < 1 ) {
    fNtracks = 2;  fType = 1; // Vertex TYPE = 1  ( -> mu mu or BG from B), 2 TRACKS
  } else {
    fNtracks = 3;  fType = 3; // Vertex TYPE = 3  ( -> mu mu K ), 3 TRACKS
  }

  // -- Primary Vertex
  int npv = primaryVertex(iEvent);

  if ( npv && fStuff->BmmRecTracks.size() == fNtracks )  {

    bmmVertex(iEvent, iSetup, fType, fNtracks); 
  }
  else {

    cout << "Not enough tracks: " << fStuff->BmmRecTracks.size() << ", needed " << fNtracks
	 << ". And/Or no primary Vtx found: "  << npv << " vertices found."<< endl;
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

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>fillGeneratorBlock> Starting to fill generator block, event: " << fNevt << endl;

  TGenCand *pCand;

  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel(fSourceLabel.c_str(), evt);
  const HepMC::GenEvent *genEvent = evt->GetEvent();
  

  int gcnt = genEvent->particles_size();
  cout << "Counted " << gcnt  << " genTracks in generator block" << endl;
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
//        cout << " SCREWED INDEED --> 1: DAU1 "
//  	   << tmpDau1 << " != " << pCand->fDau1 << " !!!!!!!!!!" << endl;
       scr0++;
     }
     if ( pCand->fDau1 == tmpDau1 && pCand->fDau2 != tmpDau2 ) {
//        cout << " SCREWED INDEED --> 1: DAU2 "
//  	   << tmpDau2 << " != " << pCand->fDau2 << " !!!!!!!!!!" << endl;
       scr1++;
     }
     if ( pCand->fDau1 != tmpDau1 && pCand->fDau2 != tmpDau2 ) {
//        cout << " SCREWED INDEED --> 2: DAU1 "
//  	   << tmpDau1 << " != " << pCand->fDau1 << " and DAU2 "
//  	   << tmpDau2 << " != " << pCand->fDau2 << " !!!!!!!!!!" << endl;
       scr2++;
     }

    if (abs(pCand->fID) == 531) {

      cout << "Daughters: #" << (*p)->beginDaughters()->barcode() - 1
	   << " - #" << (*p)->endDaughters()->barcode() - 1 << endl;
      cout << "1. Muon (fDau1)    " << pCand->fDau1 << endl;
      cout << "2. Muon (fDau2)    " << pCand->fDau2 << endl;
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
    //       cout << "  ----> Found a b quark ..." << endl;
    //     }

    if (531 == aid) {
      cout << "  ----> Found a Bs -- event " << fNgen << " should be kept!! pT = " << (*p)->Momentum().perp() << endl;
    }
  }

  fEff->AddBinContent(10, gcnt);
  fEff->AddBinContent(11, scr0);
  fEff->AddBinContent(12, scr1);
  fEff->AddBinContent(13, scr2);

//   cout << "============================================================================================" << endl;
//   fEvent->dumpGenBlock();
//   cout << "============================================================================================" << endl;
//   cout << "==>Bs2MuMu> Done with generator block." << endl;

}


// ----------------------------------------------------------------------
void Bs2MuMu::fillRecTracks(const edm::Event &iEvent) {
  
  fNrec++;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>fillRecTracks> Starting to fill reconstructed tracks, event: " << fNevt << endl;

  TAnaTrack *pTrack;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product());

  int rcnt = tracks->size();
  cout << "Counted " << rcnt  << " reconstructed tracks." << endl;
  rcnt = 0;

  // -- get the collection of TrackingParticles 
  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
  iEvent.getByLabel("trackingtruth","TrackTruth",TPCollectionH);
  const TrackingParticleCollection tPC   = *(TPCollectionH.product()); 

  int scnt = TPCollectionH->size();
  cout << "Counted " << scnt  << " simulated tracks." << endl;
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

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>bmmTracks1> Matching muons, truth-matched to be decay products of Bs, event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product()); 
 
 // -- track association
  reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection); 

 // -- perform association by either hits or chi2
//   reco::RecoToSimCollection recSimColl; 
//   if ( !strcmp( (fAssocLabel.c_str()), "TrackAssociatorByHits") ) {

//     cout << "==>bmmTracks1> -- Associator by hits --" << endl;  
//     recSimColl = associatorByHits->associateRecoToSim (tracks, TPCollectionH, &iEvent);
//   }
//   if ( !strcmp( (fAssocLabel.c_str()), "TrackAssociatorByChi2") ) {

//     cout << "==>bmmTracks1> -- Associator by chi2 --" << endl;  
//     recSimColl = associatorByChi2->associateRecoToSim (tracks, TPCollectionH, &iEvent );
//   }

  const reco::Track* tt = 0;
  const HepMC::GenParticle* genMom = 0;
  const HepMC::GenParticle* genPar = 0; 

  int ncand(0);
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

      cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-." <<endl;
      
      gen_cnt = 0; gen_pdg_id = -99999; gen_id = -99999;

      for ( TrackingParticle::genp_iterator pit=tpr->genParticle_begin(); pit!=tpr->genParticle_end(); pit++ ){
	
	genPar=pit->get();
	
	cout << ".. Rec. Track #" << track.index() << " pT = " << track->pt() << endl;
	cout << ".. Sim. Particle #" << tpr.index() << " pT = " << tpr->pt()  
	     << " PDG ID: " << tpr->pdgId() << " (NShared: "  << assocChi2 << ")" << endl;
	
	gen_pdg_id = genPar->ParticleID();
	gen_id     = genPar->barcode()-1;
	
	cout << ".. Gen. Particle #" << gen_id << " pT = " << genPar->Momentum().perp() << " ---> PDG ID: " << gen_pdg_id << endl;
	
	
	if ( abs(gen_pdg_id) == 13 ) { 
	  
	  genMom = genPar->mother();
	  mom_pdg_id = genMom->ParticleID();
	  mom_id     = genMom->barcode()-1;
	  
	  
	  cout << endl << ".. and Mother Particle #" << mom_id << " pT = " << genMom->Momentum().perp() 
	               << " ---> PDG ID: " << mom_pdg_id;
	  
	  if ( abs(mom_pdg_id) == 531) { 
	    
	    ncand++;

	    tt = &(*track);
	    
	    fStuff->BmmRecTracks.push_back(tt);
	    fStuff->BmmRecTracksIndex.push_back(track.index());
	    
	    cout << "    *** Bs -> mu- mu+ (" << ncand << ") *** ";
	  }
	  
	  cout << endl; 
	}

	gen_cnt++;	 
      }

      if ( gen_cnt == 0 ) {
	
	cout << "%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
      }

      if ( gen_cnt > 1 ) {
	
	cout << " ==> !!! More than one gen. particle for sim. particle #" << tpr.index() << " !!!" << endl;
      }
    } catch (Exception event) {

      cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	   << setprecision(2) << setw(6) << track->pt() 
	   <<  " matched to 0 MC Tracks" << endl;
    } 
    
  } 

  trimBmmTracks(iEvent);
}


// ----------------------------------------------------------------------
void Bs2MuMu::bmmTracks2(const edm::Event &iEvent) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>bmmTracks2> Matching muons, truth matched to be muons in the generator block, event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product());

 // -- track association
  reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection); 

  const reco::Track* tt = 0;
  const HepMC::GenParticle* genPar = 0; 
  const HepMC::GenParticle* genMom = 0;

  int ncand(0);
  int gen_pdg_id(-99999), gen_id(-99999);
  int mom_pdg_id(-99999), mom_id(-99999);
  int gen_cnt(0);

  for(TrackCollection::size_type i=0; i<recTC.size(); ++i) {

    TrackRef track(tracks, i);

    try{ 

      std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
//       for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
// 	     it != tp.end(); ++it) {

      TrackingParticleRef tpr = tp.begin()->first;  // ??? This takes the first associated simulated track ???	    
      double assocChi2 = tp.begin()->second;

      cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-." <<endl;
      
      gen_cnt = 0; gen_pdg_id = -99999; gen_id = -99999;

      for ( TrackingParticle::genp_iterator pit=tpr->genParticle_begin(); pit!=tpr->genParticle_end(); pit++ ){
	
	genPar=pit->get();
		
	gen_pdg_id = genPar->ParticleID();
	gen_id     = genPar->barcode()-1;

	if ( abs(gen_pdg_id) == 13 ) { 
	  
	  ncand++;
	  tt = &(*track);

	  cout << ".. Rec. Track #" << track.index() << " pT = " << track->pt() << endl;
	  cout << ".. Sim. Particle #" << tpr.index() << " pT = " << tpr->pt() 
	       << " PDG ID: " << tpr->pdgId() << " (NShared: "  << assocChi2 << ")" << endl;
	  cout << ".. Gen. Particle #" << gen_id << " pT = " << genPar->Momentum().perp() 
	       << " ---> PDG ID: " << gen_pdg_id << endl;
	  
	  
	  fStuff->BmmRecTracks.push_back(tt);
	  fStuff->BmmRecTracksIndex.push_back(track.index());

	  genMom = genPar->mother();
	  mom_pdg_id = genMom->ParticleID();
	  mom_id     = genMom->barcode()-1;
		  
	  cout << endl << ".. Mother Particle #" << mom_id << " pT = " << genMom->Momentum().perp() 
	               << " ---> PDG ID: " << mom_pdg_id;
	    
	  if ( abs(mom_pdg_id) == 13 ) { 
	    
	    cout << "    *** Bs -> mu- mu+ (" << ncand << ") *** ";
	  }
	  
	  cout << endl; 
	}

	gen_cnt++;	 
      }

      if ( gen_cnt == 0 ) {
	
	cout << "%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
      }

      if ( gen_cnt > 1 ) {
	
	cout << " ==> !!! More than one gen. particle for sim. particle #" << tpr.index() << " !!!" << endl;
      }
    } catch (Exception event) {

      cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	   << setprecision(2) << setw(6) << track->pt() 
	   <<  " matched to 0 MC Tracks" << endl;
    }
  }

  trimBmmTracks(iEvent);
}


// ----------------------------------------------------------------------
void Bs2MuMu::bmmTracks3(const edm::Event &iEvent) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>bmmTracks3> Matching muons, rec PID as muons, event: " << fNevt << endl;

  if (0 == fStuff->theMuonCollection ) {
    
    cout << "==>bmmTrack3> getting the muon collection ..." << endl;
    
    edm::Handle<reco::MuonCollection> MuCollection;
    iEvent.getByLabel(fMuonLabel.c_str(), MuCollection);
    fStuff->theMuonCollection   = new reco::MuonCollection(*(MuCollection.product()));
  }

  int nMuon  = fStuff->theMuonCollection->size();
  int nTrack = fStuff->theTkCollection->size();

  cout << "==>bmmTracks3>  tracks found: " << nTrack << endl
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

      fStuff->BmmRecTracks.push_back(tt);
      fStuff->BmmRecTracksIndex.push_back(idrec);
    } else {
      cout << "==>bmmTracks3> Coulnd't find rec. track for muon #" << mcnt << endl;  
    }
    
    mcnt++;
//     glbMuons1.push_back(MuonRef(MuCollection,mcnt-1));
//     glbMuons2.push_back(MuonRef(MuCollection,mcnt-1));
  }

  trimBmmTracks(iEvent);
}


// ----------------------------------------------------------------------
void Bs2MuMu::bmmkTracks1(const edm::Event &iEvent, const edm::EventSetup& iSetup) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>bmmkTracks1> Matching J/Psi K, truth-matched to be decay products of B+, event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product()); 
 
 // -- track association
  reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection); 

  const reco::Track* tt = 0;
  const HepMC::GenParticle* genGmo = 0;
  const HepMC::GenParticle* genMom = 0;
  const HepMC::GenParticle* genPar = 0; 

  int mcand(0), kcand(0);
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

      cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-." <<endl;
      
      gen_cnt = 0; gen_pdg_id = -99999; gen_id = -99999;
      mom_pdg_id = -99999; mom_id = -99999;
      gmo_pdg_id = -99999; gmo_id = -99999;

      for ( TrackingParticle::genp_iterator pit=tpr->genParticle_begin(); pit!=tpr->genParticle_end(); pit++ ){
	
	genPar=pit->get();
	
	cout << ".. Rec. Track #" << track.index() << " pT = " << track->pt() << endl;
	cout << ".. Sim. Particle #" << tpr.index() << " pT = " << tpr->pt() 
	                                            << " (NShared: "  << assocChi2 << ")" << endl;
	
	gen_pdg_id = genPar->ParticleID();
	gen_id     = genPar->barcode()-1;
	
	cout << ".. Gen. Particle #" << gen_id << " pT = " << genPar->Momentum().perp() << " ---> PDG ID: " << gen_pdg_id << endl;
	
	if ( abs(gen_pdg_id) == 13 ) { 
	  
	  genMom = genPar->mother();
	  mom_pdg_id = genMom->ParticleID();
	  mom_id     = genMom->barcode()-1;
	  
	  
 	  cout << endl << ".. and Mother Particle #" << mom_id << " pT = " << genMom->Momentum().perp() 
	       << " ---> PDG ID: " << mom_pdg_id;
	  
	  if ( abs(mom_pdg_id) == 443) { 

	    genGmo = genMom->mother();
	    gmo_pdg_id = genGmo->ParticleID();
	    gmo_id     = genGmo->barcode()-1;

	    cout << endl << ".. and Grandmother Particle #" << gmo_id << " pT = " << genGmo->Momentum().perp() 
		 << " ---> PDG ID: " << gmo_pdg_id;
	    
	    if ( abs(gmo_pdg_id) == 521) { 
	    
	      mcand++;
	      
	      tt = &(*track);
	      
	      fCand->BmmRecTracks.push_back(tt);
	      fCand->BmmRecTracksIndex.push_back(track.index());

	      fCand->parIndex.push_back(gen_id);
	      fCand->parPdgId.push_back(gen_pdg_id);

	      fCand->momIndex.push_back(mom_id);
	      fCand->momPdgId.push_back(mom_pdg_id);

	      fCand->gmoIndex.push_back(gmo_id);
	      fCand->gmoPdgId.push_back(gmo_pdg_id);
	      
	      cout << "    *** B+ -> (J/Psi->mu-mu+) K+ (" << mcand << ". muon-cand.) *** ";
	    }
	  }
	  
	  cout << endl; 
	}

	if ( abs(gen_pdg_id) == 321 ) { 

	  genMom = genPar->mother();
	  mom_pdg_id = genMom->ParticleID();
	  mom_id     = genMom->barcode()-1;

  	  cout << endl << ".. and Mother Particle #" << mom_id << " pT = " << genMom->Momentum().perp() 
 	       << " ---> PDG ID: " << mom_pdg_id;
	  
	  if ( abs(mom_pdg_id) == 521) { 

	    gmo_pdg_id = mom_pdg_id;
	    gmo_id     = mom_id;

	    cout << endl << ".. and Grandmother Particle was set to #" << gmo_id  
		 << " ---> PDG ID: " << gmo_pdg_id;

	    kcand++;

	    tt = &(*track);
	    
	    fCand->BmmRecTracks.push_back(tt);
	    fCand->BmmRecTracksIndex.push_back(track.index());

	    fCand->parIndex.push_back(gen_id);
	    fCand->parPdgId.push_back(gen_pdg_id);

	    fCand->momIndex.push_back(mom_id);
	    fCand->momPdgId.push_back(mom_pdg_id);
	    
	    fCand->gmoIndex.push_back(gmo_id);
	    fCand->gmoPdgId.push_back(gmo_pdg_id);
	    
	    cout << "    *** B+ -> (J/Psi->mu-mu+) K+ (" << kcand << ". kaon-cand.) *** ";
	  }
	  
	  cout << endl; 
	}

	gen_cnt++;	 
      }

      if ( gen_cnt == 0 ) {
	
	cout << "%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
      }

      if ( gen_cnt > 1 ) {
	
	cout << " ==> !!! More than one gen. particle for sim. particle #" << tpr.index() << " !!!" << endl;
      }
    } catch (Exception event) {

      cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	   << setprecision(2) << setw(6) << track->pt() 
	   <<  " matched to 0 MC Tracks" << endl;
    } 
    
  }

  if ( kcand < 9 ) { fEff->Fill(60.1 + kcand); } else { fEff->Fill(69.1); }
  if ( mcand < 9 ) { fEff->Fill(70.1 + mcand); } else { fEff->Fill(79.1); }

  if ( kcand >= 1 && mcand >= 2 ) { 

    fEff->Fill(80.1); 

  } else { 

    fEff->Fill(81.1); 

  }

  // -- Finding the daughters from same B+
  int bs_index(-1), k_index(-1), dcnt(0);
  std::vector<int> mu_index;

  for ( unsigned int i = 0; i < fCand->parIndex.size(); i++) {

    mu_index.clear();
    bs_index = -1;
    k_index  = -1;
    dcnt = 0;

    if ( abs(fCand->parPdgId[i]) == 321 ) {

      bs_index = fCand->momIndex[i];
      k_index  = i;  // Kaon Index

      for (unsigned int j = 0; j < fCand->parIndex.size(); j++ ) {

	if ( (fCand->gmoIndex[j] == bs_index) && (i != j)  ) {
	 
	  if ( abs(fCand->parPdgId[j]) == 321 ) {

	    cout << "==> !!! Found another kaon from same B+ !!!" << endl;
	  }
	  if ( (abs(fCand->parPdgId[j]) == 321) && ((*fCand->BmmRecTracks[j]).pt() > (*fCand->BmmRecTracks[i]).pt()) ) {

	    cout << "==> !!! Found another kaon from same B+ with higher pT !!!" << endl;
	    k_index = j;  // Kaon Index
	  }
	  if ( abs(fCand->parPdgId[j]) == 13 ) {

	    mu_index.push_back(j);  // Muon Index
      	  }
	}
      }

      if ( mu_index.size() >= 2 ) {
	
	int id(-1);

	// -- Fill in muons --
	for (unsigned int k = 0; k < mu_index.size(); k++) {
	  
	  id = mu_index[k];
	  fStuff->BmmRecTracks.push_back(fCand->BmmRecTracks[id]);
	  fStuff->BmmRecTracksIndex.push_back(fCand->BmmRecTracksIndex[id]);
	}
	
	// ... trim and add 2 muons to signal tracks block in ntuple
	trimBmmTracks(iEvent);            // including fillTracks (adds 2 muon signal tracks)
	bmmVertex(iEvent, iSetup, 2, 2);  // including fillVertex (adds J/Psi candidates) 
	                                  // Vertex TYPE = 2 (-> mu mu from J/Psi)
	
	// -- Fill in kaon --
	id = k_index;
	fStuff->BmmRecTracks.push_back(fCand->BmmRecTracks[id]);
	fStuff->BmmRecTracksIndex.push_back(fCand->BmmRecTracksIndex[id]);  
	
	// ... and add to signal tracks block in ntuple
	reco::Track tt(*fCand->BmmRecTracks[id]);
	cout << " Adding to signal tracks, track #" << fCand->BmmRecTracksIndex[id]
	     << " (index in RecTracks), pT = " << tt.pt()
	     << endl;

	TAnaTrack *pTrack;
	pTrack = fEvent->addSigTrack();
	pTrack->fIndex = fCand->BmmRecTracksIndex[id];
	
	fillTrack(iEvent, pTrack, &tt, fCand->BmmRecTracksIndex[id], 1);
	
	break;
      }
    }
  }
}
// ----------------------------------------------------------------------
void Bs2MuMu::bmmkTracks3(const edm::Event &iEvent, const edm::EventSetup& iSetup) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>bmmkTracks3> Matching muons from J/Psi, rec PID as muons and searching for kaon, event: " << fNevt << endl;

  if (0 == fStuff->theMuonCollection ) {
    
    cout << "==>bmmkTracks3> getting the muon collection ..." << endl;
    
    edm::Handle<reco::MuonCollection> MuCollection;
    iEvent.getByLabel(fMuonLabel.c_str(), MuCollection);
    fStuff->theMuonCollection   = new reco::MuonCollection(*(MuCollection.product()));
  }

  int nMuon  = fStuff->theMuonCollection->size();
  int nTrack = fStuff->theTkCollection->size();

  cout << "==>bmmkTracks3>  tracks found: " << nTrack << endl
       << "==>bmmkTracks3>  muons  found: " << nMuon  << endl;  

  const reco::Track* tt = 0;

  int mcnt(0);
  int mcand(0), kcand(0);
  for (MuonCollection::const_iterator glbMuon = (*fStuff->theMuonCollection).begin(); 
       glbMuon != (*fStuff->theMuonCollection).end(); 
       ++glbMuon){

    TrackRef glbTrack = glbMuon->combinedMuon();
    tt = &(*glbTrack);

    int idrec = idRecTrack(tt);

    //    if ( (idrec > -1) && (tt->pt() > 3.) ) {
    if ( idrec > -1 ) {  // pT Cut removed !!!!!

      mcand++;
      fStuff->BmmRecTracks.push_back(tt);
      fStuff->BmmRecTracksIndex.push_back(idrec);
    } else {
      cout << "==>bmmkTracks3> Coulnd't find rec. track for muon #" << mcnt << endl;  
    }

    mcnt++;
  }
  
  // ... trim and add 2 muons to signal tracks block in ntuple
  if ( fStuff->BmmRecTracks.size() > 1 ) {
    if ( jpsiCandidate(iEvent) > 0 ) {      // two muons in J/Psi mass window
      
      trimBmmTracks(iEvent);                // including fillTracks (adds 2 muon signal tracks)
      bmmVertex(iEvent, iSetup, 2, 2);      // including fillVertex (adds J/Psi candidates)
                                            // Vertex TYPE = 2 (-> mu mu from J/Psi)  

      // ... find kaon and add to signal tracks block in ntuple
      int id = kaonCandidate(iEvent, iSetup);  // Vertex all candidate tracks, 
                                               // cut on dR & select best chi2 
      if ( id >= 0  ) {
	
	kcand++;
	fStuff->BmmRecTracks.push_back(fCand->BmmRecTracks[id]);
	fStuff->BmmRecTracksIndex.push_back(fCand->BmmRecTracksIndex[id]);  

	reco::Track kt(*fCand->BmmRecTracks[id]);
	
	cout << "==> Adding to signal tracks, track #" << fCand->BmmRecTracksIndex[id]
	     << " (index in RecTracks), pT = " << kt.pt()
	     << endl;
	
	TAnaTrack *pTrack;
	pTrack = fEvent->addSigTrack();
	pTrack->fIndex = fCand->BmmRecTracksIndex[id];
	
	fillTrack(iEvent, pTrack, &kt, fCand->BmmRecTracksIndex[id], 1);
      } else {
	cout << "==>bmmkTracks3> Aborting, couldn't find kaon !!!" << endl;  
      }
    }
    else {
      fStuff->BmmRecTracks.clear();
      fStuff->BmmRecTracksIndex.clear();
      cout << "==>bmmkTracks3> Aborting, couldn't find two muons in J/Psi mass window !!!" << endl;  
    }
  } else if ( mcnt == 0 ){
    cout << "==>bmmkTracks3> Aborting, no global muons !!!" << endl;  
  } else {
    cout << "==>bmmkTracks3> Aborting, not enough global muons !!!" << endl;
  }

  if ( kcand < 9 ) { fEff->Fill(60.1 + kcand); } else { fEff->Fill(69.1); }
  if ( mcand < 9 ) { fEff->Fill(70.1 + kcand); } else { fEff->Fill(79.1); }

  if ( kcand >= 1 && mcand >= 2 ) { 

    fEff->Fill(80.1); 

  } else { 

    fEff->Fill(81.1); 

  }
}

// ----------------------------------------------------------------------
void Bs2MuMu::trimBmmTracks(const edm::Event &iEvent) {  

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>trimBmmTracks> Trimming bmm-track vector of size " << fStuff->BmmRecTracks.size() 
       << ", event: " << fNevt << endl;

  TAnaTrack *pTrack;

  // -- get the two highest-pT muons 
  int m0(-1), m1(-1);
  double ptMu(0.), pT0(0.), pT1(0.);
  for (unsigned int i = 0; i < fStuff->BmmRecTracks.size(); ++i) {
    
    ptMu = (*fStuff->BmmRecTracks[i]).pt();
    cout << "  " << i << " with pT = " << ptMu << endl;
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
    const reco::Track* tM0 = fStuff->BmmRecTracks[m0];
    const reco::Track* tM1 = fStuff->BmmRecTracks[m1];
    
    fStuff->BmmRecTracks.clear();
    
    fStuff->BmmRecTracks.push_back(tM0);
    fStuff->BmmRecTracks.push_back(tM1);
    
    // -- Index into bmm-track vector 
    m0 = fStuff->BmmRecTracksIndex[m0];
    m1 = fStuff->BmmRecTracksIndex[m1];
    
    fStuff->BmmRecTracksIndex.clear();
    
    fStuff->BmmRecTracksIndex.push_back(m0);
    fStuff->BmmRecTracksIndex.push_back(m1);
    
    //    fSel = 1;
    
    for (unsigned int i = 0; i < fStuff->BmmRecTracks.size(); ++i) {


      // -- Add to signal tracks block in ntuple
      reco::Track tt(*fStuff->BmmRecTracks[i]);
      cout << "==>trimBmmTracks>  Adding to signal tracks, track # " << fStuff->BmmRecTracksIndex[i]
 	   << " (index in RecTracks), pT = " << tt.pt()
 	   << endl;

      pTrack = fEvent->addSigTrack();
      pTrack->fIndex = fStuff->BmmRecTracksIndex[i];

      fillTrack(iEvent, pTrack, &tt, fStuff->BmmRecTracksIndex[i], 1);
    }

    //    fNrec++;

  } else {
    cout << "========> Not enough tracks for signal block" << endl;  
      
    for (unsigned int i = 0; i < fStuff->BmmRecTracks.size(); ++i) {

      const reco::Track tt(*fStuff->BmmRecTracks[i]);

      cout << " NOT added to signal tracks, track #" << fStuff->BmmRecTracksIndex[i]
 	   << " (index in RecTracks), pT = " << tt.pt()
 	   << endl;
    }

    //    fSel = 0;
  }
}


// ----------------------------------------------------------------------
int Bs2MuMu::primaryVertex(const edm::Event &iEvent) {  


  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>primaryVertex> List of primary vertices, event: " << fNevt << endl;

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

  cout << endl << "====> Primary Vertices from CTF Tracks: " << nvtx << " vertices found." << endl;

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
  }

  return nvtx;
}


// ----------------------------------------------------------------------
void Bs2MuMu::bmmVertex(const edm::Event &iEvent, const edm::EventSetup& iSetup, int type, unsigned int ntracks) {


  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>bmmVertex> Starting with secondary vertex, event: " << fNevt << endl;

  // -- MC truth vertices
//   Handle<SimVertexContainer> simVertexCollection;
//   iEvent.getByLabel("g4SimHits", simVertexCollection);
//   const SimVertexContainer simVC = *(simVertexCollection.product());

  // -- Kalman Vertex Fit
  std::vector<reco::TransientTrack> RecoTransientTrack;
  RecoTransientTrack.clear();
  
  
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  cout << "==>bmmVertex> Vertexing with " << fStuff->BmmRecTracks.size() << " reco. tracks." << endl;
   
  //get transient track 
  for ( unsigned int i=0; i<fStuff->BmmRecTracks.size(); i++ ) { 

    RecoTransientTrack.push_back( (theB->build( &(*fStuff->BmmRecTracks[i]) )) );
  }

  //get secondary vertex 
  cout << "==>bmmVertex> Vertexing with " << RecoTransientTrack.size() << " trans. tracks " << endl;
  KalmanVertexFitter theFitter(true);
  
  cout << "==>bmmVertex> Starting fitter TransSecVtx" << endl;

  // -- Does not work for CMSSW_1_2_0
  TransientVertex TransSecVtx = theFitter.vertex(RecoTransientTrack); 
  
  if ( TransSecVtx.isValid() ) {
    cout << "==>bmmVertex> KVF successful!" << endl; 
  } else {
    cout << "==>bmmVertex> KVF failed!" << endl;
    cout << "==>bmmVertex> Aborting... !" << endl;
    return;
  }

  cout << "==>bmmVertex> Filling vector SecVtx" << endl;
  
  TVector3 SecVtx; 

  SecVtx.SetX(TransSecVtx.position().x()*10);  //*10 to get mm (same unit as gen info)
  SecVtx.SetY(TransSecVtx.position().y()*10);
  SecVtx.SetZ(TransSecVtx.position().z()*10);
  
  printf ("RECO SecVtx (x, y, z) = ( %5.4f, %5.4f, %5.4f)\n", SecVtx.X(), SecVtx.Y(), SecVtx.Z());

  ChiSquared chi(TransSecVtx.totalChiSquared(), TransSecVtx.degreesOfFreedom());
  cout << "Chi2 of SecVtx-Fit: " << chi.value() << endl;

  if ( type ) {

    fStuff->secondaryVertex = SecVtx;
    fillVertex(iEvent, iSetup, &TransSecVtx, type, ntracks);
  } 
  else { 

    kaonDeltaR(iEvent, &TransSecVtx, type, ntracks);
  }
}


// ----------------------------------------------------------------------
void Bs2MuMu::fillTrack(const edm::Event &iEvent, TAnaTrack *pTrack, reco::Track *it, int idx, int verb) {

  if ( verb == 1 ) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "==>fillTrack> Filling track #" << idx << ", event: " << fNevt << endl;
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
      
      cout << "--> Track #" << idx 
	   << " (pT = " << it->pt()
	   << ", eta = " << it->eta() << ", " << it->charge() << ")"
	   << " matched to Gen. Part. #" << index 
	   << " (" << tpr->pt()
	   << ", " << tpr->eta() << ")"
	   << ", PDG: " << type << ", NShared: " << assocChi2 << endl;
    }
    else if (verb == 1) {
      
      cout << "\t%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
    }
  } catch (Exception event) {

    if (verb == 1) {

      cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
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
	cout << " --> Track #" << idx 
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

  // TYPE = 0: B->mu mu K (all candidates)
  //      = 1: B->mu mu
  //      = 2: J/Psi->mu mu
  //      = 3: B->mu mu K

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>fillVertex> Filling vertex, event: " << fNevt << endl;
  
  using namespace edm;
  using namespace reco;

  std::vector<reco::TransientTrack> refTT = v->refittedTracks();
  std::vector<reco::TransientTrack> orgTT = v->originalTracks();

  std::vector<reco::Track> refitted;
  std::vector<reco::Track> original;

  for(vector<reco::TransientTrack>::const_iterator i = refTT.begin(); 
      i != refTT.end(); i++) {

    const Track & ftt = i->track();
    // cout << "mmmmmmm : " << ftt->innermomentum() << endl;
    refitted.push_back(ftt);
  }

  for(vector<reco::TransientTrack>::const_iterator i = orgTT.begin(); 
      i != orgTT.end(); i++) {

    const Track & ftt = i->track();
    original.push_back(ftt);
  }
  
  TLorentzVector m0, m1, kp;
  TLorentzVector bs;

  int i1(0), i2(0);

  fEff->Fill(30.1);

  if ( refitted.size() == ntracks ) {

    fEff->Fill(31.1);

    if ( ntracks == 2 ) {
      
      cout << "--- Vertex with number of refitted tracks: " << refitted.size() << " ---" << endl;
      
      m0.SetXYZM(refitted[0].px(),
		 refitted[0].py(),
		 refitted[0].pz(),
		 0.1056583);
    
      m1.SetXYZM(refitted[1].px(),
		 refitted[1].py(),
		 refitted[1].pz(),
		 0.1056583);
  
      bs = m0 + m1;
   
      double mass(0.);
      mass = bs.M();
      if ( mass ) {
	cout << "fillVertex> Inv. Mass = " << mass << endl;
	fM100->Fill(mass);
      }
      else {
	cout << "fillVertex> no mass" << endl;
      }

      i1 = fEvent->getSigTrack(0)->fIndex;
      i2 = fEvent->getSigTrack(1)->fIndex;
    }
      
    if ( ntracks == 3 ) {
      
      cout << "--- Vertex with number of refitted tracks: " << refitted.size() << " ---" << endl;
      
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
      

      double mass(0.);
      mass = bs.M();
      if ( mass ) {
	cout << "fillVertex> Inv. Mass = " << mass << endl;
	fM200->Fill(mass);
      }
      else {
	cout << "fillVertex> no mass" << endl;
      }

      i1 = fEvent->getSigTrack(0)->fIndex;
      i2 = fEvent->getSigTrack(2)->fIndex;  // This should be the K+
      
    }
  } else if ( original.size() == ntracks ) {

    cout << "--- Not enough refitted tracks for this Vertex: " << refitted.size()  << " ---" << endl;
    cout << "--- Use original tracks instead ---" << endl;
    
    fEff->Fill(32.1);

    if ( ntracks == 2 ) {
      
      cout << "--- Vertex with number of original tracks: " << original.size() << " ---" << endl;
      
      m0.SetXYZM(original[0].px(),
		 original[0].py(),
		 original[0].pz(),
		 0.1056583);
      
      m1.SetXYZM(original[1].px(),
		 original[1].py(),
		 original[1].pz(),
		 0.1056583);
      
      bs = m0 + m1;

      double mass(0.);
      mass = bs.M();
      if ( mass ) {
	cout << "fillVertex> Inv. Mass = " << mass << endl;
	fM100->Fill(mass);
      }
      else {
	cout << "fillVertex> no mass" << endl;
      }
      
      i1 = fEvent->getSigTrack(0)->fIndex;
      i2 = fEvent->getSigTrack(1)->fIndex;
    }
      
    if ( ntracks == 3 ) {
      
      cout << "--- Vertex with number of original tracks: " << original.size() << " ---" << endl;
      
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

      double mass(0.);
      mass = bs.M();
      if ( mass ) {
	cout << "fillVertex> Inv. Mass = " << mass << endl;
	fM200->Fill(mass);
      }
      else {
	cout << "fillVertex> no mass" << endl;
      }
      
      i1 = fEvent->getSigTrack(0)->fIndex;
      i2 = fEvent->getSigTrack(2)->fIndex;  // This should be the K+
      
    }
  } else {
      fEff->Fill(33.1);
      cout << "!!! Not enough tracks for this Vertex: " << refitted.size() 
	   << " refitted and " << original.size() << "original tracks !!!" << endl;
  }
  
  ChiSquared chi(v->totalChiSquared(), v->degreesOfFreedom());

  // -- Adding vertex to ntuple ...
  TAnaVertex *pVtx = new TAnaVertex();

    //    pVtx->setInfo(chi2.value(), chi2.degreesOfFreedom(), chi2.probability(), 0, 2); ??? change ndof to double ???
  pVtx->setInfo(chi.value(), int(chi.degreesOfFreedom()), chi.probability(), 1, type);
  pVtx->fPoint.SetXYZ(v->position().x(), v->position().y(), v->position().z());

  pVtx->addTrack(i1);
  pVtx->addTrack(i2);

  VertexDistanceXY axy;
  double dXY      = axy.distance(fStuff->primaryVertex2, *v).value();
  double dXYE     = axy.distance(fStuff->primaryVertex2, *v).error();
  double compXY   = axy.compatibility(fStuff->primaryVertex2, *v);

  VertexDistance3D a3d;
  double d3d      = a3d.distance(fStuff->primaryVertex2, *v).value();
  double d3dE     = a3d.distance(fStuff->primaryVertex2, *v).error();
  double comp3d   = a3d.compatibility(fStuff->primaryVertex2, *v);

  pVtx->fDxy  = dXY; 
  pVtx->fDxyE = dXYE; 
  pVtx->fCxy  = compXY; 

  pVtx->fD3d  = d3d; 
  pVtx->fD3dE = d3dE; 
  pVtx->fC3d  = comp3d; 


  // -- Adding Candidate to ntuple...
  TAnaCand  *pCand = fEvent->addCand();

  pCand->fPlab =  bs.Vect();
  pCand->fMass = bs.M();

  pCand->fSig1 = i1;  
  pCand->fSig2 = i2;  
  pCand->fType = type;
  
  pCand->fVtx  = *pVtx;

  cout << "   KVF Vertex/cand of type " << type << " points to signal tracks: " << i1 << " and " << i2 << endl;
  cout << "      mass: " << pCand->fMass << endl;
  cout << "      PVF xy: " << dXY << " +/- " << dXYE << "   comp: " << compXY << endl;
  cout << "      PVF 3d: " << d3d << " +/- " << d3dE << "   comp: " << comp3d << endl;


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

//**       cout << " Qual. = " << qual
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

//   cout << "      sim position: pos  = " 
//        << " " << pVtx->fSimPoint.X()
//        << " " << pVtx->fSimPoint.Y()
//        << " " << pVtx->fSimPoint.Z()
//        << endl;//
//

}
// ----------------------------------------------------------------------
int Bs2MuMu::jpsiCandidate(const edm::Event &iEvent) {  

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>jpsiCandidate> Looking for two muons in J/Psi-mass-window, event: " << fNevt << endl;

  TLorentzVector m0, m1, jpsi; 
  int i0(-1), i1(-1);
  double mJpsi(-1.); //double dmJpsi(9999.);
  int done(0);
  
  for (unsigned int i = 0; i < fStuff->BmmRecTracks.size(); ++i) {
      
    m0.SetXYZM((*fStuff->BmmRecTracks[i]).px(),
	       (*fStuff->BmmRecTracks[i]).py(),
	       (*fStuff->BmmRecTracks[i]).pz(),
	       0.1056583);

    for (unsigned int j = 0; j < fStuff->BmmRecTracks.size(); ++j) {
      
      if ( i == j ) { continue; }
      
      m1.SetXYZM((*fStuff->BmmRecTracks[j]).px(),
		 (*fStuff->BmmRecTracks[j]).py(),
		 (*fStuff->BmmRecTracks[j]).pz(),
		 0.1056583);
    
      jpsi = m0 + m1;
      mJpsi = jpsi.M();

      if ( mJpsi > 2.9 && mJpsi < 3.2 ) { 

	if ( fStuff->BmmRecTracks[i]->charge() != fStuff->BmmRecTracks[j]->charge() ) {

	  //	if ( abs(5.29 - mJpsi) < dmJpsi ) {	  
	  //	  dmJpsi = abs(5.29 - mJpsi);

	  i0 = i;
	  i1 = j;	
	  
	  cout << endl << "==>jpsiCandidate> : Found J/Psi candidates with muons \t #"
	       <<  fStuff->BmmRecTracksIndex[i0] << " \t #" <<  fStuff->BmmRecTracksIndex[i1]
	       << "\t --->  Inv. Mass: "  << mJpsi
	       << "( q1 = " << fStuff->BmmRecTracks[i0]->charge()
	       << ", q2 = " << fStuff->BmmRecTracks[i1]->charge()
	       << ")" << endl;

	  done = 1;
	  
	}
      }

      if ( done ) { break; }
    }

    if ( done ) { break; }
  }

  if ((i0 > -1) && (i1 > -1)) {
    
    // -- RecTrack vector
    const reco::Track* tM0 = fStuff->BmmRecTracks[i0];
    const reco::Track* tM1 = fStuff->BmmRecTracks[i1];
    
    fStuff->BmmRecTracks.clear();
    
    fStuff->BmmRecTracks.push_back(tM0);
    fStuff->BmmRecTracks.push_back(tM1);
    
    // -- Index into bmm-track vector 
    i0 = fStuff->BmmRecTracksIndex[i0];
    i1 = fStuff->BmmRecTracksIndex[i1];
    
    fStuff->BmmRecTracksIndex.clear();
    
    fStuff->BmmRecTracksIndex.push_back(i0);
    fStuff->BmmRecTracksIndex.push_back(i1);

    return 1;
  
  } else {
    return 0;
  }
}


// ----------------------------------------------------------------------
int Bs2MuMu::kaonCandidate(const edm::Event &iEvent, const edm::EventSetup& iSetup) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>kaonCandidate> Checking all tracks as potential kaon candidates, event: " << fNevt << endl;

  const reco::Track* tt = 0;

  int index(0);
  for (TrackCollection::const_iterator it = (*fStuff->theTkCollection).begin(); 
       it != (*fStuff->theTkCollection).end(); 
       ++it){
    
    if (index == fStuff->BmmRecTracksIndex[0]) { continue; }
    if (index == fStuff->BmmRecTracksIndex[1]) { continue; }

    tt = &(*it);

    fStuff->BmmRecTracks.push_back(tt);
    fStuff->BmmRecTracksIndex.push_back(index);

    bmmVertex(iEvent, iSetup, 0, 3);  //  Vertex TYPE = 0, B+ -> mu mu K with kaonDeltaR() 
                                      //   = list of candidates and filling in all vertices (with TYPE = 0)
    fStuff->BmmRecTracks.pop_back();
    fStuff->BmmRecTracksIndex.pop_back();
    index++;
  }

	
  double chi_min(99999.);
  int id(-1);
  for (unsigned int i = 0; i < fCand->SecVtxChi2.size(); ++i) {
    
    if ( fCand->SecVtxChi2[i] < chi_min ) {
      
      id = i;
      chi_min = fCand->SecVtxChi2[i];
    }
   }
  
  if ( fCand->BmmRecTracks.size() > 0 ) {
    cout << "======================================================================" << endl;
    cout << "==>kaonCandidate> Selected kaon candidate : Candidate #" << id 
 	 << " which corresponds to track #" << fCand->BmmRecTracksIndex[id] << endl;
    cout << " Inv. Mass m = " << fCand->InvMass[id] << " and Chi2 = " << fCand->SecVtxChi2[id] << endl;
    cout << "======================================================================" << endl;
    
    fK200->Fill( fCand->SecVtxChi2[id], fCand->InvMass[id] );
    
    return id;
  
  } else {
    return -1;
  }	
}


// ----------------------------------------------------------------------
void Bs2MuMu::kaonDeltaR(const edm::Event &iEvent, TransientVertex *v, int type, unsigned int ntracks) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>kaonDeltaR> Calculating deltaR of kaon candidate, event: " << fNevt << endl;

  std::vector<reco::TransientTrack> refTT = v->refittedTracks();
  std::vector<reco::Track> refitted;

  for(vector<reco::TransientTrack>::const_iterator i = refTT.begin(); 
      i != refTT.end(); i++) {

    const Track & ftt = i->track();
    refitted.push_back(ftt);
  }
  
  TLorentzVector m0, m1, kp;
  TLorentzVector jpsi, bs;

  fEff->Fill(30.1);

  if ( refitted.size() == ntracks ) {

    fEff->Fill(31.1);

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

    jpsi = m0 + m1;
    bs   = m0 + m1 + kp;

    TVector3 jpsi_plab = jpsi.Vect();
    TVector3 kp_plab   = kp.Vect();
    
    ChiSquared chi(v->totalChiSquared(), v->degreesOfFreedom());  
    
    int i1 = fStuff->BmmRecTracksIndex[0];
    int i2 = fStuff->BmmRecTracksIndex[2];

    // -- Adding vertex to ntuple ...
    TAnaVertex *pVtx = new TAnaVertex();
    
    //    pVtx->setInfo(chi2.value(), chi2.degreesOfFreedom(), chi2.probability(), 0, 2); ??? change ndof to double ???
    pVtx->setInfo(chi.value(), int(chi.degreesOfFreedom()), chi.probability(), 1, type);
    pVtx->fPoint.SetXYZ(v->position().x(), v->position().y(), v->position().z());
    
    pVtx->addTrack(i1);
    pVtx->addTrack(i2);
    
    VertexDistanceXY axy;
    double dXY      = axy.distance(fStuff->primaryVertex2, *v).value();
    double dXYE     = axy.distance(fStuff->primaryVertex2, *v).error();
    double compXY   = axy.compatibility(fStuff->primaryVertex2, *v);
    
    VertexDistance3D a3d;
    double d3d      = a3d.distance(fStuff->primaryVertex2, *v).value();
    double d3dE     = a3d.distance(fStuff->primaryVertex2, *v).error();
    double comp3d   = a3d.compatibility(fStuff->primaryVertex2, *v);
    
    pVtx->fDxy  = dXY; 
    pVtx->fDxyE = dXYE; 
    pVtx->fCxy  = compXY; 
    
    pVtx->fD3d  = d3d; 
    pVtx->fD3dE = d3dE; 
    pVtx->fC3d  = comp3d; 
  

    // -- Adding Candidate to ntuple...
    TAnaCand  *pCand = fEvent->addCand();
    
    pCand->fPlab =  bs.Vect();
    pCand->fMass = bs.M();
    
    pCand->fSig1 = i1;  
    pCand->fSig2 = i2;  
    pCand->fType = type;
    
    pCand->fVtx  = *pVtx;
    
    
    double dphi = jpsi_plab.DeltaPhi(kp_plab);
    double deta = kp_plab.Eta() - jpsi_plab.Eta();
    double dr   = TMath::Sqrt(dphi*dphi + deta*deta);

    if ( dr < 1.5 ) {
      cout << "==>kaonDeltaR> Adding kaon candidate #" << fCand->SecVtxChi2.size()+1 
	   <<", corresponding to track #" << fStuff->BmmRecTracksIndex[2] 
	   << " with dr = " << dr << ", Bs mass " << bs.M() << " and chi2 " << chi.value() << endl;
      
      fCand->InvMass.push_back(bs.M());
      fCand->SecVtxChi2.push_back(chi.value());

      fCand->BmmRecTracks.push_back(fStuff->BmmRecTracks[2]);
      fCand->BmmRecTracksIndex.push_back(fStuff->BmmRecTracksIndex[2]);

      fK100->Fill(chi.value(), bs.M());

    }
    else {
      cout << "==>kaonDeltaR> !!! Kaon candidate rejected, ouside cone:  dr = " << dr
	   << ", Bs mass " << bs.M() << " and chi2 " << chi.value() << " !!!" <<endl;
    }
    
  } else {
      fEff->Fill(33.1);
      cout << "==>kaonDeltaR> !!! Not enough tracks for this Vertex: " << refitted.size() 
	   << " refitted tracks !!!" << endl;
  }  
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

  // cout << mdpt << " " << mdphi << " " << mdeta << " " << found << endl;

  if ((mdpt < ept)
      && (mdphi < ephi)
      && (mdeta < eeta)
      ) {
    return found;
  } else {
    return -1;
  }

}


// ==============================================================================
//                           PRINT OUT
// ==============================================================================

void Bs2MuMu::printGenTracks(const edm::Event &iEvent) {
  
  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>printGenTracks> Starting to print generator block, event: " << fNevt << endl;


  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel(fSourceLabel.c_str(), evt);  
  const HepMC::GenEvent *genEvent = evt->GetEvent();

  int gcnt = genEvent->particles_size();
  cout << "Counted " << gcnt  << " genTracks in generator block" << endl;
  gcnt = 0;
  
  for (HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin();
       p != genEvent->particles_end();
       ++p) {
      
    cout << "GenTrack " <<  gcnt
	 << "   Particle ID: " << (*p)->pdg_id()
	 << ", Status "        << (*p)->StatusCode()
	 << ", pT "            << (*p)->Momentum().perp()
	 << ", eta "           << (*p)->Momentum().eta() << endl;
    
    gcnt++;
  }
}

// ----------------------------------------------------------------------
void Bs2MuMu::printSimTracks(const edm::Event &iEvent) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>printSimTracks> Starting to print simulated tracks, event: " << fNevt << endl;


  int scnt = (*fStuff->theTPCollection).size();
  cout << "Counted " << scnt  << " simulated tracks." << endl;
  scnt = 0;

  const HepMC::GenParticle* genPar = 0; 
  int gen_pdg_id(-99999), gen_id(-99999); 
  int gen_cnt(0);

  for (TrackingParticleCollection::const_iterator tpr = (*fStuff->theTPCollection).begin(); 
       tpr != (*fStuff->theTPCollection).end(); 
       ++tpr){

    cout << "SimTrack #" <<  scnt
	 << ": pT "     <<  (*tpr).pt()
	 << ", eta "    <<  (*tpr).eta()
	 << ", charge " <<  (*tpr).charge()
	 << ", PDG ID " <<  (*tpr).pdgId() << endl;
    

    // -- Using fEvent->getGenIndex to find matching track in generator block
    int index = fEvent->getGenIndex((*tpr).momentum().x(), (*tpr).momentum().y(), (*tpr).momentum().z(), 
				    (*tpr).pdgId());
    if ( index < 0 ) {
      
      cout << "%%> Did not find a match in generator block (fEvent->getGenIndex)" << endl;
    }
    else {
      cout << "==> matches to generator block particle #" << index << " (fEvent->getGenIndex)" <<endl;
    }


    // -- Using CMSSW's genp_iterator to find matching track in generator block
    gen_cnt = 0; gen_pdg_id = -99999; gen_id = -99999;
    for ( TrackingParticle::genp_iterator pit=tpr->genParticle_begin(); pit!=tpr->genParticle_end(); pit++ ){

      genPar=pit->get();

      gen_pdg_id = (*genPar).ParticleID();
      gen_id     = (*genPar).barcode()-1;

      if ( gen_id == index ) {

	cout << "==> matches to generator block particle #" << gen_id
	     << " with pT = " << (*genPar).Momentum().perp() 
	     << " and (gen) PDG ID: " << gen_pdg_id << endl;
      }
      else {
	cout << "==> matches to generator block particle #" << gen_id << " (disagreement!)"
	     << " with pT = " << (*genPar).Momentum().perp() 
	     << " and (gen) PDG ID: " << gen_pdg_id << endl;
      }

      gen_cnt++;

    }
    

    if ( gen_cnt == 0 ) {
	
      cout << "%%> no match in gen. block for sim. particle #" << scnt << endl;
    }
    else if ( gen_cnt > 1 ) {

      cout << " ==> !!! More than one gen. particle for sim. particle #" << scnt << " !!!" << endl;
    }

    cout << endl;
    scnt++;
  }
}


// ----------------------------------------------------------------------
void Bs2MuMu::printRecTracks(const edm::Event &iEvent) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>printRecTracks> Starting to print reconstructed tracks, event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );

 // -- track association
  reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection); 

  int index(-99999), type(-99999);

  int rcnt = (*fStuff->theTkCollection).size();
  cout << "Counted " << rcnt  << " reconstructed tracks." << endl;
  rcnt = 0;  

  for (TrackCollection::const_iterator itTrack = (*fStuff->theTkCollection).begin(); 
       itTrack != (*fStuff->theTkCollection).end(); 
       ++itTrack){

    // -- Track parameters
    cout << "RecTrack #"  <<  rcnt
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

	  cout << "\t ----> Found gen. MC Track #" << index << " with PDG ID " << type 
	       << " (NShared: " << assocChi2 << ")" << endl;
	}
	else {

	  cout << "\t%%> no match in gen. block for sim. particle #" << tpr.index() << endl;
	}
      }   

    } catch (Exception event) {

      cout << "\t%%> no MC match found for rec. Track #" << rcnt << endl;
    }

    rcnt++;
  }    
}

// ----------------------------------------------------------------------
void Bs2MuMu::printMuonTracks(const edm::Event &iEvent) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>printMuonTracks> Starting to print tracks reconstructed as muons, event: " << fNevt << endl;


  int mcnt = (*fStuff->theMuonCollection).size();
  cout << "Counted " << mcnt  << " reconstructed tracks." << endl;
  mcnt = 0;  

  for (MuonCollection::const_iterator glbMuon = (*fStuff->theMuonCollection).begin(); 
       glbMuon != (*fStuff->theMuonCollection).end(); 
       ++glbMuon){

    TrackRef glbTrack = glbMuon->combinedMuon();

    cout << "MuonTrack #" <<  mcnt
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
 

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==>printReco2Sim> Associating recontructed tracks to simulated tracks, event: " << fNevt << endl;

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

    cout << "-- Associator by hits --" << endl;  
    reco::RecoToSimCollection recSimColl = 
      associatorByHits->associateRecoToSim (tracks, TPCollectionH, &iEvent);
    
    double ptRes(1000.);
    for(TrackCollection::size_type i=0; i<recTC.size(); ++i) {
      TrackRef track(tracks, i);
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
	cout << "Reco Track " << setw(2) << track.index() << " pT: "  << setw(6) << track->pt() 
	     <<  " matched to " << tp.size() << " MC Tracks" << std::endl;

	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  double assocChi2 = it->second;
	  cout << "\t\tMCTrack " << setw(2) << tpr.index() << " pdgId: " << setw(2) << tpr->pdgId() 
	       << " pT: " << setw(6) << tpr->pt() 
	       << " NShared: " << assocChi2 << endl;
	  
 	  ptRes=track->pt()-sqrt(tpr->momentum().Perp2());
 	  cout << "  res = " << ptRes << endl;
 	  fPT310->Fill(ptRes);
	}
      } catch (Exception event) {

	cout << "%%> Rec. Track #" << setw(2) << track.index() << " pT: " 
	     << setprecision(2) << setw(6) << track->pt() 
	     <<  " matched to 0 MC Tracks" << endl;
      }
    }
  }
  
  // -- perform association by chi2
  else if ( !strcmp(option,"chi2") ) {

    cout << "-- Associator by chi2 --" << endl;  
    reco::RecoToSimCollection recSimColl = 
      associatorByChi2->associateRecoToSim (tracks, TPCollectionH, &iEvent );

    double ptRes(1000.);
    for(TrackCollection::size_type i=0; i<recTC.size(); ++i) {
      TrackRef track(tracks, i);
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
	cout << "Reco Track " << setw(2) << track.index() << " pT: "  << setw(6) << track->pt() 
	     <<  " matched to " << tp.size() << " MC Tracks" << std::endl;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  double assocChi2 = it->second;
	  cout << "\t\tMCTrack " << setw(2) << tpr.index() << " pdgId: " << setw(2) << tpr->pdgId() 
	       << " pT: " << setw(6) << tpr->pt() 
	       << " NShared: " << assocChi2 << endl;
	  
	  ptRes=track->pt()-sqrt(tpr->momentum().Perp2());
	  cout << "  res = " << ptRes << endl;
	  fPT320->Fill(ptRes);
	}
      } catch (Exception event) {

	cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
	     << setprecision(2) << setw(6) << track->pt() 
	     <<  " matched to 0 MC Tracks" << endl;
      }
    }
  }
  else {

    cout << "-- " << fAssocLabel.c_str() << " --" << endl;  

    reco::RecoToSimCollection recSimColl = (*fStuff->recSimCollection);

    double ptRes(1000.);

    for(TrackCollection::size_type i=0; i<(*fStuff->theTkCollection).size(); ++i) {

      TrackRef track(tracks, i);
      
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
	cout << "Reco Track " << setw(2) << track.index() << " pT: "  << setw(6) << track->pt() 
	     <<  " matched to " << tp.size() << " MC Tracks" << std::endl;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  double assocChi2 = it->second;
	  cout << "\t\tMCTrack " << setw(2) << tpr.index() << " pdgId: " << setw(2) << tpr->pdgId() 
	       << " pT: " << setw(6) << tpr->pt() 
	       << " NShared: " << assocChi2 << endl;
	  
	  ptRes=track->pt()-sqrt(tpr->momentum().Perp2());
	  cout << "  res = " << ptRes << endl;
	  fPT300->Fill(ptRes);
	}
      } catch (Exception event) {

	cout << "%%>   Rec. Track #" << setw(2) << track.index() << " pT: " 
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
