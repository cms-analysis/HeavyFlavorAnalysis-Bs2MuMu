#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpCandidates.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" 
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "SimTracker/Records/interface/VertexAssociatorRecord.h"
#include "SimTracker/VertexAssociation/interface/VertexAssociatorBase.h" 
#include "SimTracker/VertexAssociation/interface/VertexAssociatorByTracks.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"


// -- Yikes!
extern TAna00Event *gBSEvent;
extern TFile       *gBSFile;


using namespace std;
using namespace reco;
using namespace edm;



// ----------------------------------------------------------------------
// ======================================================================
BSDumpCandidates::BSDumpCandidates(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fBmmSel(iConfig.getUntrackedParameter<int>("bmmsel", -1)),
  fChannel(iConfig.getUntrackedParameter<string>("channel","default")),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<string>("PrimaryVertexLabel"
							    , string("offlinePrimaryVerticesFromCTFTracks"))),
  fGenEventLabel(iConfig.getUntrackedParameter<string>("generatorEventLabel", string("source"))),
  fTrackingParticlesLabel(iConfig.getUntrackedParameter<string>("trackingParticlesLabel", string("trackingParticles"))),
  fTrackingVertexLabel(iConfig.getUntrackedParameter<string>("trackingVertexLabel", string("trackingParticles"))),
  fTracksLabel(iConfig.getUntrackedParameter<string>("tracksLabel", string("ctfWithMaterialTracks"))),
  fMuonsLabel1(iConfig.getUntrackedParameter<InputTag>("muonsLabel1")),
  fMuonsLabel2(iConfig.getUntrackedParameter<InputTag>("muonsLabel2")),
  fVtxAssociatorLabel(iConfig.getUntrackedParameter<string>("vertexAssociatorLabel", string("VertexAssociatorByTracks"))), 
  fAssociatorLabel(iConfig.getUntrackedParameter<string>("associatorLabel", string("TrackAssociatorByChi2"))), 
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1))  {

  using namespace std;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- BSDumpCandidates constructor" << endl;
  cout << "--- Verbose                : " << fVerbose << endl;
  cout << "--- BmmSel                 : " << fBmmSel << endl;
  cout << "--- Channel                : " << fChannel.c_str() << endl;
  cout << "--- prim. Vertex Label     : " << fPrimaryVertexLabel.c_str() << endl;
  cout << "--- generatorEventLabel    : " << fGenEventLabel.c_str() << endl;
  cout << "--- trackingParticlesLabel : " << fTrackingParticlesLabel.c_str() << endl;
  cout << "--- tracksLabel            : " << fTracksLabel.c_str() << endl;
  cout << "--- muonsLabel1            : " << fMuonsLabel1 << endl;
  cout << "--- muonsLabel2            : " << fMuonsLabel2 << endl;
  cout << "--- associatorLabel        : " << fAssociatorLabel.c_str() << endl;
  cout << "--- doTruthMatching        : " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS
  cout << "----------------------------------------------------------------------" << endl;

  cout << "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl;
  cout << "===>> BSDumpCandidates >>> ctor, decay mode" << endl;


  // -- Counters
  fNevt = 0; 
  fNgen = 0;
  fNrec = 0;


  if ( strcmp( (fAssociatorLabel.c_str()), "TrackAssociatorByHits") && 
       strcmp( (fAssociatorLabel.c_str()), "TrackAssociatorByChi2") )  {
  
    cout << " Please set your track associator option to either \"TrackAssociatorByHits\"                                                              or \"TrackAssociatorByChi2\" .... Abort!" << endl;
    return;
  }

  // -- Decay mode
  decayChannel(fChannel.c_str());

  cout << "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl << endl << endl;

}

// ======================================================================
// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// ======================================================================

BSDumpCandidates::~BSDumpCandidates() {

  cout << endl << "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl;
  cout << "===>> BSDumpCandidates >>> dtor, bye ..." << endl;
  cout << "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl << endl << endl;

}

// ======================================================================
// ----------------------------------------------------------------------



// ------------ method called to for each event  ------------
void BSDumpCandidates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
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
  cout << "===>> BSDumpCandidates >>> Start with event: " << fNevt << endl;

  using reco::TrackCollection;
  using reco::MuonCollection;

  // === Initialize event record ===

  theGlobalMuonCollection = 0;
  theTrackerMuonCollection = 0;
  theTkCollection = 0;
  theGenCollection.clear();
  //  recSimCollection = 0;
  //  recSimCollectionVertex = 0;

  clearTracks();
  clearCandidateTracks();

  RecTracks.clear();
  RecTracksIndex.clear();

  // primaryVertex2 = 0;

  // === Fill event record ===

  // -- Get generator block directly
  Handle<HepMCProduct> evt;
  iEvent.getByLabel(fGenEventLabel.c_str(), evt);
  const HepMC::GenEvent *genEvent = evt->GetEvent();

  for (HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin();
       p != genEvent->particles_end();
       ++p) { 

    theGenCollection.push_back((*p));
  }

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  theTkCollection  = tracks.product();
  // theTkCollection  = new reco::TrackCollection(*(tracks.product()));
   
  // -- get the collection of globalMuonTracks 
  edm::Handle<reco::MuonCollection> MuCollection;
  iEvent.getByLabel(fMuonsLabel1, MuCollection);
  theGlobalMuonCollection   = MuCollection.product();
  // theGlobalMuonCollection   = new reco::MuonCollection(*(MuCollection.product()));  
   
  // -- get the collection of trackerMuonTracks 
  edm::Handle<reco::MuonCollection> tkMuCollection;
  iEvent.getByLabel(fMuonsLabel2, tkMuCollection);
  theTrackerMuonCollection   = tkMuCollection.product();
  // theTrackerMuonCollection   = new reco::MuonCollection(*(tkMuCollection.product()));  

  // -- get the collection of sim. Tracks
//   Handle<SimTrackContainer> simTrackCollection;
//   iEvent.getByLabel("g4SimHits", simTrackCollection);
//   const SimTrackContainer simTC = *(simTrackCollection.product());

  // -- get the collection of sim. Vertices
//   Handle<SimVertexContainer> simVertexCollection;
//   iEvent.getByLabel("g4SimHits", simVertexCollection);
//   const SimVertexContainer simVC = *(simVertexCollection.product());

  // -- get the collection of TrackingParticles 
  edm::Handle<TrackingParticleCollection> trackingParticles;
  iEvent.getByLabel(fTrackingParticlesLabel.c_str(), trackingParticles);

  recSimCollection =  new 
      reco::RecoToSimCollection(fAssociator->associateRecoToSim(tracks, trackingParticles, &iEvent)); 

  // -- get the collection of TrackingVertices
  edm::Handle<TrackingVertexCollection>  TVCollectionH ;
  iEvent.getByLabel(fTrackingVertexLabel.c_str(), TVCollectionH);


  // -- get the collection of primary Vertices
//  edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
//  iEvent.getByLabel(fPrimaryVertexLabel.c_str(), recoPrimaryVertexCollection);

//  recSimCollectionVertex =  new 
//    reco::VertexRecoToSimCollection(fVtxAssociator->associateRecoToSim(recoPrimaryVertexCollection, 
// 									  TVCollectionH, iEvent,
//  									  (*recSimCollection)));
  // === Start analysis ===

  // -- Primary Vertex
  int npv = primaryVertex(iEvent);
  if (!npv) {
    cout << "*********** No primary vertex found !!!" << endl;  
  }

  // -- (Rec.) Signal Tracks
  if (fBmmSel == 1) {

    bmmTracks1(iEvent);
    truthCandTracks(iEvent, iSetup);
    secondaryVertex(iEvent, iSetup);

  } else if (fBmmSel == 2) {

    bmmTracks2(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

  } else if (fBmmSel == 3) {

    bmmTracks3(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

  } else if (fBmmSel == 4) {

    bmmTracks4(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

  } else {

    fBmmSel = 1;
    bmmTracks1(iEvent);
    truthCandTracks(iEvent, iSetup);
    secondaryVertex(iEvent, iSetup);

    fBmmSel = 2;
    bmmTracks2(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

    fBmmSel = 3;
    bmmTracks3(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

    fBmmSel = 4;
    bmmTracks4(iEvent);
    muonCandTracks(iEvent, iSetup, fMass, fMass2);
    kaonCandTracks(iEvent, iSetup, 1.5);
    secondaryVertex(iEvent, iSetup);

    fBmmSel = -1;
  }

  // -- Dump tree
//   fTree->Fill();

  delete recSimCollection;
  // delete recSimCollectionVertex;

  cout << endl << "===>> BSDumpCandidates >>> Done with event: " << fNevt << endl;
  cout << "*********************************************************************" << endl;
    
}


// ----------------------------------------------------------------------
// ======================================================================



// ----------------------------------------------------------------------
void BSDumpCandidates::bmmTracks1(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>bmmTracks1> Matching particle, truth-matched to be particle PDG #" 
       << fTruthMC_I << " and a decay products of particle PDG #" << fTruthMC_mom << ", event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product()); 
 
 // -- track association
  reco::RecoToSimCollection recSimColl = (*recSimCollection); 


  // -- Clear tracks
  clearTracks();
  clearCandidateTracks();

  const reco::Track* tt = 0;
  const HepMC::GenParticle* genGmo = 0;
  const HepMC::GenParticle* genMom = 0;
  const HepMC::GenParticle* genPar = 0; 

  int mcand1(0), mcand2(0), kcand(0);
  int gmo_pdg_id(-99999), gmo_id(-99999);
  int mom_pdg_id(-99999), mom_id(-99999);
  int gen_pdg_id(-99999), gen_id(-99999);
  int gen_cnt(0);

  int motherBarcode(-99999), grandmotherBarcode(-99999);

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

	gen_pdg_id = genPar->pdg_id();
	gen_id     = genPar->barcode()-1;
	
	if (fVerbose) cout << ".. Gen. Particle #" << gen_id << " pT = " << genPar->momentum().perp() << " ---> PDG ID: " << gen_pdg_id << endl;
	


	// -- Two body decays
	if ( (abs(gen_pdg_id) == fTruthMC_I) ||  
	     (abs(gen_pdg_id) == fTruthMC_II) ) { 
	  
	  motherBarcode = -99999;

	  motherBarcode = genPar->production_vertex() && 
	    genPar->production_vertex()->particles_in_const_begin() !=
	    genPar->production_vertex()->particles_in_const_end() ?
	    (*(genPar->production_vertex()->particles_in_const_begin()))->barcode()-1 : 0;
	  
	  if ( motherBarcode < 0 ) {

	    if ( fVerbose ) {
	      cout << " --> No good mother barcode " << motherBarcode << " <--- " << endl;
	      continue;
	    }
	  }

	  if ( motherBarcode > theGenCollection.size() ) {
	    
	    if ( fVerbose ) {
	      cout << " --> Mother barcode " << motherBarcode << " outside GenCollection <--- " << endl;
	      continue;
	    }
	  }

    	  genMom = theGenCollection[motherBarcode];
	  
	  mom_pdg_id = genMom->pdg_id();
	  mom_id     = genMom->barcode()-1;
	  
	  
	  if (fVerbose) cout << endl << ".. and Mother Particle #" << mom_id << " pT = " << genMom->momentum().perp() 
			     << " ---> PDG ID: " << mom_pdg_id;
	  
	  // -- Cand 1 (Bs)
	  if ( abs(mom_pdg_id) == fTruthMC_mom ) { 
	    
	    mcand1++;
	    
	    tt = &(*track);
	    
	    BmmRecTracks.push_back(tt);
	    BmmRecTracksIndex.push_back(track.index());
	    BmmRecTracksB.push_back(mom_id);

	    MuonRecTracks.push_back(tt);
	    MuonRecTracksIndex.push_back(track.index());
	    
	    if (fVerbose) cout << "    *** " << fPrintChannel << " (" << mcand1 << ") *** ";
	  }	  
	  
	  // -- Cand 2 (from J/Psi or from B directly)
	  if ( abs(mom_pdg_id) == fTruthMC2_mom ) { 

	    if ( fTruthMC2_gmo > 0 ) {
	      
	      grandmotherBarcode = -99999;
	  
	      grandmotherBarcode = genMom->production_vertex() && 
		genMom->production_vertex()->particles_in_const_begin() !=
		genMom->production_vertex()->particles_in_const_end() ?
		(*(genMom->production_vertex()->particles_in_const_begin()))->barcode()-1 : 0;
	      	  
	      
	      if ( grandmotherBarcode < 0 ) {
		
		if ( fVerbose ) {
		  cout << " --> No good grandmother barcode " << grandmotherBarcode << " <--- " << endl;
		  continue;
		}
	      }
	      
	      if ( grandmotherBarcode > theGenCollection.size() ) {
		
		if ( fVerbose ) {
		  cout << " --> Grandmother barcode " << grandmotherBarcode << " outside GenCollection <--- " << endl;
		  continue;
		}
	      }

	      genGmo = theGenCollection[grandmotherBarcode];
	      
	      gmo_pdg_id = genGmo->pdg_id();
	      gmo_id     = genGmo->barcode()-1;
	      
	      if (fVerbose) cout << endl << ".. and Grandmother Particle #" << gmo_id << " pT = " 
				 << genGmo->momentum().perp() << " ---> PDG ID: " << gmo_pdg_id;
	      
	      if ( abs(gmo_pdg_id) == fTruthMC2_gmo ) { 
		
		mcand2++;
		
		tt = &(*track);
		
		JpsiRecTracks.push_back(tt);
		JpsiRecTracksIndex.push_back(track.index());
		JpsiRecTracksB.push_back(gmo_id);
		
		MuonRecTracks.push_back(tt);
		MuonRecTracksIndex.push_back(track.index());
		
		if (fVerbose) cout << "    *** " << fPrintChannel2 << " (" << mcand2 << ". mu) *** ";
	      }
	    
	    } else { 
		
		mcand2++;
		
		tt = &(*track);
		
		JpsiRecTracks.push_back(tt);
		JpsiRecTracksIndex.push_back(track.index());
		JpsiRecTracksB.push_back(mom_id);
		
		MuonRecTracks.push_back(tt);
		MuonRecTracksIndex.push_back(track.index());
		
		if (fVerbose) cout << "    *** " << fPrintChannel2 << " (" << mcand2 << ". mu) *** ";


	    }
	    
	    if (fVerbose) cout << endl; 
	  }
	}
	
	// -- Cand 3 (only Kaon or 3rd track)
	if ( abs(gen_pdg_id) == fTruthMC2  ) { 
	  
	  motherBarcode = -99999;	  

	  motherBarcode = genPar->production_vertex() && 
	    genPar->production_vertex()->particles_in_const_begin() !=
	    genPar->production_vertex()->particles_in_const_end() ?
	    (*(genPar->production_vertex()->particles_in_const_begin()))->barcode()-1 : 0;
	  	  
	  if ( motherBarcode < 0 ) {

	    if ( fVerbose ) {
	      cout << " --> No good mother barcode <--- " << endl;
	      continue;
	    }
	  }
	  
	  if ( motherBarcode > theGenCollection.size() ) {
	    
	    if ( fVerbose ) {
	      cout << " --> Mother barcode outside GenCollection <--- " << endl;
	      continue;
	    }
	  }

    	  genMom = theGenCollection[motherBarcode];
	  mom_pdg_id = genMom->pdg_id();
	  mom_id     = genMom->barcode()-1;
	  
	  
	  if (fVerbose) cout << endl << ".. and Mother Particle #" << mom_id << " pT = " << genMom->momentum().perp() 
			     << " ---> PDG ID: " << mom_pdg_id;
	  
	  if ( fTruthMC2_gmo > 0 ) {

	    if ( abs(mom_pdg_id) == fTruthMC2_gmo ) { 
	      
	      kcand++;
	      
	      tt = &(*track);
	    
	      KaonRecTracks.push_back(tt);
	      KaonRecTracksIndex.push_back(track.index());
	      KaonRecTracksB.push_back(mom_id);
	      
	      if (fVerbose) cout << "    *** " << fPrintChannel2 << " (" << kcand << ". K or 3rd track) *** ";
	    }
	  } else {

	    if ( abs(mom_pdg_id) == fTruthMC2_mom ) { 
	      
	      kcand++;
	      
	      tt = &(*track);
	    
	      KaonRecTracks.push_back(tt);
	      KaonRecTracksIndex.push_back(track.index());
	      KaonRecTracksB.push_back(mom_id);
	      
	      if (fVerbose) cout << "    *** " << fPrintChannel2 << " (" << kcand << ". K or 3rd track) *** ";
	    }
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
void BSDumpCandidates::bmmTracks2(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>bmmTracks2> Matching particle, truth matched to be particle PDG #" 
       << fTruthMC_I << " in the generator block, event: " << fNevt << endl;

  // -- get the collection of RecoTracks 
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks );  
  const  reco::TrackCollection  recTC = *(tracks.product());

 // -- track association
  reco::RecoToSimCollection recSimColl = (*recSimCollection);

  // -- Clear tracks
  clearTracks();
  clearCandidateTracks();

  const reco::Track* tt = 0;
  const HepMC::GenParticle* genPar = 0; 

  int mcand(0);
  int gen_pdg_id(-99999), gen_id(-99999);
  int gen_cnt(0);

  TH1D *h1 = (TH1D*)gBSFile->Get("eff");

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
		
	gen_pdg_id = genPar->pdg_id();
	gen_id     = genPar->barcode()-1;

	if ( abs(gen_pdg_id) == 13 ) {
	  
	  mcand++;
	  tt = &(*track);
	    
	  MuonRecTracks.push_back(tt);
	  MuonRecTracksIndex.push_back(track.index());

	  if (fVerbose) {
	    
	    cout << ".. Rec. Track #" << track.index() << " pT = " << track->pt() << endl;
	    cout << ".. Sim. Particle #" << tpr.index() << " pT = " << tpr->pt() 
		 << " PDG ID: " << tpr->pdgId() << " (NShared: "  << assocChi2 << ")" << endl;
	    cout << ".. Gen. Particle #" << gen_id << " pT = " << genPar->momentum().perp() 
		 << " ---> PDG ID: " << gen_pdg_id << "    *** " << fPrintChannel << " (" << mcand << ") *** " << endl;
	  }
	}
	
	gen_cnt++;
      }

      // -- eff. histogram ---------------------------

      if ( MuonRecTracks.size() >= 0 ) {

	if ( MuonRecTracks.size() < 9 ) {
	  
	  h1->Fill(60.1 + MuonRecTracks.size());

	} else {
	  
	  h1->Fill(69.1);
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
void BSDumpCandidates::bmmTracks3(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>bmmTracks3> Get global muon tracks, event: " << fNevt << endl;

  // -- Clear tracks
  clearTracks();
  clearCandidateTracks();

  int globalMuon  = theGlobalMuonCollection->size();
  int trackerMuon      = theTrackerMuonCollection->size();

  if (fVerbose) cout << "==>bmmTracks3>  muons  found: " << globalMuon  << " global muons,  "
		     << trackerMuon << " tracker muons"  << endl;  

  const reco::Track* tt = 0;

  // -- Get global muons
  int mcnt(0), mcand(0);
  for (MuonCollection::const_iterator glbMuon = (*theGlobalMuonCollection).begin(); 
       glbMuon != (*theGlobalMuonCollection).end(); 
       ++glbMuon){

    TrackRef glbTrack = glbMuon->combinedMuon();
    tt = &(*glbTrack);

    int idrec = (glbMuon->track()).index();

    if ( idrec > -1 ) {
	    
      mcand++;
      MuonRecTracks.push_back(tt);
      MuonRecTracksIndex.push_back(idrec);

    } else {

      if (fVerbose) cout << "==>bmmTracks3> Coulnd't find rec. track for global muon #" << mcnt << endl;  
    }
    
    mcnt++;

  }


  // -- Find tracker muons, that were not found as global muons
  int tcnt(0), tcand(0);
  for (MuonCollection::const_iterator tkMuon = (*theTrackerMuonCollection).begin(); 
       tkMuon != (*theTrackerMuonCollection).end(); 
       ++tkMuon){

    TrackRef tkTrack = tkMuon->track();
    tt = &(*tkTrack);

    int idrec = (tkMuon->track()).index();

    int isGlobal(0);
    for ( unsigned int i=0; i<MuonRecTracksIndex.size(); i++ ) { 
      isGlobal = 0;
      if (idrec == MuonRecTracksIndex[i]) {
	isGlobal = 1;
	break;
      }
    }

    if ( !isGlobal && (idrec > -1) ) {
	  
      tcand++;  
      MuonRecTracks.push_back(tt);
      MuonRecTracksIndex.push_back(idrec);

    } else {

      if (fVerbose) {
	if ( idrec < 0 ) cout << "==>bmmTracks3> Coulnd't find rec. track for tracker muon #" << tcnt << endl;  
      }
    }
    
    tcnt++;
  }


  if (fVerbose) cout << "==>bmmTracks3> Taking " << mcand << "/" << mcnt << " global and " 
		     << tcand  << "/" << tcnt <<  " tracker muons " <<  endl;


  // -- eff. histogram ---------------------------

  TH1D *h1 = (TH1D*)gBSFile->Get("eff");

  if ( MuonRecTracks.size() >= 0 ) {
    if ( MuonRecTracks.size() < 9 ) {
      h1->Fill(70.1 + MuonRecTracks.size());
    } else {
      h1->Fill(79.1);
    }
  }

  if ( mcnt < 9 ) {
    h1->Fill(170.1 + mcnt);
  } else {
    h1->Fill(179.1);
  }

  if ( tcnt < 9 ) {
    h1->Fill(180.1 + tcnt);
  } else {
    h1->Fill(189.1);
  }
    
  if ( mcand < 9 ) {
    h1->Fill(270.1 + mcand);
  } else {
    h1->Fill(279.1);
  }

  if ( tcand < 9 ) {
    h1->Fill(280.1 + tcand);
  } else {
    h1->Fill(289.1);
  }

  // ----------------------------------------------
}


// ----------------------------------------------------------------------
void BSDumpCandidates::bmmTracks4(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>bmmTracks4> Matching tracks, rec PID as muons, event: " << fNevt << endl;

  // -- Clear tracks
  clearTracks();
  clearCandidateTracks();

  int globalMuon   = theGlobalMuonCollection->size();
  int trackerMuon = theTrackerMuonCollection->size();
  int nTrack       = theTkCollection->size();

  if (fVerbose) cout << "==>bmmTracks4>  tracks found: " << nTrack << endl
		     << "==>bmmTracks4>  muons  found: " << globalMuon  << " global muons,  "
		     << trackerMuon << " tracker muons"  << endl;  

 
  const reco::Track* glt = 0;

  // -- Get tracks matched to global muons
  int mcnt(0), mcand(0);
  for (MuonCollection::const_iterator glbMuon = (*theGlobalMuonCollection).begin(); 
       glbMuon != (*theGlobalMuonCollection).end(); 
       ++glbMuon){

    TrackRef glbTrack = glbMuon->combinedMuon();
    glt = &(*glbTrack);

    int idrec = (glbMuon->track()).index();

    if ( (idrec > -1) && (idrec < (*theTkCollection).size()) ) {
	    
      int index(0);
      for (TrackCollection::const_iterator it = (*theTkCollection).begin(); 
	   it != (*theTkCollection).end(); 
	   ++it){
	
	if (index == idrec ) {	
	  glt = &(*it);
	  break;
	}

	index++;
      }
	
      mcand++;
      MuonRecTracks.push_back(glt);
      MuonRecTracksIndex.push_back(idrec);

    } else {
      if (fVerbose) cout << "==>bmmTracks4> Coulnd't find rec. track for global muon #" << mcnt << endl;  
    }
    
    mcnt++;

  }


  const reco::Track* tkt = 0;

  // -- Get tracks matched to tracker muons, that were not found as global muons
  int tcnt(0), tcand(0);
  for (MuonCollection::const_iterator tkMuon = (*theTrackerMuonCollection).begin(); 
       tkMuon != (*theTrackerMuonCollection).end(); 
       ++tkMuon){

    TrackRef tkTrack = tkMuon->track();
    tkt = &(*tkTrack);

    int idrec = (tkMuon->track()).index();

    int isGlobal(0);
    for ( unsigned int i=0; i<MuonRecTracksIndex.size(); i++ ) { 
      isGlobal = 0;
      if (idrec == MuonRecTracksIndex[i]) {
	isGlobal = 1;
	break;
      }
    }

    if ( !isGlobal && (idrec > -1) && (idrec < (*theTkCollection).size()) ) {
 	    
      int index(0);
      for (TrackCollection::const_iterator it = (*theTkCollection).begin(); 
	   it != (*theTkCollection).end(); 
	   ++it){
	
	if (index == idrec ) {	
	  tkt = &(*it);
	  break;
	}

	index++;
      }
	
      tcand++;	
      MuonRecTracks.push_back(tkt);
      MuonRecTracksIndex.push_back(idrec);

    } else {

      if (fVerbose) {
	if ( idrec < 0 ) cout << "==>bmmTracks4> Coulnd't find rec. track for tracker muon #" << tcnt << endl;  
      }
    }
    
    tcnt++;

  }
  
  if (fVerbose) cout << "==>bmmTracks4> Taking " << mcand << "/" << mcnt << " global and " 
		     << tcand  << "/" << tcnt <<  " tracker muons " <<  endl;
  
}


// ----------------------------------------------------------------------
void BSDumpCandidates::truthCandTracks(const edm::Event &iEvent, const edm::EventSetup& iSetup) {

  if (fVerbose) cout << endl << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>truthCandTracks> Sort muon and kaon candidates depending on ID of B-mother, event: " 
		     << fNevt << endl;

  if ( BmmRecTracks.size() > 1 ) {

    for ( unsigned int i=0; i<BmmRecTracks.size(); i++ ) { 

      for ( unsigned int j=i; j<BmmRecTracks.size(); j++ ) { 

	if (i == j) { continue; }

	if ( BmmRecTracksB[i] == BmmRecTracksB[j] ) {
	  
	  if ((*BmmRecTracks[j]).pt() > (*BmmRecTracks[i]).pt() ) { 
	    BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (BmmRecTracks[j],BmmRecTracks[i]));
	    BmmPairTracksIndex.push_back(pair<int, int>
						 (BmmRecTracksIndex[j],BmmRecTracksIndex[i]));
	  } else {
	    BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (BmmRecTracks[i],BmmRecTracks[j]));
	    BmmPairTracksIndex.push_back(pair<int, int>
						 (BmmRecTracksIndex[i],BmmRecTracksIndex[j]));
	  }

	}
      }
    }
  }  

  
  if ( JpsiRecTracks.size() > 1 ) {
    
    for ( unsigned int i=0; i<JpsiRecTracks.size(); i++ ) { 
      
      for ( unsigned int j=i; j<JpsiRecTracks.size(); j++ ) { 
	
	if (i == j) { continue; }
	
	if (JpsiRecTracksB[i] == JpsiRecTracksB[j]) {
	  

	  if ( KaonRecTracks.size() > 0 ) {

	    for ( unsigned int k=0; k<KaonRecTracks.size(); k++ ) {

	      if ( KaonRecTracksIndex[k] == JpsiRecTracksIndex[i] ||
		   KaonRecTracksIndex[k] == JpsiRecTracksIndex[j] ) { continue; }

	      if (JpsiRecTracksB[i] == KaonRecTracksB[k]) {
		
		if ((*JpsiRecTracks[j]).pt() > (*JpsiRecTracks[i]).pt() ) { 
		  JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
						   (JpsiRecTracks[j],JpsiRecTracks[i]));
		  JpsiPairTracksIndex.push_back(pair<int, int>
							(JpsiRecTracksIndex[j],JpsiRecTracksIndex[i]));
		} else {
		  JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
						   (JpsiRecTracks[i],JpsiRecTracks[j]));
		  JpsiPairTracksIndex.push_back(pair<int, int>
							(JpsiRecTracksIndex[i],JpsiRecTracksIndex[j]));
		}
		
		KaonTrack.push_back(KaonRecTracks[k]);
		KaonTrackIndex.push_back(KaonRecTracksIndex[k]);
	      }
	    }
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
void BSDumpCandidates::muonCandTracks(const edm::Event &iEvent, const edm::EventSetup& iSetup, double m_cand1, double m_cand2) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>muonCandTracks> Sort opposite charge muons into pairs depending on inv. mass ("
		     << m_cand1 << " GeV) or (" << m_cand2 << " GeV), event: " << fNevt << endl;

  int type(0);

  if ( MuonRecTracks.size() > 1 ) {

    for ( unsigned int i=0; i<MuonRecTracks.size(); i++ ) { 

      for ( unsigned int j=i+1; j<MuonRecTracks.size(); j++ ) {

	RecTracks.clear();
	RecTracksIndex.clear();

	RecTracks.push_back(MuonRecTracks[i]);
	RecTracks.push_back(MuonRecTracks[j]);

	RecTracksIndex.push_back(MuonRecTracksIndex[i]);
	RecTracksIndex.push_back(MuonRecTracksIndex[j]);

	type = massMuonCand(iEvent, iSetup, m_cand1, m_cand2);   //  < ------- here

	if ( type == 1 ) {

	  if ((*MuonRecTracks[j]).pt() > (*MuonRecTracks[i]).pt() ) { 
	    BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (MuonRecTracks[j],MuonRecTracks[i]));
	    BmmPairTracksIndex.push_back(pair<int, int>
						 (MuonRecTracksIndex[j],MuonRecTracksIndex[i]));
	  } else {
	    BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (MuonRecTracks[i],MuonRecTracks[j]));
	    BmmPairTracksIndex.push_back(pair<int, int>
						 (MuonRecTracksIndex[i],MuonRecTracksIndex[j]));
	  }
	}

	    
	if ( type == 2 ) {
	      
	  if ((*MuonRecTracks[j]).pt() > (*MuonRecTracks[i]).pt() ) { 
	    JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					     (MuonRecTracks[j],MuonRecTracks[i]));
	    JpsiPairTracksIndex.push_back(pair<int, int>
						  (MuonRecTracksIndex[j],MuonRecTracksIndex[i]));
	  } else {
	    JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					     (MuonRecTracks[i],MuonRecTracks[j]));
	    JpsiPairTracksIndex.push_back(pair<int, int>
						  (MuonRecTracksIndex[i],MuonRecTracksIndex[j]));
	  }
	}
	    
	if ( type == 3 ) {

	  if ((*MuonRecTracks[j]).pt() > (*MuonRecTracks[i]).pt() ) { 
	    BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (MuonRecTracks[j],MuonRecTracks[i]));
	    BmmPairTracksIndex.push_back(pair<int, int>
						 (MuonRecTracksIndex[j],MuonRecTracksIndex[i]));
	  } else {
	    BmmPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					    (MuonRecTracks[i],MuonRecTracks[j]));
	    BmmPairTracksIndex.push_back(pair<int, int>
						 (MuonRecTracksIndex[i],MuonRecTracksIndex[j]));
	  }
	      
	  if ((*MuonRecTracks[j]).pt() > (*MuonRecTracks[i]).pt() ) { 
	    JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					     (MuonRecTracks[j],MuonRecTracks[i]));
	    JpsiPairTracksIndex.push_back(pair<int, int>
						  (MuonRecTracksIndex[j],MuonRecTracksIndex[i]));
	  } else {
	    JpsiPairTracks.push_back(pair<const reco::Track*, const reco::Track*>
					     (MuonRecTracks[i],MuonRecTracks[j]));
	    JpsiPairTracksIndex.push_back(pair<int, int>
						  (MuonRecTracksIndex[i],MuonRecTracksIndex[j]));
	  }
	}
      }
    }
  }
}


// ----------------------------------------------------------------------
void BSDumpCandidates::kaonCandTracks(const edm::Event &iEvent, const edm::EventSetup& iSetup, double cone) { 
   
  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>kaonCandTracks> Find " << JpsiPairTracks.size() 
		     << " kaon candidate(s) with best chi2, from all tracks in cone "
		     << cone << " around reconstructed J/Psi , event: " << fNevt << endl;

  const reco::Track* tt = 0;

  int index(0), trk(-1), ncand(0), cand(-1);
  double chi2_min(99999.), chi2(-1.);

  if ( JpsiPairTracks.size() > 0 ) {
    
    for ( unsigned int i=0; i<JpsiPairTracks.size(); i++ ) {

      // -- Muons from J/Psi
      RecTracks.clear();
      RecTracksIndex.clear();

      RecTracks.push_back(JpsiPairTracks[i].first);
      RecTracksIndex.push_back(JpsiPairTracksIndex[i].first);

      RecTracks.push_back(JpsiPairTracks[i].second); 
      RecTracksIndex.push_back(JpsiPairTracksIndex[i].second);  

      // ------------------------------------------------------------------------
      // -- Find kaon for given muon pair (if no kaon cand. found, index remains -1)
      KaonTrack.push_back(tt);
      KaonTrackIndex.push_back(-1);

      index = 0; trk = -1; ncand = 0; cand = -1;

      SecVtxChi2.clear();
      InvMass.clear();

      for (TrackCollection::const_iterator it = (*theTkCollection).begin(); 
	   it != (*theTkCollection).end(); 
	   ++it){
	
	if (index == RecTracksIndex[0]) { index++;  continue; }
	if (index == RecTracksIndex[1]) { index++;  continue; }
	
	tt = &(*it);
	
	RecTracks.push_back(tt);
	RecTracksIndex.push_back(index);
     
	chi2 =  rmmKaonCand(iEvent, iSetup, cone);   //  < ------- here: chi2 = -1, if outside cone!

	RecTracks.pop_back();
	RecTracksIndex.pop_back();
	
	if (chi2 > 0) {

	  if (chi2 < chi2_min) {
	    
	    chi2_min  = chi2;
	    
	    KaonTrack.pop_back();
	    KaonTrackIndex.pop_back();
	    
	    KaonTrack.push_back(tt);
	    KaonTrackIndex.push_back(index);
	    
	    trk  = index;
	    cand = ncand;
	  }

	  ncand++;
	}
	
	index++;	
      }

      // ------------------------------------------------------------------------  

      if (fVerbose) {

	if (ncand) { 
	  
	  cout << "==>kaonCandTracks> Fount " << ncand << " kaon candidate tracks for J/Psi candidate #"
	       << i << ". Selected track #" << trk << ", with chi2 = " << chi2_min << "   (chi2=" 
	       << SecVtxChi2[cand] << ", m=" << InvMass[cand] << ")." << endl; 

	} else {	
	  
	  cout << "==>kaonCandTracks> No kaon candidates found for J/Psi candidate #" 
	       << i << "." << endl;
	}
      }
 
      // -- eff. histogram ---------------------------
      TH1D *h1 = (TH1D*)gBSFile->Get("eff");

      if ( ncand < 9 ) {
	
	h1->Fill(80.1 + ncand );
	
      } else {
	
	h1->Fill(89.1);
      }

      // ----------------------------------------------  
 
    }
  }
}
// ----------------------------------------------------------------------
int BSDumpCandidates::massMuonCand(const edm::Event &iEvent, const edm::EventSetup& iSetup, double m_cand1, double m_cand2) {  

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>massMuonCand> Looking for two muons in mass-windows around "
		     << m_cand1 << " and " << m_cand2 << ", event: " << fNevt << endl;

  TLorentzVector m0, m1, cand; 
  double mass(-1.);
  double dmass_1 = fMassRange;
  double dmass_2 = fMassRange2;

  int cnd1(0), cnd2(0);
  
  m0.SetXYZM((*RecTracks[0]).px(),
	     (*RecTracks[0]).py(),
	     (*RecTracks[0]).pz(),
	     0.1056583);
      
  m1.SetXYZM((*RecTracks[1]).px(),
	     (*RecTracks[1]).py(),
	     (*RecTracks[1]).pz(),
	     0.1056583);
    
  cand  = m0 + m1;
  mass  = cand.M();

  if ( m_cand1 > 0 && mass > (m_cand1 - dmass_1)  && mass < (m_cand1 + dmass_1) &&
       RecTracks[0]->charge() != RecTracks[1]->charge() ) { 


    if (fVerbose) cout << endl << "==>massCandidate> : Found B candidates with muons \t #"
		       << RecTracksIndex[0] << " \t #" <<  RecTracksIndex[1]
		       << "\t --->  Inv. Mass: "  << mass
		       << "( q1 = " << RecTracks[0]->charge()
		       << ", q2 = " << RecTracks[1]->charge()
		       << ")" << endl;
    cnd1 =  1;
  }
  
  if ( m_cand2 > 0 && mass > (m_cand2 - dmass_2)  && mass < (m_cand2 + dmass_2) &&
	      RecTracks[0]->charge() != RecTracks[1]->charge() ) { 

    if (fVerbose) cout << endl << "==>massCandidate> : Found J/Psi candidates with muons \t #"
		       << RecTracksIndex[0] << " \t #" <<  RecTracksIndex[1]
		       << "\t --->  Inv. Mass: "  << mass
		       << "( q1 = " << RecTracks[0]->charge()
		       << ", q2 = " << RecTracks[1]->charge()
		       << ")" << endl;
    cnd2 = 2;
  }

  if ( !cnd1 && !cnd2 ) {

    if (fVerbose) cout << endl << "==>massCandidate> : Muons not in any of the mass windows \t #"
		       << RecTracksIndex[0] << " \t #" <<  RecTracksIndex[1]
		       << "\t --->  Inv. Mass: "  << mass
		       << "( q1 = " << RecTracks[0]->charge()
		       << ", q2 = " << RecTracks[1]->charge()
		       << ")" << endl;
  }

  return cnd1+cnd2;
}


// ----------------------------------------------------------------------
double BSDumpCandidates::rmmKaonCand(const edm::Event &iEvent, const edm::EventSetup& iSetup, double cone) {  

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>rmmKaonCand> Looking for kaon in cone of "
		     << cone << " around J/Psi, event: " << fNevt << endl;

  double chi2 = kalmanVertexFit(iEvent, iSetup, 0, 3);  //  Vertex TYPE = 0, B+ -> mu mu K (all K-candidates)
  
  TLorentzVector m0, m1, kp;
  TLorentzVector jpsi, bs;
  
  m0.SetXYZM(RefittedTracks[0].px(),
	     RefittedTracks[0].py(),
	     RefittedTracks[0].pz(),
	     0.1056583);
  
  m1.SetXYZM(RefittedTracks[1].px(),
	     RefittedTracks[1].py(),
	     RefittedTracks[1].pz(),
	     0.1056583);
  
  kp.SetXYZM(RefittedTracks[2].px(),
	     RefittedTracks[2].py(),
	     RefittedTracks[2].pz(),
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
		       << RecTracksIndex[2] 	   << " with dr = " << dr << ", B mass " 
		       << bs.M() << " and chi2 " << chi2 << endl;
    
    
    InvMass.push_back(bs.M());
    SecVtxChi2.push_back(chi2);

    return chi2;
  
  } else {
    
    if (fVerbose) cout << "==>rmmKaonCand> !!! Kaon candidate rejected, ouside cone:  dr = " << dr
		       << ", B mass " << bs.M() << " and chi2 " << chi2 << " !!!" <<endl;

    return -1.;
  }
}


// ----------------------------------------------------------------------
int BSDumpCandidates::primaryVertex(const edm::Event &iEvent) {

  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>primaryVertex> List of primary vertices, event: " << fNevt << endl;

  int nvtx(0);

  TH1D *h1 = (TH1D*)gBSFile->Get("eff");

  // -- Primary vertex
  try {

    edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
    iEvent.getByLabel(fPrimaryVertexLabel.c_str(), recoPrimaryVertexCollection);
    const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());

    // const reco::Vertex* primVtx = 0;
    for(reco::VertexCollection::const_iterator v=recoPrimaryVertexCollection->begin(); 
	v!=recoPrimaryVertexCollection->end(); 
	++v){    
      
      nvtx++;

      if ( fVerbose > 0 ) {
	printf ("==>BSDumpStuff>  %i. Primary Vertex (x, y, z) = ( %5.4f, %5.4f, %5.4f)\n", 
		nvtx, v->x(), v->y(), v->z());
      }

      //     primVtx = &(*v);
      
      //     primVtx.SetX(v->x()*10);   //*10 to get mm (same unit as gen info)
      //     primVtx.SetY(v->y()*10);
      //     primVtx.SetZ(v->z()*10);
    }

    if ( nvtx == 0 ) {

      if ( fVerbose > 0 ) {
	cout << "==>BSDumpStuff>  no primary vertex in recoPrimaryVertexCollection" << endl;
      }

      return nvtx;
    }

    const reco::Vertex pV = vertices[0]; // ???? 

    if ( fVerbose ) cout << endl << "====> Primary Vertices from CTF Tracks: " << nvtx << " vertices found." << endl;
    
    if ( nvtx < 9 )  { h1->Fill(20.1 + nvtx); } else { h1->Fill(29.1); }
        
    if (nvtx > 0) {
      
      //    printf ("Taking Primary Vertex (x, y, z) = ( %5.4f, %5.4f, %5.4f)\n", pV.x(), pV.y(), pV.z());
      if ( fVerbose > 0 ) {
	printf ("Taking Primary Vertex (x, y, z)/10. = ( %5.4f, %5.4f, %5.4f)\n", 10*pV.x(), 10*pV.y(), 10*pV.z());
      }
      
      ChiSquared chi2(pV.chi2(),pV.ndof());
      
      TAnaVertex *pVtx;
      pVtx = new TAnaVertex();
      
      //    pVtx->setInfo(chi2.value(), chi2.degreesOfFreedom(), chi2.probability(), 0, 2); ??? change ndof to double ???
      pVtx->setInfo(chi2.value(), int(chi2.degreesOfFreedom()), chi2.probability(), 0, 2);
      pVtx->fPoint.SetXYZ(pV.position().x(),
			  pV.position().y(),
			  pV.position().z());
      
      primaryVertex2 = pV;
      gBSEvent->fPrimaryVertex2 = *pVtx;
    }
      
  } catch (cms::Exception &ex) {
   
    //    cout << ex.explainSelf() << endl;
    
    if ( fVerbose > 0 ) {
      cout << "==>primaryVertex> primaryVertex " << fPrimaryVertexLabel.c_str() << " not found " << endl;  
    }
  
    TAnaVertex *pVtx;
    pVtx = new TAnaVertex();

    pVtx->setInfo(-1., -1, -1., -1, -1);
    pVtx->fPoint.SetXYZ(-100., -100., -100.);

    gBSEvent->fPrimaryVertex2 = *pVtx;
  } 

  return nvtx;
}

// ----------------------------------------------------------------------
void BSDumpCandidates::secondaryVertex(const edm::Event &iEvent, const edm::EventSetup& iSetup) { 
  
  if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
  if (fVerbose) cout << "==>secondaryVertex> Fitting secondary vertices, event: " << fNevt << endl;

  if (fVerbose) {

    cout << "==>secondaryVertex> " << BmmPairTracks.size() << " Bs -> mu+ mu- pair(s), event: " 
	 << fNevt << ", sel = " << fBmmSel << endl; 
    cout << "==>secondaryVertex> " << JpsiPairTracks.size() << " J/Psi -> mu+ mu- pair(s), event: " 
	 << fNevt << ", sel = " << fBmmSel << endl; 
    cout << "==>secondaryVertex> " << KaonTrack.size() << " Kaon(s) for B+ -> J/Psi K+, event: " 
	 << fNevt << ", sel = " << fBmmSel << endl; 
  }

  double chi2(-1.);

  int nbs0  = BmmPairTracks.size();
  int njpsi = JpsiPairTracks.size();
  int nkaon = KaonTrack.size();

  if ( BmmPairTracks.size() > 0 ) {
    
    for ( unsigned int i=0; i<BmmPairTracks.size(); i++ ) { 
      
      RecTracks.clear();
      RecTracksIndex.clear();

      RecTracks.push_back(BmmPairTracks[i].first);
      RecTracks.push_back(BmmPairTracks[i].second);

      RecTracksIndex.push_back(BmmPairTracksIndex[i].first);
      RecTracksIndex.push_back(BmmPairTracksIndex[i].second);

      chi2 = kalmanVertexFit(iEvent, iSetup, 5310+fBmmSel, 2);  // Vertex TYPE = 531x  
                                                                        // -> mu mu or BG from B, 2 TRACKS => Cand 1
    }
  }
 
  if ( JpsiPairTracks.size() > 0 ) {
    
    if ( JpsiPairTracks.size() == KaonTrack.size() ) {

      for ( unsigned int i=0; i<JpsiPairTracks.size(); i++ ) { 
	
	if ( KaonTrackIndex[i] > 0 ) {
	  
	  RecTracks.clear();
	  RecTracksIndex.clear();
	  
	  RecTracks.push_back(JpsiPairTracks[i].first);
	  RecTracksIndex.push_back(JpsiPairTracksIndex[i].first);
	  
	  RecTracks.push_back(JpsiPairTracks[i].second);
	  RecTracksIndex.push_back(JpsiPairTracksIndex[i].second);
	  
	  chi2 = kalmanVertexFit(iEvent, iSetup, 4430+fBmmSel , 2);  // Vertex TYPE = 443x 
	                                                                     // -> mu mu from J/Psi, 2 TRACKS => Cand 2
	  
	  RecTracks.push_back(KaonTrack[i]);
	  RecTracksIndex.push_back(KaonTrackIndex[i]);
	  
	  chi2 = kalmanVertexFit(iEvent, iSetup, 5210+fBmmSel, 3);  // Vertex TYPE = 521x  
	                                                                    // -> mu mu K, 3 TRACKS => Cand 3
	}
      }
    } else {

      cout << "****ERROR: Number of kaon candidates: " << KaonTrack.size()
	   << ", is not equal to number of J/Psi candidates: " <<  JpsiPairTracks.size() << "!!!! ****" <<endl;
      nkaon = 9;
      njpsi = 9;
    }
  }

  TH1D *h1 = (TH1D*)gBSFile->Get("eff");
  if (nbs0  >= 0 && nbs0  < 10) h1->Fill(100.1 + (fBmmSel-1)*100 + nbs0);
  if (njpsi >= 0 && njpsi < 10) h1->Fill(120.1 + (fBmmSel-1)*100 + njpsi);
  if (nkaon >= 0 && nkaon < 10) h1->Fill(140.1 + (fBmmSel-1)*100 + nkaon);
}


// ----------------------------------------------------------------------
double BSDumpCandidates::kalmanVertexFit(const edm::Event &iEvent, const edm::EventSetup& iSetup
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

  if (fVerbose) cout << "==>kalmanVertexFit> Vertexing with " << RecTracks.size() << " reco. tracks." << endl;
   
  //get transient track 
  for ( unsigned int i=0; i<RecTracks.size(); i++ ) { 

    RecoTransientTrack.push_back( (theB->build( &(*RecTracks[i]) )) );
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
void BSDumpCandidates::fillVertex(const edm::Event &iEvent, const edm::EventSetup& iSetup
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

  RefittedTracks.clear();

  for(vector<reco::TransientTrack>::const_iterator i = refTT.begin(); 
      i != refTT.end(); i++) {

    const Track & ftt = i->track();
    // if (fVerbose) cout << "mmmmmmm : " << ftt->innermomentum() << endl;
    refitted.push_back(ftt);
    RefittedTracks.push_back(ftt);
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

  TH1D *h1 = (TH1D*)gBSFile->Get("eff");

  h1->Fill(30.1);

  // ==============================================================================================
  // ----------------------------- REFITTED TRACKS ------------------------------------------------

  if ( refitted.size() == ntracks ) {

    h1->Fill(31.1);

    if ( ntracks == 2 ) {
      
      if (fVerbose) cout << "--- Vertex with number of refitted tracks: " << refitted.size() << " of 2 ---" << endl;
      
      i1 = RecTracksIndex[0];
      i2 = RecTracksIndex[1];

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
      
      i1 = RecTracksIndex[0];
      i2 = RecTracksIndex[2];  // This should be the K+

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
		 fMassKaon);
      
      bs = m0 + m1 + kp;
      
      mass = bs.M();     
    }
  }
  
  // ----------------------------- ORIGINAL TRACKS --------------------------------------------------
  
  else if ( original.size() == ntracks ) {

    if (fVerbose) cout << "--- Not enough refitted tracks for this Vertex: " << refitted.size()  << " ---" << endl;
    if (fVerbose) cout << "--- Use original tracks instead ---" << endl;
    
    h1->Fill(32.1);

    if ( ntracks == 2 ) {
      
      if (fVerbose) cout << "--- Vertex with number of original tracks: " << original.size() << " of 2 ---" << endl;
      
      i1 = RecTracksIndex[0];
      i2 = RecTracksIndex[1];
      
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
      
      i1 = RecTracksIndex[0];
      i2 = RecTracksIndex[2];  // This should be the K+

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

    h1->Fill(33.1);
    if (fVerbose) cout << "!!! Not enough tracks for this Vertex: " << refitted.size() 
		       << " refitted and " << original.size() << "original tracks !!!" << endl;
  }
  
  // ---------------------------------------------------------------------------------------------
  // ==============================================================================================

  if ( type > -1 ) {

    TH1D *h1 = (TH1D*)gBSFile->Get(Form("m%.3i_%i", int(type/10), fBmmSel));
    h1->Fill(mass);
  }

  ChiSquared chi(v->totalChiSquared(), v->degreesOfFreedom());

  VertexDistanceXY axy;
  double dXY      = axy.distance(primaryVertex2, *v).value();
  double dXYE     = axy.distance(primaryVertex2, *v).error();
  double compXY   = axy.compatibility(primaryVertex2, *v);

  VertexDistance3D a3d;
  double d3d      = a3d.distance(primaryVertex2, *v).value();
  double d3dE     = a3d.distance(primaryVertex2, *v).error();
  double comp3d   = a3d.compatibility(primaryVertex2, *v);


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
    TAnaCand  *pCand = gBSEvent->addCand();
    
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

//**   reco::VertexRecoToSimCollection recSimCollVertex = (*recSimCollectionVertex); 


  // -------------------------------not used---------------------------------------------------
//   reco::RecoToSimCollection recSimColl             = (*recSimCollection); 

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
void BSDumpCandidates::decayChannel(const char *fileName) {

  // -- sel = 1 settings
  fPrintChannel = "N/D";              fPrintChannel2 = "N/D";
  fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = -1;
  fTruthMC_mom = 531;                 fTruthMC2_mom = -1;
  fTruthMC_gmo = -1;                  fTruthMC2_gmo = -1;

  // -- sel = 2, 3 settings
  fMass = 5.367;                      fMass2 = -1;
  fMassRange = 5.;                    fMassRange2 = 0.2; 
  fMassKaon = 0.493677;


  // ****************************************** Bd decays *********************************************************
  if ( !strcmp("Bd2PiKp", fileName) ) { 

    fTruthMC_I = 211;  fTruthMC_II = 321;
    fTruthMC_mom = 511;
    fMass = 5.367; //    fMass = 5.279;

    fPrintChannel = "Bd -> pi- K+";

    cout << "Selected decay mode: Bd (" << fTruthMC_mom << ") to K+ (" 	 << fTruthMC_II 
	 << ") pi- (" << fTruthMC_I << ")" << endl;

  } else if ( !strcmp("Bd2PiMuNu", fileName) ) {

    fTruthMC_I = 211;  fTruthMC_II = 13;
    fTruthMC_mom = 511;
    fMass = 5.367; //    fMass = 5.279;

    fPrintChannel = "Bd -> pi- mu+ nu_mu";

    cout << "Selected decay mode: Bd (" << fTruthMC_mom << ") to mu+ (" 	 << fTruthMC_II 
	 << ") pi- (" << fTruthMC_I << ") nu_mu" << endl;

  } else if ( !strcmp("Bd2PiPi", fileName) ) {

    fTruthMC_I = 211;  fTruthMC_II = 211;
    fTruthMC_mom = 511;
    fMass = 5.367; //    fMass = 5.279;

    fPrintChannel = "Bd -> pi+ pi-";

    cout << "Selected decay mode: Bd (" << fTruthMC_mom << ") to pi+ pi- (" << fTruthMC_I << ")" << endl;

  } else if ( !strcmp("Bd2MuMu", fileName) ) {

    fTruthMC_I = 13;  fTruthMC_II = 13;
    fTruthMC_mom = 511;
    fMass = 5.367; //    fMass = 5.279;

    fPrintChannel = "Bd -> mu+ mu-";

    cout << "Selected decay mode: Bd (" << fTruthMC_mom << ") to mu+ mu- (" << fTruthMC_I << ")" << endl;

  } else if ( !strcmp("Bd2MuMuPi0", fileName) ) {

    fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = 111;
    fTruthMC_mom = 511;                 fTruthMC2_mom = 511;
    fTruthMC_gmo = -1;                  fTruthMC2_gmo = -1;
    fMass = 5.367;  /* fMass = 5.279;*/ fMass2 = -1;
    fMassKaon = 0.1349766;

    fPrintChannel = "Bd -> mu+ mu-";
    fPrintChannel2 = "Bd -> mu+ mu- pi0 ";

    cout << "Selected decay modes: Bd (" << fTruthMC_mom << ") to mu+ mu- (" << fTruthMC_I << ")" << endl
	 << "     and Bd (" << fTruthMC2_mom  << ") to mu+ mu- (" << fTruthMC_I << ") and pi0 (" 
	 << fTruthMC2  << ")" << endl;

  // ****************************************** Bs decays *********************************************************

  }  else if ( !strcmp("Bs2KmMuNu", fileName) ) {   

    fTruthMC_I = 321;  fTruthMC_II = 13;
    fTruthMC_mom = 531;
    fMass = 5.367;

    fPrintChannel = "Bs -> K- mu+ nu_mu";

    cout << "Selected decay mode: Bs (" << fTruthMC_mom << ") to mu+ (" 	 << fTruthMC_II 
	 << ") K- (" << fTruthMC_I << ") nu_mu" << endl;

  } else if ( !strcmp("Bs2KmPi", fileName) ) {

    fTruthMC_I = 321;  fTruthMC_II = 211;
    fTruthMC_mom = 531;
    fMass = 5.367;

    fPrintChannel = "Bs -> K- pi+";

    cout << "Selected decay mode: Bs (" << fTruthMC_mom << ") to pi+ (" 	 << fTruthMC_II 
	 << ") K- (" << fTruthMC_I << ")" << endl;

  } else if ( !strcmp("Bs2KpKm", fileName) ) {

    fTruthMC_I = 321;  fTruthMC_II = 321;
    fTruthMC_mom = 531;
    fMass = 5.367;

    fPrintChannel = "Bs -> K+ K-";

    cout << "Selected decay mode: Bs (" << fTruthMC_mom << ") to K+ K- (" << fTruthMC_I << ")" << endl;

  }  else if ( !strcmp("Bs2PiPi", fileName) ) {

    fTruthMC_I = 211;  fTruthMC_II = 211;
    fTruthMC_mom = 531;
    fMass = 5.367;

    fPrintChannel = "Bs -> pi+ pi-";

    cout << "Selected decay mode: Bs (" << fTruthMC_mom << ") to pi+ pi- (" << fTruthMC_I << ")" << endl;

  } else if ( !strcmp("Bs2MuMuGamma", fileName) ) {

    fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = 22;
    fTruthMC_mom = 531;                 fTruthMC2_mom = 531;
    fTruthMC_gmo = -1;                  fTruthMC2_gmo = -1;
    fMass = 5.367;                      fMass2 = -1;
    fMassKaon = 0.;

    fPrintChannel = "Bs -> mu+ mu-";
    fPrintChannel2 = "Bs -> mu+ mu- gamma ";

    cout << "Selected decay modes: Bs (" << fTruthMC_mom << ") to mu+ mu- (" << fTruthMC_I << ")" << endl
	 << "     and Bs (" << fTruthMC2_mom  << ") to mu+ mu- (" << fTruthMC_I << ") and gamma (" 
	 << fTruthMC2  << ")" << endl;

  // ****************************************** Bu decays *********************************************************

  } else if ( !strcmp("Bu2MuMuMuNu", fileName) ) {

    fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = 13;
    fTruthMC_mom = 521;                 fTruthMC2_mom = 521;
    fTruthMC_gmo = -1;                  fTruthMC2_gmo = -1;
    fMass = 5.367; /* fMass = 5.279; */ fMass2 = -1;
    fMassKaon = 0.1056583;

    fPrintChannel = "Bu -> mu+ mu-";
    fPrintChannel2 = "Bu -> mu+ mu- mu+ nu_mu ";

    cout << "Selected decay modes: Bu (" << fTruthMC_mom << ") to mu+ mu- (" << fTruthMC_I << ")" << endl
	 << "     and Bu (" << fTruthMC2_mom  << ") to mu+ mu- (" << fTruthMC_I << ") and mu+ (" 
	 << fTruthMC2  << ") nu_mu" << endl;

  // ****************************************** Bc decays *********************************************************

  } else if ( !strcmp("Bc2MuMuMuNu", fileName) ) {

    fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = 13;
    fTruthMC_mom = 521;                 fTruthMC2_mom = 521;
    fTruthMC_gmo = -1;                  fTruthMC2_gmo = -1;
    fMass = 5.367; /* fMass = 6.286; */ fMass2 = -1;
    fMassKaon = 0.1056583;

    fPrintChannel = "Bu -> mu+ mu-";
    fPrintChannel2 = "Bu -> mu+ mu- mu+ nu_mu ";

    cout << "Selected decay modes: Bc (" << fTruthMC_mom << ") to mu+ mu- (" << fTruthMC_I << ")" << endl
	 << "     and Bc (" << fTruthMC2_mom  << ") to mu+ mu- (" << fTruthMC_I << ") and mu+ (" 
	 << fTruthMC2  << ") nu_mu" << endl;

  } else if ( !strcmp("Bc2JpsiMuNu", fileName) ){
    
    fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = 13;
    fTruthMC_mom = -1;                  fTruthMC2_mom = 443;
    fTruthMC_gmo = -1;                  fTruthMC2_gmo = 521;
    fMass = 5.367; /* fMass = 6.286;*/  fMass2 = 3.096;
    fMassKaon = 0.1056583;

    fPrintChannel  = "Bc -> J/Psi mu+ nu_mu";
    fPrintChannel2 = "Bc -> J/Psi mu+ nu_mu";

    cout << "Selected decay modes: Bc (" << fTruthMC_mom << ") to mu+ mu- (" << fTruthMC_I << ")" << endl
	 << "     and Bc (" << fTruthMC2_gmo << ") to J/Psi (" << fTruthMC2_mom << ") mu+ (" << fTruthMC2
	 << ") where J/Psi to mu+ mu- (" << fTruthMC_I  << ")" << endl;

  // *************************************** Lambda_b decays ******************************************************

  } else if ( !strcmp("Lb2PKm", fileName) ) { 

    fTruthMC_I = 321;  fTruthMC_II = 2212;
    fTruthMC_mom = 5122; 
    fMass = 5.367; //    fMass = 5.624;                 

    fPrintChannel = "Lb -> p+ K-";

    cout << "Selected decay mode: Lb (" << fTruthMC_mom << ") to p+ (" 	 << fTruthMC_II 
	 << ") K- (" << fTruthMC_I << ")" << endl;

  } else if ( !strcmp("Lb2PPi", fileName) ) {

    fTruthMC_I = 211;  fTruthMC_II = 2212;
    fTruthMC_mom = 5122; 
    fMass = 5.367; //    fMass = 5.624;                 

    fPrintChannel = "Lb -> p+ pi-";

    cout << "Selected decay mode: Lb (" << fTruthMC_mom << ") to p+ (" 	 << fTruthMC_II 
	 << ") pi- (" << fTruthMC_I << ")" << endl;

  // ******************************* Signal & normalization channel *************************************************
      
  }  else {  

    fTruthMC_I = 13;  fTruthMC_II = 13; fTruthMC2 = 321;
    fTruthMC_mom = 531;                 fTruthMC2_mom = 443;
    fTruthMC_gmo = -1;                  fTruthMC2_gmo = 521;
    fMass = 5.367;                      fMass2 = 3.096;

    fPrintChannel = "Bs -> mu+ mu-";
    fPrintChannel2 = "B+ -> J/Psi K+";

    cout << "Selected decay modes: Bs (" << fTruthMC_mom << ") to mu+ mu- (" << fTruthMC_I << ")" << endl
	 << "     and B+ (" << fTruthMC2_gmo << ") to J/Psi (" << fTruthMC2_mom << ") K+ (" << fTruthMC2
	 << ") where J/Psi to mu+ mu- (" << fTruthMC_I  << ")" << endl;
  }          

  // --------------------------------------------------------------------------------------------------------------
  cout << "   ---> Channels (" << fChannel << "): " << fPrintChannel << ", " << fPrintChannel2 << endl;
  cout << "   ---> Mass windows: " << fMass << " +/- " << fMassRange 
       << ", " << fMass2  << " +/- " << fMassRange2 << endl;
  // --------------------------------------------------------------------------------------------------------------

}

// ----------------------------------------------------------------------
void BSDumpCandidates::clearTracks() {

  BmmRecTracks.clear();
  BmmRecTracksIndex.clear();
  BmmRecTracksB.clear();

  JpsiRecTracks.clear();
  JpsiRecTracksIndex.clear();
  JpsiRecTracksB.clear();

  KaonRecTracks.clear();
  KaonRecTracksIndex.clear();
  KaonRecTracksB.clear();

  MuonRecTracks.clear();
  MuonRecTracksIndex.clear(); 
}

// ----------------------------------------------------------------------
void BSDumpCandidates::clearCandidateTracks() {


  BmmPairTracks.clear();
  JpsiPairTracks.clear();
  KaonTrack.clear();

  BmmPairTracksIndex.clear();
  JpsiPairTracksIndex.clear();
  KaonTrackIndex.clear();
}



// ------------ method called once each job just before starting event loop  ------------
void BSDumpCandidates::beginJob(const edm::EventSetup& setup) {  
 
   gBSFile->cd();

   edm::ESHandle<MagneticField> theMF;
   setup.get<IdealMagneticFieldRecord>().get(theMF);
   
   edm::ESHandle<TrackAssociatorBase> theAssociator;
   setup.get<TrackAssociatorRecord>().get(fAssociatorLabel.c_str(), theAssociator);
   fAssociator = (TrackAssociatorBase*)theAssociator.product();
   
//    edm::ESHandle<VertexAssociatorBase> theTracksAssociator;
//    setup.get<VertexAssociatorRecord>().get(fVtxAssociatorLabel.c_str(), theTracksAssociator);
//    fVtxAssociator = (VertexAssociatorBase *) theTracksAssociator.product();

}

// ------------ method called once each job just after ending the event loop  ------------
void BSDumpCandidates::endJob() {  

}

//define this as a plug-in
// DEFINE_FWK_MODULE(BSDumpCandidates);
