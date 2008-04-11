#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpTracks.h"

#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" 
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h" 
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include <TFile.h>
#include <TH1.h>

// -- Yikes!
extern TAna00Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;


// ----------------------------------------------------------------------
HFDumpTracks::HFDumpTracks(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fGenEventLabel(iConfig.getUntrackedParameter<string>("generatorEventLabel", string("source"))),
  fSimTracksLabel(iConfig.getUntrackedParameter<string>("simTracksLabel", string("famosSimHits"))),
  fTrackingParticlesLabel(iConfig.getUntrackedParameter<string>("trackingParticlesLabel", string("trackingParticles"))),
  fTracksLabel(iConfig.getUntrackedParameter<string>("tracksLabel", string("ctfWithMaterialTracks"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fAssociatorLabel(iConfig.getUntrackedParameter<string>("associatorLabel", string("TrackAssociatorByChi2"))), 
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)) {

  using namespace std;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTracks constructor" << endl;
  cout << "--- Verbose                : " << fVerbose << endl;
  cout << "--- generatorEventLabel    : " << fGenEventLabel.c_str() << endl;
  cout << "--- simTracksLabel         : " << fSimTracksLabel.c_str() << endl;
  cout << "--- trackingParticlesLabel : " << fTrackingParticlesLabel.c_str() << endl;
  cout << "--- tracksLabel            : " << fTracksLabel.c_str() << endl;
  cout << "--- muonsLabel             : " << fMuonsLabel << endl;
  cout << "--- associatorLabel        : " << fAssociatorLabel.c_str() << endl;
  cout << "--- doTruthMatching        : " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS
  cout << "----------------------------------------------------------------------" << endl;

  fNevt = 0;
}


// ----------------------------------------------------------------------
HFDumpTracks::~HFDumpTracks() {
  
}


// ----------------------------------------------------------------------
void HFDumpTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  fNevt++; 
  // -- get the collection of tracks
  edm::Handle<TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks);  

  if (fVerbose > 0) cout << "==>HFDumpTracks> nTracks = " << tracks->size() 
			 << ", event: " << fNevt << endl;

  // -- get the collection of muons
  Handle<MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);

  if (fVerbose > 0) cout << "==>HFDumpTracks> nMuons = " << hMuons->size() << endl;

  // -- store their muon track indices
  vector<int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    TrackRef track = muon->track();
    muonIndices.push_back((muon->track()).index());
  }
 
  // -- get the tracking particle collection needed for truth matching. Only on RECO data tier!
  RecoToSimCollection recSimColl;
  if (1 == fDoTruthMatching) {
  if (fVerbose > 0) cout << "==>HFDumpTracks> Get tracking particles for TrackAssociator" << endl;
    try {
      edm::Handle<TrackingParticleCollection> trackingParticles;
      iEvent.getByLabel(fTrackingParticlesLabel.c_str(), trackingParticles);
      recSimColl = fAssociator->associateRecoToSim(tracks, trackingParticles, &iEvent); 
    } catch (cms::Exception &ex) {
      cout << ex.explainSelf() << endl;
    }
  }

  // -- Get the stuff needed for FAMOS truth matching
  Handle<HepMCProduct> hepmc;
  const HepMC::GenEvent *genEvent = 0;
  edm::Handle<std::vector<SimTrack> > simTracks;
  if (2 == fDoTruthMatching) {
  if (fVerbose > 0) cout << "==>HFDumpTracks> Get sim. track for FAMOS truth matching" << endl;
    iEvent.getByLabel(fGenEventLabel.c_str(), hepmc);
    genEvent = hepmc->GetEvent();
    iEvent.getByLabel(fSimTracksLabel.c_str(), simTracks); 
  }      

  TAnaTrack *pTrack; 
//   TH1D *h1 = (TH1D*)gHFFile->Get("h1");
//   TH1D *h2 = (TH1D*)gHFFile->Get("h2");
//   h1->Fill(tracks->size());
//   h2->Fill(tracks->size());

  for (unsigned int i = 0; i < tracks->size(); ++i){

    TrackRef rTrack(tracks, i);
    Track track(*rTrack);

    pTrack = gHFEvent->addRecTrack();
    pTrack->fIndex = rTrack.index();
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
    pTrack->fMuID = 0.; 

    for (unsigned int im = 0; im < muonIndices.size(); ++im) {
      if (i == muonIndices[im]) {
	pTrack->fMuID = 1.; 
      }
    }

    int gen_pdg_id(-99999), gen_id(-99999), gen_cnt(0);

    // -- RECO truth matching with TrackingParticle
    if (1 == fDoTruthMatching) {
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[rTrack];
	TrackingParticleRef tpr = tp.begin()->first;  
	const HepMC::GenParticle *genPar = 0; 
	for (TrackingParticle::genp_iterator pit = tpr->genParticle_begin(); 
	     pit != tpr->genParticle_end(); 
	     ++pit){
	  genPar     = pit->get();
	  gen_pdg_id = (*genPar).pdg_id();
	  gen_id     = (*genPar).barcode()-1;
	  //gen_pt     = (*genPar).momentum().perp();
	  //gen_phi    = (*genPar).momentum().phi();
	  //gen_eta    = (*genPar).momentum().pseudoRapidity();
	  gen_cnt++;
	}
      } catch (Exception event) {
	// 	cout << "%%>   Rec. Track #" << setw(2) << rTrack.index() << " pT: " 
	// 	     << setprecision(2) << setw(6) << track.pt() 
	// 	     <<  " matched to 0 MC Tracks" << endl;
      }
    }

    // -- FAMOS truth matching via SimHit
    if (2 == fDoTruthMatching) {
      for (trackingRecHit_iterator it = track.recHitsBegin();  it != track.recHitsEnd(); it++) {
	if ((*it)->isValid()) {

	  int currentId(-1);
	  if (const SiTrackerGSRecHit2D *rechit = dynamic_cast<const SiTrackerGSRecHit2D *> (it->get())) {
            currentId = rechit->simtrackId();          
	  }
	  
	  for (SimTrackContainer::const_iterator simTrack = simTracks->begin(); 
	       simTrack != simTracks->end(); 
	       simTrack++)   { 

	    if (simTrack->trackId() == currentId) {
	      int igen = simTrack->genpartIndex();
	      HepMC::GenParticle *genPar = genEvent->barcode_to_particle(igen);
	      if (genPar) {
		gen_pdg_id = (*genPar).pdg_id();
		gen_id     = (*genPar).barcode()-1;
		goto done;
	      }
	    }
	  }
	}
      }
    done:;
    }            
      
    pTrack->fGenIndex = gen_id;
    pTrack->fMCID     = gen_pdg_id;
    if (fVerbose > 0) {
      cout << "%%>   Rec. Track #" << setw(2)  << i << ": "; 
      pTrack->dump(); 
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpTracks::beginJob(const edm::EventSetup& setup) {
  edm::ESHandle<TrackAssociatorBase> theAssociator;
  setup.get<TrackAssociatorRecord>().get(fAssociatorLabel.c_str(), theAssociator);
  fAssociator = (TrackAssociatorBase*)theAssociator.product();

  gHFFile->cd();

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpTracks::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(HFDumpTracks);
