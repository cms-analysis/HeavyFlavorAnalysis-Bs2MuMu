#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpTracks.h"

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

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"


#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include <TFile.h>
#include <TH1.h>

// -- Yikes!
extern TAna00Event *gBSEvent;
extern TFile       *gBSFile;

using namespace std;
using namespace edm;
using namespace reco;
using namespace l1extra ;



// ----------------------------------------------------------------------
BSDumpTracks::BSDumpTracks(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fGenEventLabel(iConfig.getUntrackedParameter<string>("generatorEventLabel", string("source"))),
  fSimTracksLabel(iConfig.getUntrackedParameter<string>("simTracksLabel", string("famosSimHits"))),
  fTrackingParticlesLabel(iConfig.getUntrackedParameter<string>("trackingParticlesLabel", string("trackingParticles"))),
  fTracksLabel(iConfig.getUntrackedParameter<string>("tracksLabel", string("ctfWithMaterialTracks"))),
  fAssociatorLabel(iConfig.getUntrackedParameter<string>("associatorLabel", string("TrackAssociatorByChi2"))), 
  fL1MuLabel(iConfig.getUntrackedParameter<string> ("L1MuLabel")),
  fMuonsLabel1(iConfig.getUntrackedParameter<InputTag>("muonsLabel1")),
  fMuonsLabel2(iConfig.getUntrackedParameter<InputTag>("muonsLabel2")),
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)) {

  using namespace std;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- BSDumpTracks constructor" << endl;
  cout << "--- Verbose                : " << fVerbose << endl;
  cout << "--- generatorEventLabel    : " << fGenEventLabel.c_str() << endl;
  cout << "--- simTracksLabel         : " << fSimTracksLabel.c_str() << endl;
  cout << "--- trackingParticlesLabel : " << fTrackingParticlesLabel.c_str() << endl;
  cout << "--- tracksLabel            : " << fTracksLabel.c_str() << endl;
  cout << "--- muonsLabel1            : " << fMuonsLabel1 << endl;
  cout << "--- muonsLabel2            : " << fMuonsLabel2 << endl;
  cout << "--- L1MuonLabel            : " << fL1MuLabel << endl;
  cout << "--- associatorLabel        : " << fAssociatorLabel.c_str() << endl;
  cout << "--- doTruthMatching        : " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS
  cout << "----------------------------------------------------------------------" << endl;

  fNevt = 0;
}


// ----------------------------------------------------------------------
BSDumpTracks::~BSDumpTracks() {
  
}


// ----------------------------------------------------------------------
void BSDumpTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  fNevt++; 

  using reco::TrackCollection;

  theTkCollection = 0;

  // -- get the collection of tracks
  edm::Handle<TrackCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks);
  theTkCollection  = tracks.product();  

  if (fVerbose > 0) cout << "==>BSDumpTracks> nTracks = " << tracks->size() 
			 << ", event: " << fNevt << endl;

  // -------- get the collection of global muons -----------
  vector<int> muonIndices;

  try {

    Handle<MuonCollection> hMuons;
    iEvent.getByLabel(fMuonsLabel1, hMuons);
    
    if (fVerbose > 0) cout << "==>BSDumpTracks> globalMuons = " << hMuons->size() << endl;
    
        for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
      TrackRef track = muon->track();
      muonIndices.push_back((muon->track()).index());
    }
    
  } catch (Exception event) {
    cout << "%% -- No MuonCollection with label " << fMuonsLabel1 << endl;
  }

   // ------- get the collection of tracker muons -----------
  vector<int> tkMuonIndices;

  try {

    Handle<MuonCollection> tkMuons;
    iEvent.getByLabel(fMuonsLabel2, tkMuons);
    
    if (fVerbose > 0) cout << "==>BSDumpTracks> trackerMuons = " << tkMuons->size() << endl;
    
    for (MuonCollection::const_iterator tkmu = tkMuons->begin(); tkmu != tkMuons->end(); ++tkmu) {
      TrackRef track = tkmu->track();
      tkMuonIndices.push_back((tkmu->track()).index());
    }
    
  } catch (Exception event) {
    cout << "%% -- No MuonCollection with label " << fMuonsLabel2 << endl;
  }

  // ------- get the collection of L1 muons -----------
  vector<int> l1MuonIndices;

  try {

    edm::Handle<l1extra::L1MuonParticleCollection> l1extmu;
    iEvent.getByLabel(fL1MuLabel,l1extmu);
    
    const l1extra::L1MuonParticleCollection& L1ExtMu = *l1extmu;
    
    if ( fVerbose ) cout << "==>BSDumpTracks> L1Muons = " << l1extmu->size() << endl;
    
    
    if (&L1ExtMu) {
      
      for (l1extra::L1MuonParticleCollection::const_iterator muItr = L1ExtMu.begin(); muItr != L1ExtMu.end(); ++muItr) {
	
	int idrec = idRecTrack(muItr->pt(), muItr->eta(), muItr->phi(), 100., 0.4, 0.9);
	l1MuonIndices.push_back(idrec);
      }
    }
    
  } catch (Exception event) {
    cout << "%% -- No MuonCollection with label " << fL1MuLabel << endl;
  }
  
  // -- get the tracking particle collection needed for truth matching. Only on RECO data tier!
  RecoToSimCollection recSimColl;
  if (1 == fDoTruthMatching) {
    if (fVerbose > 0) cout << "==>BSDumpTracks> Get tracking particles for TrackAssociator" << endl;
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
    
    if (fVerbose > 0) cout << "==>BSDumpTracks> Get sim. track for FAMOS truth matching" << endl;
    
    try { 
      iEvent.getByLabel(fGenEventLabel.c_str(), hepmc);
      genEvent = hepmc->GetEvent();
    } catch (Exception event) {
      cout << "%% -- No HepMCProduct with label " << fGenEventLabel.c_str() << endl;
    }
    
    try {
      iEvent.getByLabel(fSimTracksLabel.c_str(), simTracks); 
    } catch (Exception event) {
      cout << "%% -- No SimTrack  with label " << fSimTracksLabel.c_str() << endl;
    }
  }   
  
  TAnaTrack *pTrack; 

  for (unsigned int i = 0; i < tracks->size(); ++i){

    TrackRef rTrack(tracks, i);
    Track track(*rTrack);

    pTrack = gBSEvent->addRecTrack();
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
    pTrack->fMuType = 0;  
    pTrack->fMuID = -1.;  // -- index in global muons (type = 33)
    pTrack->fKaID = -1.;  // -- index in tracker muons (type = 42)
    pTrack->fElID = -1.;  // -- index in l1 muons (type = 11)

    for (unsigned int im = 0; im < muonIndices.size(); ++im) {
      if (i == muonIndices[im]) {
	pTrack->fMuID = im; 
	pTrack->fMuType |= (0x1 << 0);
      }
    }

    for (unsigned int itm = 0; itm < tkMuonIndices.size(); ++itm) {
      if (i == tkMuonIndices[itm]) {
	pTrack->fKaID = itm; 
	pTrack->fMuType |= (0x1 << 1);
      }
    }

    for (unsigned int ilm = 0; ilm < l1MuonIndices.size(); ++ilm) {
      if (i == l1MuonIndices[ilm]) {
	pTrack->fElID = ilm; 
	pTrack->fMuType |= (0x1 << 2);
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


// ----------------------------------------------------------------------
int BSDumpTracks::idRecTrack(double pt, double eta, double phi, double ept, double eeta, double ephi) {
  
  int found(-1), index(0);

  double dpt(0.),  dphi(0.),  deta(0.);
  double mdpt(9999.), mdphi(9999.), mdeta(9999.);
  // double ept(0.2), ephi(0.01), eeta(0.01);

  for (TrackCollection::const_iterator it = (*theTkCollection).begin(); 
       it != (*theTkCollection).end(); 
       ++it){


    dpt  = fabs(pt - it->pt());
    dphi = phi - it->phi();
    while (dphi >= M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    dphi = fabs(dphi);
    deta = fabs(eta - it->eta());

//     if ((dpt < mdpt)
// 	&& (dphi < mdphi)
// 	&& (deta < mdeta)
// 	) {

    if ( 
	((dpt < mdpt) && (deta < mdeta) && (dphi < mdphi))  ||
	
	((dpt < mdpt) && (deta < eeta/2.) && (dphi < ephi/2.)) ||
	((deta < mdeta) && (dpt < ept/2.) && (dphi < ephi/2.)) ||
	((dphi < mdphi) && (dpt < ept/2.) && (deta < eeta/2.)) 
	
	) {
	 

      mdpt  = dpt;
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

// ------------ method called once each job just before starting event loop  ------------
void  BSDumpTracks::beginJob(const edm::EventSetup& setup) {

  edm::ESHandle<TrackAssociatorBase> theAssociator;
  setup.get<TrackAssociatorRecord>().get(fAssociatorLabel.c_str(), theAssociator);
  fAssociator = (TrackAssociatorBase*)theAssociator.product();

  gBSFile->cd();

}

// ------------ method called once each job just after ending the event loop  ------------
void  BSDumpTracks::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(BSDumpTracks);
