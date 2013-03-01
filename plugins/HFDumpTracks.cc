#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpTracks.h"
#include "HFDumpMuons.h"

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
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TSimpleTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaMuon.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TTrgObj.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpUtilities.hh"

#include <TFile.h>
#include <TH1.h>

// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;


// ----------------------------------------------------------------------
HFDumpTracks::HFDumpTracks(const edm::ParameterSet& iConfig):
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("ctfWithMaterialTracks"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("primaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fBeamSpotLabel(iConfig.getUntrackedParameter<InputTag>("beamSpotLabel", InputTag("offlineBeamSpot"))),
  fGenEventLabel(iConfig.getUntrackedParameter<InputTag>("generatorEventLabel", InputTag("source"))),
  fSimTracksLabel(iConfig.getUntrackedParameter<InputTag>("simTracksLabel", InputTag("famosSimHits"))),
  fAssociatorLabel(iConfig.getUntrackedParameter<InputTag>("associatorLabel", InputTag("TrackAssociatorByChi2"))), 
  fTrackingParticlesLabel(iConfig.getUntrackedParameter<InputTag>("trackingParticlesLabel", InputTag("trackingParticles"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)),
  fDumpSimpleTracks(iConfig.getUntrackedParameter<bool>("dumpSimpleTracks", true)),
  fDumpRecTracks(iConfig.getUntrackedParameter<bool>("dumpRecTracks", false)),
  fPropMuon(iConfig.getParameter<edm::ParameterSet>("propMuon"))
    {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTracks constructor  " << endl;
  cout << "---  tracksLabel:             " << fTracksLabel << endl;
  cout << "---  primaryVertexLabel:      " << fPrimaryVertexLabel << endl;
  cout << "---  beamSpotLabel:           " << fBeamSpotLabel << endl;
  cout << "---  muonsLabel:              " << fMuonsLabel << endl;
  cout << "---  generatorEventLabel:     " << fGenEventLabel << endl;
  cout << "---  simTracksLabel:          " << fSimTracksLabel << endl;
  cout << "---  associatorLabel:         " << fAssociatorLabel << endl;
  cout << "---  trackingParticlesLabel:  " << fTrackingParticlesLabel << endl;
  cout << "---  dumpSimpleTracks:        " << fDumpSimpleTracks << endl;
  cout << "---  dumpRecTracks:           " << fDumpRecTracks << endl;
  cout << "---  doTruthMatching:         " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS, 3 = TAna01Event
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDumpTracks::~HFDumpTracks() {
  
}


// ----------------------------------------------------------------------
void HFDumpTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (fVerbose > 0) cout << "==>HFDumpTracks> new event " << endl;
  // -- get the collection of RecoTracks 
  edm::Handle<edm::View<reco::Track> > tracksView;
  iEvent.getByLabel(fTracksLabel, tracksView);

  // -- get the primary vertex
  edm::Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFDumpTracks> No primary vertex collection found, skipping" << endl;
    return;
  }
  const reco::VertexCollection *vc = recoPrimaryVertexCollection.product();
  if (vc->size() == 0) {
    cout << "==>HFDumpTracks> No primary vertex found, skipping" << endl;
    return;
  }



  // -- get the beam spot
  const reco::BeamSpot *bs;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(fBeamSpotLabel, beamSpotHandle);

  if (beamSpotHandle.isValid()) {
    bs = beamSpotHandle.product();
  } else {
    bs = 0; 
    cout << "==>HFDumpTracks> No beam spot available from EventSetup" << endl;
  }


  // -- get the collection of muons and store their corresponding track indices
  vector<unsigned int> muonIndices, muonCollectionIndices; 
  Handle<MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  const reco::MuonCollection *mc = (MuonCollection*)hMuons.product();

  int muonIndex(0); 
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    TrackRef track = muon->innerTrack();
    muonIndices.push_back(track.index());
    muonCollectionIndices.push_back(muonIndex); 
    ++muonIndex; 
  }
 
  if (fVerbose > 0) cout << "==>HFDumpTracks> nMuons = " << hMuons->size() << endl;
  TH1D *h2 = (TH1D*)gHFFile->Get("h2");
  if (h2) h2->Fill(hMuons->size());
 
  if (fVerbose > 0) cout << "===> Tracks " << tracksView->size() << endl;

  TH1D *h1 = (TH1D*)gHFFile->Get("h1");
  if (h1) h1->Fill(tracksView->size());

  // -- fill association map between tracks and primary vertices
  //  tracksAndPv(iEvent); 
  //   if (fVerbose > 20) {
  //     for (unsigned int i = 0; i < tracksView->size(); ++i) {
  //       TrackBaseRef rTrackView(tracksView,i);
  //       Track trackView(*rTrackView);
  //       cout << Form("%4d: %d", i, fTrack2Pv[i]) << " pt = " << trackView.pt() <<endl; 
  //     }
  //   }
  
  int genIdx(-1); 
  for (unsigned int i = 0; i < tracksView->size(); ++i){    

    TrackBaseRef rTrackView(tracksView,i);
    const Track trackView(*rTrackView);

    // -- Muon?
    int mid = 0; 
    for (unsigned int im = 0; im < muonIndices.size(); ++im) {
      if (i == muonIndices[im]) {
	mid = 1; 
	break;
      }
    }
    
    // -- truth matching with deltaR comparision
    genIdx = -1; 
    if (3 == fDoTruthMatching) {    
      genIdx  = gHFEvent->getGenIndexWithDeltaR(trackView.pt(), trackView.eta(), trackView.phi(), trackView.charge()); 
    }
    
    // -- fill the tracks
    if (fDumpSimpleTracks) {
      TSimpleTrack *st = gHFEvent->addSimpleTrack();
      fillSimpleTrack(st, trackView, i, mid, genIdx, vc); 
    } 

    if (fDumpRecTracks) {
      TAnaTrack *at = gHFEvent->addRecTrack();
      fillAnaTrack(at, trackView, i, genIdx, vc, mc, bs); 
    }
  }


  ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  cleanupTruthMatching(tracksView, magfield);


  if (fVerbose > 0) {
    if (fDumpSimpleTracks) {
      for (int it = 0; it < gHFEvent->nSimpleTracks(); ++it) {
	gHFEvent->getSimpleTrack(it)->dump();
      }
    }

    if (fDumpRecTracks) {
      for (int it = 0; it < gHFEvent->nRecTracks(); ++it) {
	gHFEvent->getRecTrack(it)->dump();
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpTracks::beginJob() {
  gHFFile->cd();
  //  TH1D *h1 = new TH1D("h2", "h2", 20, 0., 20.);

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpTracks::endJob() {
}

void HFDumpTracks::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
	fPropMuon.init(iSetup);
} // beginRun()


//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTracks);
