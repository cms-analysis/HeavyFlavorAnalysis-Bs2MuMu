#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpStuff.h"

#include <iostream>

#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaMuon.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TTrgObj.hh"

#include <TVector3.h>

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace edm;
using namespace reco;


// ----------------------------------------------------------------------
HFDumpStuff::HFDumpStuff(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 1)),
  fGenEventScaleLabel(iConfig.getUntrackedParameter<string>("GenEventScaleLabel", string("genEventScale"))),
  fCandidates1Label(iConfig.getUntrackedParameter<string>("Candidates1Label", string("JPsiToMuMu"))),
  fCandidates2Label(iConfig.getUntrackedParameter<string>("Candidates2Label", string("JPsiToMuMu"))),
  fCandidates3Label(iConfig.getUntrackedParameter<string>("Candidates3Label", string("JPsiToMuMu"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fPrimaryVertexTracksLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexTracksLabel", InputTag("generalTracks")))
{
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpStuff constructor" << endl;
  cout << "---  verbose:                    " << fVerbose << endl;
  cout << "---  PrimaryVertexLabel:         " << fPrimaryVertexLabel << endl;
  cout << "---  PrimaryVertexTracksLabel:   " << fPrimaryVertexTracksLabel << endl;
  cout << "---  GenEventScaleLabel:         " << fGenEventScaleLabel.c_str() << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
HFDumpStuff::~HFDumpStuff() {
  
}


// ----------------------------------------------------------------------
void HFDumpStuff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  gHFEvent->fRunNumber   = iEvent.id().run();
  gHFEvent->fEventNumber = iEvent.id().event();
  gHFEvent->fBx          = iEvent.bunchCrossing();
  gHFEvent->fOrbit       = iEvent.orbitNumber();
  gHFEvent->fLumiSection = iEvent.luminosityBlock();

  const edm::Timestamp timeStamp = iEvent.time();
  unsigned int high = timeStamp.value() >> 32;       // seconds
  unsigned int low = 0xffffffff & timeStamp.value(); // microseconds  
  
  gHFEvent->fTimeLo      = low;
  gHFEvent->fTimeHi      = high;

  // -- Primary vertex
  int bestPV(-1), bestN(-1), cnt(0); 
  try {
    // -- get the collection of RecoTracks 
    edm::Handle<edm::View<Track> > tracksView;
    iEvent.getByLabel(fPrimaryVertexTracksLabel, tracksView);

    edm::Handle<VertexCollection> recoPrimaryVertexCollection;
    iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
    //    const VertexCollection vertices = *(recoPrimaryVertexCollection.product());

    int isFake(-1); 
    double cov[9]; 
    for (VertexCollection::const_iterator iv = recoPrimaryVertexCollection->begin(); iv != recoPrimaryVertexCollection->end(); ++iv) {
      TAnaVertex *pVtx = gHFEvent->addPV();
      ChiSquared chi2(iv->chi2(), iv->ndof());
      if (iv->isFake()) {
	isFake = 1; 
      } else {
	isFake = 0; 
      }
      
      int ntracks = iv->tracksSize(); 
      if (0 == isFake) {
	if (ntracks > bestN) {
	  bestN  = ntracks; 
	  bestPV = cnt;
	}
      }

      if (fVerbose > 1) {
	cout << "PV found: isFake = " << isFake << " with " << ntracks << " tracks" << endl;
      }
      pVtx->setInfo(chi2.value(), iv->ndof(), chi2.probability(), isFake, ntracks);
      pVtx->fPoint.SetXYZ(iv->x(),
			  iv->y(),
			  iv->z());

      for (int i = 0; i < 3; ++i) {
	for (int j = 0; j < 3; ++j) {
	  cov[i*3+j] = iv->covariance(i,j);
	}
      }
      pVtx->setCovXX(cov);

      Vertex::trackRef_iterator v1TrackIter;
      Vertex::trackRef_iterator v1TrackBegin = iv->tracks_begin();
      Vertex::trackRef_iterator v1TrackEnd   = iv->tracks_end();
      for (v1TrackIter = v1TrackBegin; v1TrackIter != v1TrackEnd; v1TrackIter++) {
	if (fVerbose > 10) {
	  cout << "vtx trk: " << " " << v1TrackIter->key()
	       << "  " << (**v1TrackIter).pt() 
	       << "  " << (**v1TrackIter).phi() 
	       << endl;
	}
	TrackBaseRef rTrackView(tracksView, v1TrackIter->key());
	Track track(*rTrackView);
	if (fVerbose > 10) {cout << " -> track " << "  " << track.pt() << "  " << track.phi() << endl;}
	pVtx->addTrack(v1TrackIter->key()); 
      }
      ++cnt; 
    }
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "==>HFDumpStuff> primaryVertex " << fPrimaryVertexLabel << " not found " << endl;
  } 

  if (bestPV > -1) {
    gHFEvent->fBestPV = bestPV;
  } else {
    gHFEvent->fBestPV = 0;
  }    
  if (fVerbose > 0) cout << "The best pV is at position: " << bestPV  << " and has " << bestN << " tracks" << endl;


  // -- pthat 
  gHFEvent->fPtHat = -1.; 
  try {
    edm::Handle<double> genEventScaleHandle;
    iEvent.getByLabel(fGenEventScaleLabel.c_str(), genEventScaleHandle);
    gHFEvent->fPtHat = *genEventScaleHandle;
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "genEventScale " << fGenEventScaleLabel.c_str() << " not found " << endl;
  }



}

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpStuff::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpStuff::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpStuff);
