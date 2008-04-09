#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpMuons.h"

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

using namespace std;
using namespace reco;
using namespace edm;

// ----------------------------------------------------------------------
HFDumpMuons::HFDumpMuons(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel"))

{
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpMuons constructor" << endl;
  cout << "--- muonsLabel:             " << fMuonsLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDumpMuons::~HFDumpMuons() {
  
}


// ----------------------------------------------------------------------
void HFDumpMuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // -- get the collection of muons
  Handle<reco::MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);

  if (fVerbose > 0) cout << "==>HFDumpMuons> nMuons = " << hMuons->size() << endl;

  std::vector<const reco::Track*> recTracks; 

  TAnaTrack *pTrack;

  // -- Look at muons
  const reco::Track* at = 0;
  const reco::Track* tt = 0;
  const reco::Track* ct = 0;

  int mcnt(0), type(0);
  for (reco::MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    
    mcnt++;

    // -- standalone muon
    type = 1;
    TrackRef staTrack = muon->standAloneMuon();
    at = &(*staTrack);

    pTrack   = gHFEvent->addSigTrack(); 
    pTrack->fMuType   = type;
    pTrack->fMCID     = at->charge()*-13; 
    pTrack->fGenIndex = -1; 
    pTrack->fQ        = at->charge();
    pTrack->fPlab.SetPtEtaPhi(at->pt(),
			      at->eta(),
			      at->phi()
			      );

    pTrack->fIndex  = (muon->track()).index();

    if (fVerbose > 0) {
      cout << "%%> Type-" << type <<  "  Muon #" << setw(2)  << mcnt << ": "; 
      pTrack->dump(); 
    }
    
    // -- tracker muon
    type = 2;
    TrackRef trkTrack = muon->track();
    tt = &(*trkTrack);

    pTrack   = gHFEvent->addSigTrack(); 
    pTrack->fMuType   = type;
    pTrack->fMCID     = tt->charge()*-13; 
    pTrack->fGenIndex = -1; 
    pTrack->fQ        = tt->charge();
    pTrack->fPlab.SetPtEtaPhi(tt->pt(),
			      tt->eta(),
			      tt->phi()
			      );

    pTrack->fIndex  = (muon->track()).index();

    if (fVerbose > 0) {
      cout << "%%>  Type-" << type <<  "  Muon #" << setw(2)  << mcnt << ": "; 
      pTrack->dump(); 
    }

    // -- combined muon
    type = 3;
    TrackRef glbTrack = muon->combinedMuon();
    ct = &(*glbTrack);


    pTrack   = gHFEvent->addSigTrack();
    pTrack->fMuType   = type;
    pTrack->fMCID     = ct->charge()*-13; 
    pTrack->fGenIndex = -1; 
    pTrack->fQ        = ct->charge();
    pTrack->fPlab.SetPtEtaPhi(ct->pt(),
			      ct->eta(),
			      ct->phi()
			      );

    pTrack->fIndex  = (muon->track()).index();

    if (fVerbose > 0) {
      cout << "%%>  Type-" << type <<  " Muon #" << setw(2)  << mcnt << ": "; 
      pTrack->dump(); 
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpMuons::beginJob(const edm::EventSetup& setup) {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpMuons::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(HFDumpMuons);
